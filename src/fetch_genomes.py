"""Stage 1 -- retrieve complete HSV genomes from NCBI Nucleotide and build a manifest.

Design notes
------------
* The taxon is configurable: HSV-1 is taxid 10298 (Human alphaherpesvirus 1),
  HSV-2 is taxid 10310 (Human alphaherpesvirus 2).
* Retrieval is deterministic: accessions are sorted before any --limit is applied,
  so `--limit 5` always yields the same five genomes for a given NCBI snapshot.
* Downloads are cached per-accession under data/raw/. Re-runs never re-download.
* The exact Entrez query string and the retrieval timestamp are written into the
  manifest (one column each, repeated per row) so a results table can always be
  traced back to the query that produced it.
* We NEVER pad or invent records. If NCBI returns fewer genomes than the query
  count suggested, the manifest contains exactly what was actually downloaded and
  the log says so loudly.

PRIVACY: the contact email comes from $NCBI_EMAIL only. Nothing is hardcoded.
"""

from __future__ import annotations

import argparse
import io
import json
import logging
import sys
from pathlib import Path

import pandas as pd

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.common import (  # noqa: E402
    DATA_DIR,
    DEFAULT_MANIFEST,
    META_DIR,
    RAW_DIR,
    RateLimiter,
    configure_entrez,
    ensure_dirs,
    entrez_rate_limit_interval,
    setup_logging,
    sha256_file,
    utc_now_iso,
    with_retries,
)

LOG = logging.getLogger("fetch")

KNOWN_TAXA = {
    10298: "Human alphaherpesvirus 1 (HSV-1)",
    10310: "Human alphaherpesvirus 2 (HSV-2)",
}

MANIFEST_COLUMNS = [
    "accession",
    "accession_base",
    "organism",
    "taxid",
    "length_bp",
    "n_ambiguous",
    "ambiguous_fraction",
    "strain",
    "isolate",
    "country",
    "collection_date",
    "host",
    "isolation_source",
    "definition",
    "completeness_label",
    "is_refseq",
    "fasta_path",
    "fasta_sha256",
    "entrez_query",
    "retrieval_date_utc",
]


# --------------------------------------------------------------------------------------
# Query construction
# --------------------------------------------------------------------------------------


def build_query(taxid: int, min_len: int, max_len: int, include_partial: bool) -> str:
    """Build the Entrez Nucleotide query.

    Rationale for each clause:
      txid<N>[Organism:exp]  -- the taxon and everything below it
      <min>:<max>[SLEN]      -- plausible full-length herpesvirus genome
      biomol_genomic[PROP]   -- genomic DNA, not mRNA/cDNA constructs
      NOT patent[PROP]       -- exclude patent sequences (often engineered)
      "complete genome"[Title] -- submitter asserts a complete genome

    NOTE: many high-quality clinical isolate assemblies (~152 kb) are deposited with
    the title "... , partial genome" because the terminal/internal repeats were not
    fully resolved. --include-partial drops the title clause and admits them. That
    roughly triples the HSV-1 dataset but weakens the "complete genome" claim, so it
    is OFF by default and is recorded per-record in the completeness_label column.
    """
    clauses = [
        f"txid{taxid}[Organism:exp]",
        f"{min_len}:{max_len}[SLEN]",
        "biomol_genomic[PROP]",
    ]
    query = " AND ".join(clauses)
    if not include_partial:
        query += ' AND "complete genome"[Title]'
    query += " NOT patent[PROP]"
    return query


# --------------------------------------------------------------------------------------
# Entrez interaction
# --------------------------------------------------------------------------------------


def esearch_all_uids(Entrez, query: str, limiter: RateLimiter,
                     page_size: int = 500) -> tuple[list[str], int]:
    """Return (uids, reported_total). Pages through the full result set."""

    def _count() -> int:
        handle = Entrez.esearch(db="nuccore", term=query, retmax=0)
        try:
            return int(Entrez.read(handle)["Count"])
        finally:
            handle.close()

    total = with_retries(_count, limiter=limiter, description="esearch (count)", logger=LOG)
    LOG.info("Entrez reports %d matching records for query: %s", total, query)

    uids: list[str] = []
    for start in range(0, total, page_size):
        def _page(start=start) -> list[str]:
            handle = Entrez.esearch(db="nuccore", term=query, retstart=start,
                                    retmax=page_size, idtype="acc")
            try:
                return list(Entrez.read(handle)["IdList"])
            finally:
                handle.close()

        page = with_retries(_page, limiter=limiter,
                            description=f"esearch (offset {start})", logger=LOG)
        if not page:
            LOG.warning("esearch returned an empty page at offset %d; stopping paging.", start)
            break
        uids.extend(page)

    # De-duplicate while keeping determinism, then sort.
    unique = sorted(set(uids))
    if len(unique) != total:
        LOG.warning(
            "NCBI reported %d records but paging yielded %d unique accessions. "
            "Proceeding with the %d actually returned.", total, len(unique), len(unique)
        )
    return unique, total


def fetch_fasta_batch(Entrez, accessions: list[str], limiter: RateLimiter) -> str:
    def _fetch() -> str:
        handle = Entrez.efetch(db="nuccore", id=",".join(accessions),
                               rettype="fasta", retmode="text")
        try:
            return handle.read()
        finally:
            handle.close()

    return with_retries(_fetch, limiter=limiter,
                        description=f"efetch fasta x{len(accessions)}", logger=LOG)


def fetch_metadata_batch(Entrez, accessions: list[str], limiter: RateLimiter) -> str:
    """Fetch only the GenBank header + source feature.

    Trick: requesting a 1 bp slice (seq_start=1, seq_stop=1) returns the full record
    header and the source feature qualifiers while transferring ~2 kB instead of the
    ~400 kB of a full annotated GenBank flat file.
    """

    def _fetch() -> str:
        handle = Entrez.efetch(db="nuccore", id=",".join(accessions), rettype="gb",
                               retmode="text", seq_start=1, seq_stop=1)
        try:
            return handle.read()
        finally:
            handle.close()

    return with_retries(_fetch, limiter=limiter,
                        description=f"efetch metadata x{len(accessions)}", logger=LOG)


# --------------------------------------------------------------------------------------
# Caching helpers
# --------------------------------------------------------------------------------------


def fasta_path_for(accession: str) -> Path:
    return RAW_DIR / f"{accession}.fasta"


def meta_path_for(accession: str) -> Path:
    return META_DIR / f"{accession}.gb"


def _cached_and_valid(path: Path, min_bytes: int = 100) -> bool:
    return path.is_file() and path.stat().st_size >= min_bytes


def split_fasta_records(text: str) -> dict[str, str]:
    """Split a multi-FASTA blob into {accession.version: record_text}."""
    out: dict[str, str] = {}
    current_id: str | None = None
    buf: list[str] = []
    for line in text.splitlines(keepends=True):
        if line.startswith(">"):
            if current_id is not None:
                out[current_id] = "".join(buf)
            header = line[1:].strip()
            current_id = header.split()[0] if header else None
            buf = [line]
        elif current_id is not None:
            buf.append(line)
    if current_id is not None:
        out[current_id] = "".join(buf)
    return out


def download_genomes(Entrez, accessions: list[str], limiter: RateLimiter,
                     batch_size: int, refresh: bool) -> list[str]:
    """Download FASTA + metadata for `accessions`, caching per accession.

    Returns the list of accessions for which BOTH files are present on disk.
    """
    missing_fasta = [a for a in accessions
                     if refresh or not _cached_and_valid(fasta_path_for(a), 1000)]
    missing_meta = [a for a in accessions
                    if refresh or not _cached_and_valid(meta_path_for(a))]

    LOG.info("Cache status: %d/%d genomes already downloaded, %d to fetch.",
             len(accessions) - len(missing_fasta), len(accessions), len(missing_fasta))

    for i in range(0, len(missing_fasta), batch_size):
        batch = missing_fasta[i:i + batch_size]
        LOG.info("Fetching FASTA %d-%d of %d ...", i + 1, i + len(batch), len(missing_fasta))
        blob = fetch_fasta_batch(Entrez, batch, limiter)
        records = split_fasta_records(blob)
        by_base = {rid.split(".")[0]: rid for rid in records}
        for acc in batch:
            key = acc if acc in records else by_base.get(acc.split(".")[0])
            if key is None:
                LOG.error("NCBI returned no FASTA for %s -- it will be omitted.", acc)
                continue
            fasta_path_for(acc).write_text(records[key], encoding="utf-8")

    meta_batch = max(batch_size, 50)
    for i in range(0, len(missing_meta), meta_batch):
        batch = missing_meta[i:i + meta_batch]
        LOG.info("Fetching metadata %d-%d of %d ...", i + 1, i + len(batch), len(missing_meta))
        blob = fetch_metadata_batch(Entrez, batch, limiter)
        # Records are separated by the GenBank terminator "//".
        chunks = [c for c in blob.split("\n//\n") if c.strip()]
        parsed: dict[str, str] = {}
        for chunk in chunks:
            text = chunk.rstrip() + "\n//\n"
            for line in text.splitlines():
                if line.startswith("VERSION"):
                    parsed[line.split()[1]] = text
                    break
        by_base = {k.split(".")[0]: k for k in parsed}
        for acc in batch:
            key = acc if acc in parsed else by_base.get(acc.split(".")[0])
            if key is None:
                LOG.error("NCBI returned no metadata for %s.", acc)
                continue
            meta_path_for(acc).write_text(parsed[key], encoding="utf-8")

    ok = [a for a in accessions
          if _cached_and_valid(fasta_path_for(a), 1000) and _cached_and_valid(meta_path_for(a))]
    dropped = sorted(set(accessions) - set(ok))
    if dropped:
        LOG.warning("%d accession(s) could not be fully retrieved and are excluded: %s",
                    len(dropped), ", ".join(dropped[:20]) + ("..." if len(dropped) > 20 else ""))
    return ok


# --------------------------------------------------------------------------------------
# Manifest construction
# --------------------------------------------------------------------------------------


def _first(qualifiers: dict, *keys: str) -> str:
    for key in keys:
        vals = qualifiers.get(key)
        if vals:
            return str(vals[0]).strip()
    return ""


def read_sequence(path: Path) -> str:
    """Read a single-record FASTA and return the uppercase sequence string."""
    lines = path.read_text(encoding="utf-8").splitlines()
    return "".join(ln.strip() for ln in lines if not ln.startswith(">")).upper()


def build_manifest(accessions: list[str], query: str, taxid: int,
                   retrieval_date: str) -> pd.DataFrame:
    from Bio import SeqIO

    rows = []
    for acc in accessions:
        fpath = fasta_path_for(acc)
        seq = read_sequence(fpath)
        n_ambig = len(seq) - sum(seq.count(b) for b in "ACGT")

        meta_text = meta_path_for(acc).read_text(encoding="utf-8")
        try:
            rec = next(SeqIO.parse(io.StringIO(meta_text), "genbank"))
        except (StopIteration, ValueError) as exc:
            LOG.error("Could not parse metadata for %s (%s); recording minimal row.", acc, exc)
            rec = None

        src_q: dict = {}
        definition = ""
        organism = ""
        if rec is not None:
            definition = (rec.description or "").strip()
            organism = rec.annotations.get("organism", "")
            for feat in rec.features:
                if feat.type == "source":
                    src_q = feat.qualifiers
                    break

        low_def = definition.lower()
        if "complete genome" in low_def:
            completeness = "complete"
        elif "partial genome" in low_def:
            completeness = "partial"
        else:
            completeness = "unstated"

        rows.append({
            "accession": acc,
            "accession_base": acc.split(".")[0],
            "organism": organism,
            "taxid": taxid,
            "length_bp": len(seq),
            "n_ambiguous": n_ambig,
            "ambiguous_fraction": round(n_ambig / len(seq), 8) if seq else 1.0,
            "strain": _first(src_q, "strain"),
            "isolate": _first(src_q, "isolate"),
            "country": _first(src_q, "geo_loc_name", "country"),
            "collection_date": _first(src_q, "collection_date"),
            "host": _first(src_q, "host"),
            "isolation_source": _first(src_q, "isolation_source"),
            "definition": definition,
            "completeness_label": completeness,
            "is_refseq": acc.startswith(("NC_", "NG_", "AC_")),
            "fasta_path": str(fpath.relative_to(DATA_DIR.parent)).replace("\\", "/"),
            "fasta_sha256": sha256_file(fpath),
            "entrez_query": query,
            "retrieval_date_utc": retrieval_date,
        })

    df = pd.DataFrame(rows, columns=MANIFEST_COLUMNS)
    return df.sort_values("accession", kind="stable").reset_index(drop=True)


# --------------------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------------------


def add_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--taxid", type=int, default=10298,
                        help="NCBI taxonomy id. 10298 = HSV-1 (default), 10310 = HSV-2.")
    parser.add_argument("--min-length", type=int, default=145_000,
                        help="Minimum genome length in bp (default 145000).")
    parser.add_argument("--max-length", type=int, default=160_000,
                        help="Maximum genome length in bp (default 160000).")
    parser.add_argument("--include-partial", action="store_true",
                        help="Also admit records titled 'partial genome' (near-full-length "
                             "clinical isolates with unresolved repeats). Off by default.")
    parser.add_argument("--query", default=None,
                        help="Override the whole Entrez query string.")
    parser.add_argument("--limit", type=int, default=None,
                        help="Use only the first N accessions (sorted) -- for smoke tests.")
    parser.add_argument("--batch-size", type=int, default=10,
                        help="Accessions per efetch request (default 10).")
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST,
                        help="Output manifest TSV path.")
    parser.add_argument("--refresh", action="store_true",
                        help="Ignore the download cache and re-fetch everything.")


def run(args: argparse.Namespace) -> pd.DataFrame:
    ensure_dirs()
    Entrez = configure_entrez()
    limiter = RateLimiter(entrez_rate_limit_interval())

    query = args.query or build_query(args.taxid, args.min_length, args.max_length,
                                      args.include_partial)
    retrieval_date = utc_now_iso()
    LOG.info("Taxon: %s (taxid %d)", KNOWN_TAXA.get(args.taxid, "user-specified"), args.taxid)

    accessions, reported_total = esearch_all_uids(Entrez, query, limiter)
    LOG.info("Retrieved %d unique accessions from esearch.", len(accessions))

    if args.limit is not None and args.limit < len(accessions):
        LOG.warning("--limit %d: using only the first %d of %d accessions (sorted). "
                    "Conservation fractions are NOT publication-grade in this mode.",
                    args.limit, args.limit, len(accessions))
        accessions = accessions[:args.limit]

    if not accessions:
        raise SystemExit("NCBI returned no accessions for this query. Nothing to do.")

    ok = download_genomes(Entrez, accessions, limiter, args.batch_size, args.refresh)
    if not ok:
        raise SystemExit("No genomes could be downloaded.")
    if len(ok) < len(accessions):
        LOG.warning("REAL COUNT: %d of %d requested genomes are usable.", len(ok), len(accessions))

    manifest = build_manifest(ok, query, args.taxid, retrieval_date)
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(args.manifest, sep="\t", index=False)
    LOG.info("Wrote manifest with %d genomes -> %s", len(manifest), args.manifest)

    summary = {
        "stage": "fetch_genomes",
        "taxid": args.taxid,
        "entrez_query": query,
        "retrieval_date_utc": retrieval_date,
        "esearch_reported_total": reported_total,
        "accessions_requested": len(accessions),
        "genomes_downloaded": len(manifest),
        "limit_applied": args.limit,
        "include_partial": bool(args.include_partial),
        "length_min_bp": int(manifest["length_bp"].min()),
        "length_max_bp": int(manifest["length_bp"].max()),
        "genomes_with_ambiguous_bases": int((manifest["n_ambiguous"] > 0).sum()),
    }
    (DATA_DIR / "fetch_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    add_arguments(parser)
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args(argv)
    setup_logging(args.verbose)
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
