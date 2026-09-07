"""Stage 2 -- enumerate every SpCas9 target site in the annotated target genes of a
reference genome.

Definitions used throughout the project
---------------------------------------
* protospacer : 20 nt of genomic DNA, 5'->3' on the strand the guide RNA matches.
* PAM         : the immediately 3' NGG trinucleotide on the same strand.
* target_23mer: protospacer + PAM, i.e. the 23 nt string whose EXACT presence in
                another genome (on either strand) defines conservation.

Coordinate/strand conventions
-----------------------------
Everything is reported against the plus strand of the reference in 1-based inclusive
coordinates (`ref_start`/`ref_end` span the whole 23-mer footprint). `strand` is the
strand the protospacer lies on. For a minus-strand guide the PAM is at the LOWER plus-
strand coordinates of that footprint.

Gene coordinates are parsed from the GenBank annotation of the reference -- no
coordinates are hardcoded anywhere in this file.

Spliced genes (RL2/ICP0 has three exons) are handled by iterating over each
contiguous genomic segment of the CDS separately, so every reported 23-mer is a
contiguous stretch of genomic DNA (which is what Cas9 actually cuts) that also lies
inside coding sequence.

Repeat-region genes (RL2 is present twice, in TRL and IRL) yield identical 23-mers
from both copies. These are collapsed into one guide record; every reference
position is retained in `ref_all_positions` and counted in `n_reference_copies`.
"""

from __future__ import annotations

import argparse
import logging
import sys
from collections import OrderedDict
from pathlib import Path

import pandas as pd

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.common import (  # noqa: E402
    DEFAULT_CANDIDATES,
    REF_DIR,
    RateLimiter,
    UNAMBIGUOUS,
    configure_entrez,
    ensure_dirs,
    entrez_rate_limit_interval,
    gc_fraction,
    has_polyt,
    max_homopolymer_run,
    revcomp,
    setup_logging,
    with_retries,
)

LOG = logging.getLogger("guides")

# Essential HSV-1 genes that are established or plausible anti-viral CRISPR targets.
DEFAULT_GENES = ["UL30", "UL19", "UL5", "UL52", "UL29", "RL2", "UL54"]

GENE_DESCRIPTIONS = {
    "UL30": "DNA polymerase catalytic subunit",
    "UL19": "major capsid protein VP5",
    "UL5": "helicase-primase helicase subunit",
    "UL52": "helicase-primase primase subunit",
    "UL29": "ICP8 single-stranded DNA binding protein",
    "RL2": "ICP0 ubiquitin E3 ligase (immediate early)",
    "UL54": "ICP27 multifunctional expression regulator",
}

DEFAULT_REFERENCE = "NC_001806"

CANDIDATE_COLUMNS = [
    "guide_id",
    "gene",
    "gene_product",
    "locus_tag",
    "strand",
    "ref_start",
    "ref_end",
    "cut_site_ref",
    "protospacer",
    "pam",
    "target_23mer",
    "gc_fraction",
    "max_homopolymer_run",
    "has_polyT",
    "n_reference_copies",
    "ref_all_positions",
    "reference_accession",
]


# --------------------------------------------------------------------------------------
# Reference retrieval (cached)
# --------------------------------------------------------------------------------------


def fetch_reference(accession: str, refresh: bool = False) -> Path:
    """Download the fully annotated GenBank record for `accession`, cached on disk."""
    ensure_dirs()
    path = REF_DIR / f"{accession}.gb"
    if path.is_file() and path.stat().st_size > 10_000 and not refresh:
        LOG.info("Using cached reference annotation: %s", path)
        return path

    Entrez = configure_entrez()
    limiter = RateLimiter(entrez_rate_limit_interval())

    def _fetch() -> str:
        handle = Entrez.efetch(db="nuccore", id=accession, rettype="gbwithparts",
                               retmode="text")
        try:
            return handle.read()
        finally:
            handle.close()

    LOG.info("Downloading reference annotation for %s ...", accession)
    text = with_retries(_fetch, limiter=limiter,
                        description=f"efetch reference {accession}", logger=LOG)
    if "LOCUS" not in text[:200]:
        raise RuntimeError(f"Unexpected GenBank payload for {accession}")
    path.write_text(text, encoding="utf-8")
    return path


# --------------------------------------------------------------------------------------
# Annotation parsing
# --------------------------------------------------------------------------------------


def collect_target_segments(record, genes: list[str], feature_type: str = "CDS"):
    """Return a list of (gene, product, locus_tag, seg_start0, seg_end, feature_strand).

    Segments are half-open plus-strand genomic intervals [seg_start0, seg_end).
    A spliced CDS contributes one segment per exon.
    """
    wanted = {g.upper() for g in genes}
    segments = []
    seen_genes = set()

    for feat in record.features:
        if feat.type != feature_type:
            continue
        names = [n.upper() for n in feat.qualifiers.get("gene", [])]
        names += [n.upper() for n in feat.qualifiers.get("gene_synonym", [])]
        hit = next((n for n in names if n in wanted), None)
        if hit is None:
            continue
        seen_genes.add(hit)
        product = (feat.qualifiers.get("product") or [""])[0]
        locus_tag = (feat.qualifiers.get("locus_tag") or [""])[0]
        strand = feat.location.strand
        for part in feat.location.parts:
            segments.append((hit, product, locus_tag,
                             int(part.start), int(part.end),
                             "+" if strand == 1 else "-"))

    missing = sorted(wanted - seen_genes)
    if missing:
        LOG.warning("Gene(s) not found as %s features in the reference annotation: %s",
                    feature_type, ", ".join(missing))
    segments.sort(key=lambda s: (s[0], s[3], s[4]))
    LOG.info("Found %d annotated %s segment(s) covering %d/%d requested genes.",
             len(segments), feature_type, len(seen_genes), len(wanted))
    return segments


# --------------------------------------------------------------------------------------
# Guide enumeration
# --------------------------------------------------------------------------------------


def _valid(seq: str) -> bool:
    return all(ch in UNAMBIGUOUS for ch in seq)


def enumerate_guides_in_segment(genome: str, start0: int, end: int,
                                protospacer_len: int = 20):
    """Yield (strand, ref_start1, ref_end1, protospacer, pam, target_23mer).

    Only 23-mers lying entirely inside [start0, end) are emitted.
    """
    total = protospacer_len + 3
    lo = max(0, start0)
    hi = min(len(genome), end)
    for i in range(lo, hi - total + 1):
        window = genome[i:i + total]
        if not _valid(window):
            continue
        # Plus-strand protospacer: window == 20 nt protospacer + NGG
        if window[protospacer_len + 1] == "G" and window[protospacer_len + 2] == "G":
            yield ("+", i + 1, i + total,
                   window[:protospacer_len], window[protospacer_len:], window)
        # Minus-strand protospacer: reverse complement of window ends in NGG,
        # which means the plus-strand window starts with CC.
        if window[0] == "C" and window[1] == "C":
            rc = revcomp(window)
            yield ("-", i + 1, i + total,
                   rc[:protospacer_len], rc[protospacer_len:], rc)


def build_guide_table(genome: str, segments, reference_accession: str,
                      protospacer_len: int = 20) -> pd.DataFrame:
    """Enumerate guides across all segments and collapse duplicate 23-mers."""
    collapsed: "OrderedDict[str, dict]" = OrderedDict()

    for gene, product, locus_tag, seg_start0, seg_end, _fstrand in segments:
        for strand, s1, e1, proto, pam, t23 in enumerate_guides_in_segment(
                genome, seg_start0, seg_end, protospacer_len):
            pos_label = f"{s1}-{e1}({strand})"
            existing = collapsed.get(t23)
            if existing is not None:
                existing["n_reference_copies"] += 1
                existing["_positions"].append(pos_label)
                if gene not in existing["_genes"]:
                    existing["_genes"].append(gene)
                continue

            # Cas9 cuts 3 bp 5' of the PAM. On the plus strand that is between
            # plus-coordinates (e1-5) and (e1-4); on the minus strand between
            # (s1+4) and (s1+5).
            cut = (e1 - 5) if strand == "+" else (s1 + 4)

            collapsed[t23] = {
                "guide_id": f"{gene}_{s1}{strand}",
                "gene": gene,
                "gene_product": product or GENE_DESCRIPTIONS.get(gene, ""),
                "locus_tag": locus_tag,
                "strand": strand,
                "ref_start": s1,
                "ref_end": e1,
                "cut_site_ref": cut,
                "protospacer": proto,
                "pam": pam,
                "target_23mer": t23,
                "gc_fraction": round(gc_fraction(proto), 4),
                "max_homopolymer_run": max_homopolymer_run(proto),
                "has_polyT": has_polyt(proto),
                "n_reference_copies": 1,
                "_positions": [pos_label],
                "_genes": [gene],
                "reference_accession": reference_accession,
            }

    rows = []
    for entry in collapsed.values():
        entry = dict(entry)
        entry["ref_all_positions"] = ";".join(entry.pop("_positions"))
        genes = entry.pop("_genes")
        if len(genes) > 1:
            entry["gene"] = "|".join(genes)
        rows.append(entry)

    df = pd.DataFrame(rows, columns=CANDIDATE_COLUMNS)
    return df.sort_values(["gene", "ref_start", "strand"], kind="stable").reset_index(drop=True)


# --------------------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------------------


def add_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--reference", default=DEFAULT_REFERENCE,
                        help=f"Reference accession (default {DEFAULT_REFERENCE}, "
                             "HSV-1 strain 17). For HSV-2 use e.g. NC_001798.")
    parser.add_argument("--genes", default=",".join(DEFAULT_GENES),
                        help="Comma-separated gene names to target "
                             f"(default: {','.join(DEFAULT_GENES)}).")
    parser.add_argument("--feature-type", default="CDS", choices=["CDS", "gene"],
                        help="Annotation feature type to take coordinates from "
                             "(default CDS: coding exons only).")
    parser.add_argument("--protospacer-length", type=int, default=20,
                        help="Protospacer length in nt (default 20 for SpCas9).")
    parser.add_argument("--candidates", type=Path, default=DEFAULT_CANDIDATES,
                        help="Output TSV of candidate guides.")
    parser.add_argument("--refresh-reference", action="store_true",
                        help="Re-download the reference GenBank record.")


def run(args: argparse.Namespace) -> pd.DataFrame:
    from Bio import SeqIO

    ensure_dirs()
    path = fetch_reference(args.reference, refresh=args.refresh_reference)
    record = SeqIO.read(path, "genbank")
    genome = str(record.seq).upper()
    LOG.info("Reference %s: %d bp, %d ambiguous base(s).",
             record.id, len(genome), len(genome) - sum(genome.count(b) for b in "ACGT"))

    genes = [g.strip() for g in args.genes.split(",") if g.strip()]
    segments = collect_target_segments(record, genes, args.feature_type)
    if not segments:
        raise SystemExit(f"No {args.feature_type} features matched genes: {genes}")

    df = build_guide_table(genome, segments, record.id, args.protospacer_length)
    args.candidates.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.candidates, sep="\t", index=False)

    LOG.info("Extracted %d unique candidate guides -> %s", len(df), args.candidates)
    if not df.empty:
        per_gene = df["gene"].value_counts().sort_index()
        for gene, count in per_gene.items():
            LOG.info("  %-12s %5d guides", gene, count)
    return df


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
