"""Stage 3 -- alignment-free conservation scoring by exact 23-mer presence.

Method (deliberately alignment-free)
------------------------------------
A guide is "present" in a genome if its 23-mer (protospacer + PAM) occurs VERBATIM
in that genome on either strand. No multiple sequence alignment is performed or
required: a mismatch anywhere in the protospacer or in the GG of the PAM makes the
site a different site, and for the purpose of "will this guide cut this isolate?"
exact presence is the biologically relevant test. This removes every external
binary dependency (MAFFT/MUSCLE/Clustal) from the pipeline.

conservation_fraction = n_strains_present / n_strains_total

Algorithm and complexity
------------------------
The naive approach -- searching each of G guides in each of S genomes of length L --
is O(G x S x L) and is far too slow (thousands of guides x hundreds of 152 kb
genomes). Neither is it practical to store every 23-mer of every genome in a set
(S x L entries, tens of millions of strings, gigabytes of RAM).

Instead we invert the problem and scan each genome exactly once:

  1. Build two small dictionaries keyed by the guide 23-mers themselves:
       fwd[target_23mer]              -> guide_id   (guide on the genome's plus strand)
       rev[revcomp(target_23mer)]     -> guide_id   (guide on the genome's minus strand)
     Together these cover both strands for every guide regardless of which strand
     the guide came from in the reference. Memory is O(G), a few thousand entries.

  2. Every SpCas9 target 23-mer ends in NGG, so on the plus strand it must contain
     "GG" at offset 21, and its reverse complement must contain "CC" at offset 0.
     We therefore do not test all L positions: we jump between occurrences of "GG"
     and of "CC" using str.find, which runs in optimised C. Only those candidate
     positions produce a slice and a dict lookup.

Total work is O(S x L) with a small constant, and O(G) memory. In practice this
scores thousands of guides against hundreds of genomes in seconds.

Ambiguity caveat
----------------
A genome containing N (or any non-ACGT symbol) inside a target site cannot match
exactly and is scored as absent. This is a conservative, honest behaviour -- it can
only under-state conservation, never over-state it. Per-genome ambiguous base counts
are in the manifest, and --max-ambiguous-fraction can exclude poorly resolved
assemblies from the denominator.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.common import (  # noqa: E402
    DEFAULT_CANDIDATES,
    DEFAULT_CONSERVATION,
    DEFAULT_MANIFEST,
    PROJECT_ROOT,
    ensure_dirs,
    revcomp,
    setup_logging,
)

LOG = logging.getLogger("conservation")

CONSERVATION_COLUMNS = [
    "guide_id",
    "target_23mer",
    "n_strains_present",
    "n_strains_total",
    "conservation_fraction",
    "absent_in",
]


def load_genome(path: Path) -> str:
    lines = path.read_text(encoding="utf-8").splitlines()
    return "".join(ln.strip() for ln in lines if not ln.startswith(">")).upper()


def build_lookup(target_23mers: dict[str, str]) -> tuple[dict[str, str], dict[str, str]]:
    """target_23mers maps guide_id -> 23-mer. Returns (fwd_map, rev_map)."""
    fwd: dict[str, str] = {}
    rev: dict[str, str] = {}
    for guide_id, t23 in target_23mers.items():
        fwd[t23] = guide_id
        rev[revcomp(t23)] = guide_id
    return fwd, rev


def scan_genome(seq: str, fwd: dict[str, str], rev: dict[str, str],
                k: int = 23, pam_gg_offset: int = 21) -> set[str]:
    """Return the set of guide_ids whose 23-mer occurs in `seq` on either strand."""
    hits: set[str] = set()
    n = len(seq)

    # Plus-strand occurrences: the 23-mer ends in NGG, so "GG" sits at offset 21.
    pos = seq.find("GG")
    while pos != -1:
        start = pos - pam_gg_offset
        if start >= 0:
            gid = fwd.get(seq[start:start + k])
            if gid is not None:
                hits.add(gid)
        pos = seq.find("GG", pos + 1)

    # Minus-strand occurrences: the reverse complement of the 23-mer starts with "CC".
    pos = seq.find("CC")
    while pos != -1:
        if pos + k <= n:
            gid = rev.get(seq[pos:pos + k])
            if gid is not None:
                hits.add(gid)
        pos = seq.find("CC", pos + 1)

    return hits


def score_conservation(candidates: pd.DataFrame, manifest: pd.DataFrame,
                       max_ambiguous_fraction: float | None = None,
                       record_absences: int = 25) -> pd.DataFrame:
    used = manifest
    if max_ambiguous_fraction is not None:
        before = len(used)
        used = used[used["ambiguous_fraction"] <= max_ambiguous_fraction]
        if len(used) < before:
            LOG.warning("Excluded %d genome(s) with ambiguous_fraction > %g.",
                        before - len(used), max_ambiguous_fraction)
    used = used.sort_values("accession", kind="stable")

    if used.empty:
        raise SystemExit("No genomes left after filtering; cannot score conservation.")

    target_map = dict(zip(candidates["guide_id"], candidates["target_23mer"]))
    fwd, rev = build_lookup(target_map)
    LOG.info("Scoring %d guides against %d genomes (%d lookup keys).",
             len(target_map), len(used), len(fwd) + len(rev))

    present_counts = {gid: 0 for gid in target_map}
    absent_lists: dict[str, list[str]] = {gid: [] for gid in target_map}

    for i, row in enumerate(used.itertuples(index=False), start=1):
        path = PROJECT_ROOT / str(row.fasta_path)
        if not path.is_file():
            raise SystemExit(f"Manifest references a missing FASTA: {path}")
        seq = load_genome(path)
        hits = scan_genome(seq, fwd, rev)
        for gid in target_map:
            if gid in hits:
                present_counts[gid] += 1
            elif len(absent_lists[gid]) < record_absences:
                absent_lists[gid].append(row.accession)
        if i % 25 == 0 or i == len(used):
            LOG.info("  scanned %d/%d genomes", i, len(used))

    total = len(used)
    rows = []
    for gid, t23 in target_map.items():
        present = present_counts[gid]
        absent = absent_lists[gid]
        rows.append({
            "guide_id": gid,
            "target_23mer": t23,
            "n_strains_present": present,
            "n_strains_total": total,
            "conservation_fraction": round(present / total, 6),
            "absent_in": ";".join(absent) if present < total else "",
        })

    df = pd.DataFrame(rows, columns=CONSERVATION_COLUMNS)
    return df.sort_values("guide_id", kind="stable").reset_index(drop=True)


def sanity_check_reference(candidates: pd.DataFrame, conservation: pd.DataFrame,
                           manifest: pd.DataFrame) -> None:
    """If the reference genome is itself in the downloaded set, every guide must be
    found in it. A failure here means a coordinate/strand bug, not biology."""
    if candidates.empty:
        return
    ref = str(candidates["reference_accession"].iloc[0])
    ref_base = ref.split(".")[0]
    if ref_base not in set(manifest["accession_base"]):
        LOG.info("Reference %s is not part of the downloaded genome set "
                 "(it is not counted in the conservation denominator).", ref)
        return
    zero = conservation[conservation["n_strains_present"] == 0]
    if not zero.empty:
        LOG.error("SANITY CHECK FAILED: %d guides were not found in ANY genome even "
                  "though the reference %s is in the set.", len(zero), ref)
    else:
        LOG.info("Sanity check: reference %s is in the genome set and no guide has "
                 "zero hits.", ref)


def add_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--candidates", type=Path, default=DEFAULT_CANDIDATES)
    parser.add_argument("--conservation", type=Path, default=DEFAULT_CONSERVATION)
    parser.add_argument("--max-ambiguous-fraction", type=float, default=None,
                        help="Exclude genomes whose fraction of non-ACGT bases exceeds "
                             "this value (default: keep all, but see the caveat in the "
                             "module docstring).")


def run(args: argparse.Namespace) -> pd.DataFrame:
    ensure_dirs()
    if not args.manifest.is_file():
        raise SystemExit(f"Manifest not found: {args.manifest} (run fetch_genomes first)")
    if not args.candidates.is_file():
        raise SystemExit(f"Candidates not found: {args.candidates} (run extract_guides first)")

    manifest = pd.read_csv(args.manifest, sep="\t")
    candidates = pd.read_csv(args.candidates, sep="\t")
    if candidates.empty:
        raise SystemExit("Candidate guide table is empty.")

    df = score_conservation(candidates, manifest, args.max_ambiguous_fraction)
    sanity_check_reference(candidates, df, manifest)

    args.conservation.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.conservation, sep="\t", index=False)
    LOG.info("Wrote conservation table (%d guides) -> %s", len(df), args.conservation)

    perfect = int((df["conservation_fraction"] >= 1.0).sum())
    LOG.info("Perfectly conserved guides (present in all %d genomes): %d / %d",
             int(df["n_strains_total"].iloc[0]), perfect, len(df))
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
