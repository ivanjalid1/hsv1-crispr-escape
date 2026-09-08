"""Stage 3 -- alignment-free conservation scoring by exact target-site presence.

Method (deliberately alignment-free)
------------------------------------
A guide is "present" in a genome if its target site (protospacer + PAM; 23 nt for
SpCas9, 26-27 nt for SaCas9) occurs VERBATIM in that genome on either strand. No
multiple sequence alignment is performed or
required: a mismatch anywhere in the protospacer or in the constrained positions of
the PAM makes the
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

  1. Build two small dictionaries keyed by the guide target sites themselves
     (23 nt for SpCas9, 26-27 nt for SaCas9):
       fwd[target_23mer]              -> guide_id   (guide on the genome's plus strand)
       rev[revcomp(target_23mer)]     -> guide_id   (guide on the genome's minus strand)
     Together these cover both strands for every guide regardless of which strand
     the guide came from in the reference. Memory is O(G), a few thousand entries.

  2. Every target site of a given nuclease carries the LITERAL positions of its PAM
     pattern at fixed offsets: an SpCas9 (NGG) 23-mer always has "GG" at offset 21,
     and a SaCas9 (NNGRRT) 27-mer always has "G" at offset 23 and "T" at offset 26.
     We therefore do not test all L positions: we jump between occurrences of one
     such literal anchor using str.find, which runs in optimised C. Only those
     candidate positions produce a slice and a dict lookup. The anchors come from
     src/nuclease.py, so no PAM literal is hardcoded here; when several anchors are
     equally long the scanner picks, per sequence, the one that occurs least often
     (for a 68% GC herpesvirus genome the "T" of NNGRRT is ~3x rarer than the "G").
     A PAM with no literal position at all yields no anchors and the scanner falls
     back to testing every position -- slower, still exact.

     tests/test_nuclease.py proves set-identity between this scan and a naive
     both-strand substring search for every supported PAM, on randomised sequences
     and on the real HSV-1 reference genome.

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

from src.nuclease import SPCAS9, Nuclease, get_nuclease  # noqa: E402
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
                nuclease: Nuclease = SPCAS9) -> set[str]:
    """Return the set of guide_ids whose target site occurs in `seq` on either strand.

    Exact by construction: an occurrence of a target site necessarily contains the
    nuclease's literal PAM anchor at a fixed offset, so jumping between anchor
    occurrences cannot miss a hit. Falls back to an exhaustive positional scan when
    the PAM has no literal position at all.
    """
    hits: set[str] = set()
    n = len(seq)
    k = nuclease.target_length
    if n < k:
        return hits
    anchors = nuclease.scan_anchors

    def _sweep(table: dict[str, str], literal: str, offset: int) -> None:
        """Probe every position where `literal` occurs at `offset` inside the site."""
        if not literal:
            for start in range(0, n - k + 1):
                gid = table.get(seq[start:start + k])
                if gid is not None:
                    hits.add(gid)
            return
        pos = seq.find(literal)
        while pos != -1:
            start = pos - offset
            if start >= 0 and start + k <= n:
                gid = table.get(seq[start:start + k])
                if gid is not None:
                    hits.add(gid)
            pos = seq.find(literal, pos + 1)

    if not anchors:
        # No literal PAM position anywhere: exhaustive scan against both maps.
        _sweep(fwd, "", 0)
        _sweep(rev, "", 0)
        return hits

    # Plus-strand occurrences of the target site: anchor as-is, at its own offset.
    plus_literal, plus_offset = min(anchors, key=lambda a: (seq.count(a[0]), a[1]))
    _sweep(fwd, plus_literal, plus_offset)

    # Minus-strand occurrences: the rev map is keyed on revcomp(target), in which the
    # anchor appears reverse-complemented and mirrored to offset k - offset - len.
    rc_anchors = [(revcomp(lit), k - off - len(lit)) for lit, off in anchors]
    rc_literal, rc_offset = min(rc_anchors, key=lambda a: (seq.count(a[0]), a[1]))
    _sweep(rev, rc_literal, rc_offset)

    return hits


def scan_genome_naive(seq: str, fwd: dict[str, str], rev: dict[str, str],
                      nuclease: Nuclease = SPCAS9) -> set[str]:
    """Reference implementation of `scan_genome`: test every position, no anchors.

    Kept in the shipped module rather than only in the tests, because it is what the
    fast scanner is validated against and it documents exactly what "present" means.
    O(L) slices per strand; far too slow for production use.
    """
    hits: set[str] = set()
    k = nuclease.target_length
    for start in range(0, len(seq) - k + 1):
        site = seq[start:start + k]
        gid = fwd.get(site)
        if gid is not None:
            hits.add(gid)
        gid = rev.get(site)
        if gid is not None:
            hits.add(gid)
    return hits


def score_conservation(candidates: pd.DataFrame, manifest: pd.DataFrame,
                       max_ambiguous_fraction: float | None = None,
                       record_absences: int = 25,
                       nuclease: Nuclease = SPCAS9) -> pd.DataFrame:
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
    LOG.info("Scoring %d %s guides against %d genomes (%d lookup keys).",
             len(target_map), nuclease.label, len(used), len(fwd) + len(rev))

    present_counts = {gid: 0 for gid in target_map}
    absent_lists: dict[str, list[str]] = {gid: [] for gid in target_map}

    for i, row in enumerate(used.itertuples(index=False), start=1):
        path = PROJECT_ROOT / str(row.fasta_path)
        if not path.is_file():
            raise SystemExit(f"Manifest references a missing FASTA: {path}")
        seq = load_genome(path)
        hits = scan_genome(seq, fwd, rev, nuclease)
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


def nuclease_from_candidates(candidates: pd.DataFrame,
                             args: "argparse.Namespace | None" = None) -> Nuclease:
    """Recover the nuclease that produced a candidate table.

    Stage 2 writes `nuclease` and `pam_pattern` columns; they are authoritative, so a
    conservation run can never be scored under a different grammar than the one that
    enumerated the sites. Candidate tables written before those columns existed are
    SpCas9 by definition, and fall back to the CLI arguments.
    """
    if "pam_pattern" in candidates.columns and len(candidates):
        pam = str(candidates["pam_pattern"].iloc[0])
        name = str(candidates["nuclease"].iloc[0]).split("-")[0]
        spacer = len(str(candidates["protospacer"].iloc[0]))
        return get_nuclease(name, pam, spacer)
    if args is not None:
        return get_nuclease(getattr(args, "nuclease", SPCAS9.name),
                            getattr(args, "pam", None),
                            getattr(args, "spacer_length", None))
    return SPCAS9


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

    nuclease = nuclease_from_candidates(candidates, args)
    df = score_conservation(candidates, manifest, args.max_ambiguous_fraction,
                            nuclease=nuclease)
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
