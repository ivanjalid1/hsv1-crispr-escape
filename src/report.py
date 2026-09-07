"""Stage 4 -- join candidates with conservation scores, apply expression/QC flags,
rank, and write results/guides_ranked.tsv.

Flags applied here
------------------
has_polyT              TTTT or longer inside the protospacer. RNA polymerase III
                       (U6/H1 promoters) terminates at a run of >=4 T, so such a
                       guide is transcribed as a truncated, non-functional sgRNA.
                       This is a hard exclusion criterion, not a soft penalty.
max_homopolymer_run    Longest single-base run in the protospacer. Runs of >=5 are
                       associated with poor synthesis/activity and are flagged.
gc_in_range            GC fraction of the protospacer within [--gc-min, --gc-max].
                       Very low GC gives weak binding, very high GC gives off-target
                       promiscuity and poor unwinding.
passes_filters         Conservative AND of: conservation >= --min-conservation,
                       not has_polyT, max_homopolymer_run <= --max-homopolymer,
                       and gc_in_range.

Ranking is deterministic: perfectly conserved, expression-competent, GC-balanced
guides first; ties broken by gene name and reference coordinate so that the same
inputs always produce byte-identical output.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import pandas as pd

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.common import (  # noqa: E402
    DEFAULT_CANDIDATES,
    DEFAULT_CONSERVATION,
    DEFAULT_RANKED,
    RESULTS_DIR,
    ensure_dirs,
    setup_logging,
)

LOG = logging.getLogger("report")

REPORT_COLUMNS = [
    "rank",
    "guide_id",
    "gene",
    "gene_product",
    "strand",
    "position_in_reference",
    "ref_start",
    "ref_end",
    "cut_site_ref",
    "protospacer",
    "pam",
    "target_23mer",
    "gc_content",
    "gc_in_range",
    "max_homopolymer_run",
    "has_polyT",
    "conservation_fraction",
    "n_strains_present",
    "n_strains_total",
    "passes_filters",
    "n_reference_copies",
    "ref_all_positions",
    "absent_in",
    "reference_accession",
]


def build_report(candidates: pd.DataFrame, conservation: pd.DataFrame,
                 gc_min: float, gc_max: float, max_homopolymer: int,
                 min_conservation: float) -> pd.DataFrame:
    merged = candidates.merge(
        conservation.drop(columns=["target_23mer"]), on="guide_id", how="left",
        validate="one_to_one",
    )
    if merged["n_strains_present"].isna().any():
        missing = int(merged["n_strains_present"].isna().sum())
        raise SystemExit(f"{missing} guides have no conservation score; tables are out of sync.")

    merged["gc_content"] = merged["gc_fraction"].astype(float).round(4)
    merged["gc_in_range"] = merged["gc_content"].between(gc_min, gc_max)
    merged["position_in_reference"] = merged["ref_start"].astype(int)
    merged["absent_in"] = merged["absent_in"].fillna("")
    merged["has_polyT"] = merged["has_polyT"].astype(bool)
    merged["passes_filters"] = (
        (merged["conservation_fraction"] >= min_conservation)
        & (~merged["has_polyT"])
        & (merged["max_homopolymer_run"] <= max_homopolymer)
        & merged["gc_in_range"]
    )

    merged["_gc_center_dist"] = (merged["gc_content"] - (gc_min + gc_max) / 2.0).abs()
    merged = merged.sort_values(
        by=["conservation_fraction", "passes_filters", "has_polyT",
            "max_homopolymer_run", "_gc_center_dist", "gene", "ref_start", "strand"],
        ascending=[False, False, True, True, True, True, True, True],
        kind="stable",
    ).reset_index(drop=True)
    merged["rank"] = merged.index + 1

    for col in REPORT_COLUMNS:
        if col not in merged.columns:
            merged[col] = ""
    return merged[REPORT_COLUMNS]


def summarise(df: pd.DataFrame) -> dict:
    total_strains = int(df["n_strains_total"].iloc[0]) if len(df) else 0
    frac = df["conservation_fraction"]
    bins = {
        "== 1.00 (perfect)": int((frac >= 1.0).sum()),
        ">= 0.99": int((frac >= 0.99).sum()),
        ">= 0.95": int((frac >= 0.95).sum()),
        ">= 0.90": int((frac >= 0.90).sum()),
        ">= 0.50": int((frac >= 0.50).sum()),
        "== 0.00": int((frac <= 0.0).sum()),
    }
    per_gene = (
        df[df["conservation_fraction"] >= 1.0]["gene"].value_counts().sort_index().to_dict()
    )
    return {
        "n_guides": len(df),
        "n_strains_total": total_strains,
        "conservation_cumulative_counts": bins,
        "conservation_mean": round(float(frac.mean()), 6) if len(df) else 0.0,
        "conservation_median": round(float(frac.median()), 6) if len(df) else 0.0,
        "conservation_min": round(float(frac.min()), 6) if len(df) else 0.0,
        "n_has_polyT": int(df["has_polyT"].sum()),
        "n_gc_out_of_range": int((~df["gc_in_range"]).sum()),
        "n_passing_all_filters": int(df["passes_filters"].sum()),
        "perfectly_conserved_per_gene": per_gene,
    }


def add_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--candidates", type=Path, default=DEFAULT_CANDIDATES)
    parser.add_argument("--conservation", type=Path, default=DEFAULT_CONSERVATION)
    parser.add_argument("--ranked", type=Path, default=DEFAULT_RANKED)
    parser.add_argument("--gc-min", type=float, default=0.35)
    parser.add_argument("--gc-max", type=float, default=0.75)
    parser.add_argument("--max-homopolymer", type=int, default=4,
                        help="Maximum tolerated single-base run in the protospacer.")
    parser.add_argument("--min-conservation", type=float, default=1.0,
                        help="Conservation fraction required for passes_filters "
                             "(default 1.0 = present in every genome).")


def run(args: argparse.Namespace) -> pd.DataFrame:
    ensure_dirs()
    candidates = pd.read_csv(args.candidates, sep="\t")
    conservation = pd.read_csv(args.conservation, sep="\t")
    df = build_report(candidates, conservation, args.gc_min, args.gc_max,
                      args.max_homopolymer, args.min_conservation)

    args.ranked.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.ranked, sep="\t", index=False)
    LOG.info("Wrote ranked guide table (%d rows) -> %s", len(df), args.ranked)

    summary = summarise(df)
    (RESULTS_DIR / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    LOG.info("Conservation distribution over %d genomes:", summary["n_strains_total"])
    for label, count in summary["conservation_cumulative_counts"].items():
        LOG.info("  %-18s %5d / %d guides", label, count, summary["n_guides"])
    LOG.info("poly-T (TTTT) guides: %d | GC out of range: %d | passing all filters: %d",
             summary["n_has_polyT"], summary["n_gc_out_of_range"],
             summary["n_passing_all_filters"])
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
