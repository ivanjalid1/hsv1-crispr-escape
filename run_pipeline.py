#!/usr/bin/env python
"""Single entrypoint: fetch genomes -> extract guides -> score conservation -> report
-> (opt-in) robustness analysis.

Every stage caches its output, so re-running is cheap and re-runnable. Stages are
skipped when their output already exists and is newer than its inputs, unless
--force is given.

Examples
--------
    # fast smoke test on 5 genomes
    python run_pipeline.py --limit 5

    # full run plus the stage-5 robustness / denominator-sensitivity analysis
    python run_pipeline.py --robustness

    # full HSV-1 run
    python run_pipeline.py

    # include near-full-length "partial genome" clinical isolates
    python run_pipeline.py --include-partial

    # HSV-2 instead
    python run_pipeline.py --taxid 10310 --reference NC_001798

NCBI_EMAIL must be set in the environment. Nothing is hardcoded.
"""

from __future__ import annotations

import argparse
import json
import logging
import platform
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from src import conservation as stage_conservation  # noqa: E402
from src import extract_guides as stage_guides  # noqa: E402
from src import fetch_genomes as stage_fetch  # noqa: E402
from src import offtarget  # noqa: E402
from src import report as stage_report  # noqa: E402
from src import robustness as stage_robustness  # noqa: E402
from src.common import (  # noqa: E402
    DEFAULT_RUNLOG,
    MissingCredentialsError,
    ensure_dirs,
    setup_logging,
    utc_now_iso,
)

LOG = logging.getLogger("pipeline")


def _needs_rebuild(output: Path, inputs: list[Path], force: bool) -> bool:
    if force or not output.is_file():
        return True
    out_mtime = output.stat().st_mtime
    return any(p.is_file() and p.stat().st_mtime > out_mtime for p in inputs)


def build_parser() -> argparse.ArgumentParser:
    # The four stages intentionally share option names (--manifest, --candidates,
    # --conservation) with identical defaults so that paths stay consistent across
    # the pipeline. conflict_handler="resolve" lets the later stage re-register the
    # same flag instead of raising; because the defaults are defined once in
    # src/common.py, the resolved option is equivalent either way.
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        conflict_handler="resolve",
    )
    stage_fetch.add_arguments(parser)
    stage_guides.add_arguments(parser)
    stage_conservation.add_arguments(parser)
    stage_report.add_arguments(parser)
    stage_robustness.add_arguments(parser)

    parser.add_argument("--force", action="store_true",
                        help="Recompute every stage even if outputs exist.")
    parser.add_argument("--skip-fetch", action="store_true",
                        help="Reuse the existing manifest; do not contact NCBI for genomes.")
    parser.add_argument("-v", "--verbose", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    setup_logging(args.verbose)
    ensure_dirs()

    started = utc_now_iso()
    t0 = time.time()
    stages_run: list[str] = []

    LOG.info("=" * 78)
    LOG.info("HSV CRISPR-Cas9 conserved guide pipeline")
    LOG.info("Python %s on %s", platform.python_version(), platform.platform())
    LOG.info("=" * 78)

    # ---- Stage 1: genomes -------------------------------------------------------
    if args.skip_fetch:
        if not args.manifest.is_file():
            raise SystemExit(f"--skip-fetch given but no manifest at {args.manifest}")
        LOG.info("[1/4] fetch_genomes SKIPPED (--skip-fetch); reusing %s", args.manifest)
    elif _needs_rebuild(args.manifest, [], args.force):
        LOG.info("[1/4] fetch_genomes ...")
        stage_fetch.run(args)
        stages_run.append("fetch_genomes")
    else:
        LOG.info("[1/4] fetch_genomes: manifest present, reusing %s "
                 "(use --force to rebuild)", args.manifest)

    import pandas as pd
    manifest = pd.read_csv(args.manifest, sep="\t")
    LOG.info("      genomes in manifest: %d", len(manifest))
    if args.limit is not None and len(manifest) != args.limit:
        LOG.warning("      --limit %d was requested but the manifest holds %d genome(s). "
                    "Either the cached manifest is from a different run (use --force) or "
                    "NCBI simply returned fewer records.", args.limit, len(manifest))

    # ---- Stage 2: candidate guides ---------------------------------------------
    rebuild_guides = _needs_rebuild(args.candidates, [], args.force)
    if not rebuild_guides:
        # Guard against silently reusing a table built for a different reference or
        # a different gene set.
        cached = pd.read_csv(args.candidates, sep="\t")
        cached_ref = str(cached["reference_accession"].iloc[0]).split(".")[0] if len(cached) else ""
        cached_genes = {g for cell in cached["gene"].astype(str) for g in cell.split("|")}
        wanted_genes = {g.strip().upper() for g in args.genes.split(",") if g.strip()}
        if cached_ref != args.reference.split(".")[0] or not wanted_genes <= cached_genes:
            LOG.warning("      cached candidate table does not match --reference/--genes; "
                        "rebuilding.")
            rebuild_guides = True

    if rebuild_guides:
        LOG.info("[2/4] extract_guides ...")
        stage_guides.run(args)
        stages_run.append("extract_guides")
    else:
        LOG.info("[2/4] extract_guides: reusing %s", args.candidates)
    candidates = pd.read_csv(args.candidates, sep="\t")

    # ---- Stage 3: conservation --------------------------------------------------
    if _needs_rebuild(args.conservation, [args.manifest, args.candidates], args.force):
        LOG.info("[3/4] conservation ...")
        stage_conservation.run(args)
        stages_run.append("conservation")
    else:
        LOG.info("[3/4] conservation: reusing %s", args.conservation)

    # ---- Stage 4: report --------------------------------------------------------
    LOG.info("[4/4] report ...")
    ranked = stage_report.run(args)
    stages_run.append("report")

    # ---- Stage 5 (opt-in): robustness / sensitivity ------------------------------
    robustness_summary = None
    if args.robustness:
        LOG.info("[5/5] robustness (opt-in) ...")
        robustness_summary = stage_robustness.run(args)
        stages_run.append("robustness")
    else:
        LOG.info("[5/5] robustness: SKIPPED (pass --robustness to run it). The headline "
                 "conservation numbers above are measured over %d complete genomes only; "
                 "see results/robustness_report.md for what that denominator is worth.",
                 len(manifest))

    # ---- Run log ----------------------------------------------------------------
    import numpy
    import Bio
    run_log = {
        "started_utc": started,
        "finished_utc": utc_now_iso(),
        "elapsed_seconds": round(time.time() - t0, 2),
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "package_versions": {
            "biopython": Bio.__version__,
            "pandas": pd.__version__,
            "numpy": numpy.__version__,
        },
        "command_line": ["run_pipeline.py"] + (argv if argv is not None else sys.argv[1:]),
        "stages_executed": stages_run,
        "taxid": args.taxid,
        "reference_accession": str(candidates["reference_accession"].iloc[0])
        if len(candidates) else None,
        "target_genes": args.genes,
        "entrez_query": str(manifest["entrez_query"].iloc[0]) if len(manifest) else None,
        "retrieval_date_utc": str(manifest["retrieval_date_utc"].iloc[0])
        if len(manifest) else None,
        "n_genomes_used": int(len(manifest)),
        "accessions_used": sorted(manifest["accession"].astype(str).tolist()),
        "n_candidate_guides": int(len(candidates)),
        "n_ranked_guides": int(len(ranked)),
        "summary": stage_report.summarise(ranked),
        "offtarget_screening": {
            "implemented": offtarget.OFFTARGET_IMPLEMENTED,
            "caveat": offtarget.caveat(),
        },
        "robustness": robustness_summary,
    }
    DEFAULT_RUNLOG.parent.mkdir(parents=True, exist_ok=True)
    DEFAULT_RUNLOG.write_text(json.dumps(run_log, indent=2), encoding="utf-8")

    LOG.info("-" * 78)
    LOG.info("DONE in %.1fs. Ranked guides: %s", run_log["elapsed_seconds"], args.ranked)
    LOG.info("Run log (accessions, query, versions): %s", DEFAULT_RUNLOG)
    LOG.warning("REMINDER: %s", offtarget.caveat())
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except MissingCredentialsError as exc:
        print(f"\nERROR: {exc}\n", file=sys.stderr)
        raise SystemExit(2) from None
