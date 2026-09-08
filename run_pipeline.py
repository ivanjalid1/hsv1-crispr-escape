#!/usr/bin/env python
"""Single entrypoint: fetch genomes -> extract guides -> score conservation -> report
-> (opt-in) robustness, SaCas9 benchmark, escape model and human off-target screen.

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

    # stage 6: the SaCas9 head-to-head against Amrani et al. 2024 (offline)
    python run_pipeline.py --skip-fetch --benchmark-sacas9

    # stage 7: the multiplex escape-probability model (offline, needs stages 1 and 5)
    python run_pipeline.py --skip-fetch --escape

    # stage 8: human GRCh38 off-target screening. SKIPPED unless the ~4 GB genome
    # cache already exists -- it is never downloaded implicitly. To build it once:
    #     python src/offtarget.py --stage fetch
    python run_pipeline.py --skip-fetch --offtarget

    # SaCas9 (21 nt spacer, NNGRRT PAM) instead of SpCas9. Outputs are namespaced
    # under results/sacas9/ so the SpCas9 results are never overwritten.
    python run_pipeline.py --nuclease sacas9 --skip-fetch

    # the more permissive SaCas9 PAM, and the 20 nt spacers of Amrani et al. 2024
    python run_pipeline.py --nuclease sacas9 --pam NNGRRN --skip-fetch
    python run_pipeline.py --nuclease sacas9 --spacer-length 20 --skip-fetch

NCBI_EMAIL must be set in the environment. Nothing is hardcoded.
"""

from __future__ import annotations

import argparse
import json
import logging
import platform
import shlex
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from src import benchmark_sacas9 as stage_benchmark  # noqa: E402
from src import conservation as stage_conservation  # noqa: E402
from src import escape as stage_escape  # noqa: E402
from src import extract_guides as stage_guides  # noqa: E402
from src import fetch_genomes as stage_fetch  # noqa: E402
from src import offtarget as stage_offtarget  # noqa: E402
from src import report as stage_report  # noqa: E402
from src import robustness as stage_robustness  # noqa: E402
from src.common import (  # noqa: E402
    DEFAULT_CANDIDATES,
    DEFAULT_CONSERVATION,
    DEFAULT_RANKED,
    DEFAULT_ROBUSTNESS_REPORT,
    MissingCredentialsError,
    ensure_dirs,
    namespaced,
    setup_logging,
    utc_now_iso,
)
from src.nuclease import get_nuclease  # noqa: E402

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
    stage_benchmark.add_arguments(parser)

    # Stages 7 and 8 are NOT registered with add_arguments() here, unlike stages 1-6.
    # Their option names collide with the pipeline's own on flags whose defaults
    # differ (--pam, --nuclease, --spacer-length, --stage), and conflict_handler
    # "resolve" would silently let the last registration win -- which would change
    # the default SpCas9 behaviour of stages 2-4. Both stages are nuclease-fixed by
    # construction (they exist to compare against Amrani et al.'s SaCas9 guides), so
    # instead they are opted into with one flag each and configured through a
    # pass-through argument string parsed by the module's own parser.
    parser.add_argument("--escape", action="store_true",
                        help="Stage 7: run the multiplex escape-probability model "
                             "(offline; needs the stage-1 manifest and genome cache).")
    parser.add_argument("--escape-args", default="",
                        help="Extra arguments passed verbatim to src/escape.py, "
                             "e.g. --escape-args=\"--repeat-model all-copies\".")
    parser.add_argument("--offtarget", action="store_true",
                        help="Stage 8: run human GRCh38 off-target screening. SKIPPED "
                             "with a message unless the ~4 GB genome cache already "
                             "exists; it is never downloaded implicitly.")
    parser.add_argument("--offtarget-download", action="store_true",
                        help="Permit stage 8 to download and encode GRCh38 (~4.1 GB on "
                             "disk) if the cache is absent. Explicit by design.")
    parser.add_argument("--offtarget-args", default="",
                        help="Extra arguments passed verbatim to src/offtarget.py, "
                             "e.g. --offtarget-args=\"--chrom 21\".")

    parser.add_argument("--force", action="store_true",
                        help="Recompute every stage even if outputs exist.")
    parser.add_argument("--skip-fetch", action="store_true",
                        help="Reuse the existing manifest; do not contact NCBI for genomes.")
    parser.add_argument("-v", "--verbose", action="store_true")
    return parser


def _substage_args(module, extra: str, verbose: bool) -> argparse.Namespace:
    """Build a namespace for an opt-in stage using that module's OWN parser.

    Defaults therefore come from the stage itself, not from the pipeline, which is
    what keeps stages 7 and 8 nuclease-fixed and keeps the pipeline's SpCas9 defaults
    untouched.
    """
    sub = argparse.ArgumentParser(prog=module.__name__, add_help=False)
    module.add_arguments(sub)
    ns = sub.parse_args(shlex.split(extra))
    ns.verbose = verbose
    return ns


def escape_preconditions(args) -> str | None:
    """Reason stage 7 cannot run, or None. Never downloads anything to find out."""
    if not args.manifest.is_file():
        return f"no manifest at {args.manifest} (run stage 1)"
    try:
        gb = stage_robustness.reference_genbank_path(args.reference)
    except SystemExit as exc:
        # reference_genbank_path exits when nothing is cached for the accession.
        # Inside the pipeline that is a reason to skip an opt-in stage, not a reason
        # to abort a run whose stages 1-4 have already succeeded.
        return str(exc)
    if not gb.is_file():
        return (f"the reference GenBank record {gb} is not cached "
                f"(run stage 2, or python src/extract_guides.py)")
    return None


def offtarget_preconditions(sub_args, allow_download: bool) -> str | None:
    """Reason stage 8 cannot run, or None.

    The GRCh38 cache is ~4.1 GB on disk. It is NEVER built as a side effect of a
    default pipeline run: absence of the cache is a skip, not a download, unless
    --offtarget-download says otherwise.
    """
    pool = Path(sub_args.pool)
    if not pool.is_file():
        return f"the stage-6 candidate pool {pool} is absent (run --benchmark-sacas9)"
    if allow_download:
        return None
    missing = [p for p in (stage_offtarget.CACHE_BIN, stage_offtarget.CACHE_INDEX)
               if not p.is_file()]
    if missing:
        return ("the GRCh38 genome cache is absent "
                f"({', '.join(str(p) for p in missing)}). Building it downloads "
                "~1.0 GB and occupies ~4.1 GB on disk, so it is never done "
                "implicitly. Run `python src/offtarget.py --stage fetch` once, or "
                "pass --offtarget-download.")
    if not sub_args.no_gtf and not stage_offtarget.ANNOTATION_CACHE.is_file():
        return (f"the Ensembl annotation cache {stage_offtarget.ANNOTATION_CACHE} is "
                "absent. Run `python src/offtarget.py --stage fetch`, pass "
                "--offtarget-download, or add --offtarget-args=\"--no-gtf\".")
    return None


def apply_nuclease_namespace(args: argparse.Namespace):
    """Resolve the nuclease and, for anything other than SpCas9, move the outputs
    into results/<tag>/ unless the user gave an explicit path.

    Outputs are namespaced rather than overwritten because the SpCas9 result set is
    published in the README: a SaCas9 run must never silently replace it.
    """
    nuclease = get_nuclease(args.nuclease, args.pam, args.spacer_length)
    tag = nuclease.tag
    defaults = {
        "candidates": DEFAULT_CANDIDATES,
        "conservation": DEFAULT_CONSERVATION,
        "ranked": DEFAULT_RANKED,
        "robustness_report": DEFAULT_ROBUSTNESS_REPORT,
    }
    for attr, default in defaults.items():
        if getattr(args, attr) == default:
            setattr(args, attr, namespaced(default, tag))
    args.ranked.parent.mkdir(parents=True, exist_ok=True)
    args.run_log = args.ranked.parent / "run_log.json"
    return nuclease


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    setup_logging(args.verbose)
    ensure_dirs()
    nuclease = apply_nuclease_namespace(args)

    started = utc_now_iso()
    t0 = time.time()
    stages_run: list[str] = []

    LOG.info("=" * 78)
    LOG.info("HSV CRISPR-Cas9 conserved guide pipeline")
    LOG.info("Python %s on %s", platform.python_version(), platform.platform())
    LOG.info("Nuclease: %s; target footprint %d nt; outputs -> %s",
             nuclease.label, nuclease.target_length, args.ranked.parent)
    LOG.info("=" * 78)

    # ---- Stage 1: genomes -------------------------------------------------------
    if args.skip_fetch:
        if not args.manifest.is_file():
            raise SystemExit(f"--skip-fetch given but no manifest at {args.manifest}")
        LOG.info("[1/8] fetch_genomes SKIPPED (--skip-fetch); reusing %s", args.manifest)
    elif _needs_rebuild(args.manifest, [], args.force):
        LOG.info("[1/8] fetch_genomes ...")
        stage_fetch.run(args)
        stages_run.append("fetch_genomes")
    else:
        LOG.info("[1/8] fetch_genomes: manifest present, reusing %s "
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
        cached_tag = (str(cached["nuclease"].iloc[0]) if "nuclease" in cached.columns
                      and len(cached) else "spcas9")
        if (cached_ref != args.reference.split(".")[0] or not wanted_genes <= cached_genes
                or cached_tag != nuclease.tag):
            LOG.warning("      cached candidate table does not match --reference/--genes; "
                        "rebuilding.")
            rebuild_guides = True

    if rebuild_guides:
        LOG.info("[2/8] extract_guides ...")
        stage_guides.run(args)
        stages_run.append("extract_guides")
    else:
        LOG.info("[2/8] extract_guides: reusing %s", args.candidates)
    candidates = pd.read_csv(args.candidates, sep="\t")

    # ---- Stage 3: conservation --------------------------------------------------
    if _needs_rebuild(args.conservation, [args.manifest, args.candidates], args.force):
        LOG.info("[3/8] conservation ...")
        stage_conservation.run(args)
        stages_run.append("conservation")
    else:
        LOG.info("[3/8] conservation: reusing %s", args.conservation)

    # ---- Stage 4: report --------------------------------------------------------
    LOG.info("[4/8] report ...")
    ranked = stage_report.run(args)
    stages_run.append("report")

    # ---- Stage 5 (opt-in): robustness / sensitivity ------------------------------
    robustness_summary = None
    if args.robustness:
        LOG.info("[5/8] robustness (opt-in) ...")
        robustness_summary = stage_robustness.run(args)
        stages_run.append("robustness")
    else:
        LOG.info("[5/8] robustness: SKIPPED (pass --robustness to run it). The headline "
                 "conservation numbers above are measured over %d complete genomes only; "
                 "see results/robustness_report.md for what that denominator is worth.",
                 len(manifest))

    # ---- Stage 6 (opt-in): SaCas9 head-to-head -----------------------------------
    benchmark_summary = None
    if args.benchmark_sacas9:
        LOG.info("[6/8] SaCas9 benchmark vs Amrani et al. 2024 (opt-in) ...")
        benchmark_summary = stage_benchmark.run(args)
        stages_run.append("benchmark_sacas9")
    else:
        LOG.info("[6/8] SaCas9 benchmark: SKIPPED (pass --benchmark-sacas9).")

    # ---- Stage 7 (opt-in): multiplex escape-probability model ---------------------
    escape_summary = None
    escape_skipped = None
    if args.escape:
        escape_skipped = escape_preconditions(args)
        if escape_skipped:
            LOG.warning("[7/8] escape model: SKIPPED -- %s", escape_skipped)
        else:
            LOG.info("[7/8] escape-probability model (opt-in) ...")
            esc_args = _substage_args(stage_escape, args.escape_args, args.verbose)
            esc_args.manifest = args.manifest
            escape_summary = stage_escape.run(esc_args)
            stages_run.append("escape")
    else:
        LOG.info("[7/8] escape model: SKIPPED (pass --escape).")

    # ---- Stage 8 (opt-in): human GRCh38 off-target screening ---------------------
    offtarget_summary = None
    offtarget_skipped = None
    if args.offtarget:
        ot_args = _substage_args(stage_offtarget, args.offtarget_args, args.verbose)
        if ot_args.stage != "all":
            raise SystemExit(
                "--offtarget runs the full screen; --stage "
                f"{ot_args.stage!r} is only meaningful when src/offtarget.py is "
                "invoked directly (python src/offtarget.py --stage "
                f"{ot_args.stage}).")
        offtarget_skipped = offtarget_preconditions(ot_args, args.offtarget_download)
        if offtarget_skipped:
            LOG.warning("[8/8] off-target screening: SKIPPED -- %s", offtarget_skipped)
        else:
            LOG.info("[8/8] human GRCh38 off-target screening (opt-in) ...")
            offtarget_summary = stage_offtarget.run_screen(ot_args)
            stages_run.append("offtarget")
    else:
        LOG.info("[8/8] off-target screening: SKIPPED (pass --offtarget; needs the "
                 "~4 GB GRCh38 cache, which is never downloaded implicitly).")

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
        "nuclease": nuclease.as_dict(),
        "entrez_query": str(manifest["entrez_query"].iloc[0]) if len(manifest) else None,
        "retrieval_date_utc": str(manifest["retrieval_date_utc"].iloc[0])
        if len(manifest) else None,
        "n_genomes_used": int(len(manifest)),
        "accessions_used": sorted(manifest["accession"].astype(str).tolist()),
        "n_candidate_guides": int(len(candidates)),
        "n_ranked_guides": int(len(ranked)),
        "summary": stage_report.summarise(ranked),
        "offtarget_screening": {
            "implemented": stage_offtarget.OFFTARGET_IMPLEMENTED,
            "ran_in_this_invocation": offtarget_summary is not None,
            "skipped_because": offtarget_skipped,
            "caveat": (stage_offtarget.caveat() if offtarget_summary is not None
                       else stage_offtarget.unscreened_caveat()),
        },
        "robustness": robustness_summary,
        "sacas9_benchmark": benchmark_summary,
        "escape": escape_summary,
        "escape_skipped_because": escape_skipped,
        "offtarget": offtarget_summary,
    }
    args.run_log.parent.mkdir(parents=True, exist_ok=True)
    args.run_log.write_text(json.dumps(run_log, indent=2, default=str),
                            encoding="utf-8")

    LOG.info("-" * 78)
    LOG.info("DONE in %.1fs. Ranked guides: %s", run_log["elapsed_seconds"], args.ranked)
    LOG.info("Run log (accessions, query, versions): %s", args.run_log)
    if offtarget_summary is not None:
        LOG.warning("REMINDER: %s", stage_offtarget.caveat())
    else:
        LOG.warning("REMINDER: %s", stage_offtarget.unscreened_caveat())
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except MissingCredentialsError as exc:
        print(f"\nERROR: {exc}\n", file=sys.stderr)
        raise SystemExit(2) from None
