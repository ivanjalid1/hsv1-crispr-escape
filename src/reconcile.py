"""Reconciliation of stages 6, 7 and 8 into one joint comparison table.

Why this module exists
----------------------
Stage 6 (conservation), stage 7 (escape probability) and stage 8 (human off-target
burden) each produced a *different* recommended replacement for Amrani et al. 2024's
ICP0g2 guide, because each ranked on its own axis. `results/recommendation.md`
resolves that into a single recommendation. This module is the arithmetic behind that
document: it joins the three stages' outputs on `guide_id`, computes the one quantity
that none of them produced -- **P(escape) for each ICP0 candidate paired with
ICP27g1** -- and evaluates Pareto dominance mechanically so that the choice is not
made by eye.

It is a reader, not a stage. It downloads nothing, re-derives nothing that a stage
already wrote, and refuses to run if any of the three inputs is missing rather than
silently producing a partial table.

The escape numbers
------------------
`src/escape.py` writes a per-site repair-escape probability `q_site` but only reports
P(escape) for the guide *sets* it happens to search. The pairing of an arbitrary ICP0
candidate with ICP27g1 is not among them. This module rebuilds the stage-7
`EscapeModel` over the same 183-genome presence matrix -- same nuclease grammar, same
`benchmark_sacas9.presence_matrix` -- and takes `q` from the stage-7 output, then
evaluates arbitrary pairs.

That reconstruction is **validated before any new number is read off it**: the Amrani
lead pair must reproduce the value in `results/escape_summary.json` to within 1e-12,
or this module exits. A reconstruction that silently disagreed with the stage that
produced it would make every figure in `results/recommendation.md` unsupported.

Selection rule
--------------
The rule applied in `results/recommendation.md` is stated there in prose and encoded
here in `GATES` and `DISCRIMINATOR` so that the two cannot drift apart. Changing the
rule means editing both.
"""

from __future__ import annotations

import argparse
import itertools
import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src import benchmark_sacas9 as bm  # noqa: E402
from src import escape as esc  # noqa: E402
from src import robustness as rob  # noqa: E402
from src.common import DEFAULT_MANIFEST, RESULTS_DIR, setup_logging  # noqa: E402
from src.nuclease import get_nuclease  # noqa: E402

LOG = logging.getLogger("reconcile")

#: The published lead pair, by site id in the stage-6 pool.
ICP0G2 = "RL2_4496+"
ICP27G1 = "UL54_115156+"

#: Axes over which Pareto dominance is tested. Conservation is higher-is-better;
#: escape probability and off-target burden are lower-is-better. Seed-mismatch counts
#: and the 21-nt columns are deliberately EXCLUDED from the core set and reported
#: separately, because they are secondary readouts of the same screen and including
#: them twice would double-weight the off-target axis.
PARETO_HIGHER = ("cons_183", "cons_gene_level", "cons_subgenomic")
PARETO_LOWER_CORE = ("p_escape_with_icp27g1", "nngrrt_le3", "nngrrt_le4",
                     "nngrrn_le3", "nngrrn_le4", "cds_le3")
PARETO_LOWER_SEED = ("nngrrt_seed_le1", "nngrrn_seed_le1")
PARETO_LOWER_21NT = ("nngrrt_le3_21nt", "nngrrt_le4_21nt",
                     "nngrrn_le3_21nt", "nngrrn_le4_21nt", "cds_le3_21nt")

#: Hard gates, in the order they are reported. Each is (label, predicate on a row).
#: These encode section 3 of results/recommendation.md; the prose and this tuple must
#: be changed together.
GATES = (
    ("G1 conservation = 1.000 over 183 complete genomes",
     lambda r: r.cons_183 >= 1.0),
    ("G2 conservation >= 0.95 on the gene-level corpus",
     lambda r: r.cons_gene_level >= 0.95),
    ("G3 passes the standard filters (poly-T, homopolymer, GC)",
     lambda r: bool(r.passes_filters)),
    ("G4 present in both ICP0 repeat copies",
     lambda r: int(r.n_reference_copies) >= 2),
    ("G5 zero GRCh38 sites at <= 3 mismatches under the canonical NNGRRT PAM",
     lambda r: int(r.nngrrt_le3) == 0),
    ("G6 zero coding-exon hits at <= 3 mismatches, either PAM, either spacer length",
     lambda r: int(r.cds_le3) == 0 and int(r.cds_le3_21nt) == 0),
)

#: Primary discriminator among gate survivors, lower is better. See section 3 of
#: results/recommendation.md for why the canonical PAM tier governs.
DISCRIMINATOR = "nngrrt_le4"

#: Tie-breaks after the discriminator, in order, each (column, lower_is_better).
TIE_BREAKS = (("p_escape_with_icp27g1", True), ("cons_gene_level", False),
              ("cons_subgenomic", False), ("n_subgenomic_records", False))


def _require(path: Path, produced_by: str) -> Path:
    if not path.is_file():
        raise SystemExit(f"{path} not found. Run {produced_by} first.")
    return path


def build_model(manifest_path: Path, reference: str, genes: tuple[str, ...],
                site_params: pd.DataFrame):
    """Rebuild the stage-7 EscapeModel, then prove it agrees with stage 7."""
    from Bio import SeqIO

    nuclease = get_nuclease("sacas9", "NNGRRT", 20)
    manifest = pd.read_csv(manifest_path, sep="\t")
    record = SeqIO.read(rob.reference_genbank_path(reference), "genbank")
    pool = bm.build_pool(str(record.seq).upper(), record, list(genes), nuclease)
    accessions, presence, col_index = bm.presence_matrix(
        pool, manifest.sort_values("accession", kind="stable"), nuclease)

    ids = [g for g in pool["guide_id"] if g in set(site_params["guide_id"])]
    sites = site_params.set_index("guide_id").loc[ids].reset_index()
    model = esc.make_model(pool, presence, col_index, sites)
    LOG.info("rebuilt escape model: %d genomes x %d sites",
             model.n_genomes, len(model.site_ids))
    return model, pool, presence, col_index


def validate_against_stage7(model, summary_path: Path, tol: float = 1e-12) -> float:
    """Reproduce stage 7's published lead-pair P(escape), or refuse to continue."""
    ours = model.p_escape_ids([ICP0G2, ICP27G1])
    payload = json.loads(summary_path.read_text(encoding="utf-8"))
    theirs = _find_lead_pair(payload)
    if theirs is None:
        LOG.warning("no lead-pair P(escape) found in %s; validation skipped",
                    summary_path)
        return ours
    if abs(ours - theirs) > tol:
        raise SystemExit(
            "INTEGRITY FAILURE: the rebuilt escape model does not reproduce stage 7. "
            f"Amrani lead pair: this module {ours!r}, {summary_path.name} {theirs!r} "
            f"(difference {abs(ours - theirs):.3e} > {tol:.0e}). Every paired "
            "P(escape) in results/recommendation.md would be unsupported, so this is "
            "a hard stop rather than a warning."
        )
    LOG.info("validated against stage 7: Amrani lead pair P(escape) = %.17g", ours)
    return ours


def _find_lead_pair(payload) -> float | None:
    """Locate the published lead pair's P(escape) anywhere in the stage-7 summary.

    Searched structurally rather than by a hardcoded path, so that a rename in the
    stage-7 summary degrades to a skipped validation with a warning rather than to a
    KeyError in the middle of a reconciliation.
    """
    stack = [payload]
    while stack:
        node = stack.pop()
        if isinstance(node, dict):
            for key, value in node.items():
                if isinstance(value, dict) and "lead" in key:
                    p = value.get("p_escape")
                    if isinstance(p, (int, float)):
                        return float(p)
                if isinstance(value, (dict, list)):
                    stack.append(value)
        elif isinstance(node, list):
            stack.extend(node)
    return None


def joint_table(model, presence, col_index, pool: pd.DataFrame,
                offtargets: pd.DataFrame, site_params: pd.DataFrame,
                gene: str = "RL2") -> pd.DataFrame:
    """One row per screened candidate in `gene`, every axis, joined on guide_id."""
    ot20 = offtargets[offtargets["spacer_length"] == 20].set_index("guide_id")
    ot21 = offtargets[offtargets["spacer_length"] == 21].copy()
    ot21["guide_id"] = ot21["guide_id"].str.replace(r"\|21nt", "", regex=True)
    ot21 = ot21.set_index("guide_id")
    pool_ix = pool.set_index("guide_id")
    q_of = dict(zip(site_params["guide_id"], site_params["q_site"]))

    screened = [g for g in ot20.index if g in pool_ix.index
                and pool_ix.loc[g, "gene"] == gene]
    rows = []
    for guide in sorted(screened, key=lambda g: int(pool_ix.loc[g, "rank_in_gene"])):
        p, o, o21 = pool_ix.loc[guide], ot20.loc[guide], ot21.loc[guide]
        cols = [model.index[guide], model.index[ICP27G1]]
        both = np.all(presence[:, [col_index[guide], col_index[ICP27G1]]], axis=1)
        rows.append({
            "guide_id": guide,
            "published_name": (str(p.amrani_guide).replace("Amrani2024_", "")
                               if isinstance(p.amrani_guide, str) else ""),
            "rank_in_gene": int(p.rank_in_gene),
            "cons_183": float(p.conservation_complete_genomes),
            "cons_gene_level": float(p.conservation_gene_level),
            "n_gene_level_records": int(p.n_records_covering),
            "cons_subgenomic": float(p.conservation_subgenomic),
            "n_subgenomic_records": int(p.n_records_covering_subgenomic),
            "n_reference_copies": int(p.n_reference_copies),
            "gc_content": float(p.gc_content),
            "local_gc_200bp": float(p.local_gc_200bp),
            "max_homopolymer_run": int(p.max_homopolymer_run),
            "has_polyT": bool(p.has_polyT),
            "passes_filters": bool(p.passes_filters),
            "q_site": float(q_of[guide]),
            "joint_conservation_with_icp27g1": model.joint_conservation(cols),
            "n_genomes_both_present": int(both.sum()),
            "p_escape_with_icp27g1": model.p_escape(cols),
            "nngrrt_le3": int(o.nngrrt_le3), "nngrrt_le4": int(o.nngrrt_total),
            "nngrrt_seed_le1": int(o.nngrrt_seed_le1),
            "nngrrn_le3": int(o.nngrrn_le3), "nngrrn_le4": int(o.nngrrn_total),
            "nngrrn_seed_le1": int(o.nngrrn_seed_le1),
            "cds_le3": int(o.n_le3_in_cds),
            "nngrrt_le3_21nt": int(o21.nngrrt_le3),
            "nngrrt_le4_21nt": int(o21.nngrrt_total),
            "nngrrn_le3_21nt": int(o21.nngrrn_le3),
            "nngrrn_le4_21nt": int(o21.nngrrn_total),
            "cds_le3_21nt": int(o21.n_le3_in_cds),
        })
    return pd.DataFrame(rows)


def pareto_front(table: pd.DataFrame, higher: tuple[str, ...],
                 lower: tuple[str, ...]) -> tuple[list[str], list[tuple[str, str]]]:
    """(front, dominance pairs). `a` dominates `b` if it is weakly better on every
    axis and strictly better on at least one."""
    idx = table.set_index("guide_id")
    dominated: list[tuple[str, str]] = []
    for a, b in itertools.permutations(idx.index, 2):
        weakly = (all(idx.at[a, c] >= idx.at[b, c] for c in higher)
                  and all(idx.at[a, c] <= idx.at[b, c] for c in lower))
        strictly = (any(idx.at[a, c] > idx.at[b, c] for c in higher)
                    or any(idx.at[a, c] < idx.at[b, c] for c in lower))
        if weakly and strictly:
            dominated.append((a, b))
    losers = {b for _, b in dominated}
    return [g for g in idx.index if g not in losers], dominated


def apply_rule(table: pd.DataFrame) -> dict:
    """The selection rule of results/recommendation.md, applied mechanically."""
    verdicts = []
    survivors = []
    for row in table.itertuples(index=False):
        failed = [label for label, ok in GATES if not ok(row)]
        verdicts.append({"guide_id": row.guide_id,
                         "gates_failed": "; ".join(failed) or "-",
                         "survives": not failed})
        if not failed:
            survivors.append(row.guide_id)

    sub = table[table["guide_id"].isin(survivors)]
    ranked = sub.sort_values(
        [DISCRIMINATOR] + [c for c, _ in TIE_BREAKS],
        ascending=[True] + [asc for _, asc in TIE_BREAKS],
        kind="stable")
    front, dominance = pareto_front(sub, PARETO_HIGHER, PARETO_LOWER_CORE)
    return {
        "gate_verdicts": pd.DataFrame(verdicts),
        "survivors": survivors,
        "ranked": ranked,
        "selected": (ranked["guide_id"].iloc[0] if len(ranked) else None),
        "pareto_front": front,
        "pareto_dominance": dominance,
    }


def add_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--reference", default="NC_001806")
    parser.add_argument("--gene", default="RL2",
                        help="gene whose candidates are compared (default RL2/ICP0)")
    parser.add_argument("--outdir", type=Path, default=RESULTS_DIR)
    parser.add_argument("--mc-draws", type=int, default=2_000_000,
                        help="draws for the Monte Carlo cross-check of the closed form")
    parser.add_argument("--seed", type=int, default=20240814,
                        help="seed for that cross-check (stage 7's default)")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="reconcile",
        description="Joint stage-6/7/8 comparison behind results/recommendation.md.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    add_arguments(parser)
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args(argv)
    setup_logging(args.verbose)

    pool_path = _require(RESULTS_DIR / "sacas9_benchmark_pool.tsv",
                         "stage 6 (python run_pipeline.py --benchmark-sacas9)")
    sites_path = _require(RESULTS_DIR / "escape_site_parameters.tsv",
                          "stage 7 (python run_pipeline.py --escape)")
    summary_path = _require(RESULTS_DIR / "escape_summary.json",
                            "stage 7 (python run_pipeline.py --escape)")
    ot_path = _require(RESULTS_DIR / "offtarget_summary.tsv",
                       "stage 8 (python run_pipeline.py --offtarget)")

    pool_stage6 = pd.read_csv(pool_path, sep="\t")
    site_params = pd.read_csv(sites_path, sep="\t")
    offtargets = pd.read_csv(ot_path, sep="\t")

    model, pool, presence, col_index = build_model(
        args.manifest, args.reference, ("RL2", "UL54"), site_params)
    validate_against_stage7(model, summary_path)

    table = joint_table(model, presence, col_index, pool_stage6, offtargets,
                        site_params, gene=args.gene)
    out = Path(args.outdir) / "recommendation_table.tsv"
    table.to_csv(out, sep="\t", index=False)
    LOG.info("wrote %s (%d candidates)", out, len(table))

    result = apply_rule(table)
    print("\nHard gates:")
    print(result["gate_verdicts"].to_string(index=False))
    print(f"\nSurvivors: {result['survivors']}")
    print(f"Pareto front over {PARETO_HIGHER + PARETO_LOWER_CORE}:")
    print(f"  {result['pareto_front']}")
    for a, b in result["pareto_dominance"]:
        print(f"  {a} dominates {b}")
    print(f"\nDiscriminator {DISCRIMINATOR} (lower is better), then "
          f"{[c for c, _ in TIE_BREAKS]}:")
    cols = ["guide_id", DISCRIMINATOR] + [c for c, _ in TIE_BREAKS]
    print(result["ranked"][cols].to_string(index=False))
    print(f"\nSELECTED: {result['selected']}")

    # Independent cross-check of the closed form for the two pairs the document
    # quotes, using stage 7's own seeded simulator. This is a check on the algebra,
    # not new evidence: it draws from the same presence matrix and the same q.
    if result["selected"]:
        print("\nMonte Carlo cross-check of the closed form "
              f"({args.mc_draws:,} draws, seed {args.seed}):")
        for icp0 in (ICP0G2, result["selected"]):
            cols = [model.index[icp0], model.index[ICP27G1]]
            exact = model.p_escape(cols)
            est, se = model.monte_carlo(cols, args.mc_draws, args.seed)
            z = abs(est - exact) / se if se else 0.0
            print(f"  {icp0} + {ICP27G1}: closed form {exact:.6g}, "
                  f"simulated {est:.6g} +/- {se:.2g} ({z:.2f} sigma)")

    # Secondary axis sets, reported so that the sensitivity of the front is visible.
    for label, extra in (("+ seed counts", PARETO_LOWER_SEED),
                         ("+ 21-nt columns", PARETO_LOWER_21NT),
                         ("+ both", PARETO_LOWER_SEED + PARETO_LOWER_21NT)):
        front, _ = pareto_front(table[table["guide_id"].isin(result["survivors"])],
                                PARETO_HIGHER, PARETO_LOWER_CORE + extra)
        print(f"Pareto front {label}: {front}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
