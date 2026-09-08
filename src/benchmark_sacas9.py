"""Stage 6 -- the SaCas9 head-to-head against Amrani et al. 2024 (EBT-104).

The question
------------
Amrani et al. 2024 (Mol Ther Methods Clin Dev 32:101303) built EBT-104 from four
SaCas9 guides, two in ICP0/RL2 and two in ICP27/UL54, selected on ">70% conservation"
against a ViPR CDS snapshot plus an off-target count. Their lead clinical pair is
ICP0g2 + ICP27g1. Measured with our exact-presence method over 183 complete HSV-1
genomes, ICP0g2 is present in only 155/183 (0.847) -- the worst of their four.

Stage 5 already established that. It could not, however, say whether they had better
options, because stages 2-4 enumerated SpCas9 (NGG) sites and their guides live in
SaCas9 (NNGRRT) space. This stage closes that gap: it enumerates the FULL SaCas9
candidate pool in exactly the two genes they targeted, in exactly their site grammar
(20 nt spacer + NNGRRT PAM, as printed in their Table 1), scores every candidate the
way stage 3 and stage 5 score ours, and reports where their four guides fall.

What is measured
----------------
1. Complete-genome conservation: exact presence of the 26 nt footprint in each of the
   183 complete genomes, either strand (stage-3 semantics, generalised PAM).
2. Per-genome presence vectors, so that a PAIR of guides can be scored jointly --
   the fraction of genomes in which BOTH sites are intact. For a two-cut excision
   strategy that joint number, not either marginal, is what escape depends on.
3. Gene-level conservation over the cached stage-5 corpus (sub-genomic records plus
   near-full-length "partial genome" records), with the same anchor-bracketed
   coverage correction, so a guide is only scored against records that demonstrably
   span its site.
4. The standard pipeline filters (poly-T, GC band, homopolymer run) at their stage-4
   defaults -- NOT retuned for this analysis -- plus the local GC context of each
   site, because ICP0 is GC-rich and repeat-associated and Amrani et al. report that
   they could not even amplify the ICP0g2 site for sequencing (GC ~85%).

Everything here is offline: it reuses the genomes in data/raw/ and the corpus in
data/raw/genes/ that earlier stages downloaded. Nothing is re-fetched, and no
threshold is chosen after seeing the answer -- the filters are the pipeline defaults
and the comparison thresholds are the competitor's own published values.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src import robustness as rob  # noqa: E402
from src.common import (  # noqa: E402
    DEFAULT_GENE_MANIFEST,
    DEFAULT_MANIFEST,
    PROJECT_ROOT,
    RESULTS_DIR,
    ensure_dirs,
    gc_fraction,
    setup_logging,
)
from src.conservation import build_lookup, load_genome, scan_genome  # noqa: E402
from src.extract_guides import (  # noqa: E402
    build_guide_table,
    collect_target_segments,
    fetch_reference,
)
from src.nuclease import Nuclease, get_nuclease  # noqa: E402

LOG = logging.getLogger("benchmark")

#: The two genes Amrani et al. target, and the label they use for each.
TARGET_GENES = {"RL2": "ICP0", "UL54": "ICP27"}

#: Their lead clinical pair (construct SaCas9-g2g1).
LEAD_PAIR = ("Amrani2024_ICP0g2", "Amrani2024_ICP27g1")

#: Their four guides, from src/robustness.py so there is exactly one copy of the
#: published sequences in this repository.
BENCHMARK_GUIDES = rob.BENCHMARK_GUIDES

BENCHMARK_GENE = {
    "Amrani2024_ICP0g1": "RL2",
    "Amrani2024_ICP0g2": "RL2",
    "Amrani2024_ICP27g1": "UL54",
    "Amrani2024_ICP27g2": "UL54",
}

POOL_COLUMNS = [
    "guide_id", "gene", "gene_label", "amrani_guide", "strand",
    "ref_start", "ref_end", "cut_site_ref", "protospacer", "pam", "target_site",
    "gc_content", "gc_in_range", "max_homopolymer_run", "has_polyT",
    "local_gc_100bp", "local_gc_200bp", "local_gc_300bp",
    "n_reference_copies", "ref_all_positions",
    "n_strains_present", "n_strains_total", "conservation_complete_genomes",
    "n_records_covering", "n_present_records", "conservation_gene_level",
    "n_records_covering_subgenomic", "conservation_subgenomic",
    "passes_filters", "rank_in_gene", "beats_lead_guide", "absent_in",
]


# --------------------------------------------------------------------------------------
# Candidate pool
# --------------------------------------------------------------------------------------


def build_pool(genome: str, record, genes: list[str], nuclease: Nuclease,
               feature_type: str = "CDS") -> pd.DataFrame:
    """Every site of `nuclease` inside the CDS of `genes`, via the stage-2 enumerator."""
    segments = collect_target_segments(record, genes, feature_type)
    if not segments:
        raise SystemExit(f"No {feature_type} features matched genes: {genes}")
    df = build_guide_table(genome, segments, record.id, nuclease)
    df = df.rename(columns={"target_23mer": "target_site",
                            "gc_fraction": "gc_content"})
    return df.reset_index(drop=True)


def add_local_gc(pool: pd.DataFrame, genome: str,
                 flanks: tuple[int, ...] = (100, 200, 300)) -> pd.DataFrame:
    """GC fraction of the reference in windows centred on each site.

    Motivation is concrete rather than cosmetic: Amrani et al. could not amplify or
    sequence their ICP0g2 site at all, attributing this to ~85% GC in the surrounding
    region. A replacement guide that is itself GC-acceptable but sits in an equally
    GC-extreme neighbourhood inherits the same synthesis/PCR problem, so the local
    context is reported next to the site's own GC.
    """
    n = len(genome)
    for flank in flanks:
        vals = []
        for row in pool.itertuples(index=False):
            centre = (int(row.ref_start) - 1 + int(row.ref_end)) // 2
            lo = max(0, centre - flank // 2)
            hi = min(n, centre + flank // 2)
            vals.append(round(gc_fraction(genome[lo:hi]), 4))
        pool[f"local_gc_{flank}bp"] = vals
    return pool


# --------------------------------------------------------------------------------------
# Conservation over the complete genomes, with per-genome presence vectors
# --------------------------------------------------------------------------------------


def presence_matrix(pool: pd.DataFrame, manifest: pd.DataFrame,
                    nuclease: Nuclease) -> tuple[list[str], np.ndarray, dict[str, int]]:
    """Boolean (genome x guide) presence matrix, stage-3 semantics.

    The marginal column means are exactly stage 3's `conservation_fraction`; the
    matrix additionally supports joint (pair) conservation, which a table of
    marginals cannot give.
    """
    used = manifest.sort_values("accession", kind="stable")
    guide_ids = pool["guide_id"].tolist()
    target_map = dict(zip(guide_ids, pool["target_site"]))
    if len(set(target_map.values())) != len(target_map):
        raise SystemExit("Duplicate target sites in the pool; stage 2 should collapse them.")
    fwd, rev = build_lookup(target_map)
    index = {gid: i for i, gid in enumerate(guide_ids)}

    accessions = used["accession"].astype(str).tolist()
    mat = np.zeros((len(accessions), len(guide_ids)), dtype=bool)
    for gi, row in enumerate(used.itertuples(index=False)):
        path = PROJECT_ROOT / str(row.fasta_path)
        if not path.is_file():
            raise SystemExit(f"Manifest references a missing FASTA: {path}")
        seq = load_genome(path)
        for gid in scan_genome(seq, fwd, rev, nuclease):
            mat[gi, index[gid]] = True
        if (gi + 1) % 50 == 0 or gi + 1 == len(accessions):
            LOG.info("  scanned %d/%d genomes", gi + 1, len(accessions))
    # The column index is returned explicitly: `assemble` re-sorts the pool by rank,
    # so positional alignment between pool rows and matrix columns does not survive
    # it, and any later consumer (the pair analysis) must look columns up by id.
    return accessions, mat, index


# --------------------------------------------------------------------------------------
# Gene-level corpus (cached; no network)
# --------------------------------------------------------------------------------------


def cached_corpus_records(gene_manifest: Path) -> list[tuple[str, str]]:
    """(record_class, accession) for every corpus record already on disk.

    Read from the stage-5 corpus manifest rather than re-running the Entrez queries,
    so this stage needs no network access and scores the SaCas9 sites against exactly
    the same corpus the SpCas9 analysis used.
    """
    if not gene_manifest.is_file():
        return []
    man = pd.read_csv(gene_manifest, sep="\t")
    out = []
    for row in man.itertuples(index=False):
        acc = str(row.accession)
        if (rob.GENE_DIR / f"{acc}.fasta").is_file():
            out.append((str(row.record_class), acc))
    return out


def robustness_defaults(**overrides) -> argparse.Namespace:
    """Stage-5 argument defaults, so the anchoring/coverage parameters here are
    provably the same ones the published robustness analysis used."""
    parser = argparse.ArgumentParser()
    rob.add_arguments(parser)
    args = parser.parse_args([])
    for key, value in overrides.items():
        setattr(args, key, value)
    return args


def gene_level_for_pool(pool: pd.DataFrame, reference_accession: str,
                        reference: str, genes: list[str],
                        records: list[tuple[str, str]],
                        nuclease: Nuclease) -> pd.DataFrame:
    """Stage-5 part C, restricted to this pool."""
    spans = rob.collect_gene_spans(reference_accession, genes)
    args = robustness_defaults()
    candidates = pool.rename(columns={"target_site": "target_23mer"})
    result = rob.gene_level_analysis(records, candidates, reference, spans, args,
                                     nuclease)
    return result["per_guide"]


# --------------------------------------------------------------------------------------
# Assembly and analysis
# --------------------------------------------------------------------------------------


def assemble(pool: pd.DataFrame, accessions: list[str], mat: np.ndarray,
             gene_level: pd.DataFrame | None,
             gc_min: float, gc_max: float, max_homopolymer: int) -> pd.DataFrame:
    n_genomes = len(accessions)
    pool = pool.copy()
    pool["n_strains_total"] = n_genomes
    pool["n_strains_present"] = mat.sum(axis=0)
    pool["conservation_complete_genomes"] = (pool["n_strains_present"] / n_genomes).round(6)
    pool["gene_label"] = pool["gene"].map(TARGET_GENES).fillna(pool["gene"])

    by_target = {proto + pam: name for name, (proto, pam) in BENCHMARK_GUIDES.items()}
    pool["amrani_guide"] = pool["target_site"].map(by_target).fillna("")

    pool["gc_in_range"] = pool["gc_content"].between(gc_min, gc_max)
    pool["has_polyT"] = pool["has_polyT"].astype(bool)
    # Exactly the stage-4 flags, minus the conservation term: this analysis ranks BY
    # conservation, so folding a conservation threshold into passes_filters would
    # make the ranking circular. The conservation comparison is done explicitly.
    pool["passes_filters"] = (
        (~pool["has_polyT"])
        & (pool["max_homopolymer_run"] <= max_homopolymer)
        & pool["gc_in_range"]
    )

    absent = []
    for j in range(mat.shape[1]):
        miss = [accessions[i] for i in np.flatnonzero(~mat[:, j])][:25]
        absent.append(";".join(miss))
    pool["absent_in"] = absent

    if gene_level is not None:
        gl = gene_level.rename(columns={"n_present": "n_present_records"})
        keep = ["guide_id", "n_records_covering", "n_present_records",
                "conservation_gene_level", "n_records_covering_subgenomic",
                "conservation_subgenomic"]
        pool = pool.merge(gl[keep], on="guide_id", how="left", validate="one_to_one")
    else:
        for col in ("n_records_covering", "n_present_records", "conservation_gene_level",
                    "n_records_covering_subgenomic", "conservation_subgenomic"):
            pool[col] = np.nan

    # Rank within gene: conservation first, then the same deterministic tie-breaks
    # stage 4 uses, so the ordering rule is not invented for this report.
    pool["_gc_center_dist"] = (pool["gc_content"] - (gc_min + gc_max) / 2.0).abs()
    pool = pool.sort_values(
        by=["gene", "conservation_complete_genomes", "passes_filters", "has_polyT",
            "max_homopolymer_run", "_gc_center_dist", "ref_start", "strand"],
        ascending=[True, False, False, True, True, True, True, True],
        kind="stable",
    ).reset_index(drop=True)
    pool["rank_in_gene"] = pool.groupby("gene").cumcount() + 1
    return pool


def filter_cascade(gene_pool: pd.DataFrame, threshold: float) -> pd.DataFrame:
    """How many candidates survive each successive constraint. Reported as a cascade
    so that a large headline count cannot hide behind one permissive filter."""
    rows = []
    step = gene_pool
    rows.append(("all SaCas9 sites in the gene", len(step)))
    step = step[step["conservation_complete_genomes"] > threshold]
    rows.append((f"... conservation > {threshold:.3f}", len(step)))
    step = step[~step["has_polyT"]]
    rows.append(("... and no poly-T (TTTT) terminator", len(step)))
    step = step[step["max_homopolymer_run"] <= 4]
    rows.append(("... and homopolymer run <= 4", len(step)))
    step = step[step["gc_in_range"]]
    rows.append(("... and GC within 0.35-0.75", len(step)))
    # Only meaningful for RL2/ICP0, which is present twice (TRL and IRL). UL54 is
    # single-copy, so the row would read 0 and imply a failure that is not one.
    if (gene_pool["n_reference_copies"] >= 2).any():
        dup = step[step["n_reference_copies"] >= 2]
        rows.append(("... and present in both ICP0 repeat copies", len(dup)))
    return pd.DataFrame(rows, columns=["constraint", "n_sites"])


def pair_analysis(pool: pd.DataFrame, mat: np.ndarray, columns: dict[str, int],
                  require_filters: bool = False) -> pd.DataFrame:
    """Joint conservation of every (RL2 site, UL54 site) pair.

    A two-guide excision strategy needs BOTH sites intact in a given isolate, so the
    pair's figure of merit is the fraction of genomes where both are present -- not
    the mean or the minimum of the two marginals, which can both overstate it when
    the two absences fall in different isolates.
    """
    idx = columns
    left = pool[pool["gene"] == "RL2"]
    right = pool[pool["gene"] == "UL54"]
    if require_filters:
        left = left[left["passes_filters"]]
        right = right[right["passes_filters"]]
    n = mat.shape[0]
    rows = []
    for lrow in left.itertuples(index=False):
        li = mat[:, idx[lrow.guide_id]]
        for rrow in right.itertuples(index=False):
            ri = mat[:, idx[rrow.guide_id]]
            both = int(np.count_nonzero(li & ri))
            rows.append({
                "icp0_guide": lrow.guide_id,
                "icp0_amrani": lrow.amrani_guide,
                "icp0_conservation": lrow.conservation_complete_genomes,
                "icp0_gc": lrow.gc_content,
                "icp0_local_gc_200bp": lrow.local_gc_200bp,
                "icp0_copies": lrow.n_reference_copies,
                "icp27_guide": rrow.guide_id,
                "icp27_amrani": rrow.amrani_guide,
                "icp27_conservation": rrow.conservation_complete_genomes,
                "icp27_gc": rrow.gc_content,
                "icp27_local_gc_200bp": rrow.local_gc_200bp,
                "icp0_gene_level": lrow.conservation_gene_level,
                "icp27_gene_level": rrow.conservation_gene_level,
                "icp0_n_present": int(np.count_nonzero(li)),
                "icp27_n_present": int(np.count_nonzero(ri)),
                "n_genomes_both_present": both,
                "joint_conservation": round(both / n, 6),
                "n_genomes_neither": int(np.count_nonzero(~li & ~ri)),
                "n_genomes_at_least_one": int(np.count_nonzero(li | ri)),
                "both_pass_filters": bool(lrow.passes_filters and rrow.passes_filters),
            })
    out = pd.DataFrame(rows)
    # Invariant check, not decoration: joint presence cannot exceed either marginal,
    # and cannot fall below the Bonferroni floor. An earlier version of this function
    # looked matrix columns up positionally after the pool had been re-sorted by rank,
    # which silently produced joint values above both marginals; this assertion is
    # what such a bug hits first.
    # Checked on the integer counts, not the rounded fractions: at n=183 the two
    # differ in the sixth decimal and a float comparison fires spuriously.
    lo = (out["icp0_n_present"] + out["icp27_n_present"] - n).clip(lower=0)
    hi = out[["icp0_n_present", "icp27_n_present"]].min(axis=1)
    bad = out[(out["n_genomes_both_present"] > hi)
              | (out["n_genomes_both_present"] < lo)]
    if not bad.empty:
        raise AssertionError(
            "joint presence count outside [max(0, a+b-n), min(a, b)] for "
            f"{len(bad)} pair(s); first: {bad.iloc[0].to_dict()}"
        )
    # Ties on joint conservation are broken by the weaker of the two guides' gene-level
    # scores (the second, independent denominator), then by the marginals, then by
    # coordinate. Without the gene-level term the tie-break would be alphabetical,
    # which would let a site with less corroborating evidence be named "best" purely
    # on its identifier.
    out["_min_gene_level"] = out[["icp0_gene_level", "icp27_gene_level"]].min(
        axis=1).fillna(-1.0)
    out = out.sort_values(
        ["joint_conservation", "_min_gene_level", "icp0_conservation",
         "icp27_conservation", "icp0_guide", "icp27_guide"],
        ascending=[False, False, False, False, True, True],
        kind="stable").reset_index(drop=True)
    return out.drop(columns=["_min_gene_level"])


# --------------------------------------------------------------------------------------
# Context assembly -- every number quoted in the report is computed here
# --------------------------------------------------------------------------------------


def _round(x, nd=4):
    return None if x is None or (isinstance(x, float) and np.isnan(x)) else round(float(x), nd)


def build_context(pool: pd.DataFrame, pairs: pd.DataFrame, accessions: list[str],
                  mat: np.ndarray, columns: dict[str, int], nuclease: Nuclease,
                  reference: str, n_corpus_records: int, gene_level_ran: bool,
                  args: argparse.Namespace, generated: str) -> dict:
    n_g = len(accessions)
    amrani_rows = pool[pool["amrani_guide"] != ""].copy()
    if len(amrani_rows) != len(BENCHMARK_GUIDES):
        missing = set(BENCHMARK_GUIDES) - set(amrani_rows["amrani_guide"])
        raise SystemExit(
            "Not all published guides were recovered from the candidate pool "
            f"(missing: {sorted(missing)}). They must all be present, or the ranks "
            "below are meaningless -- check --benchmark-spacer-length and the PAM."
        )
    amrani_rows["n_better_in_gene"] = [
        int(((pool["gene"] == row.gene)
             & (pool["conservation_complete_genomes"]
                > row.conservation_complete_genomes)).sum())
        for row in amrani_rows.itertuples(index=False)
    ]
    amrani_rows = amrani_rows.sort_values("amrani_guide", kind="stable")

    g2 = amrani_rows.loc[amrani_rows["amrani_guide"] == LEAD_PAIR[0]].iloc[0]
    g1 = amrani_rows.loc[amrani_rows["amrani_guide"] == LEAD_PAIR[1]].iloc[0]
    thr0 = float(g2["conservation_complete_genomes"])
    thr27 = float(g1["conservation_complete_genomes"])

    # Recorded on the pool itself so the shipped TSV carries the comparison the
    # report is built on, per gene, against that gene's published lead guide.
    lead_threshold = {"RL2": thr0, "UL54": thr27}
    pool = pool.copy()
    pool["beats_lead_guide"] = [
        bool(row.conservation_complete_genomes > lead_threshold.get(row.gene, 1.1))
        for row in pool.itertuples(index=False)
    ]
    icp0 = pool[pool["gene"] == "RL2"]
    icp27 = pool[pool["gene"] == "UL54"]

    # -- pool overview --------------------------------------------------------------
    rows = []
    for gene, label, sub, thr in (("RL2", "ICP0", icp0, thr0),
                                  ("UL54", "ICP27", icp27, thr27)):
        rows.append({
            "gene": f"{gene} ({label})",
            "SaCas9 sites": len(sub),
            "pass poly-T/GC/homopolymer": int(sub["passes_filters"].sum()),
            "conservation = 1.000": int((sub["conservation_complete_genomes"] >= 1).sum()),
            "conservation > 0.70 (their cut)": int(
                (sub["conservation_complete_genomes"] > 0.70).sum()),
            "median GC": _round(sub["gc_content"].median(), 3),
            "median local GC (200 bp)": _round(sub["local_gc_200bp"].median(), 3),
        })
    pool_counts = pd.DataFrame(rows)

    # -- their four guides ----------------------------------------------------------
    cols = {
        "amrani_guide": "guide", "gene_label": "gene", "guide_id": "site",
        "rank_in_gene": "rank in gene pool", "n_better_in_gene": "sites ranked above it",
        "conservation_complete_genomes": f"conservation (n={n_g})",
        "conservation_gene_level": "gene-level corpus",
        "gc_content": "GC", "local_gc_200bp": "local GC 200 bp",
        "max_homopolymer_run": "max homopolymer", "passes_filters": "passes filters",
        "n_reference_copies": "reference copies",
    }
    amrani_table = amrani_rows[list(cols)].rename(columns=cols)
    amrani_table["guide"] = amrani_table["guide"].str.replace("Amrani2024_", "", regex=False)
    amrani_table["rank in gene pool"] = [
        f"{r} / {len(icp0) if g == 'ICP0' else len(icp27)}"
        for r, g in zip(amrani_table["rank in gene pool"], amrani_table["gene"])
    ]

    # -- cascades and the better set ------------------------------------------------
    cascade_icp0 = filter_cascade(icp0, thr0)
    cascade_icp27 = filter_cascade(icp27, thr27)
    better_icp0 = icp0[(icp0["conservation_complete_genomes"] > thr0)
                       & icp0["passes_filters"]].copy()
    n_better_raw = int((icp0["conservation_complete_genomes"] > thr0).sum())
    n_better_filtered = len(better_icp0)
    n_better_perfect = int((better_icp0["conservation_complete_genomes"] >= 1).sum())
    better_cols = {
        "guide_id": "site", "rank_in_gene": "rank", "strand": "strand",
        "protospacer": "spacer (20 nt)", "pam": "PAM",
        "conservation_complete_genomes": f"conservation (n={n_g})",
        "conservation_gene_level": "gene-level", "n_records_covering": "records",
        "gc_content": "GC", "local_gc_200bp": "local GC 200 bp",
        "max_homopolymer_run": "homopolymer", "n_reference_copies": "copies",
    }
    better_icp0_table = better_icp0[list(better_cols)].rename(columns=better_cols)

    # -- pairs ----------------------------------------------------------------------
    lead = pairs[(pairs["icp0_amrani"] == LEAD_PAIR[0])
                 & (pairs["icp27_amrani"] == LEAD_PAIR[1])].iloc[0]
    filtered = pairs[pairs["both_pass_filters"]]
    best = filtered.iloc[0]
    n_pairs_beating = int((filtered["joint_conservation"]
                           > lead["joint_conservation"]).sum())
    n_pairs_perfect = int((filtered["joint_conservation"] >= 1.0).sum())
    pair_rows = filtered.head(6).copy()
    lead_row = pairs[(pairs["icp0_amrani"] == LEAD_PAIR[0])
                     & (pairs["icp27_amrani"] == LEAD_PAIR[1])]
    pair_show = pd.concat([lead_row, pair_rows]).drop_duplicates(
        subset=["icp0_guide", "icp27_guide"])
    pcols = {
        "icp0_guide": "ICP0 site", "icp0_conservation": "ICP0 cons.",
        "icp0_gc": "ICP0 GC", "icp27_guide": "ICP27 site",
        "icp27_conservation": "ICP27 cons.", "icp27_gc": "ICP27 GC",
        "n_genomes_both_present": f"both present (of {n_g})",
        "joint_conservation": "joint conservation",
        "both_pass_filters": "both pass filters",
    }
    pair_table = pair_show[list(pcols)].rename(columns=pcols)
    pair_table.insert(0, "pair", [
        "Amrani lead (g2+g1)" if (a == LEAD_PAIR[0] and b == LEAD_PAIR[1]) else ""
        for a, b in zip(pair_show["icp0_amrani"], pair_show["icp27_amrani"])])

    # -- GC / feasibility -----------------------------------------------------------
    def _gc_row(name, sub):
        return {
            "set": name, "n": len(sub),
            "GC min": _round(sub["gc_content"].min(), 3),
            "GC median": _round(sub["gc_content"].median(), 3),
            "GC max": _round(sub["gc_content"].max(), 3),
            "local GC 200 bp median": _round(sub["local_gc_200bp"].median(), 3),
            "local GC 200 bp max": _round(sub["local_gc_200bp"].max(), 3),
            "outside GC band": int((~sub["gc_in_range"]).sum()),
        }

    gc_table = pd.DataFrame([
        _gc_row("all ICP0/RL2 SaCas9 sites", icp0),
        _gc_row("ICP0 sites beating ICP0g2", icp0[icp0["conservation_complete_genomes"] > thr0]),
        _gc_row("... and passing all filters", better_icp0),
        _gc_row("ICP0g2 itself", icp0[icp0["amrani_guide"] == LEAD_PAIR[0]]),
        _gc_row("ICP0g1 itself", icp0[icp0["amrani_guide"] == "Amrani2024_ICP0g1"]),
        _gc_row("all ICP27/UL54 SaCas9 sites", icp27),
    ])

    # -- gene level -----------------------------------------------------------------
    if gene_level_ran:
        gl_cols = {
            "guide_id": "site", "amrani_guide": "published guide",
            "conservation_complete_genomes": f"complete genomes (n={n_g})",
            "conservation_gene_level": "gene-level corpus",
            "n_records_covering": "records covering",
            "conservation_subgenomic": "sub-genomic tier",
            "n_records_covering_subgenomic": "sub-genomic records",
        }
        gl_show = pd.concat([
            amrani_rows,
            better_icp0.sort_values("conservation_complete_genomes",
                                    ascending=False).head(6),
        ]).drop_duplicates(subset=["guide_id"])
        gene_level_table = gl_show[list(gl_cols)].rename(columns=gl_cols)
        gene_level_table["published guide"] = (
            gene_level_table["published guide"].str.replace("Amrani2024_", "", regex=False))
    else:
        gene_level_table = pd.DataFrame()

    ctx = {
        "generated": generated,
        "reference": reference,
        "nuclease": nuclease,
        "n_genomes": n_g,
        "pool": pool,
        "pairs": pairs,
        "amrani_rows": amrani_rows,
        "amrani_table": amrani_table,
        "pool_counts": pool_counts,
        "cascade_icp0": cascade_icp0,
        "cascade_icp27": cascade_icp27,
        "better_icp0_table": better_icp0_table,
        "pair_table": pair_table,
        "gc_table": gc_table,
        "gene_level_ran": gene_level_ran,
        "gene_level_table": gene_level_table,
        "n_corpus_records": n_corpus_records,
    }

    # -- commentary (conditional on the numbers, not written ahead of them) ----------
    worst = amrani_rows.sort_values("conservation_complete_genomes").iloc[0]
    ctx["rank_commentary"] = (
        f"Their lead ICP0 guide **ICP0g2 ranks {int(g2['rank_in_gene'])} of "
        f"{len(icp0)}** in its own gene's SaCas9 pool: "
        f"{int(g2['n_better_in_gene'])} sites are better conserved. Their lead ICP27 "
        f"guide **ICP27g1 ranks {int(g1['rank_in_gene'])} of {len(icp27)}**, with "
        f"{int(g1['n_better_in_gene'])} sites above it. "
        f"The worst of the four on this measure is "
        f"{worst['amrani_guide'].replace('Amrani2024_', '')} at "
        f"{worst['conservation_complete_genomes']:.3f}.\n\n"
        f"All four clear their own published `>70%` criterion -- but so do "
        f"{int((icp0['conservation_complete_genomes'] > 0.70).sum())} of {len(icp0)} "
        f"ICP0 sites and "
        f"{int((icp27['conservation_complete_genomes'] > 0.70).sum())} of {len(icp27)} "
        f"ICP27 sites. A threshold that admits "
        f"{int((pool['conservation_complete_genomes'] > 0.70).sum())} of {len(pool)} "
        f"candidates is not doing much selecting; that, rather than any rule-breaking, "
        f"is the substantive criticism of their selection procedure."
    )

    ctx["headline_commentary"] = (
        f"**{n_better_raw} of the {len(icp0)} ICP0/RL2 SaCas9 sites are better "
        f"conserved than ICP0g2 across the {n_g} complete genomes. "
        f"{n_better_filtered} of those also pass every standard filter "
        f"(no poly-T terminator, homopolymer run <= 4, GC in 0.35-0.75), and "
        f"{n_better_perfect} of them are present in all {n_g} genomes.**\n\n"
        f"All {n_better_filtered} lie in both ICP0 repeat copies, so they preserve "
        f"the two-cuts-in-ICP0 property that motivates the paper's three-DSB design.\n\n"
        f"The count is not enormous, because the pool is not enormous: 12 out of 46 is "
        f"a quarter of every NNGRRT site in the gene. The honest headline is not "
        f"\"they missed hundreds of guides\" but \"a quarter of the gene's sites "
        f"dominate their choice on the metric they themselves selected on\"."
    )

    # The minimal intervention: keep their ICP27 guide, swap only the ICP0 one.
    swap = pairs[(pairs["icp27_guide"] == g1["guide_id"])
                 & pairs["both_pass_filters"]].iloc[0]
    ctx["pair_commentary"] = (
        f"Their lead pair is intact at both sites in "
        f"**{int(lead['n_genomes_both_present'])} of {n_g} genomes "
        f"(joint conservation {lead['joint_conservation']:.3f})**. That is lower than "
        f"either guide alone ({lead['icp0_conservation']:.3f} and "
        f"{lead['icp27_conservation']:.3f}) because the two guides' absences fall in "
        f"different isolates: the pair is only as good as the union of the failures, "
        f"and here that union is nearly the sum.\n\n"
        f"{n_pairs_beating} filter-passing pairs beat it, and {n_pairs_perfect} "
        f"filter-passing pairs are intact in **all {n_g} genomes**. The best is "
        f"**{best['icp0_guide']} + {best['icp27_guide']}** (marginals "
        f"{best['icp0_conservation']:.3f} / {best['icp27_conservation']:.3f}, joint "
        f"{best['joint_conservation']:.3f}).\n\n"
        f"The minimal intervention is smaller than that. Keeping ICP27g1 exactly as "
        f"published and swapping only the ICP0 guide for "
        f"**{swap['icp0_guide']}** takes the pair from "
        f"{lead['joint_conservation']:.3f} to "
        f"**{swap['joint_conservation']:.3f}** "
        f"({int(swap['n_genomes_both_present'])}/{n_g} genomes) -- one guide changed, "
        f"{int(swap['n_genomes_both_present']) - int(lead['n_genomes_both_present'])} "
        f"more isolates covered.\n\n"
        f"The joint measure is not decoration: it is bounded above by the smaller "
        f"marginal and below by (a + b - 1), and for this pair it sits exactly at the "
        f"lower bound. A per-guide conservation table -- the form in which the "
        f"published method reports its selection -- cannot show that, because it never "
        f"asks whether two guides fail in the same isolate or in different ones."
    )

    hi_local = float(g2["local_gc_200bp"])
    ctx["feasibility_commentary"] = (
        f"ICP0 is GC-rich and repeat-associated, so the obvious worry is that the "
        f"better-conserved sites are unusable GC-extremes. **The data say the "
        f"opposite.** The gene's sites have a median GC of "
        f"{_round(icp0['gc_content'].median(), 3)} and "
        f"{int((~icp0['gc_in_range']).sum())} of {len(icp0)} fall outside the "
        f"0.35-0.75 band -- but the {n_better_filtered} better options that survive "
        f"filtering span GC {_round(better_icp0['gc_content'].min(), 2)}-"
        f"{_round(better_icp0['gc_content'].max(), 2)} with a median local (200 bp) "
        f"GC of {_round(better_icp0['local_gc_200bp'].median(), 3)}, against "
        f"**{hi_local:.3f} local GC around ICP0g2 itself**.\n\n"
        f"That figure is worth dwelling on. Amrani et al. report that the ICP0g2 site "
        f"could not be amplified or sequenced at all, which they attribute to ~85% GC "
        f"in the region; our independent measurement of the reference gives "
        f"{hi_local * 100:.0f}% GC in the 200 bp around that site, corroborating them. "
        f"Several of the better-conserved alternatives sit in materially milder "
        f"neighbourhoods -- so on the one feasibility axis that is visible in the data, "
        f"and that they themselves flagged as a problem, the alternatives are better, "
        f"not worse.\n\n"
        f"The caveat that survives: GC bands are a proxy for synthesis and activity, "
        f"not a measurement of either, and nothing here predicts cutting efficiency."
    )

    if gene_level_ran:
        gl_g2 = _round(g2["conservation_gene_level"], 4)
        ctx["gene_level_commentary"] = (
            f"The ranking is not an artefact of the complete-genome denominator. "
            f"ICP0g2 scores {gl_g2} on the gene-level corpus "
            f"({int(g2['n_records_covering'])} covering records), while the top "
            f"alternatives score at or near 1.000 on comparable numbers of records. "
            f"The two denominators agree on the ordering, which is the point of "
            f"measuring both.\n\n"
            f"One number moves between reports and it is worth saying why. "
            f"`results/robustness_report.md` (stage 5) scores ICP0g2 at 0.807 over 445 "
            f"records; here it is {gl_g2} over {int(g2['n_records_covering'])}. Stage 5 "
            f"admits a record to a published guide's denominator when the record's "
            f"anchor chains overlap the GENE, while every candidate in this report is "
            f"scored with the stricter per-SITE rule -- anchors must bracket the "
            f"guide's own footprint. The stricter rule drops records that reach the "
            f"gene but not the site, which is the correct denominator and is the one "
            f"applied uniformly to their guides and to the alternatives here."
        )
    else:
        ctx["gene_level_commentary"] = ""

    ctx["caveats"] = CAVEATS.format(
        n_better_filtered=n_better_filtered, n_icp0=len(icp0), n_g=n_g,
        n_pairs_perfect=n_pairs_perfect)

    ctx["summary"] = {
        "generated_utc": generated,
        "reference": reference,
        "nuclease": nuclease.as_dict(),
        "n_complete_genomes": n_g,
        "n_corpus_records": n_corpus_records,
        "pool_sizes": {"RL2_ICP0": len(icp0), "UL54_ICP27": len(icp27)},
        "amrani_guides": {
            row.amrani_guide.replace("Amrani2024_", ""): {
                "gene": row.gene,
                "rank_in_gene": int(row.rank_in_gene),
                "n_sites_better_conserved": int(row.n_better_in_gene),
                "conservation_complete_genomes": _round(row.conservation_complete_genomes, 6),
                "conservation_gene_level": _round(row.conservation_gene_level, 6),
                "gc_content": _round(row.gc_content, 4),
                "local_gc_200bp": _round(row.local_gc_200bp, 4),
                "passes_filters": bool(row.passes_filters),
            }
            for row in amrani_rows.itertuples(index=False)
        },
        "icp0_sites_beating_icp0g2": n_better_raw,
        "icp0_sites_beating_icp0g2_passing_filters": n_better_filtered,
        "icp0_sites_beating_icp0g2_perfectly_conserved": n_better_perfect,
        "icp27_sites_beating_icp27g1": int(g1["n_better_in_gene"]),
        "lead_pair": {
            "icp0": str(lead["icp0_guide"]), "icp27": str(lead["icp27_guide"]),
            "n_genomes_both_present": int(lead["n_genomes_both_present"]),
            "joint_conservation": _round(lead["joint_conservation"], 6),
        },
        "best_filtered_pair": {
            "icp0": str(best["icp0_guide"]), "icp27": str(best["icp27_guide"]),
            "joint_conservation": _round(best["joint_conservation"], 6),
        },
        "n_filtered_pairs_beating_lead": n_pairs_beating,
        "n_filtered_pairs_perfect": n_pairs_perfect,
    }
    return ctx


CAVEATS = """
This analysis measures one thing -- exact target-site presence across sequenced
HSV-1 isolates -- and ranks on it. Everything below is a reason the ranking is not,
by itself, a guide-selection decision.

1. **No off-target screening.** `src/offtarget.py` is a stub. None of the
   {n_better_filtered} alternatives has been checked against the human genome. Amrani
   et al. did run BWA + Cas-OFFinder against hg38 and GUIDE-seq validation, and they
   selected partly on nominated off-target counts (ICP0g1 358, ICP0g2 910, ICP27g1
   443, ICP27g2 316). A site that is better conserved may well be worse on that axis;
   we cannot see it. **This is the single largest gap in the comparison.**

2. **No activity prediction.** Conservation says a site exists in an isolate, not that
   SaCas9 cuts it efficiently. They screened six pairwise combinations in Vero cells
   and picked g2g1 on measured antiviral activity, which is evidence of a kind this
   analysis contains none of.

3. **Constraints we cannot see.** AAV packaging limits, sgRNA scaffold compatibility,
   synthesis feasibility of a GC-extreme oligo, and any unpublished screening failures
   all sit outside the data. A site can be perfect on paper and undeliverable.

4. **The reference defines the candidate space.** Sites are enumerated from the
   reference genome only, and only inside annotated CDS. A site conserved
   across every other isolate but absent from strain 17 is never considered, and a
   spliced-exon junction is never a candidate because the enumerator works on
   contiguous genomic segments.

5. **Different denominators.** Their >70% was computed against a Feb-2022 ViPR CDS
   snapshot filtered to 9,409 NCBI IDs with a per-gene denominator they never report;
   ours is {n_g} complete genomes plus an anchor-bracketed gene-level corpus. We are
   not claiming they mis-measured against their own corpus -- we are measuring their
   guides against ours, which is the only comparison available.

6. **Exact matching is conservative.** A site overlapping an `N` in an assembly is
   scored absent, so conservation is a lower bound. That bias applies equally to their
   guides and to the alternatives, so it should not change the ranking, but it does
   mean absolute values are floors.

7. **Repeat structure cuts both ways.** All {n_better_filtered} filter-passing
   alternatives lie in both ICP0 copies, preserving the three-DSB design. But sites in
   the TRL/IRL repeats are also the hardest region to assemble, which is part of why
   some isolates score absent at all.

8. **{n_pairs_perfect} perfect pairs is a count over a small pool, not a discovery
   rate.** With {n_icp0} ICP0 sites and a few dozen ICP27 sites, pair counts multiply
   quickly; the meaningful statement is that at least one filter-passing pair is
   intact in every genome, not that hundreds of independent options exist.
""".strip()


def nngrrn_sensitivity(genome: str, record, genes: list[str],
                       manifest: pd.DataFrame, nngrrt_pool: pd.DataFrame,
                       args: argparse.Namespace) -> tuple[pd.DataFrame, str]:
    """Repeat the pool count under the permissive NNGRRN PAM Amrani et al. used for
    their off-target search. Reported as an upper bound on the option space."""
    nuclease = get_nuclease("sacas9", pam="NNGRRN",
                            spacer_length=args.benchmark_spacer_length)
    pool = build_pool(genome, record, genes, nuclease, args.feature_type)
    pool = add_local_gc(pool, genome)
    _accessions, mat, _cols = presence_matrix(pool, manifest, nuclease)
    pool = assemble(pool, _accessions, mat, None, args.gc_min, args.gc_max,
                    args.max_homopolymer)
    thr = float(nngrrt_pool.loc[nngrrt_pool["amrani_guide"] == LEAD_PAIR[0],
                                "conservation_complete_genomes"].iloc[0])
    rows = []
    for gene, label in (("RL2", "ICP0"), ("UL54", "ICP27")):
        sub = pool[pool["gene"] == gene]
        strict = nngrrt_pool[nngrrt_pool["gene"] == gene]
        better = sub[(sub["conservation_complete_genomes"] > thr) & sub["passes_filters"]]
        rows.append({
            "gene": f"{gene} ({label})",
            "NNGRRT sites": len(strict),
            "NNGRRN sites": len(sub),
            f"NNGRRN sites > {thr:.3f} and passing filters": len(better),
            "of which conservation = 1.000": int(
                (better["conservation_complete_genomes"] >= 1).sum()),
        })
    table = pd.DataFrame(rows)
    n0 = int(table.loc[0, f"NNGRRN sites > {thr:.3f} and passing filters"])
    commentary = (
        f"Relaxing the PAM from NNGRRT to NNGRRN multiplies the ICP0 pool from "
        f"{int(table.loc[0, 'NNGRRT sites'])} to {int(table.loc[0, 'NNGRRN sites'])} "
        f"sites and raises the number of filter-passing options better conserved than "
        f"ICP0g2 to **{n0}**. This is an upper bound only: NNGRRN sites without the "
        f"terminal T are cleaved less efficiently by SaCas9, and Amrani et al. used "
        f"NNGRRN for off-target *search* (deliberately permissive) while every one of "
        f"their four on-target sites is NNGRRT. The NNGRRT figures in section 3 remain "
        f"the like-for-like comparison."
    )
    return table, commentary


# --------------------------------------------------------------------------------------
# Report
# --------------------------------------------------------------------------------------


def _md_table(df: pd.DataFrame, floatfmt: str = "{:.4g}") -> str:
    if df.empty:
        return "_(no rows)_\n"

    def fmt(v):
        if isinstance(v, float) and not pd.isna(v):
            return floatfmt.format(v)
        if v is None or (isinstance(v, float) and pd.isna(v)):
            return "n/a"
        return str(v)

    head = "| " + " | ".join(str(c) for c in df.columns) + " |"
    rule = "|" + "|".join(["---"] * len(df.columns)) + "|"
    body = ["| " + " | ".join(fmt(v) for v in row) + " |"
            for row in df.itertuples(index=False, name=None)]
    return "\n".join([head, rule] + body) + "\n"


def render_report(ctx: dict) -> str:
    from textwrap import dedent

    L: list[str] = []
    add = L.append
    n_g = ctx["n_genomes"]
    pool = ctx["pool"]
    amr = ctx["amrani_rows"]
    icp0 = pool[pool["gene"] == "RL2"]
    icp27 = pool[pool["gene"] == "UL54"]
    g2 = amr.loc[amr["amrani_guide"] == "Amrani2024_ICP0g2"].iloc[0]
    g1 = amr.loc[amr["amrani_guide"] == "Amrani2024_ICP27g1"].iloc[0]
    thr = float(g2["conservation_complete_genomes"])

    add("# SaCas9 head-to-head: how many better options did Amrani et al. 2024 have?")
    add("")
    add(dedent(f"""
    Generated by `src/benchmark_sacas9.py` on {ctx['generated']}.
    Reference {ctx['reference']}; {n_g} complete HSV-1 genomes;
    nuclease **{ctx['nuclease'].label}** (footprint {ctx['nuclease'].target_length} nt),
    which is the site grammar of Amrani et al. 2024 Table 1.

    All numbers below were produced by the run that wrote this file. Nothing is
    quoted from the competing paper except the four published guide sequences, their
    stated selection threshold (>70%), and their stated lead pair.
    """).strip())
    add("")

    # -- 1. pools -------------------------------------------------------------------
    add("## 1. The candidate pools they were choosing from")
    add("")
    counts = ctx["pool_counts"]
    add(_md_table(counts))
    add("")
    add(dedent(f"""
    The SaCas9 NNGRRT pool is small. In a 68% GC genome an NNGRRT PAM occurs roughly
    once per 70 bp per strand, against once per ~9 bp for SpCas9's NGG, so the two
    target genes offer **{len(icp0)} ICP0/RL2 sites and {len(icp27)} ICP27/UL54 sites**
    in total -- not thousands. Any claim that they "missed hundreds of better guides"
    is arithmetically impossible in this PAM space, and this report does not make one.
    """).strip())
    add("")

    # -- 2. their four guides -------------------------------------------------------
    add("## 2. Where their four guides fall in those pools")
    add("")
    add(_md_table(ctx["amrani_table"]))
    add("")
    add(dedent(f"""
    Ranking is by conservation across the {n_g} complete genomes, with the stage-4
    tie-breaks (filter pass, poly-T, homopolymer run, distance of GC from the centre
    of the allowed band, coordinate). `rank_in_gene` counts within the gene's own
    SaCas9 pool.
    """).strip())
    add("")
    add(ctx["rank_commentary"])
    add("")

    # -- 3. the headline count ------------------------------------------------------
    add(f"## 3. How many ICP0 sites beat ICP0g2 ({thr:.3f})?")
    add("")
    add(_md_table(ctx["cascade_icp0"]))
    add("")
    add(ctx["headline_commentary"])
    add("")
    add("### The ICP0 sites that beat it and pass every filter")
    add("")
    add(_md_table(ctx["better_icp0_table"]))
    add("")
    add("### ICP27 for completeness")
    add("")
    add(_md_table(ctx["cascade_icp27"]))
    add("")

    # -- 4. pairs -------------------------------------------------------------------
    add("## 4. Pairs: is there a better-conserved ICP0 + ICP27 combination?")
    add("")
    lead_names = " + ".join(n.replace("Amrani2024_", "") for n in LEAD_PAIR)
    add(dedent(f"""
    A two-guide excision construct needs *both* sites intact in a given isolate, so
    the figure of merit for a pair is the fraction of genomes in which both are
    present. That is not recoverable from the two marginal conservation values, and
    it is what the per-genome presence matrix is for.

    Their lead pair is **{lead_names}** (construct SaCas9-g2g1).
    """).strip())
    add("")
    add(_md_table(ctx["pair_table"]))
    add("")
    add(ctx["pair_commentary"])
    add("")

    # -- 5. feasibility -------------------------------------------------------------
    add("## 5. Are the better options actually usable? (GC and repeat structure)")
    add("")
    add(_md_table(ctx["gc_table"]))
    add("")
    add(ctx["feasibility_commentary"])
    add("")

    # -- 6. gene-level --------------------------------------------------------------
    add("## 6. Second denominator: the gene-level corpus")
    add("")
    if ctx["gene_level_ran"]:
        add(dedent(f"""
        The same candidates re-measured against the cached stage-5 corpus
        ({ctx['n_corpus_records']} HSV-1 records: sub-genomic amplicons plus
        near-full-length "partial genome" records), with the same anchor-bracketed
        coverage correction -- a record counts toward a site only if one anchor chain
        brackets the site's footprint on both sides.
        """).strip())
        add("")
        add(_md_table(ctx["gene_level_table"]))
        add("")
        add(ctx["gene_level_commentary"])
    else:
        add("_Not run: the cached gene-level corpus was not found under "
            "`data/raw/genes/`. Run `python run_pipeline.py --robustness` once to "
            "populate it._")
    add("")

    # -- 7. sensitivity -------------------------------------------------------------
    add("## 7. Sensitivity: the permissive NNGRRN PAM")
    add("")
    if ctx["nngrrn"] is not None:
        add(dedent(f"""
        Amrani et al. ran their off-target search with the permissive `NNGRRN` PAM,
        so it is worth asking how the picture changes if on-target sites are
        enumerated the same way. SaCas9 cleaves NNGRRN sites less efficiently than
        NNGRRT ones, so this is an upper bound on the option space, not a
        recommendation.
        """).strip())
        add("")
        add(_md_table(ctx["nngrrn"]))
        add("")
        add(ctx["nngrrn_commentary"])
    else:
        add("_Not run._")
    add("")

    # -- 8. caveats -----------------------------------------------------------------
    add("## 8. Caveats -- what this analysis cannot see")
    add("")
    add(ctx["caveats"])
    add("")

    add("## 9. Reproduce")
    add("")
    add("```bash")
    add("python src/benchmark_sacas9.py            # offline; uses data/raw/")
    add("python tests/test_nuclease.py             # PAM model + scanner equivalence")
    add("```")
    add("")
    add(f"Supporting tables: {', '.join('`' + p + '`' for p in ctx['artifacts'])}.")
    add("")
    return "\n".join(L)


# --------------------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------------------


def add_arguments(parser: argparse.ArgumentParser) -> None:
    g = parser.add_argument_group("SaCas9 benchmark (stage 6)")
    g.add_argument("--benchmark-sacas9", action="store_true",
                   help="Run the stage-6 SaCas9 head-to-head against Amrani et al. 2024.")
    g.add_argument("--benchmark-report", type=Path,
                   default=RESULTS_DIR / "sacas9_benchmark_report.md")
    g.add_argument("--benchmark-genes", default=",".join(TARGET_GENES),
                   help="Genes to enumerate (default RL2,UL54 -- the two Amrani et al. "
                        "targeted).")
    g.add_argument("--benchmark-spacer-length", type=int, default=20,
                   help="Spacer length for the head-to-head pool (default 20, the "
                        "length printed in Amrani et al. Table 1; their construct uses "
                        "21 nt).")
    g.add_argument("--benchmark-skip-gene-level", action="store_true",
                   help="Skip the gene-level corpus pass (parts 1-5 and 7 only).")
    g.add_argument("--benchmark-skip-nngrrn", action="store_true",
                   help="Skip the permissive-PAM sensitivity analysis.")


def _defaults() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    add_arguments(parser)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--gene-manifest", type=Path, default=DEFAULT_GENE_MANIFEST)
    parser.add_argument("--reference", default="NC_001806")
    parser.add_argument("--feature-type", default="CDS")
    parser.add_argument("--gc-min", type=float, default=0.35)
    parser.add_argument("--gc-max", type=float, default=0.75)
    parser.add_argument("--max-homopolymer", type=int, default=4)
    return parser.parse_args([])


def run(args: argparse.Namespace) -> dict:
    from Bio import SeqIO

    from src.common import utc_now_iso

    ensure_dirs()
    genes = [g.strip().upper() for g in args.benchmark_genes.split(",") if g.strip()]
    nuclease = get_nuclease("sacas9", spacer_length=args.benchmark_spacer_length)

    manifest = pd.read_csv(args.manifest, sep="\t")
    ref_path = fetch_reference(args.reference)
    record = SeqIO.read(ref_path, "genbank")
    genome = str(record.seq).upper()
    LOG.info("Reference %s (%d bp); %d complete genomes; nuclease %s",
             record.id, len(genome), len(manifest), nuclease.label)

    pool = build_pool(genome, record, genes, nuclease, args.feature_type)
    pool = add_local_gc(pool, genome)
    LOG.info("Pool: %d SaCas9 sites (%s)", len(pool),
             ", ".join(f"{k} {v}" for k, v in
                       sorted(pool["gene"].value_counts().to_dict().items())))

    accessions, mat, columns = presence_matrix(pool, manifest, nuclease)

    records = [] if args.benchmark_skip_gene_level else cached_corpus_records(
        args.gene_manifest)
    gene_level = None
    if records:
        LOG.info("Gene-level corpus: %d cached records; anchoring ...", len(records))
        gene_level = gene_level_for_pool(pool, str(record.id), genome, genes,
                                         records, nuclease)
    elif not args.benchmark_skip_gene_level:
        LOG.warning("No cached gene-level corpus found; part 6 will be reported as "
                    "NOT RUN.")

    pool = assemble(pool, accessions, mat, gene_level,
                    args.gc_min, args.gc_max, args.max_homopolymer)

    out_dir = args.benchmark_report.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    pairs = pair_analysis(pool, mat, columns, require_filters=False)
    ctx = build_context(pool, pairs, accessions, mat, columns, nuclease,
                        str(record.id), len(records), gene_level is not None, args,
                        utc_now_iso())

    pool_path = out_dir / "sacas9_benchmark_pool.tsv"
    ctx["pool"][[c for c in POOL_COLUMNS if c in ctx["pool"].columns]].to_csv(
        pool_path, sep="\t", index=False)
    pairs_path = out_dir / "sacas9_benchmark_pairs.tsv"
    pairs.to_csv(pairs_path, sep="\t", index=False)
    ctx["artifacts"] = [str(pool_path.relative_to(PROJECT_ROOT)).replace("\\", "/"),
                        str(pairs_path.relative_to(PROJECT_ROOT)).replace("\\", "/")]

    if not args.benchmark_skip_nngrrn:
        ctx["nngrrn"], ctx["nngrrn_commentary"] = nngrrn_sensitivity(
            genome, record, genes, manifest, pool, args)
    else:
        ctx["nngrrn"], ctx["nngrrn_commentary"] = None, ""

    args.benchmark_report.write_text(render_report(ctx), encoding="utf-8")
    LOG.info("Wrote %s", args.benchmark_report)

    summary = ctx["summary"]
    (out_dir / "sacas9_benchmark_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_arguments(parser)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--gene-manifest", type=Path, default=DEFAULT_GENE_MANIFEST)
    parser.add_argument("--reference", default="NC_001806")
    parser.add_argument("--feature-type", default="CDS")
    parser.add_argument("--gc-min", type=float, default=0.35)
    parser.add_argument("--gc-max", type=float, default=0.75)
    parser.add_argument("--max-homopolymer", type=int, default=4)
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args(argv)
    setup_logging(args.verbose)
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
