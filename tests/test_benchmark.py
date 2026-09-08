"""Offline unit tests for stage 6 (`src/benchmark_sacas9.py`).

Run with:  python tests/test_benchmark.py      (or: pytest tests/)

Two of these encode bugs that actually occurred while the stage was being written,
which is why they are phrased as invariants rather than as spot values:

* `test_pair_joint_survives_reordering` -- the pair analysis originally looked matrix
  columns up by *position* in the pool, but the pool is re-sorted by rank before the
  pairs are formed. The result was joint conservation values HIGHER than either
  marginal, which is impossible.
* `test_pair_bounds` -- joint presence is bounded by min(a, b) above and by
  (a + b - n) below. Checking it on integer counts rather than rounded fractions
  matters: at n=183 the two disagree in the sixth decimal.

The real-data tests skip if data/raw is not populated.
"""

from __future__ import annotations

import random
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.benchmark_sacas9 import (
    BENCHMARK_GUIDES,
    add_local_gc,
    assemble,
    build_pool,
    filter_cascade,
    pair_analysis,
)
from src.common import DEFAULT_MANIFEST, REF_DIR, gc_fraction
from src.nuclease import get_nuclease


def _toy(n_genomes: int = 10, seed: int = 3):
    """A tiny pool of 3 RL2 + 2 UL54 sites with a hand-made presence matrix."""
    rng = random.Random(seed)
    ids = ["RL2_1+", "RL2_2+", "RL2_3-", "UL54_1+", "UL54_2-"]
    genes = ["RL2", "RL2", "RL2", "UL54", "UL54"]
    mat = np.array([[rng.random() > 0.25 for _ in ids] for _ in range(n_genomes)])
    # Force a known structure: RL2_1+ and UL54_1+ fail in DIFFERENT genomes, so the
    # joint must be strictly below both marginals.
    mat[:, 0] = True
    mat[:, 3] = True
    mat[0, 0] = False
    mat[1, 3] = False
    pool = pd.DataFrame({
        "guide_id": ids,
        "gene": genes,
        "gc_content": [0.5, 0.6, 0.7, 0.55, 0.65],
        "max_homopolymer_run": [2, 2, 2, 2, 2],
        "has_polyT": [False] * 5,
        "n_reference_copies": [2, 2, 2, 1, 1],
        "ref_start": [100, 200, 300, 400, 500],
        "ref_end": [125, 225, 325, 425, 525],
        "strand": ["+", "+", "-", "+", "-"],
        "local_gc_200bp": [0.7] * 5,
        "target_site": [f"ACGT{i}" for i in range(5)],
        "amrani_guide": [""] * 5,
        "conservation_gene_level": [np.nan] * 5,
    })
    columns = {gid: i for i, gid in enumerate(ids)}
    accessions = [f"ACC{i}" for i in range(n_genomes)]
    return pool, accessions, mat, columns


def test_assemble_marginals_match_matrix():
    pool, accessions, mat, _cols = _toy()
    out = assemble(pool.copy(), accessions, mat, None, 0.35, 0.75, 4)
    for gid, count in zip(out["guide_id"], out["n_strains_present"]):
        col = list(pool["guide_id"]).index(gid)
        assert int(count) == int(mat[:, col].sum()), gid
    assert (out["conservation_complete_genomes"]
            == (out["n_strains_present"] / len(accessions)).round(6)).all()
    print("ok  assemble marginals survive the rank re-sort")


def test_pair_joint_survives_reordering():
    """Regression: pool row order must not affect the pair result."""
    pool, accessions, mat, cols = _toy()
    assembled = assemble(pool.copy(), accessions, mat, None, 0.35, 0.75, 4)
    pairs = pair_analysis(assembled, mat, cols)

    # Same computation from the unsorted pool must give the same table.
    unsorted = pool.copy()
    unsorted["conservation_complete_genomes"] = mat.mean(axis=0).round(6)
    unsorted["passes_filters"] = True
    unsorted["local_gc_200bp"] = 0.7
    ref = pair_analysis(unsorted, mat, cols)
    key = ["icp0_guide", "icp27_guide"]
    a = pairs.sort_values(key).reset_index(drop=True)
    b = ref.sort_values(key).reset_index(drop=True)
    assert (a["n_genomes_both_present"] == b["n_genomes_both_present"]).all()

    # And the hand-made structure must come out right.
    row = a[(a["icp0_guide"] == "RL2_1+") & (a["icp27_guide"] == "UL54_1+")].iloc[0]
    n = mat.shape[0]
    assert int(row["icp0_n_present"]) == n - 1
    assert int(row["icp27_n_present"]) == n - 1
    assert int(row["n_genomes_both_present"]) == n - 2, row.to_dict()
    print("ok  pair joint conservation is index-safe and strictly sub-marginal")


def test_pair_bounds():
    """min(a, b) >= both >= a + b - n, on integer counts, for every pair."""
    for seed in range(6):
        pool, accessions, mat, cols = _toy(n_genomes=12, seed=seed)
        pool["conservation_complete_genomes"] = mat.mean(axis=0).round(6)
        pool["passes_filters"] = True
        pairs = pair_analysis(pool, mat, cols)   # raises if the invariant breaks
        n = mat.shape[0]
        for row in pairs.itertuples(index=False):
            hi = min(row.icp0_n_present, row.icp27_n_present)
            lo = max(0, row.icp0_n_present + row.icp27_n_present - n)
            assert lo <= row.n_genomes_both_present <= hi, row
            assert row.n_genomes_at_least_one + row.n_genomes_neither == n
    print("ok  pair bounds hold on 6 randomised matrices")


def test_filter_cascade_is_monotone():
    pool, accessions, mat, _cols = _toy()
    out = assemble(pool.copy(), accessions, mat, None, 0.35, 0.75, 4)
    for gene in ("RL2", "UL54"):
        casc = filter_cascade(out[out["gene"] == gene], 0.0)
        counts = casc["n_sites"].tolist()
        assert counts == sorted(counts, reverse=True), (gene, counts)
        has_repeat_row = casc["constraint"].str.contains("repeat copies").any()
        # The repeat-copy row is meaningful only for the duplicated gene.
        assert has_repeat_row == (gene == "RL2"), gene
    print("ok  filter cascade is monotone and gene-appropriate")


def test_local_gc_windows():
    genome = "A" * 500 + "G" * 100 + "A" * 500
    pool = pd.DataFrame({"ref_start": [551], "ref_end": [576]})  # inside the GC block
    out = add_local_gc(pool.copy(), genome, flanks=(50, 200))
    assert out["local_gc_50bp"].iloc[0] == 1.0
    # The 200 bp window reaches well outside the 100 bp GC block.
    assert 0.4 < out["local_gc_200bp"].iloc[0] < 0.6
    assert gc_fraction("GGCC") == 1.0
    print("ok  local GC windows")


def test_published_guides_are_in_the_real_pool():
    """Every Amrani et al. Table 1 guide must be recoverable as a SaCas9 20 nt +
    NNGRRT site inside the annotated RL2/UL54 CDS of the reference.

    If this fails, every rank in the stage-6 report is meaningless -- their guides
    would be outside the enumerated space rather than low in it.
    """
    gb = REF_DIR / "NC_001806.gb"
    if not gb.is_file():
        print("SKIP  published-guide recovery (data/raw/reference absent)")
        return
    from Bio import SeqIO

    record = SeqIO.read(gb, "genbank")
    genome = str(record.seq).upper()
    nuclease = get_nuclease("sacas9", spacer_length=20)
    pool = build_pool(genome, record, ["RL2", "UL54"], nuclease)
    sites = set(pool["target_site"])
    for name, (spacer, pam) in BENCHMARK_GUIDES.items():
        assert spacer + pam in sites, f"{name} is not in the enumerated SaCas9 pool"
    assert not DEFAULT_MANIFEST.is_file() or len(pool) > 0
    print(f"ok  all {len(BENCHMARK_GUIDES)} published guides recovered from the "
          f"real reference pool ({len(pool)} sites)")


def main() -> int:
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in tests:
        fn()
    print(f"\nAll {len(tests)} tests passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
