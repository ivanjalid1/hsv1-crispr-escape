"""Offline unit tests for `src/reconcile.py`, the stage-6/7/8 reconciliation.

Run with:  python tests/test_reconcile.py      (or: pytest tests/)

`results/recommendation.md` states a selection rule in prose and `src/reconcile.py`
encodes the same rule in `GATES`, `DISCRIMINATOR` and `TIE_BREAKS`. The danger is not
that the code is wrong in isolation -- it is short -- but that it quietly stops
meaning what the document says. These tests pin the properties the document actually
relies on:

* Pareto dominance is antisymmetric, irreflexive, and insensitive to row order, so
  "the front has three members" is a fact about the data and not about a sort;
* a candidate that is weakly better everywhere and strictly better somewhere really
  does dominate, and one that trades an axis really does not;
* every hard gate can individually eliminate a candidate, so none of the six is
  decorative;
* the discriminator is applied only *after* the gates, which is the whole reason
  `RL2_3441+` is excluded rather than merely ranked low.

Everything here runs on synthetic frames. No pipeline output is required, so these
tests pass on a clean checkout.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.reconcile import (
    DISCRIMINATOR,
    GATES,
    PARETO_HIGHER,
    PARETO_LOWER_CORE,
    TIE_BREAKS,
    apply_rule,
    pareto_front,
)

#: A candidate that passes every gate. Individual tests degrade one field at a time.
CLEAN = {
    "guide_id": "X", "published_name": "", "rank_in_gene": 1,
    "cons_183": 1.0, "cons_gene_level": 1.0, "n_gene_level_records": 429,
    "cons_subgenomic": 1.0, "n_subgenomic_records": 36, "n_reference_copies": 2,
    "gc_content": 0.65, "local_gc_200bp": 0.72, "max_homopolymer_run": 2,
    "has_polyT": False, "passes_filters": True, "q_site": 0.009486,
    "joint_conservation_with_icp27g1": 0.9781, "n_genomes_both_present": 179,
    "p_escape_with_icp27g1": 2.5145e-04,
    "nngrrt_le3": 0, "nngrrt_le4": 3, "nngrrt_seed_le1": 1,
    "nngrrn_le3": 2, "nngrrn_le4": 41, "nngrrn_seed_le1": 19,
    "cds_le3": 0, "nngrrt_le3_21nt": 0, "nngrrt_le4_21nt": 0,
    "nngrrn_le3_21nt": 0, "nngrrn_le4_21nt": 6, "cds_le3_21nt": 0,
}


def frame(*overrides: dict) -> pd.DataFrame:
    rows = []
    for i, over in enumerate(overrides):
        row = dict(CLEAN)
        row["guide_id"] = over.pop("guide_id", f"G{i}")
        row.update(over)
        rows.append(row)
    return pd.DataFrame(rows)


def test_dominance_is_irreflexive_and_antisymmetric() -> None:
    df = frame({"guide_id": "A"}, {"guide_id": "B", "nngrrt_le4": 9},
               {"guide_id": "C", "nngrrn_le4": 90})
    front, dominance = pareto_front(df, PARETO_HIGHER, PARETO_LOWER_CORE)
    assert all(a != b for a, b in dominance), "a candidate dominated itself"
    pairs = set(dominance)
    assert not any((b, a) in pairs for a, b in pairs), "dominance is not antisymmetric"
    assert set(front) <= set(df["guide_id"])
    print("ok  Pareto dominance is irreflexive and antisymmetric")


def test_front_is_independent_of_row_order() -> None:
    df = frame({"guide_id": "A"},
               {"guide_id": "B", "nngrrt_le4": 5, "nngrrn_le3": 1},
               {"guide_id": "C", "nngrrt_le4": 10, "cons_gene_level": 0.99})
    a = pareto_front(df, PARETO_HIGHER, PARETO_LOWER_CORE)[0]
    b = pareto_front(df.iloc[::-1].reset_index(drop=True),
                     PARETO_HIGHER, PARETO_LOWER_CORE)[0]
    assert sorted(a) == sorted(b), (a, b)
    print("ok  the Pareto front does not depend on row order")


def test_strict_dominance_removes_the_dominated() -> None:
    # B is weakly worse everywhere and strictly worse on two off-target axes.
    df = frame({"guide_id": "A"},
               {"guide_id": "B", "nngrrt_le4": 19, "nngrrn_le4": 78})
    front, dominance = pareto_front(df, PARETO_HIGHER, PARETO_LOWER_CORE)
    assert ("A", "B") in dominance and ("B", "A") not in dominance
    assert front == ["A"], front
    print("ok  a weakly-worse-everywhere candidate is dominated and drops off the front")


def test_a_traded_axis_blocks_dominance() -> None:
    # This is the real situation: A is better on NNGRRT, B is better on NNGRRN <= 3.
    # Neither may dominate, or the document's "the front has three members" is wrong.
    df = frame({"guide_id": "A", "nngrrt_le4": 3, "nngrrn_le3": 2},
               {"guide_id": "B", "nngrrt_le4": 10, "nngrrn_le3": 1})
    front, dominance = pareto_front(df, PARETO_HIGHER, PARETO_LOWER_CORE)
    assert dominance == [], dominance
    assert sorted(front) == ["A", "B"]
    print("ok  a single traded axis blocks dominance in both directions")


def test_every_gate_can_eliminate_on_its_own() -> None:
    """No gate is decorative: each one, violated alone, removes a candidate."""
    breakers = [
        ("G1", {"cons_183": 0.99}),
        ("G2", {"cons_gene_level": 0.90}),
        ("G3", {"passes_filters": False}),
        ("G4", {"n_reference_copies": 1}),
        ("G5", {"nngrrt_le3": 6}),
        ("G6", {"cds_le3": 1}),
    ]
    for label, override in breakers:
        df = frame({"guide_id": "good"}, {"guide_id": "bad", **override})
        result = apply_rule(df)
        assert result["survivors"] == ["good"], (label, result["survivors"])
        failed = result["gate_verdicts"].set_index("guide_id").at["bad", "gates_failed"]
        assert failed.startswith(label), (label, failed)
    # G6 must also see the 21-nt column, which is easy to forget.
    df = frame({"guide_id": "good"}, {"guide_id": "bad", "cds_le3_21nt": 1})
    assert apply_rule(df)["survivors"] == ["good"]
    assert len(GATES) == len(breakers)
    print(f"ok  each of the {len(GATES)} hard gates eliminates on its own "
          f"(including the 21-nt CDS column)")


def test_gates_are_applied_before_the_discriminator() -> None:
    """The cleanest guide on the discriminator must still be gated out if it fails.

    This is exactly the RL2_3441+ case: it is perfectly conserved and passes every
    filter, so a rule that ranked before gating could rescue it.
    """
    df = frame({"guide_id": "gated_out", DISCRIMINATOR: 0, "nngrrt_le3": 6},
               {"guide_id": "survivor", DISCRIMINATOR: 99})
    result = apply_rule(df)
    assert result["survivors"] == ["survivor"]
    assert result["selected"] == "survivor", result["selected"]
    print("ok  hard gates are applied before the discriminator, not after")


def test_discriminator_then_tie_breaks_order_the_survivors() -> None:
    df = frame({"guide_id": "second", DISCRIMINATOR: 5},
               {"guide_id": "first", DISCRIMINATOR: 3},
               {"guide_id": "third", DISCRIMINATOR: 10})
    result = apply_rule(df)
    assert list(result["ranked"]["guide_id"]) == ["first", "second", "third"]
    assert result["selected"] == "first"

    # With the discriminator tied, the first tie-break decides -- and it must be
    # applied in the direction declared in TIE_BREAKS.
    col, lower_is_better = TIE_BREAKS[0]
    better, worse = (1.0, 2.0) if lower_is_better else (2.0, 1.0)
    df = frame({"guide_id": "worse", col: worse}, {"guide_id": "better", col: better})
    assert apply_rule(df)["selected"] == "better", col
    print("ok  survivors are ordered by the discriminator, then by the tie-breaks")


def test_empty_survivor_set_is_not_a_crash() -> None:
    df = frame({"guide_id": "a", "cons_183": 0.5}, {"guide_id": "b", "cons_183": 0.5})
    result = apply_rule(df)
    assert result["survivors"] == []
    assert result["selected"] is None
    print("ok  a fully eliminated pool returns no selection instead of raising")


def main() -> int:
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in tests:
        fn()
    print(f"\nAll {len(tests)} tests passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
