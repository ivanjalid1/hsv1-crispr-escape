"""Tests for `src/figures.py`, the five publication figures.

Run with:  python tests/test_figures.py      (or: pytest tests/)

These tests check the DATA, not the pixels. A figure is a claim about numbers, and
the failure mode that matters is a figure that draws something the results files do
not say -- a stale constant, a column renamed upstream, a bound recomputed with the
wrong denominator. Pixel comparison would catch none of that and would fail on every
cosmetic tweak.

So each `figN_data()` is checked against the file that is authoritative for its
numbers:

* fig 1  point estimate  -> results/escape_k_curve.tsv, value for value;
         certifiable bound -> the table stage 7 printed in section 7.1 of
         results/escape_model_report.md, and the crossings against its section 7.2;
* fig 2  marginals -> results/sacas9_benchmark_pool.tsv;
         joint     -> results/sacas9_benchmark_pairs.tsv;
* fig 3  counts and medians -> results/summary.json, results/sacas9/summary.json;
* fig 4  every plotted column -> results/recommendation_table.tsv;
         Pareto front -> recomputed by src.reconcile.apply_rule;
* fig 5  the variant partitions the corpus exactly (carriers + reference + neither
         = 183) and agrees with the pool's published presence count; the metadata
         missingness matches data/manifest.tsv.

Plus one end-to-end test that all five figures really are written, as both a
non-empty 300 dpi PNG of plausible pixel dimensions and a real PDF, and one that
re-rendering is byte-identical.

These need the pipeline's outputs and the genome FASTAs, so unlike
tests/test_reconcile.py they do NOT pass on a clean checkout; they skip with a clear
message instead.
"""

from __future__ import annotations

import json
import re
import struct
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.common import DEFAULT_MANIFEST, PROJECT_ROOT, RESULTS_DIR  # noqa: E402
from src import figures  # noqa: E402
from src import reconcile  # noqa: E402


class Skip(Exception):
    """Raised when a pipeline output this test needs has not been produced."""


REQUIRED = [
    RESULTS_DIR / "escape_k_curve.tsv",
    RESULTS_DIR / "escape_site_parameters.tsv",
    RESULTS_DIR / "escape_model_report.md",
    RESULTS_DIR / "escape_summary.json",
    RESULTS_DIR / "sacas9-20nt" / "conservation.tsv",
    RESULTS_DIR / "sacas9" / "conservation.tsv",
    RESULTS_DIR / "sacas9" / "summary.json",
    RESULTS_DIR / "summary.json",
    RESULTS_DIR / "conservation.tsv",
    RESULTS_DIR / "sacas9_benchmark_pool.tsv",
    RESULTS_DIR / "sacas9_benchmark_pairs.tsv",
    RESULTS_DIR / "recommendation_table.tsv",
    DEFAULT_MANIFEST,
]


def require_outputs() -> None:
    missing = [str(p.relative_to(PROJECT_ROOT)) for p in REQUIRED if not p.is_file()]
    if missing:
        raise Skip("pipeline outputs not present: " + ", ".join(missing))


_CACHE: dict[str, dict] = {}


def data(n: int) -> dict:
    """figN_data(), computed once per process (fig 2 and 5 scan 183 genomes)."""
    key = f"fig{n}"
    if key not in _CACHE:
        require_outputs()
        _CACHE[key] = figures.FIGURES[n][0]()
    return _CACHE[key]


# --------------------------------------------------------------------------------------
# Helpers that read what the reports actually printed
# --------------------------------------------------------------------------------------


def _report_section(start: str, end: str) -> str:
    """The slice of escape_model_report.md between two headings.

    Necessary because several sections tabulate rows keyed `best k=N`; matching on
    the row label alone silently picks up whichever table comes last.
    """
    text = (RESULTS_DIR / "escape_model_report.md").read_text(encoding="utf-8")
    i = text.find(start)
    if i < 0:
        return ""
    j = text.find(end, i + len(start))
    return text[i:j if j > 0 else len(text)]


def resolution_floor_table() -> dict[str, tuple[float, float]]:
    """Section 7.1: set name -> (P(escape) measured, P(escape) 95% upper bound)."""
    out: dict[str, tuple[float, float]] = {}
    for line in _report_section("### 7.1", "### 7.2").splitlines():
        m = re.match(r"^\|\s*(Amrani lead pair|best k=\d+)\s*\|\s*\d+\s*\|"
                     r"\s*([0-9.e+-]+)\s*\|\s*([0-9.e+-]+)\s*\|", line)
        if m:
            out[m.group(1)] = (float(m.group(2)), float(m.group(3)))
    return out


def certifiable_min_k_table() -> dict[float, tuple[int, int]]:
    """Section 7.2, second table: threshold -> (min k point estimate, min k certifiable)."""
    out: dict[float, tuple[int, int]] = {}
    for line in _report_section("### 7.2", "### 7.3").splitlines():
        m = re.match(r"^\|\s*([0-9.]+e[+-]\d+)\s*\|\s*(\d+)\s*\|\s*(\d+)\s*\|"
                     r"\s*(-?\d+)\s*\|\s*$", line)
        if m:
            out[float(m.group(1))] = (int(m.group(2)), int(m.group(3)))
    return out


def close(a: float, b: float, rel: float = 5e-4) -> bool:
    """Agreement to the 4 significant figures the reports print."""
    return abs(a - b) <= rel * max(abs(a), abs(b), 1e-300)


def png_size(path: Path) -> tuple[int, int]:
    """(width, height) in pixels, read from the IHDR chunk -- no image library."""
    raw = path.read_bytes()
    if raw[:8] != b"\x89PNG\r\n\x1a\n":
        raise AssertionError(f"{path} is not a PNG")
    w, h = struct.unpack(">II", raw[16:24])
    return int(w), int(h)


# --------------------------------------------------------------------------------------
# Figure 1
# --------------------------------------------------------------------------------------


def test_fig1_point_estimate_is_the_k_curve_file_verbatim() -> None:
    d = data(1)
    curve = pd.read_csv(RESULTS_DIR / "escape_k_curve.tsv", sep="\t")
    assert list(d["k"]) == [int(k) for k in curve["k"]]
    assert np.array_equal(d["p_point"], curve["p_escape_best"].to_numpy())
    assert d["p_point"][0] > d["p_point"][-1], "the curve must fall with k"
    print("ok  fig1 point estimate is escape_k_curve.tsv, value for value")


def test_fig1_certifiable_bound_reproduces_stage7s_own_report() -> None:
    """The bound is regenerated here; stage 7 printed it. They must agree."""
    d = data(1)
    published = resolution_floor_table()
    if not published:
        raise Skip("section 7.1 of escape_model_report.md could not be parsed")
    checked = 0
    for k, point, bound in zip(d["k"], d["p_point"], d["p_bound"]):
        row = published.get(f"best k={k}")
        if row is None:
            continue
        assert close(point, row[0]), (k, point, row[0])
        assert close(bound, row[1]), (
            f"k={k}: regenerated bound {bound:.4e} != report's {row[1]:.4e}")
        checked += 1
    assert checked >= 6, f"only {checked} k values cross-checked against the report"
    print(f"ok  fig1 certifiable bound matches stage 7's section 7.1 for {checked} k")


def test_fig1_bound_is_above_the_point_estimate_and_uses_the_rule_of_three() -> None:
    d = data(1)
    assert np.all(d["p_bound"] > d["p_point"]), "a 95% upper bound cannot be lower"
    # Every site in the best k-sets is present in all 183 genomes, so the per-site
    # bound is the rule of three -- the 0.0162 the figure annotates.
    assert all(a == 0 for row in d["absences_per_k"] for a in row), \
        "a best-set site is absent somewhere; the annotation would be wrong"
    assert 0.0160 < d["rule_of_three"] < 0.0165, d["rule_of_three"]
    assert d["n_genomes"] == 183
    print("ok  fig1 bound exceeds the point estimate; rule of three is "
          f"{d['rule_of_three']:.4f} on 0/183")


def test_fig1_crossings_match_the_reports_minimum_k_tables() -> None:
    d = data(1)
    summary = json.loads((RESULTS_DIR / "escape_summary.json").read_text())
    for t, expected in (("0.001", 1e-3), ("1e-06", 1e-6)):
        assert d["crossings"][expected]["point"] == summary["minimum_k"][t], t
    published = certifiable_min_k_table()
    if not published:
        raise Skip("section 7.2 of escape_model_report.md could not be parsed")
    for threshold, (point_k, certifiable_k) in published.items():
        if threshold not in d["crossings"]:
            continue
        assert d["crossings"][threshold]["point"] == point_k, threshold
        assert d["crossings"][threshold]["bound"] == certifiable_k, threshold
    # The figure's headline: the bound needs one more guide at 1e-6.
    assert d["crossings"][1e-6]["point"] == 3
    assert d["crossings"][1e-6]["bound"] == 4
    print("ok  fig1 crossings agree with escape_summary.json and section 7.2 "
          "(1e-06: k=3 point, k=4 certifiable)")


# --------------------------------------------------------------------------------------
# Figure 2
# --------------------------------------------------------------------------------------


def test_fig2_marginals_and_joint_match_the_published_tables() -> None:
    d = data(2)
    for gid in (figures.ICP0G2, figures.ICP27G1):
        assert close(d["marginals"][gid], d["published_marginals"][gid], 1e-5), gid
    assert close(d["joint"], d["published_joint"], 1e-5)
    assert close(d["product_of_marginals"],
                 d["marginals"][figures.ICP0G2] * d["marginals"][figures.ICP27G1])
    print("ok  fig2 marginals match sacas9_benchmark_pool.tsv and the joint matches "
          "sacas9_benchmark_pairs.tsv")


def test_fig2_matrix_is_the_whole_corpus_and_the_counts_add_up() -> None:
    d = data(2)
    assert d["n_genomes"] == 183
    assert len(d["accessions"]) == 183 == len(set(d["accessions"]))
    assert d["icp0_present"].shape == (183,)
    total = (d["n_both"] + d["n_icp0_only_absent"] + d["n_icp27_only_absent"]
             + d["n_neither"])
    assert total == 183, total
    assert d["n_both"] == int(round(d["joint"] * 183))
    assert sorted(np.asarray(d["order"])) == list(range(183))
    print(f"ok  fig2 presence matrix covers all 183 genomes "
          f"({d['n_both']} both present)")


def test_fig2_the_two_guides_never_fail_in_the_same_isolate() -> None:
    """The visual claim of panel A, asserted as arithmetic."""
    d = data(2)
    assert d["n_neither"] == 0, "some isolate lacks both sites; panel A would mislead"
    assert d["joint"] < d["product_of_marginals"], \
        "joint >= product would make the anti-correlation annotation false"
    print("ok  fig2 failures are disjoint (0 isolates lack both) and joint < product")


def test_fig2_absent_in_column_is_truncated_so_the_matrix_is_necessary() -> None:
    """Why fig 2 rescans the genomes instead of reading `absent_in`."""
    pool = pd.read_csv(RESULTS_DIR / "sacas9_benchmark_pool.tsv", sep="\t")
    row = pool[pool["guide_id"] == figures.ICP0G2].iloc[0]
    listed = len(str(row["absent_in"]).split(";"))
    true_absent = int(row["n_strains_total"]) - int(row["n_strains_present"])
    assert listed < true_absent, (
        "absent_in is no longer truncated; the docstring in fig2_data must be updated")
    assert data(2)["n_icp0_only_absent"] == true_absent
    print(f"ok  fig2 rebuilds presence because absent_in lists {listed} of "
          f"{true_absent} ICP0g2 absences")


# --------------------------------------------------------------------------------------
# Figure 3
# --------------------------------------------------------------------------------------


def test_fig3_counts_and_medians_match_the_two_summary_files() -> None:
    d = data(3)
    for label, s in d["spaces"].items():
        assert s["n"] == s["published_n"], label
        assert s["n_perfect"] == s["published_n_perfect"], label
        assert close(s["median"], s["published_median"], 1e-5), label
        assert len(s["conservation"]) == s["n"]
        assert s["conservation"].min() >= 0.0 and s["conservation"].max() <= 1.0
    print("ok  fig3 candidate counts, perfect counts and medians match "
          "summary.json and sacas9/summary.json")


def test_fig3_shows_a_tenfold_contraction_and_a_small_fraction_change() -> None:
    d = data(3)
    ngg, nngrrt = (d["spaces"][k] for k in d["labels"])
    assert ngg["n"] == 4777 and nngrrt["n"] == 448, (ngg["n"], nngrrt["n"])
    assert 9.0 < d["contraction_factor"] < 12.0, d["contraction_factor"]
    assert abs(d["pct_perfect_delta"]) < 5.0, (
        "the conserved fraction moved a lot; panel B's title would be wrong")
    print(f"ok  fig3 {d['contraction_factor']:.1f}x fewer sites, conserved fraction "
          f"moves {d['pct_perfect_delta']:+.1f} points")


# --------------------------------------------------------------------------------------
# Figure 4
# --------------------------------------------------------------------------------------


def test_fig4_plots_only_columns_that_exist_in_the_recommendation_table() -> None:
    d = data(4)
    table = pd.read_csv(RESULTS_DIR / "recommendation_table.tsv", sep="\t")
    for col, _, _ in d["axes"]:
        assert col in table.columns, col
        assert d["table"][col].notna().all(), col
        pd.testing.assert_series_equal(
            d["table"][col].reset_index(drop=True),
            table[col].reset_index(drop=True), check_names=False)
    assert d["n_candidates"] == len(table) == 14
    print(f"ok  fig4 plots {len(d['axes'])} columns straight from "
          f"recommendation_table.tsv ({d['n_candidates']} candidates)")


def test_fig4_marks_the_guides_the_argument_needs() -> None:
    d = data(4)
    ids = set(d["table"]["guide_id"])
    for gid in ("RL2_4496+", "RL2_5335+", "RL2_3441+"):
        assert gid in ids, gid
        assert gid in d["marked"], f"{gid} must be drawn as a named, visible line"
    assert figures.ICP0G2 == reconcile.ICP0G2
    assert figures.ICP27G1 == reconcile.ICP27G1
    # The withdrawn candidate is the whole point: it must be in the pool and it must
    # NOT survive the gates, or figure 4 is telling a different story from section 7.1
    # of results/recommendation.md.
    assert "RL2_3441+" not in d["survivors"]
    assert d["selected"] == "RL2_5335+", d["selected"]
    print("ok  fig4 marks ICP0g2, RL2_5335+ (selected) and RL2_3441+ (eliminated)")


def test_fig4_pareto_front_is_recomputed_not_asserted() -> None:
    d = data(4)
    table = pd.read_csv(RESULTS_DIR / "recommendation_table.tsv", sep="\t")
    rule = reconcile.apply_rule(table)
    assert sorted(d["pareto_front"]) == sorted(rule["pareto_front"])
    assert sorted(d["pareto_front"]) == ["RL2_3364+", "RL2_5080+", "RL2_5335+"], \
        d["pareto_front"]
    # A front member must not be dominated by anything else on the front.
    front, dominance = reconcile.pareto_front(
        table[table["guide_id"].isin(d["survivors"])],
        reconcile.PARETO_HIGHER, reconcile.PARETO_LOWER_CORE)
    assert not [(a, b) for a, b in dominance if b in front]
    print("ok  fig4 Pareto front is src.reconcile's, and no member is dominated")


def test_fig4_withdrawn_candidate_is_worse_where_the_figure_says_it_is() -> None:
    d = data(4)
    t = d["table"].set_index("guide_id")
    assert t.at["RL2_3441+", "cons_183"] == 1.0, \
        "3441+ was withdrawn on off-target, not conservation"
    assert t.at["RL2_3441+", "nngrrt_le3"] > t.at["RL2_4496+", "nngrrt_le3"]
    assert t.at["RL2_3441+", "nngrrt_le4"] > t.at["RL2_5335+", "nngrrt_le4"]
    print("ok  fig4 RL2_3441+ is perfectly conserved yet worse than ICP0g2 on "
          "canonical-PAM off-targets")


# --------------------------------------------------------------------------------------
# Figure 5
# --------------------------------------------------------------------------------------


def test_fig5_variant_is_one_substitution_inside_the_spacer() -> None:
    d = data(5)
    assert len(d["variant"]) == len(d["target"])
    diffs = [i for i, (a, b) in enumerate(zip(d["target"], d["variant"])) if a != b]
    assert diffs == [d["variant_index0"]]
    assert 1 <= d["spacer_position"] <= len(d["spacer"]), \
        "the variant fell in the PAM, not the spacer"
    assert d["ref_base"] == d["target"][d["variant_index0"]]
    assert d["alt_base"] == d["variant"][d["variant_index0"]]
    assert (d["ref_base"], d["alt_base"], d["spacer_position"]) == ("G", "A", 9), \
        "the recurrent variant is no longer G->A at spacer position 9"
    print(f"ok  fig5 variant is {d['ref_base']}->{d['alt_base']} at spacer "
          f"position {d['spacer_position']}")


def test_fig5_variant_partitions_the_corpus_and_agrees_with_the_pool() -> None:
    d = data(5)
    assert d["n_genomes"] == 183
    assert d["n_present"] == d["published_n_present"], (
        "the rescan disagrees with sacas9_benchmark_pool.tsv's n_strains_present")
    other = d["n_genomes"] - d["n_present"] - d["n_carriers"]
    assert other == len(d["other_absent"]) >= 0
    assert d["n_absent"] == d["n_carriers"] + other
    assert set(d["carriers"]).issubset(set(d["absent"])), \
        "a carrier also contains the reference site; the panel-A counts would double-count"
    assert d["n_carriers"] == 25, d["n_carriers"]
    print(f"ok  fig5 {d['n_present']} reference + {d['n_carriers']} variant + "
          f"{other} neither = {d['n_genomes']}")


def test_fig5_dominant_variant_really_is_dominant() -> None:
    d = data(5)
    counts = list(d["all_variant_counts"].values())
    assert counts[0] == d["n_carriers"]
    if len(counts) > 1:
        assert counts[0] > counts[1] * 3, (
            "no single variant dominates; 'the recurrent SNP' would be a mis-title")
    print(f"ok  fig5 dominant variant carried by {counts[0]} isolates; "
          f"next most common by {counts[1] if len(counts) > 1 else 0}")


def test_fig5_metadata_missingness_matches_the_manifest() -> None:
    d = data(5)
    manifest = pd.read_csv(DEFAULT_MANIFEST, sep="\t", dtype=str).set_index("accession")
    dated = sum(1 for a in d["carriers"]
                if pd.notna(manifest.at[a, "collection_date"]))
    placed = sum(1 for a in d["carriers"] if pd.notna(manifest.at[a, "country"]))
    assert d["n_with_year"] == dated, (d["n_with_year"], dated)
    assert d["n_with_country"] == placed, (d["n_with_country"], placed)
    assert d["n_with_year"] + d["n_missing_year"] == d["n_carriers"]
    assert d["n_with_country"] + d["n_missing_country"] == d["n_carriers"]
    assert d["n_missing_year"] > 0 and d["n_missing_country"] > 0, (
        "no missing metadata; the honesty panels would be empty")
    for acc, year in d["years"].items():
        assert year is None or 1900 <= year <= 2100, (acc, year)
    print(f"ok  fig5 {d['n_with_year']}/{d['n_carriers']} carriers dated, "
          f"{d['n_with_country']}/{d['n_carriers']} placed -- missingness kept, not dropped")


# --------------------------------------------------------------------------------------
# The files themselves
# --------------------------------------------------------------------------------------


def test_all_five_figures_are_written_as_png_and_pdf() -> None:
    require_outputs()
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp)
        written = figures.build(outdir=out)
        assert sorted(written) == [1, 2, 3, 4, 5]
        for n, paths in written.items():
            png, pdf = paths
            assert png.name == f"fig{n}.png" and pdf.name == f"fig{n}.pdf"
            assert png.stat().st_size > 10_000, f"{png} is suspiciously small"
            assert pdf.stat().st_size > 5_000, f"{pdf} is suspiciously small"
            assert pdf.read_bytes()[:5] == b"%PDF-", f"{pdf} is not a PDF"
            w, h = png_size(png)
            # 85 mm at 300 dpi, +-2 px of rounding.
            assert abs(w - round(figures.COL_IN * 300)) <= 2, (n, w)
            assert 700 <= h <= 1700, (n, h)
        print("ok  five PNGs (300 dpi, 85 mm wide) and five PDFs written")


def test_rendering_is_deterministic() -> None:
    """Two renders of the same figure must be byte-identical, in both formats."""
    require_outputs()
    cheap = (1, 3)  # 2 and 5 rescan 183 genomes; the writer path is shared
    with tempfile.TemporaryDirectory() as a, tempfile.TemporaryDirectory() as b:
        first = figures.build(cheap, Path(a))
        second = figures.build(cheap, Path(b))
        for n in cheap:
            for p, q in zip(first[n], second[n]):
                assert p.read_bytes() == q.read_bytes(), f"{p.name} is not reproducible"
    print("ok  re-rendering produces byte-identical PNG and PDF")


def main() -> int:
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    skipped = 0
    for fn in tests:
        try:
            fn()
        except Skip as exc:
            skipped += 1
            print(f"SKIP {fn.__name__}: {exc}")
    print(f"\n{len(tests) - skipped} passed, {skipped} skipped.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
