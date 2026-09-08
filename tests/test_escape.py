"""Offline unit tests for stage 7, the multiplex escape-probability model.

Run with:  python tests/test_escape.py      (or: pytest tests/)

The load-bearing tests are the LIMITING CASES. A model of this kind is only
checkable if it collapses onto quantities that are already known:

  * `test_k1_reduces_to_single_site_conservation` -- with no repair escape, a
    one-guide set must give exactly 1 - conservation_fraction. That is the number
    stage 3 already publishes, so the escape model cannot quietly disagree with the
    rest of the repository.
  * `test_perfectly_conserved_independent_sites_multiply` -- sites present in every
    genome must give exactly the product of their per-site probabilities, which is
    the textbook independent-multiplex answer the model must reproduce in the regime
    where it applies.
  * `test_correlated_absence_is_not_the_product_of_marginals` -- and it must NOT
    reproduce it anywhere else. This is the whole reason the empirical joint matrix
    is carried around, and the test is built on the real, published stage-6 shape:
    two guides with high marginals whose absences fall in different isolates.
  * `test_monte_carlo_agrees_with_closed_form` -- the seeded simulation of the same
    generative process must land within a few standard errors of the closed form.

The remaining tests pin the indel spectrum arithmetic, the codon-tolerance
extraction, the anchor-based coordinate map (against synthetic SNPs, indels and N
blocks), the repeat models and the determinism of the search.

No network access, and no dependence on any pipeline output: every test builds its
own inputs. The two tests that use the real HSV-1 reference skip loudly if the cache
is absent, so the file runs on a fresh clone.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src import robustness as rob
from src.common import REF_DIR
from src.escape import (
    ESSENTIAL_GENES,
    EscapeModel,
    IndelSpectrum,
    Params,
    binomial_upper_bound,
    build_offset_map,
    builtin_spectrum,
    collect_cds_copies,
    inframe_viability,
    offset_intervals,
    separation_matrix,
    site_repair_escape,
    translate,
)

import pandas as pd

RNG = np.random.default_rng(12345)


def _rand_seq(n: int, rng=RNG) -> str:
    return "".join(rng.choice(list("ACGT"), size=n))


# ======================================================================================
# Limiting cases -- the reason this model is checkable at all
# ======================================================================================


def test_k1_reduces_to_single_site_conservation():
    """k=1 with no repair escape MUST equal 1 - conservation_fraction, exactly.

    Checked on every attainable conservation value for n=20, and then with a non-zero
    q to confirm the full one-site identity (1-c) + c*q.
    """
    n = 20
    for present in range(0, n + 1):
        col = np.zeros((n, 1), dtype=bool)
        col[:present, 0] = True
        conservation = present / n
        m0 = EscapeModel(col, np.array([0.0]), ["g"])
        assert abs(m0.p_escape([0]) - (1 - conservation)) < 1e-15, present
        for q in (0.0, 0.01, 0.5, 1.0):
            m = EscapeModel(col, np.array([q]), ["g"])
            expect = (1 - conservation) + conservation * q
            assert abs(m.p_escape([0]) - expect) < 1e-15, (present, q)
    print("ok  k=1 reduces exactly to single-site conservation")


def test_perfectly_conserved_independent_sites_multiply():
    """Sites present in every genome give exactly PROD q_i, for every subset."""
    n, k = 37, 5
    presence = np.ones((n, k), dtype=bool)
    q = np.array([0.5, 0.1, 0.02, 0.3, 0.007])
    model = EscapeModel(presence, q, [f"g{i}" for i in range(k)])
    for r in range(1, k + 1):
        for cols in ([0], [0, 1], [1, 3, 4], [0, 1, 2, 3], list(range(k)))[:r + 1]:
            expect = float(np.prod(q[cols]))
            assert abs(model.p_escape(cols) - expect) < 1e-15, cols
            # With everything present, the independence shortcut must also agree.
            assert abs(model.p_escape_independent(cols) - expect) < 1e-12, cols
    print("ok  perfectly conserved sites multiply, exactly")


def test_correlated_absence_is_not_the_product_of_marginals():
    """The model must reproduce, not assume away, correlated failure across sites.

    The construction is the published stage-6 shape: two guides whose absences fall
    in DIFFERENT isolates, so joint conservation sits BELOW the product of marginals.
    That has two consequences with opposite signs and the model must get both right:
    coverage is worse than the marginals suggest, while escape (which needs BOTH to
    fail in the SAME genome) is BETTER than an independence calculation would say.
    """
    n = 183
    a = np.ones(n, dtype=bool)
    b = np.ones(n, dtype=bool)
    a[:28] = False              # 155/183 = 0.847, as measured for ICP0g2
    b[28:32] = False            # 179/183 = 0.978, as measured for ICP27g1
    presence = np.column_stack([a, b])
    q = np.array([0.01, 0.01])
    model = EscapeModel(presence, q, ["ICP0g2", "ICP27g1"])

    marg_a, marg_b = a.mean(), b.mean()
    joint = model.joint_conservation([0, 1])
    assert joint < marg_a * marg_b, (joint, marg_a * marg_b)
    assert abs(joint - 151 / 183) < 1e-12          # the published 0.825

    emp = model.p_escape([0, 1])
    ind = model.p_escape_independent([0, 1])
    assert ind > emp, (ind, emp)                    # independence over-states escape
    # No genome has both sites absent here, so the empirical value is exactly the
    # single-escape term from the 32 half-covered isolates plus the double term.
    expect = (28 * q[1] + 4 * q[0] + 151 * q[0] * q[1]) / n
    assert abs(emp - expect) < 1e-15, (emp, expect)
    print("ok  correlated absence is measured, not assumed independent")


def test_completely_correlated_absence_is_worse_than_independence():
    """The opposite sign must also come out right when absences coincide."""
    n = 100
    a = np.ones(n, dtype=bool)
    b = np.ones(n, dtype=bool)
    a[:10] = False
    b[:10] = False               # the SAME ten isolates
    model = EscapeModel(np.column_stack([a, b]), np.array([0.01, 0.01]),
                        ["x", "y"])
    assert model.joint_conservation([0, 1]) > a.mean() * b.mean()
    assert model.p_escape([0, 1]) > model.p_escape_independent([0, 1])
    assert abs(model.p_escape([0, 1]) - (10 + 90 * 1e-4) / 100) < 1e-15
    print("ok  coincident absence makes escape worse than independence predicts")


def test_escape_is_monotone_in_k_and_in_q():
    """Adding a guide can never raise escape; raising any q can never lower it."""
    rng = np.random.default_rng(7)
    presence = rng.random((60, 8)) < 0.9
    q = rng.random(8) * 0.2
    model = EscapeModel(presence, q, [f"s{i}" for i in range(8)])
    prev = 1.0
    for k in range(1, 6):
        cols, val = model.best_set(k, list(range(8)))
        assert val <= prev + 1e-15, (k, val, prev)
        prev = val
    hotter = EscapeModel(presence, np.clip(q * 2, 0, 1), [f"s{i}" for i in range(8)])
    for cols in ([0], [1, 2], [0, 3, 5]):
        assert hotter.p_escape(cols) >= model.p_escape(cols) - 1e-15
    print("ok  escape is monotone in k and in q")


def test_monte_carlo_agrees_with_closed_form():
    """The seeded simulation must land within 4 SE of the closed form."""
    rng = np.random.default_rng(3)
    presence = rng.random((40, 4)) < 0.85
    q = np.array([0.2, 0.3, 0.25, 0.4])
    model = EscapeModel(presence, q, list("abcd"))
    for cols in ([0], [0, 1], [1, 2, 3], [0, 1, 2, 3]):
        closed = model.p_escape(cols)
        est, se = model.monte_carlo(cols, 400_000, seed=99)
        assert abs(est - closed) < 4 * se + 1e-9, (cols, closed, est, se)
    # Determinism: the same seed must give the same answer.
    assert model.monte_carlo([0, 1], 10_000, seed=5) == \
        model.monte_carlo([0, 1], 10_000, seed=5)
    print("ok  Monte Carlo agrees with the closed form and is seed-deterministic")


# ======================================================================================
# Indel spectrum
# ======================================================================================


def test_spectrum_normalisation_and_inframe_fraction():
    for name in ("default", "deletion-heavy", "insertion-heavy", "short-indels",
                 "uniform-1to20", "neuronal_nhej", "ipsc_dividing"):
        spec = builtin_spectrum(name)
        cond = spec.conditional
        assert abs(sum(cond.values()) - 1.0) < 1e-12, name
        assert all(L != 0 for L in cond), name
        assert 0.0 < spec.inframe_fraction < 0.6, (name, spec.inframe_fraction)
    # A spectrum with only multiples of three is entirely in frame, and one with none
    # is entirely out of frame -- the two degenerate ends.
    assert IndelSpectrum("x", {0: 0.5, -3: 0.25, 3: 0.25}).inframe_fraction == 1.0
    assert IndelSpectrum("y", {0: 0.5, -1: 0.25, 1: 0.25}).inframe_fraction == 0.0
    print("ok  indel spectra normalise and report a sane in-frame fraction")


def test_spectrum_rescaling_hits_the_target_exactly():
    spec = builtin_spectrum("default")
    for target in (0.05, 0.1, 0.25, 0.3333, 0.5, 0.9):
        r = spec.rescaled_to_inframe(target)
        assert abs(r.inframe_fraction - target) < 1e-12, target
        assert abs(sum(r.conditional.values()) - 1.0) < 1e-12
        # Relative weights WITHIN each frame class must be untouched.
        base, new = spec.conditional, r.conditional
        assert abs(base[-3] / base[-6] - new[-3] / new[-6]) < 1e-12
    print("ok  in-frame rescaling hits its target and preserves within-class shape")


def test_neuronal_nhej_scenario_is_pinned_to_the_archived_histogram():
    """Pin the measured post-mitotic-neuron scenario, end to end.

    This scenario is the only built-in whose SHAPE comes from published data
    (Ramadoss et al. 2025, Nat Commun 16:9883, Fig. 1d deposited source data, pooled
    over six replicates of SpCas9 RNP editing in human iPSC-derived neurons). If the
    archived histogram in refs/ ever changes, or the loader stops reading the neuron
    column, these numbers move and the test fails -- which is the point: a number
    attributed to a paper must not be silently editable.

    The values below were computed once from
    refs/ramadoss2025_fig1d_indel_histogram.tsv and are documented, with the
    extraction method, in refs/ramadoss2025_notes.md.
    """
    from src.escape import RAMADOSS2025_HISTOGRAM, neuronal_nhej_spectrum

    assert RAMADOSS2025_HISTOGRAM.is_file(), RAMADOSS2025_HISTOGRAM

    neu = builtin_spectrum("neuronal_nhej")
    assert neu.name == "neuronal_nhej"
    assert "MEASURED" in neu.provenance and "ASSUMPTION" in neu.provenance, \
        "the scenario must state BOTH halves of its provenance"
    assert abs(sum(neu.conditional.values()) - 1.0) < 1e-12

    # Pooled measured summary statistics.
    assert abs(neu.p_wt - 0.537735) < 1e-5
    assert abs(neu.inframe_fraction - 0.090341) < 1e-5
    assert abs(neu.mean_deletion_length - 6.598602) < 1e-5

    cond = neu.conditional
    # +1 is the single most probable outcome; +-1 and +-2 dominate and all frameshift.
    assert max(cond, key=cond.get) == 1
    assert abs(cond[1] - 0.294036) < 1e-5
    assert sum(p for L, p in cond.items() if abs(L) <= 2) > 0.75
    assert all(L % 3 != 0 for L in (1, -1, 2, -2))

    # The direction that matters for an essential gene: the measured post-mitotic
    # spectrum is MORE frameshifting than the assumed default, so escape must fall.
    default = builtin_spectrum("default")
    assert neu.inframe_fraction < default.inframe_fraction
    # ... and below the whole 0.10-0.50 band the inframe_fraction group sweeps.
    assert neu.inframe_fraction < 0.10

    # The isogenic dividing-cell arm of the same experiment is the control, and it
    # must sit on the other side of the default.
    ipsc = builtin_spectrum("ipsc_dividing")
    assert abs(ipsc.inframe_fraction - 0.289293) < 1e-5
    assert ipsc.inframe_fraction > default.inframe_fraction > neu.inframe_fraction
    assert ipsc.mean_deletion_length > neu.mean_deletion_length

    # Escape at a site must actually be lower under the neuronal scenario, holding
    # everything else fixed. Frameshifts are lethal in an essential gene, so a more
    # frameshifting spectrum is a less escapable one.
    tol = np.full(200, 0.5)
    p = Params()
    q_neu = site_repair_escape(neu, tol, 100, "UL30", 1, p)["q_site"]
    q_def = site_repair_escape(default, tol, 100, "UL30", 1, p)["q_site"]
    q_ips = site_repair_escape(ipsc, tol, 100, "UL30", 1, p)["q_site"]
    assert q_neu < q_def, (q_neu, q_def)
    # NOTE the iPSC arm deliberately does NOT complete an ordering here. It has the
    # highest in-frame fraction of the three but also much longer deletions (mean
    # 12.7 nt against 6.6 nt), and a long in-frame deletion must survive a product of
    # per-codon tolerances over every codon it removes. The two effects oppose each
    # other, so q is not monotone in the in-frame fraction alone. That is worth
    # pinning rather than glossing: the report's claim is that the spectrum acts
    # ALMOST entirely through its in-frame fraction, not entirely.
    assert q_ips < q_def, (q_ips, q_def)
    assert ipsc.inframe_fraction > default.inframe_fraction

    # It is a SCENARIO, not the default. The default must be untouched.
    from src.escape import load_spectrum
    assert load_spectrum("default").name == "parametric-default"
    assert abs(load_spectrum("default").inframe_fraction - 0.198571) < 1e-5
    print("ok  neuronal_nhej is pinned to the archived Ramadoss et al. 2025 histogram")


def test_neuronal_nhej_is_swept_but_is_not_the_default():
    """The scenario must appear in the sweep and must not displace the baseline."""
    import inspect

    from src import escape as esc

    src = inspect.getsource(esc.sensitivity)
    assert "neuronal_nhej" in src, "the scenario is not in the sensitivity sweep"
    # The sweep's own baseline row uses the spectrum passed in, and the CLI default
    # for that is 'default'. Both must remain so, or every previously published
    # absolute number in this repository silently changes.
    parser_src = inspect.getsource(esc.add_arguments)
    assert '"--indel-spectrum", default="default"' in parser_src
    print("ok  neuronal_nhej is swept as a scenario and is not the default")


def test_user_supplied_spectrum_round_trips():
    import tempfile
    from src.escape import load_spectrum
    with tempfile.TemporaryDirectory() as d:
        path = Path(d) / "spec.tsv"
        path.write_text("length\tprobability\n0\t2\n-3\t1\n1\t1\n", encoding="utf-8")
        spec = load_spectrum(str(path))
        assert spec.provenance == "USER-SUPPLIED"
        assert abs(spec.p_wt - 0.5) < 1e-12
        assert abs(spec.inframe_fraction - 0.5) < 1e-12
    print("ok  user-supplied spectra load, normalise and are labelled as such")


# ======================================================================================
# Per-site repair escape
# ======================================================================================


def test_frameshift_in_an_essential_gene_is_never_escape():
    """A frameshift must contribute exactly zero when frameshift viability is 0."""
    spec = builtin_spectrum("default")
    tol = np.full(100, 1.0)                      # maximally tolerant protein
    d = site_repair_escape(spec, tol, 50, "UL30", 1, Params())
    assert d["q_from_frameshift"] == 0.0
    # With a fully tolerant protein the only surviving route is the in-frame mass.
    assert abs(d["q_site"] - spec.inframe_fraction) < 1e-12
    # Declaring the gene dispensable recovers the whole spectrum.
    d2 = site_repair_escape(spec, tol, 50, "RL2", 1,
                            Params(frameshift_viability=(("RL2", 1.0),)))
    assert abs(d2["q_site"] - 1.0) < 1e-12
    print("ok  frameshifts are lethal in essential genes and only there")


def test_invariant_codons_suppress_escape_relative_to_variable_ones():
    spec = builtin_spectrum("default")
    params = Params(tol_variable=0.5, tol_invariant=0.05)
    variable = np.full(100, params.tol_variable)
    invariant = np.full(100, params.tol_invariant)
    qv = site_repair_escape(spec, variable, 50, "UL30", 1, params)["q_site"]
    qi = site_repair_escape(spec, invariant, 50, "UL30", 1, params)["q_site"]
    assert qi < qv, (qi, qv)
    # A single-codon indel is the dominant in-frame term, so the ratio is close to
    # the tolerance ratio but strictly smaller (longer indels are punished harder).
    assert qi / qv < params.tol_invariant / params.tol_variable + 1e-9
    print("ok  invariant codons suppress escape relative to variable ones")


def test_inframe_viability_window_and_edges():
    tol = np.array([0.1, 0.9, 0.9, 0.1, 0.1])
    # A one-codon indel at a tolerant position is exactly that position's tolerance.
    assert abs(inframe_viability(tol, 1, 1) - 0.9) < 1e-12
    assert abs(inframe_viability(tol, 0, 1) - 0.1) < 1e-12
    # A two-codon indel averages over the windows containing the cut codon.
    got = inframe_viability(tol, 2, 2)
    assert abs(got - (0.9 * 0.9 + 0.9 * 0.1) / 2) < 1e-12, got
    # Widths larger than the protein are clipped rather than crashing.
    assert 0.0 <= inframe_viability(tol, 4, 99) <= 1.0
    assert inframe_viability(np.array([]), 0, 1) == 0.0
    # Monotone: a longer in-frame indel is never more survivable on a uniform profile.
    flat = np.full(50, 0.5)
    vals = [inframe_viability(flat, 25, m) for m in range(1, 6)]
    assert all(vals[i] >= vals[i + 1] for i in range(len(vals) - 1)), vals
    print("ok  in-frame viability window behaves at the edges and is monotone")


def test_repeat_models_have_the_documented_ordering():
    """all-copies <= single <= redundant, and each matches its closed form."""
    spec = builtin_spectrum("default")
    tol = np.full(50, 0.3)
    q1 = site_repair_escape(spec, tol, 25, "RL2", 1, Params())["q_site"]
    out = {m: site_repair_escape(spec, tol, 25, "RL2", 2,
                                 Params(repeat_model=m))["q_site"]
           for m in ("redundant", "single", "all-copies")}
    assert abs(out["single"] - q1) < 1e-15
    assert abs(out["redundant"] - (1 - (1 - q1) ** 2)) < 1e-15
    assert abs(out["all-copies"] - q1 ** 2) < 1e-15
    assert out["all-copies"] < out["single"] < out["redundant"]
    print("ok  the three repeat models match their closed forms and ordering")


def test_essential_gene_set_records_icp0_as_dispensable():
    """A documentation test: RL2/ICP0 must NOT be silently listed as essential."""
    assert "RL2" not in ESSENTIAL_GENES
    for g in ("UL30", "UL19", "UL5", "UL52", "UL29", "UL54"):
        assert g in ESSENTIAL_GENES
    print("ok  ICP0 is recorded as not strictly essential")


# ======================================================================================
# Set search
# ======================================================================================


def test_best_pair_is_exhaustive_and_matches_brute_force():
    rng = np.random.default_rng(11)
    presence = rng.random((50, 12)) < 0.9
    q = rng.random(12) * 0.3
    model = EscapeModel(presence, q, [f"s{i}" for i in range(12)])
    cols, val = model.best_set(2, list(range(12)))
    ex_cols, ex_val = model.exhaustive(2, list(range(12)))
    assert abs(val - ex_val) < 1e-15, (val, ex_val)
    assert cols == ex_cols, (cols, ex_cols)
    print("ok  the k=2 search is exhaustive and agrees with brute force")


def test_greedy_plus_swap_matches_brute_force_on_a_small_pool():
    rng = np.random.default_rng(13)
    presence = rng.random((40, 14)) < 0.88
    q = rng.random(14) * 0.4
    model = EscapeModel(presence, q, [f"s{i}" for i in range(14)])
    for k in (3, 4):
        _, val = model.best_set(k, list(range(14)))
        _, ex_val = model.exhaustive(k, list(range(14)))
        assert val <= ex_val * (1 + 1e-9), (k, val, ex_val)
    print("ok  greedy + local swap reaches the brute-force optimum on a small pool")


def test_search_is_deterministic():
    rng = np.random.default_rng(17)
    presence = rng.random((60, 20)) < 0.9
    q = rng.random(20) * 0.2
    model = EscapeModel(presence, q, [f"s{i}" for i in range(20)])
    for k in (2, 3, 4):
        first = model.best_set(k, list(range(20)))
        for _ in range(3):
            assert model.best_set(k, list(range(20))) == first, k
    print("ok  the set search is deterministic across repeated calls")


def test_constraints_are_respected():
    """Separation and per-gene caps must actually bind."""
    sites = pd.DataFrame({
        "guide_id": [f"g{i}" for i in range(6)],
        "gene": ["A", "A", "A", "B", "B", "C"],
        "ref_start": [100, 120, 5000, 200, 9000, 300],
        "ref_end": [126, 146, 5026, 226, 9026, 326],
    })
    sep = separation_matrix(sites)
    assert not np.isfinite(sep[0, 3])            # different genes -> never constrained
    assert sep[0, 1] == 20.0
    presence = np.ones((10, 6), dtype=bool)
    q = np.array([0.01, 0.011, 0.5, 0.5, 0.5, 0.5])
    model = EscapeModel(presence, q, sites["guide_id"].tolist())
    # Unconstrained the best pair is the two adjacent A sites.
    assert model.best_set(2, list(range(6)))[0] == [0, 1]
    # With a separation floor of 150 they become illegal.
    cols, _ = model.best_set(2, list(range(6)), separation=sep, min_separation=150)
    assert cols != [0, 1] and 0 in cols
    # With at most one guide per gene, no two picks share a gene.
    gene_of = sites["gene"].to_numpy()
    cols, _ = model.best_set(3, list(range(6)), gene_of=gene_of, max_per_gene=1)
    assert len(cols) == 3 and len({gene_of[c] for c in cols}) == 3
    print("ok  separation and per-gene constraints bind in the search")


# ======================================================================================
# Statistical resolution floor
# ======================================================================================


def test_binomial_upper_bound():
    # Rule of three: the 95% upper limit on 0/n is approximately 3/n.
    for n in (100, 183, 500):
        u = binomial_upper_bound(0, n)
        assert abs(u - 3.0 / n) < 0.5 / n, (n, u)
        # It is an exact CP bound: P(X = 0 | p = u) must equal alpha.
        assert abs((1 - u) ** n - 0.05) < 1e-6, (n, u)
    assert binomial_upper_bound(5, 10) > 0.5
    assert binomial_upper_bound(10, 10) == 1.0
    # Monotone in the observed count.
    vals = [binomial_upper_bound(a, 50) for a in range(0, 20)]
    assert all(vals[i] < vals[i + 1] for i in range(len(vals) - 1))
    print("ok  the Clopper-Pearson upper bound is exact, monotone and rule-of-three")


# ======================================================================================
# Coordinate map and codon tolerance extraction
# ======================================================================================


def _map_query(ref: str, query: str, k: int = 25, step: int = 3,
               max_unanchored: int = 60) -> np.ndarray:
    index = rob.build_reference_index(ref, k)
    anchors = rob.find_anchors(query, index, k, step)
    chains = rob.chain_anchors(anchors, k)
    return build_offset_map(len(ref), chains, k, max_unanchored)


def test_offset_map_is_identity_on_an_identical_sequence():
    ref = _rand_seq(3000)
    qmap = _map_query(ref, ref)
    mapped = qmap >= 0
    assert mapped.mean() > 0.95, mapped.mean()
    assert np.all(qmap[mapped] == np.flatnonzero(mapped))
    print("ok  the coordinate map is the identity on an identical sequence")


def test_offset_map_bridges_substitutions():
    """A SNP destroys every k-mer overlapping it; the map must bridge that gap.

    This is the property the whole tolerance profile depends on -- variant positions
    are exactly the positions that lose their anchors, so failing to bridge would
    discard precisely the data the model is built from.
    """
    ref = _rand_seq(3000)
    pos = 1500
    alt = "A" if ref[pos] != "A" else "C"
    query = ref[:pos] + alt + ref[pos + 1:]
    qmap = _map_query(ref, query)
    assert qmap[pos] == pos, qmap[pos]
    assert np.all(qmap[pos - 40:pos + 40] == np.arange(pos - 40, pos + 40))
    print("ok  the coordinate map bridges substitutions rather than dropping them")


def test_offset_map_tracks_an_insertion_and_a_deletion():
    ref = _rand_seq(4000)
    ins = ref[:2000] + "GGGCCC" + ref[2000:]
    qmap = _map_query(ref, ins)
    assert qmap[500] == 500
    assert qmap[3500] == 3506, qmap[3500]         # shifted by the 6 nt insertion
    dele = ref[:2000] + ref[2009:]
    qmap = _map_query(ref, dele)
    assert qmap[500] == 500
    assert qmap[3500] == 3491, qmap[3500]         # shifted by the 9 nt deletion
    print("ok  the coordinate map tracks insertions and deletions by offset")


def test_offset_map_reports_unmapped_regions_as_negative():
    ref = _rand_seq(4000)
    query = ref[:1500] + "N" * 400 + ref[1900:]
    qmap = _map_query(ref, query, max_unanchored=60)
    assert np.all(qmap[1600:1800] < 0), "an N block must not be silently bridged"
    assert qmap[500] == 500 and qmap[3000] == 3000
    print("ok  unresolved regions are reported as unmapped, not bridged")


def test_offset_intervals_split_on_indels():
    ref = _rand_seq(3000)
    query = ref[:1500] + "GGG" + ref[1500:]
    index = rob.build_reference_index(ref, 25)
    chains = rob.chain_anchors(rob.find_anchors(query, index, 25, 3), 25)
    runs = offset_intervals(chains, 25, 60)
    offsets = {r[2] for r in runs if r[3] >= 3}
    assert 0 in offsets and 3 in offsets, offsets
    print("ok  constant-offset runs split at an indel")


def test_translate_matches_biopython():
    from Bio.Seq import Seq
    for _ in range(20):
        cds = _rand_seq(300)
        assert translate(cds) == str(Seq(cds).translate())
    print("ok  the codon table agrees with Biopython")


def test_reference_cds_copies_translate_cleanly():
    """On the real reference every target CDS must translate M...* with one stop.

    This is what pins the reading frame the tolerance profile is computed in. If it
    were wrong, every codon would look variable and the model would be nonsense.
    """
    gb = REF_DIR / "NC_001806.gb"
    if not gb.is_file():
        print("skip  reference GenBank not cached; run stage 1 first")
        return
    from Bio import SeqIO
    record = SeqIO.read(gb, "genbank")
    copies = collect_cds_copies(record, ["RL2", "UL54", "UL30", "UL19"])
    assert {c.gene for c in copies} == {"RL2", "UL54", "UL30", "UL19"}
    assert sum(1 for c in copies if c.gene == "RL2") == 2, "ICP0 is duplicated"
    for c in copies:
        assert c.protein.startswith("M") and c.protein.endswith("*")
        assert c.protein.count("*") == 1
        assert len(c.positions) == 3 * c.n_codons
        assert len(set(c.positions.tolist())) == len(c.positions)
    # The two RL2 copies must encode the same protein.
    rl2 = [c for c in copies if c.gene == "RL2"]
    assert rl2[0].protein == rl2[1].protein
    print("ok  reference CDS copies translate cleanly and ICP0 is duplicated")


def test_tolerance_extraction_finds_a_planted_substitution():
    """End-to-end: plant a codon change in a synthetic isolate and recover it."""
    from src.escape import measure_tolerance, CdsCopy
    import tempfile

    rng = np.random.default_rng(101)
    backbone = _rand_seq(6000, rng)
    # A clean synthetic CDS: ATG ... TAA with no internal stop.
    codons = ["ATG"]
    pool = [c for c in ("GCT", "TGC", "GAT", "TTC", "CAC", "AAA", "CTG", "AAC",
                        "CCG", "CGT", "TCT", "ACC", "GTT", "TAC")]
    codons += [pool[int(i)] for i in rng.integers(0, len(pool), 199)]
    codons.append("TAA")
    cds = "".join(codons)
    ref = backbone[:2000] + cds + backbone[2000:]
    start = 2000

    class FakeRecord:
        pass
    from Bio.Seq import Seq
    rec = FakeRecord()
    rec.seq = Seq(ref)

    copy = CdsCopy("TEST", 0, 1,
                   np.arange(start, start + len(cds), dtype=np.int64),
                   cds, translate(cds))

    with tempfile.TemporaryDirectory() as d:
        rows = []
        from src.common import PROJECT_ROOT
        for i in range(10):
            variant = ref
            if i < 3:
                # Change codon 50 from its reference to a different amino acid.
                j = start + 50 * 3
                variant = ref[:j] + "TGG" + ref[j + 3:]     # Trp
            path = Path(d) / f"iso{i}.fasta"
            path.write_text(">x\n" + variant + "\n", encoding="utf-8")
            try:
                rel = path.relative_to(PROJECT_ROOT)
            except ValueError:
                rel = path
            rows.append({"accession": f"ISO{i}", "fasta_path": str(rel)})
        manifest = pd.DataFrame(rows)
        profiles, qc = measure_tolerance(rec, [copy], manifest, anchor_step=3,
                                         min_resolved_fraction=0.9)

    prof = profiles["TEST"]
    assert prof.n_genomes_used == 10, (prof.n_genomes_used, prof.n_genomes_dropped)
    assert prof.n_variant[50] == 3, prof.n_variant[50]
    assert prof.n_variant.sum() == 3, "no other codon may look variable"
    v = prof.varies(min_variant_strains=1)
    assert v[50] and v.sum() == 1
    assert not prof.varies(min_variant_strains=5)[50]
    tol = prof.tolerance_vector(Params(tol_variable=0.5, tol_invariant=0.05))
    assert abs(tol[50] - 0.5) < 1e-12 and abs(tol[49] - 0.05) < 1e-12
    print("ok  tolerance extraction recovers a planted amino-acid substitution")


def test_tolerance_extraction_rejects_a_scrambled_isolate():
    """An isolate that cannot be mapped must be dropped, not counted as variable."""
    from src.escape import measure_tolerance, CdsCopy
    import tempfile
    from Bio.Seq import Seq
    from src.common import PROJECT_ROOT

    rng = np.random.default_rng(202)
    ref = _rand_seq(4000, rng)
    # A random (not low-complexity) CDS: repeats would be dropped as unanchorable
    # k-mers and the mapping would fail for legitimate isolates too.
    pool = ["GCT", "TGC", "GAT", "TTC", "CAC", "AAA", "CTG", "AAC",
            "CCG", "CGT", "TCT", "ACC", "GTT", "TAC"]
    cds = "ATG" + "".join(pool[int(i)] for i in rng.integers(0, len(pool), 150)) + "TAA"
    ref = ref[:1500] + cds + ref[1500:]
    copy = CdsCopy("TEST", 0, 1, np.arange(1500, 1500 + len(cds), dtype=np.int64),
                   cds, translate(cds))

    class FakeRecord:
        pass
    rec = FakeRecord()
    rec.seq = Seq(ref)

    with tempfile.TemporaryDirectory() as d:
        rows = []
        for i, seq in enumerate([ref, ref, _rand_seq(4000, rng)]):
            path = Path(d) / f"iso{i}.fasta"
            path.write_text(">x\n" + seq + "\n", encoding="utf-8")
            try:
                rel = path.relative_to(PROJECT_ROOT)
            except ValueError:
                rel = path
            rows.append({"accession": f"ISO{i}", "fasta_path": str(rel)})
        profiles, qc = measure_tolerance(rec, [copy], pd.DataFrame(rows),
                                         anchor_step=3)
    prof = profiles["TEST"]
    assert prof.n_genomes_used == 2 and prof.n_genomes_dropped == 1
    assert prof.n_variant.sum() == 0, "a rejected isolate must contribute no variation"
    print("ok  unmappable isolates are dropped rather than scored as variable")


# ======================================================================================
# End-to-end smoke test on the real data, if it is cached
# ======================================================================================


def test_cli_runs_end_to_end():
    from src.common import DEFAULT_MANIFEST
    gb = REF_DIR / "NC_001806.gb"
    if not (gb.is_file() and DEFAULT_MANIFEST.is_file()):
        print("skip  cached genomes/reference absent; run stage 1 first")
        return
    import tempfile
    from src.escape import main
    with tempfile.TemporaryDirectory() as d:
        rc = main(["--max-genomes", "6", "--max-k", "3", "--skip-sensitivity",
                   "--mc-draws", "20000", "--genes", "UL54,RL2",
                   "--outdir", d])
        assert rc == 0
        out = Path(d)
        for name in ("escape_model_report.md", "escape_site_parameters.tsv",
                     "escape_k_curve.tsv", "escape_summary.json"):
            assert (out / name).is_file(), name
            assert (out / name).stat().st_size > 0, name
        import json
        summary = json.loads((out / "escape_summary.json").read_text(encoding="utf-8"))
        assert summary["n_genomes"] == 6
        assert summary["amrani_lead_pair"]["p_escape"] > 0
    print("ok  the CLI runs end to end and writes a complete result set")


def main() -> int:
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in tests:
        fn()
    print(f"\nAll {len(tests)} tests passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
