"""Offline unit tests for stage 5 (src/robustness.py). No network needed.

Run with:  python tests/test_robustness.py     (or: pytest tests/)

The load-bearing claims tested here are the two that decide whether the robustness
numbers mean anything:

  * the MinHash sketch is a real similarity estimator AND is reproducible (it must
    not depend on Python's salted str hash);
  * anchor chaining distinguishes "this record does not contain the site" from
    "this record does not reach the site", and distinguishes a real mismatch from
    an assembly gap. Getting that wrong would silently manufacture either a
    flattering or a catastrophic result in part C.
"""

from __future__ import annotations

import random
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.common import revcomp
from src.robustness import (
    ST_ABSENT,
    ST_AMBIGUOUS,
    ST_NOT_COVERED,
    ST_PRESENT,
    ST_WIDE_BRACKET,
    ambiguous_positions,
    bracket_in_chain,
    build_reference_index,
    canonical_kmer_hashes,
    chain_anchors,
    choose_representatives,
    classify_footprint,
    compare_denominators,
    find_anchors,
    mash_identity,
    minhash_sketch,
    orient_and_anchor,
    parse_footprints,
    score_benchmark,
    single_linkage_clusters,
    sketch_jaccard,
)

import pandas as pd

RNG = random.Random(20260907)


def _seq(n: int, rng: random.Random | None = None) -> str:
    r = rng or RNG
    return "".join(r.choice("ACGT") for _ in range(n))


def _mutate(seq: str, rate: float, rng: random.Random) -> str:
    out = list(seq)
    for i in range(len(out)):
        if rng.random() < rate:
            out[i] = rng.choice([b for b in "ACGT" if b != out[i]])
    return "".join(out)


# --------------------------------------------------------------------------------------
# MinHash
# --------------------------------------------------------------------------------------


def test_kmer_hashes_are_strand_canonical():
    """A sequence and its reverse complement must give the identical hash multiset."""
    s = _seq(500)
    a = np.sort(canonical_kmer_hashes(s, 21))
    b = np.sort(canonical_kmer_hashes(revcomp(s), 21))
    assert a.size == b.size > 0
    assert np.array_equal(a, b)
    print("ok  canonical k-mer hashing is strand-independent")


def test_kmer_hashes_skip_ambiguous_windows():
    s = "ACGT" * 10
    clean = canonical_kmer_hashes(s, 21)
    dirty = canonical_kmer_hashes(s[:20] + "N" + s[21:], 21)
    # 21 windows of the 40-mer contain position 20.
    assert dirty.size == clean.size - min(21, clean.size)
    print("ok  ambiguous windows are excluded from sketches")


def test_sketch_jaccard_tracks_divergence():
    rng = random.Random(11)
    base = _seq(20_000, rng)
    identical = base
    close = _mutate(base, 0.0005, rng)
    far = _mutate(base, 0.05, rng)
    unrelated = _seq(20_000, rng)

    s = 2000
    sk = [minhash_sketch(x, 21, s) for x in (base, identical, close, far, unrelated)]
    j_ident = sketch_jaccard(sk[0], sk[1], s)
    j_close = sketch_jaccard(sk[0], sk[2], s)
    j_far = sketch_jaccard(sk[0], sk[3], s)
    j_none = sketch_jaccard(sk[0], sk[4], s)

    assert j_ident == 1.0
    assert 0.9 < j_close < 1.0, j_close
    assert j_far < j_close
    assert j_none < 0.01, j_none
    assert mash_identity(j_ident, 21) == 1.0
    assert mash_identity(j_close, 21) > mash_identity(j_far, 21)
    print(f"ok  sketch Jaccard ordering (identical=1.0 close={j_close:.3f} "
          f"far={j_far:.3f} unrelated={j_none:.4f})")


def test_sketch_is_reproducible_across_processes():
    """Python's hash() of a str is salted per process; ours must not be."""
    code = (
        "import sys; sys.path.insert(0, r'%s');"
        "from src.robustness import minhash_sketch;"
        "print(int(minhash_sketch('ACGTTGCA'*400, 21, 64)[0]))"
        % str(Path(__file__).resolve().parent.parent)
    )
    outs = set()
    for seed in ("0", "1", "12345"):
        res = subprocess.run([sys.executable, "-c", code], capture_output=True,
                             text=True, env={**_env(), "PYTHONHASHSEED": seed})
        assert res.returncode == 0, res.stderr
        outs.add(res.stdout.strip())
    assert len(outs) == 1, f"sketch changed with PYTHONHASHSEED: {outs}"
    print("ok  sketch is independent of PYTHONHASHSEED")


def _env() -> dict:
    import os
    return dict(os.environ)


# --------------------------------------------------------------------------------------
# Clustering / representatives
# --------------------------------------------------------------------------------------


def test_single_linkage_clustering():
    labels = ["a", "b", "c", "d"]
    ident = np.array([
        [1.00, 0.999, 0.90, 0.50],
        [0.999, 1.00, 0.999, 0.50],   # b bridges a and c
        [0.90, 0.999, 1.00, 0.50],
        [0.50, 0.50, 0.50, 1.00],
    ])
    clusters = single_linkage_clusters(labels, ident, 0.999)
    assert clusters[0] == clusters[1] == clusters[2], clusters
    assert clusters[3] != clusters[0]
    assert clusters == [0, 0, 0, 1], clusters   # deterministic numbering
    assert single_linkage_clusters(labels, ident, 0.9999) == [0, 1, 2, 3]
    print("ok  single-linkage clustering (transitive, deterministic)")


def test_representative_choice_is_deterministic_and_prefers_resolved():
    man = pd.DataFrame([
        {"accession": "Z9", "ambiguous_fraction": 0.0, "length_bp": 100, "is_refseq": False},
        {"accession": "A1", "ambiguous_fraction": 0.0, "length_bp": 100, "is_refseq": True},
        {"accession": "B2", "ambiguous_fraction": 0.01, "length_bp": 999, "is_refseq": True},
    ])
    reps = choose_representatives(man, {"Z9": 0, "A1": 0, "B2": 0})
    # ambiguity first, then length, then RefSeq, then accession -> A1 beats Z9 on RefSeq.
    assert reps[0] == "A1", reps
    print("ok  representative choice prefers the best-resolved record, deterministically")


# --------------------------------------------------------------------------------------
# Anchoring, chaining, coverage
# --------------------------------------------------------------------------------------


def test_anchor_chaining_on_a_truncated_record():
    ref = _seq(5000)
    index = build_reference_index(ref, 25, 4)
    sub = ref[1000:2000]                       # a record covering ref[1000:2000)
    oriented, strand, anchors = orient_and_anchor(sub, index, 25, 4)
    chains = chain_anchors(anchors, 25)
    assert strand == "+" and oriented == sub
    assert len(chains) == 1, [len(c) for c in chains]
    assert chains[0].ref_min == 1000
    # Last sampled anchor starts at query 972 (step 4), so the chain reaches 1997;
    # it can never exceed the true right edge.
    assert 1990 <= chains[0].ref_max <= 2000, chains[0].ref_max
    print("ok  anchor chain recovers the exact reference span of a truncated record")


def test_reverse_complemented_record_is_reoriented():
    ref = _seq(5000)
    index = build_reference_index(ref, 25, 4)
    sub = revcomp(ref[1000:2000])
    oriented, strand, anchors = orient_and_anchor(sub, index, 25, 4)
    assert strand == "-"
    assert oriented == ref[1000:2000]
    chains = chain_anchors(anchors, 25)
    assert chains and chains[0].ref_min == 1000
    print("ok  reverse-complemented record is detected and reoriented")


def test_not_covered_is_not_scored_as_absent():
    """The whole point of part C: truncation must never look like variation."""
    ref = _seq(6000)
    index = build_reference_index(ref, 25, 4)
    record = ref[0:1000]                        # does not reach the guide at 3000
    _o, _s, anchors = orient_and_anchor(record, index, 25, 4)
    chains = chain_anchors(anchors, 25)
    amb = ambiguous_positions(record)
    status = classify_footprint([(3000, 3023)], chains, amb, 25, 2000)
    assert status == ST_NOT_COVERED, status
    print("ok  a record that stops short of the site is NOT_COVERED, not ABSENT")


def test_real_mismatch_is_scored_absent():
    ref = _seq(6000)
    index = build_reference_index(ref, 25, 4)
    record = list(ref[2000:4000])
    record[1000] = "A" if record[1000] != "A" else "C"     # SNP at ref 3000
    record = "".join(record)
    _o, _s, anchors = orient_and_anchor(record, index, 25, 4)
    chains = chain_anchors(anchors, 25)
    amb = ambiguous_positions(record)
    status = classify_footprint([(2990, 3013)], chains, amb, 25, 2000)
    assert status == ST_ABSENT, status
    print("ok  a single substitution inside a covered site is ABSENT (real variation)")


def test_assembly_gap_is_scored_unknown_not_absent():
    ref = _seq(6000)
    index = build_reference_index(ref, 25, 4)
    record = ref[2000:4000]
    record = record[:950] + "N" * 100 + record[1050:]       # N block over ref 2950-3050
    _o, _s, anchors = orient_and_anchor(record, index, 25, 4)
    chains = chain_anchors(anchors, 25)
    amb = ambiguous_positions(record)
    status = classify_footprint([(2990, 3013)], chains, amb, 25, 2000)
    assert status == ST_AMBIGUOUS, status
    print("ok  an N block over the site is UNKNOWN (assembly gap), not ABSENT")


def test_divergent_but_bracketed_region_is_still_covered():
    """A guide in a variable stretch must stay in the denominator.

    Dropping variable regions for lack of anchors would bias conservation upward,
    which is exactly the failure mode this analysis is meant to expose.
    """
    rng = random.Random(7)
    ref = _seq(6000, rng)
    index = build_reference_index(ref, 25, 4)
    record = list(ref[2000:4000])
    for off in range(980, 1030):                 # 50 bp of heavy divergence
        record[off] = rng.choice([b for b in "ACGT" if b != record[off]])
    record = "".join(record)
    _o, _s, anchors = orient_and_anchor(record, index, 25, 4)
    chains = chain_anchors(anchors, 25)
    amb = ambiguous_positions(record)
    status = classify_footprint([(2990, 3013)], chains, amb, 25, 2000)
    assert status == ST_ABSENT, status
    print("ok  a divergent-but-bracketed site stays in the denominator as ABSENT")


def test_large_assembly_gap_is_unknown_not_truncation():
    """An N block too big to chain across is a GAP, not a missing region.

    This is the distinction that decides whether "absent" means "this virus differs"
    or "this submitter could not resolve the region".
    """
    ref = _seq(12_000)
    index = build_reference_index(ref, 25, 4)
    record = ref[2000:4000] + "N" * 4000 + ref[8000:10_000]
    _o, _s, anchors = orient_and_anchor(record, index, 25, 4)
    chains = chain_anchors(anchors, 25)
    amb = ambiguous_positions(record)
    status = classify_footprint([(5000, 5023)], chains, amb, 25, 2000)
    assert status == ST_AMBIGUOUS, status
    print("ok  a chain-breaking N block over the site is UNKNOWN (gap), not NOT_COVERED")


def test_structural_difference_is_unknown_not_absent():
    """Record reaches both sides with no ambiguity, but nothing colinear crosses."""
    ref = _seq(12_000)
    index = build_reference_index(ref, 25, 4)
    record = ref[2000:4000] + ref[8000:10_000]        # the middle is simply not there
    _o, _s, anchors = orient_and_anchor(record, index, 25, 4)
    chains = chain_anchors(anchors, 25)
    amb = ambiguous_positions(record)
    status = classify_footprint([(5000, 5023)], chains, amb, 25, 2000)
    assert status == ST_WIDE_BRACKET, status
    print("ok  a structural difference over the site is UNKNOWN, not ABSENT")


def test_wide_bracket_is_unknown():
    ref = _seq(8000)
    index = build_reference_index(ref, 25, 4)
    # A record built from two distant reference blocks: nothing anchors in between,
    # so the site at ref 3000 is bracketed only across a huge span.
    record = ref[1000:2000] + ref[5000:6000]
    _o, _s, anchors = orient_and_anchor(record, index, 25, 4)
    chains = chain_anchors(anchors, 25, max_gap=6000, max_drift=60)
    amb = ambiguous_positions(record)
    status = classify_footprint([(3000, 3023)], chains, amb, 25, 200)
    assert status in (ST_WIDE_BRACKET, ST_NOT_COVERED), status
    print("ok  an implausibly wide bracket is UNKNOWN, not ABSENT")


def test_chain_handles_two_repeat_copies():
    """A duplicated segment must produce two chains, not one corrupt one."""
    unit = _seq(1500)
    ref = _seq(1000) + unit + _seq(1000) + unit + _seq(1000)
    index = build_reference_index(ref, 25, 4)
    record = unit
    _o, _s, anchors = orient_and_anchor(record, index, 25, 4)
    chains = chain_anchors(anchors, 25)
    spans = sorted((c.ref_min, c.ref_max) for c in chains if len(c) > 10)
    assert len(spans) == 2, spans
    assert spans[0][0] == 1000 and spans[1][0] == 3500, spans
    print("ok  a duplicated repeat yields one chain per copy")


def test_footprint_in_either_repeat_copy_counts():
    unit = _seq(1500)
    ref = _seq(1000) + unit + _seq(1000) + unit + _seq(1000)
    index = build_reference_index(ref, 25, 4)
    record = unit                                  # only anchors both copies
    _o, _s, anchors = orient_and_anchor(record, index, 25, 4)
    chains = chain_anchors(anchors, 25)
    amb = ambiguous_positions(record)
    # Guide recorded at the second copy only; the record must still count.
    status = classify_footprint([(3700, 3723)], chains, amb, 25, 2000)
    assert status == ST_ABSENT, status           # covered and resolved but not matching
    print("ok  a guide is covered if EITHER repeat copy is covered")


def test_bracket_requires_anchors_on_both_sides():
    ref = _seq(3000)
    index = build_reference_index(ref, 25, 4)
    record = ref[500:1500]
    _o, _s, anchors = orient_and_anchor(record, index, 25, 4)
    chain = chain_anchors(anchors, 25)[0]
    assert bracket_in_chain(chain, 900, 923, 25) is not None
    assert bracket_in_chain(chain, 480, 503, 25) is None      # left edge, no left anchor
    assert bracket_in_chain(chain, 1490, 1513, 25) is None    # right edge
    print("ok  bracketing requires anchors strictly on both sides")


# --------------------------------------------------------------------------------------
# Bookkeeping helpers
# --------------------------------------------------------------------------------------


def test_parse_footprints_handles_repeat_copies():
    cand = pd.DataFrame([
        {"guide_id": "RL2_10+", "ref_start": 10, "ref_end": 32,
         "ref_all_positions": "10-32(+);5000-5022(-)"},
        {"guide_id": "UL30_7-", "ref_start": 7, "ref_end": 29, "ref_all_positions": ""},
    ])
    fps = parse_footprints(cand)
    assert fps[0] == [(9, 32), (4999, 5022)], fps[0]
    assert fps[1] == [(6, 29)], fps[1]
    print("ok  reference footprints parsed, including repeat copies")


def test_score_benchmark_both_strands():
    target = "ACGTACGTACGTACGTACGT" + "TGGAGT"
    genome = "TTTT" + target + "AAAA"
    assert score_benchmark(genome, {"g": target}) == {"g"}
    assert score_benchmark(revcomp(genome), {"g": target}) == {"g"}
    assert score_benchmark("ACGT" * 20, {"g": target}) == set()
    print("ok  benchmark guides are matched on either strand")


def test_compare_denominators_counts_the_headline():
    genome = pd.DataFrame({
        "guide_id": ["g1", "g2", "g3", "g4"],
        "conservation_strict": [1.0, 1.0, 1.0, 0.5],
    })
    gene = pd.DataFrame({
        "guide_id": ["g1", "g2", "g3", "g4"],
        "n_records_covering": [100, 100, 3, 100],   # g3 has too few records
        "conservation_gene_level": [1.0, 0.60, 1.0, 0.40],
    })
    out = compare_denominators(genome, gene, min_records=10)
    assert out["n_guides_with_enough_records"] == 3, out
    assert out["n_perfect_at_183"] == 2, out       # g1, g2 (g3 excluded, too few records)
    assert out["n_perfect_at_183_below_0.95_gene_level"] == 1, out
    assert out["n_perfect_at_183_below_0.70_gene_level"] == 1, out
    assert out["n_perfect_at_183_still_perfect"] == 1, out
    print("ok  denominator comparison counts the headline correctly")


def test_compare_denominators_can_select_a_corpus_tier():
    """Pooling record classes can hide which one carries the result."""
    genome = pd.DataFrame({"guide_id": ["g1", "g2"], "conservation_strict": [1.0, 1.0]})
    gene = pd.DataFrame({
        "guide_id": ["g1", "g2"],
        "n_records_covering": [100, 100],
        "conservation_gene_level": [0.99, 0.99],
        "n_records_covering_subgenomic": [40, 2],     # g2 has almost no independent data
        "conservation_subgenomic": [0.50, 1.00],
    })
    pooled = compare_denominators(genome, gene, 10)
    assert pooled["n_perfect_at_183"] == 2
    assert pooled["n_perfect_at_183_below_0.95_gene_level"] == 0

    sub = compare_denominators(genome, gene, 10,
                               "n_records_covering_subgenomic",
                               "conservation_subgenomic")
    assert sub["n_guides_with_enough_records"] == 1, sub   # g2 dropped, too few records
    assert sub["n_perfect_at_183_below_0.95_gene_level"] == 1, sub
    print("ok  per-tier comparison exposes what the pooled corpus hides")


def test_status_codes_order_is_load_bearing():
    """classify_footprint returns the most informative status; ABSENT must win."""
    assert ST_PRESENT < ST_ABSENT < ST_AMBIGUOUS < ST_WIDE_BRACKET < ST_NOT_COVERED
    print("ok  status ordering (ABSENT beats UNKNOWN beats NOT_COVERED)")


def main() -> int:
    tests = [v for k, v in sorted(globals().items())
             if k.startswith("test_") and callable(v)]
    for fn in tests:
        fn()
    print(f"\nAll {len(tests)} robustness tests passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
