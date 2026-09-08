"""Offline unit tests for stage 7 (human off-target screening).

Run with:  python tests/test_offtarget.py      (or: pytest tests/)

The load-bearing test is `test_fast_matches_bruteforce_chr21`. `src/offtarget.py`
replaces an exhaustive mismatch scan with a PAM prefilter plus a pigeonhole chunk
index; that is legitimate only if it returns *exactly* the same hit set. This file
proves it three ways, each implementation structurally independent of the one below
it:

    naive_scan_python        pure Python, character by character, no numpy at all
        validates
    bruteforce_scan_sequence numpy, mismatch counts at EVERY genomic position by one
                             whole-array slice comparison per protospacer base, PAM
                             applied afterwards, no chunk index, no gathers
        validates
    scan_sequence            the fast path actually used for the genome-wide screen

The chromosome-21 test runs the middle implementation over all 46.7 Mb of chr21, both
strands, for every guide the screen reports on, and asserts set identity. It is the
same standard `tests/test_nuclease.py::test_scan_matches_bruteforce_every_pam` applies
to the stage-3 PAM scanner.

No network access. Tests that need the cached human genome skip *loudly* if
`data/genome/` has not been populated, so this file is runnable on a fresh clone.
"""

from __future__ import annotations

import os
import random
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.common import revcomp
from src.nuclease import IUPAC_CODES
from src.offtarget import (
    BASE_CODE,
    CACHE_INDEX,
    GenomeCache,
    Guide,
    MAX_CHUNK_WIDTH,
    IntervalIndex,
    bruteforce_scan_sequence,
    chunk_layout,
    decode_sequence,
    encode_sequence,
    load_guides,
    materialise_hits,
    naive_scan_python,
    pam_position_mask,
    revcomp_codes,
    scan_sequence,
)

RNG = random.Random(20260908)


def _rand_seq(n: int, gc: float = 0.5, rng: random.Random | None = None) -> str:
    rng = rng or RNG
    at, gcb = "AT", "GC"
    return "".join(rng.choice(gcb) if rng.random() < gc else rng.choice(at)
                   for _ in range(n))


def _guides(seqs: list[str]) -> list[Guide]:
    return [Guide(guide_id=f"g{i}", label=f"g{i}", group="test", protospacer=s)
            for i, s in enumerate(seqs)]


# --------------------------------------------------------------------------------------
# encoding
# --------------------------------------------------------------------------------------


def test_encoding_roundtrip():
    seq = "ACGTNNACGTacgtRYKM"
    codes = encode_sequence(seq)
    assert codes.tolist() == [0, 1, 2, 3, 4, 4, 0, 1, 2, 3, 0, 1, 2, 3, 4, 4, 4, 4]
    assert decode_sequence(codes) == "ACGTNNACGTACGTNNNN"
    print("ok test_encoding_roundtrip")


def test_revcomp_codes_agrees_with_string_revcomp():
    for _ in range(200):
        seq = _rand_seq(RNG.randint(1, 60))
        assert decode_sequence(revcomp_codes(encode_sequence(seq))) == revcomp(seq)
    # N must survive reverse complementation as N, not silently become a base.
    assert decode_sequence(revcomp_codes(encode_sequence("ACNGT"))) == "ACNGT"
    print("ok test_revcomp_codes_agrees_with_string_revcomp")


# --------------------------------------------------------------------------------------
# the pigeonhole property itself
# --------------------------------------------------------------------------------------


def test_chunk_layout_tiles_the_protospacer():
    for L in range(5, 40):
        for m in range(0, min(6, L - 1)):
            layout = chunk_layout(L, m)
            # At least m+1 chunks (pigeonhole), possibly more so that every chunk
            # stays narrow enough for a dense 4**width lookup table.
            assert len(layout) >= m + 1
            assert max(w for _, w in layout) <= MAX_CHUNK_WIDTH
            assert layout[0][0] == 0
            covered = 0
            for off, w in layout:
                assert off == covered
                covered += w
            assert covered == L, (L, m, layout)
    print("ok test_chunk_layout_tiles_the_protospacer")


def test_pigeonhole_guarantee_holds_by_exhaustion():
    """With <= m mismatches, at least one of the m+1 chunks must match exactly.

    Asserted directly rather than argued: mutate a random protospacer at every
    possible combination of up to m positions and check the invariant each time.
    """
    from itertools import combinations

    for L, m in ((20, 4), (21, 4), (20, 3), (12, 2)):
        layout = chunk_layout(L, m)
        for _ in range(30):
            proto = _rand_seq(L)
            for k in range(m + 1):
                for positions in combinations(range(L), k):
                    mutant = list(proto)
                    for p in positions:
                        mutant[p] = RNG.choice([b for b in "ACGT" if b != proto[p]])
                    mutant = "".join(mutant)
                    exact = [c for c, (off, w) in enumerate(layout)
                             if mutant[off:off + w] == proto[off:off + w]]
                    assert exact, (L, m, positions)
                if k >= 2:      # the full combinatorial sweep is only cheap for small k
                    break
    print("ok test_pigeonhole_guarantee_holds_by_exhaustion")


# --------------------------------------------------------------------------------------
# PAM matching
# --------------------------------------------------------------------------------------


def _naive_pam_positions(seq: str, spacer_length: int, pam: str) -> set[int]:
    """Independent PAM matcher written straight from the IUPAC table."""
    sets = [set(IUPAC_CODES[c]) for c in pam.upper()]
    span = spacer_length + len(pam)
    out = set()
    for s in range(0, len(seq) - span + 1):
        window = seq[s + spacer_length:s + span]
        if all(ch in st for ch, st in zip(window, sets)):
            out.add(s)
    return out


def test_pam_position_mask_matches_naive():
    for pam in ("NNGRRT", "NNGRRN", "NGG", "TTTV", "ACGT"):
        for gc in (0.3, 0.5, 0.75):
            seq = _rand_seq(4000, gc=gc)
            # Sprinkle Ns: an N must never satisfy any IUPAC code, N included.
            s = list(seq)
            for _ in range(120):
                s[RNG.randrange(len(s))] = "N"
            seq = "".join(s)
            arr = encode_sequence(seq)
            fast = set(np.flatnonzero(pam_position_mask(arr, 20, pam)).tolist())
            assert fast == _naive_pam_positions(seq, 20, pam), (pam, gc)
    print("ok test_pam_position_mask_matches_naive")


def test_ambiguous_bases_never_match():
    """An N in the protospacer is a mismatch; an N in the PAM disqualifies the site."""
    proto = "ACGTACGTACGTACGTACGT"
    guides = _guides([proto])
    # Perfect site, then the same site with one N in the protospacer and one in the PAM.
    ok = "T" * 5 + proto + "CAGAGT" + "T" * 5
    hits = scan_sequence(encode_sequence(ok), guides, pam="NNGRRT", max_mismatches=0)
    assert [(0, 5, 0)] == hits, hits

    n_in_proto = "T" * 5 + "N" + proto[1:] + "CAGAGT" + "T" * 5
    hits0 = scan_sequence(encode_sequence(n_in_proto), guides, pam="NNGRRT",
                          max_mismatches=0)
    hits1 = scan_sequence(encode_sequence(n_in_proto), guides, pam="NNGRRT",
                          max_mismatches=1)
    assert hits0 == [] and hits1 == [(0, 5, 1)], (hits0, hits1)

    n_in_pam = "T" * 5 + proto + "CAGAGN" + "T" * 5
    for pam in ("NNGRRT", "NNGRRN"):
        assert scan_sequence(encode_sequence(n_in_pam), guides, pam=pam,
                             max_mismatches=4) == [], pam
    print("ok test_ambiguous_bases_never_match")


# --------------------------------------------------------------------------------------
# three-level equivalence on synthetic sequence
# --------------------------------------------------------------------------------------


def test_python_naive_matches_numpy_bruteforce():
    for _ in range(12):
        seq = _rand_seq(6000, gc=RNG.choice([0.35, 0.5, 0.68]))
        guides = _guides([_rand_seq(20) for _ in range(3)])
        # Plant a few near-matches so the test is not just asserting "no hits".
        s = list(seq)
        for gi, g in enumerate(guides):
            for k in (0, 1, 3, 4, 5):
                start = 200 + 700 * gi + 120 * k
                site = list(g.protospacer)
                for p in RNG.sample(range(20), k):
                    site[p] = RNG.choice([b for b in "ACGT" if b != site[p]])
                s[start:start + 26] = list("".join(site) + "CAGAGT")
        seq = "".join(s)
        arr = encode_sequence(seq)
        for pam in ("NNGRRT", "NNGRRN"):
            a = set(naive_scan_python(seq, guides, pam=pam, max_mismatches=4))
            b = set(bruteforce_scan_sequence(arr, guides, pam=pam, max_mismatches=4))
            assert a == b, sorted(a ^ b)[:10]
    print("ok test_python_naive_matches_numpy_bruteforce")


def test_fast_matches_bruteforce_synthetic():
    """Set identity on randomised sequence, over PAMs, spacer lengths and thresholds."""
    for pam in ("NNGRRT", "NNGRRN"):
        for L in (20, 21):
            for max_mm in (0, 2, 4):
                for gc in (0.35, 0.5, 0.72):
                    seq = _rand_seq(30000, gc=gc)
                    guides = _guides([_rand_seq(L) for _ in range(6)])
                    s = list(seq)
                    for gi, g in enumerate(guides):
                        for k in range(0, 6):
                            start = 100 + 900 * gi + 130 * k
                            site = list(g.protospacer)
                            for p in RNG.sample(range(L), k):
                                site[p] = RNG.choice([b for b in "ACGT" if b != site[p]])
                            s[start:start + L + 6] = list("".join(site) + "CAGAGT")
                    for _ in range(300):
                        s[RNG.randrange(len(s))] = "N"
                    arr = encode_sequence("".join(s))
                    fast = set(scan_sequence(arr, guides, pam=pam,
                                             max_mismatches=max_mm, block=4096))
                    slow = set(bruteforce_scan_sequence(arr, guides, pam=pam,
                                                        max_mismatches=max_mm))
                    assert fast == slow, (pam, L, max_mm, gc, sorted(fast ^ slow)[:10])
    print("ok test_fast_matches_bruteforce_synthetic")


def test_fast_matches_bruteforce_low_complexity():
    """Homopolymer / microsatellite sequence is the adversarial case for a chunk index."""
    for motif in ("A", "AT", "CAG", "GGCGC", "TTTTA"):
        seq = (motif * (30000 // len(motif)))[:30000]
        s = list(seq)
        for _ in range(400):
            s[RNG.randrange(len(s))] = RNG.choice("ACGTN")
        arr = encode_sequence("".join(s))
        guides = _guides([(motif * 20)[:20], _rand_seq(20), "A" * 20,
                          "CAGCAGCAGCAGCAGCAGCA"])
        for pam in ("NNGRRT", "NNGRRN"):
            fast = set(scan_sequence(arr, guides, pam=pam, max_mismatches=4, block=1024))
            slow = set(bruteforce_scan_sequence(arr, guides, pam=pam, max_mismatches=4))
            assert fast == slow, (motif, pam, len(fast), len(slow))
    print("ok test_fast_matches_bruteforce_low_complexity")


def test_block_size_does_not_change_the_answer():
    seq = _rand_seq(50000, gc=0.55)
    arr = encode_sequence(seq)
    guides = _guides([_rand_seq(20) for _ in range(4)])
    ref = None
    for block in (17, 1000, 1 << 16, 1 << 23):
        got = set(scan_sequence(arr, guides, pam="NNGRRN", max_mismatches=4, block=block))
        if ref is None:
            ref = got
        assert got == ref, block
    print("ok test_block_size_does_not_change_the_answer")


# --------------------------------------------------------------------------------------
# coordinate mapping
# --------------------------------------------------------------------------------------


def test_minus_strand_coordinates_round_trip():
    """A site planted on the minus strand must be reported at the right plus-strand
    coordinate, with the protospacer read 5'->3' on the guide's own strand."""
    proto = "GTACCCGACGGCCCCCGCGT"
    pam = "CGGAGT"
    footprint = proto + pam
    left = _rand_seq(500)
    right = _rand_seq(500)
    plus = left + revcomp(footprint) + right          # site is on the minus strand
    arr = encode_sequence(plus)
    rc = revcomp_codes(arr)
    guides = _guides([proto])
    hits = scan_sequence(rc, guides, pam="NNGRRT", max_mismatches=0)
    assert len(hits) == 1, hits
    rows = materialise_hits(hits, guides, rc, "chrTest", "-", arr.size, len(pam))
    r = rows[0]
    assert r["strand"] == "-"
    assert r["start"] == len(left) + 1, r          # 1-based plus-strand start
    assert r["end"] == len(left) + len(footprint)
    assert r["protospacer_genomic"] == proto
    assert r["pam_genomic"] == pam
    assert r["mismatches"] == 0
    assert plus[r["start"] - 1:r["end"]] == revcomp(footprint)
    print("ok test_minus_strand_coordinates_round_trip")


def test_plus_strand_coordinates_round_trip():
    proto = "AATCCTAGACACGCACCGCC"
    pam = "AGGAGT"
    left = _rand_seq(300)
    plus = left + proto + pam + _rand_seq(300)
    arr = encode_sequence(plus)
    guides = _guides([proto])
    hits = scan_sequence(arr, guides, pam="NNGRRT", max_mismatches=0)
    rows = materialise_hits(hits, guides, arr, "chrTest", "+", arr.size, len(pam))
    r = rows[0]
    assert r["start"] == len(left) + 1 and r["end"] == len(left) + 26
    assert plus[r["start"] - 1:r["end"]] == proto + pam
    assert r["seed_mismatches"] == 0
    print("ok test_plus_strand_coordinates_round_trip")


def test_mismatch_positions_are_reported_1_based():
    proto = "ACGTACGTACGTACGTACGT"
    assert proto[0] == "A" and proto[12] == "A"
    mutated = "T" + proto[1:12] + "T" + proto[13:]     # positions 1 and 13 (1-based)
    plus = _rand_seq(100) + mutated + "CAGAGT" + _rand_seq(100)
    arr = encode_sequence(plus)
    guides = _guides([proto])
    hits = scan_sequence(arr, guides, pam="NNGRRT", max_mismatches=4)
    rows = materialise_hits(hits, guides, arr, "chrTest", "+", arr.size, 6)
    assert len(rows) == 1
    assert rows[0]["mismatch_positions"] == "1,13", rows[0]
    assert rows[0]["mismatches"] == 2
    assert rows[0]["seed_mismatches"] == 1        # only position 13 is in the last 12 nt
    print("ok test_mismatch_positions_are_reported_1_based")


# --------------------------------------------------------------------------------------
# interval index
# --------------------------------------------------------------------------------------


def test_interval_index_overlaps():
    chrom = np.array(["1", "1", "1", "2", "2"], dtype="U")
    start = np.array([100, 500, 120, 10, 1_000_000], dtype=np.int64)
    end = np.array([200, 600, 100_000, 20, 1_000_050], dtype=np.int64)
    name = np.array(["a", "b", "long", "c", "d"], dtype="U")
    extra = np.array(["x"] * 5, dtype="U")
    ix = IntervalIndex(chrom, start, end, name, extra)

    def names(c, s, e):
        return sorted(str(x) for x in ix.name[ix.overlaps(c, s, e)])

    assert names("1", 150, 160) == ["a", "long"]
    assert names("1", 250, 260) == ["long"]          # bounded backward scan must find it
    assert names("1", 90, 99) == []
    assert names("1", 100, 100) == ["a"]
    assert names("1", 600, 600) == ["b", "long"]
    assert names("2", 15, 15) == ["c"]
    assert names("3", 1, 10**9) == []
    print("ok test_interval_index_overlaps")


# --------------------------------------------------------------------------------------
# guide selection
# --------------------------------------------------------------------------------------


def test_guide_selection_is_the_documented_rule():
    from src.offtarget import BENCHMARK_POOL
    if not BENCHMARK_POOL.exists():
        print("SKIP test_guide_selection_is_the_documented_rule: run stage 6 first")
        return
    guides = load_guides(spacer_lengths=(20,))
    labels = {g.label for g in guides}
    for published in ("ICP0g1", "ICP0g2", "ICP27g1", "ICP27g2"):
        assert published in labels, f"{published} missing from the screened set"
    ids = {g.guide_id for g in guides}
    assert "RL2_3441+" in ids, "the headline recommendation must be screened"
    assert all(g.spacer_length == 20 for g in guides)
    assert len({g.protospacer for g in guides}) == len(guides), "duplicate protospacers"
    n_icp0 = sum(1 for g in guides if g.group == "icp0_alternative")
    assert n_icp0 == 12, f"expected the 12 stage-6 ICP0 alternatives, got {n_icp0}"
    print(f"ok test_guide_selection_is_the_documented_rule ({len(guides)} guides)")


def test_21nt_variants_extend_the_20nt_spacer():
    from src.offtarget import BENCHMARK_POOL
    if not BENCHMARK_POOL.exists():
        print("SKIP test_21nt_variants_extend_the_20nt_spacer: run stage 6 first")
        return
    both = load_guides(spacer_lengths=(20, 21))
    by_id = {g.guide_id: g for g in both}
    n = 0
    for gid, g in by_id.items():
        if not gid.endswith("|21nt"):
            continue
        base = by_id[gid[:-5]]
        assert g.spacer_length == 21
        assert g.protospacer[1:] == base.protospacer, gid
        n += 1
    assert n > 0, "no 21-nt variants were built (is the HSV-1 reference cached?)"
    print(f"ok test_21nt_variants_extend_the_20nt_spacer ({n} variants)")


# --------------------------------------------------------------------------------------
# THE load-bearing test: whole of chr21, every guide, fast path == brute force
# --------------------------------------------------------------------------------------


def _load_chrom(name: str):
    if not CACHE_INDEX.exists():
        return None
    import json
    index = json.loads(CACHE_INDEX.read_text(encoding="utf-8"))
    if not any(s["name"] == name for s in index["sequences"]):
        return None
    return GenomeCache(index).get(name)


def test_fast_matches_bruteforce_chr21():
    """Set identity over an entire real human chromosome, for every guide screened.

    This is the test that makes the genome-wide numbers usable. The brute force here
    computes a mismatch count at *every* position of chr21 by one whole-array slice
    comparison per protospacer base, applies the PAM constraint afterwards, and knows
    nothing about chunks, bitmasks or PAM-position gathers. If the fast path and the
    brute force disagree even once, the fast path is wrong.

    Set OFFTARGET_TEST_CHROM to run it on a different sequence (e.g. 22).
    """
    name = os.environ.get("OFFTARGET_TEST_CHROM", "21")
    arr = _load_chrom(name)
    if arr is None:
        print(f"SKIP test_fast_matches_bruteforce_chr21: chromosome {name} is not in "
              "data/genome/. Populate it with `python src/offtarget.py --stage fetch`.")
        return
    from src.offtarget import BENCHMARK_POOL
    if not BENCHMARK_POOL.exists():
        print("SKIP test_fast_matches_bruteforce_chr21: run stage 6 first")
        return

    guides = load_guides(spacer_lengths=(20, 21))
    by_len: dict[int, list[Guide]] = {}
    for g in guides:
        by_len.setdefault(g.spacer_length, []).append(g)

    rc = revcomp_codes(arr)
    total_hits = 0
    for L, gl in sorted(by_len.items()):
        for pam in ("NNGRRN",):
            for strand, a in (("+", arr), ("-", rc)):
                t0 = time.time()
                fast = set(scan_sequence(a, gl, pam=pam, max_mismatches=4))
                t_fast = time.time() - t0
                t0 = time.time()
                slow = set(bruteforce_scan_sequence(a, gl, pam=pam, max_mismatches=4))
                t_slow = time.time() - t0
                assert fast == slow, (
                    f"chr{name} {strand} L={L} {pam}: fast and brute force differ on "
                    f"{len(fast ^ slow)} tuples, e.g. {sorted(fast ^ slow)[:5]}")
                total_hits += len(fast)
                print(f"   chr{name} {strand} L={L:2d} {pam}: {len(fast):6d} hits "
                      f"identical  (fast {t_fast:6.1f}s, brute force {t_slow:7.1f}s, "
                      f"{t_slow / max(t_fast, 1e-9):5.1f}x)")
    print(f"ok test_fast_matches_bruteforce_chr21 "
          f"({arr.size:,} bp, {len(guides)} guides, {total_hits} hits, set-identical)")


# --------------------------------------------------------------------------------------


def main() -> int:
    tests = [v for k, v in sorted(globals().items())
             if k.startswith("test_") and callable(v)]
    t0 = time.time()
    for fn in tests:
        fn()
    print(f"\n{len(tests)} tests passed in {time.time() - t0:.1f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
