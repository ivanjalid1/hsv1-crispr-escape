"""Offline unit tests for the load-bearing correctness claims. No network needed.

Run with:  python -m tests.test_core      (or: pytest tests/)

The most important test here is `test_scan_matches_bruteforce`: conservation.py
replaces a naive substring search with a "GG"/"CC" jump scan for speed, and that
optimisation is only legitimate if it returns *exactly* the same answer. This test
proves the equivalence on randomised sequences.
"""

from __future__ import annotations

import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.common import gc_fraction, has_polyt, max_homopolymer_run, revcomp
from src.conservation import build_lookup, scan_genome
from src.extract_guides import enumerate_guides_in_segment


def test_revcomp():
    assert revcomp("ATGC") == "GCAT"
    assert revcomp("AAAGGG") == "CCCTTT"
    assert revcomp(revcomp("ACGTACGTNN")) == "ACGTACGTNN"
    print("ok  revcomp")


def test_flags():
    assert gc_fraction("GGCC") == 1.0
    assert gc_fraction("ATAT") == 0.0
    assert max_homopolymer_run("AATTTTGC") == 4
    assert max_homopolymer_run("ACGT") == 1
    assert has_polyt("ACGTTTTACG") is True
    assert has_polyt("ACGTTTACG") is False
    print("ok  gc / homopolymer / polyT flags")


def test_enumerate_plus_and_minus():
    # 20 nt protospacer + AGG PAM on the plus strand.
    proto = "ACGTACGTACGTACGTACGT"
    seq = "TT" + proto + "AGG" + "TT"
    hits = list(enumerate_guides_in_segment(seq, 0, len(seq)))
    plus = [h for h in hits if h[0] == "+"]
    assert len(plus) == 1, plus
    strand, s1, e1, p, pam, t23 = plus[0]
    assert p == proto and pam == "AGG" and t23 == proto + "AGG"
    assert s1 == 3 and e1 == 25          # 1-based inclusive over the 23-mer
    assert seq[s1 - 1:e1] == t23
    print("ok  plus-strand enumeration")

    # Same site written on the minus strand: the plus strand now starts with CCT.
    rc = revcomp(proto + "AGG")
    seq2 = "TT" + rc + "TT"
    minus = [h for h in enumerate_guides_in_segment(seq2, 0, len(seq2)) if h[0] == "-"]
    assert len(minus) == 1, minus
    _, s1, e1, p, pam, t23 = minus[0]
    assert t23 == proto + "AGG" and pam == "AGG"
    assert revcomp(seq2[s1 - 1:e1]) == t23
    print("ok  minus-strand enumeration")


def test_segment_bounds_are_respected():
    proto = "ACGTACGTACGTACGTACGT"
    seq = "TT" + proto + "AGG" + "TT"
    # Window that cannot contain the whole 23-mer must yield nothing.
    assert list(enumerate_guides_in_segment(seq, 5, 20)) == []
    print("ok  segment bounds respected")


def test_ambiguous_bases_rejected():
    proto = "ACGTACGTACGTACGTACGN"
    seq = "TT" + proto + "AGG" + "TT"
    assert list(enumerate_guides_in_segment(seq, 0, len(seq))) == []
    print("ok  ambiguous bases rejected")


def test_scan_matches_bruteforce():
    """The GG/CC jump scan must equal a naive both-strand substring search."""
    rng = random.Random(20260907)
    for trial in range(8):
        genome = "".join(rng.choice("ACGT") for _ in range(60_000))

        # Guides: a mix of real substrings of the genome and random decoys.
        targets = {}
        for i in range(150):
            start = rng.randrange(0, len(genome) - 23)
            cand = genome[start:start + 23]
            if cand[21:23] == "GG":
                targets[f"g{i}"] = cand
            elif revcomp(cand)[21:23] == "GG":
                targets[f"g{i}"] = revcomp(cand)
        for i in range(50):
            decoy = "".join(rng.choice("ACGT") for _ in range(21)) + "GG"
            targets[f"d{i}"] = decoy

        fwd, rev = build_lookup(targets)
        fast = scan_genome(genome, fwd, rev)
        brute = {gid for gid, t in targets.items()
                 if (t in genome) or (revcomp(t) in genome)}
        assert fast == brute, (
            f"trial {trial}: jump scan and brute force disagree; "
            f"only-fast={sorted(fast - brute)} only-brute={sorted(brute - fast)}"
        )
    print("ok  jump scan == brute force (8 randomised genomes)")


def test_scan_finds_guide_on_reverse_strand():
    """A guide present only as a reverse-complement must still be detected."""
    proto = "ACGTTGCAACGTTGCAACGT"
    t23 = proto + "TGG"
    genome = "AAAA" + revcomp(t23) + "TTTT"
    fwd, rev = build_lookup({"x": t23})
    assert scan_genome(genome, fwd, rev) == {"x"}
    print("ok  reverse-strand presence detected")


def main() -> int:
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in tests:
        fn()
    print(f"\nAll {len(tests)} tests passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
