"""Offline unit tests for the nuclease/PAM model and the generalised scanner.

Run with:  python tests/test_nuclease.py      (or: pytest tests/)

The load-bearing test here is `test_scan_matches_bruteforce_every_pam`, the
generalisation of `tests/test_core.py::test_scan_matches_bruteforce`. Stage 3
replaces a naive both-strand substring search with a jump scan between literal PAM
anchors; for SpCas9 those anchors are "GG"/"CC", for SaCas9 NNGRRT they are the "G"
and the "T" of the PAM. The optimisation is legitimate only if it returns *exactly*
the same set for every supported PAM, which is proved here on randomised sequences,
on adversarial hand-built sequences, and -- when the cached reference is present --
on the real HSV-1 genome.

No network access. The real-genome tests skip (loudly) if data/raw/reference is
absent, so the file is runnable on a fresh clone.
"""

from __future__ import annotations

import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.common import REF_DIR, revcomp
from src.conservation import build_lookup, scan_genome, scan_genome_naive
from src.extract_guides import enumerate_guides_in_segment
from src.nuclease import (
    IUPAC_CODES,
    NUCLEASES,
    SACAS9,
    SPCAS9,
    Nuclease,
    get_nuclease,
    iupac_revcomp,
)

# Every PAM grammar the tests exercise. "NNN" and "NNNN" have no literal position at
# all, which forces the scanner's exhaustive fallback path; "TTTV" is a 3'-PAM stand-in
# with a leading literal run, exercising a non-zero anchor at the low end of the PAM.
TEST_SPECS = [
    Nuclease(name="spcas9", spacer_length=20, pam="NGG"),
    Nuclease(name="spcas9", spacer_length=20, pam="NAG"),
    Nuclease(name="sacas9", spacer_length=21, pam="NNGRRT"),
    Nuclease(name="sacas9", spacer_length=21, pam="NNGRRN"),
    Nuclease(name="sacas9", spacer_length=20, pam="NNGRRT"),
    Nuclease(name="test", spacer_length=8, pam="NNN"),
    Nuclease(name="test", spacer_length=8, pam="TTTV"),
    Nuclease(name="test", spacer_length=6, pam="ACGT"),
]


# --------------------------------------------------------------------------------------
# helpers -- independent reference implementations
# --------------------------------------------------------------------------------------


def _naive_pam_ok(seq: str, pattern: str) -> bool:
    """Independent IUPAC matcher, written from the code table rather than reusing
    Nuclease's precomputed sets."""
    if len(seq) != len(pattern):
        return False
    return all(base in IUPAC_CODES[code] for base, code in zip(seq, pattern))


def _naive_enumerate(genome: str, spec: Nuclease) -> set[tuple[str, int, str]]:
    """(strand, 1-based start, target site) by brute force over both strands.

    Deliberately written the slow, obvious way: build the reverse complement of the
    whole window and test the PAM at its 3' end, rather than testing a mirrored
    pattern at the 5' end of the plus-strand window.
    """
    out: set[tuple[str, int, str]] = set()
    k = spec.target_length
    for i in range(0, len(genome) - k + 1):
        window = genome[i:i + k]
        if any(ch not in "ACGT" for ch in window):
            continue
        if _naive_pam_ok(window[spec.spacer_length:], spec.pam):
            out.add(("+", i + 1, window))
        rc = revcomp(window)
        if _naive_pam_ok(rc[spec.spacer_length:], spec.pam):
            out.add(("-", i + 1, rc))
    return out


def _naive_presence(seq: str, targets: dict[str, str]) -> set[str]:
    """Plain `in` test on both strands -- the definition of 'present' in this project."""
    rc = revcomp(seq)
    return {gid for gid, t in targets.items() if t in seq or t in rc}


def _random_seq(rng: random.Random, n: int, alphabet: str = "ACGT") -> str:
    return "".join(rng.choice(alphabet) for _ in range(n))


def _reference_sequence() -> str | None:
    """The cached HSV-1 reference, or None when data/raw is not populated."""
    path = REF_DIR / "NC_001806.gb"
    if not path.is_file():
        return None
    from Bio import SeqIO

    return str(SeqIO.read(path, "genbank").seq).upper()


# --------------------------------------------------------------------------------------
# IUPAC model
# --------------------------------------------------------------------------------------


def test_iupac_revcomp():
    assert iupac_revcomp("NGG") == "CCN"
    assert iupac_revcomp("NNGRRT") == "AYYCNN"
    assert iupac_revcomp("NNGRRN") == "NYYCNN"
    for pattern in ("NGG", "NNGRRT", "TTTV", "ACGT", "NNN"):
        assert iupac_revcomp(iupac_revcomp(pattern)) == pattern
    # The reverse-complement pattern must admit exactly the reverse complements of
    # the sequences the original admits.
    for pattern in ("NGG", "NNGRRT", "TTTV"):
        rcp = iupac_revcomp(pattern)
        rng = random.Random(1)
        for _ in range(400):
            s = _random_seq(rng, len(pattern))
            assert _naive_pam_ok(s, pattern) == _naive_pam_ok(revcomp(s), rcp)
    print("ok  IUPAC reverse complement")


def test_pam_matching():
    assert SPCAS9.pam_matches("AGG") and SPCAS9.pam_matches("TGG")
    assert not SPCAS9.pam_matches("AGA") and not SPCAS9.pam_matches("AG")
    # NNGRRT: R = A/G.
    for pam in ("AAGAAT", "CCGAGT", "TTGGGT", "GGGGGT"):
        assert SACAS9.pam_matches(pam), pam
    for pam in ("AAGACT", "AACAAT", "AAGAAC", "AAGAA"):
        assert not SACAS9.pam_matches(pam), pam
    # NNGRRN admits the NNGRRT sites plus any final base.
    permissive = get_nuclease("sacas9", pam="NNGRRN")
    assert permissive.pam_matches("AAGAAC") and permissive.pam_matches("GGGGGT")
    assert not permissive.pam_matches("AACAAC")
    # Amrani et al. 2024 Table 1 PAMs must all satisfy NNGRRT.
    for pam in ("CGGAGT", "CAGAGT", "AGGAGT", "GGGGGT"):
        assert SACAS9.pam_matches(pam), pam
    print("ok  IUPAC PAM matching (NGG, NNGRRT, NNGRRN)")


def test_registry_and_geometry():
    assert set(NUCLEASES) == {"spcas9", "sacas9"}
    assert (SPCAS9.spacer_length, SPCAS9.pam, SPCAS9.target_length) == (20, "NGG", 23)
    assert (SACAS9.spacer_length, SACAS9.pam, SACAS9.target_length) == (21, "NNGRRT", 27)
    assert get_nuclease("spcas9").tag == "spcas9"
    assert get_nuclease("sacas9").tag == "sacas9"
    assert get_nuclease("sacas9", pam="NNGRRN").tag == "sacas9-nngrrn"
    assert get_nuclease("sacas9", spacer_length=20).tag == "sacas9-20nt"
    for bad in ("NXG", "", "N-G"):
        try:
            Nuclease(name="x", spacer_length=20, pam=bad)
        except ValueError:
            pass
        else:
            raise AssertionError(f"PAM {bad!r} should have been rejected")
    try:
        Nuclease(name="cas12a", spacer_length=23, pam="TTTV", pam_side="5prime")
    except NotImplementedError:
        pass
    else:
        raise AssertionError("5'-PAM nucleases must be rejected, not half-supported")
    print("ok  registry, geometry and input validation")


def test_scan_anchors_are_invariants():
    """Every enumerated target site must contain each anchor at its stated offset."""
    rng = random.Random(4242)
    genome = _random_seq(rng, 30_000)
    for spec in TEST_SPECS:
        sites = _naive_enumerate(genome, spec)
        assert sites, spec.pam
        for literal, offset in spec.scan_anchors:
            for _strand, _start, target in sites:
                assert target[offset:offset + len(literal)] == literal, (spec.pam, target)
    print("ok  PAM anchors hold for every enumerated site")


def test_cut_site_matches_historical_spcas9_formula():
    """SpCas9 cut coordinates must be unchanged from the SpCas9-only implementation."""
    for s1, e1 in ((1, 23), (1000, 1022), (99_999, 100_021)):
        assert SPCAS9.cut_site(s1, e1, "+") == e1 - 5
        assert SPCAS9.cut_site(s1, e1, "-") == s1 + 4
    # SaCas9's 6 nt PAM shifts the same geometry by three bases.
    assert SACAS9.cut_site(100, 126, "+") == 126 - 8
    assert SACAS9.cut_site(100, 126, "-") == 100 + 7
    print("ok  cut-site geometry (SpCas9 unchanged, SaCas9 PAM-shifted)")


# --------------------------------------------------------------------------------------
# enumeration
# --------------------------------------------------------------------------------------


def test_enumerate_sacas9_plus_and_minus():
    proto = "ACGTACGTACGTACGTACGTA"           # 21 nt
    seq = "TT" + proto + "CAGAGT" + "TT"       # NNGRRT PAM
    plus = [h for h in enumerate_guides_in_segment(seq, 0, len(seq), SACAS9)
            if h[0] == "+"]
    assert len(plus) == 1, plus
    _strand, s1, e1, p, pam, target = plus[0]
    assert (p, pam) == (proto, "CAGAGT")
    assert target == proto + "CAGAGT" and seq[s1 - 1:e1] == target
    assert e1 - s1 + 1 == 27

    rc = revcomp(proto + "CAGAGT")
    seq2 = "TT" + rc + "TT"
    minus = [h for h in enumerate_guides_in_segment(seq2, 0, len(seq2), SACAS9)
             if h[0] == "-"]
    assert len(minus) == 1, minus
    _strand, s1, e1, p, pam, target = minus[0]
    assert target == proto + "CAGAGT" and pam == "CAGAGT"
    assert revcomp(seq2[s1 - 1:e1]) == target
    print("ok  SaCas9 enumeration, both strands")


def test_enumeration_matches_bruteforce_every_pam():
    rng = random.Random(20260908)
    genomes = [_random_seq(rng, 4_000) for _ in range(3)]
    # A GC-rich genome, because HSV-1 is 68% GC and NNGRRT anchors behave differently.
    genomes.append("".join(rng.choice("GGCCAT") for _ in range(4_000)))
    # A genome with ambiguity, which must be rejected from every candidate window.
    noisy = list(_random_seq(rng, 4_000))
    for i in range(0, 4_000, 37):
        noisy[i] = "N"
    genomes.append("".join(noisy))

    for spec in TEST_SPECS:
        for gi, genome in enumerate(genomes):
            fast = {(h[0], h[1], h[5])
                    for h in enumerate_guides_in_segment(genome, 0, len(genome), spec)}
            brute = _naive_enumerate(genome, spec)
            assert fast == brute, (
                f"{spec.pam} genome {gi}: enumeration disagrees; "
                f"only-fast={sorted(fast - brute)[:3]} only-brute={sorted(brute - fast)[:3]}"
            )
    print(f"ok  enumeration == brute force ({len(TEST_SPECS)} PAMs x {len(genomes)} genomes)")


def test_enumeration_matches_bruteforce_on_reference():
    genome = _reference_sequence()
    if genome is None:
        print("SKIP  reference enumeration test (data/raw/reference/NC_001806.gb absent)")
        return
    window = genome[60_000:80_000]          # a 20 kb slice keeps the brute force honest
    for spec in TEST_SPECS:
        fast = {(h[0], h[1], h[5])
                for h in enumerate_guides_in_segment(window, 0, len(window), spec)}
        assert fast == _naive_enumerate(window, spec), spec.pam
    print("ok  enumeration == brute force on 20 kb of the real HSV-1 reference")


# --------------------------------------------------------------------------------------
# the scanner -- the load-bearing equivalence
# --------------------------------------------------------------------------------------


def _targets_from(genome: str, spec: Nuclease, rng: random.Random,
                  n_real: int = 120, n_decoy: int = 40) -> dict[str, str]:
    """Real sites drawn from the genome (both strands) plus never-present decoys.

    Distinct target STRINGS only: `build_lookup` is keyed on the target site, so two
    guide ids sharing one target string would collide in the lookup and the scan
    could return only one of them. That cannot happen in production -- stage 2
    collapses identical target sites into a single guide record (this is exactly how
    the two RL2 repeat copies are handled) -- so it is a fixture concern, not a
    scanner concern. Short synthetic PAM grammars make the collision likely enough
    to matter here.
    """
    sites = sorted(_naive_enumerate(genome, spec))
    targets: dict[str, str] = {}
    seen: set[str] = set()
    if sites:
        for _strand, _start, target in rng.sample(sites, len(sites)):
            if target in seen:
                continue
            seen.add(target)
            targets[f"real{len(targets)}"] = target
            if len(targets) >= n_real:
                break
    for i in range(n_decoy):
        while True:
            spacer = _random_seq(rng, spec.spacer_length)
            pam = "".join(rng.choice(IUPAC_CODES[c]) for c in spec.pam)
            cand = spacer + pam
            if (cand not in seen and cand not in genome
                    and revcomp(cand) not in genome):
                seen.add(cand)
                targets[f"decoy{i}"] = cand
                break
    return targets


def test_scan_matches_bruteforce_every_pam():
    """The anchor jump scan must equal a naive both-strand search, for every PAM."""
    rng = random.Random(20260908)
    for spec in TEST_SPECS:
        for trial in range(4):
            if trial == 3:
                genome = "".join(rng.choice("GGCCAT") for _ in range(40_000))
            else:
                genome = _random_seq(rng, 40_000)
            targets = _targets_from(genome, spec, rng)
            fwd, rev = build_lookup(targets)

            fast = scan_genome(genome, fwd, rev, spec)
            brute = _naive_presence(genome, targets)
            naive = scan_genome_naive(genome, fwd, rev, spec)
            assert fast == brute, (
                f"{spec.pam} trial {trial}: jump scan != substring search; "
                f"only-fast={sorted(fast - brute)[:5]} only-brute={sorted(brute - fast)[:5]}"
            )
            assert naive == brute, f"{spec.pam} trial {trial}: naive scanner disagrees"
            # The decoys must genuinely be absent, or the test proves nothing.
            assert not any(gid.startswith("decoy") for gid in fast), spec.pam
            assert any(gid.startswith("real") for gid in fast), spec.pam
    print(f"ok  jump scan == brute force ({len(TEST_SPECS)} PAMs x 4 randomised genomes)")


def test_scan_matches_bruteforce_on_reference():
    """Same equivalence on the real HSV-1 genome, where GC content is 68%."""
    genome = _reference_sequence()
    if genome is None:
        print("SKIP  reference scan test (data/raw/reference/NC_001806.gb absent)")
        return
    rng = random.Random(7)
    slice_ = genome[100_000:140_000]
    for spec in TEST_SPECS:
        targets = _targets_from(slice_, spec, rng)
        fwd, rev = build_lookup(targets)
        fast = scan_genome(genome, fwd, rev, spec)
        brute = _naive_presence(genome, targets)
        assert fast == brute, (
            f"{spec.pam}: jump scan != substring search on the reference; "
            f"only-fast={sorted(fast - brute)[:5]} only-brute={sorted(brute - fast)[:5]}"
        )
    print("ok  jump scan == brute force on the real HSV-1 reference (all PAMs)")


def test_scan_edge_cases():
    """Sites at the very start and end of a sequence, and sequences shorter than one
    site, are the classic off-by-one traps in an anchored scan."""
    for spec in (SPCAS9, SACAS9, get_nuclease("sacas9", pam="NNGRRN")):
        spacer = ("ACGT" * 8)[:spec.spacer_length]
        pam = "AGG" if spec.pam == "NGG" else "CAGAGT"
        target = spacer + pam
        fwd, rev = build_lookup({"x": target})
        assert scan_genome(target, fwd, rev, spec) == {"x"}                 # exact fit
        assert scan_genome("A" + target, fwd, rev, spec) == {"x"}           # at the end
        assert scan_genome(target + "A", fwd, rev, spec) == {"x"}           # at the start
        assert scan_genome(revcomp(target), fwd, rev, spec) == {"x"}        # minus strand
        assert scan_genome("A" + revcomp(target) + "A", fwd, rev, spec) == {"x"}
        assert scan_genome(target[:-1], fwd, rev, spec) == set()            # truncated
        assert scan_genome("", fwd, rev, spec) == set()
        # A single mismatch anywhere -- protospacer or PAM -- must lose the site.
        for i in range(len(target)):
            mutated = target[:i] + ("A" if target[i] != "A" else "C") + target[i + 1:]
            assert scan_genome(mutated, fwd, rev, spec) == set(), (spec.pam, i)
    print("ok  scanner edge cases (boundaries, truncation, single mismatches)")


def main() -> int:
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in tests:
        fn()
    print(f"\nAll {len(tests)} tests passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
