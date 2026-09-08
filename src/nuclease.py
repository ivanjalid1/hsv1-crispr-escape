"""Nuclease / PAM model -- makes the pipeline nuclease-agnostic instead of NGG-only.

Why this module exists
----------------------
Stages 2 and 3 originally hardcoded the SpCas9 `NGG` PAM: enumeration tested
`window[21] == "G" and window[22] == "G"`, and the conservation scanner jumped
between literal `GG`/`CC` occurrences. That is correct for SpCas9 and wrong for
everything else, which matters because the published competitor guide set we
benchmark against (Amrani et al. 2024, doi:10.1016/j.omtm.2024.101303) is
**SaCas9 / NNGRRT** and therefore occupies a completely different site space.

A `Nuclease` here is a small, immutable description of a target site:

    protospacer (spacer_length nt)  +  PAM (IUPAC pattern, immediately 3')

and everything else -- footprint length, PAM matching on both strands, the
literal anchors the fast scanner may jump between, the predicted cut site -- is
derived from it. No PAM literal appears anywhere else in the codebase.

IUPAC
-----
PAM patterns are IUPAC nucleotide codes, so `NNGRRT` means
"any, any, G, {A,G}, {A,G}, T". Genomic windows are only ever matched when they
are unambiguous ACGT, so ambiguity lives in the pattern, never in the sequence.

Scanning anchors
----------------
For a 3'-PAM nuclease every concrete target site shares the *literal* positions
of its PAM pattern: every SpCas9 site really does carry `GG` at target offset
`spacer_length + 1`, and every SaCas9 NNGRRT site really does carry `G` at
`spacer_length + 2` and `T` at `spacer_length + 5`. `scan_anchors` exposes those
maximal literal runs so `src/conservation.py` can jump between their occurrences
with `str.find` instead of testing every position. The optimisation is exact by
construction -- an occurrence of the target implies an occurrence of the anchor at
the fixed offset -- and `tests/test_nuclease.py` proves set-identity against a
naive both-strand search for every supported PAM, on randomised sequences and on
the real HSV-1 reference.

A PAM with no literal position at all (e.g. `NNN`) yields no anchors, and the
scanner falls back to testing every position. Correct, just slower.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass, field
from functools import cached_property
from pathlib import Path

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.common import revcomp  # noqa: E402

# IUPAC nucleotide ambiguity codes -> the concrete DNA bases they admit.
IUPAC_CODES: dict[str, str] = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG",    # puRine
    "Y": "CT",    # pYrimidine
    "S": "GC",
    "W": "AT",
    "K": "GT",    # Keto
    "M": "AC",    # aMino
    "B": "CGT",
    "D": "AGT",
    "H": "ACT",
    "V": "ACG",
    "N": "ACGT",
}

#: Codes matching exactly one base -- the positions a scanner may anchor on.
LITERAL_CODES = frozenset("ACGT")


def iupac_complement(code: str) -> str:
    """Complement of an IUPAC code (the code matching the complements of its bases)."""
    bases = IUPAC_CODES[code]
    comp = {revcomp(b) for b in bases}
    for candidate, admitted in IUPAC_CODES.items():
        if set(admitted) == comp:
            return candidate
    raise ValueError(f"No IUPAC code complements {code!r}")  # unreachable


def iupac_revcomp(pattern: str) -> str:
    """Reverse complement of an IUPAC pattern, e.g. NGG -> CCN, NNGRRT -> AYYCNN."""
    return "".join(iupac_complement(c) for c in reversed(pattern))


@dataclass(frozen=True)
class Nuclease:
    """An RNA-guided nuclease's target-site grammar.

    Attributes
    ----------
    name            short identifier, also the output namespace tag stem
    spacer_length   protospacer length in nt
    pam             IUPAC PAM pattern, 5'->3' on the protospacer strand
    pam_side        "3prime" only. 5'-PAM nucleases (Cas12a) are rejected loudly
                    rather than half-supported: their cut geometry is staggered
                    and distal to the PAM, which this project does not model.
    description     free text for reports
    """

    name: str
    spacer_length: int
    pam: str
    pam_side: str = "3prime"
    description: str = ""
    _explicit_pam: bool = field(default=False, repr=False, compare=False)
    _explicit_spacer: bool = field(default=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        if self.spacer_length < 1:
            raise ValueError(f"spacer_length must be >= 1, got {self.spacer_length}")
        pam = self.pam.upper()
        if not pam:
            raise ValueError("PAM pattern must be non-empty")
        bad = sorted(set(pam) - set(IUPAC_CODES))
        if bad:
            raise ValueError(
                f"PAM {self.pam!r} contains non-IUPAC symbol(s): {', '.join(bad)}. "
                f"Allowed: {''.join(sorted(IUPAC_CODES))}"
            )
        object.__setattr__(self, "pam", pam)
        if self.pam_side != "3prime":
            raise NotImplementedError(
                f"pam_side={self.pam_side!r} is not supported. This pipeline models "
                "3'-PAM nucleases (SpCas9, SaCas9) only; a 5'-PAM nuclease such as "
                "Cas12a also needs a different (staggered, PAM-distal) cut model, "
                "which is deliberately not faked here."
            )

    # -- geometry ---------------------------------------------------------------

    @property
    def pam_length(self) -> int:
        return len(self.pam)

    @property
    def target_length(self) -> int:
        """Length of the footprint whose exact presence defines conservation."""
        return self.spacer_length + self.pam_length

    @property
    def tag(self) -> str:
        """Output-namespace tag, e.g. 'spcas9', 'sacas9', 'sacas9-nngrrn-20nt'."""
        parts = [self.name]
        if self._explicit_pam:
            parts.append(self.pam.lower())
        if self._explicit_spacer:
            parts.append(f"{self.spacer_length}nt")
        return "-".join(parts)

    @property
    def label(self) -> str:
        return f"{self.name} ({self.spacer_length} nt spacer + {self.pam} PAM)"

    # -- PAM matching -----------------------------------------------------------

    @cached_property
    def _plus_sets(self) -> tuple[frozenset[str], ...]:
        return tuple(frozenset(IUPAC_CODES[c]) for c in self.pam)

    @cached_property
    def _minus_sets(self) -> tuple[frozenset[str], ...]:
        """Per-position base sets for the PAM as it appears on the PLUS strand when
        the protospacer is on the minus strand (i.e. the reverse complement pattern,
        which occupies the first `pam_length` bases of the plus-strand window)."""
        return tuple(frozenset(IUPAC_CODES[c]) for c in iupac_revcomp(self.pam))

    def pam_matches(self, seq: str) -> bool:
        """True if `seq` (concrete ACGT, length == pam_length) satisfies the pattern."""
        sets = self._plus_sets
        if len(seq) != len(sets):
            return False
        return all(ch in s for ch, s in zip(seq, sets))

    def window_has_plus_pam(self, window: str) -> bool:
        """`window` is a plus-strand target-length slice; is there a plus-strand PAM?"""
        sets = self._plus_sets
        off = self.spacer_length
        for i, s in enumerate(sets):
            if window[off + i] not in s:
                return False
        return True

    def window_has_minus_pam(self, window: str) -> bool:
        """Same window, protospacer on the minus strand: the PAM is the reverse
        complement pattern sitting at the LOW plus-strand coordinates."""
        for i, s in enumerate(self._minus_sets):
            if window[i] not in s:
                return False
        return True

    # -- scanner anchors --------------------------------------------------------

    @cached_property
    def scan_anchors(self) -> tuple[tuple[str, int], ...]:
        """Maximal literal runs of the PAM as (literal, offset within the target).

        Every concrete target site contains each of these literals at exactly the
        given offset, so a scanner may jump between their occurrences. Ordered
        longest-first, then leftmost, so a caller that just takes the first gets the
        most selective one; `src/conservation.py` refines the choice per sequence.
        """
        runs: list[tuple[str, int]] = []
        i = 0
        pam = self.pam
        while i < len(pam):
            if pam[i] in LITERAL_CODES:
                j = i
                while j < len(pam) and pam[j] in LITERAL_CODES:
                    j += 1
                runs.append((pam[i:j], self.spacer_length + i))
                i = j
            else:
                i += 1
        runs.sort(key=lambda r: (-len(r[0]), r[1]))
        return tuple(runs)

    # -- cut site ---------------------------------------------------------------

    def cut_site(self, ref_start: int, ref_end: int, strand: str) -> int:
        """Plus-strand coordinate reported as the predicted blunt cut position.

        Convention (unchanged from the original SpCas9-only implementation, so that
        SpCas9 output stays byte-identical): the cut is reported near the PAM-proximal
        end of the protospacer, `pam_length + 2` bases inside the footprint from the
        PAM edge. `ref_start`/`ref_end` are 1-based inclusive plus-strand coordinates
        of the whole footprint.
        """
        if strand == "+":
            return ref_end - self.pam_length - 2
        return ref_start + self.pam_length + 1

    def as_dict(self) -> dict:
        return {
            "nuclease": self.name,
            "spacer_length": self.spacer_length,
            "pam": self.pam,
            "pam_side": self.pam_side,
            "target_length": self.target_length,
            "tag": self.tag,
            "description": self.description,
        }


# --------------------------------------------------------------------------------------
# Registry
# --------------------------------------------------------------------------------------

SPCAS9 = Nuclease(
    name="spcas9",
    spacer_length=20,
    pam="NGG",
    description="Streptococcus pyogenes Cas9; canonical 20 nt spacer, 3' NGG PAM.",
)

SACAS9 = Nuclease(
    name="sacas9",
    spacer_length=21,
    pam="NNGRRT",
    description=(
        "Staphylococcus aureus Cas9; 21 nt spacer, 3' NNGRRT PAM. Used by "
        "Amrani et al. 2024 (EBT-104), whose Table 1 reports 20 nt DNA spacers -- "
        "pass --spacer-length 20 to reproduce that exact site space."
    ),
)

NUCLEASES: dict[str, Nuclease] = {n.name: n for n in (SPCAS9, SACAS9)}

#: PAM variants that may be selected with --pam for a given nuclease. The values are
#: not enforced (any IUPAC string is accepted) -- this is documentation for --help.
KNOWN_PAM_VARIANTS = {
    "spcas9": ["NGG", "NAG", "NGA"],
    "sacas9": ["NNGRRT", "NNGRRN"],
}

DEFAULT_NUCLEASE = "spcas9"


def get_nuclease(name: str = DEFAULT_NUCLEASE, pam: str | None = None,
                 spacer_length: int | None = None) -> Nuclease:
    """Resolve a nuclease from CLI arguments, applying PAM/spacer overrides.

    Overrides are recorded so that `Nuclease.tag` -- and therefore the output
    namespace -- distinguishes e.g. sacas9/NNGRRT from sacas9/NNGRRN.
    """
    key = (name or DEFAULT_NUCLEASE).strip().lower()
    if key not in NUCLEASES:
        raise SystemExit(
            f"Unknown nuclease {name!r}. Available: {', '.join(sorted(NUCLEASES))}"
        )
    base = NUCLEASES[key]
    new_pam = (pam or base.pam).strip().upper()
    new_spacer = int(spacer_length) if spacer_length is not None else base.spacer_length
    return Nuclease(
        name=base.name,
        spacer_length=new_spacer,
        pam=new_pam,
        pam_side=base.pam_side,
        description=base.description,
        _explicit_pam=new_pam != base.pam,
        _explicit_spacer=new_spacer != base.spacer_length,
    )
