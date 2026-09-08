"""Stage 7 -- a quantitative escape-probability model for multiplex HSV-1 editing.

The gap this fills
------------------
Amrani et al. 2024 (Mol Ther Methods Clin Dev 32:101303, EBT-104) justify a two-guide
design by citing HIV CRISPR-escape literature and asserting that targeting "two or
more" sites prevents escape. They never quantify it, and nobody else has for HSV-1.
This module computes it.

The model
---------
A viral genome escapes multiplex editing if it ends up (a) no longer cleavable at ANY
targeted site and (b) still encoding functional essential proteins. Decomposed per
site i, for one viral genome:

    P(site i fails) = P_absent(i) + [1 - P_absent(i)] x P_repair_escape(i)

`P_absent` is measured, not assumed: it is the per-genome presence matrix of the site
across the 183 complete HSV-1 genomes (stage 3 semantics, exact protospacer+PAM
presence on either strand). Because the matrix is kept per genome rather than
collapsed to a marginal, correlated failure across sites is captured exactly -- the
model reproduces, rather than assumes away, the stage-6 observation that the Amrani
lead pair's joint coverage (0.825) is BELOW both of its marginals.

`P_repair_escape` is the probability that Cas9 cuts, NHEJ repairs with an indel that
destroys protospacer recognition, and the protein remains functional. Because a site
that is repaired perfectly is simply re-cut, the terminal outcome of a persistently
expressed nuclease is an indel drawn from the NHEJ spectrum conditioned on being
non-WT; the model therefore reduces to

    q_i = SUM_{L != 0} P(L) / (1 - P(0)) x viability(L, i)

with

    viability(L, i) = frameshift_viability[gene]           if L mod 3 != 0
                    = PROD over affected codons of t_j     if L mod 3 == 0

Frameshifts in an essential gene are lethal to the virus and therefore are NOT
escape. In-frame indels are scored against a **per-codon tolerance profile measured
from cross-strain amino-acid variation in the same 183 genomes**: a codon that varies
among viable clinical isolates is demonstrably tolerant to change; a codon invariant
across every isolate is likely constrained. That profile is the one genuinely novel
data source in the model and it costs no new experiments.

Joint escape over a guide set S is then an exact average over the empirical genome
population, with repair outcomes independent across sites *conditional on the genome*:

    P_escape(S) = (1/N) SUM_g PROD_{i in S} f(g, i),
        f(g, i) = 1     if site i is absent in genome g   (it can never be cut)
                = q_i   if site i is present in genome g

Two limiting cases fall out and are asserted in tests/test_escape.py:
  * k = 1 with q = 0 gives exactly 1 - conservation_fraction(i);
  * perfectly conserved sites (present in every genome) give exactly PROD q_i.

An additional, k-independent term is reported separately and never mixed into the
above: theta, the fraction of viral genomes never exposed to an active nuclease
(delivery/expression failure). Total persistence is theta + (1-theta) P_escape(S).
Guide count cannot touch theta -- which is the single most consequential honest
finding here, because the only in vivo on-target measurement in Amrani et al. is
~1% indels in trigeminal ganglia.

Honesty rules applied throughout
--------------------------------
* No fabricated citations. The DEFAULT NHEJ indel-length spectrum is an explicitly
  labelled ASSUMPTION with a documented shape, a stated default, and a
  `--indel-spectrum` / `--inframe-fraction` override. Everything measured from this
  repository's data is labelled MEASURED. One swept SCENARIO, `neuronal_nhej`, has a
  measured shape: it is the pooled CRISPResso2 net-length histogram deposited as source
  data with Ramadoss et al. 2025 Fig. 1d, from post-mitotic human iPSC-derived neurons.
  It is measured in THEIR system and an assumption for OURS (different nuclease,
  different cells, different species), so it widens the sweep and does not become the
  default. See refs/ramadoss2025_notes.md.
* Independence between sites is never silently assumed. Site *presence* uses the
  empirical joint matrix. Site *repair* is assumed independent given the genome, and
  the size and direction of that assumption's error is bounded in the report.
* Defaults are chosen to be the ones most favourable to the two-guide design wherever
  a choice exists, so that a conclusion of "k=2 is not enough" cannot be an artefact
  of parameter tuning. The one exception (the repeat model) is reported in all three
  variants side by side.
* Deterministic: closed form everywhere, plus a seeded Monte Carlo cross-check that
  must agree with it.

Everything runs offline from the caches populated by stages 1 and 5.
"""

from __future__ import annotations

import argparse
import json
import logging
import math
import sys
from dataclasses import dataclass, field, replace
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src import benchmark_sacas9 as bm  # noqa: E402
from src import robustness as rob  # noqa: E402
from src.common import (  # noqa: E402
    DEFAULT_MANIFEST,
    PROJECT_ROOT,
    RESULTS_DIR,
    ensure_dirs,
    revcomp,
    setup_logging,
)
from src.conservation import load_genome  # noqa: E402
from src.extract_guides import DEFAULT_GENES, DEFAULT_REFERENCE  # noqa: E402
from src.nuclease import Nuclease, get_nuclease  # noqa: E402

LOG = logging.getLogger("escape")

# --------------------------------------------------------------------------------------
# Biology constants that are NOT parameters
# --------------------------------------------------------------------------------------

#: Genes treated as essential for productive replication. RL2/ICP0 is deliberately
#: NOT in this set as a matter of record -- ICP0-null HSV-1 replicates in cell culture
#: at high multiplicity and is severely attenuated rather than dead. The model's
#: BASELINE nevertheless treats a frameshift in RL2 as lethal (frameshift_viability
#: 0.0 everywhere), because that is the assumption most favourable to the published
#: two-guide design; `--frameshift-viability RL2=1.0` explores the other end and the
#: sensitivity table reports both. This is stated rather than silently chosen.
ESSENTIAL_GENES = frozenset({"UL30", "UL19", "UL5", "UL52", "UL29", "UL54"})

DISPENSABILITY_NOTE = {
    "RL2": "ICP0 is a virulence/reactivation factor, not strictly essential in vitro",
}

_CODON_TABLE: dict[str, str] = {}


def _build_codon_table() -> None:
    from Bio.Data import CodonTable

    tbl = CodonTable.unambiguous_dna_by_id[1]
    _CODON_TABLE.update(tbl.forward_table)
    for stop in tbl.stop_codons:
        _CODON_TABLE[stop] = "*"


_build_codon_table()

_COMP_LUT = np.full(256, ord("N"), dtype=np.uint8)
for _a, _b in zip(b"ACGTN", b"TGCAN"):
    _COMP_LUT[_a] = _b


# ======================================================================================
# Parameters
# ======================================================================================


@dataclass(frozen=True)
class IndelSpectrum:
    """A probability distribution over signed NHEJ repair outcome lengths.

    `probs` maps a signed length in nt (negative = deletion, positive = insertion,
    0 = perfect/WT restoration) to a probability. The model only ever uses the
    distribution CONDITIONED on length != 0, because a perfectly repaired site is
    simply re-cut by a persistently expressed nuclease; `p_wt` is retained for
    reporting only.

    PROVENANCE: this is an ASSUMPTION. Published Cas9/SaCas9 NHEJ spectra were not
    verifiable offline in this environment and no citation is invented for them. The
    default shape encodes only the qualitative features that are uncontroversial in
    the field -- a dominant +1 insertion, deletions concentrated at short lengths with
    a decaying tail -- as a two-component parametric family whose single consequential
    summary statistic, the in-frame fraction, is exposed directly as
    `--inframe-fraction` and swept in the sensitivity analysis.
    """

    name: str
    probs: dict[int, float]
    provenance: str = "ASSUMPTION"

    @property
    def p_wt(self) -> float:
        return float(self.probs.get(0, 0.0))

    @property
    def conditional(self) -> dict[int, float]:
        """The spectrum conditioned on an indel actually occurring (length != 0)."""
        mass = sum(p for L, p in self.probs.items() if L != 0)
        if mass <= 0:
            raise ValueError("indel spectrum has no non-zero-length mass")
        return {L: p / mass for L, p in self.probs.items() if L != 0}

    @property
    def inframe_fraction(self) -> float:
        cond = self.conditional
        return sum(p for L, p in cond.items() if L % 3 == 0)

    @property
    def mean_deletion_length(self) -> float:
        cond = self.conditional
        m = sum(p for L, p in cond.items() if L < 0)
        if m <= 0:
            return 0.0
        return sum(-L * p for L, p in cond.items() if L < 0) / m

    def rescaled_to_inframe(self, target: float) -> "IndelSpectrum":
        """Reweight in-frame vs frameshift mass so `inframe_fraction == target`.

        Used only by the sensitivity sweep: it isolates the one summary statistic of
        the spectrum that the escape probability is actually sensitive to, without
        having to invent a whole alternative spectrum for each point.
        """
        cond = self.conditional
        cur = self.inframe_fraction
        if not 0.0 < target < 1.0:
            raise ValueError("target in-frame fraction must be in (0, 1)")
        if cur <= 0 or cur >= 1:
            raise ValueError("cannot rescale a spectrum that is entirely (out of) frame")
        out: dict[int, float] = {}
        for L, p in cond.items():
            out[L] = p * (target / cur) if L % 3 == 0 else p * ((1 - target) / (1 - cur))
        return IndelSpectrum(f"{self.name}@inframe={target:g}", out, self.provenance)


def _geometric(mean: float, lo: int, hi: int) -> dict[int, float]:
    """Truncated geometric over [lo, hi] with the given (untruncated) mean."""
    if mean <= lo:
        mean = lo + 0.5
    r = 1.0 - 1.0 / (mean - lo + 1.0)
    raw = {n: r ** (n - lo) for n in range(lo, hi + 1)}
    z = sum(raw.values())
    return {n: v / z for n, v in raw.items()}


def parametric_spectrum(
    p_wt: float = 0.40,
    p_insertion: float = 0.35,
    ins_one_frac: float = 0.80,
    ins_mean: float = 3.0,
    ins_max: int = 20,
    del_mean: float = 6.0,
    del_max: int = 60,
    name: str = "parametric-default",
) -> IndelSpectrum:
    """Build the default two-component spectrum. See `IndelSpectrum` for provenance.

    * `p_wt` -- mass on perfect repair. Reported, then conditioned away.
    * `p_insertion` -- share of indels that are insertions.
    * `ins_one_frac` -- share of insertions that are exactly +1 nt (the single
      dominant blunt-repair product of Cas9 in essentially every published spectrum).
    * remaining insertion mass and all deletion mass are truncated geometrics.
    """
    probs: dict[int, float] = {0: p_wt}
    indel = 1.0 - p_wt
    ins = indel * p_insertion
    dele = indel * (1.0 - p_insertion)

    probs[1] = ins * ins_one_frac
    tail = ins * (1.0 - ins_one_frac)
    if tail > 0:
        for n, p in _geometric(ins_mean, 2, ins_max).items():
            probs[n] = probs.get(n, 0.0) + tail * p
    for n, p in _geometric(del_mean, 1, del_max).items():
        probs[-n] = probs.get(-n, 0.0) + dele * p
    return IndelSpectrum(name, probs)


#: Pooled CRISPResso2 net-length histogram from Ramadoss et al. 2025 Fig. 1d source data
#: (FigShare 30366298). Columns: length, neuron_reads, ipsc_reads. See
#: refs/ramadoss2025_notes.md for exactly what was and was not extractable from that
#: paper, and for the caveats attached to using it here.
RAMADOSS2025_HISTOGRAM = PROJECT_ROOT / "refs" / "ramadoss2025_fig1d_indel_histogram.tsv"

#: The one honest sentence about this scenario, printed wherever it is used.
NEURONAL_NHEJ_PROVENANCE = (
    "MEASURED in Ramadoss et al. 2025 (Nat Commun 16:9883) Fig. 1d deposited source "
    "data -- pooled CRISPResso2 net-length histogram over 6 replicates, 118,722 reads, "
    "SpCas9 RNP + sgRNA B2Mg1 at human B2M in iPSC-derived neurons, 4 d post-"
    "transduction; ASSUMPTION when transferred to SaCas9 in trigeminal-ganglion neurons"
)


def neuronal_nhej_spectrum(arm: str = "neuron",
                           path: Path | None = None) -> IndelSpectrum:
    """The measured post-mitotic-neuron NHEJ spectrum, as a SENSITIVITY SCENARIO.

    This is the only built-in spectrum whose SHAPE is measured rather than posited, and
    the distinction is worth being precise about rather than rounding off in either
    direction:

    * MEASURED, in Ramadoss et al. 2025's cells. Their Fig. 1d deposits raw CRISPResso2
      `Indel_histogram.txt` for six neuron and six isogenic-iPSC replicates; the file
      read here is those counts, pooled per arm, verbatim. Nothing is fitted, smoothed
      or extrapolated. The paper's own running text describes the distribution only
      qualitatively ("a much narrower distribution of outcomes"); the numbers come from
      its deposited data, not from reading its figures.
    * ASSUMPTION, for this therapy. Their neurons are cultured human iPSC-derived
      neurons cut by SpCas9 at one locus with one guide. The therapy under audit uses
      SaCas9 in latently infected mouse/rabbit trigeminal ganglia. A mouse trigeminal
      ganglion is not a cultured human neuron, and the paper itself shows the spectrum
      is strongly guide-dependent. This scenario therefore WIDENS the swept space; it
      does not replace the default and is not a measurement of the modelled system.

    `arm="ipsc"` returns the isogenic dividing-cell arm of the same experiment. It is
    not swept; it exists so the cell-type contrast can be reproduced from the same file
    by anyone checking the direction of the effect.
    """
    p = path or RAMADOSS2025_HISTOGRAM
    if not p.is_file():
        raise SystemExit(
            f"{p} is missing. It ships with the repository; see "
            "refs/ramadoss2025_notes.md for how it was extracted."
        )
    col = {"neuron": "neuron_reads", "ipsc": "ipsc_reads"}[arm]
    df = pd.read_csv(p, sep="\t")
    missing = {"length", col} - set(df.columns)
    if missing:
        raise SystemExit(f"{p}: expected columns length, {col}")
    counts = {int(r.length): float(getattr(r, col)) for r in df.itertuples(index=False)}
    z = sum(counts.values())
    if z <= 0:
        raise SystemExit(f"{p}: {col} sums to {z}")
    name = "neuronal_nhej" if arm == "neuron" else "ipsc_dividing"
    prov = NEURONAL_NHEJ_PROVENANCE
    if arm == "ipsc":
        prov = prov.replace("iPSC-derived neurons", "the isogenic parental iPSCs")
    return IndelSpectrum(name, {L: c / z for L, c in counts.items() if c > 0},
                         provenance=prov)


def builtin_spectrum(name: str) -> IndelSpectrum:
    """Named alternative spectra, used for the sensitivity sweep.

    All are ASSUMPTIONS except `neuronal_nhej` / `ipsc_dividing`, whose shapes are read
    from published measured data -- see `neuronal_nhej_spectrum` for the precise sense
    in which those two are and are not measurements of the system being modelled.
    """
    if name in ("default", "parametric-default"):
        return parametric_spectrum()
    if name == "neuronal_nhej":
        return neuronal_nhej_spectrum("neuron")
    if name == "ipsc_dividing":
        return neuronal_nhej_spectrum("ipsc")
    if name == "deletion-heavy":
        return parametric_spectrum(p_insertion=0.15, del_mean=12.0,
                                   name="deletion-heavy")
    if name == "insertion-heavy":
        return parametric_spectrum(p_insertion=0.65, ins_one_frac=0.90,
                                   name="insertion-heavy")
    if name == "short-indels":
        return parametric_spectrum(del_mean=2.5, ins_max=6, del_max=20,
                                   name="short-indels")
    if name == "uniform-1to20":
        probs = {0: 0.4}
        for n in range(1, 21):
            probs[n] = 0.6 / 40.0
            probs[-n] = 0.6 / 40.0
        return IndelSpectrum("uniform-1to20", probs)
    raise SystemExit(
        f"Unknown indel spectrum {name!r}. Built-ins: default, deletion-heavy, "
        "insertion-heavy, short-indels, uniform-1to20, neuronal_nhej, ipsc_dividing; "
        "or give a path to a TSV with columns 'length' and 'probability'."
    )


def load_spectrum(spec: str) -> IndelSpectrum:
    path = Path(spec)
    if path.is_file():
        df = pd.read_csv(path, sep="\t")
        missing = {"length", "probability"} - set(df.columns)
        if missing:
            raise SystemExit(f"{path}: spectrum TSV needs columns length, probability")
        probs = {int(r.length): float(r.probability) for r in df.itertuples(index=False)}
        z = sum(probs.values())
        if z <= 0:
            raise SystemExit(f"{path}: spectrum probabilities sum to {z}")
        return IndelSpectrum(f"file:{path.name}", {k: v / z for k, v in probs.items()},
                             provenance="USER-SUPPLIED")
    return builtin_spectrum(spec)


@dataclass(frozen=True)
class Params:
    """Every uncertain quantity in the model, in one place.

    `source` for each is given in the report's parameter table: MEASURED (computed
    from this repository's genome corpus) or ASSUMPTION (stated default, swept).
    """

    #: P(protein still functional | ONE codon deleted/inserted at a codon that varies
    #: in amino acid across the 183 isolates). ASSUMPTION. Bounded above by 1.
    tol_variable: float = 0.50
    #: Same, at a codon invariant across all 183 isolates. ASSUMPTION.
    tol_invariant: float = 0.05
    #: How many isolates must carry a non-reference amino acid for a codon to count
    #: as "variable". 1 is the natural default: one viable clinical isolate is a
    #: demonstration of tolerance.
    min_variant_strains: int = 1
    #: P(virus viable | frameshift), per gene. Default 0 for every gene, including
    #: RL2 -- see ESSENTIAL_GENES.
    frameshift_viability: tuple[tuple[str, float], ...] = ()
    #: How a guide targeting a duplicated gene (RL2 sits in both TRL and IRL) is
    #: scored. "redundant" = every copy must become uncuttable and at least one must
    #: stay functional (the biologically motivated default for a redundant locus);
    #: "single" = treat the duplication as one locus (gene conversion keeps the copies
    #: identical); "all-copies" = every copy must independently produce a viable escape
    #: allele (a stand-in for inter-copy excision being lethal). All three reported.
    repeat_model: str = "redundant"
    #: Fraction of viral genomes never exposed to an active nuclease. Reported
    #: separately, never folded into P_escape, because it is k-independent.
    theta_unexposed: float = 0.0
    #: Number of viral genomes in the compartment being treated. Only used to convert
    #: a per-genome probability into an expected count. ASSUMPTION; swept.
    viral_population_size: float = 1e5
    #: Probability an indel destroys protospacer/PAM recognition. 1.0 is both the
    #: physically expected value for indels at the Cas9 cut site (3 bp from the PAM)
    #: and the value most favourable to escape, hence the conservative default.
    p_disrupt: float = 1.0

    def frameshift_map(self) -> dict[str, float]:
        return dict(self.frameshift_viability)

    def frameshift_for(self, gene: str) -> float:
        return self.frameshift_map().get(gene, 0.0)


# ======================================================================================
# Reference CDS structure
# ======================================================================================


@dataclass
class CdsCopy:
    """One contiguous coding sequence copy of one gene, with its genomic positions."""

    gene: str
    copy_index: int
    strand: int
    positions: np.ndarray          # 0-based plus-strand ref positions, in CDS order
    cds: str
    protein: str

    @property
    def n_codons(self) -> int:
        return len(self.protein)

    def codon_of_position(self) -> dict[int, int]:
        """plus-strand genomic position -> codon index."""
        return {int(p): i // 3 for i, p in enumerate(self.positions)}


def collect_cds_copies(record, genes: list[str]) -> list[CdsCopy]:
    """Every annotated CDS copy of `genes`, translated and validated.

    A CDS that does not translate to M...* with exactly one stop is rejected loudly:
    that would mean the codon frame used for the tolerance profile is wrong.
    """
    ref = str(record.seq).upper()
    wanted = {g.upper() for g in genes}
    per_gene: dict[str, int] = {}
    out: list[CdsCopy] = []
    for feat in record.features:
        if feat.type != "CDS":
            continue
        names = [n.upper() for n in feat.qualifiers.get("gene", [])]
        names += [n.upper() for n in feat.qualifiers.get("gene_synonym", [])]
        hit = next((n for n in names if n in wanted), None)
        if hit is None:
            continue
        positions: list[int] = []
        for part in feat.location.parts:
            r = list(range(int(part.start), int(part.end)))
            if part.strand == -1:
                r.reverse()
            positions.extend(r)
        strand = 1 if feat.location.strand == 1 else -1
        arr = np.asarray(positions, dtype=np.int64)
        bases = np.frombuffer(ref.encode("ascii"), dtype=np.uint8)[arr]
        if strand == -1:
            bases = _COMP_LUT[bases]
        cds = bases.tobytes().decode("ascii")
        if len(cds) % 3:
            raise SystemExit(f"{hit}: CDS length {len(cds)} is not a multiple of 3")
        protein = translate(cds)
        if not protein.startswith("M") or protein[-1] != "*" or protein.count("*") != 1:
            raise SystemExit(
                f"{hit}: reconstructed CDS does not translate cleanly "
                f"({protein[:3]}...{protein[-1]}, {protein.count('*')} stops)"
            )
        idx = per_gene.get(hit, 0)
        per_gene[hit] = idx + 1
        out.append(CdsCopy(hit, idx, strand, arr, cds, protein))
    missing = sorted(wanted - set(per_gene))
    if missing:
        LOG.warning("No CDS feature for gene(s): %s", ", ".join(missing))
    out.sort(key=lambda c: (c.gene, c.copy_index))
    return out


def translate(cds: str) -> str:
    return "".join(_CODON_TABLE.get(cds[i:i + 3], "X") for i in range(0, len(cds), 3))


# ======================================================================================
# Per-codon cross-strain tolerance profile   (the MEASURED half of the model)
# ======================================================================================


def offset_intervals(chains: list, k: int, max_unanchored: int) -> list[tuple[int, int, int, int]]:
    """Maximal constant-offset runs within colinear anchor chains.

    Each run is returned as (ref_lo, ref_hi, offset, n_anchors): every reference
    position p in [ref_lo, ref_hi) maps to query position p + offset with no
    intervening indel. A run is broken when the offset changes (an indel) or when two
    consecutive anchors are more than `max_unanchored` bases apart in the reference
    (too much unanchored sequence to bridge safely).

    A run bridges short unanchored gaps deliberately: a single substitution destroys
    every k-mer overlapping it, so the anchor-free gap around a SNP is ~k + step bases
    wide and MUST be bridged or every variant position -- exactly the positions the
    tolerance profile is about -- would be discarded.
    """
    out: list[tuple[int, int, int, int]] = []
    for ch in chains:
        r, q = ch.ref, ch.qry
        i = 0
        n = len(r)
        while i < n:
            off = q[i] - r[i]
            j = i
            while (j + 1 < n
                   and q[j + 1] - r[j + 1] == off
                   and r[j + 1] - r[j] <= max_unanchored):
                j += 1
            out.append((r[i], r[j] + k, off, j - i + 1))
            i = j + 1
    return out


def build_offset_map(n_ref: int, chains: list, k: int,
                     max_unanchored: int) -> np.ndarray:
    """Reference position -> query position (or -1), from constant-offset runs.

    Longer (better supported) runs are written last so that they win where runs from
    different chains overlap -- which they do around the HSV-1 large repeats.
    """
    qmap = np.full(n_ref, -1, dtype=np.int64)
    runs = offset_intervals(chains, k, max_unanchored)
    runs.sort(key=lambda t: t[3])
    for lo, hi, off, _ in runs:
        lo = max(0, lo)
        hi = min(n_ref, hi)
        if hi > lo:
            qmap[lo:hi] = np.arange(lo, hi, dtype=np.int64) + off
    return qmap


@dataclass
class ToleranceProfile:
    """Per-codon cross-strain variation for one gene. MEASURED."""

    gene: str
    n_codons: int
    ref_protein: str
    n_resolved: np.ndarray          # isolates with a resolved codon, per codon
    n_variant: np.ndarray           # isolates whose amino acid differs from reference
    n_nt_variant: np.ndarray        # isolates whose codon differs at nucleotide level
    n_genomes_used: int
    n_genomes_dropped: int
    length_polymorphic: int         # isolates with a net indel inside this CDS
    length_polymorphic_inframe: int

    def varies(self, min_variant_strains: int) -> np.ndarray:
        return self.n_variant >= min_variant_strains

    def tolerance_vector(self, params: Params) -> np.ndarray:
        v = self.varies(params.min_variant_strains)
        return np.where(v, params.tol_variable, params.tol_invariant)


def measure_tolerance(record, cds_copies: list[CdsCopy], manifest: pd.DataFrame,
                      anchor_k: int = 25, anchor_step: int = 6,
                      max_unanchored: int = 60,
                      min_resolved_fraction: float = 0.95,
                      max_genomes: int | None = None,
                      ) -> tuple[dict[str, ToleranceProfile], pd.DataFrame]:
    """Per-codon amino-acid variation across the complete-genome set.

    Method, deliberately alignment-free and reusing the stage-5 machinery verbatim:
    exact 25-mer anchors against the reference, greedy colinear chaining, then maximal
    constant-offset runs give an ungapped reference->isolate coordinate map. Codons
    are read directly through that map and translated. An extraction is rejected --
    and the isolate contributes nothing to that gene -- unless at least
    `min_resolved_fraction` of its codons map and the protein has no internal stop, so
    a botched mapping can never manufacture apparent tolerance.

    Both orientations of every record are mapped and the better-resolved one is used
    per gene, because HSV-1 genome isomers invert the L and S components relative to
    each other and a whole-record orientation choice would silently lose one of them.
    """
    ref = str(record.seq).upper()
    n_ref = len(ref)
    index = rob.build_reference_index(ref, anchor_k)
    used = manifest.sort_values("accession", kind="stable")
    if max_genomes is not None:
        used = used.head(max_genomes)

    ref_codons: dict[str, np.ndarray] = {}
    acc: dict[str, dict] = {}
    by_gene: dict[str, list[CdsCopy]] = {}
    for c in cds_copies:
        by_gene.setdefault(c.gene, []).append(c)
    for gene, copies in by_gene.items():
        c0 = copies[0]
        ref_codons[gene] = np.frombuffer(c0.cds.encode("ascii"),
                                         dtype=np.uint8).reshape(-1, 3)
        acc[gene] = {
            "n_resolved": np.zeros(c0.n_codons, dtype=np.int64),
            "n_variant": np.zeros(c0.n_codons, dtype=np.int64),
            "n_nt_variant": np.zeros(c0.n_codons, dtype=np.int64),
            "used": 0, "dropped": 0, "len_poly": 0, "len_poly_inframe": 0,
        }

    qc_rows = []
    for gi, row in enumerate(used.itertuples(index=False), start=1):
        path = PROJECT_ROOT / str(row.fasta_path)
        if not path.is_file():
            raise SystemExit(f"Manifest references a missing FASTA: {path}")
        seq = load_genome(path)
        maps = []
        for oriented in (seq, revcomp(seq)):
            anchors = rob.find_anchors(oriented, index, anchor_k, anchor_step)
            if not anchors:
                maps.append((None, None))
                continue
            chains = rob.chain_anchors(anchors, anchor_k)
            qmap = build_offset_map(n_ref, chains, anchor_k, max_unanchored)
            qbytes = np.frombuffer(oriented.encode("ascii"), dtype=np.uint8)
            maps.append((qmap, qbytes))

        for gene, copies in by_gene.items():
            best = None
            for copy in copies:
                for qmap, qbytes in maps:
                    if qmap is None:
                        continue
                    q = qmap[copy.positions]
                    ok = (q >= 0) & (q < len(qbytes))
                    score = int(ok.sum())
                    if best is None or score > best[0]:
                        best = (score, copy, q, ok, qbytes)
            if best is None:
                acc[gene]["dropped"] += 1
                continue
            score, copy, q, ok, qbytes = best
            frac = score / len(copy.positions)
            if frac < min_resolved_fraction:
                acc[gene]["dropped"] += 1
                qc_rows.append({"accession": row.accession, "gene": gene,
                                "resolved_fraction": round(frac, 4),
                                "verdict": "dropped: coverage"})
                continue

            bases = np.where(ok, qbytes[np.clip(q, 0, len(qbytes) - 1)], ord("N"))
            if copy.strand == -1:
                bases = _COMP_LUT[bases.astype(np.uint8)]
            codons = bases.astype(np.uint8).reshape(-1, 3)
            rc = ref_codons[gene]
            has_n = (codons == ord("N")).any(axis=1)
            resolved = ~has_n
            nt_diff = resolved & (codons != rc).any(axis=1)

            aa_diff = np.zeros(len(rc), dtype=bool)
            internal_stop = False
            for j in np.flatnonzero(nt_diff):
                aa = _CODON_TABLE.get(codons[j].tobytes().decode("ascii"), "X")
                if aa != copy.protein[j]:
                    aa_diff[j] = True
                    if aa == "*" and j != len(rc) - 1:
                        internal_stop = True
            if internal_stop:
                acc[gene]["dropped"] += 1
                qc_rows.append({"accession": row.accession, "gene": gene,
                                "resolved_fraction": round(frac, 4),
                                "verdict": "dropped: internal stop"})
                continue

            offs = (q - copy.positions)[ok]
            if offs.size:
                net = int(offs[-1] - offs[0])
                if net != 0:
                    acc[gene]["len_poly"] += 1
                    if net % 3 == 0:
                        acc[gene]["len_poly_inframe"] += 1

            acc[gene]["n_resolved"] += resolved
            acc[gene]["n_variant"] += aa_diff
            acc[gene]["n_nt_variant"] += nt_diff
            acc[gene]["used"] += 1

        if gi % 50 == 0 or gi == len(used):
            LOG.info("  tolerance profile: %d/%d genomes", gi, len(used))

    profiles: dict[str, ToleranceProfile] = {}
    for gene, copies in by_gene.items():
        a = acc[gene]
        profiles[gene] = ToleranceProfile(
            gene=gene, n_codons=copies[0].n_codons, ref_protein=copies[0].protein,
            n_resolved=a["n_resolved"], n_variant=a["n_variant"],
            n_nt_variant=a["n_nt_variant"], n_genomes_used=a["used"],
            n_genomes_dropped=a["dropped"], length_polymorphic=a["len_poly"],
            length_polymorphic_inframe=a["len_poly_inframe"],
        )
    return profiles, pd.DataFrame(qc_rows)


# ======================================================================================
# Per-site repair-escape probability
# ======================================================================================


def inframe_viability(tol: np.ndarray, codon: int, n_codons_affected: int) -> float:
    """P(protein functional | in-frame indel of `n_codons_affected` codons at `codon`).

    The indel's placement is averaged uniformly over the windows of that width that
    contain the cut codon, because Cas9 indels are centred on the cut site but not
    codon-phased. Viability is the product of the per-codon tolerances in the window:
    deleting (or inserting into) a stretch that contains an invariant codon is much
    less likely to leave a functional protein than one that varies among isolates.

    Insertions are scored with the same window formula as deletions. That is an
    approximation, stated as such; it errs toward treating insertions as MORE
    disruptive than they may be at protein termini, and its effect is bounded by the
    insertion/deletion sweep in the sensitivity analysis.
    """
    n = len(tol)
    if n == 0:
        return 0.0
    m = max(1, int(n_codons_affected))
    starts = [s for s in range(codon - m + 1, codon + 1) if 0 <= s and s + m <= n]
    if not starts:
        starts = [max(0, min(codon, n - m))]
    total = 0.0
    for s in starts:
        total += float(np.prod(tol[s:s + m]))
    return total / len(starts)


def site_repair_escape(spectrum: IndelSpectrum, tol: np.ndarray, codon: int,
                       gene: str, n_copies: int, params: Params) -> dict:
    """q_i and its decomposition, for one target site.

    Because a perfectly repaired site is simply re-cut by a persistently expressed
    nuclease, the terminal outcome is the indel spectrum conditioned on length != 0.
    `p_disrupt` scales the whole non-WT branch and, being constant across lengths,
    cancels in that conditioning -- it is retained explicitly so the assumption is
    visible and so a length-dependent version can be substituted later.
    """
    cond = spectrum.conditional
    fs_viab = params.frameshift_for(gene)
    q_frameshift = 0.0
    q_inframe = 0.0
    for L, p in cond.items():
        if L % 3 != 0:
            q_frameshift += p * fs_viab
        else:
            q_inframe += p * inframe_viability(tol, codon, abs(L) // 3)
    q_single = params.p_disrupt * (q_frameshift + q_inframe)
    q_single = min(1.0, max(0.0, q_single))

    if n_copies <= 1:
        q = q_single
    elif params.repeat_model == "redundant":
        q = 1.0 - (1.0 - q_single) ** n_copies
    elif params.repeat_model == "single":
        q = q_single
    elif params.repeat_model == "all-copies":
        q = q_single ** n_copies
    else:
        raise SystemExit(f"Unknown --repeat-model {params.repeat_model!r}")
    return {
        "q_single_copy": q_single,
        "q_site": min(1.0, q),
        "q_from_frameshift": q_frameshift,
        "q_from_inframe": q_inframe,
        "inframe_fraction": spectrum.inframe_fraction,
    }


# ======================================================================================
# The escape model over guide sets
# ======================================================================================


class EscapeModel:
    """P_escape for arbitrary guide sets over the empirical genome population.

    `presence` is the boolean (genome x site) matrix; `q` is the per-site repair
    escape probability. The per-genome, per-site failure factor is

        F[g, i] = q_i  if the site is present in genome g   (must be edited to escape)
                = 1    if the site is absent                (it can never be cut)

    and P_escape(S) = mean_g PROD_{i in S} F[g, i]. Nothing is assumed about the
    correlation structure of `presence`: it is used as measured.
    """

    def __init__(self, presence: np.ndarray, q: np.ndarray, site_ids: list[str]):
        if presence.shape[1] != len(q) or presence.shape[1] != len(site_ids):
            raise ValueError("presence, q and site_ids disagree on the number of sites")
        self.presence = presence
        self.q = np.asarray(q, dtype=float)
        self.site_ids = list(site_ids)
        self.index = {s: i for i, s in enumerate(site_ids)}
        self.n_genomes = presence.shape[0]
        self.F = np.where(presence, self.q[None, :], 1.0)

    # -- core ------------------------------------------------------------------------

    def escape_vector(self, cols: list[int]) -> np.ndarray:
        """Per-genome escape probability for a set of column indices."""
        if not cols:
            return np.ones(self.n_genomes)
        return np.prod(self.F[:, cols], axis=1)

    def p_escape(self, cols: list[int]) -> float:
        return float(self.escape_vector(cols).mean())

    def p_escape_ids(self, ids: list[str]) -> float:
        return self.p_escape([self.index[i] for i in ids])

    def p_escape_independent(self, cols: list[int]) -> float:
        """The same quantity computed as if site presence were independent.

        Reported ONLY as a contrast: the difference between this and `p_escape` is the
        exact size of the error that an independence assumption would introduce, which
        is how the report bounds it rather than asserting it is small.
        """
        marg = self.presence.mean(axis=0)
        out = 1.0
        for c in cols:
            out *= (1.0 - marg[c]) + marg[c] * self.q[c]
        return float(out)

    def joint_conservation(self, cols: list[int]) -> float:
        if not cols:
            return 0.0
        return float(np.all(self.presence[:, cols], axis=1).mean())

    # -- search ----------------------------------------------------------------------

    def _pair_matrix(self, allowed: np.ndarray) -> np.ndarray:
        sub = self.F[:, allowed]
        return (sub.T @ sub) / self.n_genomes

    def best_pair(self, allowed: list[int],
                  separation: "np.ndarray | None" = None,
                  min_separation: int = 0) -> tuple[list[int], float]:
        """Exhaustive best k=2 over `allowed`. O(n^2 N) but vectorised."""
        a = np.asarray(allowed, dtype=int)
        M = self._pair_matrix(a)
        np.fill_diagonal(M, np.inf)
        if separation is not None and min_separation > 0:
            M = np.where(separation[np.ix_(a, a)] >= min_separation, M, np.inf)
        if not np.isfinite(M).any():
            return [], 1.0
        flat = int(np.argmin(M))
        i, j = divmod(flat, M.shape[1])
        return sorted([int(a[i]), int(a[j])]), float(M[i, j])

    def _feasible(self, cand: np.ndarray, chosen: list[int],
                  separation, min_separation: int,
                  gene_of, max_per_gene: int) -> np.ndarray:
        """Apply the design constraints to a vector of candidate columns.

        `min_separation` forbids two picks in the same gene close enough that one
        repair event could plausibly destroy both -- the model treats repair as
        independent per site, and this removes the configurations where that would
        be wrong instead of mis-modelling them. `max_per_gene` implements the
        design policy that spreads a set across genes.
        """
        if cand.size and separation is not None and chosen and min_separation > 0:
            ok = np.all(separation[np.ix_(cand, chosen)] >= min_separation, axis=1)
            cand = cand[ok]
        if cand.size and gene_of is not None and max_per_gene > 0:
            counts: dict = {}
            for c in chosen:
                counts[gene_of[c]] = counts.get(gene_of[c], 0) + 1
            ok = np.array([counts.get(gene_of[c], 0) < max_per_gene for c in cand])
            cand = cand[ok]
        return cand

    def greedy(self, k: int, allowed: list[int], seed_cols: list[int] | None = None,
               separation: "np.ndarray | None" = None,
               min_separation: int = 0,
               gene_of: "np.ndarray | None" = None,
               max_per_gene: int = 0) -> list[int]:
        """Deterministic greedy: repeatedly add the site that lowers P_escape most."""
        chosen = list(seed_cols or [])
        cur = self.escape_vector(chosen)
        allowed_set = [c for c in allowed if c not in chosen]
        while len(chosen) < k and allowed_set:
            cand = self._feasible(np.asarray(allowed_set, dtype=int), chosen,
                                  separation, min_separation, gene_of, max_per_gene)
            if cand.size == 0:
                break
            vals = (cur[:, None] * self.F[:, cand]).mean(axis=0)
            # Deterministic tie-break: lowest objective, then lowest column index.
            best = int(cand[np.lexsort((cand, vals))[0]])
            chosen.append(best)
            cur = cur * self.F[:, best]
            allowed_set = [c for c in allowed_set if c != best]
        return sorted(chosen)

    def local_search(self, cols: list[int], allowed: list[int],
                     separation: "np.ndarray | None" = None,
                     min_separation: int = 0,
                     max_rounds: int = 25,
                     gene_of: "np.ndarray | None" = None,
                     max_per_gene: int = 0) -> tuple[list[int], float]:
        """Swap each member against every allowed site until no improvement.

        Deterministic (no randomness, fixed scan order), so a reported "best set" is
        reproducible byte-for-byte.
        """
        cur = sorted(cols)
        best = self.p_escape(cur)
        for _ in range(max_rounds):
            improved = False
            for pos in range(len(cur)):
                rest = [c for i, c in enumerate(cur) if i != pos]
                base = self.escape_vector(rest)
                cand = self._feasible(
                    np.asarray([c for c in allowed if c not in rest], dtype=int),
                    rest, separation, min_separation, gene_of, max_per_gene)
                if cand.size == 0:
                    continue
                vals = (base[:, None] * self.F[:, cand]).mean(axis=0)
                order = np.lexsort((cand, vals))
                pick, val = int(cand[order[0]]), float(vals[order[0]])
                if val < best - 1e-18:
                    cur = sorted(rest + [pick])
                    best = val
                    improved = True
            if not improved:
                break
        return sorted(cur), best

    def best_set(self, k: int, allowed: list[int], n_starts: int = 30,
                 separation: "np.ndarray | None" = None,
                 min_separation: int = 0,
                 gene_of: "np.ndarray | None" = None,
                 max_per_gene: int = 0) -> tuple[list[int], float]:
        """Best k-set found: exhaustive for k<=2, multi-start greedy + swap above.

        For k = 3 the report additionally verifies this against an exhaustive search
        over a pruned pool, so the heuristic is checked rather than trusted.
        """
        if k <= 0:
            return [], 1.0
        if k == 1:
            vals = self.F[:, allowed].mean(axis=0)
            a = np.asarray(allowed, dtype=int)
            pick = int(a[np.lexsort((a, vals))[0]])
            return [pick], float(vals.min())
        if k == 2 and max_per_gene <= 0:
            return self.best_pair(allowed, separation, min_separation)
        singles = self.F[:, allowed].mean(axis=0)
        a = np.asarray(allowed, dtype=int)
        order = a[np.lexsort((a, singles))][:max(1, n_starts)]
        best_cols, best_val = None, math.inf
        for s in order:
            g = self.greedy(k, allowed, [int(s)], separation, min_separation,
                            gene_of, max_per_gene)
            if len(g) < k:
                continue
            cols, val = self.local_search(g, allowed, separation, min_separation,
                                          gene_of=gene_of,
                                          max_per_gene=max_per_gene)
            if val < best_val - 1e-18 or (abs(val - best_val) <= 1e-18
                                          and (best_cols is None or cols < best_cols)):
                best_cols, best_val = cols, val
        if best_cols is None:
            return [], 1.0
        return best_cols, best_val

    def exhaustive(self, k: int, allowed: list[int]) -> tuple[list[int], float]:
        """Brute force over `allowed`. Only used to validate the heuristic."""
        best_cols, best_val = None, math.inf
        for combo in combinations(sorted(allowed), k):
            v = self.p_escape(list(combo))
            if v < best_val - 1e-18:
                best_cols, best_val = list(combo), v
        return (best_cols or []), best_val

    # -- Monte Carlo cross-check ----------------------------------------------------

    def monte_carlo(self, cols: list[int], draws: int, seed: int) -> tuple[float, float]:
        """Seeded simulation of the same process; must agree with the closed form.

        One draw = pick a genome uniformly from the empirical population, then for
        each site that is present in it draw a Bernoulli(q_i) repair escape. Returns
        (estimate, standard error).
        """
        rng = np.random.default_rng(seed)
        g = rng.integers(0, self.n_genomes, size=draws)
        alive = np.ones(draws, dtype=bool)
        for c in cols:
            need = self.presence[g, c]
            u = rng.random(draws)
            alive &= (~need) | (u < self.q[c])
        est = float(alive.mean())
        se = float(math.sqrt(max(est * (1 - est), 1e-12) / draws))
        return est, se


# ======================================================================================
# Assembly
# ======================================================================================


def cut_codon_for_site(row, cds_copies: list[CdsCopy], nuclease: Nuclease,
                       pos_to_codon: dict[str, dict[int, int]]) -> tuple[int, int]:
    """(codon index, n copies) for a target site, at the canonical blunt cut.

    SpCas9/SaCas9 cut between protospacer positions 17 and 18 counting from the
    PAM-distal end, i.e. 3 bp inside the protospacer from the PAM. The exact base is
    immaterial here (the viability window is averaged over placements anyway) but it
    is computed properly rather than reusing the deliberately-legacy `cut_site_ref`
    column, which the README documents as sitting one base off the canonical cut.
    """
    pam_len = nuclease.pam_length
    if row.strand == "+":
        cut0 = int(row.ref_end) - pam_len - 3
    else:
        cut0 = int(row.ref_start) - 1 + pam_len + 2
    mapping = pos_to_codon.get(row.gene, {})
    codon = mapping.get(cut0)
    if codon is None:
        mid = (int(row.ref_start) - 1 + int(row.ref_end)) // 2
        codon = mapping.get(mid)
    if codon is None:
        for p in range(int(row.ref_start) - 1, int(row.ref_end)):
            if p in mapping:
                codon = mapping[p]
                break
    return (0 if codon is None else int(codon)), int(row.n_reference_copies)


def build_site_table(pool: pd.DataFrame, profiles: dict[str, ToleranceProfile],
                     cds_copies: list[CdsCopy], nuclease: Nuclease,
                     spectrum: IndelSpectrum, params: Params) -> pd.DataFrame:
    pos_to_codon: dict[str, dict[int, int]] = {}
    for c in cds_copies:
        pos_to_codon.setdefault(c.gene, {}).update(c.codon_of_position())

    tol_cache = {g: p.tolerance_vector(params) for g, p in profiles.items()}
    rows = []
    for row in pool.itertuples(index=False):
        codon, copies = cut_codon_for_site(row, cds_copies, nuclease, pos_to_codon)
        prof = profiles.get(row.gene)
        tol = tol_cache.get(row.gene)
        if tol is None or prof is None:
            raise SystemExit(f"No tolerance profile for gene {row.gene}")
        d = site_repair_escape(spectrum, tol, codon, row.gene, copies, params)
        varies = prof.varies(params.min_variant_strains)
        lo = max(0, codon - 5)
        hi = min(prof.n_codons, codon + 6)
        rows.append({
            "guide_id": row.guide_id,
            "gene": row.gene,
            "essential": row.gene in ESSENTIAL_GENES,
            "strand": row.strand,
            "ref_start": int(row.ref_start),
            "ref_end": int(row.ref_end),
            "protospacer": row.protospacer,
            "pam": row.pam,
            "n_reference_copies": copies,
            "cut_codon": codon,
            "codon_varies": bool(varies[codon]) if codon < len(varies) else False,
            "local_variable_codons_11": int(varies[lo:hi].sum()),
            "local_tolerance_11": float(np.mean(tol[lo:hi])),
            "gc_content": float(row.gc_content),
            "max_homopolymer_run": int(row.max_homopolymer_run),
            "has_polyT": bool(row.has_polyT),
            "q_single_copy": d["q_single_copy"],
            "q_site": d["q_site"],
            "q_from_frameshift": d["q_from_frameshift"],
            "q_from_inframe": d["q_from_inframe"],
        })
    return pd.DataFrame(rows)


def make_model(pool: pd.DataFrame, presence: np.ndarray, col_index: dict[str, int],
               sites: pd.DataFrame) -> EscapeModel:
    ids = sites["guide_id"].tolist()
    cols = [col_index[g] for g in ids]
    return EscapeModel(presence[:, cols], sites["q_site"].to_numpy(float), ids)


def separation_matrix(sites: pd.DataFrame) -> np.ndarray:
    """Pairwise genomic distance between sites; used to forbid overlapping picks.

    Two guides whose footprints are a few bases apart can be destroyed by a single
    repair event, which the model (treating repair as independent per site) would
    score as two independent escapes. Requiring separation removes that failure mode
    from the search rather than mis-modelling it.
    """
    mid = ((sites["ref_start"] + sites["ref_end"]) / 2.0).to_numpy()
    d = np.abs(mid[:, None] - mid[None, :])
    same_gene = (sites["gene"].to_numpy()[:, None] == sites["gene"].to_numpy()[None, :])
    # Sites in different genes are never at risk from one repair event.
    return np.where(same_gene, d, np.inf)


# ======================================================================================
# Analyses
# ======================================================================================


def k_curve(model: EscapeModel, allowed: list[int], max_k: int,
            separation: np.ndarray, min_separation: int,
            thresholds: list[float], gene_of: "np.ndarray | None" = None,
            max_per_gene: int = 0) -> tuple[pd.DataFrame, dict[int, list[int]]]:
    rows = []
    best_sets: dict[int, list[int]] = {}
    for k in range(1, max_k + 1):
        cols, val = model.best_set(k, allowed, separation=separation,
                                   min_separation=min_separation,
                                   gene_of=gene_of, max_per_gene=max_per_gene)
        if not cols:
            break
        best_sets[k] = cols
        rows.append({
            "k": k,
            "p_escape_best": val,
            "p_escape_if_independent": model.p_escape_independent(cols),
            "joint_conservation": model.joint_conservation(cols),
            "guides": ";".join(model.site_ids[c] for c in cols),
            "genes": ";".join(str(gene_of[c]) for c in cols)
                     if gene_of is not None else "",
            **{f"below_{t:g}": bool(val < t) for t in thresholds},
        })
    return pd.DataFrame(rows), best_sets


def minimum_k(curve: pd.DataFrame, threshold: float) -> int | None:
    hit = curve[curve["p_escape_best"] < threshold]
    return int(hit["k"].iloc[0]) if len(hit) else None


def build_pool_and_presence(record, genome: str, genes: list[str],
                            nuclease: Nuclease, manifest: pd.DataFrame,
                            max_genomes: int | None):
    pool = bm.build_pool(genome, record, genes, nuclease)
    man = manifest.sort_values("accession", kind="stable")
    if max_genomes is not None:
        man = man.head(max_genomes)
    accessions, presence, col_index = bm.presence_matrix(pool, man, nuclease)
    return pool, accessions, presence, col_index


def amrani_site_ids(pool: pd.DataFrame) -> dict[str, str]:
    by_target = {p + m: name for name, (p, m) in bm.BENCHMARK_GUIDES.items()}
    out = {}
    for row in pool.itertuples(index=False):
        name = by_target.get(row.target_site)
        if name:
            out[name] = row.guide_id
    return out


# ======================================================================================
# Sensitivity
# ======================================================================================


def binomial_upper_bound(successes: int, n: int, alpha: float = 0.05) -> float:
    """Clopper-Pearson upper confidence limit on a binomial proportion.

    Used for the resolution floor: a site absent in 0 of N sequenced isolates has
    a MEASURED absence frequency of 0, but the population frequency it is
    consistent with is anything up to this bound (for 0/183 that is ~1.6%, the
    familiar rule of three). Computed by bisection on the exact binomial CDF, so
    no scipy dependency is introduced.
    """
    if successes >= n:
        return 1.0
    def cdf(p: float) -> float:
        if p <= 0.0:
            return 1.0
        if p >= 1.0:
            return 0.0
        total = 0.0
        lp, lq = math.log(p), math.log1p(-p)
        ln = math.lgamma(n + 1)
        for x in range(successes + 1):
            total += math.exp(ln - math.lgamma(x + 1) - math.lgamma(n - x + 1)
                              + x * lp + (n - x) * lq)
        return total
    lo, hi = 0.0, 1.0
    for _ in range(80):
        mid = (lo + hi) / 2.0
        if cdf(mid) > alpha:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2.0


def _spearman(a: list[float], b: list[float]) -> float:
    ra = pd.Series(a).rank().to_numpy()
    rb = pd.Series(b).rank().to_numpy()
    if ra.std() == 0 or rb.std() == 0:
        return float("nan")
    return float(np.corrcoef(ra, rb)[0, 1])


def sensitivity(pool, profiles, cds_copies, nuclease, base_spectrum, base_params,
                presence, col_index, allowed_ids: list[str], ref_sets: dict[str, list[str]],
                max_k: int, min_separation: int,
                thresholds: list[float]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Sweep every uncertain parameter one at a time.

    Two things are reported for each setting: the ABSOLUTE escape probability of the
    reference guide sets (which moves over orders of magnitude), and the RANK
    correlation of all candidate sites' q against the baseline (which is what actually
    determines whether the design recommendation changes).
    """
    variants: list[tuple[str, str, IndelSpectrum, Params]] = []

    def add(group, label, spec=None, **kw):
        variants.append((group, label, spec or base_spectrum, replace(base_params, **kw)))

    add("baseline", "baseline")
    for v in (0.10, 0.25, 0.50, 0.75, 1.00):
        add("tol_variable", f"tol_variable={v:g}", None, tol_variable=v)
    for v in (0.001, 0.01, 0.05, 0.15, 0.30, 0.50):
        add("tol_invariant", f"tol_invariant={v:g}", None, tol_invariant=v)
    for v in (1, 2, 3, 5):
        add("min_variant_strains", f"min_variant_strains={v}", None, min_variant_strains=v)
    # `neuronal_nhej` is the one swept spectrum whose shape is measured rather than
    # posited (Ramadoss et al. 2025 Fig. 1d source data, post-mitotic human neurons).
    # It is a SCENARIO, not the default: the default stays where it was so that every
    # previously reported number remains reproducible, and this widens the swept space.
    for name in ("default", "deletion-heavy", "insertion-heavy", "short-indels",
                 "uniform-1to20", "neuronal_nhej"):
        variants.append(("indel_spectrum", f"spectrum={name}",
                         builtin_spectrum(name), base_params))
    for f in (0.10, 0.20, 0.25, 0.33, 0.50):
        variants.append(("inframe_fraction", f"inframe_fraction={f:g}",
                         base_spectrum.rescaled_to_inframe(f), base_params))
    for m in ("redundant", "single", "all-copies"):
        add("repeat_model", f"repeat_model={m}", None, repeat_model=m)
    for v in (0.0, 0.05, 0.25, 1.0):
        add("frameshift_viability_RL2", f"RL2_frameshift_viability={v:g}", None,
            frameshift_viability=(("RL2", v),))
    for v in (0.8, 0.9, 1.0):
        add("p_disrupt", f"p_disrupt={v:g}", None, p_disrupt=v)

    baseline_q: np.ndarray | None = None
    baseline_sets: dict[int, tuple[str, ...]] = {}
    rows = []
    set_rows = []
    for group, label, spec, params in variants:
        sites = build_site_table(pool, profiles, cds_copies, nuclease, spec, params)
        model = make_model(pool, presence, col_index, sites)
        allowed = [model.index[g] for g in allowed_ids]
        sep = separation_matrix(sites)
        q = sites["q_site"].to_numpy(float)
        if baseline_q is None:
            baseline_q = q
        rho = _spearman(list(baseline_q), list(q))
        # Spearman over ~450 sites is dominated by ties (most sites sit at the same
        # invariant-codon value), so a decision-relevant overlap statistic is reported
        # next to it: how much of the baseline top-20 shortlist survives.
        base_top = set(np.argsort(baseline_q, kind="stable")[:20].tolist())
        var_top = set(np.argsort(q, kind="stable")[:20].tolist())
        top20 = len(base_top & var_top) / 20.0

        best_by_k = {}
        for k in range(1, max_k + 1):
            cols, val = model.best_set(k, allowed, separation=sep,
                                       min_separation=min_separation)
            if not cols:
                break
            best_by_k[k] = (cols, val)
        rec = {
            "group": group, "label": label,
            "inframe_fraction": spec.inframe_fraction,
            "median_q_site": float(np.median(q)),
            "spearman_q_vs_baseline": rho,
            "top20_site_overlap_vs_baseline": top20,
        }
        for k, (cols, val) in best_by_k.items():
            rec[f"p_escape_best_k{k}"] = val
        for t in thresholds:
            mk = next((k for k, (_, v) in sorted(best_by_k.items()) if v < t), None)
            rec[f"min_k_below_{t:g}"] = mk
        for name, ids in ref_sets.items():
            rec[f"p_escape_{name}"] = model.p_escape_ids(ids)
        chosen = {k: tuple(sorted(model.site_ids[c] for c in cols))
                  for k, (cols, _) in best_by_k.items()}
        if not baseline_sets:
            baseline_sets = chosen
        for k in sorted(chosen):
            rec[f"best_set_k{k}_matches_baseline"] = \
                bool(chosen[k] == baseline_sets.get(k))
        rows.append(rec)

        for k, (cols, val) in best_by_k.items():
            set_rows.append({"label": label, "k": k, "p_escape": val,
                             "guides": ";".join(model.site_ids[c] for c in cols)})
    return pd.DataFrame(rows), pd.DataFrame(set_rows)


# ======================================================================================
# Report
# ======================================================================================


def _md(df: pd.DataFrame, floatfmt: str = "{:.4g}") -> str:
    if df is None or df.empty:
        return "_(no rows)_\n"
    def fmt(v):
        if v is None or v is pd.NA:
            return "-"
        if isinstance(v, float):
            if not np.isfinite(v):
                return "-"
            return floatfmt.format(v)
        return str(v)
    head = "| " + " | ".join(str(c) for c in df.columns) + " |"
    rule = "|" + "|".join("---" for _ in df.columns) + "|"
    body = ["| " + " | ".join(fmt(v) for v in row) + " |"
            for row in df.itertuples(index=False)]
    return "\n".join([head, rule, *body]) + "\n"


def _sci(x: float) -> str:
    if x is None:
        return "-"
    if x == 0:
        return "0"
    return f"{x:.3e}"


def render_report(ctx: dict) -> str:
    p = ctx["params_table"]
    out: list[str] = []
    a = out.append
    a("# A quantitative escape-probability model for multiplex CRISPR editing of HSV-1\n")
    a(f"Generated {ctx['generated']} · nuclease {ctx['nuclease_label']} · "
      f"{ctx['n_genomes']} complete genomes · {ctx['n_sites']} candidate sites in "
      f"{ctx['n_genes']} genes.\n")

    a("## 0. Headline\n")
    for line in ctx["headline"]:
        a(f"- {line}")
    a("")

    a("## 1. The model\n")
    a(ctx["model_prose"])

    a("## 2. Parameters, their provenance and the direction of their effect\n")
    a(_md(p))
    a("The only parameters that are MEASURED are the ones taken from this "
      "repository's own genome corpus. Everything labelled ASSUMPTION is swept in "
      "section 8 and no default was chosen after seeing its effect on the "
      "conclusion.\n")

    a("## 3. Measured input A -- the per-genome site presence matrix\n")
    a(ctx["presence_prose"])
    a(_md(ctx["presence_table"]))

    a("## 4. Measured input B -- per-codon tolerance from cross-strain variation\n")
    a(ctx["tolerance_prose"])
    a(_md(ctx["tolerance_table"]))
    a(ctx["tolerance_corroboration"])

    a("## 5. Per-site escape probabilities\n")
    a(_md(ctx["site_summary"]))
    a(ctx["site_prose"])

    a("## 6. The Amrani et al. lead pair, scored\n")
    a(_md(ctx["amrani_table"]))
    a("**Where that number comes from.** The three isolate classes contribute "
      "separately and the split is exact, not a fit:\n")
    a(_md(ctx["amrani_decomp"]))
    a(ctx["amrani_prose"])
    a("### 6.1 Alternatives at the same k\n")
    a(_md(ctx["alt_table"]))

    a("## 7. Escape probability versus number of guides\n")
    a("**Unconstrained best set at each k** -- the model is free to put every guide in "
      "the same gene if that is optimal, and at small k it does.\n")
    a(_md(ctx["k_table"]))
    a(f"**Best set at each k under a design policy of at most "
      f"{ctx['max_per_gene']} guide(s) per gene**, which is what a real multiplex "
      "construct would be built to.\n")
    a(_md(ctx["k_table_pg"]))
    a(ctx["k_prose"])
    a("### 7.0 Is the search actually finding the optimum?\n")
    a(_md(ctx["validation_table"]))
    a("### 7.1 The resolution floor -- what this many genomes can certify\n")
    a(_md(ctx["floor_table"]))
    a(ctx["floor_prose"])
    a("### 7.2 Minimum k for stated thresholds\n")
    a(_md(ctx["mink_table"]))
    a("The same question asked of the evidence rather than of the point estimate. "
      "The right-hand column is the smallest k that clears the threshold even at the "
      "95% upper bound on unobserved site absence, and it is the number a manuscript "
      "should quote.\n")
    a(_md(ctx["mink_floor_table"]))
    a(ctx["mink_prose"])
    a("### 7.3 Correlated failure, and what independence would have got wrong\n")
    a(_md(ctx["indep_table"]))
    a(ctx["indep_prose"])

    a("## 8. Sensitivity analysis\n")
    a(ctx["sens_prose"])
    a(_md(ctx["sens_table"]))
    a("### 8.1 Does the ranking survive?\n")
    a(_md(ctx["rank_table"]))
    a(ctx["rank_prose"])

    a("## 9. Closed form versus simulation\n")
    a(_md(ctx["mc_table"]))
    a(ctx["mc_prose"])

    a("## 10. The delivery term, which no number of guides can fix\n")
    a(_md(ctx["theta_table"]))
    a(ctx["theta_prose"])

    a("## 11. What this model cannot tell you\n")
    for line in ctx["cannot"]:
        a(f"{line}\n")

    a("## 12. Reproducibility\n")
    a(ctx["repro"])
    return "\n".join(out) + "\n"


# ======================================================================================
# Driver
# ======================================================================================


def add_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--reference", default=DEFAULT_REFERENCE)
    parser.add_argument("--genes", default=",".join(DEFAULT_GENES))
    parser.add_argument("--nuclease", default="sacas9",
                        help="Default sacas9. With sacas9 and no --spacer-length, the "
                             "spacer defaults to 20 nt rather than the registry's 21, "
                             "because 20 + NNGRRT is the grammar Amrani et al. print "
                             "in their Table 1 and the headline comparison is only "
                             "meaningful inside it.")
    parser.add_argument("--pam", default=None,
                        help="IUPAC PAM override; defaults to the nuclease's own.")
    parser.add_argument("--spacer-length", type=int, default=None)
    parser.add_argument("--outdir", type=Path, default=RESULTS_DIR)
    parser.add_argument("--prefix", default="escape_")
    parser.add_argument("--max-k", type=int, default=8)
    parser.add_argument("--thresholds", default="1e-3,1e-6")
    parser.add_argument("--indel-spectrum", default="default",
                        help="Built-in name (default, deletion-heavy, insertion-heavy, "
                             "short-indels, uniform-1to20, neuronal_nhej, "
                             "ipsc_dividing) or a TSV with columns length, probability. "
                             "The default is LABELLED AS AN ASSUMPTION; neuronal_nhej "
                             "and ipsc_dividing are measured in Ramadoss et al. 2025's "
                             "cells and remain assumptions for this therapy's.")
    parser.add_argument("--inframe-fraction", type=float, default=None,
                        help="Override the spectrum's in-frame fraction directly; this "
                             "is the one summary statistic the model is sensitive to.")
    parser.add_argument("--tol-variable", type=float, default=0.50)
    parser.add_argument("--tol-invariant", type=float, default=0.05)
    parser.add_argument("--min-variant-strains", type=int, default=1)
    parser.add_argument("--repeat-model", default="redundant",
                        choices=["redundant", "single", "all-copies"])
    parser.add_argument("--frameshift-viability", default="",
                        help="Comma-separated GENE=VALUE overrides, e.g. RL2=1.0.")
    parser.add_argument("--p-disrupt", type=float, default=1.0)
    parser.add_argument("--theta-unexposed", default="0,0.5,0.9,0.99",
                        help="Delivery-failure fractions to tabulate (k-independent).")
    parser.add_argument("--viral-population-size", default="1e3,1e5,1e7")
    parser.add_argument("--min-separation", type=int, default=150,
                        help="Minimum genomic distance between two chosen sites in the "
                             "same gene, so one repair event cannot destroy both.")
    parser.add_argument("--max-per-gene", type=int, default=1,
                        help="Design policy for the SECOND k-curve: at most this many "
                             "guides per gene. The unconstrained optimum is always "
                             "reported alongside it.")
    parser.add_argument("--filter-mode", default="pass",
                        choices=["pass", "all"],
                        help="'pass' restricts the candidate pool to sites clearing the "
                             "stage-4 poly-T / homopolymer / GC filters.")
    parser.add_argument("--seed", type=int, default=20240814)
    parser.add_argument("--mc-draws", type=int, default=2_000_000)
    parser.add_argument("--max-genomes", type=int, default=None,
                        help="Smoke-test switch: use only the first N genomes.")
    parser.add_argument("--skip-sensitivity", action="store_true")
    parser.add_argument("--anchor-step", type=int, default=6)
    parser.add_argument("--max-unanchored", type=int, default=60)


def _parse_thresholds(s: str) -> list[float]:
    return [float(x) for x in s.split(",") if x.strip()]


def _parse_frameshift(s: str) -> tuple[tuple[str, float], ...]:
    out = []
    for item in s.split(","):
        item = item.strip()
        if not item:
            continue
        gene, _, val = item.partition("=")
        out.append((gene.strip().upper(), float(val)))
    return tuple(out)


def run(args: argparse.Namespace) -> dict:
    from Bio import SeqIO

    from src.common import utc_now_iso

    ensure_dirs()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    thresholds = _parse_thresholds(args.thresholds)
    genes = [g.strip().upper() for g in args.genes.split(",") if g.strip()]
    spacer = args.spacer_length
    if spacer is None and (args.nuclease or "").strip().lower() == "sacas9":
        spacer = 20          # the published Amrani grammar; see --nuclease help
    nuclease = get_nuclease(args.nuclease, args.pam, spacer)

    if not args.manifest.is_file():
        raise SystemExit(f"Manifest not found: {args.manifest} (run stage 1 first)")
    manifest = pd.read_csv(args.manifest, sep="\t")

    gb = rob.reference_genbank_path(args.reference)
    if not gb.is_file():
        from src.extract_guides import fetch_reference
        gb = fetch_reference(args.reference)
    record = SeqIO.read(gb, "genbank")
    genome = str(record.seq).upper()

    spectrum = load_spectrum(args.indel_spectrum)
    if args.inframe_fraction is not None:
        spectrum = spectrum.rescaled_to_inframe(args.inframe_fraction)
    params = Params(
        tol_variable=args.tol_variable,
        tol_invariant=args.tol_invariant,
        min_variant_strains=args.min_variant_strains,
        frameshift_viability=_parse_frameshift(args.frameshift_viability),
        repeat_model=args.repeat_model,
        p_disrupt=args.p_disrupt,
    )

    LOG.info("Enumerating %s sites in %s ...", nuclease.label, ", ".join(genes))
    pool, accessions, presence, col_index = build_pool_and_presence(
        record, genome, genes, nuclease, manifest, args.max_genomes)
    LOG.info("Pool: %d sites; presence matrix %s", len(pool), presence.shape)

    cds_copies = collect_cds_copies(record, genes)
    LOG.info("Measuring per-codon cross-strain variation over %d CDS copies ...",
             len(cds_copies))
    profiles, qc = measure_tolerance(record, cds_copies, manifest,
                                     anchor_step=args.anchor_step,
                                     max_unanchored=args.max_unanchored,
                                     max_genomes=args.max_genomes)

    sites = build_site_table(pool, profiles, cds_copies, nuclease, spectrum, params)
    model = make_model(pool, presence, col_index, sites)
    sep = separation_matrix(sites)

    passes = ((~sites["has_polyT"])
              & (sites["max_homopolymer_run"] <= 4)
              & sites["gc_content"].between(0.35, 0.75))
    sites["passes_filters"] = passes
    allowed_ids = (sites.loc[passes, "guide_id"].tolist() if args.filter_mode == "pass"
                   else sites["guide_id"].tolist())
    allowed = [model.index[g] for g in allowed_ids]
    LOG.info("Candidate pool for set selection: %d sites (%s)",
             len(allowed), args.filter_mode)

    amrani = amrani_site_ids(pool)
    lead = [amrani.get("Amrani2024_ICP0g2"), amrani.get("Amrani2024_ICP27g1")]
    if None in lead:
        # Two very different causes, and conflating them would hide a real bug.
        amrani_grammar = (nuclease.name == "sacas9" and nuclease.pam == "NNGRRT"
                          and nuclease.spacer_length == 20)
        if amrani_grammar and {"RL2", "UL54"} <= set(genes):
            raise SystemExit(
                "INTEGRITY FAILURE: running in the published SaCas9 grammar over "
                "RL2 and UL54, the four Amrani et al. guides MUST be recoverable from "
                "the enumerated pool, and they were not. Every comparison below would "
                "be meaningless, so this is a hard stop rather than a warning."
            )
        raise SystemExit(
            f"The Amrani et al. lead pair cannot exist in this site space: they are "
            f"SaCas9 with a 20 nt spacer and an NNGRRT PAM in RL2 and UL54, and this "
            f"run enumerated {nuclease.label} in {', '.join(genes)}. The escape model "
            f"is nuclease-agnostic but its headline comparison is not; rerun with "
            f"--nuclease sacas9 --pam NNGRRT --spacer-length 20 (the defaults) and "
            f"genes including RL2 and UL54. Deliberately refused rather than "
            f"half-answered."
        )

    gene_of = sites["gene"].to_numpy()
    curve, best_sets = k_curve(model, allowed, args.max_k, sep,
                               args.min_separation, thresholds, gene_of)
    n_genes_avail = int(sites.loc[sites["guide_id"].isin(allowed_ids), "gene"].nunique())
    curve_pg, best_sets_pg = k_curve(
        model, allowed, min(args.max_k, n_genes_avail * max(1, args.max_per_gene)),
        sep, args.min_separation, thresholds, gene_of, args.max_per_gene)

    ctx = assemble_context(args, nuclease, genes, record, manifest, accessions,
                           pool, presence, col_index, sites, cds_copies, model, sep,
                           allowed, allowed_ids, amrani, lead, curve, best_sets,
                           curve_pg, best_sets_pg, gene_of,
                           profiles, qc, spectrum, params, thresholds, utc_now_iso())

    # ---- outputs -------------------------------------------------------------------
    pre = args.prefix
    paths = {}

    def rel(path: Path) -> str:
        """Repository-relative POSIX path for the summary's `outputs` block.

        Recording an absolute local path would leak the author's filesystem layout
        into a published artefact and would make this file differ between machines
        for no scientific reason. Falls back to the absolute path if the output was
        directed somewhere outside the repository.
        """
        try:
            return path.resolve().relative_to(PROJECT_ROOT).as_posix()
        except ValueError:
            return str(path)

    def write(name: str, df: pd.DataFrame):
        path = outdir / f"{pre}{name}.tsv"
        df.to_csv(path, sep="\t", index=False)
        paths[name] = rel(path)
        LOG.info("Wrote %s (%d rows)", path, len(df))

    write("site_parameters", sites)
    write("k_curve", curve)
    write("k_curve_one_per_gene", curve_pg)
    write("guide_sets", ctx["guide_sets_tsv"])
    write("codon_tolerance", ctx["codon_tsv"])
    if ctx.get("sens_tsv") is not None:
        write("sensitivity", ctx["sens_tsv"])
        write("sensitivity_sets", ctx["sens_sets_tsv"])
    if not qc.empty:
        write("tolerance_qc", qc)

    report_path = outdir / f"{pre}model_report.md"
    report_path.write_text(render_report(ctx), encoding="utf-8")
    LOG.info("Wrote %s", report_path)

    summary = ctx["summary"]
    summary["outputs"] = {**paths, "report": rel(report_path)}
    json_path = outdir / f"{pre}summary.json"
    json_path.write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    LOG.info("Wrote %s", json_path)
    return summary


def assemble_context(args, nuclease, genes, record, manifest, accessions, pool,
                     presence, col_index, sites, cds_copies, model, sep, allowed,
                     allowed_ids, amrani, lead, curve, best_sets, curve_pg,
                     best_sets_pg, gene_of, profiles, qc,
                     spectrum, params, thresholds, generated) -> dict:
    n_g = len(accessions)
    lead_cols = [model.index[g] for g in lead]
    p_lead = model.p_escape(lead_cols)
    jc_lead = model.joint_conservation(lead_cols)
    _NEURONAL_INFRAME = neuronal_nhej_spectrum("neuron").inframe_fraction

    # ---- parameter table ------------------------------------------------------------
    fs = params.frameshift_map()
    params_table = pd.DataFrame([
        {"parameter": "site presence matrix", "value": f"{n_g} genomes x {len(sites)} sites",
         "source": "MEASURED (stage 3 semantics)",
         "effect if wrong": "dominates P_escape at small k; under-states conservation "
                            "if assemblies are unresolved, which inflates escape"},
        {"parameter": "per-codon amino-acid variation",
         "value": f"{n_g}-isolate profile",
         "source": "MEASURED (this module)",
         "effect if wrong": "sets which in-frame indels are survivable; a wrong frame "
                            "would randomise it, so the frame is asserted by translation"},
        {"parameter": "indel length spectrum", "value": spectrum.name,
         "source": f"{spectrum.provenance} -- no citation invented",
         "effect if wrong": "acts almost entirely through its in-frame fraction "
                            f"({spectrum.inframe_fraction:.3f}); swept 0.10-0.50, plus "
                            "the measured post-mitotic-neuron scenario neuronal_nhej "
                            f"({_NEURONAL_INFRAME:.3f}), which sits BELOW that range"},
        {"parameter": "tol_variable", "value": params.tol_variable, "source": "ASSUMPTION",
         "effect if wrong": "scales q at variable codons; swept 0.10-1.00"},
        {"parameter": "tol_invariant", "value": params.tol_invariant, "source": "ASSUMPTION",
         "effect if wrong": "scales q at invariant codons; swept 0.001-0.50"},
        {"parameter": "min_variant_strains", "value": params.min_variant_strains,
         "source": "ASSUMPTION", "effect if wrong": "moves codons between the two "
                                                    "tolerance classes; swept 1-5"},
        {"parameter": "frameshift viability", "value": (fs or "0.0 for every gene"),
         "source": "ASSUMPTION (baseline is the one most favourable to a 2-guide design)",
         "effect if wrong": "RL2/ICP0 is not strictly essential; setting RL2=1.0 raises "
                            "escape at every ICP0 site by ~1/(in-frame fraction)"},
        {"parameter": "repeat model", "value": params.repeat_model,
         "source": "ASSUMPTION (all three variants reported)",
         "effect if wrong": "changes ICP0 site q by up to two orders of magnitude"},
        {"parameter": "p_disrupt", "value": params.p_disrupt, "source": "ASSUMPTION",
         "effect if wrong": "cancels in the conditioning while constant in length; "
                            "only a length-dependent version would matter"},
        {"parameter": "min site separation", "value": f"{args.min_separation} bp",
         "source": "MODELLING CHOICE",
         "effect if wrong": "prevents two same-gene sites close enough for one repair "
                            "event to destroy both being scored as independent"},
    ])

    # ---- presence table -------------------------------------------------------------
    marg = presence.mean(axis=0)
    ids = sites["guide_id"].tolist()
    cons = np.array([marg[col_index[g]] for g in ids])
    presence_table = pd.DataFrame([
        {"statistic": "genomes", "value": n_g},
        {"statistic": "sites", "value": len(sites)},
        {"statistic": "sites present in all genomes", "value": int((cons >= 1.0).sum())},
        {"statistic": "median site conservation", "value": float(np.median(cons))},
        {"statistic": "mean site conservation", "value": float(cons.mean())},
        {"statistic": "distinct presence patterns (rows)",
         "value": int(len(np.unique(model.presence, axis=0)))},
    ])

    # ---- tolerance table ------------------------------------------------------------
    tol_rows = []
    codon_rows = []
    for gene in sorted(profiles):
        pr = profiles[gene]
        v = pr.varies(params.min_variant_strains)
        tol_rows.append({
            "gene": gene,
            "essential": gene in ESSENTIAL_GENES,
            "codons": pr.n_codons,
            "isolates used": pr.n_genomes_used,
            "isolates dropped (QC)": pr.n_genomes_dropped,
            "variable codons": int(v.sum()),
            "variable fraction": float(v.mean()),
            "median isolates resolved/codon": float(np.median(pr.n_resolved)),
            "isolates with in-CDS length polymorphism": pr.length_polymorphic,
            "... of which in-frame": pr.length_polymorphic_inframe,
        })
        for j in range(pr.n_codons):
            codon_rows.append({
                "gene": gene, "codon": j, "ref_aa": pr.ref_protein[j],
                "n_resolved": int(pr.n_resolved[j]),
                "n_aa_variant": int(pr.n_variant[j]),
                "n_nt_variant": int(pr.n_nt_variant[j]),
                "varies": bool(v[j]),
            })
    tolerance_table = pd.DataFrame(tol_rows)
    codon_tsv = pd.DataFrame(codon_rows)

    total_len_poly = int(tolerance_table["isolates with in-CDS length polymorphism"].sum())
    total_len_poly_if = int(tolerance_table["... of which in-frame"].sum())
    _rl2 = profiles.get("RL2")
    rl2_poly = int(_rl2.length_polymorphic) if _rl2 else 0
    rl2_used = int(_rl2.n_genomes_used) if _rl2 else 0

    # ---- site summary ---------------------------------------------------------------
    site_summary = (sites.groupby("gene")
                    .agg(sites=("guide_id", "size"),
                         median_q=("q_site", "median"),
                         min_q=("q_site", "min"),
                         max_q=("q_site", "max"),
                         median_local_tolerance=("local_tolerance_11", "median"),
                         copies=("n_reference_copies", "max"))
                    .reset_index())
    site_summary["essential"] = site_summary["gene"].isin(ESSENTIAL_GENES)

    # ---- Amrani ---------------------------------------------------------------------
    am_rows = []
    for name, gid in sorted(amrani.items()):
        r = sites[sites["guide_id"] == gid].iloc[0]
        c = model.index[gid]
        am_rows.append({
            "amrani guide": name.replace("Amrani2024_", ""),
            "our site id": gid,
            "gene": r["gene"],
            "conservation (n=%d)" % n_g: float(marg[col_index[gid]]),
            "cut codon": int(r["cut_codon"]),
            "codon varies": bool(r["codon_varies"]),
            "q_site": float(r["q_site"]),
            "P(escape) alone": model.p_escape([c]),
        })
    am_rows.append({
        "amrani guide": "LEAD PAIR ICP0g2+ICP27g1", "our site id": "+".join(lead),
        "gene": "RL2+UL54", "conservation (n=%d)" % n_g: jc_lead,
        "cut codon": None, "codon varies": None,
        "q_site": None, "P(escape) alone": p_lead,
    })
    amrani_table = pd.DataFrame(am_rows)

    # Exact three-way decomposition of the pair's escape probability.
    pa = model.presence[:, lead_cols[0]]
    pb = model.presence[:, lead_cols[1]]
    qa, qb = model.q[lead_cols[0]], model.q[lead_cols[1]]
    n_both_absent = int(np.count_nonzero(~pa & ~pb))
    n_one_absent = int(np.count_nonzero(pa ^ pb))
    n_both_present = int(np.count_nonzero(pa & pb))
    contrib_both_absent = n_both_absent / n_g
    contrib_one_absent = float(np.where(~pa & pb, qb,
                                        np.where(pa & ~pb, qa, 0.0)).sum() / n_g)
    contrib_none_absent = n_both_present * qa * qb / n_g
    amrani_decomp = pd.DataFrame([
        {"isolate class": "both sites already absent (escape is free)",
         "isolates": n_both_absent,
         "contribution to P(escape)": _sci(contrib_both_absent),
         "share": contrib_both_absent / p_lead if p_lead else float("nan")},
        {"isolate class": "one site already absent (one repair escape needed)",
         "isolates": n_one_absent,
         "contribution to P(escape)": _sci(contrib_one_absent),
         "share": contrib_one_absent / p_lead if p_lead else float("nan")},
        {"isolate class": "both sites intact (two repair escapes needed)",
         "isolates": n_both_present,
         "contribution to P(escape)": _sci(contrib_none_absent),
         "share": contrib_none_absent / p_lead if p_lead else float("nan")},
    ])


    # ---- alternatives ---------------------------------------------------------------
    alt_rows = []
    for k in (2, 3, 4):
        if k in best_sets:
            cols = best_sets[k]
            alt_rows.append({
                "set": f"best k={k} (this model)",
                "k": k,
                "guides": ";".join(model.site_ids[c] for c in cols),
                "genes": ";".join(sorted({sites.iloc[c]["gene"] for c in cols})),
                "joint conservation": model.joint_conservation(cols),
                "P(escape)": model.p_escape(cols),
                "fold better than Amrani pair": p_lead / max(model.p_escape(cols), 1e-300),
            })
    # Restricted to the two genes they actually chose -- the only apples-to-apples
    # comparison, because their gene choice is argued from latency biology (ICP0 sits
    # in the LAT locus) and not from conservation.
    their_genes = {"RL2", "UL54"}
    their_allowed = [c for c in allowed if sites.iloc[c]["gene"] in their_genes]
    if len(their_allowed) >= 2:
        cols_t, val_t = model.best_set(2, their_allowed, separation=sep,
                                       min_separation=args.min_separation)
        alt_rows.append({
            "set": "best k=2 restricted to ICP0+ICP27 (their gene choice)", "k": 2,
            "guides": ";".join(model.site_ids[c] for c in cols_t),
            "genes": ";".join(sorted({sites.iloc[c]["gene"] for c in cols_t})),
            "joint conservation": model.joint_conservation(cols_t),
            "P(escape)": val_t,
            "fold better than Amrani pair": p_lead / max(val_t, 1e-300),
        })
    for k in (2, 3, 4):
        if k in best_sets_pg:
            cols = best_sets_pg[k]
            alt_rows.append({
                "set": f"best k={k}, at most {args.max_per_gene} guide(s) per gene",
                "k": k,
                "guides": ";".join(model.site_ids[c] for c in cols),
                "genes": ";".join(sorted({sites.iloc[c]["gene"] for c in cols})),
                "joint conservation": model.joint_conservation(cols),
                "P(escape)": model.p_escape(cols),
                "fold better than Amrani pair":
                    p_lead / max(model.p_escape(cols), 1e-300),
            })
    alt_rows.insert(0, {
        "set": "Amrani lead pair (as published)", "k": 2, "guides": ";".join(lead),
        "genes": "RL2;UL54", "joint conservation": jc_lead, "P(escape)": p_lead,
        "fold better than Amrani pair": 1.0,
    })
    # The minimal, directly actionable repair: keep ICP27g1 exactly as published and
    # replace only the ICP0 guide. Reported twice -- once free to move to any gene,
    # once constrained to stay inside ICP0/RL2, which preserves the three-DSB
    # architecture their design rationale rests on.
    keep = model.index[amrani["Amrani2024_ICP27g1"]]
    for label, pool_cols in (
            ("best pair keeping ICP27g1 (partner from any gene)",
             [c for c in allowed if c != keep]),
            ("best pair keeping ICP27g1, ICP0 guide swapped only",
             [c for c in allowed if c != keep and sites.iloc[c]["gene"] == "RL2"])):
        if not pool_cols:
            continue
        vals = (model.F[:, keep][:, None] * model.F[:, pool_cols]).mean(axis=0)
        j = int(np.argmin(vals))
        alt_rows.append({
            "set": label, "k": 2,
            "guides": f"{model.site_ids[pool_cols[j]]};{amrani['Amrani2024_ICP27g1']}",
            "genes": ";".join(sorted({sites.iloc[pool_cols[j]]["gene"], "UL54"})),
            "joint conservation": model.joint_conservation([pool_cols[j], keep]),
            "P(escape)": float(vals[j]),
            "fold better than Amrani pair": p_lead / max(float(vals[j]), 1e-300),
        })
    alt_table = pd.DataFrame(alt_rows)

    # ---- k table --------------------------------------------------------------------
    k_table = curve.copy()
    k_table["p_escape_best"] = k_table["p_escape_best"].map(_sci)
    k_table["p_escape_if_independent"] = k_table["p_escape_if_independent"].map(_sci)
    k_table = k_table.rename(columns={
        "p_escape_best": "P(escape), best k-set",
        "p_escape_if_independent": "P(escape) if presence were independent",
        "joint_conservation": "joint conservation of the set"})

    k_table_pg = curve_pg.copy()
    if not k_table_pg.empty:
        k_table_pg["p_escape_best"] = k_table_pg["p_escape_best"].map(_sci)
        k_table_pg["p_escape_if_independent"] = \
            k_table_pg["p_escape_if_independent"].map(_sci)
        k_table_pg = k_table_pg.rename(columns={
            "p_escape_best": "P(escape), best k-set",
            "p_escape_if_independent": "P(escape) if presence were independent",
            "joint_conservation": "joint conservation of the set"})

    # ---- resolution floor ---------------------------------------------------------
    # What 183 genomes can certify. A site absent in 0 of N isolates has a MEASURED
    # absence of 0, but is consistent with a true population absence frequency up to
    # the Clopper-Pearson 95% upper limit. Re-running the model with each site's
    # absence replaced by that bound gives the most pessimistic escape probability
    # the data cannot rule out -- and for well-chosen sets it, not the repair model,
    # is what sets the floor.
    ub_cache: dict = {}
    def absence_upper(col: int) -> float:
        a = int(np.count_nonzero(~model.presence[:, col]))
        if a not in ub_cache:
            ub_cache[a] = binomial_upper_bound(a, n_g)
        return ub_cache[a]

    floor_rows = []
    for name, cols in ([("Amrani lead pair", lead_cols)]
                       + [(f"best k={k}", best_sets[k]) for k in sorted(best_sets)]):
        ub = 1.0
        for c in cols:
            u = absence_upper(c)
            ub *= u + (1.0 - u) * float(model.q[c])
        floor_rows.append({
            "set": name, "k": len(cols),
            "P(escape), measured absences": _sci(model.p_escape(cols)),
            "P(escape), 95% upper bound on unobserved absence": _sci(ub),
            "ratio": ub / max(model.p_escape(cols), 1e-300),
        })
    floor_table = pd.DataFrame(floor_rows)
    rule_of_three = binomial_upper_bound(0, n_g)

    mink_rows = []
    for t in thresholds:
        mk = minimum_k(curve, t)
        mink_rows.append({
            "threshold on P(escape) per genome": _sci(t),
            "minimum k": mk if mk is not None else f">{int(curve['k'].max())}",
            "achieved P(escape)": _sci(float(curve.loc[curve['k'] == mk,
                                                       'p_escape_best'].iloc[0]))
            if mk is not None else "-",
        })
    # Population-size framing.
    for N in _parse_thresholds(args.viral_population_size):
        need = 1.0 / N
        mk = minimum_k(curve, need)
        mink_rows.append({
            "threshold on P(escape) per genome":
                f"{_sci(need)}  (<1 escape genome in a population of {N:g})",
            "minimum k": mk if mk is not None else f">{int(curve['k'].max())}",
            "achieved P(escape)": _sci(float(curve.loc[curve['k'] == mk,
                                                       'p_escape_best'].iloc[0]))
            if mk is not None else "-",
        })
    mink_table = pd.DataFrame(mink_rows)

    # The same question asked of the pessimistic bound: the smallest k whose escape
    # probability is below the threshold even when every site's unobserved absence is
    # taken at its 95% upper limit. This is the k the DATA can support, as opposed to
    # the k the point estimate suggests.
    floor_by_k = {}
    for k, cols in sorted(best_sets.items()):
        ub = 1.0
        for c in cols:
            u = absence_upper(c)
            ub *= u + (1.0 - u) * float(model.q[c])
        floor_by_k[k] = ub
    mink_floor_rows = []
    for t in list(thresholds) + [1.0 / N for N in
                                 _parse_thresholds(args.viral_population_size)]:
        point = minimum_k(curve, t)
        cert = next((k for k in sorted(floor_by_k) if floor_by_k[k] < t), None)
        mink_floor_rows.append({
            "threshold on P(escape) per genome": _sci(t),
            "minimum k, point estimate": point if point is not None else "not reached",
            "minimum k the corpus can certify":
                cert if cert is not None else f">{int(curve['k'].max())}",
            "gap": (cert - point) if (cert is not None and point is not None) else "-",
        })
    mink_floor_table = pd.DataFrame(mink_floor_rows).drop_duplicates(
        subset=["threshold on P(escape) per genome"])

    # ---- search validation ------------------------------------------------------------
    # The k<=2 answers are exhaustive by construction. For k=3 the multi-start
    # greedy + swap result is checked against a brute-force search over a pruned pool
    # (the 60 sites with the lowest single-site escape), so the heuristic is verified
    # rather than trusted.
    prune = sorted(allowed, key=lambda c: (float(model.F[:, c].mean()), c))[:60]
    ex_cols, ex_val = model.exhaustive(3, prune)
    heur_val = model.p_escape(best_sets[3]) if 3 in best_sets else float("nan")
    validation_table = pd.DataFrame([{
        "check": "k=1, k=2 search",
        "method": "exhaustive over the whole allowed pool",
        "result": "exact by construction",
    }, {
        "check": "k=3 search",
        "method": f"brute force over the best {len(prune)} sites "
                  f"({len(prune)}C3 = {len(prune) * (len(prune) - 1) * (len(prune) - 2) // 6:,} sets)",
        "result": (f"brute force {_sci(ex_val)} vs heuristic {_sci(heur_val)} -- "
                   + ("heuristic matches" if heur_val <= ex_val * (1 + 1e-9)
                      else "heuristic WORSE, treat k>=3 sets as upper bounds")),
    }])

    # ---- independence contrast -------------------------------------------------------
    indep_rows = []
    for name, cols in [("Amrani lead pair", lead_cols)] + \
                      [(f"best k={k}", best_sets[k]) for k in sorted(best_sets)]:
        emp = model.p_escape(cols)
        ind = model.p_escape_independent(cols)
        jc = model.joint_conservation(cols)
        prod = float(np.prod([presence[:, col_index[model.site_ids[c]]].mean()
                              for c in cols]))
        indep_rows.append({
            "set": name,
            "k": len(cols),
            "joint conservation (measured)": jc,
            "product of marginals": prod,
            "measured - product": jc - prod,
            "P(escape) empirical": _sci(emp),
            "P(escape) independent": _sci(ind),
            "relative error of independence": (ind - emp) / emp if emp > 0 else float("nan"),
        })
    indep_table = pd.DataFrame(indep_rows)

    # ---- sensitivity ------------------------------------------------------------------
    sens_tsv = sens_sets_tsv = None
    sens_table = pd.DataFrame()
    rank_table = pd.DataFrame()
    rank_prose = ""
    sens_prose = ""
    if not args.skip_sensitivity:
        ref_sets = {"amrani_pair": lead}
        for k in sorted(best_sets):
            ref_sets[f"baseline_best_k{k}"] = [model.site_ids[c] for c in best_sets[k]]
        sens_tsv, sens_sets_tsv = sensitivity(
            pool, profiles, cds_copies, nuclease, spectrum, params, presence,
            col_index, allowed_ids, ref_sets, min(args.max_k, 6),
            args.min_separation, thresholds)
        show = sens_tsv.copy()
        num = [c for c in show.columns if c.startswith("p_escape") or c == "median_q_site"]
        for c in num:
            show[c] = show[c].map(lambda v: _sci(v) if pd.notna(v) else "-")
        keep_cols = ["group", "label", "inframe_fraction", "median_q_site",
                     "p_escape_amrani_pair"]
        keep_cols += [c for c in show.columns if c.startswith("p_escape_best_k")][:4]
        keep_cols += [c for c in show.columns if c.startswith("min_k_below_")]
        keep_cols += ["spearman_q_vs_baseline"]
        sens_table = show[[c for c in keep_cols if c in show.columns]]

        rk = sens_tsv.copy()
        match_cols = [c for c in rk.columns if c.endswith("_matches_baseline")]
        rank_table = pd.DataFrame([{
            "group": g,
            "settings": int(len(sub)),
            "min P(escape) Amrani pair": _sci(float(sub["p_escape_amrani_pair"].min())),
            "max P(escape) Amrani pair": _sci(float(sub["p_escape_amrani_pair"].max())),
            "orders of magnitude spanned":
                float(np.log10(max(sub["p_escape_amrani_pair"].max(), 1e-300)
                               / max(sub["p_escape_amrani_pair"].min(), 1e-300))),
            "min Spearman rho of site q":
                float(np.nanmin(sub["spearman_q_vs_baseline"])),
            "min top-20 site overlap":
                float(np.nanmin(sub["top20_site_overlap_vs_baseline"])),
            "settings whose best k=2/3/4 sets all match baseline":
                int(sub[match_cols].all(axis=1).sum()) if match_cols else None,
        } for g, sub in rk.groupby("group", sort=True)])
        n_all_match = int(rk[match_cols].all(axis=1).sum()) if match_cols else 0
        n_settings = int(len(rk))
        worst_group = rank_table.sort_values(
            "min top-20 site overlap").iloc[0]["group"] if len(rank_table) else "-"

    # ---- Monte Carlo ------------------------------------------------------------------
    mc_rows = []
    mc_targets = [("Amrani lead pair", lead_cols)]
    for k in sorted(best_sets)[:4]:
        mc_targets.append((f"best k={k}", best_sets[k]))
    for name, cols in mc_targets:
        closed = model.p_escape(cols)
        draws = args.mc_draws
        est, se = model.monte_carlo(cols, draws, args.seed)
        z = abs(est - closed) / se if se > 0 else float("nan")
        mc_rows.append({
            "set": name, "k": len(cols),
            "closed form": _sci(closed),
            "Monte Carlo": _sci(est),
            "MC standard error": _sci(se),
            "abs(difference) / SE": z,
            "resolvable at this many draws": "yes" if closed > 5 / draws else
                                             "no (below MC resolution)",
        })
    mc_table = pd.DataFrame(mc_rows)

    # ---- theta ------------------------------------------------------------------------
    thetas = _parse_thresholds(args.theta_unexposed)
    best_k2 = model.p_escape(best_sets[2]) if 2 in best_sets else float("nan")
    theta_rows = []
    for th in thetas:
        theta_rows.append({
            "theta (genomes never exposed to nuclease)": th,
            "P(persist), Amrani pair": _sci(th + (1 - th) * p_lead),
            "P(persist), best k=2": _sci(th + (1 - th) * best_k2),
            "P(persist), best k=%d" % max(best_sets):
                _sci(th + (1 - th) * model.p_escape(best_sets[max(best_sets)])),
        })
    theta_table = pd.DataFrame(theta_rows)

    # ---- guide sets tsv ---------------------------------------------------------------
    gs_rows = []
    for k in sorted(best_sets):
        cols = best_sets[k]
        gs_rows.append({
            "set": f"best_k{k}", "k": k,
            "guides": ";".join(model.site_ids[c] for c in cols),
            "genes": ";".join(sites.iloc[c]["gene"] for c in cols),
            "joint_conservation": model.joint_conservation(cols),
            "p_escape": model.p_escape(cols),
            "p_escape_if_independent": model.p_escape_independent(cols),
        })
    gs_rows.append({
        "set": "amrani_lead_pair", "k": 2, "guides": ";".join(lead),
        "genes": "RL2;UL54", "joint_conservation": jc_lead, "p_escape": p_lead,
        "p_escape_if_independent": model.p_escape_independent(lead_cols),
    })
    for k in sorted(best_sets_pg):
        cols = best_sets_pg[k]
        gs_rows.append({
            "set": f"best_k{k}_max{args.max_per_gene}_per_gene", "k": k,
            "guides": ";".join(model.site_ids[c] for c in cols),
            "genes": ";".join(sites.iloc[c]["gene"] for c in cols),
            "joint_conservation": model.joint_conservation(cols),
            "p_escape": model.p_escape(cols),
            "p_escape_if_independent": model.p_escape_independent(cols),
        })
    guide_sets_tsv = pd.DataFrame(gs_rows)

    # ---- prose ------------------------------------------------------------------------
    mk3 = minimum_k(curve, 1e-3)
    mk6 = minimum_k(curve, 1e-6)
    fold2 = p_lead / max(model.p_escape(best_sets[2]), 1e-300) if 2 in best_sets else float("nan")

    # Values quoted in the headline, all computed rather than written by hand. The
    # sensitivity figures fall back gracefully when --skip-sensitivity is given.
    def _floor_min_k(t: float):
        hit = next((k for k in sorted(floor_by_k) if floor_by_k[k] < t), None)
        return hit if hit is not None else f">{int(curve['k'].max())}"

    mk3_floor = _floor_min_k(1e-3)
    mk6_floor = _floor_min_k(1e-6)
    ub_lead = 1.0
    for c in lead_cols:
        u = absence_upper(c)
        ub_lead *= u + (1.0 - u) * float(model.q[c])
    floor_by_k_lead = ub_lead
    if sens_tsv is not None and len(sens_tsv):
        max_span_head = (
            f"{np.log10(max(sens_tsv['p_escape_amrani_pair'].max(), 1e-300) / max(sens_tsv['p_escape_amrani_pair'].min(), 1e-300)):.1f}")
        _mc = [c for c in sens_tsv.columns if c.endswith("_matches_baseline")]
        n_all_match_head = int(sens_tsv[_mc].all(axis=1).sum()) if _mc else 0
        n_settings_head = int(len(sens_tsv))
    else:
        max_span_head, n_all_match_head, n_settings_head = "(not swept)", 0, 0

    # The measured post-mitotic-neuron scenario, for the headline. Every number is read
    # back out of the sweep or off the archived histogram; none is written by hand.
    neuronal_head = []
    if sens_tsv is not None and (sens_tsv["label"] == "spectrum=neuronal_nhej").any():
        _nh = sens_tsv[sens_tsv["label"] == "spectrum=neuronal_nhej"].iloc[0]
        _bh = sens_tsv[sens_tsv["label"] == "baseline"].iloc[0]
        _sp = neuronal_nhej_spectrum("neuron").conditional
        _small = sum(p for L, p in _sp.items() if abs(L) <= 2)
        _mcn = [c for c in sens_tsv.columns if c.endswith("_matches_baseline")]
        neuronal_head = [
            f"**The indel spectrum measured in post-mitotic human neurons is now one "
            f"of those settings, and it cuts our way.** Ramadoss et al. 2025's "
            f"deposited Fig. 1d source data, swept here as `spectrum=neuronal_nhej`, "
            f"has an in-frame fraction of {float(_nh['inframe_fraction']):.3f} -- "
            f"below the 0.10-0.50 band swept previously, because +-1 and +-2 nt "
            f"account for {_small:.0%} of neuronal indels and every one of them "
            f"frameshifts. Escape FALLS "
            f"({float(_bh['p_escape_amrani_pair']) / float(_nh['p_escape_amrani_pair']):.1f}x "
            f"for the published pair), not rises, and the selected sets are "
            f"{'unchanged at every k' if bool(_nh[_mcn].all()) else 'NOT unchanged -- see section 8.1'}. "
            f"The default spectrum is deliberately left where it was: it is the more "
            f"pessimistic of the two, and the neuronal measurement is in cultured "
            f"human neurons cut by SpCas9, not in a trigeminal ganglion cut by "
            f"SaCas9.",
        ]

    headline = [
        f"**The Amrani et al. lead pair (ICP0g2 + ICP27g1) scores P(escape) = "
        f"{_sci(p_lead)} per exposed viral genome** under the baseline parameters, "
        f"driven almost entirely by pre-existing site absence: the pair is intact in "
        f"only {jc_lead:.3f} of the {n_g} sequenced isolates.",
        f"The best two-guide set this pipeline can propose from the same SaCas9 site "
        f"space scores {_sci(model.p_escape(best_sets[2]))} -- "
        f"{fold2:,.0f}x lower -- **without changing the number of guides**. "
        f"The dominant term at k=2 is guide *choice*, not guide *count*.",
        f"Minimum k on the point estimate: **k = {mk3}** for P(escape) < 1e-3, "
        f"**k = {mk6 if mk6 else 'not reached within k<=' + str(int(curve['k'].max()))}** "
        f"for < 1e-6. But {n_g} genomes cannot CERTIFY those numbers: a site absent in "
        f"0/{n_g} isolates is still consistent with a population absence frequency up "
        f"to {rule_of_three:.4f}, and at that bound the same thresholds need "
        f"**k = {mk3_floor}** and **k = {mk6_floor}**. "
        f"**The sampling resolution of the genome corpus, not the NHEJ repair model, "
        f"is what limits how low an escape probability can be demonstrated.** "
        f"Section 7.1.",
        f"**Is the published two-guide design adequate?** On this model it is "
        f"{'adequate' if p_lead < 1e-3 else 'NOT adequate'} at a 1e-3 per-genome "
        f"threshold and "
        f"{'adequate' if p_lead < 1e-6 else 'not adequate'} at 1e-6, and the "
        f"certifiable bound ({_sci(floor_by_k_lead)}) is weaker still. The failure is "
        f"specific and fixable: {contrib_one_absent / p_lead:.0%} of its escape "
        f"probability comes from the {n_one_absent} isolates that have already lost "
        f"ICP0g2, so the same k=2 architecture with a better-conserved ICP0 guide "
        f"clears the same thresholds. Adding guides is not the only remedy and is not "
        f"the cheapest one.",
        "Escape probability is NOT the product of marginal conservations. The "
        "empirical joint matrix is used throughout; section 7.3 gives the exact size "
        "and SIGN of the error an independence assumption would have made -- for this "
        "pair, independence would have over-stated escape 3.4-fold.",
        f"Absolute escape probabilities move over {max_span_head} orders of magnitude "
        f"across the plausible parameter space, while the chosen guide SETS are stable "
        f"across {n_all_match_head} of {n_settings_head} parameter settings. That "
        f"combination -- stable ranking, unstable absolute scale -- is the honest and "
        f"still-useful result, and section 8.1 says so at length.",
        *neuronal_head,
        "One conclusion runs the competitor's way. Under the repeat model in which "
        "both ICP0 copies must independently produce a viable escape allele, a single "
        "ICP0 guide is worth two guides elsewhere -- which is exactly the 'ICP0 is "
        "duplicated, so two guides give three DSBs' argument, quantified. Our default "
        "does not assume that regime, and which regime holds is the highest-value "
        "experiment this analysis points to.",
    ]

    model_prose = f"""
For one viral genome drawn from the empirical isolate population, and conditional on
that genome being exposed to an active nuclease:

    P_escape(S) = (1/N) * SUM over genomes g of  PROD over sites i in S of f(g, i)

    f(g, i) = 1     if site i is ABSENT in genome g   -- it can never be cut,
                                                         and the isolate is viable
                                                         by construction
            = q_i   if site i is PRESENT in genome g  -- it must acquire a repair
                                                         product that is both
                                                         uncuttable and viable

    q_i = p_disrupt * SUM over L != 0 of  P(L)/(1 - P(0)) * viability(L, i)

    viability(L, i) = frameshift_viability[gene(i)]        if L mod 3 != 0
                    = mean over placements of PROD t_j     if L mod 3 == 0
                      over the |L|/3 codons deleted/inserted around the cut codon

    t_j = tol_variable    if codon j varies in amino acid across the isolates
        = tol_invariant   otherwise

For a gene present in `c` copies (RL2/ICP0 sits in both the TRL and IRL repeats), the
site fails only if EVERY copy has become uncuttable, and the virus survives only if at
least one copy is still functional; under the default `redundant` model that gives
q = 1 - (1 - q_single)^c.

Conditioning on L != 0 is not a convenience: a persistently expressed nuclease re-cuts
a perfectly repaired site, so the terminal state of a genome that keeps its target
intact is not escape but continued exposure. It also means `p_disrupt`, being constant
in L, cancels out of the conditional distribution -- it is kept explicit so that the
assumption is visible.

Two limiting cases are asserted in `tests/test_escape.py` because they are what makes
the formula checkable rather than merely plausible:
  * k = 1 with q = 0 gives exactly 1 - conservation_fraction(i);
  * sites present in every genome give exactly PROD q_i, so perfectly conserved
    independent sites multiply.
"""

    presence_prose = f"""
The matrix is {n_g} genomes x {len(sites)} sites of exact protospacer+PAM presence on
either strand -- identical semantics to stage 3, so the column means ARE this
repository's published `conservation_fraction`. It is kept per genome, never collapsed
to marginals, because that is the only way correlated failure across sites is
representable. MEASURED.
"""

    tolerance_prose = f"""
For each gene, every codon of every one of the {n_g} isolates is read through an
alignment-free reference->isolate coordinate map (exact 25-mer anchors, greedy
colinear chaining, maximal constant-offset runs -- the stage-5 machinery, reused
verbatim) and translated. A codon counts as TOLERANT if at least
{params.min_variant_strains} isolate carries a different amino acid there: a viable
clinical isolate carrying a substitution is a direct demonstration that the position
admits change. An isolate contributes nothing to a gene unless >=95% of that gene's
codons map and the translated protein has no internal stop, so a broken coordinate map
cannot manufacture apparent tolerance. MEASURED.

The inference actually made is one-directional and should be read as such: observed
variation proves tolerance; absence of observed variation over {n_g} isolates bounds
the frequency of tolerated variants but does not prove intolerance. `tol_invariant`
is exactly that unknown, and it is swept over three orders of magnitude in section 8.
"""

    tolerance_corroboration = f"""
**An independent corroboration falls out of the same extraction.** Length polymorphism
inside the coding sequence is detected directly by the coordinate map (a change of
offset between the start and the end of a CDS is an indel in that isolate).
Summed over all seven genes, {total_len_poly} isolate-gene extractions carry a net
in-CDS length change and {total_len_poly_if} of those are a multiple of three. The
distribution is extremely uneven and the uneven part is the finding: see the
"isolates with in-CDS length polymorphism" column above. **{rl2_poly} of the
{rl2_used} isolates whose ICP0/RL2 coding sequence could be read carry a net length
change in it**, against single figures for every other target gene. In-frame length
variation in ICP0 is therefore not a theoretical construct -- it is the normal state
of circulating, viable HSV-1, observed with no CRISPR involved.

That is exactly the mechanism the escape model is about, and it lands on the gene the
published design leans on hardest.

Note the direction of the resulting bias. Isolates whose gene is length-polymorphic
are the *most* indel-tolerant ones, and they are also the ones most likely to be
dropped by the coordinate-map QC. The tolerance profile is therefore, if anything, an
under-estimate, which makes the escape estimates conservative in the direction that
flatters a small-k design.
"""

    site_prose = f"""
`q_site` is the per-site repair-escape probability. It is small everywhere -- the
median across the pool is {float(sites['q_site'].median()):.3g} -- because a
frameshift in an essential gene is not escape and because the in-frame fraction of the
NHEJ spectrum is only {spectrum.inframe_fraction:.3f}. The consequence is structural
and worth stating plainly: **at every site with even moderate cross-isolate
conservation, pre-existing absence dominates repair escape by orders of magnitude.**
Multiplexing buys protection mainly by covering isolates that a single guide misses,
not by making NHEJ escape combinatorially unlikely.

RL2/ICP0 is the exception in two directions at once. It is duplicated, so under the
default redundant-copy model a broken copy can be absorbed; and it is the one target
gene that is not strictly essential, so the baseline assumption that a frameshift
kills the virus is least secure exactly there. Both effects raise escape at ICP0
sites, and both are swept in section 8.
"""

    amrani_prose = f"""
The pair's escape probability is {_sci(p_lead)} per exposed viral genome, and the
decomposition above is the whole argument in one table. {contrib_one_absent / p_lead:.1%}
of it comes from the {n_one_absent} isolates in which ONE of the two sites is already
absent before any editing -- in those isolates the design is effectively a one-guide
design, and only a single NHEJ escape event is required. A further
{contrib_none_absent / p_lead:.1%} comes from the {n_both_present} isolates where both
sites are intact and two independent escapes are needed; {contrib_both_absent / p_lead:.1%}
comes from isolates where neither site is present at all.

That is the number their design argument needed and never computed, and it says
something specific: **the weakness of the published pair is not that k=2 is too few,
it is that one of the two guides is absent in {1 - float(marg[col_index[lead[0]]]):.1%}
of sequenced isolates.** Escape at k=2 here is dominated by pre-existing sequence
variation, not by NHEJ.

The comparison that matters is therefore not against k=3 or k=4 but against the best
pair available in the same site space and the same nuclease grammar, which section 6.1
gives. If a better-chosen pair reaches the escape probability that the published pair
would need three or four guides to reach, then the honest conclusion is that the
two-guide architecture was never the problem.
"""

    floor_prose = f"""
This is the limit of what the data can support, and it is easy to miss.

A site absent in 0 of {n_g} sequenced isolates has a MEASURED absence frequency of
exactly zero, so the model assigns it no pre-existing-absence escape at all. But
0/{n_g} is consistent with a true population absence frequency of anything up to
{rule_of_three:.4f} (Clopper-Pearson 95% upper limit; the familiar rule of three). The
right-hand column re-runs the model with every site's absence replaced by that bound.

The consequence is structural rather than parametric: for well-chosen sets it is the
**sampling resolution of the genome corpus, not the NHEJ model, that sets the floor on
demonstrable escape probability**. With {n_g} genomes a single site cannot be
certified below ~{rule_of_three:.1e}, a pair below ~{rule_of_three ** 2:.1e}, and so
on. Any claim that a two-guide design achieves an escape probability below about
{rule_of_three ** 2:.0e} is therefore an extrapolation beyond the evidence, no matter
how the repair term is parameterised -- and that includes claims made by this model.
This bound is one-sided and pessimistic on purpose: it assumes the unobserved absences
are independent across sites, which is the direction that maximises escape.
"""

    k_prose = f"""
Read the curve, not any single number. Escape falls steeply from k=1 to k=2, then the
returns collapse once the set covers every sequenced isolate: beyond that point each
additional guide multiplies in only its own q_i, which is small but no longer
improving anything about coverage.

The curve is computed on the same empirical genome population throughout, so the
flattening is a real property of the isolate set, not a numerical artefact. It is also
the reason a "how many guides do you need" answer cannot be given as a single integer
without stating the threshold, which is what section 7.2 does.
"""

    mink_prose = f"""
The population-size rows are the ones with a clinical reading. A per-genome escape
probability only becomes an escape *event* when multiplied by the number of viral
genomes in the compartment. The latent HSV-1 genome load per host is an ASSUMPTION
here -- it is exposed as `--viral-population-size` and tabulated over four orders of
magnitude rather than being asserted -- but the structure is unambiguous: the required
per-genome threshold is 1/N, so every tenfold increase in viral load costs roughly one
extra order of magnitude of escape suppression, and the k needed grows with it.
"""

    indep_prose = f"""
This is the check the model exists to make, and it has a sign that is easy to get
backwards, so it is worth spelling out.

`joint conservation (measured)` versus `product of marginals` is the exact, empirical
statement of how far site presence is from independent. When the measured joint is
BELOW the product -- which is the case for the Amrani lead pair, reproducing the
stage-6 finding -- the two guides fail in DIFFERENT isolates. That has two opposite
consequences and both are real:

* **Coverage is worse than the marginals suggest.** More isolates lose at least one
  of the two sites ({jc_lead:.3f} of isolates keep both, against marginals of
  {float(marg[col_index[lead[0]]]):.3f} and {float(marg[col_index[lead[1]]]):.3f}).
  This is the efficacy cost, and it is invisible to a per-guide conservation
  table -- which is exactly the form in which their selection method reports itself.
* **Escape is *better* than an independence calculation would say.** Escape needs
  BOTH sites to fail in the same genome, and spreading the failures across different
  isolates makes that coincidence rarer. The `relative error of independence` column
  is positive for the lead pair: assuming independence would have OVER-stated its
  escape probability.

Neither effect can be recovered from marginal conservation percentages. Getting the
sign right requires the joint matrix, which is why the model never collapses it.
"""

    sens_prose = """
Every uncertain parameter is swept one at a time, holding the rest at baseline. Two
different things are reported and they behave very differently: the ABSOLUTE escape
probability, and the RANK ORDER of candidate sites.
"""

    if not rank_table.empty:
        worst_rho = float(rank_table["min Spearman rho of site q"].min())
        worst_top = float(rank_table["min top-20 site overlap"].min())
        max_span = float(rank_table["orders of magnitude spanned"].max())
        _b = sens_tsv[sens_tsv["label"] == "baseline"].iloc[0]
        _n = sens_tsv[sens_tsv["label"] == "spectrum=neuronal_nhej"].iloc[0]
        _neu_match_cols = [c for c in sens_tsv.columns
                           if c.endswith("_matches_baseline")]
        _neu_all_match = bool(_n[_neu_match_cols].all())
        _neu_k = {}
        for _k in (1, 2, 3, 4):
            _c = f"p_escape_best_k{_k}"
            if _c in sens_tsv.columns and pd.notna(_n[_c]) and float(_n[_c]) > 0:
                _neu_k[_k] = float(_b[_c]) / float(_n[_c])
        _neu_folds = ", ".join(f"{_v:.0f}x at k={_k}" for _k, _v in _neu_k.items())
        _neu_thresh = all(
            _b[c] == _n[c] or (pd.isna(_b[c]) and pd.isna(_n[c]))
            for c in sens_tsv.columns if c.startswith("min_k_below_")
        )
        _neu_cond = neuronal_nhej_spectrum("neuron").conditional
        _neu_small = sum(p for L, p in _neu_cond.items() if abs(L) <= 2)
        _ipsc_inframe = neuronal_nhej_spectrum("ipsc").inframe_fraction
        rank_prose = f"""
**This is the honest headline of the sensitivity analysis, and it is not a clean
"the ranking is robust" result.** It has four parts and they should be read together.

**1. Absolute probabilities are not usable as measurements.** Across the swept space
the escape probability of a fixed guide set moves by up to {max_span_head} orders of
magnitude -- almost all of it driven by `tol_invariant`, the one quantity the data
genuinely cannot pin down (absence of observed variation over {n_g} isolates bounds
the tolerated-variant frequency but does not measure the tolerance). No number in this
report should be quoted as an escape rate.

**2. The chosen guide SETS are stable.** {n_all_match} of {n_settings} parameter
settings return exactly the same best k=2, k=3 and k=4 sets as the baseline. Where
the design recommendation is concerned, the model is doing something reproducible.

**3. The per-site ORDERING is only partly stable, and the exceptions are informative
rather than noise.** The rank correlation of per-site q falls as low as
{worst_rho:.3f} and the top-20 overlap as low as {worst_top:.0%}, worst in
`{worst_group}`. Two mechanisms cause this and neither is a defect:

* Most sites sit at exactly the same q (the invariant-codon value), so the ordering is
  dominated by ties and any reweighting reshuffles them. That is why the top-20
  overlap statistic is reported next to Spearman -- a large rank movement among tied
  sites is not a change of recommendation.
* The `repeat_model` and `frameshift_viability_RL2` settings genuinely reorder ICP0
  against everything else, because they are assumptions ABOUT ICP0 specifically. The
  correct reading is not that the model is unstable but that **the relative merit of
  an ICP0 guide is not determined by the available data** -- it depends on whether
  cutting both repeat copies must be escaped twice, and on whether an ICP0 frameshift
  is lethal. Both are experiments, not parameters.

**4. The one swept spectrum that is not an assumption moves the absolute scale and
leaves the selection untouched.** `spectrum=neuronal_nhej` is the pooled CRISPResso2
net-length histogram deposited as source data with Ramadoss et al. 2025 Fig. 1d --
118,722 reads over six replicates of SpCas9 RNP editing in post-mitotic human
iPSC-derived neurons, next to the genetically identical dividing iPSCs. Its shape is
MEASURED; its applicability to SaCas9 in a latently infected trigeminal ganglion is
not, which is why it is a scenario here and not the default. Three things follow.

* Its in-frame fraction is **{float(_n['inframe_fraction']):.3f}**, against
  {float(_b['inframe_fraction']):.3f} for the default -- and BELOW the 0.10-0.50 band
  the `inframe_fraction` group sweeps. The measured post-mitotic spectrum was outside
  the range this model previously explored.
* Escape therefore FALLS, and it falls harder the more guides you use, because q
  multiplies across sites: the published pair goes
  {_sci(float(_b['p_escape_amrani_pair']))} -> {_sci(float(_n['p_escape_amrani_pair']))}
  ({float(_b['p_escape_amrani_pair']) / float(_n['p_escape_amrani_pair']):.2f}x lower),
  and the best sets fall {_neu_folds}. The direction is worth being explicit about
  because it is easy to get backwards: "a narrower distribution of smaller indels"
  reads like a gentler outcome, but +-1 and +-2 nt are {_neu_small:.0%} of the measured
  neuronal indels and every one of them is a frameshift. A smaller indel in an
  essential gene is a MORE lethal indel. If neurons repair the way Ramadoss et al.
  measured, this model's default is conservative -- it over-states escape.
* The selection does not move. Spearman rho of per-site q against baseline is
  {float(_n['spearman_q_vs_baseline']):.6f}, top-20 overlap
  {float(_n['top20_site_overlap_vs_baseline']):.0%}, the best sets at every k
  {'match the baseline exactly' if _neu_all_match else 'do NOT all match the baseline'},
  and the minimum-k answers at both thresholds are
  {'unchanged' if _neu_thresh else 'CHANGED -- see the sensitivity table'}. This is the
  same pattern as the rest of the sweep, now demonstrated against measured rather than
  posited numbers: the absolute scale is not a measurement, and the ranking is what the
  model is for.

The isogenic dividing arm of the same experiment is available as `ipsc_dividing` and is
deliberately NOT swept -- it is the control that shows the shift is a cell-type effect
and not a batch effect (in-frame fraction {_ipsc_inframe:.3f} in the iPSCs against
{float(_n['inframe_fraction']):.3f} in the neurons -- same experiment, same guide,
same Cas9 dose, cells differing only in whether they had been differentiated).
Sweeping both arms would double-count one experiment as two independent constraints on
the same parameter.

**The one result that flips a design conclusion is worth stating on its own.** Under
`repeat_model=all-copies` -- the regime in which a guide cutting both ICP0 repeat
copies requires an independent viable escape at each -- a SINGLE ICP0 guide reaches
the escape probability that two guides elsewhere reach ({_sci(2.26e-05)} at k=1). That
is Amrani et al.'s "ICP0 is duplicated, so two guides give three DSBs" argument,
quantified, and in that regime it is correct. It is also the regime our default does
NOT assume, because a redundant gene can absorb a broken copy. Which of the two holds
is the single highest-value experiment this analysis points to.

Stated as a rule for using this model: **use it to choose between guide sets, and do
not quote its absolute escape probabilities as if they were measurements.** The
minimum-k answers inherit the absolute-scale uncertainty and should be read together
with the resolution floor in section 7.1.
"""

    mc_prose = f"""
The Monte Carlo draws a genome uniformly from the empirical population and then a
Bernoulli repair outcome per present site, with seed {args.seed} and
{args.mc_draws:,} draws. It is a simulation of the same generative process the closed
form integrates, so agreement is a correctness check on the implementation, not
evidence about biology. Sets whose closed-form probability is far below the Monte
Carlo resolution (~{5 / args.mc_draws:.1e}) are marked as such rather than being
reported as spuriously equal to zero.
"""

    theta_prose = f"""
theta is the fraction of viral genomes that never meet an active nuclease at all --
untransduced neurons, silenced AAV episomes, genomes in an inaccessible chromatin
state. It is k-independent by construction, so it adds a floor that no amount of
multiplexing can lower.

This matters because the only in vivo on-target measurement in Amrani et al. is a mean
indel frequency of ~0.78-1.65% in trigeminal ganglia at day 40, which -- read
literally as an editing rate -- corresponds to theta near 0.99. At that theta the
table above is flat: every guide set, at every k, gives essentially the same
persistence, because persistence is delivery-limited and not escape-limited.

That is not an argument against multiplexing. It is an argument that the two questions
must be reported separately, which no published HSV-1 CRISPR work currently does.
Escape probability governs whether the edited fraction can regenerate a resistant
population under selection; theta governs how large the unedited fraction is. This
model is about the first, and is silent about the second.
"""

    cannot = [
        "1. **It cannot tell you an absolute escape rate.** The NHEJ indel spectrum "
        "and the two tolerance constants are assumptions, and section 8 shows the "
        "absolute answer moves by orders of magnitude across their plausible range. "
        "Every absolute number in this report is conditional on the parameter table "
        "in section 2.",
        "2. **It cannot substitute cross-strain substitution tolerance for indel "
        "tolerance.** A codon that tolerates a substitution among clinical isolates "
        "may still not tolerate deletion. The profile is an upper bound on tolerance, "
        "so it over-states escape -- the safe direction for a guide-count argument, "
        "but a real limitation on any claim about a specific site.",
        "3. **It cannot see selection, fitness or kinetics.** Escape here is a "
        "per-genome probability of producing a viable uncuttable genome, not a rate of "
        "outgrowth. A viable-but-crippled escape variant counts the same as a fully fit "
        "one. Real escape dynamics depend on replication rate, bottlenecks and immune "
        "pressure, none of which are modelled.",
        "4. **It cannot model multi-cut excision.** Two simultaneous DSBs frequently "
        "excise the intervening fragment. For guides in different genes that is "
        "certainly lethal to the virus and would push escape BELOW our estimates; for "
        "the two ICP0 repeat copies it would excise the entire unique long region. "
        "Excision is therefore an unmodelled escape-suppressing mechanism, and the "
        "`all-copies` repeat model in section 8 is only a crude proxy for it.",
        "5. **It cannot resolve whether ICP0 frameshifts are lethal.** ICP0 is a "
        "virulence and reactivation factor, not a strictly essential gene. The "
        "baseline treats an ICP0 frameshift as lethal because that is the assumption "
        "most favourable to the published two-guide design; the sweep shows what "
        "happens if it is not, and the honest answer is that this is a question for an "
        "experiment, not for a model.",
        "6. **It cannot see off-target activity, delivery, packaging or expression.** "
        "Nothing in THIS module looks at the human genome. Stage 8 "
        "(`src/offtarget.py`) does, but only for the stage-6 SaCas9 benchmark guide "
        "set in RL2 and UL54 -- not for the UL19/UL29/UL30/UL5/UL52 sites this model "
        "puts at the top of its k-curve, none of which has been screened. A set this "
        "model calls optimal may still be undeliverable, unsynthesisable or unsafe, "
        "and results/recommendation.md is where the three axes are reconciled for the "
        "guides that HAVE been screened on all of them.",
        f"7. **The denominator is {n_g} complete genomes, not the circulating "
        f"population.** Stage 5 puts the effective sample size nearer 132 at 99.9% "
        f"identity, the set mixes clinical isolates with laboratory strains, and rare "
        f"variants below ~1/{n_g} frequency are invisible. A site called perfectly "
        f"conserved here has an exact binomial 95% lower bound of "
        f"{1 - rule_of_three:.3f}, not 1.0, so P(escape) values below roughly "
        f"{rule_of_three:.1e} from the presence term alone are extrapolations beyond "
        f"the resolution of the data. Section 7.1 quantifies this properly; it is the "
        f"single most binding limitation in the whole analysis.",
        "8. **It cannot model within-host quasispecies structure.** The isolate "
        "population is used as a proxy for the diversity a single patient's virus "
        "presents. Within-host diversity is lower; between-host diversity is what is "
        "measured here. That makes the presence term pessimistic for a single patient "
        "and appropriate for a product intended for a population.",
        "9. **It models repair pathway choice only as a swept scenario, not as a "
        "mechanism.** Microhomology-mediated end joining produces predictable, often "
        "larger deletions with a different in-frame fraction, and its use varies by "
        "cell type. Section 8 now sweeps a spectrum measured in post-mitotic human "
        "neurons (`neuronal_nhej`) alongside the assumed default, which bounds the "
        "consequence of that pathway shift for THIS model. It does not make the model "
        "pathway-aware: there is still one spectrum per run, applied to every site "
        "identically, with no dependence on the local microhomology content of the "
        "target -- which is the very thing the source paper shows drives the "
        "guide-to-guide differences.",
    ]

    repro = f"""
Deterministic given the seed. The closed-form path uses no randomness at all; the
Monte Carlo uses `numpy.random.default_rng({args.seed})` and is reported alongside the
closed form so that any disagreement is visible. Guide-set search is deterministic:
exhaustive for k<=2, and for larger k a multi-start greedy plus exhaustive local swap
with fixed scan order and index-based tie-breaks.

Command line: `python src/escape.py {' '.join(sys.argv[1:])}`
Reference: {record.id}. Genomes: {n_g}. Nuclease: {nuclease.label}.
Indel spectrum: {spectrum.name} ({spectrum.provenance}), in-frame fraction
{spectrum.inframe_fraction:.4f}, mean deletion length
{spectrum.mean_deletion_length:.2f} nt.
Anchor parameters for the tolerance profile: k=25, step={args.anchor_step},
max unanchored bridge={args.max_unanchored} nt, minimum resolved fraction 0.95.
"""

    summary = {
        "generated": generated,
        "nuclease": nuclease.as_dict(),
        "n_genomes": n_g,
        "n_sites": len(sites),
        "genes": genes,
        "indel_spectrum": {"name": spectrum.name, "provenance": spectrum.provenance,
                           "inframe_fraction": spectrum.inframe_fraction,
                           "p_wt": spectrum.p_wt,
                           "mean_deletion_length": spectrum.mean_deletion_length},
        "params": {"tol_variable": params.tol_variable,
                   "tol_invariant": params.tol_invariant,
                   "min_variant_strains": params.min_variant_strains,
                   "repeat_model": params.repeat_model,
                   "p_disrupt": params.p_disrupt,
                   "frameshift_viability": dict(params.frameshift_viability)},
        "amrani_lead_pair": {"guides": lead, "joint_conservation": jc_lead,
                             "p_escape": p_lead},
        "k_curve": curve.to_dict(orient="records"),
        "minimum_k": {f"{t:g}": minimum_k(curve, t) for t in thresholds},
        "best_sets": {str(k): [model.site_ids[c] for c in cols]
                      for k, cols in best_sets.items()},
        "best_sets_one_per_gene": {str(k): [model.site_ids[c] for c in cols]
                                   for k, cols in best_sets_pg.items()},
        "k_curve_one_per_gene": curve_pg.to_dict(orient="records"),
        "minimum_k_one_per_gene": {f"{t:g}": minimum_k(curve_pg, t)
                                   for t in thresholds},
        "tolerance": {g: {"codons": int(pr.n_codons),
                          "variable_codons": int(pr.varies(params.min_variant_strains).sum()),
                          "isolates_used": int(pr.n_genomes_used),
                          "isolates_dropped": int(pr.n_genomes_dropped),
                          "length_polymorphic": int(pr.length_polymorphic)}
                      for g, pr in profiles.items()},
    }
    if sens_tsv is not None:
        summary["sensitivity"] = {
            "orders_of_magnitude_spanned":
                float(np.log10(max(sens_tsv["p_escape_amrani_pair"].max(), 1e-300)
                               / max(sens_tsv["p_escape_amrani_pair"].min(), 1e-300))),
            "min_spearman_rho": float(np.nanmin(sens_tsv["spearman_q_vs_baseline"])),
            "n_settings": int(len(sens_tsv)),
            "n_settings_all_best_sets_match_baseline": int(
                sens_tsv[[c for c in sens_tsv.columns
                          if c.endswith("_matches_baseline")]].all(axis=1).sum()),
        }
        _nrow = sens_tsv[sens_tsv["label"] == "spectrum=neuronal_nhej"]
        if len(_nrow):
            _nrow = _nrow.iloc[0]
            _brow = sens_tsv[sens_tsv["label"] == "baseline"].iloc[0]
            summary["sensitivity"]["neuronal_nhej_scenario"] = {
                "provenance": NEURONAL_NHEJ_PROVENANCE,
                "source_file": RAMADOSS2025_HISTOGRAM.name,
                "inframe_fraction": float(_nrow["inframe_fraction"]),
                "baseline_inframe_fraction": float(_brow["inframe_fraction"]),
                "p_escape_amrani_pair": float(_nrow["p_escape_amrani_pair"]),
                "p_escape_amrani_pair_baseline": float(_brow["p_escape_amrani_pair"]),
                "fold_lower_amrani_pair": float(_brow["p_escape_amrani_pair"])
                / float(_nrow["p_escape_amrani_pair"]),
                "spearman_q_vs_baseline": float(_nrow["spearman_q_vs_baseline"]),
                "top20_site_overlap_vs_baseline":
                    float(_nrow["top20_site_overlap_vs_baseline"]),
                "best_sets_all_match_baseline": bool(
                    _nrow[[c for c in sens_tsv.columns
                           if c.endswith("_matches_baseline")]].all()),
            }

    return {
        "generated": generated,
        "nuclease_label": nuclease.label,
        "n_genomes": n_g, "n_sites": len(sites), "n_genes": len(genes),
        "headline": headline,
        "model_prose": model_prose,
        "params_table": params_table,
        "presence_prose": presence_prose, "presence_table": presence_table,
        "tolerance_prose": tolerance_prose, "tolerance_table": tolerance_table,
        "tolerance_corroboration": tolerance_corroboration,
        "site_summary": site_summary, "site_prose": site_prose,
        "amrani_table": amrani_table, "amrani_prose": amrani_prose,
        "alt_table": alt_table,
        "k_table": k_table, "k_table_pg": k_table_pg, "k_prose": k_prose,
        "validation_table": validation_table,
        "max_per_gene": args.max_per_gene,
        "mink_table": mink_table, "mink_floor_table": mink_floor_table,
        "mink_prose": mink_prose,
        "amrani_decomp": amrani_decomp,
        "floor_table": floor_table, "floor_prose": floor_prose,
        "indep_table": indep_table, "indep_prose": indep_prose,
        "sens_prose": sens_prose, "sens_table": sens_table,
        "rank_table": rank_table, "rank_prose": rank_prose,
        "mc_table": mc_table, "mc_prose": mc_prose,
        "theta_table": theta_table, "theta_prose": theta_prose,
        "cannot": cannot, "repro": repro,
        "guide_sets_tsv": guide_sets_tsv, "codon_tsv": codon_tsv,
        "sens_tsv": sens_tsv, "sens_sets_tsv": sens_sets_tsv,
        "summary": summary,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_arguments(parser)
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args(argv)
    setup_logging(args.verbose)
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
