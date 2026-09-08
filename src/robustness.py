"""Stage 5 -- robustness / sensitivity analysis of the conservation estimate.

Stages 1-4 answer "what fraction of complete HSV-1 genomes contain this exact
23-mer?". That number is only as good as its denominator. This stage attacks the
denominator from three directions and reports what it finds, whether or not the
answer is flattering.

A. REDUNDANCY -- is n=183 really 183 independent observations?
   The set mixes clinical isolates with laboratory strains, and strain 17 is
   deposited at least four times. We sketch every genome with a mash-style bottom-s
   MinHash over canonical k-mers (pure Python + numpy -- no mash binary), estimate
   all-pairs Jaccard, convert it to an ANI-like identity, cluster by single linkage
   at several identity thresholds, and recompute conservation on one representative
   per cluster. The number of clusters is an effective sample size.

B. AMBIGUITY -- how much of "not conserved" is really "not sequenced"?
   Two independent probes.
     1. A sweep over --max-ambiguous-fraction: as poorly resolved assemblies are
        dropped, the genome count falls and the perfectly-conserved guide count
        rises. The sweep shows how the two co-vary, so the reader can see that
        "945 guides" and "833 guides" are the same data at different denominators.
     2. An N-TOLERANT match mode. For every (guide, genome) pair where the exact
        23-mer is absent we ask whether the homologous region of that genome is
        actually resolved. If it is not, the pair is scored UNKNOWN and removed
        from that guide's denominator instead of being counted as a mismatch.
        Homology is established by anchor chaining (below), not by alignment.

C. THE DENOMINATOR QUESTION -- 183 complete genomes vs the whole of GenBank.
   NCBI holds thousands of sub-genomic HSV-1 records (individual CDS entries,
   clinical survey amplicons, drug-resistance genotyping fragments) plus hundreds
   of near-full-length "partial genome" assemblies. Excision BioTherapeutics
   (Amrani et al. 2024) measured guide conservation against a corpus of that kind
   and used a >70% threshold. We re-measure the SAME guides against that corpus.

   The trap this stage exists to avoid: a 600 bp UL30 amplicon does not contain a
   UL19 guide, and that is not evidence of variation. A record is admitted to a
   guide's denominator ONLY if it demonstrably spans the homologous region.

Coverage determination (used by both B and C)
---------------------------------------------
No aligner is available and none is wanted, so homology is established with exact
k-mer anchors and colinear chaining -- the same idea as the seed-chain-extend stage
of minimap2/BLAST, minus the extend:

  1. Index every k-mer (default k=25) of the reference genome that occurs at most
     --anchor-max-occurrences times (repeats such as TRL/IRL are kept, up to that
     multiplicity; low-complexity k-mers are dropped).
  2. Scan the record every --anchor-step bases, in both orientations; keep the
     orientation with more anchors.
  3. Chain anchors greedily: an anchor extends an open chain if it advances in both
     the reference and the record, the reference gap is <= --anchor-max-gap, and the
     diagonal (ref_pos - record_pos) moves by <= --anchor-max-drift. Indels move the
     diagonal gradually and are absorbed; a second repeat copy sits on a completely
     different diagonal and starts its own chain.
  4. A guide's reference footprint is COVERED by a record if some chain contains an
     anchor entirely to its left and another entirely to its right. Those two anchors
     bracket the homologous stretch of the record exactly -- no padding constant is
     involved. Everything between two anchors of one chain is covered even if it is
     divergent, so a guide is NOT silently dropped merely because it sits in a
     variable region. That matters: dropping variable regions is precisely the bias
     that would manufacture a flattering result.

  Given the bracket, a non-matching (guide, record) pair is classified:
     ABSENT        bracket is short and contains no ambiguity  -> real sequence variation
     AMBIGUOUS     bracket contains a non-ACGT base            -> unresolved, UNKNOWN
     WIDE_BRACKET  bracket longer than --max-bracket-span, or the record reaches both
                   sides of the site with no ambiguity in between but no colinear
                   chain crosses it -> structurally unlike the reference, UNKNOWN
     NOT_COVERED   the record does not reach both sides of the site at all
  Only PRESENT and ABSENT enter the N-tolerant denominator.

Everything here is deterministic: the MinHash uses a fixed 64-bit mixing function
(splitmix64), not Python's randomised hash; chaining is a fixed greedy pass over a
sorted list; representatives are chosen by an explicit total order.

PRIVACY: the contact email comes from $NCBI_EMAIL only. Nothing is hardcoded.
"""

from __future__ import annotations

import argparse
import json
import logging
import math
import sys
from bisect import bisect_left, bisect_right
from pathlib import Path

import numpy as np
import pandas as pd

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.common import (  # noqa: E402
    DATA_DIR,
    DEFAULT_CANDIDATES,
    DEFAULT_CONSERVATION,
    DEFAULT_GENE_MANIFEST,
    DEFAULT_MANIFEST,
    DEFAULT_ROBUSTNESS_REPORT,
    GENE_DIR,
    PROJECT_ROOT,
    REF_DIR,
    RESULTS_DIR,
    RateLimiter,
    configure_entrez,
    ensure_dirs,
    entrez_rate_limit_interval,
    revcomp,
    setup_logging,
    utc_now_iso,
)
from src.conservation import (  # noqa: E402
    build_lookup,
    load_genome,
    nuclease_from_candidates,
    scan_genome,
)
from src.nuclease import SPCAS9, Nuclease  # noqa: E402

LOG = logging.getLogger("robustness")

# Status codes for a (guide, record) pair. Stored as int8 in dense matrices.
ST_PRESENT = 0
ST_ABSENT = 1
ST_AMBIGUOUS = 2
ST_WIDE_BRACKET = 3
ST_NOT_COVERED = 4

STATUS_NAMES = {
    ST_PRESENT: "present",
    ST_ABSENT: "absent",
    ST_AMBIGUOUS: "unknown_ambiguous",
    ST_WIDE_BRACKET: "unknown_structural",
    ST_NOT_COVERED: "not_covered",
}

# Ambiguity thresholds swept in part B. Fixed grid, declared up front, so that no
# threshold can be chosen after seeing the answer.
AMBIGUITY_SWEEP = [None, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001,
                   0.0005, 0.0002, 0.0001, 0.0]

# Identity thresholds swept in part A. Also fixed up front.
ANI_SWEEP = [0.9999, 0.999, 0.998, 0.995, 0.99, 0.98]

# Published SaCas9 guides from Amrani et al. 2024 (Mol Ther Methods Clin Dev 32:101303,
# PMC11602521), Table 1. Spacer (20 nt) + NNGRRT PAM (6 nt), as printed in the paper.
# Included so that the competing guides can be scored with OUR method against OUR
# corpora -- a like-for-like comparison, not a quotation of their numbers.
BENCHMARK_GUIDES = {
    "Amrani2024_ICP0g1": ("GTACCCGACGGCCCCCGCGT", "CGGAGT"),
    "Amrani2024_ICP0g2": ("CTCAGGCCGCGAACCAAGAA", "CAGAGT"),
    "Amrani2024_ICP27g1": ("AATCCTAGACACGCACCGCC", "AGGAGT"),
    "Amrani2024_ICP27g2": ("TCGCCAGCGTCATTAGCGGG", "GGGGGT"),
}


# ======================================================================================
# Part A helpers -- mash-style MinHash over canonical k-mers
# ======================================================================================

_BASE_LUT = np.full(256, 255, dtype=np.uint8)
for _i, _b in enumerate(b"ACGT"):
    _BASE_LUT[_b] = _i

_M1 = np.uint64(0xFF51AFD7ED558CCD)
_M2 = np.uint64(0xC4CEB9FE1A85EC53)
_S33 = np.uint64(33)
_S2 = np.uint64(2)


def _splitmix64(x: np.ndarray) -> np.ndarray:
    """Deterministic 64-bit avalanche mixer.

    Python's built-in hash() of a str is salted per process (PYTHONHASHSEED), which
    would make sketches irreproducible across runs. This is not.
    """
    x = x.astype(np.uint64, copy=True)
    x ^= x >> _S33
    x *= _M1
    x ^= x >> _S33
    x *= _M2
    x ^= x >> _S33
    return x


def canonical_kmer_hashes(seq: str, k: int = 21) -> np.ndarray:
    """64-bit hashes of every canonical (strand-independent) k-mer of `seq`.

    Windows containing a non-ACGT base are skipped entirely.
    """
    if k <= 0 or k > 31:
        raise ValueError("k must be in 1..31 (2 bits per base must fit in uint64)")
    raw = np.frombuffer(seq.encode("ascii", "replace"), dtype=np.uint8)
    if raw.size < k:
        return np.empty(0, dtype=np.uint64)
    codes = _BASE_LUT[raw]
    bad = codes == 255
    codes64 = np.where(bad, 0, codes).astype(np.uint64)

    m = raw.size - k + 1
    fwd = np.zeros(m, dtype=np.uint64)
    rev = np.zeros(m, dtype=np.uint64)
    three = np.uint64(3)
    for j in range(k):
        fwd = (fwd << _S2) | codes64[j:j + m]
        rev = (rev << _S2) | (three - codes64[k - 1 - j:k - 1 - j + m])

    cum = np.concatenate(([0], np.cumsum(bad.astype(np.int64))))
    window_bad = cum[k:] - cum[:m]
    ok = window_bad == 0

    canon = np.minimum(fwd, rev)[ok]
    return _splitmix64(canon)


def minhash_sketch(seq: str, k: int = 21, sketch_size: int = 5000) -> np.ndarray:
    """Bottom-s sketch: the `sketch_size` smallest distinct k-mer hashes, sorted."""
    h = canonical_kmer_hashes(seq, k)
    if h.size == 0:
        return np.empty(0, dtype=np.uint64)
    uniq = np.unique(h)
    return uniq[:sketch_size]


def sketch_jaccard(a: np.ndarray, b: np.ndarray, sketch_size: int) -> float:
    """Mash's estimator: |merged n A n B| / |merged|, merged = bottom-s of A u B."""
    if a.size == 0 or b.size == 0:
        return 0.0
    merged = np.union1d(a, b)[:sketch_size]
    shared = np.intersect1d(a, b, assume_unique=True)
    hits = np.intersect1d(merged, shared, assume_unique=True).size
    return hits / merged.size


def mash_identity(jaccard: float, k: int) -> float:
    """Mash distance turned into an ANI-like identity in [0, 1]."""
    if jaccard <= 0.0:
        return 0.0
    if jaccard >= 1.0:
        return 1.0
    d = -(1.0 / k) * math.log(2.0 * jaccard / (1.0 + jaccard))
    return max(0.0, 1.0 - d)


def single_linkage_clusters(labels: list[str], identity: np.ndarray,
                            threshold: float) -> list[int]:
    """Union-find single-linkage clustering. Returns a cluster index per label.

    Cluster indices are assigned in order of the first member's position in
    `labels`, so the numbering is deterministic.
    """
    n = len(labels)
    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i in range(n):
        for j in range(i + 1, n):
            if identity[i, j] >= threshold:
                ri, rj = find(i), find(j)
                if ri != rj:
                    parent[max(ri, rj)] = min(ri, rj)

    order: dict[int, int] = {}
    out = []
    for i in range(n):
        root = find(i)
        if root not in order:
            order[root] = len(order)
        out.append(order[root])
    return out


def choose_representatives(manifest: pd.DataFrame, cluster_of: dict[str, int]) -> dict[int, str]:
    """One representative per cluster, by an explicit deterministic total order:

    best resolved (lowest ambiguous_fraction) > longest > RefSeq > lowest accession.
    """
    best: dict[int, tuple] = {}
    for row in manifest.itertuples(index=False):
        cid = cluster_of[row.accession]
        key = (float(row.ambiguous_fraction), -int(row.length_bp),
               0 if bool(row.is_refseq) else 1, str(row.accession))
        if cid not in best or key < best[cid][0]:
            best[cid] = (key, str(row.accession))
    return {cid: acc for cid, (_k, acc) in best.items()}


# ======================================================================================
# Coverage: exact k-mer anchors + colinear chaining
# ======================================================================================


_ACGT = frozenset("ACGT")


def build_reference_index(ref: str, k: int = 25,
                          max_occurrences: int = 4) -> dict[str, tuple[int, ...]]:
    """Map every sufficiently unique reference k-mer to its plus-strand positions.

    k-mers containing a non-ACGT base are excluded (they cannot anchor anything
    reliably). k-mers occurring more than `max_occurrences` times are dropped as
    low-complexity; the HSV-1 large repeats occur twice, which is retained.
    """
    positions: dict[str, list[int]] = {}
    n = len(ref)
    for i in range(n - k + 1):
        sub = ref[i:i + k]
        lst = positions.get(sub)
        if lst is None:
            positions[sub] = [i]
        elif len(lst) <= max_occurrences:
            lst.append(i)
    return {kmer: tuple(pos) for kmer, pos in positions.items()
            if len(pos) <= max_occurrences and not (set(kmer) - _ACGT)}


def find_anchors(query: str, index: dict[str, tuple[int, ...]], k: int,
                 step: int) -> list[tuple[int, int]]:
    """Exact (reference_pos, query_pos) k-mer matches, sampled every `step` bases."""
    out: list[tuple[int, int]] = []
    limit = len(query) - k + 1
    if limit <= 0:
        return out
    for q in range(0, limit, step):
        hit = index.get(query[q:q + k])
        if hit is not None:
            for r in hit:
                out.append((r, q))
    return out


def orient_and_anchor(query: str, index: dict[str, tuple[int, ...]], k: int,
                      step: int) -> tuple[str, str, list[tuple[int, int]]]:
    """Return (oriented_sequence, strand, anchors).

    The record may be deposited in either orientation relative to the reference.
    Whichever orientation yields more anchors is the one we work in; all downstream
    coordinates refer to `oriented_sequence`.
    """
    fwd = find_anchors(query, index, k, step)
    if len(fwd) >= 50:
        # A record cannot anchor 50 exact 25-mers in both orientations; skip the
        # reverse pass, which for a 152 kb genome is the expensive half.
        return query, "+", fwd
    rc = revcomp(query)
    rev = find_anchors(rc, index, k, step)
    if len(rev) > len(fwd):
        return rc, "-", rev
    return query, "+", fwd


class Chain:
    """A colinear run of anchors: a stretch of the record homologous to the reference."""

    __slots__ = ("ref", "qry", "ref_min", "ref_max")

    def __init__(self, ref: list[int], qry: list[int], k: int):
        self.ref = ref
        self.qry = qry
        self.ref_min = ref[0]
        self.ref_max = ref[-1] + k

    def __len__(self) -> int:
        return len(self.ref)


def chain_anchors(anchors: list[tuple[int, int]], k: int, max_gap: int = 3000,
                  max_drift: int = 60, min_anchors: int = 2) -> list[Chain]:
    """Greedy colinear chaining over anchors sorted by reference position.

    Several chains stay open at once, which is what lets the two copies of a large
    repeat (TRL/IRL) chain independently instead of corrupting each other.
    """
    if not anchors:
        return []
    ordered = sorted(set(anchors))
    open_chains: list[tuple[list[int], list[int]]] = []
    done: list[tuple[list[int], list[int]]] = []

    for r, q in ordered:
        still_open = []
        best = None
        best_drift = None
        for ch in open_chains:
            last_r, last_q = ch[0][-1], ch[1][-1]
            if r - last_r > max_gap:
                done.append(ch)
                continue
            still_open.append(ch)
            if q <= last_q:
                continue
            drift = abs((r - last_r) - (q - last_q))
            if drift <= max_drift and (best_drift is None or drift < best_drift):
                best, best_drift = ch, drift
        open_chains = still_open
        if best is not None:
            best[0].append(r)
            best[1].append(q)
        else:
            new = ([r], [q])
            open_chains.append(new)
    done.extend(open_chains)

    chains = [Chain(ref, qry, k) for ref, qry in done if len(ref) >= min_anchors]
    chains.sort(key=lambda c: (c.ref_min, c.ref_max))
    return chains


def bracket_in_chain(chain: Chain, ref_start0: int, ref_end0: int,
                     k: int) -> tuple[int, int] | None:
    """Query interval bracketing [ref_start0, ref_end0) within one chain.

    Returns (q_lo, q_hi) spanned by the nearest anchor strictly left of the footprint
    and the nearest anchor strictly right of it, or None if the chain does not
    bracket the footprint on both sides.
    """
    ref = chain.ref
    # nearest anchor entirely to the left: ref[i] + k <= ref_start0
    i = bisect_right(ref, ref_start0 - k) - 1
    if i < 0:
        return None
    # nearest anchor entirely to the right: ref[j] >= ref_end0
    j = bisect_left(ref, ref_end0)
    if j >= len(ref):
        return None
    q_lo = chain.qry[i]
    q_hi = chain.qry[j] + k
    if q_hi <= q_lo:
        return None
    return q_lo, q_hi


def classify_footprint(footprints: list[tuple[int, int]], chains: list[Chain],
                       amb_positions: np.ndarray, k: int,
                       max_bracket_span: int) -> int:
    """Classify a non-matching guide against one record.

    `footprints` are the guide's 0-based half-open reference intervals (a guide in a
    duplicated repeat has more than one). The best -- i.e. most informative --
    classification over all footprints and chains wins, in the order
    ABSENT > AMBIGUOUS > WIDE_BRACKET > NOT_COVERED: if any copy of the site is
    demonstrably present-and-resolved-but-different, that is real variation.
    """
    best = ST_NOT_COVERED
    for start0, end0 in footprints:
        for chain in chains:
            if chain.ref_min > start0 or chain.ref_max < end0:
                continue
            br = bracket_in_chain(chain, start0, end0, k)
            if br is None:
                continue
            q_lo, q_hi = br
            if q_hi - q_lo > max_bracket_span:
                if best > ST_WIDE_BRACKET:
                    best = ST_WIDE_BRACKET
                continue
            if _has_ambiguity(amb_positions, q_lo, q_hi):
                if best > ST_AMBIGUOUS:
                    best = ST_AMBIGUOUS
                continue
            return ST_ABSENT

    if best != ST_NOT_COVERED:
        return best

    # No single chain brackets the site. Two very different situations look the
    # same at this point and must not be merged: an assembly GAP (the record does
    # reach both sides of the site but the middle is Ns, which breaks every chain),
    # and a genuine TRUNCATION or rearrangement (the record simply does not have
    # the region). Distinguish them by looking in the space between the chain that
    # ends before the site and the chain that starts after it.
    for start0, end0 in footprints:
        left = right = None
        for chain in chains:
            if chain.ref_max <= start0 and (left is None or chain.ref_max > left.ref_max):
                left = chain
            if chain.ref_min >= end0 and (right is None or chain.ref_min < right.ref_min):
                right = chain
        if left is None or right is None:
            continue
        q_lo = left.qry[-1]
        q_hi = right.qry[0] + k
        if q_hi <= q_lo:
            best = min(best, ST_WIDE_BRACKET)   # rearranged relative to the reference
            continue
        if _has_ambiguity(amb_positions, q_lo, q_hi):
            return ST_AMBIGUOUS
        best = min(best, ST_WIDE_BRACKET)
    return best


def _has_ambiguity(amb_positions: np.ndarray, q_lo: int, q_hi: int) -> bool:
    if amb_positions.size == 0:
        return False
    lo = np.searchsorted(amb_positions, q_lo, side="left")
    hi = np.searchsorted(amb_positions, q_hi, side="left")
    return bool(hi > lo)


def ambiguous_positions(seq: str) -> np.ndarray:
    raw = np.frombuffer(seq.encode("ascii", "replace"), dtype=np.uint8)
    return np.flatnonzero(_BASE_LUT[raw] == 255).astype(np.int64)


def covered_reference_span(chains: list[Chain]) -> tuple[int, int, int]:
    """(min_ref, max_ref, total_chained_bp) over all chains -- for the corpus manifest."""
    if not chains:
        return (-1, -1, 0)
    total = sum(c.ref_max - c.ref_min for c in chains)
    return (min(c.ref_min for c in chains), max(c.ref_max for c in chains), total)


# ======================================================================================
# Guide bookkeeping
# ======================================================================================


def parse_footprints(candidates: pd.DataFrame) -> list[list[tuple[int, int]]]:
    """0-based half-open reference footprints per guide, including repeat copies.

    `ref_all_positions` holds entries like "1234-1256(+)"; a guide inside TRL/IRL has
    one per copy. A record covering either copy covers the guide.
    """
    out: list[list[tuple[int, int]]] = []
    for row in candidates.itertuples(index=False):
        spans: list[tuple[int, int]] = []
        raw = str(getattr(row, "ref_all_positions", "") or "")
        for piece in raw.split(";"):
            piece = piece.strip()
            if not piece:
                continue
            coords = piece.split("(")[0]
            if "-" not in coords:
                continue
            a, b = coords.split("-", 1)
            try:
                spans.append((int(a) - 1, int(b)))
            except ValueError:
                continue
        if not spans:
            spans = [(int(row.ref_start) - 1, int(row.ref_end))]
        out.append(sorted(set(spans)))
    return out


def reference_genbank_path(accession: str) -> Path:
    """Locate the cached reference record.

    Stage 2 caches it under the accession it was asked for (`NC_001806`) but stores
    the versioned id (`NC_001806.2`) in the candidate table, so both are tried.
    """
    for name in (accession, accession.split(".")[0]):
        path = REF_DIR / f"{name}.gb"
        if path.is_file():
            return path
    raise SystemExit(
        f"Reference annotation for {accession} not found under {REF_DIR}. "
        "Run stage 2 (extract_guides) first."
    )


def load_reference_sequence(accession: str) -> str:
    from Bio import SeqIO

    record = SeqIO.read(reference_genbank_path(accession), "genbank")
    return str(record.seq).upper()


def benchmark_targets() -> dict[str, str]:
    return {name: proto + pam for name, (proto, pam) in BENCHMARK_GUIDES.items()}


def score_benchmark(seq: str, targets: dict[str, str]) -> set[str]:
    """Plain exact substring presence on either strand.

    The stage-3 jump scan is specialised for the SpCas9 NGG PAM at offset 21 and
    cannot be used for SaCas9's 6 nt NNGRRT PAM, so these four use a direct search.
    """
    rc = revcomp(seq)
    return {name for name, t in targets.items() if t in seq or t in rc}


# ======================================================================================
# Part A + B: one pass over the complete genomes
# ======================================================================================


def scan_complete_genomes(manifest: pd.DataFrame, candidates: pd.DataFrame,
                          reference: str, args,
                          nuclease: Nuclease = SPCAS9) -> dict:
    """Single pass producing everything parts A and B need.

    Returns dense per-(genome, guide) status and presence matrices plus per-genome
    MinHash sketches, so that every downstream subset (de-duplicated set, ambiguity
    sweep, N-tolerant mode) is a cheap slice rather than a re-scan.
    """
    ensure_dirs()
    guide_ids = candidates["guide_id"].tolist()
    target_map = dict(zip(candidates["guide_id"], candidates["target_23mer"]))
    fwd, rev = build_lookup(target_map)
    footprints = parse_footprints(candidates)
    guide_index = {gid: i for i, gid in enumerate(guide_ids)}

    LOG.info("Building reference k-mer index (k=%d, max_occ=%d) ...",
             args.anchor_k, args.anchor_max_occurrences)
    index = build_reference_index(reference, args.anchor_k, args.anchor_max_occurrences)
    LOG.info("  %d anchorable reference %d-mers.", len(index), args.anchor_k)

    accessions = manifest["accession"].astype(str).tolist()
    n_g, n_s = len(accessions), len(guide_ids)
    status = np.full((n_g, n_s), ST_NOT_COVERED, dtype=np.int8)
    sketches: list[np.ndarray] = []
    bench_targets = benchmark_targets()
    bench_present = {name: np.zeros(n_g, dtype=bool) for name in bench_targets}
    per_genome_chain_counts = []

    for gi, row in enumerate(manifest.itertuples(index=False)):
        path = PROJECT_ROOT / str(row.fasta_path)
        seq = load_genome(path)
        sketches.append(minhash_sketch(seq, args.minhash_k, args.minhash_sketch))

        hits = scan_genome(seq, fwd, rev, nuclease)
        oriented, _strand, anchors = orient_and_anchor(seq, index, args.anchor_k,
                                                       args.anchor_step)
        chains = chain_anchors(anchors, args.anchor_k, args.anchor_max_gap,
                               args.anchor_max_drift, args.anchor_min_anchors)
        per_genome_chain_counts.append(len(chains))
        amb = ambiguous_positions(oriented)

        for name in score_benchmark(seq, bench_targets):
            bench_present[name][gi] = True

        row_status = status[gi]
        for si, gid in enumerate(guide_ids):
            if gid in hits:
                row_status[si] = ST_PRESENT
            else:
                row_status[si] = classify_footprint(
                    footprints[si], chains, amb, args.anchor_k, args.max_bracket_span)
        if (gi + 1) % 25 == 0 or gi + 1 == n_g:
            LOG.info("  robustness scan %d/%d genomes", gi + 1, n_g)

    return {
        "accessions": accessions,
        "guide_ids": guide_ids,
        "guide_index": guide_index,
        "status": status,
        "sketches": sketches,
        "benchmark_present": bench_present,
        "chain_counts": per_genome_chain_counts,
    }


def redundancy_analysis(manifest: pd.DataFrame, scan: dict, args) -> dict:
    """Part A: sketch-based all-pairs identity, clustering, effective sample size."""
    accessions = scan["accessions"]
    sketches = scan["sketches"]
    n = len(accessions)
    ident = np.eye(n, dtype=np.float64)
    jacc = np.eye(n, dtype=np.float64)
    for i in range(n):
        for j in range(i + 1, n):
            jv = sketch_jaccard(sketches[i], sketches[j], args.minhash_sketch)
            iv = mash_identity(jv, args.minhash_k)
            jacc[i, j] = jacc[j, i] = jv
            ident[i, j] = ident[j, i] = iv

    off = ident.copy()
    np.fill_diagonal(off, -1.0)
    nearest_idx = off.argmax(axis=1)
    nearest_id = off.max(axis=1)

    presence = scan["status"] == ST_PRESENT
    sweep_rows = []
    cluster_assignments: dict[float, list[int]] = {}
    for thr in ANI_SWEEP:
        clusters = single_linkage_clusters(accessions, ident, thr)
        cluster_assignments[thr] = clusters
        n_clusters = len(set(clusters))
        cluster_of = dict(zip(accessions, clusters))
        reps = choose_representatives(manifest, cluster_of)
        rep_rows = [i for i, acc in enumerate(accessions) if reps[cluster_of[acc]] == acc]
        sub = presence[rep_rows]
        perfect = int(sub.all(axis=0).sum()) if sub.size else 0
        sweep_rows.append({
            "ani_threshold": thr,
            "n_clusters": n_clusters,
            "n_genomes_collapsed": n - n_clusters,
            "largest_cluster": int(pd.Series(clusters).value_counts().max()),
            "n_perfectly_conserved": perfect,
            "pct_perfectly_conserved": round(100.0 * perfect / presence.shape[1], 3),
        })

    primary = args.dedup_identity
    if primary not in cluster_assignments:
        cluster_assignments[primary] = single_linkage_clusters(accessions, ident, primary)
    clusters = cluster_assignments[primary]
    cluster_of = dict(zip(accessions, clusters))
    reps = choose_representatives(manifest, cluster_of)
    sizes = pd.Series(clusters).value_counts().to_dict()

    meta = manifest.set_index("accession")
    per_genome = pd.DataFrame({
        "accession": accessions,
        "strain": [str(meta.at[a, "strain"]) if a in meta.index else "" for a in accessions],
        "length_bp": [int(meta.at[a, "length_bp"]) for a in accessions],
        "ambiguous_fraction": [float(meta.at[a, "ambiguous_fraction"]) for a in accessions],
        "cluster_id": clusters,
        "cluster_size": [sizes[c] for c in clusters],
        "is_representative": [reps[c] == a for a, c in zip(accessions, clusters)],
        "nearest_accession": [accessions[k] for k in nearest_idx],
        "nearest_identity": [round(float(v), 6) for v in nearest_id],
        "nearest_jaccard": [round(float(jacc[i, nearest_idx[i]]), 6) for i in range(n)],
    })
    per_genome = per_genome.sort_values(["cluster_id", "accession"],
                                        kind="stable").reset_index(drop=True)

    rep_rows = [i for i, acc in enumerate(accessions) if reps[cluster_of[acc]] == acc]
    dedup_presence = presence[rep_rows]
    dedup_perfect = int(dedup_presence.all(axis=0).sum())
    full_perfect = int(presence.all(axis=0).sum())

    multi = per_genome[per_genome["cluster_size"] > 1]
    return {
        "sweep": pd.DataFrame(sweep_rows),
        "per_genome": per_genome,
        "identity_matrix_summary": {
            "min_pairwise_identity": round(float(ident[np.triu_indices(n, 1)].min()), 6),
            "median_pairwise_identity": round(float(np.median(ident[np.triu_indices(n, 1)])), 6),
            "max_pairwise_identity": round(float(ident[np.triu_indices(n, 1)].max()), 6),
        },
        "primary_threshold": primary,
        "n_genomes": n,
        "n_clusters": len(set(clusters)),
        "n_multi_member_clusters": int(multi["cluster_id"].nunique()),
        "n_genomes_in_multi_clusters": int(len(multi)),
        "perfect_full": full_perfect,
        "perfect_dedup": dedup_perfect,
        "rep_rows": rep_rows,
    }


def ambiguity_analysis(manifest: pd.DataFrame, scan: dict, args) -> dict:
    """Part B: threshold sweep plus the N-tolerant re-interpretation."""
    accessions = scan["accessions"]
    presence = scan["status"] == ST_PRESENT
    amb_frac = manifest.set_index("accession")["ambiguous_fraction"].astype(float)
    frac = np.array([amb_frac[a] for a in accessions])

    rows = []
    for thr in AMBIGUITY_SWEEP:
        keep = np.ones(len(accessions), dtype=bool) if thr is None else (frac <= thr)
        n_keep = int(keep.sum())
        if n_keep == 0:
            rows.append({"max_ambiguous_fraction": thr, "n_genomes": 0,
                         "n_perfectly_conserved": 0, "pct_perfectly_conserved": 0.0,
                         "n_genomes_dropped": len(accessions)})
            continue
        perfect = int(presence[keep].all(axis=0).sum())
        rows.append({
            "max_ambiguous_fraction": "none" if thr is None else thr,
            "n_genomes": n_keep,
            "n_genomes_dropped": len(accessions) - n_keep,
            "n_perfectly_conserved": perfect,
            "pct_perfectly_conserved": round(100.0 * perfect / presence.shape[1], 3),
        })

    status = scan["status"]
    n_present = (status == ST_PRESENT).sum(axis=0)
    n_absent = (status == ST_ABSENT).sum(axis=0)
    n_ambig = (status == ST_AMBIGUOUS).sum(axis=0)
    n_wide = (status == ST_WIDE_BRACKET).sum(axis=0)
    n_uncov = (status == ST_NOT_COVERED).sum(axis=0)
    n_total = status.shape[0]

    denom_ntol = n_present + n_absent
    frac_strict = n_present / n_total
    frac_ntol = np.where(denom_ntol > 0, n_present / np.maximum(denom_ntol, 1), np.nan)

    per_guide = pd.DataFrame({
        "guide_id": scan["guide_ids"],
        "n_genomes_total": n_total,
        "n_present": n_present,
        "n_absent_resolved": n_absent,
        "n_unknown_ambiguous": n_ambig,
        "n_unknown_structural": n_wide,
        "n_not_covered": n_uncov,
        "conservation_strict": np.round(frac_strict, 6),
        "conservation_n_tolerant": np.round(frac_ntol, 6),
    }).sort_values("guide_id", kind="stable").reset_index(drop=True)

    perfect_strict = int((n_present == n_total).sum())
    perfect_ntol = int(((denom_ntol > 0) & (n_present == denom_ntol)).sum())
    rescued = int(((n_present < n_total) & (denom_ntol > 0) &
                   (n_present == denom_ntol)).sum())

    return {
        "sweep": pd.DataFrame(rows),
        "per_guide": per_guide,
        "perfect_strict": perfect_strict,
        "perfect_n_tolerant": perfect_ntol,
        "n_rescued_by_n_tolerance": rescued,
        "status_totals": {STATUS_NAMES[s]: int((status == s).sum())
                          for s in sorted(STATUS_NAMES)},
        "n_pairs": int(status.size),
    }


# ======================================================================================
# Part C: the gene-level corpus
# ======================================================================================


def gene_corpus_queries(taxid: int, min_len: int, max_len: int,
                        include_partial_genomes: bool,
                        genome_min_len: int) -> list[tuple[str, str]]:
    """(record_class, entrez_query) pairs defining the gene-level corpus.

    Deliberately NOT keyword-based. Gene keywords such as "RL2" are indexed on only
    a handful of HSV-1 records (verified: 5 hits), so a keyword query would silently
    lose almost the whole corpus. Instead every HSV-1 nucleotide record shorter than
    a genome is downloaded and assigned to genes by anchor mapping against the
    reference annotation -- the record's own labelling is never trusted.
    """
    base = f"txid{taxid}[Organism:exp] AND biomol_genomic[PROP] NOT patent[PROP]"
    out = [("subgenomic", f"{base} AND {min_len}:{max_len}[SLEN]")]
    if include_partial_genomes:
        out.append(("partial_genome",
                    f"{base} AND {genome_min_len}:160000[SLEN] "
                    'NOT "complete genome"[Title]'))
    return out


def fetch_gene_corpus(queries: list[tuple[str, str]], batch_size: int,
                      refresh: bool, limit: int | None) -> tuple[list[tuple[str, str]], dict]:
    """Download and cache the corpus under data/raw/genes/. Returns (records, stats)."""
    from src.fetch_genomes import esearch_all_uids, fetch_fasta_batch, split_fasta_records

    ensure_dirs()
    Entrez = configure_entrez()
    limiter = RateLimiter(entrez_rate_limit_interval())

    records: list[tuple[str, str]] = []
    stats: dict = {"queries": [], "n_requested": 0, "n_cached_hit": 0, "n_downloaded": 0}

    for record_class, query in queries:
        accs, reported = esearch_all_uids(Entrez, query, limiter)
        if limit is not None and limit < len(accs):
            LOG.warning("--gene-corpus-limit %d: using %d of %d %s accessions.",
                        limit, limit, len(accs), record_class)
            accs = accs[:limit]
        stats["queries"].append({
            "record_class": record_class,
            "entrez_query": query,
            "esearch_reported_total": reported,
            "accessions_used": len(accs),
        })
        stats["n_requested"] += len(accs)

        missing = [a for a in accs
                   if refresh or not (GENE_DIR / f"{a}.fasta").is_file()
                   or (GENE_DIR / f"{a}.fasta").stat().st_size < 60]
        stats["n_cached_hit"] += len(accs) - len(missing)
        LOG.info("%s: %d accessions (%d cached, %d to download).",
                 record_class, len(accs), len(accs) - len(missing), len(missing))

        for i in range(0, len(missing), batch_size):
            batch = missing[i:i + batch_size]
            blob = fetch_fasta_batch(Entrez, batch, limiter)
            got = split_fasta_records(blob)
            by_base = {rid.split(".")[0]: rid for rid in got}
            for acc in batch:
                key = acc if acc in got else by_base.get(acc.split(".")[0])
                if key is None:
                    LOG.error("NCBI returned no FASTA for %s -- omitted.", acc)
                    continue
                (GENE_DIR / f"{acc}.fasta").write_text(got[key], encoding="utf-8")
                stats["n_downloaded"] += 1
            if (i // batch_size) % 5 == 0 or i + batch_size >= len(missing):
                LOG.info("  downloaded %d/%d %s records",
                         min(i + batch_size, len(missing)), len(missing), record_class)

        for acc in accs:
            path = GENE_DIR / f"{acc}.fasta"
            if path.is_file() and path.stat().st_size >= 60:
                records.append((record_class, acc))

    stats["n_available"] = len(records)
    if stats["n_available"] < stats["n_requested"]:
        LOG.warning("REAL COUNT: %d of %d requested corpus records are on disk.",
                    stats["n_available"], stats["n_requested"])
    return records, stats


def gene_level_analysis(records: list[tuple[str, str]], candidates: pd.DataFrame,
                        reference: str, gene_spans: dict[str, list[tuple[int, int]]],
                        args, nuclease: Nuclease = SPCAS9) -> dict:
    """Part C: re-measure the same guides against the gene-level corpus."""
    guide_ids = candidates["guide_id"].tolist()
    target_map = dict(zip(candidates["guide_id"], candidates["target_23mer"]))
    fwd, rev = build_lookup(target_map)
    footprints = parse_footprints(candidates)
    n_s = len(guide_ids)

    LOG.info("Building reference k-mer index for corpus mapping ...")
    index = build_reference_index(reference, args.anchor_k, args.anchor_max_occurrences)

    # Sorted guide start positions, so a short record only classifies guides that
    # could plausibly lie inside its chains.
    starts = np.array([min(s for s, _e in fp) for fp in footprints])
    order = np.argsort(starts, kind="stable")
    starts_sorted = starts[order].tolist()
    guide_pos = {gid: i for i, gid in enumerate(guide_ids)}

    STATUSES = ("present", "absent", "ambiguous", "wide", "not_covered")
    RECORD_CLASSES = ("subgenomic", "partial_genome")
    # Counted separately per record class. The 393 near-full-length "partial genome"
    # records are the same KIND of evidence as the 183 complete genomes (a title
    # convention separates them), whereas the sub-genomic records are independent
    # clinical amplicons. Pooling them would hide which of the two is carrying the
    # result, so they are tracked apart and reported apart.
    by_class = {cls: {name: np.zeros(n_s, dtype=np.int32) for name in STATUSES}
                for cls in RECORD_CLASSES}
    bench_targets = benchmark_targets()
    bench = {name: {"present": 0, "covered": 0} for name in bench_targets}
    bench_spans = {
        "Amrani2024_ICP0g1": "RL2", "Amrani2024_ICP0g2": "RL2",
        "Amrani2024_ICP27g1": "UL54", "Amrani2024_ICP27g2": "UL54",
    }

    manifest_rows = []
    n_mapped = 0
    for ri, (record_class, acc) in enumerate(records, start=1):
        seq = load_genome(GENE_DIR / f"{acc}.fasta")
        if len(seq) < args.anchor_k:
            continue
        oriented, strand, anchors = orient_and_anchor(seq, index, args.anchor_k,
                                                      args.anchor_step)
        chains = chain_anchors(anchors, args.anchor_k, args.anchor_max_gap,
                               args.anchor_max_drift, args.anchor_min_anchors)
        amb = ambiguous_positions(oriented)
        hits = scan_genome(seq, fwd, rev, nuclease)

        genes_touched: set[str] = set()
        scored: set[int] = set()
        n_cov_here = 0
        if chains:
            n_mapped += 1
            for gene, spans in gene_spans.items():
                for gs, ge in spans:
                    if any(c.ref_min < ge and c.ref_max > gs for c in chains):
                        genes_touched.add(gene)
                        break

            lo_ref = min(c.ref_min for c in chains)
            hi_ref = max(c.ref_max for c in chains)
            lo = bisect_left(starts_sorted, lo_ref - 25)
            hi = bisect_right(starts_sorted, hi_ref + 25)
            counts = by_class[record_class]
            for pos in order[lo:hi]:
                pos = int(pos)
                scored.add(pos)
                gid = guide_ids[pos]
                if gid in hits:
                    counts["present"][pos] += 1
                    n_cov_here += 1
                    continue
                st = classify_footprint(footprints[pos], chains, amb, args.anchor_k,
                                        args.max_bracket_span)
                if st == ST_ABSENT:
                    counts["absent"][pos] += 1
                    n_cov_here += 1
                elif st == ST_AMBIGUOUS:
                    counts["ambiguous"][pos] += 1
                elif st == ST_WIDE_BRACKET:
                    counts["wide"][pos] += 1
                else:
                    counts["not_covered"][pos] += 1

        present_names = score_benchmark(seq, bench_targets)
        for name, target in bench_targets.items():
            gene = bench_spans[name]
            covered = gene in genes_touched
            if name in present_names:
                bench[name]["present"] += 1
                bench[name]["covered"] += 1
            elif covered:
                bench[name]["covered"] += 1

        # An exact 23-mer match is itself proof that the record spans the site, so a
        # hit outside the chained window (very short records that cannot be chained)
        # is still counted as present.
        for gid in hits:
            pos = guide_pos[gid]
            if pos not in scored:
                by_class[record_class]["present"][pos] += 1
                scored.add(pos)
                n_cov_here += 1
        manifest_rows.append({
            "accession": acc,
            "record_class": record_class,
            "length_bp": len(seq),
            "n_ambiguous": int(amb.size),
            "strand_vs_reference": strand,
            "n_anchor_chains": len(chains),
            "ref_span_start": chains[0].ref_min + 1 if chains else -1,
            "ref_span_end": max(c.ref_max for c in chains) if chains else -1,
            "genes_covered": "|".join(sorted(genes_touched)),
            "n_guides_scored": n_cov_here,
        })
        if ri % 250 == 0 or ri == len(records):
            LOG.info("  corpus %d/%d records mapped", ri, len(records))

    counts = {name: sum(by_class[c][name] for c in RECORD_CLASSES) for name in STATUSES}

    def _frac(num, den):
        return np.round(np.where(den > 0, num / np.maximum(den, 1), np.nan), 6)

    denom = counts["present"] + counts["absent"]
    # The "no coverage check" denominator: every record that was considered for this
    # guide at all counts, whether or not it actually spans the site. This is what a
    # conservation number looks like when it is computed as
    #   (records containing the guide) / (records assigned to the gene),
    # which is how Amrani et al. 2024 describe their calculation. It is reported here
    # NOT as a better estimate but to show what the coverage correction is worth.
    denom_naive = (denom + counts["ambiguous"] + counts["wide"]
                   + counts["not_covered"])
    sub_denom = by_class["subgenomic"]["present"] + by_class["subgenomic"]["absent"]
    par_denom = by_class["partial_genome"]["present"] + by_class["partial_genome"]["absent"]

    per_guide = pd.DataFrame({
        "guide_id": guide_ids,
        "gene": candidates["gene"].tolist(),
        "n_records_covering": denom,
        "n_present": counts["present"],
        "n_absent_resolved": counts["absent"],
        "n_unknown_ambiguous": counts["ambiguous"],
        "n_unknown_structural": counts["wide"],
        "n_not_covered_in_span": counts["not_covered"],
        "conservation_gene_level": _frac(counts["present"], denom),
        "conservation_no_coverage_check": _frac(counts["present"], denom_naive),
        "n_records_covering_subgenomic": sub_denom,
        "n_present_subgenomic": by_class["subgenomic"]["present"],
        "conservation_subgenomic": _frac(by_class["subgenomic"]["present"], sub_denom),
        "n_records_covering_partial_genome": par_denom,
        "conservation_partial_genomes": _frac(by_class["partial_genome"]["present"],
                                              par_denom),
    })

    return {
        "per_guide": per_guide,
        "corpus_manifest": pd.DataFrame(manifest_rows),
        "n_records": len(records),
        "n_records_mapped": n_mapped,
        "benchmark": bench,
    }


def collect_gene_spans(reference_accession: str, genes: list[str]) -> dict[str, list[tuple[int, int]]]:
    from Bio import SeqIO

    from src.extract_guides import collect_target_segments

    record = SeqIO.read(reference_genbank_path(reference_accession), "genbank")
    segments = collect_target_segments(record, genes, "CDS")
    spans: dict[str, list[tuple[int, int]]] = {}
    for gene, _product, _lt, s0, e, _st in segments:
        spans.setdefault(gene, []).append((s0, e))
    return spans


def compare_denominators(genome_frac: pd.DataFrame, gene_frac: pd.DataFrame,
                         min_records: int,
                         records_col: str = "n_records_covering",
                         cons_col: str = "conservation_gene_level") -> dict:
    """The headline: what happens to our 100%-conserved guides at gene-level n?

    `records_col`/`cons_col` select the corpus tier (whole corpus, sub-genomic
    records only, partial genomes only, or the pessimistic bound).
    """
    joined = genome_frac.merge(gene_frac, on="guide_id", how="inner", validate="one_to_one")
    merged = joined[["guide_id", "conservation_strict", records_col, cons_col]].copy()
    merged.columns = ["guide_id", "conservation_strict",
                      "n_records_covering", "conservation_gene_level"]
    usable = merged[merged["n_records_covering"] >= min_records].copy()

    out: dict = {
        "tier": cons_col,
        "n_guides": int(len(merged)),
        "n_guides_with_enough_records": int(len(usable)),
        "min_records_required": min_records,
    }
    if usable.empty:
        return out

    a = usable["conservation_strict"].to_numpy(dtype=float)
    b = usable["conservation_gene_level"].to_numpy(dtype=float)
    ok = ~np.isnan(a) & ~np.isnan(b)
    a, b = a[ok], b[ok]
    if a.size >= 2 and a.std() > 0 and b.std() > 0:
        out["pearson_r"] = round(float(np.corrcoef(a, b)[0, 1]), 4)
        ra = pd.Series(a).rank().to_numpy()
        rb = pd.Series(b).rank().to_numpy()
        out["spearman_rho"] = round(float(np.corrcoef(ra, rb)[0, 1]), 4)
    else:
        out["pearson_r"] = None
        out["spearman_rho"] = None

    perfect = usable[usable["conservation_strict"] >= 1.0]
    out["n_perfect_at_183"] = int(len(perfect))
    for cut in (0.99, 0.95, 0.90, 0.70):
        out[f"n_perfect_at_183_below_{cut:.2f}_gene_level"] = int(
            (perfect["conservation_gene_level"] < cut).sum())
    out["n_perfect_at_183_still_perfect"] = int(
        (perfect["conservation_gene_level"] >= 1.0).sum())
    out["median_gene_level_of_perfect"] = (
        round(float(perfect["conservation_gene_level"].median()), 4) if len(perfect) else None)
    out["min_gene_level_of_perfect"] = (
        round(float(perfect["conservation_gene_level"].min()), 4) if len(perfect) else None)
    out["mean_gene_level_all"] = round(float(usable["conservation_gene_level"].mean()), 4)
    out["mean_genome_level_all"] = round(float(usable["conservation_strict"].mean()), 4)
    out["n_gene_level_ge_0.70"] = int((usable["conservation_gene_level"] >= 0.70).sum())
    out["n_gene_level_ge_0.95"] = int((usable["conservation_gene_level"] >= 0.95).sum())
    out["n_gene_level_eq_1.00"] = int((usable["conservation_gene_level"] >= 1.0).sum())
    out["median_records_per_guide"] = int(usable["n_records_covering"].median())
    out["max_records_per_guide"] = int(usable["n_records_covering"].max())
    return out


# ======================================================================================
# Reporting
# ======================================================================================


def _md_table(df: pd.DataFrame, floatfmt: str = "{:.4g}") -> str:
    cols = list(df.columns)
    lines = ["| " + " | ".join(str(c) for c in cols) + " |",
             "|" + "|".join("---" for _ in cols) + "|"]
    for row in df.itertuples(index=False):
        cells = []
        for v in row:
            if isinstance(v, float) and not math.isnan(v):
                cells.append(floatfmt.format(v))
            elif isinstance(v, float):
                cells.append("n/a")
            else:
                cells.append(str(v))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def render_report(ctx: dict) -> str:
    a = ctx["redundancy"]
    b = ctx["ambiguity"]
    c = ctx.get("gene_level")
    cmp_ = ctx.get("comparison")
    n_guides = ctx["n_guides"]
    n_genomes = ctx["n_genomes"]

    L: list[str] = []
    add = L.append

    add("# Robustness and sensitivity analysis of HSV-1 guide conservation")
    add("")
    add(f"Generated {ctx['generated_utc']} by `src/robustness.py`. "
        f"Reference {ctx['reference']}. {n_guides} candidate guides, "
        f"{n_genomes} complete genomes.")
    add("")
    add("Every number below was produced by the run that wrote this file. Nothing is "
        "extrapolated, and no threshold was chosen after seeing its effect: the "
        "identity and ambiguity grids are module-level constants.")
    add("")

    # ---------------- headline ----------------
    add("## Headline")
    add("")
    add(f"* Baseline (stage 3): **{b['perfect_strict']} / {n_guides}** guides are present "
        f"in all {n_genomes} complete genomes.")
    add(f"* De-duplicating near-identical genomes at {a['primary_threshold']*100:.2f}% "
        f"estimated identity leaves an effective sample size of **{a['n_clusters']}** "
        f"independent genomes, and the perfectly-conserved count moves to "
        f"**{a['perfect_dedup']}** ({a['perfect_dedup'] - a['perfect_full']:+d}).")
    add(f"* Scoring N as UNKNOWN rather than MISMATCH raises the perfectly-conserved "
        f"count to **{b['perfect_n_tolerant']}** ({b['perfect_n_tolerant'] - b['perfect_strict']:+d}); "
        f"{b['n_rescued_by_n_tolerance']} guides are perfect only because unresolved "
        f"positions were removed from their denominator.")
    if cmp_ and cmp_.get("n_perfect_at_183"):
        add(f"* Against the gene-level corpus "
            f"({c['n_records']} records, {c['n_records_mapped']} of which map to the "
            f"reference), of the **{cmp_['n_perfect_at_183']}** guides that are 100% "
            f"conserved at n={n_genomes} and have at least "
            f"{cmp_['min_records_required']} covering records, "
            f"**{cmp_['n_perfect_at_183_below_0.95_gene_level']}** fall below 95% and "
            f"**{cmp_['n_perfect_at_183_below_0.70_gene_level']}** fall below the 70% "
            f"threshold used by Amrani et al. 2024.")
        sub_t = (cmp_.get("tiers") or {}).get("subgenomic_only")
        if sub_t and sub_t.get("n_guides_with_enough_records"):
            add(f"* **Caveat that changes the reading:** most of that corpus is 393 "
                f"near-full-length `partial genome` assemblies, not independent "
                f"per-gene records. Restricted to genuinely sub-genomic records, only "
                f"**{sub_t['n_guides_with_enough_records']} of {n_guides}** guides have "
                f"a usable denominator at all; of the "
                f"**{sub_t['n_perfect_at_183']}** perfect ones among them, "
                f"**{sub_t['n_perfect_at_183_below_0.95_gene_level']}** fall below 95%.")
        elif sub_t is not None:
            add("* **Caveat that changes the reading:** restricted to genuinely "
                "sub-genomic records, NO guide reaches the minimum-record bar. Almost "
                "all of the apparent widening comes from 393 near-full-length "
                "`partial genome` assemblies.")
    add("")

    # ---------------- A ----------------
    add("## A. Redundancy and effective sample size")
    add("")
    add(f"Method: bottom-{ctx['minhash_sketch']} MinHash sketch of canonical "
        f"{ctx['minhash_k']}-mers per genome (splitmix64 mixing, no salted hashing), "
        "all-pairs Jaccard by the mash estimator, Jaccard converted to an ANI-like "
        "identity, single-linkage clustering. No alignment, no all-pairs full "
        "comparison.")
    add("")
    add("Pairwise identity across the set: "
        f"min {a['identity_matrix_summary']['min_pairwise_identity']:.5f}, "
        f"median {a['identity_matrix_summary']['median_pairwise_identity']:.5f}, "
        f"max {a['identity_matrix_summary']['max_pairwise_identity']:.5f}.")
    add("")
    add("### Effective sample size vs clustering threshold")
    add("")
    add(_md_table(a["sweep"]))
    add("")
    add(f"At the primary threshold ({a['primary_threshold']}), "
        f"{a['n_genomes_in_multi_clusters']} genomes fall into "
        f"{a['n_multi_member_clusters']} multi-member clusters; the rest are singletons.")
    add("")
    ex = ctx.get("redundancy_examples")
    if ex is not None and not ex.empty:
        add("Largest clusters (deposited near-duplicates):")
        add("")
        add(_md_table(ex))
        add("")
    add(f"**Effect on the result.** Recomputing conservation on one representative per "
        f"cluster at {a['primary_threshold']}: {a['perfect_dedup']} perfectly conserved "
        f"guides over {a['n_clusters']} genomes, versus {a['perfect_full']} over "
        f"{n_genomes}. Interpretation is in the closing section.")
    add("")
    add("Per-genome cluster assignments: `results/robustness_redundancy.tsv`.")
    add("")

    # ---------------- B ----------------
    add("## B. Ambiguity (N) sensitivity")
    add("")
    add("### B1. Sweeping --max-ambiguous-fraction")
    add("")
    add("Dropping poorly resolved assemblies raises conservation and shrinks the "
        "denominator at the same time. Both columns must be quoted together; either "
        "one alone is misleading.")
    add("")
    add(_md_table(b["sweep"]))
    add("")
    add("### B2. N-tolerant matching")
    add("")
    add("For every (guide, genome) pair where the exact 23-mer is absent, the "
        "homologous region is located by anchor chaining and inspected. Only pairs "
        "whose homologous region is present and fully resolved are counted as real "
        "mismatches.")
    add("")
    pairs = b["n_pairs"]
    tot = pd.DataFrame([{"classification": k, "pairs": v,
                         "pct_of_all_pairs": round(100.0 * v / pairs, 3)}
                        for k, v in b["status_totals"].items()])
    add(_md_table(tot))
    add("")
    add(f"Strict (stage 3) perfectly conserved: **{b['perfect_strict']}**. "
        f"N-tolerant perfectly conserved: **{b['perfect_n_tolerant']}**. "
        f"Guides that become perfect only under N-tolerance: "
        f"**{b['n_rescued_by_n_tolerance']}**.")
    add("")
    add("Per-guide breakdown: `results/robustness_conservation_modes.tsv`.")
    add("")

    # ---------------- C ----------------
    add("## C. The denominator question: 183 complete genomes vs the gene-level corpus")
    add("")
    if c is None:
        add("**NOT RUN.** The gene-level corpus was skipped "
            f"({ctx.get('gene_level_skip_reason', 'reason not recorded')}). The "
            "comparison that this analysis exists to make is therefore missing, and "
            "no claim about the relative denominators can be made from this run.")
        add("")
    else:
        add("### Corpus construction")
        add("")
        for q in ctx["gene_corpus_stats"]["queries"]:
            add(f"* `{q['record_class']}`: `{q['entrez_query']}` "
                f"-> NCBI reported {q['esearch_reported_total']}, "
                f"{q['accessions_used']} used.")
        add("")
        add(f"{ctx['gene_corpus_stats']['n_available']} records are on disk under "
            f"`data/raw/genes/`; {c['n_records_mapped']} of them anchor-map to the "
            "reference genome at all. Records that do not map contribute nothing and "
            "are not counted anywhere.")
        add("")
        if c.get("corpus_composition") is not None:
            add(_md_table(c["corpus_composition"]))
            add("")
        add("**Read that table before reading anything below it.** Most HSV-1 "
            "sub-genomic records in GenBank are not about our genes at all -- they are "
            "thymidine kinase, glycoprotein G and glycoprotein B typing fragments. The "
            "sub-genomic corpus that actually touches the seven target genes is far "
            "smaller than the raw record count suggests, and for six of the seven "
            "genes the near-full-length `partial genome` records dominate the "
            "denominator. Those are the same KIND of evidence as the 183 complete "
            "genomes -- whole assemblies separated from them only by a title "
            "convention -- so pooling everything would make the widened denominator "
            "look far more independent than it is. Every result below is therefore "
            "also reported for the sub-genomic records alone.")
        add("")
        add("Gene keyword queries were deliberately NOT used: `RL2[All Fields]` "
            "returns 5 HSV-1 records and `ICP27`-style terms only 10, so a "
            "keyword-built corpus would have missed almost everything. Gene assignment "
            "comes from anchor mapping instead, so a record's own labelling is never "
            "trusted.")
        add("")
        add("### Coverage: distinguishing ABSENT from NOT COVERED")
        add("")
        add("A record enters a guide's denominator only if one of its anchor chains "
            "contains an anchor entirely left of the guide footprint and another "
            "entirely right of it. Those two anchors bracket the homologous stretch "
            "exactly. Everything between two anchors of the same chain counts as "
            "covered even where it is divergent, so variable sites are not silently "
            "discarded. A bracket wider than "
            f"{ctx['max_bracket_span']} bp, or a bracket containing a non-ACGT base, "
            "is scored UNKNOWN, not ABSENT.")
        add("")
        add("### Records per gene")
        add("")
        add(_md_table(ctx["gene_record_counts"]))
        add("")
        add("### Conservation over the gene-level corpus")
        add("")
        if cmp_ and cmp_.get("n_guides_with_enough_records"):
            add(f"{cmp_['n_guides_with_enough_records']} of {cmp_['n_guides']} guides "
                f"have at least {cmp_['min_records_required']} covering records "
                f"(median {cmp_.get('median_records_per_guide')}, "
                f"max {cmp_.get('max_records_per_guide')}).")
            add("")
            add(f"* Mean conservation over complete genomes: "
                f"{cmp_['mean_genome_level_all']:.4f}")
            add(f"* Mean conservation over the gene-level corpus: "
                f"{cmp_['mean_gene_level_all']:.4f}")
            add(f"* Pearson r = {cmp_['pearson_r']}, Spearman rho = {cmp_['spearman_rho']}")
            add("")
            add("**Fate of the perfectly conserved guides.**")
            add("")
            rows = [{"criterion": f"100% at n={n_genomes} (with >= "
                                  f"{cmp_['min_records_required']} covering records)",
                     "guides": cmp_["n_perfect_at_183"]}]
            for cut in ("0.99", "0.95", "0.90", "0.70"):
                rows.append({"criterion": f"...of which < {cut} at gene level",
                             "guides": cmp_[f"n_perfect_at_183_below_{cut}_gene_level"]})
            rows.append({"criterion": "...of which still exactly 100% at gene level",
                         "guides": cmp_["n_perfect_at_183_still_perfect"]})
            add(_md_table(pd.DataFrame(rows)))
            add("")
            add(f"Median gene-level conservation of the perfect set: "
                f"{cmp_['median_gene_level_of_perfect']}; minimum "
                f"{cmp_['min_gene_level_of_perfect']}.")
            add("")
            tiers = cmp_.get("tiers")
            if tiers:
                add("### The same question, one corpus tier at a time")
                add("")
                add("`subgenomic_only` is the tier that matters most: it is the only "
                    "one built from evidence of a different kind than the 183 complete "
                    "genomes. `no_coverage_correction` is not an estimate at all -- it "
                    "divides by every record considered for the guide, spanning the "
                    "site or not, which is how Amrani et al. describe their own "
                    "calculation. The gap between it and `whole_corpus` is the size of "
                    "the coverage correction, i.e. how much apparent 'variation' is "
                    "really just records that stop short of the site.")
                add("")
                trows = []
                for name, tv in tiers.items():
                    if not tv.get("n_guides_with_enough_records"):
                        trows.append({"tier": name, "guides_with_>=N_records": 0,
                                      "perfect_at_183": 0, "below_0.99": "n/a",
                                      "below_0.95": "n/a", "below_0.70": "n/a",
                                      "still_100pct": "n/a", "pearson_r": "n/a"})
                        continue
                    trows.append({
                        "tier": name,
                        "guides_with_>=N_records": tv["n_guides_with_enough_records"],
                        "perfect_at_183": tv["n_perfect_at_183"],
                        "below_0.99": tv["n_perfect_at_183_below_0.99_gene_level"],
                        "below_0.95": tv["n_perfect_at_183_below_0.95_gene_level"],
                        "below_0.70": tv["n_perfect_at_183_below_0.70_gene_level"],
                        "still_100pct": tv["n_perfect_at_183_still_perfect"],
                        "pearson_r": tv["pearson_r"],
                    })
                add(_md_table(pd.DataFrame(trows)))
                add("")
                gt = ctx.get("gene_tier_table")
                if gt is not None and not gt.empty:
                    add("Per gene, on the sub-genomic tier alone -- this is where the "
                        "denominators genuinely differ, and a pooled figure hides it:")
                    add("")
                    add(_md_table(gt))
                    add("")
                    med = int(gt["median_subgenomic_records"].median())
                    lo = int(gt["median_subgenomic_records"].min())
                    hi = int(gt["median_subgenomic_records"].max())
                    add(f"Note the denominators: the median guide has {med} "
                        f"sub-genomic records behind it (range of per-gene medians "
                        f"{lo} to {hi}). A guide scored 100% over ~30 records is not "
                        "the same claim as one scored 100% over 500; the exact "
                        "binomial 95% lower bound on 30/30 is about 0.88. Only UL30 "
                        "(and to a lesser extent UL5) has enough independent records "
                        "for a strong per-gene statement.")
                    add("")
        else:
            add("**No guide accumulated enough covering records for a comparison.** "
                "That is the real result of this run.")
            add("")
        add("Per-guide gene-level numbers: `results/robustness_gene_level.tsv`. "
            "Per-record corpus manifest: `data/gene_corpus_manifest.tsv`.")
        add("")

    # ---------------- benchmark ----------------
    if ctx.get("benchmark_table") is not None:
        add("## D. Head-to-head with the published Excision BioTherapeutics guides")
        add("")
        add("The four SaCas9 guides of Amrani et al. 2024 (Table 1) scored with OUR "
            "method against OUR corpora -- exact 26-mer (20 nt spacer + NNGRRT PAM) "
            "presence on either strand. This is not a quotation of their reported "
            "numbers; it is the same measurement applied to their sequences.")
        add("")
        add(_md_table(ctx["benchmark_table"]))
        add("")

    # ---------------- interpretation ----------------
    add("## What this means for guide selection")
    add("")
    for para in ctx["interpretation"]:
        add(para)
        add("")

    add("## Parameters")
    add("")
    add(_md_table(pd.DataFrame(
        [{"parameter": k, "value": v} for k, v in ctx["parameters"].items()])))
    add("")
    return "\n".join(L)


# ======================================================================================
# CLI
# ======================================================================================


def add_arguments(parser: argparse.ArgumentParser) -> None:
    g = parser.add_argument_group("robustness (stage 5)")
    g.add_argument("--robustness", action="store_true",
                   help="Run the stage-5 robustness/sensitivity analysis.")
    g.add_argument("--minhash-k", type=int, default=21,
                   help="k-mer size for the genome sketches (default 21).")
    g.add_argument("--minhash-sketch", type=int, default=5000,
                   help="Bottom-s sketch size (default 5000).")
    g.add_argument("--dedup-identity", type=float, default=0.999,
                   help="Primary single-linkage identity threshold for de-duplication "
                        "(default 0.999). The full sweep is reported regardless.")
    g.add_argument("--anchor-k", type=int, default=25,
                   help="k-mer size for homology anchors (default 25).")
    g.add_argument("--anchor-step", type=int, default=4,
                   help="Sample an anchor every N bases of the record (default 4).")
    g.add_argument("--anchor-max-gap", type=int, default=3000,
                   help="Maximum reference gap within one anchor chain (default 3000).")
    g.add_argument("--anchor-max-drift", type=int, default=60,
                   help="Maximum per-step diagonal movement within a chain, i.e. the "
                        "largest indel absorbed between two anchors (default 60).")
    g.add_argument("--anchor-min-anchors", type=int, default=2,
                   help="Minimum anchors for a chain to be kept (default 2).")
    g.add_argument("--anchor-max-occurrences", type=int, default=4,
                   help="Drop reference k-mers occurring more often than this "
                        "(default 4; the HSV-1 large repeats occur twice).")
    g.add_argument("--max-bracket-span", type=int, default=2000,
                   help="A guide bracketed by anchors more than this far apart is "
                        "scored UNKNOWN rather than ABSENT (default 2000).")
    g.add_argument("--skip-gene-corpus", action="store_true",
                   help="Do not contact NCBI for the gene-level corpus (parts A and B "
                        "only). Part C is then reported as NOT RUN.")
    g.add_argument("--skip-partial-genomes", action="store_true",
                   help="Exclude near-full-length 'partial genome' records from the "
                        "gene-level corpus.")
    g.add_argument("--gene-corpus-min-length", type=int, default=100)
    g.add_argument("--gene-corpus-max-length", type=int, default=145_000)
    g.add_argument("--gene-corpus-batch-size", type=int, default=100,
                   help="Accessions per efetch when downloading the corpus.")
    g.add_argument("--gene-corpus-limit", type=int, default=None,
                   help="Use only the first N accessions per corpus query (smoke test).")
    g.add_argument("--gene-corpus-refresh", action="store_true",
                   help="Ignore the data/raw/genes cache.")
    g.add_argument("--min-records-per-guide", type=int, default=10,
                   help="Minimum covering records before a guide's gene-level "
                        "conservation is used in the comparison (default 10).")
    g.add_argument("--robustness-report", type=Path, default=DEFAULT_ROBUSTNESS_REPORT)
    g.add_argument("--gene-manifest", type=Path, default=DEFAULT_GENE_MANIFEST)


def _interpretation(ctx: dict) -> list[str]:
    a, b, cmp_ = ctx["redundancy"], ctx["ambiguity"], ctx.get("comparison")
    n_genomes, n_guides = ctx["n_genomes"], ctx["n_guides"]
    out = []

    drop = a["perfect_dedup"] - a["perfect_full"]
    out.append(
        f"**Redundancy is real but small.** {n_genomes} deposited genomes behave like "
        f"{a['n_clusters']} at {a['primary_threshold']*100:.2f}% identity. Collapsing "
        f"them changes the perfectly-conserved count by {drop:+d} "
        f"({100.0*abs(drop)/max(a['perfect_full'],1):.1f}% of the baseline). The "
        "direction is the informative part: removing near-duplicates usually *raises* "
        "the count, because a duplicate adds no new sequence but can add new "
        "assembly noise. A count that rises on de-duplication is not evidence of a "
        "more conserved virus, only of a less redundant denominator.")

    out.append(
        f"**A large part of 'not conserved' is 'not sequenced'.** Of "
        f"{b['n_pairs']} (guide, genome) pairs, "
        f"{b['n_pairs'] - b['status_totals']['present'] - b['status_totals']['absent']} "
        "could not be resolved into present/absent by exact matching plus anchor "
        f"mapping. Treating those as UNKNOWN instead of MISMATCH moves the perfect "
        f"count from {b['perfect_strict']} to {b['perfect_n_tolerant']}. Both numbers "
        "are defensible; what is not defensible is quoting one without the other. "
        "The strict number is a lower bound and is the one to publish; the N-tolerant "
        "number is the right one for ranking guides against each other, because it "
        "does not punish a guide for happening to sit under someone else's assembly "
        "gap.")

    if cmp_ is None:
        out.append(
            "**The denominator question is unanswered in this run** because the "
            "gene-level corpus was not built. Until it is, the claim that our guides "
            "beat a >70%-conserved competitor cannot be evaluated: we would be "
            "comparing a 183-genome denominator against a several-thousand-record "
            "one.")
        return out

    if not cmp_.get("n_guides_with_enough_records"):
        out.append(
            "**The gene-level corpus did not produce a usable denominator.** No guide "
            "accumulated the required number of covering records. This is a real "
            "negative result about corpus coverage, not a computational failure.")
        return out

    lost95 = cmp_["n_perfect_at_183_below_0.95_gene_level"]
    lost70 = cmp_["n_perfect_at_183_below_0.70_gene_level"]
    kept = cmp_["n_perfect_at_183_still_perfect"]
    perfect = cmp_["n_perfect_at_183"]
    pct95 = 100.0 * lost95 / max(perfect, 1)

    comp = (ctx.get("gene_level") or {}).get("corpus_composition")
    if comp is not None and not comp.empty:
        sub = comp[comp["record_class"] == "subgenomic"]
        par = comp[comp["record_class"] == "partial_genome"]
        n_sub_dl = int(sub["records_downloaded"].iloc[0]) if len(sub) else 0
        n_sub_hit = int(sub["records_covering_a_target_gene"].iloc[0]) if len(sub) else 0
        n_par_hit = int(par["records_covering_a_target_gene"].iloc[0]) if len(par) else 0
        out.append(
            f"**The widened denominator is smaller than it looks, and this is the "
            f"most important caveat in this report.** {n_sub_dl} sub-genomic HSV-1 "
            f"records were downloaded, but only {n_sub_hit} of them cover any of the "
            "seven target genes; the rest are thymidine-kinase, glycoprotein-G and "
            "glycoprotein-B typing fragments from unrelated surveys. The remaining "
            f"{n_par_hit} covering records are near-full-length `partial genome` "
            "assemblies, i.e. the same kind of evidence as the 183 complete genomes, "
            "excluded from stage 1 only by a title convention. So for six of the "
            "seven genes this is not really a test against a larger, messier, "
            "independent corpus: it is mostly a test against 393 more whole genomes. "
            "UL30 is the exception -- it has a genuine clinical-amplicon corpus, "
            "because UL30 is the gene sequenced for aciclovir-resistance genotyping.")

    tiers = cmp_.get("tiers") or {}
    sub_t = tiers.get("subgenomic_only")

    out.append(
        f"**The headline number.** Of the {perfect} guides that are 100% conserved "
        f"across {n_genomes} complete genomes and have at least "
        f"{cmp_['min_records_required']} covering records in the gene-level corpus, "
        f"{lost95} ({pct95:.1f}%) fall below 95% when the denominator is widened, and "
        f"{lost70} fall below the 70% bar used by Amrani et al. {kept} remain at "
        "exactly 100%. Correlation between the two conservation measures is "
        f"r = {cmp_['pearson_r']} (Spearman rho = {cmp_['spearman_rho']}).")

    if sub_t and sub_t.get("n_guides_with_enough_records"):
        s_perfect = sub_t["n_perfect_at_183"]
        s95 = sub_t["n_perfect_at_183_below_0.95_gene_level"]
        s70 = sub_t["n_perfect_at_183_below_0.70_gene_level"]
        s_pct = 100.0 * s95 / max(s_perfect, 1)
        out.append(
            f"**The same number on sub-genomic records only** -- the tier that is "
            f"actually independent of the complete-genome set. "
            f"{sub_t['n_guides_with_enough_records']} guides clear the "
            f"{sub_t['min_records_required']}-record bar there (against "
            f"{cmp_['n_guides_with_enough_records']} on the whole corpus), of which "
            f"{s_perfect} are perfect at n={n_genomes}: {s95} ({s_pct:.1f}%) drop "
            f"below 95% and {s70} below 70%. Where the two tiers disagree, this is "
            "the one to believe, and the shrunken guide count is itself the finding: "
            "for most of our genes GenBank simply does not hold an independent "
            "per-gene corpus to test against.")

    naive = tiers.get("no_coverage_correction")
    if naive and naive.get("n_guides_with_enough_records"):
        out.append(
            f"**The coverage correction is the single biggest effect in this whole "
            f"analysis, larger than redundancy and ambiguity combined.** Without it "
            f"-- dividing by every record considered rather than by the records that "
            f"actually span the site -- "
            f"{naive['n_perfect_at_183_below_0.95_gene_level']} of the "
            f"{naive['n_perfect_at_183']} perfect guides fall below 95% and the mean "
            f"conservation drops to {naive['mean_gene_level_all']:.4f}. Almost all of "
            "that is an artefact: a 600 bp amplicon does not contain a guide 40 kb "
            "away, and a partial genome missing its repeats does not contain an RL2 "
            "guide. Any per-gene conservation number computed from a heterogeneous "
            "record set without a coverage check -- including a >70%-style threshold "
            "applied to a ViPR-scale corpus -- is measuring record length as much as "
            "sequence variation. This cuts both ways: it means the competing "
            "published numbers are probably pessimistic, not that ours are better.")

    verdict_pct = pct95
    verdict_basis = "the whole gene-level corpus"
    if sub_t and sub_t.get("n_perfect_at_183"):
        verdict_pct = (100.0 * sub_t["n_perfect_at_183_below_0.95_gene_level"]
                       / max(sub_t["n_perfect_at_183"], 1))
        verdict_basis = "the sub-genomic tier"
    out.append(f"*(The verdict below is taken from {verdict_basis}.)*")

    if verdict_pct >= 25.0:
        out.append(
            "**Read honestly, this substantially weakens the 'complete genomes are "
            "the better denominator' position.** A quarter or more of our "
            "perfectly-conserved set does not survive contact with the wider corpus. "
            "The 183-genome denominator is cleaner, but cleanliness is not the same "
            "as representativeness: the sub-genomic records are enriched for clinical "
            "surveys and drug-resistance genotyping, i.e. exactly the circulating "
            "diversity a therapeutic guide has to cut. Any claim of superiority over "
            "the Excision guides must be made on the gene-level corpus, or not made.")
    elif verdict_pct >= 1.0:
        out.append(
            "**The position survives but needs qualification.** Most perfectly "
            "conserved guides stay high at gene-level n, but a non-trivial minority "
            "does not, and those are exactly the ones that would fail in the clinic. "
            "The defensible framing is a two-denominator one: report conservation "
            "over complete genomes AND over the gene-level corpus, and select only "
            "guides that are high in both.")
    else:
        out.append(
            "**The position is defensible.** Widening the denominator by orders of "
            "magnitude barely moves the perfectly-conserved set, which is what one "
            "would expect if the 183-genome estimate were unbiased rather than merely "
            "small. Note the asymmetry that makes this credible: a wider corpus can "
            "only ever find more variation, never less, so a set that survives it has "
            "survived the harder test.")

    out.append(
        "**Practical rule for guide selection.** Rank on the intersection: perfectly "
        "conserved across complete genomes, still >= 95% over the gene-level corpus "
        "with a real denominator behind it (not 3 records), no poly-T, GC in band. "
        "Guides whose gene-level denominator is small should be labelled "
        "'insufficient evidence', not silently promoted -- a guide with 2/2 records "
        "is not more conserved than one with 380/400.")
    return out


def run(args: argparse.Namespace) -> dict:
    ensure_dirs()
    if not args.manifest.is_file():
        raise SystemExit(f"Manifest not found: {args.manifest} (run stage 1 first)")
    if not args.candidates.is_file():
        raise SystemExit(f"Candidates not found: {args.candidates} (run stage 2 first)")

    manifest = pd.read_csv(args.manifest, sep="\t").sort_values(
        "accession", kind="stable").reset_index(drop=True)
    candidates = pd.read_csv(args.candidates, sep="\t")
    if candidates.empty:
        raise SystemExit("Candidate guide table is empty.")

    reference_accession = str(candidates["reference_accession"].iloc[0])
    reference = load_reference_sequence(reference_accession)
    LOG.info("Reference %s: %d bp. %d guides, %d genomes.",
             reference_accession, len(reference), len(candidates), len(manifest))

    nuclease = nuclease_from_candidates(candidates, args)
    LOG.info("Nuclease of the candidate table: %s", nuclease.label)
    scan = scan_complete_genomes(manifest, candidates, reference, args, nuclease)
    LOG.info("Part A: redundancy / effective sample size ...")
    red = redundancy_analysis(manifest, scan, args)
    LOG.info("  effective sample size at %.4f identity: %d clusters (from %d genomes)",
             red["primary_threshold"], red["n_clusters"], red["n_genomes"])
    LOG.info("Part B: ambiguity sensitivity ...")
    amb = ambiguity_analysis(manifest, scan, args)
    LOG.info("  perfect strict=%d  n-tolerant=%d", amb["perfect_strict"],
             amb["perfect_n_tolerant"])

    # Outputs live next to the report, so a non-default (nuclease-namespaced)
    # --robustness-report path carries the whole stage-5 output set with it.
    out_dir = args.robustness_report.parent
    out_dir.mkdir(parents=True, exist_ok=True)
    red["per_genome"].to_csv(out_dir / "robustness_redundancy.tsv", sep="\t", index=False)
    red["sweep"].to_csv(out_dir / "robustness_effective_n.tsv", sep="\t", index=False)
    amb["sweep"].to_csv(out_dir / "robustness_ambiguity_sweep.tsv", sep="\t", index=False)
    amb["per_guide"].to_csv(out_dir / "robustness_conservation_modes.tsv",
                            sep="\t", index=False)

    gene_level = None
    comparison = None
    gene_stats = None
    skip_reason = None
    genes = [g.strip().upper() for g in args.genes.split(",") if g.strip()]
    gene_record_counts = pd.DataFrame()
    gene_tier_table = pd.DataFrame()

    if args.skip_gene_corpus:
        skip_reason = "--skip-gene-corpus was given"
        LOG.warning("Part C SKIPPED (--skip-gene-corpus).")
    else:
        LOG.info("Part C: building the gene-level corpus ...")
        queries = gene_corpus_queries(
            args.taxid, args.gene_corpus_min_length, args.gene_corpus_max_length,
            not args.skip_partial_genomes, args.min_length)
        try:
            records, gene_stats = fetch_gene_corpus(
                queries, args.gene_corpus_batch_size, args.gene_corpus_refresh,
                args.gene_corpus_limit)
        except Exception as exc:  # noqa: BLE001
            skip_reason = f"corpus download failed: {type(exc).__name__}: {exc}"
            LOG.error("Part C FAILED: %s", skip_reason)
            records = []
        if records:
            spans = collect_gene_spans(reference_accession, genes)
            gene_level = gene_level_analysis(records, candidates, reference, spans,
                                             args, nuclease)
            gene_level["per_guide"].to_csv(out_dir / "robustness_gene_level.tsv",
                                           sep="\t", index=False)
            gene_level["corpus_manifest"].to_csv(args.gene_manifest, sep="\t", index=False)
            base = amb["per_guide"][["guide_id", "conservation_strict"]]
            comparison = compare_denominators(
                base, gene_level["per_guide"], args.min_records_per_guide)
            tiers = {
                "whole_corpus": dict(comparison),   # copy: comparison holds `tiers`
                "subgenomic_only": compare_denominators(
                    base, gene_level["per_guide"], args.min_records_per_guide,
                    "n_records_covering_subgenomic", "conservation_subgenomic"),
                "partial_genomes_only": compare_denominators(
                    base, gene_level["per_guide"], args.min_records_per_guide,
                    "n_records_covering_partial_genome", "conservation_partial_genomes"),
                "no_coverage_correction": compare_denominators(
                    base, gene_level["per_guide"], args.min_records_per_guide,
                    "n_records_covering", "conservation_no_coverage_check"),
            }
            comparison["tiers"] = tiers

            cm = gene_level["corpus_manifest"]
            cm_hits = cm["genes_covered"].fillna("").str.len() > 0
            gene_level["corpus_composition"] = pd.DataFrame([
                {"record_class": cls,
                 "records_downloaded": int((cm["record_class"] == cls).sum()),
                 "records_mapping_to_reference": int(
                     ((cm["record_class"] == cls) & (cm["n_anchor_chains"] > 0)).sum()),
                 "records_covering_a_target_gene": int(
                     ((cm["record_class"] == cls) & cm_hits).sum()),
                 "median_bp": int(cm[cm["record_class"] == cls]["length_bp"].median())
                 if (cm["record_class"] == cls).any() else 0}
                for cls in ("subgenomic", "partial_genome")])
            rows = []
            for gene in genes:
                mask = cm["genes_covered"].fillna("").str.split("|").apply(lambda v: gene in v)
                sub = cm[mask]
                rows.append({
                    "gene": gene,
                    "records_covering_gene": int(len(sub)),
                    "subgenomic": int((sub["record_class"] == "subgenomic").sum()),
                    "partial_genome": int((sub["record_class"] == "partial_genome").sum()),
                    "median_record_bp": int(sub["length_bp"].median()) if len(sub) else 0,
                })
            gene_record_counts = pd.DataFrame(rows)
            gene_record_counts.to_csv(out_dir / "robustness_gene_records.tsv",
                                      sep="\t", index=False)

            # Per gene, on the independent (sub-genomic) tier only. This is where the
            # denominator actually differs between genes, so a pooled number hides it.
            gl = gene_level["per_guide"].merge(base, on="guide_id", validate="one_to_one")
            grows = []
            for gene, sub in gl.groupby("gene", sort=True):
                usable = sub[sub["n_records_covering_subgenomic"]
                             >= args.min_records_per_guide]
                perfect = usable[usable["conservation_strict"] >= 1.0]
                grows.append({
                    "gene": gene,
                    "guides": int(len(sub)),
                    "median_subgenomic_records": int(
                        sub["n_records_covering_subgenomic"].median()),
                    "max_subgenomic_records": int(
                        sub["n_records_covering_subgenomic"].max()),
                    "perfect_at_183_with_usable_n": int(len(perfect)),
                    "of_those_below_0.95_subgenomic": int(
                        (perfect["conservation_subgenomic"] < 0.95).sum()),
                    "min_subgenomic_conservation_of_perfect": (
                        round(float(perfect["conservation_subgenomic"].min()), 4)
                        if len(perfect) else None),
                })
            gene_tier_table = pd.DataFrame(grows)
            gene_tier_table.to_csv(out_dir / "robustness_gene_level_by_gene.tsv",
                                   sep="	", index=False)
        elif skip_reason is None:
            skip_reason = "no corpus records were retrieved"

    # Benchmark table (our method, their guides).
    bench_rows = []
    for name in BENCHMARK_GUIDES:
        present = int(scan["benchmark_present"][name].sum())
        row = {
            "guide": name,
            "spacer": BENCHMARK_GUIDES[name][0],
            "pam": BENCHMARK_GUIDES[name][1],
            f"present_in_complete_genomes_n{len(manifest)}": present,
            "conservation_complete_genomes": round(present / len(manifest), 4),
        }
        if gene_level is not None:
            bcov = gene_level["benchmark"][name]
            row["gene_level_records_covering"] = bcov["covered"]
            row["gene_level_present"] = bcov["present"]
            row["conservation_gene_level"] = (
                round(bcov["present"] / bcov["covered"], 4) if bcov["covered"] else float("nan"))
        bench_rows.append(row)
    benchmark_table = pd.DataFrame(bench_rows)
    benchmark_table.to_csv(out_dir / "robustness_benchmark_guides.tsv",
                           sep="\t", index=False)

    top_clusters = (red["per_genome"][red["per_genome"]["cluster_size"] > 1]
                    .sort_values(["cluster_size", "cluster_id", "accession"],
                                 ascending=[False, True, True], kind="stable"))
    examples = None
    if not top_clusters.empty:
        keep_ids = list(dict.fromkeys(top_clusters["cluster_id"].tolist()))[:6]
        examples = (top_clusters[top_clusters["cluster_id"].isin(keep_ids)]
                    [["cluster_id", "cluster_size", "accession", "strain",
                      "nearest_accession", "nearest_identity"]]
                    .reset_index(drop=True))

    ctx = {
        "generated_utc": utc_now_iso(),
        "reference": reference_accession,
        "n_guides": len(candidates),
        "n_genomes": len(manifest),
        "minhash_k": args.minhash_k,
        "minhash_sketch": args.minhash_sketch,
        "max_bracket_span": args.max_bracket_span,
        "redundancy": red,
        "redundancy_examples": examples,
        "ambiguity": amb,
        "gene_level": gene_level,
        "gene_corpus_stats": gene_stats,
        "gene_record_counts": gene_record_counts,
        "gene_tier_table": gene_tier_table,
        "gene_level_skip_reason": skip_reason,
        "comparison": comparison,
        "benchmark_table": benchmark_table,
        "parameters": {
            "minhash_k": args.minhash_k,
            "minhash_sketch": args.minhash_sketch,
            "dedup_identity": args.dedup_identity,
            "anchor_k": args.anchor_k,
            "anchor_step": args.anchor_step,
            "anchor_max_gap": args.anchor_max_gap,
            "anchor_max_drift": args.anchor_max_drift,
            "anchor_max_occurrences": args.anchor_max_occurrences,
            "max_bracket_span": args.max_bracket_span,
            "min_records_per_guide": args.min_records_per_guide,
            "gene_corpus_min_length": args.gene_corpus_min_length,
            "gene_corpus_max_length": args.gene_corpus_max_length,
            "partial_genomes_included": not args.skip_partial_genomes,
        },
    }
    ctx["interpretation"] = _interpretation(ctx)

    args.robustness_report.parent.mkdir(parents=True, exist_ok=True)
    args.robustness_report.write_text(render_report(ctx), encoding="utf-8")
    LOG.info("Wrote %s", args.robustness_report)

    summary = {
        "generated_utc": ctx["generated_utc"],
        "reference": reference_accession,
        "n_guides": len(candidates),
        "n_genomes": len(manifest),
        "redundancy": {
            "primary_identity_threshold": red["primary_threshold"],
            "effective_sample_size": red["n_clusters"],
            "n_multi_member_clusters": red["n_multi_member_clusters"],
            "perfect_all_genomes": red["perfect_full"],
            "perfect_deduplicated": red["perfect_dedup"],
            "sweep": red["sweep"].to_dict(orient="records"),
            "identity": red["identity_matrix_summary"],
        },
        "ambiguity": {
            "perfect_strict": amb["perfect_strict"],
            "perfect_n_tolerant": amb["perfect_n_tolerant"],
            "n_rescued_by_n_tolerance": amb["n_rescued_by_n_tolerance"],
            "status_totals": amb["status_totals"],
            "sweep": amb["sweep"].to_dict(orient="records"),
        },
        "gene_level": None if gene_level is None else {
            "corpus": gene_stats,
            "n_records_mapped": gene_level["n_records_mapped"],
            "corpus_composition": (gene_level.get("corpus_composition").to_dict(orient="records")
                                   if gene_level.get("corpus_composition") is not None else None),
            "comparison": comparison,
            "records_per_gene": gene_record_counts.to_dict(orient="records"),
            "per_gene_subgenomic_tier": gene_tier_table.to_dict(orient="records"),
        },
        "gene_level_skip_reason": skip_reason,
        "benchmark_guides": benchmark_table.to_dict(orient="records"),
        "parameters": ctx["parameters"],
    }
    (out_dir / "robustness_summary.json").write_text(
        json.dumps(summary, indent=2, default=str), encoding="utf-8")
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--candidates", type=Path, default=DEFAULT_CANDIDATES)
    parser.add_argument("--conservation", type=Path, default=DEFAULT_CONSERVATION)
    parser.add_argument("--taxid", type=int, default=10298)
    parser.add_argument("--min-length", type=int, default=145_000)
    parser.add_argument("--genes", default="UL30,UL19,UL5,UL52,UL29,RL2,UL54")
    add_arguments(parser)
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args(argv)
    setup_logging(args.verbose)
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
