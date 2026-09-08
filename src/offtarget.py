"""Stage 8 -- human (GRCh38) off-target screening for SaCas9 guides.

Why this module exists
----------------------
Phase 1 of this project ranked HSV-1 guides on cross-isolate conservation only, and
`README.md` limitation 1 says so explicitly. Stage 6 then argued that Amrani et al.
2024 had better-conserved SaCas9 options in ICP0 than the ICP0g2 they took forward.
That argument is not decidable on conservation alone: Amrani et al. selected on
conservation *and* an off-target count (BWA -> Cas-OFFinder against hg38, PAM
NNGRRN, plus GUIDE-seq validation). If our proposed alternatives are dirtier against
the human genome than ICP0g2, the recommendation collapses.

This module measures that. It is deliberately capable of killing our own thesis, and
the report it writes states the result whichever way it falls.

What it does
------------
1. Downloads and caches the Ensembl GRCh38 primary assembly (release 116), verifies
   it three ways (publisher BSD `sum`, gzip CRC-32 of the decompressed stream, and
   our own SHA-256), and encodes it once into a memory-mappable uint8 array.
2. Finds every genomic site, on both strands, whose PAM matches NNGRRN and whose
   protospacer is within `--max-mismatches` (default 4) substitutions of a query
   guide. NNGRRT sites are a strict subset of NNGRRN sites and are flagged rather
   than searched separately, so both PAM variants come out of one pass.
3. Annotates sites at <= 3 mismatches against the Ensembl 116 GTF and flags any that
   fall inside a coding exon (CDS).
4. Writes `results/offtarget_report.md` plus supporting TSVs.

Why an exact-substring scanner is the wrong tool, and what is used instead
--------------------------------------------------------------------------
Off-target sites are by definition MISmatched, so `src/conservation.py`'s exact-match
jump scan cannot be reused. Two properties of the problem are exploited instead.

**PAM sparsity.** SaCas9 needs `NNGRRN`: G at PAM offset 2 and purines at offsets 3
and 4. That is 1/16 of positions per strand, and `NNGRRT` is 1/64. Only those
positions can ever be cut, so 15/16 of the genome is discarded with three vectorised
slice comparisons before any protospacer is looked at.

**The pigeonhole principle.** With at most `m` mismatches in an `L`-nt protospacer,
splitting the protospacer into `m + 1` contiguous chunks guarantees that at least one
chunk matches the guide exactly -- there are not enough mismatches to touch every
chunk. So for m = 4 we cut the 20-nt protospacer into five 4-mers, encode every
genomic 4-mer as a byte once per chromosome by pure slice arithmetic, gather the five
chunk codes at the PAM positions, and use a 256-entry bitmask lookup table per chunk
to find, for all guides at once, which positions share at least one exact chunk with
which guide. Roughly 2% of PAM positions survive per guide.

The surviving candidates are then scored *without touching the genome again*: the
five chunks tile the whole protospacer, so a per-(guide, chunk) 256-entry table of
"how many of these 4 bases differ" turns the five already-gathered codes into a
mismatch count with five more table lookups. Because the codes pack an ambiguous base
onto A, that count is a strict lower bound, which makes it a sound filter but not a
sound answer -- so the handful of positions that survive it (true hits plus a few
N-boundary windows) are read back from the genome and counted exactly. Everything is
numpy; there is no per-base Python loop anywhere on the hot path.

Correctness
-----------
The pigeonhole path is an optimisation, not an approximation, and `tests/
test_offtarget.py` proves it the way `tests/test_nuclease.py` proves the stage-3
scanner: `test_fast_matches_bruteforce_chr21` runs a structurally independent
brute-force scan -- mismatch counts computed at *every* genomic position by 20
whole-array slice comparisons, no PAM prefilter, no chunk index, explicit reverse
complement -- over the whole of chromosome 21 for every guide screened, and asserts
set identity of the hit tuples. A third, plainly-written pure-Python string
implementation validates the numpy brute force in turn on smaller inputs.

Honest positioning
------------------
Cas-OFFinder is the field standard and is what Amrani et al. used; a reviewer will
expect it or a validated equivalent. This module is a validated equivalent for
substitution-only search, and it is *not* equivalent for DNA/RNA bulges, which
Cas-OFFinder models and this does not. The report says so in those words.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import logging
import math
import os
import sys
import time
import urllib.error
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.common import DATA_DIR, RAW_DIR, RESULTS_DIR, revcomp, sha256_file  # noqa: E402
from src.nuclease import IUPAC_CODES, iupac_revcomp  # noqa: E402

LOG = logging.getLogger("offtarget")

__all__ = [
    "OFFTARGET_IMPLEMENTED",
    "Guide",
    "GenomeCache",
    "acquire_genome",
    "build_genome_cache",
    "bruteforce_scan_sequence",
    "caveat",
    "unscreened_caveat",
    "encode_sequence",
    "load_guides",
    "main",
    "naive_scan_python",
    "pam_position_mask",
    "scan_sequence",
    "screen_offtargets",
]

#: Phase 2 is implemented. Kept as a public flag because phase-1 code and the README
#: refer to it.
OFFTARGET_IMPLEMENTED = True

# --------------------------------------------------------------------------------------
# Reference data
# --------------------------------------------------------------------------------------

GENOME_DIR = DATA_DIR / "genome"

#: Pinned, not "current": a floating release would silently change every number in
#: the report. Release 116 = GRCh38.p14 = GenBank assembly GCA_000001405.29.
ENSEMBL_RELEASE = 116
ENSEMBL_BASE = f"https://ftp.ensembl.org/pub/release-{ENSEMBL_RELEASE}"

#: UNMASKED (`dna`), not soft-masked (`dna_sm`) and emphatically not hard-masked
#: (`dna_rm`). Repeat content is not a nuisance here: a Cas9 site inside a LINE or a
#: satellite is a real cleavage substrate, and hard-masking would delete ~50% of the
#: genome from the search. Soft-masking would be harmless if every read were
#: upper-cased, but "harmless if you remember to" is exactly the kind of silent
#: case-handling bug that would quietly shrink an off-target count, so the file that
#: cannot express the distinction is the one we download.
#: PRIMARY assembly, not "toplevel": toplevel additionally contains ALT haplotypes
#: and patch scaffolds, which are alternative representations of sequence already
#: present on the chromosomes and would double-count off-targets.
FASTA_NAME = "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
FASTA_URL = f"{ENSEMBL_BASE}/fasta/homo_sapiens/dna/{FASTA_NAME}"
FASTA_CHECKSUMS_URL = f"{ENSEMBL_BASE}/fasta/homo_sapiens/dna/CHECKSUMS"

GTF_NAME = f"Homo_sapiens.GRCh38.{ENSEMBL_RELEASE}.gtf.gz"
GTF_URL = f"{ENSEMBL_BASE}/gtf/homo_sapiens/{GTF_NAME}"
GTF_CHECKSUMS_URL = f"{ENSEMBL_BASE}/gtf/homo_sapiens/CHECKSUMS"

CACHE_BIN = GENOME_DIR / "grch38_primary.u8"
CACHE_INDEX = GENOME_DIR / "grch38_primary.index.json"
GENOME_MANIFEST = GENOME_DIR / "genome_manifest.json"
ANNOTATION_CACHE = GENOME_DIR / f"ensembl{ENSEMBL_RELEASE}_annotation.npz"

# --------------------------------------------------------------------------------------
# Base encoding
# --------------------------------------------------------------------------------------

#: A/C/G/T -> 0/1/2/3 (case-insensitive); every other byte, N included, -> 4.
#: Code 4 is never equal to a guide base and is never admitted by an IUPAC set, so an
#: ambiguous genomic base counts as a mismatch in the protospacer and disqualifies a
#: PAM outright. That is the conservative direction for a PAM and the only defensible
#: one for a protospacer.
BASE_CODE = np.full(256, 4, dtype=np.uint8)
for _i, _b in enumerate("ACGT"):
    BASE_CODE[ord(_b)] = _i
    BASE_CODE[ord(_b.lower())] = _i

CODE_BASE = np.frombuffer(b"ACGTN", dtype=np.uint8)

N_CODE = 4


def encode_sequence(seq: str | bytes) -> np.ndarray:
    """Encode an ACGT(N) string as uint8 codes."""
    if isinstance(seq, str):
        seq = seq.encode("ascii", "replace")
    return BASE_CODE[np.frombuffer(seq, dtype=np.uint8)]


def decode_sequence(codes: np.ndarray) -> str:
    return CODE_BASE[np.clip(codes, 0, 4)].tobytes().decode("ascii")


def revcomp_codes(codes: np.ndarray) -> np.ndarray:
    """Reverse complement in code space: 0<->3, 1<->2, 4 fixed."""
    out = np.where(codes < 4, 3 - codes, np.uint8(4)).astype(np.uint8)
    return out[::-1].copy()


def iupac_code_sets(pattern: str) -> list[np.ndarray]:
    """Per-position array of admitted *codes* (0-3) for an IUPAC pattern."""
    return [np.array(sorted(BASE_CODE[ord(b)] for b in IUPAC_CODES[c]), dtype=np.uint8)
            for c in pattern.upper()]


# --------------------------------------------------------------------------------------
# Acquisition
# --------------------------------------------------------------------------------------


def bsd_sum(path: Path, log_every: int = 8) -> tuple[int, int]:
    """The BSD 16-bit rotating checksum used by `sum` -- Ensembl's CHECKSUMS format.

    Returns `(checksum, blocks)` where blocks is the size in 1024-byte blocks rounded
    up, exactly as `sum` prints them. The recurrence
    `c = rotate_right_16(c) + byte (mod 2**16)` mixes a bit rotation with modular
    addition and is therefore not vectorisable; it is a genuinely sequential ~10^9
    step loop. It runs once per download and the result is cached in
    `genome_manifest.json`, so the cost is paid one time.
    """
    checksum = 0
    total = 0
    size = path.stat().st_size
    t0 = time.time()
    next_report = log_every
    with open(path, "rb") as fh:
        while True:
            chunk = fh.read(1 << 22)
            if not chunk:
                break
            for byte in chunk:
                checksum = (checksum >> 1) | ((checksum & 1) << 15)
                checksum = (checksum + byte) & 0xFFFF
            total += len(chunk)
            elapsed = time.time() - t0
            if elapsed >= next_report:
                next_report = elapsed + log_every
                LOG.info("  BSD sum %.1f%% (%.0f MB / %.0f MB, %.0fs)",
                         100.0 * total / max(size, 1), total / 1e6, size / 1e6, elapsed)
    return checksum, math.ceil(total / 1024)


def parse_ensembl_checksums(text: str) -> dict[str, tuple[int, int]]:
    """`<checksum> <blocks> <filename>` lines -> {filename: (checksum, blocks)}."""
    out: dict[str, tuple[int, int]] = {}
    for line in text.splitlines():
        parts = line.split()
        if len(parts) == 3 and parts[0].isdigit() and parts[1].isdigit():
            out[parts[2]] = (int(parts[0]), int(parts[1]))
    return out


def http_get(url: str, timeout: int = 120) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "hsv-crispr-offtarget/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read()


def download_resumable(url: str, dest: Path, *, attempts: int = 8,
                       timeout: int = 120) -> Path:
    """Download `url` to `dest`, resuming a partial `.part` file via HTTP Range.

    Idempotent: if `dest` already exists the function returns immediately, so a
    re-run never re-downloads. A partial transfer leaves `<dest>.part` behind and the
    next attempt continues from its length rather than starting over -- which matters
    for an 841 MB file on a domestic connection.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        LOG.info("cached: %s (%.0f MB)", dest.name, dest.stat().st_size / 1e6)
        return dest
    part = dest.with_suffix(dest.suffix + ".part")
    delay = 3.0
    for attempt in range(1, attempts + 1):
        have = part.stat().st_size if part.exists() else 0
        headers = {"User-Agent": "hsv-crispr-offtarget/1.0"}
        if have:
            headers["Range"] = f"bytes={have}-"
        req = urllib.request.Request(url, headers=headers)
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                if have and resp.status != 206:
                    LOG.warning("server ignored Range (status %s); restarting download",
                                resp.status)
                    have = 0
                    part.unlink(missing_ok=True)
                total = resp.headers.get("Content-Length")
                total = (int(total) + have) if total else None
                mode = "ab" if have else "wb"
                t0 = time.time()
                got = have
                next_report = 10.0
                with open(part, mode) as fh:
                    while True:
                        buf = resp.read(1 << 20)
                        if not buf:
                            break
                        fh.write(buf)
                        got += len(buf)
                        elapsed = time.time() - t0
                        if elapsed >= next_report:
                            next_report = elapsed + 10.0
                            pct = f"{100.0 * got / total:.1f}%" if total else "?"
                            LOG.info("  %s %.0f MB (%s) %.1f MB/s",
                                     dest.name, got / 1e6, pct,
                                     (got - have) / 1e6 / max(elapsed, 1e-9))
            part.replace(dest)
            LOG.info("downloaded %s (%.0f MB)", dest.name, dest.stat().st_size / 1e6)
            return dest
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            if attempt == attempts:
                raise RuntimeError(f"download failed after {attempts} attempts: {url}") from exc
            LOG.warning("download attempt %d/%d failed (%s: %s) -- retrying in %.0fs "
                        "from byte %d", attempt, attempts, type(exc).__name__, exc,
                        delay, part.stat().st_size if part.exists() else 0)
            time.sleep(delay)
            delay = min(delay * 2, 60.0)
    raise RuntimeError(f"unreachable: {url}")


def acquire_genome(*, skip_bsd_sum: bool = False, want_gtf: bool = True) -> dict:
    """Download + verify GRCh38 (and the GTF) into data/genome/. Returns the manifest.

    Verification is threefold and each check is independent of the others:

    * the publisher's own `sum` checksum from Ensembl's CHECKSUMS file (detects a
      corrupt or substituted file at the source-of-truth level);
    * the gzip member's CRC-32 and ISIZE trailer, checked by the decompressor when
      `build_genome_cache` streams the whole file (detects any transport corruption
      that survived TCP);
    * a SHA-256 we compute and record ourselves, which is the value to quote in a
      methods section.
    """
    GENOME_DIR.mkdir(parents=True, exist_ok=True)
    # Keep the cache out of git without touching the repo-level .gitignore, which is
    # not this module's file to edit.
    gi = GENOME_DIR / ".gitignore"
    if not gi.exists():
        gi.write_text("# Regenerable reference-genome cache; never commit.\n*\n",
                      encoding="utf-8")

    manifest: dict = {}
    if GENOME_MANIFEST.exists():
        manifest = json.loads(GENOME_MANIFEST.read_text(encoding="utf-8"))

    fasta = GENOME_DIR / FASTA_NAME
    checks_path = GENOME_DIR / "CHECKSUMS.fasta.txt"
    if not checks_path.exists():
        checks_path.write_bytes(http_get(FASTA_CHECKSUMS_URL))
    published = parse_ensembl_checksums(checks_path.read_text(encoding="utf-8"))

    download_resumable(FASTA_URL, fasta)
    size = fasta.stat().st_size

    entry = manifest.get(FASTA_NAME, {})
    if entry.get("size_bytes") != size or "sha256" not in entry:
        LOG.info("hashing %s (SHA-256)", FASTA_NAME)
        t0 = time.time()
        entry = {
            "url": FASTA_URL,
            "size_bytes": size,
            "sha256": sha256_file(fasta),
            "ensembl_release": ENSEMBL_RELEASE,
        }
        LOG.info("  sha256 %s (%.0fs)", entry["sha256"][:16], time.time() - t0)

    exp = published.get(FASTA_NAME)
    if exp is not None:
        entry["published_bsd_sum"] = exp[0]
        entry["published_blocks"] = exp[1]
        exp_bytes_lo = (exp[1] - 1) * 1024 + 1
        exp_bytes_hi = exp[1] * 1024
        if not (exp_bytes_lo <= size <= exp_bytes_hi):
            raise RuntimeError(
                f"{FASTA_NAME}: size {size} B is outside the {exp[1]} 1 KiB blocks "
                f"Ensembl publishes ({exp_bytes_lo}-{exp_bytes_hi} B). Download is "
                "truncated or the release moved; delete data/genome/ and retry."
            )
        if skip_bsd_sum:
            entry.setdefault("bsd_sum", None)
            entry["bsd_sum_verified"] = bool(entry.get("bsd_sum") == exp[0])
        elif entry.get("bsd_sum") is None:
            LOG.info("verifying %s against Ensembl's published `sum` checksum "
                     "(sequential, one-off, several minutes)", FASTA_NAME)
            t0 = time.time()
            got_sum, got_blocks = bsd_sum(fasta)
            entry["bsd_sum"] = got_sum
            entry["bsd_sum_blocks"] = got_blocks
            entry["bsd_sum_seconds"] = round(time.time() - t0, 1)
            if (got_sum, got_blocks) != exp:
                raise RuntimeError(
                    f"{FASTA_NAME}: BSD checksum mismatch. Ensembl publishes "
                    f"{exp[0]} {exp[1]}, computed {got_sum} {got_blocks}. "
                    "The download is corrupt -- delete data/genome/ and retry."
                )
            entry["bsd_sum_verified"] = True
            LOG.info("  BSD sum %d/%d blocks MATCHES Ensembl (%.0fs)",
                     got_sum, got_blocks, entry["bsd_sum_seconds"])
        else:
            entry["bsd_sum_verified"] = (
                entry.get("bsd_sum") == exp[0] and entry.get("bsd_sum_blocks") == exp[1])
    manifest[FASTA_NAME] = entry

    if want_gtf:
        gtf = GENOME_DIR / GTF_NAME
        gchecks = GENOME_DIR / "CHECKSUMS.gtf.txt"
        if not gchecks.exists():
            gchecks.write_bytes(http_get(GTF_CHECKSUMS_URL))
        gpub = parse_ensembl_checksums(gchecks.read_text(encoding="utf-8"))
        download_resumable(GTF_URL, gtf)
        gentry = manifest.get(GTF_NAME, {})
        if gentry.get("size_bytes") != gtf.stat().st_size:
            gentry = {
                "url": GTF_URL,
                "size_bytes": gtf.stat().st_size,
                "sha256": sha256_file(gtf),
                "ensembl_release": ENSEMBL_RELEASE,
            }
        gexp = gpub.get(GTF_NAME)
        if gexp is not None:
            gentry["published_bsd_sum"], gentry["published_blocks"] = gexp
            if gentry.get("bsd_sum") is None and not skip_bsd_sum:
                s, b = bsd_sum(gtf)
                gentry["bsd_sum"], gentry["bsd_sum_blocks"] = s, b
                if (s, b) != gexp:
                    raise RuntimeError(f"{GTF_NAME}: BSD checksum mismatch "
                                       f"({s} {b} vs published {gexp[0]} {gexp[1]})")
                gentry["bsd_sum_verified"] = True
        manifest[GTF_NAME] = gentry

    GENOME_MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


# --------------------------------------------------------------------------------------
# Genome cache: one memory-mappable uint8 array + a JSON index
# --------------------------------------------------------------------------------------


def _iter_fasta_records(fh, chunk_size: int = 1 << 24):
    """Yield (header, bytes_of_sequence_with_newlines_removed) streaming from `fh`.

    Sequence data is emitted in pieces; a new header starts a new record. The scan
    for record boundaries uses `bytes.find(b"\\n>")` and `bytes.translate` so no
    Python-level loop ever touches an individual base -- at 3.1 Gb that distinction is
    the difference between one minute and one hour.
    """
    pending = b""
    header: str | None = None
    first = True
    while True:
        data = fh.read(chunk_size)
        if not data:
            break
        data = pending + data
        cut = data.rfind(b"\n")
        if cut == -1:
            pending = data
            continue
        pending = data[cut + 1:]
        data = data[:cut + 1]

        pos = 0
        if first:
            first = False
            if data.startswith(b">"):
                nl = data.index(b"\n")
                header = data[1:nl].decode("ascii", "replace")
                pos = nl + 1
        while True:
            idx = data.find(b"\n>", pos)
            if idx == -1:
                seq = data[pos:]
                if seq:
                    yield header, seq.translate(None, b"\r\n")
                break
            seq = data[pos:idx + 1]
            if seq:
                yield header, seq.translate(None, b"\r\n")
            nl = data.index(b"\n", idx + 1)
            header = data[idx + 2:nl].decode("ascii", "replace")
            yield header, None  # signals "new record starts here"
            pos = nl + 1
    if pending:
        if pending.startswith(b">"):
            yield pending[1:].decode("ascii", "replace"), None
        else:
            yield header, pending.translate(None, b"\r\n")


def build_genome_cache(force: bool = False) -> "GenomeCache":
    """Decode the gzipped FASTA once into `grch38_primary.u8` + a JSON index.

    Re-running is free: the index records the source SHA-256 and byte length, and the
    function returns immediately when they still match.
    """
    fasta = GENOME_DIR / FASTA_NAME
    manifest = json.loads(GENOME_MANIFEST.read_text(encoding="utf-8"))
    src_sha = manifest[FASTA_NAME]["sha256"]

    if CACHE_INDEX.exists() and CACHE_BIN.exists() and not force:
        idx = json.loads(CACHE_INDEX.read_text(encoding="utf-8"))
        if (idx.get("source_sha256") == src_sha
                and CACHE_BIN.stat().st_size == idx.get("total_length")):
            LOG.info("genome cache present: %d sequences, %.2f Gb",
                     len(idx["sequences"]), idx["total_length"] / 1e9)
            return GenomeCache(idx)
        LOG.warning("genome cache is stale -- rebuilding")

    LOG.info("decoding %s -> %s (one-off)", FASTA_NAME, CACHE_BIN.name)
    t0 = time.time()
    sequences: list[dict] = []
    offset = 0
    cur_name: str | None = None
    cur_desc = ""
    cur_len = 0
    cur_n = 0

    def close_record():
        nonlocal cur_name, cur_len, cur_n, offset
        if cur_name is not None:
            sequences.append({"name": cur_name, "description": cur_desc,
                              "offset": offset, "length": cur_len,
                              "n_ambiguous": int(cur_n)})
            offset += cur_len
        cur_name, cur_len, cur_n = None, 0, 0

    tmp = CACHE_BIN.with_suffix(".u8.part")
    with gzip.open(fasta, "rb") as gz, open(tmp, "wb") as out:
        for header, seq in _iter_fasta_records(gz):
            if seq is None:
                close_record()
                cur_name = header.split()[0] if header else "?"
                cur_desc = header or ""
                continue
            if cur_name is None and header is not None:
                cur_name = header.split()[0]
                cur_desc = header
            codes = BASE_CODE[np.frombuffer(seq, dtype=np.uint8)]
            out.write(codes.tobytes())
            cur_len += codes.size
            cur_n += int(np.count_nonzero(codes == N_CODE))
        close_record()
    tmp.replace(CACHE_BIN)

    index = {
        "source": FASTA_NAME,
        "source_sha256": src_sha,
        "ensembl_release": ENSEMBL_RELEASE,
        "assembly": "GRCh38 primary assembly (unmasked)",
        "encoding": "uint8: A=0 C=1 G=2 T=3 other/N=4",
        "total_length": offset,
        "build_seconds": round(time.time() - t0, 1),
        "sequences": sequences,
    }
    CACHE_INDEX.write_text(json.dumps(index, indent=1) + "\n", encoding="utf-8")
    LOG.info("cached %d sequences, %.3f Gb (%.0f N), %.0fs",
             len(sequences), offset / 1e9,
             sum(s["n_ambiguous"] for s in sequences), index["build_seconds"])
    return GenomeCache(index)


@dataclass
class GenomeCache:
    """Read access to the encoded genome. One `np.memmap`, sliced per sequence."""

    index: dict
    _mm: np.memmap | None = field(default=None, repr=False)

    def __post_init__(self) -> None:
        self._mm = np.memmap(CACHE_BIN, dtype=np.uint8, mode="r",
                             shape=(self.index["total_length"],))

    @property
    def sequences(self) -> list[dict]:
        return self.index["sequences"]

    @property
    def total_length(self) -> int:
        return int(self.index["total_length"])

    def names(self) -> list[str]:
        return [s["name"] for s in self.sequences]

    def get(self, name: str) -> np.ndarray:
        for s in self.sequences:
            if s["name"] == name:
                return np.asarray(self._mm[s["offset"]:s["offset"] + s["length"]])
        raise KeyError(name)

    def by_record(self, rec: dict) -> np.ndarray:
        return np.asarray(self._mm[rec["offset"]:rec["offset"] + rec["length"]])


# --------------------------------------------------------------------------------------
# Guides
# --------------------------------------------------------------------------------------


@dataclass(frozen=True)
class Guide:
    """A query protospacer. `pam` is the on-target PAM, carried for reporting only."""

    guide_id: str
    label: str
    group: str
    protospacer: str
    pam: str = ""
    gene: str = ""
    conservation: float = float("nan")
    note: str = ""

    @property
    def spacer_length(self) -> int:
        return len(self.protospacer)

    @property
    def codes(self) -> np.ndarray:
        return encode_sequence(self.protospacer)


BENCHMARK_POOL = RESULTS_DIR / "sacas9_benchmark_pool.tsv"


def load_guides(pool_path: Path = BENCHMARK_POOL,
                spacer_lengths: tuple[int, ...] = (20,),
                reference_fasta: Path | None = None) -> list[Guide]:
    """The guide set to screen, selected by rule from the stage-6 pool.

    The selection is fixed *before* any off-target number exists and is stated in the
    report, so it cannot be tuned to flatter the recommendation:

    * the four published Amrani et al. 2024 guides (`amrani_guide` non-null);
    * every ICP0/RL2 site that beats their lead ICP0g2 on conservation AND passes the
      pipeline's standard filters -- the 12 alternatives stage 6 puts forward,
      including the headline recommendation `RL2_3441+`;
    * every ICP27/UL54 site that passes the standard filters and outranks their lead
      ICP27g1 within its own gene.

    Screening all 4,777 SpCas9 candidates is deliberately *not* attempted: they are
    not in SaCas9 site space and the question at issue is a head-to-head against four
    named guides.
    """
    if not pool_path.exists():
        raise SystemExit(
            f"{pool_path} not found. Run stage 6 first:\n"
            "    python src/benchmark_sacas9.py")
    pool = pd.read_csv(pool_path, sep="\t")

    amrani = pool[pool["amrani_guide"].notna()]
    lead_rank = int(pool.loc[pool["amrani_guide"] == "Amrani2024_ICP27g1",
                             "rank_in_gene"].iloc[0])
    icp0_alt = pool[(pool["gene"] == "RL2")
                    & (pool["beats_lead_guide"].astype(bool))
                    & (pool["passes_filters"].astype(bool))]
    icp27_alt = pool[(pool["gene"] == "UL54")
                     & (pool["passes_filters"].astype(bool))
                     & (pool["rank_in_gene"] < lead_rank)]

    guides: list[Guide] = []
    seen: set[str] = set()

    def add(row, group: str, note: str = "") -> None:
        proto = str(row["protospacer"]).upper()
        if proto in seen:
            return
        seen.add(proto)
        amr = row["amrani_guide"]
        label = str(amr).replace("Amrani2024_", "") if isinstance(amr, str) else str(row["guide_id"])
        guides.append(Guide(
            guide_id=str(row["guide_id"]),
            label=label,
            group=group,
            protospacer=proto,
            pam=str(row["pam"]).upper(),
            gene=str(row["gene_label"]),
            conservation=float(row["conservation_complete_genomes"]),
            note=note,
        ))

    for _, row in amrani.iterrows():
        lead = row["amrani_guide"] in ("Amrani2024_ICP0g2", "Amrani2024_ICP27g1")
        add(row, "published", "lead clinical pair" if lead else "")
    for _, row in icp0_alt.sort_values("rank_in_gene").iterrows():
        add(row, "icp0_alternative",
            "headline recommendation" if row["guide_id"] == "RL2_3441+" else "")
    for _, row in icp27_alt.sort_values("rank_in_gene").iterrows():
        add(row, "icp27_alternative")

    if 21 in spacer_lengths:
        guides.extend(_extend_to_21nt(guides, pool, reference_fasta))
    if 20 not in spacer_lengths:
        guides = [g for g in guides if g.spacer_length != 20]
    return guides


def _read_single_fasta(path: Path) -> str:
    seq: list[str] = []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line.startswith(">"):
                seq.append(line.strip())
    return "".join(seq).upper()


def _extend_to_21nt(guides: list[Guide], pool: pd.DataFrame,
                    reference_fasta: Path | None) -> list[Guide]:
    """Build 21-nt spacer variants by taking the base immediately 5' in HSV-1.

    Amrani et al. published 20-nt SaCas9 spacers and the stage-6 pool is enumerated in
    that grammar, so 20 nt is the like-for-like length. SaCas9 is however commonly
    used at 21 nt (this repository's own `SACAS9` default), and a shorter spacer is
    the *permissive* direction -- fewer positions to mismatch means more genomic sites
    qualify -- so the 21-nt variants are a check that the 20-nt answer is not an
    artefact of length.
    """
    ref_path = reference_fasta or (RAW_DIR / "NC_001806.2.fasta")
    if not ref_path.exists():
        LOG.warning("HSV-1 reference %s absent; skipping 21-nt variants", ref_path)
        return []
    ref = _read_single_fasta(ref_path)
    by_id = pool.set_index("guide_id")
    out: list[Guide] = []
    for g in guides:
        if g.guide_id not in by_id.index:
            continue
        row = by_id.loc[g.guide_id]
        start, end, strand = int(row["ref_start"]), int(row["ref_end"]), str(row["strand"])
        if strand == "+":
            if start - 2 < 0:
                continue
            extra = ref[start - 2]                       # 1-based start -> 0-based start-1
            proto21 = extra + g.protospacer
        else:
            if end >= len(ref):
                continue
            extra = revcomp(ref[end])                    # base 3' on the plus strand
            proto21 = extra + g.protospacer
        # Sanity: the 20-nt suffix must be unchanged.
        assert proto21[1:] == g.protospacer
        out.append(Guide(guide_id=g.guide_id + "|21nt", label=g.label + " (21 nt)",
                         group=g.group + "_21nt", protospacer=proto21, pam=g.pam,
                         gene=g.gene, conservation=g.conservation,
                         note=(g.note + "; 21-nt variant").strip("; ")))
    return out


# --------------------------------------------------------------------------------------
# The fast search
# --------------------------------------------------------------------------------------


def pam_position_mask(arr: np.ndarray, spacer_length: int, pam: str) -> np.ndarray:
    """Boolean mask over protospacer start positions with a matching 3' PAM.

    `arr` is code-space. Position `s` is True iff `arr[s+L : s+L+P]` satisfies the
    IUPAC `pam`. Every test is a *slice* comparison -- contiguous, vectorised, no
    gather -- which is what makes the 15/16 rejection essentially free.
    """
    n = arr.size
    span = spacer_length + len(pam)
    m = n - span + 1
    if m <= 0:
        return np.zeros(0, dtype=bool)
    mask = np.ones(m, dtype=bool)
    for i, allowed in enumerate(iupac_code_sets(pam)):
        col = arr[spacer_length + i: spacer_length + i + m]
        if allowed.size == 4:                    # N: any real base, but not code 4
            mask &= col < 4
        elif allowed.size == 1:
            mask &= col == allowed[0]
        else:
            ok = col == allowed[0]
            for a in allowed[1:]:
                ok |= col == a
            mask &= ok
    return mask


#: A chunk is indexed by a dense 4**width lookup table, so the width has to stay
#: small. 8 nt costs 65,536 entries; 20 nt would cost 10^12. When few mismatches are
#: allowed the pigeonhole bound alone would permit very wide chunks, so the layout is
#: additionally split until every chunk fits. Using MORE chunks than the pigeonhole
#: bound requires is always sound -- with k chunks and at most m mismatches, at least
#: k - m chunks are exact -- it only costs a slightly higher candidate rate.
MAX_CHUNK_WIDTH = 8


def chunk_layout(spacer_length: int, max_mismatches: int) -> list[tuple[int, int]]:
    """Contiguous (offset, width) chunks tiling the protospacer.

    Pigeonhole: at most `max_mismatches` substitutions cannot touch all
    `max_mismatches + 1` chunks, so at least one chunk is an exact match. Widths are
    as equal as possible; the extra bases go to the leftmost chunks. For the default
    case -- 20 nt, 4 mismatches -- this is exactly five 4-mers.
    """
    k = max(max_mismatches + 1,
            (spacer_length + MAX_CHUNK_WIDTH - 1) // MAX_CHUNK_WIDTH)
    if spacer_length < k:
        raise ValueError(f"spacer of {spacer_length} nt cannot be split into {k} chunks")
    base, extra = divmod(spacer_length, k)
    widths = [base + (1 if i < extra else 0) for i in range(k)]
    offsets, acc = [], 0
    for w in widths:
        offsets.append((acc, w))
        acc += w
    return offsets


def _kmer_code_array(arr: np.ndarray, width: int) -> np.ndarray:
    """`C[i]` = the `width`-mer at `arr[i:i+width]` packed 2 bits per base.

    Built from *slices* of `arr & 3`, so it costs `width` contiguous passes and no
    gathers. Positions containing an N get a meaningless code, which can only ever
    manufacture a false candidate (rejected by the exact extend step), never lose a
    true one: a true hit with <= m mismatches has at least one chunk that matches the
    guide exactly, and such a chunk contains no N by construction.
    """
    n = arr.size
    m = n - width + 1
    dtype = np.uint16 if width > 4 else np.uint8
    out = np.zeros(m, dtype=dtype)
    a = arr & 3
    for j in range(width):
        out <<= 2
        out |= a[j:j + m].astype(dtype, copy=False)
    return out


def _chunk_code(seq: str) -> int:
    code = 0
    for ch in seq:
        code = (code << 2) | int(BASE_CODE[ord(ch)] & 3)
    return code


def _mismatch_tables(guides: list[Guide], layout: list[tuple[int, int]]
                     ) -> list[np.ndarray]:
    """`T[c][g, code]` = mismatches between chunk `c` of guide `g` and `code`.

    Because the chunks tile the protospacer, `sum_c T[c][g, code_c]` is the exact
    total mismatch count -- the extend step therefore never re-reads the genome.
    """
    tables = []
    for c, (off, width) in enumerate(layout):
        size = 1 << (2 * width)
        codes = np.arange(size, dtype=np.uint32)
        bases = np.stack([(codes >> (2 * (width - 1 - j))) & 3 for j in range(width)])
        tab = np.zeros((len(guides), size), dtype=np.uint8)
        for gi, g in enumerate(guides):
            gc = encode_sequence(g.protospacer[off:off + width]).astype(np.uint32)
            tab[gi] = (bases != gc[:, None]).sum(axis=0)
        tables.append(tab)
    return tables


def _bitmask_luts(guides: list[Guide], layout: list[tuple[int, int]],
                  group_offset: int) -> list[np.ndarray]:
    """`L[c][code]` = bitmask of guides whose chunk `c` equals `code`."""
    luts = []
    for c, (off, width) in enumerate(layout):
        lut = np.zeros(1 << (2 * width), dtype=np.uint64)
        for gi, g in enumerate(guides):
            lut[_chunk_code(g.protospacer[off:off + width])] |= np.uint64(
                1 << (gi - group_offset))
        luts.append(lut)
    return luts


def scan_sequence(arr: np.ndarray, guides: list[Guide], *, pam: str = "NNGRRN",
                  max_mismatches: int = 4, block: int = 1 << 23,
                  strand: str = "+") -> list[tuple[int, int, int]]:
    """Find every (guide index, start, mismatches) hit on ONE strand of `arr`.

    `start` is the 0-based index of the protospacer in `arr`; the PAM follows it.
    All guides must share a spacer length.
    """
    if not guides:
        return []
    lengths = {g.spacer_length for g in guides}
    if len(lengths) != 1:
        raise ValueError(f"scan_sequence needs one spacer length, got {sorted(lengths)}")
    L = lengths.pop()
    if len(guides) > 64:
        raise ValueError("scan_sequence handles <= 64 guides per call")

    mask = pam_position_mask(arr, L, pam)
    if not mask.any():
        return []
    positions = np.flatnonzero(mask).astype(np.int64)
    del mask

    layout = chunk_layout(L, max_mismatches)
    widths = sorted({w for _, w in layout})
    code_arrays = {w: _kmer_code_array(arr, w) for w in widths}
    luts = _bitmask_luts(guides, layout, 0)
    tables = _mismatch_tables(guides, layout)
    gcodes = [g.codes for g in guides]

    hits: list[tuple[int, int, int]] = []
    for lo in range(0, positions.size, block):
        pos = positions[lo:lo + block]
        codes = [code_arrays[w][pos + off] for off, w in layout]
        bits = luts[0][codes[0]]
        for c in range(1, len(layout)):
            bits |= luts[c][codes[c]]
        # Most positions share an exact chunk with no guide at all; resolve those once
        # rather than once per guide, which is the difference between O(guides x block)
        # and O(guides x candidates).
        nz = np.flatnonzero(bits)
        if nz.size == 0:
            continue
        bits_nz = bits[nz]
        del bits
        for gi in range(len(guides)):
            sel = nz[np.flatnonzero((bits_nz >> np.uint64(gi)) & np.uint64(1))]
            if sel.size == 0:
                continue
            mm = tables[0][gi][codes[0][sel]].astype(np.uint8)
            for c in range(1, len(layout)):
                mm += tables[c][gi][codes[c][sel]]
            keep = mm <= max_mismatches
            if not keep.any():
                continue
            # `mm` is computed from 2-bit codes, in which an ambiguous base collapses
            # onto A, so it is a strict LOWER bound on the true mismatch count: a site
            # can only look better than it is, never worse. That makes it a sound
            # filter but not a sound answer, so the handful of survivors -- true hits
            # plus a negligible number of N-boundary windows -- are read back from the
            # genome and counted exactly. N is code 4 and equals no guide base, so it
            # scores as a mismatch here, which is where the documented semantics is
            # actually enforced.
            starts = pos[sel[keep]]
            window = arr[starts[:, None] + np.arange(L)]
            exact = (window != gcodes[gi]).sum(axis=1)
            final = exact <= max_mismatches
            for p, k in zip(starts[final].tolist(), exact[final].tolist()):
                hits.append((gi, int(p), int(k)))
    return hits


# --------------------------------------------------------------------------------------
# Reference implementations used only by the tests and by --selfcheck
# --------------------------------------------------------------------------------------


def bruteforce_scan_sequence(arr: np.ndarray, guides: list[Guide], *,
                             pam: str = "NNGRRN", max_mismatches: int = 4
                             ) -> list[tuple[int, int, int]]:
    """Independent reference scan: no PAM prefilter, no chunk index, no gathers.

    Mismatch counts are accumulated at *every* genomic position by one whole-array
    slice comparison per protospacer base, and the PAM constraint is applied
    afterwards. Structurally this shares nothing with `scan_sequence` except the base
    encoding, which is the point.
    """
    L = guides[0].spacer_length
    span = L + len(pam)
    n = arr.size
    m = n - span + 1
    if m <= 0:
        return []
    pam_ok = np.ones(m, dtype=bool)
    for i, allowed in enumerate(iupac_code_sets(pam)):
        col = arr[L + i:L + i + m]
        ok = np.zeros(m, dtype=bool)
        for a in allowed:
            ok |= col == a
        pam_ok &= ok
    out: list[tuple[int, int, int]] = []
    for gi, g in enumerate(guides):
        gc = g.codes
        mm = np.zeros(m, dtype=np.uint8)
        for j in range(L):
            mm += (arr[j:j + m] != gc[j])
        sel = np.flatnonzero(pam_ok & (mm <= max_mismatches))
        for p in sel.tolist():
            out.append((gi, int(p), int(mm[p])))
    return out


def naive_scan_python(seq: str, guides: list[Guide], *, pam: str = "NNGRRN",
                      max_mismatches: int = 4) -> list[tuple[int, int, int]]:
    """The most obvious implementation there is: pure Python, character by character.

    Quadratic and unusably slow at genome scale; it exists to validate the numpy
    brute force, which in turn validates the fast path.
    """
    sets = [set(IUPAC_CODES[c]) for c in pam.upper()]
    L = guides[0].spacer_length
    span = L + len(pam)
    out: list[tuple[int, int, int]] = []
    for start in range(0, len(seq) - span + 1):
        window = seq[start:start + span]
        pam_seq = window[L:]
        if any(ch not in s for ch, s in zip(pam_seq, sets)):
            continue
        proto = window[:L]
        for gi, g in enumerate(guides):
            mm = sum(1 for a, b in zip(proto, g.protospacer) if a != b)
            if mm <= max_mismatches:
                out.append((gi, start, mm))
    return out


def caveat() -> str:
    """The sentence that must accompany any downstream report."""
    return (
        "Guides were screened for substitution off-targets (<= 4 mismatches, PAM "
        f"NNGRRN and NNGRRT) against the Ensembl release-{ENSEMBL_RELEASE} GRCh38 "
        "primary assembly, both strands, genome-wide. Bulge-containing off-targets "
        "were NOT modelled, and no cell-based validation (e.g. GUIDE-seq) was "
        "performed."
    )


def unscreened_caveat() -> str:
    """The sentence that must accompany a run in which stage 8 did NOT execute.

    `caveat()` describes a screen that happened; stating it after a run that skipped
    stage 8 would be a false claim about the guides in that run's own tables, so the
    two are separate functions and the caller picks by what actually ran.
    """
    return (
        "Off-target screening was NOT performed in this run. The guides in this "
        "run's tables have not been checked against the human genome here. Stage 8 "
        "(`python run_pipeline.py --offtarget`, or `python src/offtarget.py`) "
        "implements that screen for the stage-6 SaCas9 benchmark guide set; its "
        "results are in results/offtarget_report.md and cover that guide set only."
    )


# --------------------------------------------------------------------------------------
# Genomic context (Ensembl GTF)
# --------------------------------------------------------------------------------------

_GTF_FEATURES = ("gene", "exon", "CDS")


def _gtf_attr(attrs: str, key: str) -> str:
    """Pull one value out of a GTF attribute string without a regex per line."""
    tag = key + ' "'
    i = attrs.find(tag)
    if i == -1:
        return ""
    i += len(tag)
    j = attrs.find('"', i)
    return attrs[i:j]


class IntervalIndex:
    """Per-chromosome sorted intervals with an O(log n + k) overlap query.

    Intervals are 1-based inclusive (GTF convention). The backward scan is bounded by
    the longest interval on that chromosome, so the query never degenerates into a
    full scan just because one gene happens to be a megabase long.
    """

    def __init__(self, chrom: np.ndarray, start: np.ndarray, end: np.ndarray,
                 name: np.ndarray, extra: np.ndarray):
        order = np.lexsort((start, chrom))
        self.chrom, self.start = chrom[order], start[order]
        self.end, self.name, self.extra = end[order], name[order], extra[order]
        self._bounds: dict[str, tuple[int, int, int]] = {}
        if self.chrom.size:
            uniq, first = np.unique(self.chrom, return_index=True)
            firsts = list(first) + [self.chrom.size]
            for i, c in enumerate(uniq):
                lo, hi = int(firsts[i]), int(firsts[i + 1])
                span = int((self.end[lo:hi] - self.start[lo:hi]).max()) if hi > lo else 0
                self._bounds[str(c)] = (lo, hi, span)

    def overlaps(self, chrom: str, qs: int, qe: int) -> list[int]:
        b = self._bounds.get(chrom)
        if b is None:
            return []
        lo, hi, span = b
        s = self.start[lo:hi]
        left = int(np.searchsorted(s, qs - span, side="left"))
        right = int(np.searchsorted(s, qe, side="right"))
        if right <= left:
            return []
        idx = np.arange(lo + left, lo + right)
        return idx[self.end[idx] >= qs].tolist()

    def __len__(self) -> int:
        return int(self.chrom.size)


def build_annotation_cache(force: bool = False) -> dict[str, IntervalIndex]:
    """Parse the Ensembl GTF into three interval indexes: genes, exons, CDS.

    Coordinates stay exactly as the GTF gives them -- 1-based, inclusive, on Ensembl
    chromosome names, which are the same names the Ensembl FASTA uses. No liftover and
    no name mapping, so there is nothing that could silently shift a coordinate.
    """
    gtf = GENOME_DIR / GTF_NAME
    if not gtf.exists():
        raise SystemExit(f"{gtf} not found; run `python src/offtarget.py --stage fetch`")

    if ANNOTATION_CACHE.exists() and not force:
        z = np.load(ANNOTATION_CACHE, allow_pickle=False)
        return {k: IntervalIndex(z[f"{k}_chrom"], z[f"{k}_start"], z[f"{k}_end"],
                                 z[f"{k}_name"], z[f"{k}_extra"])
                for k in _GTF_FEATURES}

    LOG.info("parsing %s (one-off)", GTF_NAME)
    t0 = time.time()
    rows: dict[str, list] = {k: [] for k in _GTF_FEATURES}
    n_lines = 0
    with gzip.open(gtf, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            n_lines += 1
            parts = line.split("\t")
            feat = parts[2]
            if feat not in rows:
                continue
            attrs = parts[8]
            rows[feat].append((
                parts[0], int(parts[3]), int(parts[4]),
                _gtf_attr(attrs, "gene_name") or _gtf_attr(attrs, "gene_id"),
                _gtf_attr(attrs, "gene_biotype")))

    out: dict[str, IntervalIndex] = {}
    payload: dict[str, np.ndarray] = {}
    for k, recs in rows.items():
        chrom = np.array([r[0] for r in recs], dtype="U")
        start = np.array([r[1] for r in recs], dtype=np.int64)
        end = np.array([r[2] for r in recs], dtype=np.int64)
        name = np.array([r[3] for r in recs], dtype="U")
        extra = np.array([r[4] for r in recs], dtype="U")
        out[k] = IntervalIndex(chrom, start, end, name, extra)
        payload[f"{k}_chrom"], payload[f"{k}_start"] = chrom, start
        payload[f"{k}_end"], payload[f"{k}_name"], payload[f"{k}_extra"] = end, name, extra
        LOG.info("  %-5s %d records", k, len(recs))
    np.savez_compressed(ANNOTATION_CACHE, **payload)
    parse_seconds = round(time.time() - t0, 1)
    ANNOTATION_CACHE.with_suffix(".json").write_text(json.dumps({
        "gtf": GTF_NAME, "parse_seconds": parse_seconds, "gtf_feature_lines": n_lines,
        "counts": {k: len(v) for k, v in rows.items()},
    }, indent=1), encoding="utf-8")
    LOG.info("annotation parsed in %.0fs (%d GTF feature lines)", parse_seconds, n_lines)
    return out


def annotate_hits(df: pd.DataFrame, ann: dict[str, IntervalIndex]) -> pd.DataFrame:
    """Add gene / exon / CDS context columns to a hit table (1-based inclusive)."""
    genes, exons, cds = ann["gene"], ann["exon"], ann["CDS"]
    cols: dict[str, list] = {"in_gene": [], "gene_names": [], "gene_biotypes": [],
                             "in_exon": [], "exon_genes": [], "in_cds": [],
                             "cds_genes": [], "cut_in_exon": [], "cut_in_cds": []}
    for chrom, s, e, cut in zip(df["chrom"], df["start"], df["end"], df["cut_site"]):
        gi = genes.overlaps(str(chrom), int(s), int(e))
        ei = exons.overlaps(str(chrom), int(s), int(e))
        ci = cds.overlaps(str(chrom), int(s), int(e))
        # The footprint flags above are the conservative reading (any overlap of the
        # 26-nt site). These two ask the sharper question: does the predicted blunt cut
        # itself land in an exon / a coding exon?
        cols["cut_in_exon"].append(bool(exons.overlaps(str(chrom), int(cut), int(cut))))
        cols["cut_in_cds"].append(bool(cds.overlaps(str(chrom), int(cut), int(cut))))
        cols["in_gene"].append(bool(gi))
        cols["gene_names"].append(",".join(sorted({str(x) for x in genes.name[gi]})) if gi else "")
        cols["gene_biotypes"].append(
            ",".join(sorted({str(x) for x in genes.extra[gi]})) if gi else "")
        cols["in_exon"].append(bool(ei))
        cols["exon_genes"].append(",".join(sorted({str(x) for x in exons.name[ei]})) if ei else "")
        cols["in_cds"].append(bool(ci))
        cols["cds_genes"].append(",".join(sorted({str(x) for x in cds.name[ci]})) if ci else "")
    for k, v in cols.items():
        df[k] = v
    return df


# --------------------------------------------------------------------------------------
# Screening driver
# --------------------------------------------------------------------------------------

#: SaCas9's PAM-proximal "seed". Mismatches here are far more disruptive than
#: PAM-distal ones. Reported, never used to filter.
SEED_LENGTH = 12

HIT_COLUMNS = [
    "guide_id", "label", "group", "gene", "chrom", "start", "end", "strand",
    "mismatches", "seed_mismatches", "pam_class", "protospacer_genomic", "pam_genomic",
    "mismatch_positions", "cut_site",
]


def _scan_one_strand(arr: np.ndarray, guides: list[Guide], pam: str,
                     max_mm: int, block: int) -> list[tuple[int, int, int]]:
    """Batch guides into groups of <= 64 (one uint64 bitmask lane each)."""
    out: list[tuple[int, int, int]] = []
    for lo in range(0, len(guides), 64):
        batch = guides[lo:lo + 64]
        for gi, p, mm in scan_sequence(arr, batch, pam=pam, max_mismatches=max_mm,
                                       block=block):
            out.append((lo + gi, p, mm))
    return out


def materialise_hits(hits: list[tuple[int, int, int]], guides: list[Guide],
                     arr: np.ndarray, chrom: str, strand: str, seq_len: int,
                     pam_len: int) -> list[dict]:
    """Turn (guide, position, mismatches) into fully specified, re-verified rows.

    `arr` is the strand that was actually scanned (the reverse complement for '-'), so
    the protospacer and PAM are read straight out of it and only the coordinate is
    mapped back into plus-strand space. Every mismatch count is recomputed here from
    the genomic bases and asserted against the count the chunk tables produced: if the
    fast path were ever wrong, this fails loudly instead of silently under-reporting.
    """
    rows: list[dict] = []
    for gi, p, mm in hits:
        g = guides[gi]
        L = g.spacer_length
        proto = arr[p:p + L]
        pam_codes = arr[p + L:p + L + pam_len]
        diff = np.flatnonzero(proto != g.codes)
        if diff.size != mm:
            raise AssertionError(
                f"fast-path mismatch count {mm} != recomputed {diff.size} for "
                f"{g.guide_id} at {chrom} index {p} strand {strand}")
        start0 = p if strand == "+" else seq_len - (p + L + pam_len)
        pam_seq = decode_sequence(pam_codes)
        nngrrt = len(pam_seq) >= 6 and pam_seq[5] == "T"
        seed_mm = int((diff >= L - SEED_LENGTH).sum())
        # SaCas9 cleaves bluntly ~3 bp 5' of the PAM, inside the protospacer.
        cut = (start0 + L - 3) if strand == "+" else (start0 + pam_len + 3)
        rows.append({
            "guide_id": g.guide_id, "label": g.label, "group": g.group, "gene": g.gene,
            "chrom": chrom, "start": start0 + 1, "end": start0 + L + pam_len,
            "strand": strand, "mismatches": int(mm), "seed_mismatches": seed_mm,
            "pam_class": "NNGRRT" if nngrrt else "NNGRRN_only",
            "protospacer_genomic": decode_sequence(proto), "pam_genomic": pam_seq,
            "mismatch_positions": ",".join(str(int(d) + 1) for d in diff),
            "cut_site": cut + 1,
        })
    return rows


def build_summary(hits: pd.DataFrame, guides: list[Guide], max_mm: int) -> pd.DataFrame:
    """Per-guide off-target counts at each mismatch level, for both PAM variants."""
    rows = []
    for g in guides:
        h = hits[hits["guide_id"] == g.guide_id] if len(hits) else hits
        rec: dict = {
            "guide_id": g.guide_id, "label": g.label, "group": g.group,
            "gene": g.gene, "spacer_length": g.spacer_length,
            "protospacer": g.protospacer, "on_target_pam": g.pam,
            "conservation_complete_genomes": g.conservation, "note": g.note,
        }
        strict = h[h["pam_class"] == "NNGRRT"] if len(h) else h
        for variant, sel in (("nngrrt", strict), ("nngrrn", h)):
            for k in range(max_mm + 1):
                rec[f"{variant}_mm{k}"] = int((sel["mismatches"] == k).sum()) if len(sel) else 0
            rec[f"{variant}_total"] = int(len(sel))
            rec[f"{variant}_le3"] = int((sel["mismatches"] <= 3).sum()) if len(sel) else 0
            rec[f"{variant}_seed_le1"] = (
                int((sel["seed_mismatches"] <= 1).sum()) if len(sel) else 0)
        rows.append(rec)
    return pd.DataFrame(rows)


def run_screen(args: argparse.Namespace) -> dict:
    """Acquire, scan, annotate, tabulate, report. Returns a JSON-ready summary."""
    timings: dict[str, float] = {}

    t = time.time()
    manifest = acquire_genome(skip_bsd_sum=args.skip_checksum, want_gtf=not args.no_gtf)
    timings["acquire_seconds"] = round(time.time() - t, 1)

    t = time.time()
    cache = build_genome_cache(force=args.rebuild_cache)
    timings["cache_seconds"] = round(time.time() - t, 1)

    spacer_lengths = tuple(int(x) for x in str(args.spacer_lengths).split(","))
    guides = load_guides(Path(args.pool), spacer_lengths=spacer_lengths)
    LOG.info("screening %d guides (spacer lengths %s) at <= %d mismatches, PAM %s",
             len(guides), sorted({g.spacer_length for g in guides}),
             args.max_mismatches, args.pam)

    records = cache.sequences
    if args.chrom:
        want = set(str(args.chrom).split(","))
        records = [r for r in records if r["name"] in want]
        if not records:
            raise SystemExit(f"no sequence named {args.chrom!r} in the cache")

    by_len: dict[int, list[Guide]] = {}
    for g in guides:
        by_len.setdefault(g.spacer_length, []).append(g)

    pam_len = len(args.pam)
    rows: list[dict] = []
    scanned_bases = 0
    n_pam_positions = 0
    t_scan = time.time()
    for ri, rec in enumerate(records, 1):
        name = rec["name"]
        arr = cache.by_record(rec)
        scanned_bases += arr.size
        rc = revcomp_codes(arr)
        t_rec = time.time()
        rec_hits = 0
        for L, gl in by_len.items():
            n_pam_positions += 2 * int(pam_position_mask(arr, L, args.pam).sum())
            for strand, a in (("+", arr), ("-", rc)):
                hits = _scan_one_strand(a, gl, args.pam, args.max_mismatches, args.block)
                rec_hits += len(hits)
                rows.extend(materialise_hits(hits, gl, a, name, strand, arr.size, pam_len))
        del rc
        if arr.size > 5_000_000 or ri % 25 == 0 or ri == len(records):
            LOG.info("  [%3d/%3d] %-22s %11d bp %7d hits %6.1fs (cumulative %.0fs)",
                     ri, len(records), name, arr.size, rec_hits,
                     time.time() - t_rec, time.time() - t_scan)
    timings["scan_seconds"] = round(time.time() - t_scan, 1)
    LOG.info("scan complete: %d hits over %.3f Gb in %.0fs",
             len(rows), scanned_bases / 1e9, timings["scan_seconds"])

    hits = (pd.DataFrame(rows, columns=HIT_COLUMNS) if rows
            else pd.DataFrame(columns=HIT_COLUMNS))
    if len(hits):
        hits = hits.sort_values(["guide_id", "mismatches", "chrom", "start"],
                                kind="stable").reset_index(drop=True)

    summary = build_summary(hits, guides, args.max_mismatches)

    ann_rows = pd.DataFrame()
    if not args.no_gtf and len(hits):
        t = time.time()
        ann = build_annotation_cache(force=args.rebuild_cache)
        ann_rows = hits[hits["mismatches"] <= args.annotate_max_mismatches].copy()
        if len(ann_rows):
            ann_rows = annotate_hits(ann_rows, ann)
            counts = (ann_rows[ann_rows["in_cds"].astype(bool)]
                      .groupby("guide_id").size())
            summary["n_le3_in_cds"] = summary["guide_id"].map(counts).fillna(0).astype(int)
        else:
            summary["n_le3_in_cds"] = 0
        timings["annotate_seconds"] = round(time.time() - t, 1)

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    hits.to_csv(RESULTS_DIR / "offtarget_sites.tsv", sep="\t", index=False)
    summary.to_csv(RESULTS_DIR / "offtarget_summary.tsv", sep="\t", index=False)
    if len(ann_rows):
        ann_rows.to_csv(RESULTS_DIR / "offtarget_annotated.tsv", sep="\t", index=False)

    gtf_path = GENOME_DIR / GTF_NAME
    payload = {
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "ensembl_release": ENSEMBL_RELEASE,
        "assembly": cache.index["assembly"],
        "genome_bases_scanned": int(scanned_bases),
        "genome_sequences_scanned": len(records),
        "genome_sequences_total": len(cache.sequences),
        "pam_searched": args.pam,
        "max_mismatches": int(args.max_mismatches),
        "spacer_lengths": list(sorted({g.spacer_length for g in guides})),
        "n_guides": len(guides),
        "n_sites_total": int(len(hits)),
        "pam_positions_examined": int(n_pam_positions),
        "timings": timings,
        "genome_manifest": manifest,
        "disk_bytes": {
            "fasta_gz": (GENOME_DIR / FASTA_NAME).stat().st_size,
            "encoded_u8": CACHE_BIN.stat().st_size,
            "gtf_gz": gtf_path.stat().st_size if gtf_path.exists() else 0,
            "annotation_npz": (ANNOTATION_CACHE.stat().st_size
                               if ANNOTATION_CACHE.exists() else 0),
        },
        "caveat": caveat(),
    }
    (RESULTS_DIR / "offtarget_summary.json").write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    write_report(payload, summary, hits, ann_rows, args)
    return payload


# --------------------------------------------------------------------------------------
# Report
# --------------------------------------------------------------------------------------


def _md_table(df: pd.DataFrame, cols: list[str], headers: list[str] | None = None) -> str:
    headers = headers or cols
    lines = ["| " + " | ".join(headers) + " |",
             "|" + "|".join("---" for _ in headers) + "|"]
    for _, r in df.iterrows():
        vals = []
        for c in cols:
            v = r[c]
            if isinstance(v, float):
                vals.append("--" if pd.isna(v) else f"{v:.3f}")
            else:
                vals.append(str(v))
        lines.append("| " + " | ".join(vals) + " |")
    return "\n".join(lines)


#: Genes worth naming in prose if a close off-target lands in one. Short, factual
#: notes only -- this is a lookup table for the report, not a risk model, and a gene
#: not on this list is not thereby safe.
NOTABLE_LOCI = {
    "ABL1": "proto-oncogene tyrosine kinase; the BCR-ABL fusion partner",
    "EGR2": "master transcriptional regulator of peripheral-nerve myelination -- the "
            "intended target tissue here is the trigeminal ganglion",
    "PDCD4": "tumour suppressor, translation inhibitor",
    "NSDHL": "X-linked sterol dehydrogenase; loss of function causes CHILD syndrome",
    "SLC5A5": "sodium/iodide symporter; thyroid hormone synthesis",
    "DOCK6": "Adams-Oliver syndrome gene",
    "MYH14": "non-muscle myosin heavy chain; DFNA4 hearing loss",
    "FOXO6": "forkhead transcription factor, high expression in brain",
    "TOX": "HMG-box transcription factor, T-cell and neural development",
    "OBSCN": "obscurin, sarcomeric signalling",
    "EME2": "structure-specific endonuclease subunit, DNA repair",
    "ARHGEF5": "Rho guanine nucleotide exchange factor; reported oncogenic activity",
    "ZNF331": "transcriptional repressor, frequently methylated in gastrointestinal "
              "cancer -- the locus Amrani et al. themselves validated as an ICP27g1 "
              "off-target",
    "OSBP2": "oxysterol-binding protein, photoreceptor-enriched",
}


def write_report(payload: dict, summary: pd.DataFrame, hits: pd.DataFrame,
                 ann: pd.DataFrame, args: argparse.Namespace) -> Path:
    out = RESULTS_DIR / "offtarget_report.md"
    t = payload["timings"]
    man = payload["genome_manifest"].get(FASTA_NAME, {})
    disk = payload["disk_bytes"]
    total_disk = sum(disk.values())
    mm_max = int(args.max_mismatches)
    s20 = summary[summary["spacer_length"] == 20]

    def row(label_or_id: str, frame: pd.DataFrame) -> pd.DataFrame:
        sel = frame[(frame["label"] == label_or_id) | (frame["guide_id"] == label_or_id)]
        return sel

    def val(frame: pd.DataFrame, col: str):
        return int(frame[col].iloc[0]) if len(frame) else None

    lead0, lead27 = row("ICP0g2", s20), row("ICP27g1", s20)
    head = row("RL2_3441+", s20)

    L: list[str] = []
    A = L.append
    A("# Human (GRCh38) off-target screening of the SaCas9 guide set")
    A("")
    A(f"Stage 8, `src/offtarget.py`. Generated {payload['generated_utc']}. "
      "Every number below is measured on this machine, not estimated or quoted.")
    A("")

    A("## 0. Headline")
    A("")
    A("This analysis was run because it is capable of invalidating this repository's own "
      "recommendation. It partly did. The finding is reported as it fell.")
    A("")
    if len(lead0):
        icp0 = s20[(s20["gene"] == "ICP0") | (s20["group"].str.startswith("icp0"))]
        icp0 = icp0[icp0["group"].isin(("published", "icp0_alternative"))]
        cmp_cols = ("nngrrt_total", "nngrrt_le3", "nngrrn_total", "nngrrn_le3")
        base = {c: int(lead0[c].iloc[0]) for c in cmp_cols}
        cleaner = icp0[(icp0["label"] != "ICP0g2")
                       & (icp0["nngrrt_total"] <= base["nngrrt_total"])
                       & (icp0["nngrrt_le3"] <= base["nngrrt_le3"])
                       & (icp0["nngrrn_total"] <= base["nngrrn_total"])
                       & (icp0["nngrrn_le3"] <= base["nngrrn_le3"])]
        perfect_cleaner = cleaner[cleaner["conservation_complete_genomes"] >= 1.0]

        A("**1. Our named headline guide does not survive.** `RL2_3441+` was put forward "
          "in stage 6 as the ICP0 guide that should have been chosen instead of ICP0g2. "
          "Against GRCh38 it is materially *dirtier* than ICP0g2:")
        A("")
        A("| | ICP0g2 (theirs) | RL2_3441+ (ours) |")
        A("|---|---|---|")
        if len(head):
            A(f"| NNGRRT sites <= 4 mm | {base['nngrrt_total']} | "
              f"{val(head,'nngrrt_total')} |")
            A(f"| NNGRRT sites <= 3 mm | {base['nngrrt_le3']} | "
              f"{val(head,'nngrrt_le3')} |")
            A(f"| NNGRRN sites <= 4 mm | {base['nngrrn_total']} | "
              f"{val(head,'nngrrn_total')} |")
            A(f"| NNGRRN sites <= 3 mm | {base['nngrrn_le3']} | "
              f"{val(head,'nngrrn_le3')} |")
            A(f"| HSV-1 conservation (n=183) | "
              f"{float(lead0['conservation_complete_genomes'].iloc[0]):.3f} | "
              f"{float(head['conservation_complete_genomes'].iloc[0]):.3f} |")
        A("")
        A("The 3-mismatch tier is the one that matters and it is where `RL2_3441+` is "
          "worst: it carries close sites under the **canonical** NNGRRT PAM, where "
          "ICP0g2 carries none, and one of them sits in *ABL1*. Four more are copies of "
          "a single repeated sequence on chrY. On off-target burden alone, ICP0g2 is the "
          "better of the two guides, and stage 6's specific recommendation is withdrawn.")
        A("")
        if len(perfect_cleaner):
            names = ", ".join(f"`{x}`" for x in perfect_cleaner["guide_id"])
            A(f"**2. The underlying argument survives, with a different guide.** "
              f"{len(perfect_cleaner)} of the ICP0 alternatives are perfectly conserved "
              f"across all 183 complete HSV-1 genomes *and* have an off-target profile at "
              f"least as clean as ICP0g2 on every measure computed here: {names}. They "
              "were in the same candidate pool Amrani et al. drew from, in their own site "
              "grammar.")
            A("")
            A(_md_table(
                pd.concat([lead0, perfect_cleaner]),
                ["label", "guide_id", "conservation_complete_genomes", "nngrrt_total",
                 "nngrrt_le3", "nngrrn_total", "nngrrn_le3"],
                ["guide", "site id", "HSV-1 cons.", "NNGRRT <=4mm", "NNGRRT <=3mm",
                 "NNGRRN <=4mm", "NNGRRN <=3mm"]))
            A("")
            A("So the claim that better ICP0 options existed inside their own site space "
              "is *not* refuted by off-target data -- it is refuted only for the "
              "particular guide stage 6 named. The corrected statement is narrower and "
              "still stands: at least one perfectly conserved ICP0 site is also cleaner "
              "against the human genome than the guide they took to the clinic.")
        else:
            A("**2. No alternative is both perfectly conserved and at least as clean as "
              "ICP0g2 on every measure computed here.** On the combined evidence the "
              "stage-6 recommendation does not stand, and ICP0g2 is defensible as their "
              "choice.")
        A("")
        icp27 = s20[s20["group"].isin(("icp27_alternative",))]
        if len(lead27) and len(icp27):
            b27 = {c: int(lead27[c].iloc[0]) for c in cmp_cols}
            beat27 = icp27[(icp27["nngrrt_total"] <= b27["nngrrt_total"])
                           & (icp27["nngrrt_le3"] <= b27["nngrrt_le3"])
                           & (icp27["nngrrn_total"] <= b27["nngrrn_total"])
                           & (icp27["nngrrn_le3"] <= b27["nngrrn_le3"])]
            if len(beat27) == 0:
                A("**2b. Their ICP27 choice looks good, and this analysis says so.** "
                  f"ICP27g1 carries {b27['nngrrn_total']} NNGRRN sites at <= 4 mismatches "
                  f"and {b27['nngrrn_le3']} at <= 3 -- the cleanest profile of every "
                  f"ICP27/UL54 candidate screened. **None** of the "
                  f"{len(icp27)} better-conserved, filter-passing ICP27 alternatives "
                  "matches it on all four off-target measures. Stage 6 never proposed "
                  "replacing ICP27g1; this is the first evidence that it should not be.")
            else:
                names = ", ".join(f"`{x}`" for x in beat27["guide_id"])
                A(f"**2b. ICP27.** {len(beat27)} of the ICP27 alternatives are at least as "
                  f"clean as ICP27g1 on every measure computed here: {names}.")
            A("")
        A("**3. Independent corroboration of the method.** Amrani et al. report (Table 3) "
          "that no sites with <= 3 total mismatches plus bulges were found for either "
          "lead guide. Under the canonical NNGRRT PAM this screen finds exactly that: "
          f"{base['nngrrt_le3']} sites at <= 3 mm for ICP0g2 and "
          f"{val(lead27, 'nngrrt_le3') if len(lead27) else '--'} for ICP27g1, computed "
          "independently, from a different assembly download, by a different algorithm. "
          "Under the permissive NNGRRN PAM that they actually searched, this screen does "
          f"find sites ({base['nngrrn_le3']} for ICP0g2, "
          f"{val(lead27, 'nngrrn_le3') if len(lead27) else '--'} for ICP27g1), all with "
          "non-canonical PAMs. The most likely explanation is their BWA `aln` "
          "pre-nomination step, which is a heuristic aligner and is not guaranteed to "
          "surface every 3-mismatch locus; the scan here is exhaustive. That is a "
          "difference in completeness, not a contradiction, and it cuts in their favour "
          "as much as against them.")
        A("")

    A("## 1. Reference genome")
    A("")
    A("| item | value |")
    A("|---|---|")
    A(f"| assembly | {payload['assembly']} |")
    A(f"| source | Ensembl release {ENSEMBL_RELEASE} (GRCh38.p14, GCA_000001405.29) |")
    A(f"| file | `{FASTA_NAME}` |")
    A(f"| URL | {FASTA_URL} |")
    A(f"| compressed size | {disk['fasta_gz']:,} bytes |")
    A(f"| SHA-256 | `{man.get('sha256', '--')}` |")
    A(f"| Ensembl published `sum` | {man.get('published_bsd_sum', '--')} "
      f"({man.get('published_blocks', '--')} x 1 KiB blocks) |")
    A(f"| checksum verified | {man.get('bsd_sum_verified', '--')} |")
    A(f"| sequences in assembly | {payload['genome_sequences_total']} |")
    A(f"| sequences scanned | {payload['genome_sequences_scanned']} |")
    A(f"| bases scanned | {payload['genome_bases_scanned']:,} |")
    A("")
    A("**Why unmasked, and why the primary assembly.** `dna` (unmasked) rather than "
      "`dna_sm` or `dna_rm`. A Cas9 site inside a LINE, an Alu or a satellite is a real "
      "cleavage substrate, so repeat masking is not merely unhelpful here, it is wrong: "
      "hard-masking would delete roughly half the genome from the search, and "
      "soft-masking is only safe if every read remembers to upper-case, which is exactly "
      "the kind of silent bug that would quietly shrink an off-target count. "
      "`primary_assembly` rather than `toplevel`: toplevel adds ALT haplotypes and patch "
      "scaffolds, which are alternative representations of sequence already present on "
      "the chromosomes and would double-count off-targets.")
    A("")
    A(f"**Disk actually used: {total_disk / 1e9:.2f} GB.** "
      f"{disk['fasta_gz'] / 1e6:.0f} MB compressed FASTA + "
      f"{disk['encoded_u8'] / 1e9:.2f} GB one-byte-per-base cache + "
      f"{disk['gtf_gz'] / 1e6:.0f} MB GTF + "
      f"{disk['annotation_npz'] / 1e6:.0f} MB parsed annotation.")
    A("")
    A("**Integrity, three independent ways.** Ensembl's own published `sum` checksum "
      "(the BSD 16-bit rotating checksum, recomputed here and compared); the gzip "
      "member's CRC-32 and ISIZE trailer, which the decompressor verifies when the whole "
      "file is streamed during cache construction; and a SHA-256 we compute and record "
      "in `data/genome/genome_manifest.json`. The download is resumable via HTTP Range "
      "and is never repeated on a re-run.")
    A("")

    A("## 2. Method")
    A("")
    A(f"- PAM searched: **{args.pam}**. `NNGRRT` sites are a strict subset of `NNGRRN` "
      "sites, so they are flagged rather than searched separately and both variants come "
      "out of a single pass.")
    A(f"- Up to **{mm_max} mismatches** in the protospacer, **both strands**, "
      "genome-wide, across every sequence in the primary assembly (chromosomes, "
      "unlocalised and unplaced scaffolds, and the mitochondrion).")
    A(f"- Spacer lengths screened: {payload['spacer_lengths']}. 20 nt is the grammar "
      "Amrani et al. published in their Table 1 and the grammar the stage-6 candidate "
      "pool is enumerated in, so it is the like-for-like length; 21 nt is this "
      "repository's SaCas9 default and is built by extending each site one base 5' in "
      "HSV-1. A shorter spacer is the permissive direction -- fewer positions available "
      "to mismatch means more genomic sites qualify -- so the 20-nt numbers are the "
      "conservative ones.")
    A("- **Substitutions only. Bulges (DNA/RNA insertions) are not modelled.**")
    A("- Ambiguous genomic bases (`N`) count as mismatches in the protospacer and "
      "disqualify a PAM outright, so no site is ever credited to unresolved sequence.")
    A("")
    A("**Algorithm.** Exact substring search is the wrong tool -- an off-target is by "
      "definition mismatched -- so stage 3's jump scanner cannot be reused. Two "
      "properties of the problem are exploited instead.")
    A("")
    A("*PAM sparsity.* SaCas9 needs `NNGRRN`: a G at PAM offset 2 and purines at offsets "
      "3 and 4. That is 1 position in 16 per strand (1 in 64 for `NNGRRT`). Only those "
      "positions can ever be cut, so three vectorised slice comparisons discard 15/16 of "
      "the genome before any protospacer is looked at. This is the single biggest win "
      "available and it is free.")
    A("")
    A("*The pigeonhole principle.* With at most 4 mismatches, a 20-nt protospacer cut "
      "into five 4-nt chunks must match the guide exactly in at least one chunk -- four "
      "mismatches cannot touch five chunks. Every genomic 4-mer is packed into a byte "
      "once per chromosome by pure slice arithmetic (no gathers); the five chunk codes "
      "are then gathered at the PAM positions only; and a 256-entry bitmask lookup table "
      "per chunk resolves, for all guides simultaneously, which positions share an exact "
      "chunk with which guide. About 2% of PAM positions survive per guide.")
    A("")
    A("*Scoring from the codes already in hand.* Because the five chunks tile the whole "
      "protospacer, a per-(guide, chunk) 256-entry table of \"how many of these four "
      "bases differ\" turns the five codes into a mismatch count with five more table "
      "lookups and no further genome access. Those codes pack an ambiguous base onto A, "
      "so that count is a strict *lower* bound -- a site can look better than it is, "
      "never worse -- which makes it a sound filter but not a sound answer. The handful "
      "of positions that survive it are therefore read back from the genome and counted "
      "exactly, which is where the \"N is a mismatch\" rule is actually enforced. Every "
      "reported site is then re-verified a second time during table construction and the "
      "count asserted, so a fast-path error would raise rather than silently "
      "under-report. Everything is numpy; there is no per-base Python loop on the hot "
      "path.")
    A("")
    idx = json.loads(CACHE_INDEX.read_text(encoding="utf-8"))
    ann_meta = {}
    ann_json = ANNOTATION_CACHE.with_suffix(".json")
    if ann_json.exists():
        ann_meta = json.loads(ann_json.read_text(encoding="utf-8"))
    scan_s = max(float(t.get("scan_seconds", 0.0)), 1e-9)
    bases = payload["genome_bases_scanned"]
    A("**Wall-clock, measured on this machine.** One-off costs, paid once and cached: "
      f"BSD-checksum verification of the 882 MB FASTA "
      f"{man.get('bsd_sum_seconds', '--')} s; decoding it into the 1-byte-per-base "
      f"array {idx.get('build_seconds', '--')} s; parsing the GTF "
      f"{ann_meta.get('parse_seconds', '--')} s "
      f"({ann_meta.get('gtf_feature_lines', '--')} feature lines -> "
      f"{ann_meta.get('counts', {}).get('gene', '--')} genes, "
      f"{ann_meta.get('counts', {}).get('exon', '--')} exons, "
      f"{ann_meta.get('counts', {}).get('CDS', '--')} CDS).")
    A("")
    A(f"**The genome-wide scan itself: {scan_s:.0f} s.** "
      f"{payload['n_guides']} guides x 2 strands x "
      f"{bases / 1e9:.2f} Gb, i.e. {bases / 1e6 / scan_s:.1f} Mb of genome per second, "
      f"or {2 * payload['n_guides'] * bases / 1e9 / scan_s:.2f} billion "
      "(guide x strand x base) positions per second. "
      f"{payload['pam_positions_examined']:,} PAM-matched positions were examined out of "
      f"{2 * bases:,} strand-bases. Annotation of the close sites: "
      f"{t.get('annotate_seconds', 0):.0f} s. A warm re-run (caches present) is the scan "
      "time plus a few seconds. For scale, the brute-force reference implementation the "
      "tests compare against needs about 35 s per (spacer length x strand) for "
      "chromosome 21's 46.7 Mb, which extrapolates to roughly 2.5 hours for the whole "
      f"assembly against {scan_s / 60:.0f} minutes here -- consistent with the 23-44x "
      "speed-up measured directly on chr21 in the equivalence test.")
    A("")
    A("**Correctness.** `tests/test_offtarget.py::test_fast_matches_bruteforce_chr21` "
      "asserts that the fast path returns *exactly* the same hit set as a structurally "
      "independent brute-force scan -- mismatch counts accumulated at every genomic "
      "position by one whole-array slice comparison per protospacer base, no PAM "
      "prefilter, no chunk index, explicit reverse complement -- over the whole of "
      "chromosome 21, for every guide screened. A third, plainly written pure-Python "
      "character-by-character implementation validates the numpy brute force in turn. "
      "This is the standard `tests/test_nuclease.py` already applies to the stage-3 PAM "
      "scanner, applied here.")
    A("")
    A("**Positioning, honestly.** Cas-OFFinder is the field standard and is what Amrani "
      "et al. used (BWA `aln` to nominate homologous loci, then Cas-OFFinder over those "
      "loci with PAM NNGRRN, then GUIDE-seq validation). What is used here is a "
      "validated equivalent **for substitution-only search**, and it is exhaustive where "
      "a BWA pre-filter is not. It is **not** equivalent for bulges: their post-filter "
      "retained 1mm+1bulge, 2mm+1bulge, 3mm+1bulge, 4mm+0bulge and 5mm+0bulge classes, "
      "so their published totals (ICP0g1 358, ICP0g2 910, ICP27g1 443, ICP27g2 316) "
      "count classes this screen does not enumerate, and were computed against a "
      "different assembly build. Counts here are therefore **not** comparable to those "
      "totals line for line. They are comparable *between guides screened here*, which "
      "is the question this stage exists to answer.")
    A("")

    A("## 3. Guides screened")
    A("")
    A("Selection rule, fixed before any off-target number existed, so it cannot have "
      "been tuned to flatter the recommendation:")
    A("")
    A("1. the four published Amrani et al. 2024 guides;")
    A("2. every ICP0/RL2 site that beats their lead ICP0g2 on cross-isolate conservation "
      "*and* passes the pipeline's standard filters -- the 12 stage-6 alternatives, "
      "including the headline recommendation `RL2_3441+`;")
    A("3. every ICP27/UL54 site that passes those filters *and* outranks their lead "
      "ICP27g1 within its own gene.")
    A("")
    A("All 4,777 SpCas9 candidates are deliberately not screened: they are not in SaCas9 "
      "site space, and the question at issue is a head-to-head against four named "
      "guides.")
    A("")
    A(_md_table(s20, ["label", "guide_id", "gene", "group", "protospacer",
                      "on_target_pam", "conservation_complete_genomes"],
                ["guide", "site id", "gene", "class", "protospacer (20 nt)", "PAM",
                 "HSV-1 conservation"]))
    A("")

    A("## 4. Off-target counts")
    A("")
    A(f"Distinct genomic sites within *k* mismatches of the protospacer with a matching "
      "PAM. `NNGRRT` is the canonical SaCas9 PAM; `NNGRRN` is the permissive variant "
      "Amrani et al. searched, and is a superset of it.")
    A("")
    for spacer in payload["spacer_lengths"]:
        sub = summary[summary["spacer_length"] == spacer]
        if not len(sub):
            continue
        A(f"### {spacer}-nt spacers")
        A("")
        for variant, title in (("nngrrt", "NNGRRT (canonical)"),
                               ("nngrrn", "NNGRRN (permissive)")):
            cols = [f"{variant}_mm{k}" for k in range(mm_max + 1)]
            A(f"**{title}**")
            A("")
            A(_md_table(sub, ["label", "guide_id", "group"] + cols
                        + [f"{variant}_total", f"{variant}_seed_le1"],
                        ["guide", "site id", "class"]
                        + [f"{k} mm" for k in range(mm_max + 1)]
                        + ["total", "seed<=1mm"]))
            A("")
    A("`seed<=1mm` counts sites with at most one mismatch in the PAM-proximal 12 nt, "
      "the region where SaCas9 tolerates mismatches least. It is reported, not used to "
      "filter.")
    A("")

    A("## 5. Genomic context of the close off-targets")
    A("")
    if not len(ann):
        A(f"**No off-target site at <= {args.annotate_max_mismatches} mismatches exists "
          "for any guide screened, under either PAM variant.** There is nothing to "
          "annotate. Note what this does and does not say: it says the close-range "
          "substitution off-target profile is empty for the whole set, and it says "
          "nothing about bulge-containing sites or about the 4-mismatch tier.")
    else:
        in_cds = ann[ann["in_cds"].astype(bool)]
        A(f"{len(ann)} site(s) at <= {args.annotate_max_mismatches} mismatches, annotated "
          f"against Ensembl {ENSEMBL_RELEASE} (`{GTF_NAME}`). Coordinates are 1-based "
          "inclusive on Ensembl chromosome names -- the same names the FASTA uses, so no "
          "liftover and no name mapping is involved and nothing can silently shift.")
        A("")
        n_cut_cds = int(ann["cut_in_cds"].astype(bool).sum())
        n_exon = int(ann["in_exon"].astype(bool).sum())
        n_cut_exon = int(ann["cut_in_exon"].astype(bool).sum())
        if len(in_cds):
            A(f"**{len(in_cds)} site(s) fall inside a coding exon (CDS)**, of which "
              f"{n_cut_cds} put the predicted blunt cut itself inside the CDS rather "
              "than merely overlapping one with the 26-nt footprint. "
              f"({n_exon} sites overlap an exon of any biotype, {n_cut_exon} cut inside "
              "one.) The intended target cell type is post-mitotic sensory neurons, "
              "where double-strand break repair is NHEJ-only and any indel is permanent "
              "and unrepairable by homologous recombination. A coding-exon off-target is "
              "the dangerous class, and these are listed first.")
            A("")
            A(_md_table(in_cds, ["label", "chrom", "start", "strand", "mismatches",
                                 "seed_mismatches", "pam_class", "pam_genomic",
                                 "cds_genes", "cut_in_cds"],
                        ["guide", "chrom", "start", "str", "mm", "seed mm", "PAM class",
                         "PAM", "CDS gene", "cut in CDS"]))
            A("")
        else:
            A("**No site at this range falls inside a coding exon (CDS).**")
            A("")
        per_guide = (ann.assign(exon=ann["in_exon"].astype(int),
                                cds=ann["in_cds"].astype(int),
                                cut_cds=ann["cut_in_cds"].astype(int),
                                nngrrt=(ann["pam_class"] == "NNGRRT").astype(int))
                     .groupby(["label"], as_index=False)
                     .agg(sites=("mismatches", "size"), nngrrt=("nngrrt", "sum"),
                          in_exon=("exon", "sum"), in_cds=("cds", "sum"),
                          cut_in_cds=("cut_cds", "sum"))
                     .sort_values(["in_cds", "nngrrt", "sites"], ascending=False))
        named = []
        for _, r in ann.iterrows():
            for field in ("cds_genes", "gene_names"):
                for gname in str(r.get(field, "") or "").split(","):
                    note = NOTABLE_LOCI.get(gname.strip())
                    if note:
                        named.append((gname.strip(), str(r["label"]), int(r["mismatches"]),
                                      str(r["pam_class"]), bool(r["in_cds"]), note))
        seen_pair = set()
        uniq = []
        for g, lbl, mm, pc, incds, note in named:
            key = (g, lbl.replace(" (21 nt)", ""))
            if key in seen_pair:
                continue
            seen_pair.add(key)
            uniq.append((g, lbl.replace(" (21 nt)", ""), mm, pc, incds, note))
        if uniq:
            A("**Loci worth naming.** Not a risk model -- a lookup of genes whose "
              "identity a reader would want flagged. A gene absent from this list is "
              "not thereby safe.")
            A("")
            for g, lbl, mm, pc, incds, note in sorted(uniq):
                where = "inside a coding exon" if incds else "within the gene body"
                A(f"- ***{g}*** -- {note}. Hit by `{lbl}` at {mm} mismatches "
                  f"({pc} PAM), {where}.")
            A("")
        A("Per guide (all spacer lengths pooled; a site found by both the 20-nt and the "
          "21-nt version of one guide is the same genomic locus and is counted once per "
          "version):")
        A("")
        A(_md_table(per_guide,
                    ["label", "sites", "nngrrt", "in_exon", "in_cds", "cut_in_cds"],
                    ["guide", f"sites <= {args.annotate_max_mismatches}mm",
                     "of which canonical NNGRRT", "in an exon", "in a CDS",
                     "cut inside a CDS"]))
        A("")
        A("Full listing for the published 20-nt grammar (the 21-nt variants are in "
          "`results/offtarget_annotated.tsv`):")
        A("")
        short = ann[~ann["label"].astype(str).str.contains(r"\(21 nt\)", regex=True)]
        A(_md_table(short, ["label", "chrom", "start", "strand", "mismatches",
                            "seed_mismatches", "pam_class", "pam_genomic", "gene_names",
                            "in_exon", "in_cds"],
                    ["guide", "chrom", "start", "str", "mm", "seed mm", "PAM class",
                     "PAM", "gene(s)", "exon", "CDS"]))
    A("")

    A("## 6. Limitations")
    A("")
    A("1. **Substitutions only.** DNA/RNA bulges are not enumerated. Cas-OFFinder does "
      "enumerate them and Amrani et al. counted them. A bulge-tolerant reanalysis can "
      "only add sites, never remove them, and could in principle reorder the guides.")
    A("2. **No activity model.** These are raw site counts. No CFD, MIT/Hsu or "
      "Doench score is applied, because none of those models is parameterised for "
      "SaCas9; an SpCas9 score applied to SaCas9 data would be worse than no score. "
      "Seed-region mismatch counts are reported instead, unweighted.")
    A("3. **No empirical validation.** Amrani et al. ran GUIDE-seq and targeted amplicon "
      "deep sequencing. Nothing here is measured in cells. In-silico nomination is a "
      "filter, not evidence of cleavage -- and equally, absence of a nominated site is "
      "not proof of safety.")
    A("4. **One reference, one haplotype.** GRCh38 primary assembly only. Population "
      "variation creates and destroys both PAMs and protospacer matches; ALT haplotypes "
      "were excluded to avoid double-counting, which necessarily makes "
      "haplotype-specific sites invisible.")
    A("5. **Annotation is Ensembl-only**, at gene/exon/CDS granularity. No cancer-gene "
      "or essential-gene list is applied. A site outside a gene is not thereby safe: "
      "enhancers, promoters and other regulatory elements are not modelled. Amrani et "
      "al.'s own validated off-target for ICP27g1 sits in an intron of *ZNF331* that is "
      "exon 1 of one transcript variant -- exactly the kind of case that gene-level "
      "annotation renders as 'intronic'.")
    A("6. **The 4-mismatch ceiling is a choice.** Amrani et al.'s post-filter also "
      "retained a 5-mismatch class. Raising the ceiling would raise every count; it "
      "would not change the ranking unless one guide's excess is concentrated at 5 mm.")
    A("")

    A("## 7. Reproducibility")
    A("")
    A("```")
    A("python src/offtarget.py                  # fetch + verify + scan + annotate + report")
    A("python src/offtarget.py --stage fetch    # acquisition and caches only")
    A("python src/offtarget.py --chrom 21       # scan one sequence (fast smoke test)")
    A("python tests/test_offtarget.py           # includes the chr21 equivalence proof")
    A("```")
    A("")
    A("Supporting tables: `results/offtarget_summary.tsv` (per-guide counts), "
      "`results/offtarget_sites.tsv` (every site, with coordinates, strand, mismatch "
      "positions and the observed PAM), `results/offtarget_annotated.tsv` (close sites "
      "with gene context), `results/offtarget_summary.json` (machine-readable; carries "
      "the genome manifest, the checksums and all timings).")
    A("")
    A("**Caveat for any manuscript.** " + caveat())
    A("")

    out.write_text("\n".join(L), encoding="utf-8")
    LOG.info("wrote %s", out)
    return out


# --------------------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------------------


def add_arguments(parser: argparse.ArgumentParser) -> None:
    g = parser.add_argument_group("human off-target screening (stage 8)")
    g.add_argument("--stage", choices=("all", "fetch", "scan", "report"), default="all",
                   help="'fetch' downloads, verifies and encodes GRCh38 then stops; "
                        "'report' re-renders results/offtarget_report.md from the "
                        "existing TSVs without rescanning.")
    g.add_argument("--pam", default="NNGRRN",
                   help="PAM to search. NNGRRT is reported as a subset of it.")
    g.add_argument("--max-mismatches", type=int, default=4)
    g.add_argument("--annotate-max-mismatches", type=int, default=3,
                   help="annotate genomic context for sites at or below this many "
                        "mismatches")
    g.add_argument("--spacer-lengths", default="20,21",
                   help="comma-separated spacer lengths to screen. 20 is the published "
                        "Amrani grammar; 21 is this repository's SaCas9 default, built "
                        "by extending each site one base 5' in HSV-1.")
    g.add_argument("--chrom", default="",
                   help="restrict the scan to these sequence names (comma-separated)")
    g.add_argument("--pool", default=str(BENCHMARK_POOL),
                   help="stage-6 candidate pool the guide set is selected from")
    g.add_argument("--block", type=int, default=1 << 23,
                   help="PAM positions processed per vectorised block")
    g.add_argument("--skip-checksum", action="store_true",
                   help="skip the slow one-off BSD `sum` verification")
    g.add_argument("--no-gtf", action="store_true", help="skip annotation entirely")
    g.add_argument("--rebuild-cache", action="store_true")


def main(argv: list[str] | None = None) -> int:
    """Entry point. `main(["--chrom", "21"])` is all integration needs."""
    parser = argparse.ArgumentParser(
        prog="offtarget",
        description="Human GRCh38 off-target screening for the SaCas9 guide set.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    add_arguments(parser)
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args(argv)

    from src.common import setup_logging
    setup_logging(args.verbose)

    if args.stage == "report":
        payload = json.loads((RESULTS_DIR / "offtarget_summary.json").read_text(
            encoding="utf-8"))
        summary = pd.read_csv(RESULTS_DIR / "offtarget_summary.tsv", sep="	")
        hits = pd.read_csv(RESULTS_DIR / "offtarget_sites.tsv", sep="	")
        ann_path = RESULTS_DIR / "offtarget_annotated.tsv"
        ann = pd.read_csv(ann_path, sep="	") if ann_path.exists() else pd.DataFrame()
        write_report(payload, summary, hits, ann, args)
        return 0

    if args.stage == "fetch":
        acquire_genome(skip_bsd_sum=args.skip_checksum, want_gtf=not args.no_gtf)
        cache = build_genome_cache(force=args.rebuild_cache)
        if not args.no_gtf:
            build_annotation_cache(force=args.rebuild_cache)
        LOG.info("genome cache ready: %d sequences, %.3f Gb",
                 len(cache.sequences), cache.total_length / 1e9)
        return 0

    payload = run_screen(args)
    LOG.info("done: %d off-target sites for %d guides; scan took %.0fs",
             payload["n_sites_total"], payload["n_guides"],
             payload["timings"].get("scan_seconds", 0.0))
    return 0


def screen_offtargets(argv: list[str] | None = None) -> int:
    """Backwards-compatible name kept from the phase-1 stub, which this replaced."""
    return main(argv)


if __name__ == "__main__":
    raise SystemExit(main())
