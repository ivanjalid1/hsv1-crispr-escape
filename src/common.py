"""Shared utilities: paths, logging, NCBI Entrez configuration, rate limiting, sequence helpers.

PRIVACY: no email address is ever hardcoded in this repository. The NCBI contact
address is read exclusively from the environment variable NCBI_EMAIL, and the
optional API key from NCBI_API_KEY.
"""

from __future__ import annotations

import hashlib
import logging
import os
import sys
import threading
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, TypeVar

# --------------------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
REF_DIR = RAW_DIR / "reference"
META_DIR = RAW_DIR / "metadata"
RESULTS_DIR = PROJECT_ROOT / "results"

DEFAULT_MANIFEST = DATA_DIR / "manifest.tsv"
DEFAULT_CANDIDATES = RESULTS_DIR / "guides_candidates.tsv"
DEFAULT_CONSERVATION = RESULTS_DIR / "conservation.tsv"
DEFAULT_RANKED = RESULTS_DIR / "guides_ranked.tsv"
DEFAULT_RUNLOG = RESULTS_DIR / "run_log.json"

TOOL_NAME = "hsv-crispr-conservation"


def ensure_dirs() -> None:
    for d in (DATA_DIR, RAW_DIR, REF_DIR, META_DIR, RESULTS_DIR):
        d.mkdir(parents=True, exist_ok=True)


# --------------------------------------------------------------------------------------
# Logging
# --------------------------------------------------------------------------------------


def setup_logging(verbose: bool = False) -> logging.Logger:
    level = logging.DEBUG if verbose else logging.INFO
    root = logging.getLogger()
    if not root.handlers:
        handler = logging.StreamHandler(sys.stderr)
        handler.setFormatter(
            logging.Formatter("%(asctime)s %(levelname)-7s %(name)-18s %(message)s",
                              datefmt="%H:%M:%S")
        )
        root.addHandler(handler)
    root.setLevel(level)
    return root


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


# --------------------------------------------------------------------------------------
# NCBI Entrez configuration
# --------------------------------------------------------------------------------------

_ENV_HELP = """
NCBI Entrez requires a contact email address so that NCBI can reach you if your
queries cause problems. This pipeline deliberately does NOT hardcode one.

Set it before running:

    PowerShell : $env:NCBI_EMAIL = "you@example.org"
    cmd.exe    : set NCBI_EMAIL=you@example.org
    bash       : export NCBI_EMAIL=you@example.org

Optionally, an NCBI API key raises your rate limit from 3 to 10 requests/second:

    PowerShell : $env:NCBI_API_KEY = "<your key>"

Get a key at https://www.ncbi.nlm.nih.gov/account/settings/
""".strip()


class MissingCredentialsError(RuntimeError):
    pass


def configure_entrez():
    """Configure Bio.Entrez from the environment. Raises if NCBI_EMAIL is unset."""
    from Bio import Entrez

    email = (os.environ.get("NCBI_EMAIL") or "").strip()
    if not email:
        raise MissingCredentialsError(
            "Environment variable NCBI_EMAIL is not set.\n\n" + _ENV_HELP
        )
    if "@" not in email or email.startswith("@") or email.endswith("@"):
        raise MissingCredentialsError(
            f"NCBI_EMAIL does not look like an email address: {email!r}\n\n" + _ENV_HELP
        )

    Entrez.email = email
    Entrez.tool = TOOL_NAME
    api_key = (os.environ.get("NCBI_API_KEY") or "").strip()
    if api_key:
        Entrez.api_key = api_key
    return Entrez


def entrez_rate_limit_interval() -> float:
    """NCBI permits 3 requests/second without an API key, 10 with one.

    We stay comfortably under the limit.
    """
    return 0.12 if (os.environ.get("NCBI_API_KEY") or "").strip() else 0.40


class RateLimiter:
    """Simple thread-safe minimum-interval rate limiter."""

    def __init__(self, min_interval: float):
        self.min_interval = min_interval
        self._last = 0.0
        self._lock = threading.Lock()

    def wait(self) -> None:
        with self._lock:
            now = time.monotonic()
            delta = now - self._last
            if delta < self.min_interval:
                time.sleep(self.min_interval - delta)
            self._last = time.monotonic()


T = TypeVar("T")


def with_retries(
    fn: Callable[[], T],
    *,
    attempts: int = 5,
    base_delay: float = 2.0,
    limiter: RateLimiter | None = None,
    description: str = "NCBI request",
    logger: logging.Logger | None = None,
) -> T:
    """Run `fn` with deterministic exponential backoff.

    Backoff is deterministic (no jitter) so that runs are reproducible in structure;
    only wall-clock timing varies.
    """
    log = logger or logging.getLogger("ncbi")
    last_exc: BaseException | None = None
    for attempt in range(1, attempts + 1):
        if limiter is not None:
            limiter.wait()
        try:
            return fn()
        except KeyboardInterrupt:
            raise
        except Exception as exc:  # noqa: BLE001 - network layer is genuinely broad
            last_exc = exc
            if attempt == attempts:
                break
            delay = base_delay * (2 ** (attempt - 1))
            log.warning(
                "%s failed (attempt %d/%d): %s: %s -- retrying in %.1fs",
                description, attempt, attempts, type(exc).__name__, exc, delay,
            )
            time.sleep(delay)
    raise RuntimeError(f"{description} failed after {attempts} attempts") from last_exc


# --------------------------------------------------------------------------------------
# Sequence helpers
# --------------------------------------------------------------------------------------

_COMPLEMENT = str.maketrans("ACGTUNRYSWKMBDHVacgtunryswkmbdhv",
                            "TGCAANYRSWMKVHDBtgcaanyrswmkvhdb")

UNAMBIGUOUS = frozenset("ACGT")


def revcomp(seq: str) -> str:
    """Reverse complement of an uppercase (or mixed-case) nucleotide string."""
    return seq.translate(_COMPLEMENT)[::-1]


def gc_fraction(seq: str) -> float:
    if not seq:
        return 0.0
    s = seq.upper()
    return (s.count("G") + s.count("C")) / len(s)


def max_homopolymer_run(seq: str) -> int:
    """Length of the longest run of a single identical character."""
    if not seq:
        return 0
    best = run = 1
    prev = seq[0]
    for ch in seq[1:]:
        if ch == prev:
            run += 1
            if run > best:
                best = run
        else:
            run = 1
            prev = ch
    return best


def has_polyt(seq: str, run: int = 4) -> bool:
    """TTTT (or longer) terminates RNA polymerase III transcription and kills sgRNA
    expression from a U6/H1 promoter."""
    return ("T" * run) in seq.upper()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()
