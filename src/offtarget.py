"""STUB -- human off-target screening. NOT IMPLEMENTED IN PHASE 1.

===============================================================================
 TODO (PHASE 2): human genome off-target screening
===============================================================================

This module is a deliberate, clearly marked placeholder. Phase 1 of this project
scores only *on-target conservation across HSV isolates*. No candidate guide
produced by this pipeline has been screened against the human genome, and none of
them should be described as "safe", "specific", or "therapeutically viable" until
this module exists.

Any results table emitted by report.py therefore carries an implicit caveat that
MUST be stated in any manuscript:

    "Guides were not screened for off-target activity against the human genome;
     conservation ranking reflects on-target coverage across HSV-1 isolates only."

--------------------------------------------------------------------------------
What phase 2 needs to do
--------------------------------------------------------------------------------
1. Obtain a human reference assembly (GRCh38 / T2T-CHM13). This is ~3.1 Gb and
   cannot be handled by the exact-substring approach used in conservation.py --
   off-target sites are by definition MISmatched, so exact matching is the wrong
   tool.

2. Enumerate near-matches allowing up to 3-4 mismatches in the protospacer, with
   PAM-proximal ("seed", ~12 nt adjacent to the PAM) mismatches weighted far more
   heavily than PAM-distal ones. Include NAG as a permissive PAM, not just NGG.

3. Score each near-match with a published model (e.g. CFD, or Hsu/Zhang MIT
   specificity score) and aggregate per guide into a specificity score.

4. Flag guides whose off-target sites fall in exons, cancer-associated genes, or
   essential genes.

--------------------------------------------------------------------------------
Constraint conflict to resolve before implementing
--------------------------------------------------------------------------------
The pure-Python / no-external-binary constraint that makes phase 1 clean does NOT
survive contact with a 3.1 Gb mismatch search. Realistic options, in order of
preference:

  (a) Cas-OFFinder or CRISPRitz as an external binary (breaks the constraint but is
      the field standard and is what reviewers will expect).
  (b) A pure-Python seed-and-extend index over GRCh38 built from the ~2 x 10^8
      NGG/CCN sites. Feasible but memory-hungry and slow; would need on-disk
      index shards.
  (c) A remote service (e.g. CRISPOR / CHOPCHOP) queried over HTTP. No local
      dependency, but introduces a network dependency and a rate limit, and the
      exact tool version becomes hard to pin for reproducibility.

Decision deferred to phase 2. Document whichever is chosen in the methods section.
===============================================================================
"""

from __future__ import annotations

__all__ = ["screen_offtargets", "OFFTARGET_IMPLEMENTED"]

OFFTARGET_IMPLEMENTED = False

_MESSAGE = (
    "Human off-target screening is not implemented (phase 1 scope). "
    "See src/offtarget.py for the phase-2 plan. Guides from this pipeline are "
    "ranked by cross-isolate conservation ONLY and have not been checked for "
    "off-target activity in the human genome."
)


def screen_offtargets(*_args, **_kwargs):
    """Intentionally unimplemented. See module docstring."""
    raise NotImplementedError(_MESSAGE)


def caveat() -> str:
    """Return the caveat string that must accompany any downstream report."""
    return _MESSAGE
