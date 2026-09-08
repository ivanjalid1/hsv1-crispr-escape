# The reconciled recommendation: one ICP0 guide, chosen on three axes jointly

Stages 6, 7 and 8 of this pipeline each proposed a *different* replacement for Amrani
et al. 2024's ICP0g2 guide, because each optimised a different axis:

| stage | module | axis optimised | proposed |
|---|---|---|---|
| 6 | `src/benchmark_sacas9.py` | cross-isolate conservation | `RL2_3441+` |
| 7 | `src/escape.py` | multiplex escape probability | `RL2_5080+` |
| 8 | `src/offtarget.py` | human GRCh38 off-target burden | `RL2_5335+`, and it **refuted `RL2_3441+`** |

Three recommendations is not a recommendation. This document produces one, chosen on
all three axes together, states the selection rule before applying it, and records
what is withdrawn and why.

---

## 0. The recommendation

**Keep `ICP27g1` (`UL54_115156+`) exactly as published. Replace `ICP0g2`
(`RL2_4496+`) with `RL2_5335+`.** One guide changes.

`RL2_5335+` — `GCGTCGGAGTGGAACAGCCT` + `CTGGAT` (20 nt spacer, `NNGRRT` PAM, plus
strand, reference positions 5,335–5,360 of `NC_001806`, present in **both** ICP0
repeat copies).

| | ICP0g2 + ICP27g1 (published lead pair) | **RL2_5335+ + ICP27g1** |
|---|---|---|
| conservation, 183 complete genomes | 0.847 | **1.000** |
| conservation, gene-level corpus | 0.898 (400 records) | **1.000 (429 records)** |
| conservation, sub-genomic tier only | 0.952 (21 records) | **1.000 (36 records)** |
| **joint** conservation of the pair | 0.8251 (151/183) | **0.9781 (179/183)** |
| P(escape) per exposed viral genome | 9.720e-04 | **2.5145e-04** (3.87× lower) |
| GRCh38 `NNGRRT` sites ≤ 4 mm | 5 | **3** |
| GRCh38 `NNGRRT` sites ≤ 3 mm | 0 | **0** |
| GRCh38 `NNGRRN` sites ≤ 4 mm | 56 | **41** |
| GRCh38 `NNGRRN` sites ≤ 3 mm | 5 | **2** |
| coding-exon hits ≤ 3 mm (either PAM) | 0 | **0** |
| 21-nt variant, `NNGRRT` ≤ 4 mm / `NNGRRN` ≤ 4 mm | 0 / 15 | **0 / 6** |
| GC (spacer) | 0.60 | 0.65 |
| local GC, 200 bp | **0.840** | **0.725** |
| longest homopolymer / poly-T | 2 / no | 2 / no |
| in both ICP0 repeat copies | yes | yes |

Every number above is measured in this repository. Nothing is quoted from the paper
except the identity of their guides.

**Two things about this answer are not tidy, and both are stated rather than
smoothed over.**

1. **Escape probability contributes no information whatsoever to the choice between
   the ICP0 candidates.** Section 4 shows they are *exactly* tied.
2. **No candidate dominates on every measure.** The Pareto front has three members
   and the choice needs one explicitly stated priority. Section 5.

---

## 1. The candidate pool

The pool is every ICP0/RL2 site in Amrani et al.'s own site grammar — 20 nt spacer +
`NNGRRT` PAM, inside the annotated RL2 CDS of `NC_001806` — that stage 8 actually
screened against the human genome. That is the 12 stage-6 alternatives (sites that
beat ICP0g2 on conservation *and* pass the standard filters) plus both published ICP0
guides: 14 of the gene's 46 `NNGRRT` sites. Sites outside this set have no off-target
measurement, so they cannot be compared on all three axes and are not ranked here.

`ICP27g1` is held fixed. Stage 8 examined that assumption and vindicated it:
ICP27g1 has the cleanest off-target profile of every ICP27/UL54 candidate screened
(24 `NNGRRN` sites ≤ 4 mm, 1 at ≤ 3 mm, 0 coding-exon hits), and **none** of the 9
better-conserved filter-passing ICP27 alternatives matches it on all four off-target
measures. The obvious "better" swap, `UL54_114833-` (conservation 1.000, and it would
take the pair's joint conservation to 1.000), carries 33 `NNGRRT` sites at ≤ 4 mm,
4 at ≤ 3 mm and 2 coding-exon hits. Changing both guides buys 4 more isolates at the
cost of a materially dirtier ICP27 guide, so the recommendation is the single swap.

---

## 2. The joint table

Every measure, every screened ICP0 candidate, ordered by stage-6 rank within the
gene. Conservation columns are higher-is-better; escape and off-target columns are
lower-is-better.

Column key: `cons183` = fraction of the 183 complete genomes containing the exact
site; `cons_gene` / `n_gene` = conservation over the stage-5 gene-level corpus and
the number of anchor-bracketed records covering the site; `cons_sub` / `n_sub` = the
same restricted to the sub-genomic tier, the only tier independent of the complete
genomes; `q` = the stage-7 per-site repair-escape probability; `P(esc)` and
`joint cons` = paired with ICP27g1; `T` = `NNGRRT`, `N` = `NNGRRN`, `seed` = sites
with ≤ 1 mismatch in the PAM-proximal 12 nt; `CDS` = sites at ≤ 3 mm inside a coding
exon; `21nt` columns repeat the screen for the 21-nt spacer variant.

| guide | pub. | rank | cons183 | cons_gene | n_gene | cons_sub | n_sub | copies | GC | localGC | homo | filters | q | joint cons | P(esc) | T ≤3 | T ≤4 | T seed | N ≤3 | N ≤4 | N seed | CDS | 21nt T ≤3/≤4 | 21nt N ≤3/≤4 | 21nt CDS |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **RL2_5335+** | | 1 | 1.0000 | 1.0000 | 429 | 1.0000 | 36 | 2 | 0.65 | 0.725 | 2 | pass | 0.009486 | 0.9781 | 2.515e-04 | **0** | **3** | 1 | 2 | **41** | 19 | 0 | 0 / 0 | 0 / 6 | 0 |
| RL2_3364+ | | 2 | 1.0000 | 0.9788 | 425 | 1.0000 | 32 | 2 | 0.75 | 0.710 | 2 | pass | 0.009486 | 0.9781 | 2.515e-04 | 0 | 5 | 1 | **1** | 54 | 21 | 0 | 0 / 2 | 0 / 20 | 0 |
| RL2_5080+ | | 3 | 1.0000 | 0.9907 | 429 | 0.9722 | 36 | 2 | 0.65 | 0.690 | 3 | pass | 0.009486 | 0.9781 | 2.515e-04 | 0 | 10 | 5 | **1** | 61 | 31 | 0 | 0 / 2 | 1 / 18 | 0 |
| RL2_3441+ | | 4 | 1.0000 | 1.0000 | 425 | 1.0000 | 32 | 2 | 0.50 | **0.635** | 4 | pass | 0.009486 | 0.9781 | 2.515e-04 | **6** | **19** | 7 | 9 | 78 | 17 | 0 | 1 / 11 | 2 / 28 | 0 |
| RL2_3196- | | 7 | 0.9945 | 0.9857 | 421 | 0.8929 | 28 | 2 | 0.75 | 0.740 | 2 | pass | 0.009486 | 0.9727 | 2.772e-04 | 1 | 11 | 1 | 4 | 70 | 9 | **1** | 0 / 6 | 2 / 34 | 0 |
| RL2_5167- | | 8 | 0.9945 | 1.0000 | 429 | 1.0000 | 36 | 2 | 0.55 | 0.680 | 4 | pass | 0.009486 | 0.9727 | 2.772e-04 | 0 | 23 | 4 | 3 | 90 | 16 | **1** | 0 / 8 | 2 / 26 | **1** |
| RL2_5395- | | 10 | 0.9891 | 1.0000 | 429 | 1.0000 | 36 | 2 | 0.75 | 0.695 | 3 | pass | 0.009486 | 0.9672 | 3.029e-04 | 0 | 13 | 4 | 5 | 89 | 25 | 0 | 0 / 3 | 3 / 26 | 0 |
| RL2_5319+ | ICP0g1 | 12 | 0.9836 | 0.9835 | 424 | 1.0000 | 31 | 2 | 0.80 | 0.735 | 5 | **fail** | 0.009506 | 0.9617 | 3.292e-04 | 1 | 4 | 0 | 1 | 19 | 5 | 0 | 0 / 2 | 0 / 4 | 0 |
| RL2_3568- | | 18 | 0.9672 | 0.9835 | 423 | 0.9667 | 30 | 2 | 0.70 | 0.685 | 3 | pass | 0.009506 | 0.9454 | 4.064e-04 | 0 | 28 | 1 | 38 | 349 | 8 | 0 | 0 / 8 | 3 / 77 | 0 |
| RL2_3576- | | 22 | 0.9508 | 0.9409 | 423 | 0.9667 | 30 | 2 | 0.70 | 0.675 | 2 | pass | 0.092824 | 0.9290 | 2.673e-03 | 2 | 26 | 4 | 12 | 151 | 36 | **5** | 1 / 5 | 3 / 37 | 0 |
| RL2_4810- | | 23 | 0.9508 | 0.9481 | 424 | 0.9697 | 33 | 2 | 0.75 | 0.815 | 2 | pass | 0.114751 | 0.9290 | 3.249e-03 | 0 | 18 | 2 | 4 | 97 | 7 | 0 | 0 / 4 | 2 / 32 | 0 |
| RL2_4659+ | | 25 | 0.9454 | 0.9423 | 416 | 0.9630 | 27 | 2 | 0.70 | 0.810 | 2 | pass | 0.009486 | 0.9235 | 5.088e-04 | 3 | 24 | 6 | 8 | 146 | 33 | 0 | 1 / 3 | 2 / 23 | 0 |
| RL2_3515+ | | 30 | 0.8907 | 0.8759 | 419 | 1.0000 | 26 | 2 | 0.65 | 0.635 | 3 | pass | 0.009486 | 0.8689 | 7.661e-04 | 0 | 4 | 0 | 3 | 38 | 7 | **2** | 0 / 1 | 1 / 10 | 0 |
| RL2_4496+ | **ICP0g2** | 31 | 0.8470 | 0.8975 | 400 | 0.9524 | 21 | 2 | 0.60 | **0.840** | 2 | pass | 0.009486 | 0.8251 | 9.720e-04 | 0 | 5 | 1 | 5 | 56 | 8 | 0 | 0 / 0 | 1 / 15 | 0 |

This table is **generated, not transcribed**: `python src/reconcile.py` joins the
three stages on `guide_id` and writes it to `results/recommendation_table.tsv`.
Sources are `results/sacas9_benchmark_pool.tsv` (conservation, GC, filters, repeat
copies), `results/escape_site_parameters.tsv` (`q`), `results/offtarget_summary.tsv`
and `results/offtarget_annotated.tsv` (off-target), and the stage-7 `EscapeModel`
re-run over the same presence matrix for the paired `P(escape)` and
joint-conservation columns (section 4).

---

## 3. The selection rule, stated before it is applied

The rule is fixed here, in this order, and then applied mechanically in section 5.
It was not tuned to produce a particular winner; section 6 reports what a defensible
*different* rule would have produced.

**Step 1 — hard gates.** A candidate is eliminated outright if it fails any of:

* **G1.** Conservation = 1.000 across all 183 complete genomes. *Why a gate and not a
  score: the whole criticism of the published selection is that a `>70%` threshold
  admits 69 of 81 sites and therefore barely discriminates. Anything less than
  perfect leaves known isolates uncovered when perfect options exist.*
* **G2.** Conservation ≥ 0.95 on the stage-5 gene-level corpus. *The second, partly
  independent denominator; a candidate that is perfect only on the corpus that
  enumerated it has not been checked.*
* **G3.** Passes the pipeline's standard filters — no poly-T terminator, longest
  homopolymer ≤ 4, GC in 0.35–0.75. *A U6/H1-driven sgRNA containing `TTTT` is
  transcribed truncated; this is a hard exclusion, not a penalty.*
* **G4.** Present in **both** ICP0 repeat copies. *The three-DSB architecture is the
  paper's own design rationale; a single-copy site silently changes the construct.*
* **G5.** **Zero** GRCh38 off-target sites at ≤ 3 mismatches under the canonical
  `NNGRRT` PAM. *This is the tier Amrani et al. themselves report as empty for both
  lead guides, and `NNGRRT` is the PAM SaCas9 cleaves efficiently and the PAM all
  four published on-target sites use. A candidate that introduces a near-cognate site
  where the published guide had none is not an improvement, whatever its
  conservation.*
* **G6.** Zero coding-exon hits at ≤ 3 mismatches under either PAM, at either spacer
  length.

**Step 2 — primary discriminator among survivors: total `NNGRRT` off-target burden at
≤ 4 mismatches**, lower is better. *Rationale, stated as the judgement it is: the
`NNGRRN` search space is deliberately permissive — it is what Amrani et al. searched
in order to over-nominate — while every real on-target site in the design is
`NNGRRT`, and SaCas9 cleaves `NNGRRN` sites without the terminal T substantially less
efficiently. When the two PAM tiers disagree, the canonical one governs.*

**Step 3 — tie-breaks, in order:** (a) P(escape) paired with ICP27g1; (b)
conservation on the gene-level corpus; (c) conservation on the sub-genomic tier;
(d) number of covering records on the sub-genomic tier.

**Reported but never decisive:** GC and local 200 bp GC, seed-region mismatch counts,
21-nt-variant off-target counts. These appear in the table, are discussed, and do not
enter the rule — because a rule with eleven weighted terms is a rule that can be made
to say anything.

---

## 4. Escape probability, computed for `RL2_5335+` and `RL2_3364+`

This number did not exist before this document. It was produced by rebuilding the
stage-7 presence matrix exactly as `src/escape.py` builds it (same nuclease grammar,
same 183 genomes, same `bm.presence_matrix`) and reusing the per-site `q` written by
the canonical stage-7 run. **The reconstruction is validated before any new number
is read off it**: it must reproduce stage 7's own published Amrani lead-pair
P(escape) to within 1e-12, and it does — 0.00097198790268153488 here against
`results/escape_summary.json`'s 0.0009719879026815389. A 2,000,000-draw seeded Monte
Carlo over the same presence matrix agrees with the closed form for both quoted
pairs: 9.975e-04 ± 2.2e-05 against 9.7199e-04 (1.14 sigma), and 2.595e-04 ±
1.1e-05 against 2.5145e-04 (0.71 sigma). `python src/reconcile.py` prints all of
this.

| ICP0 guide, paired with ICP27g1 | joint conservation | both present | P(escape) | fold better than the published pair |
|---|---|---|---|---|
| ICP0g2 (`RL2_4496+`) — as published | 0.8251 | 151/183 | 9.7199e-04 | 1.00 |
| **`RL2_5335+`** | 0.9781 | 179/183 | **2.5145e-04** | **3.87** |
| `RL2_3364+` | 0.9781 | 179/183 | **2.5145e-04** | **3.87** |
| `RL2_5080+` | 0.9781 | 179/183 | 2.5145e-04 | 3.87 |
| `RL2_3441+` | 0.9781 | 179/183 | 2.5145e-04 | 3.87 |

### The escape model cannot choose between these four, and it is important to say so

All four are present in all 183 genomes, so their presence vectors are identical.
All four cut at codons whose surrounding 11-codon window contains **zero** codons
that vary in amino acid across the isolates, so their tolerance vectors are identical
too. Their repair-escape probabilities are therefore the same number:

```
RL2_3364+   0.0094858626500384
RL2_3441+   0.0094858626558352
RL2_5080+   0.0094858626498880
RL2_5335+   0.0094858626500384
```

The spread is 5.9e-12 absolute, 6.3e-10 relative — floating-point summation order
inside `inframe_viability`'s placement average, not biology.

**Stage 7's preference for `RL2_5080+` was an artefact of that.** The line in
`results/escape_model_report.md` labelled "best pair keeping ICP27g1, ICP0 guide
swapped only" is computed with `np.argmin` over those four values, and `RL2_5080+`
happens to hold the smallest float by roughly one part in 1.6 billion. It is not a
finding. It is `argmin` over a tie.

Two consequences:

* The escape axis is **live for the ICP0g2-vs-alternatives comparison** — it separates
  0.847-conserved ICP0g2 (9.72e-04) from any perfectly conserved candidate
  (2.51e-04) by 3.87× — and **dead for the choice among the alternatives**. Step 3(a)
  of the selection rule is invoked and returns nothing.
* The pair's escape probability sits essentially *at* the resolution floor of the
  evidence. Stage 7's own section 7.1 shows that with 183 genomes a two-guide set
  cannot be certified below ~2.6e-04 (Clopper–Pearson 95% upper bound on an absence
  observed 0/183 times is 0.0162). The recommended pair is at 2.5145e-04. **Any
  further improvement in escape probability from guide choice alone would be an
  extrapolation past what this corpus can support** — which is an argument for
  choosing on the off-target axis, where the measurements are not at their floor.

---

## 5. Applying the rule

**Step 1, hard gates.** Of the 14 screened candidates, three survive. This table
is `src/reconcile.py`'s own output, not a hand transcription:

| candidate | G1 cons = 1.000 | G2 corpus ≥ 0.95 | G3 filters | G4 both copies | G5 `NNGRRT` ≤3 mm = 0 | G6 no CDS hits | verdict |
|---|---|---|---|---|---|---|---|
| **RL2_5335+** | ✓ 1.0000 | ✓ 1.0000 | ✓ | ✓ 2 | ✓ 0 | ✓ | **survives** |
| **RL2_3364+** | ✓ 1.0000 | ✓ 0.9788 | ✓ | ✓ 2 | ✓ 0 | ✓ | **survives** |
| **RL2_5080+** | ✓ 1.0000 | ✓ 0.9907 | ✓ | ✓ 2 | ✓ 0 | ✓ | **survives** |
| RL2_3441+ | ✓ 1.0000 | ✓ 1.0000 | ✓ | ✓ 2 | ✗ 6 | ✓ | eliminated |
| RL2_3196- | ✗ 0.9945 | ✓ 0.9857 | ✓ | ✓ 2 | ✗ 1 | ✗ 1 at 20 nt + 0 at 21 nt | eliminated |
| RL2_5167- | ✗ 0.9945 | ✓ 1.0000 | ✓ | ✓ 2 | ✓ 0 | ✗ 1 at 20 nt + 1 at 21 nt | eliminated |
| RL2_5395- | ✗ 0.9891 | ✓ 1.0000 | ✓ | ✓ 2 | ✓ 0 | ✓ | eliminated |
| RL2_5319+ (ICP0g1) | ✗ 0.9836 | ✓ 0.9835 | ✗ GC 0.8, homo 5 | ✓ 2 | ✗ 1 | ✓ | eliminated |
| RL2_3568- | ✗ 0.9672 | ✓ 0.9835 | ✓ | ✓ 2 | ✓ 0 | ✓ | eliminated |
| RL2_3576- | ✗ 0.9508 | ✗ 0.9409 | ✓ | ✓ 2 | ✗ 2 | ✗ 5 at 20 nt + 0 at 21 nt | eliminated |
| RL2_4810- | ✗ 0.9508 | ✗ 0.9481 | ✓ | ✓ 2 | ✓ 0 | ✓ | eliminated |
| RL2_4659+ | ✗ 0.9454 | ✗ 0.9423 | ✓ | ✓ 2 | ✗ 3 | ✓ | eliminated |
| RL2_3515+ | ✗ 0.8907 | ✗ 0.8759 | ✓ | ✓ 2 | ✓ 0 | ✗ 2 at 20 nt + 0 at 21 nt | eliminated |
| RL2_4496+ (ICP0g2) | ✗ 0.8470 | ✗ 0.8975 | ✓ | ✓ 2 | ✓ 0 | ✓ | eliminated |

**Step 2, primary discriminator — `NNGRRT` sites at ≤ 4 mismatches:**

```
RL2_5335+   3      <   RL2_3364+   5      <   RL2_5080+   10
```

`RL2_5335+` wins outright. Steps 3(a)–(d) are not reached. `src/reconcile.py`
applies the same gates and the same discriminator and prints `SELECTED: RL2_5335+`;
the rule lives in that module's `GATES`, `DISCRIMINATOR` and `TIE_BREAKS` constants
so that the prose above and the arithmetic cannot drift apart.

**For the record, had they been reached:** 3(a) P(escape) is an exact tie (section 4);
3(b) gene-level conservation orders `RL2_5335+` 1.000 > `RL2_5080+` 0.991 >
`RL2_3364+` 0.979 — the same winner. The result does not depend on where in the
tie-break chain the decision lands.

### The trade-off that stops this being a clean domination

`RL2_5335+` does **not** dominate the other two survivors. Testing Pareto dominance
mechanically over the axis set {`cons183`, `cons_gene`, `cons_sub`, `P(escape)`,
`NNGRRT ≤3`, `NNGRRT ≤4`, `NNGRRN ≤3`, `NNGRRN ≤4`, `CDS ≤3`}:

**Pareto front = {`RL2_5335+`, `RL2_3364+`, `RL2_5080+`}** — computed by
`src/reconcile.py`, not asserted.

The single measure that blocks domination is `NNGRRN` sites at ≤ 3 mismatches:
`RL2_5335+` has **2**, the other two have **1** each. `RL2_5335+`'s two are
chr16:75,480,285 (`CHST6`, intronic, PAM `CAGGGA`) and chr18:59,279,076 (`CPLX4`,
intronic, PAM `GAGAGC`) — both non-canonical PAMs, both non-coding, both with 1 seed
mismatch. Adding the seed-region counts to the axis set does not change the front;
adding the 21-nt columns does not either.

So the choice among these three **is** a judgement, and it is exactly the judgement
written into step 2: whether one extra `NNGRRN`-only 3-mismatch site outweighs a 2×
to 3× larger burden in the canonical `NNGRRT` tier. This document says it does not,
for the stated biological reason, and section 6 gives the answer under the opposite
view. On the 21-nt data the question does not arise at all: at 21 nt `RL2_5335+` has
0 `NNGRRT` sites and 6 `NNGRRN` sites at ≤ 4 mismatches — fewer than `RL2_3364+`
(2 / 20), fewer than `RL2_5080+` (2 / 18), and the cleanest profile of **every** guide
screened, published guides included.

Against `RL2_3441+` and against ICP0g2 there is no judgement to make: `RL2_5335+`
Pareto-dominates both outright on the axis set above.

---

## 6. What a different, defensible rule would have said

Honesty requires reporting the sensitivity of the answer to the one judgement it
rests on.

* **If `NNGRRN` sites at ≤ 3 mismatches were the primary discriminator instead**, the
  survivors tie at 1 site (`RL2_3364+`, `RL2_5080+`) and `RL2_5335+` is eliminated
  first. The next tie-break decides: on total `NNGRRN` burden it is `RL2_3364+`
  (54 < 61); on gene-level conservation it is `RL2_5080+` (0.991 > 0.979). So the
  permissive-PAM rule does not even yield a stable single answer without a further
  arbitrary choice, which is part of why it is not the rule adopted.
* **If conservation alone were the rule** (stage 6's rule), the answer is a four-way
  tie broken by identifier — which is how `RL2_3441+` got named, and it is refuted.
* **If escape probability alone were the rule** (stage 7's rule), the answer is a
  four-way tie broken by floating-point noise — which is how `RL2_5080+` got named.
* **If local GC were decisive** (the synthesis/amplification feasibility proxy),
  `RL2_3441+` (0.635) would win, then `RL2_5080+` (0.690), `RL2_3364+` (0.710),
  `RL2_5335+` (0.725). Local GC is deliberately *not* in the rule: all four sit far
  below the 0.840 measured around ICP0g2, which is the level at which Amrani et al.
  report the site could not be amplified or sequenced at all. The axis separates the
  candidates from the published guide, not from each other.

The recommendation is robust to every ordering of the tie-breaks and to the seed-count
and 21-nt data. It is *not* robust to inverting step 2. That inversion is the one
place a reasonable reviewer could land somewhere else, and it is named here rather
than buried.

---

## 7. Recommendations withdrawn, and why

### 7.1 `RL2_3441+` (stage 6) — **withdrawn, refuted**

Stage 6's ranking is not retracted: `RL2_3441+` really is perfectly conserved on both
denominators, really does pass every filter, really does sit in both ICP0 repeat
copies, and really does take the pair's joint conservation from 0.825 to 0.978.
Nothing measured in stage 6 was wrong.

What was wrong was treating a conservation rank as a guide selection. Stage 8
screened it against GRCh38:

| | ICP0g2 (theirs) | RL2_3441+ (stage 6's pick) |
|---|---|---|
| `NNGRRT` sites ≤ 4 mm | 5 | **19** |
| `NNGRRT` sites ≤ 3 mm | 0 | **6** |
| `NNGRRN` sites ≤ 4 mm | 56 | 78 |
| `NNGRRN` sites ≤ 3 mm | 5 | 9 |
| seed ≤ 1 mm, `NNGRRT` | 1 | 7 |

Six near-cognate sites under the *canonical* PAM where the published guide has none,
one of them inside ***ABL1*** (chr9:130,830,349, 3 mismatches, 1 in the seed, PAM
`ATGGGT`, within the gene body). Four more are copies of one repeated sequence on
chrY. On off-target burden alone ICP0g2 is the better of the two guides. Proposing a
guide that is *worse than the incumbent on the axis the incumbent was partly selected
on* is a straightforward failure of the recommendation, and it is withdrawn.

The generic lesson, recorded because it is the reusable part: **`RL2_3441+` was never
uniquely indicated even within stage 6.** Four ICP0 sites tie exactly on stage 6's
own figure of merit; `RL2_3441+` was named because it sorts first by identifier. A
recommendation produced by a tie-break should never have been reported as a finding
without saying that it was one. `src/benchmark_sacas9.py` now emits that caveat
inline.

### 7.2 `RL2_5080+` (stage 7) — **withdrawn, superseded**

Stage 7 named `RL2_5080+` as the ICP0 guide in "best pair keeping ICP27g1, ICP0 guide
swapped only". This is superseded, and more than that, it was never a real preference:
as section 4 shows, `RL2_5080+` was selected by `argmin` over four probabilities that
agree to ten significant figures and differ only in floating-point summation order.

`RL2_5080+` is a perfectly respectable guide — it clears every hard gate and sits on
the Pareto front. It is simply not the best of the three on the discriminator the
rule uses: 10 `NNGRRT` sites at ≤ 4 mismatches against `RL2_5335+`'s 3, with 5 of
them carrying ≤ 1 seed mismatch against `RL2_5335+`'s 1. It is also slightly worse on
the second and third conservation denominators (0.991 and 0.972 against 1.000 and
1.000).

If a reviewer rejects step 2 of the selection rule, `RL2_5080+` is the most likely
alternative answer. That is said plainly in section 6.

### 7.3 Not withdrawn: everything else stage 6 claims

The substantive stage-6 findings stand and are unaffected: 30 of the 46 `NNGRRT`
sites in ICP0 are better conserved than ICP0g2; 12 of those pass every filter; the
published `>70%` threshold admits 69 of 81 sites across both genes; their lead pair's
joint conservation (0.825) is below either marginal; and the local GC around ICP0g2
is 0.840, independently corroborating their own report that the site could not be
amplified. Stage 8 sharpened the conclusion rather than overturning it: **at least one
perfectly conserved ICP0 site is also cleaner against the human genome than the guide
they took to the clinic.** That site is `RL2_5335+`.

---

## 8. What this still does not settle

1. **No activity measurement.** Conservation and off-target counts say nothing about
   whether SaCas9 cuts `RL2_5335+` efficiently. Amrani et al. screened six pairwise
   combinations in Vero cells and chose their lead pair on measured antiviral
   activity. That is evidence of a kind this repository contains none of, and it is
   the single most likely reason a well-founded paper recommendation could still be
   the wrong drug.
2. **Bulges are not modelled.** Cas-OFFinder, which Amrani et al. used, searches
   DNA/RNA bulges; `src/offtarget.py` searches substitutions only. A bulge-tolerant
   site would be invisible here.
3. **No cell-based off-target validation.** No GUIDE-seq, no CIRCLE-seq. Nominated
   sites are not measured cleavage.
4. **The escape model's absolute scale is soft.** Its NHEJ indel-length spectrum is a
   labelled ASSUMPTION, absolute escape probabilities move over ~3 orders of
   magnitude across the plausible parameter space, and the pair recommended here sits
   at the corpus's resolution floor. What is stable across the sweep is the *ranking*,
   which is all this document uses it for — and even that ranking is silent between
   the three finalists.
5. **183 genomes is the denominator.** Effective sample size is nearer 132 at 99.9%
   identity; a site absent in 0/183 isolates is still consistent with a population
   absence frequency up to 1.6%.
6. **AAV packaging, sgRNA scaffold compatibility and synthesis feasibility** are
   outside everything measured here.

The claim this document supports, stated at the width the evidence actually covers:
**within Amrani et al.'s own SaCas9 site space, `RL2_5335+` is better than their
ICP0g2 on cross-isolate conservation, on multiplex escape probability, and on human
off-target burden simultaneously — and it is the best of the three candidates that
survive a stated set of hard gates, under one stated priority between two PAM tiers.**
It is not a claim that it is a better drug.

---

## 9. Reproduce

```bash
python run_pipeline.py --skip-fetch --benchmark-sacas9   # stage 6 pool + pairs
python run_pipeline.py --skip-fetch --escape             # stage 7 q, presence, k-curve
python src/offtarget.py --stage fetch                    # once: ~4.1 GB on disk
python run_pipeline.py --skip-fetch --offtarget          # stage 8 GRCh38 screen
python src/reconcile.py                                  # this document's arithmetic
```

`src/reconcile.py` refuses to run if any of the three stages' outputs is missing, and
**refuses to report anything if its rebuilt escape model does not reproduce stage 7's
published Amrani lead-pair P(escape) to within 1e-12** — because a reconstruction that
silently disagreed with the stage that produced it would leave every paired
`P(escape)` here unsupported. Observed on this machine: `0.00097198790268153488`
against `results/escape_summary.json`'s `0.0009719879026815389`.

Inputs to every number in this document:
`results/sacas9_benchmark_pool.tsv`, `results/sacas9_benchmark_pairs.tsv`,
`results/escape_site_parameters.tsv`, `results/escape_guide_sets.tsv`,
`results/offtarget_summary.tsv`, `results/offtarget_annotated.tsv`.

The paired `P(escape)` and joint-conservation columns are obtained by loading
`src/escape.py`'s `EscapeModel` over the stage-7 presence matrix and the `q_site`
column of `escape_site_parameters.tsv`; the reconstruction is checked by reproducing
the published Amrani lead-pair value to machine precision before any new number is
read off it.
