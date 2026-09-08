# Conserved CRISPR-Cas9 target sites in HSV-1

[![DOI](https://zenodo.org/badge/1361724608.svg)](https://doi.org/10.5281/zenodo.22664837)

**Analysis code and archived results for the preprint *"What a genome corpus can and
cannot certify about CRISPR antiviral escape: resolution floors, joint coverage, and a
worked audit of an HSV-1 guide pair"*.**
Manuscript: [`manuscript/manuscript.md`](manuscript/manuscript.md) ·
References: [`manuscript/references.md`](manuscript/references.md) ·
Reconciled recommendation: [`results/recommendation.md`](results/recommendation.md)

Multiplex CRISPR antivirals are justified by an escape argument — cut a viral genome at
several conserved sites at once and no single repair event restores a viable, uncuttable
genome — and that argument is almost always supported by a table of per-guide
conservation percentages computed against a public sequence corpus. This repository
implements the measurement end-to-end for herpes simplex virus type 1 and shows that the
evidentiary chain has two structural weaknesses that are properties of the method rather
than of any particular design. **First, the corpus imposes a resolution floor**: a site
absent in 0 of 183 complete HSV-1 genomes is still consistent with a population absence
frequency of 1.62 %, so sampling resolution — not the repair model — sets how low an
escape probability can be *demonstrated*. **Second, joint intactness is not the product
of marginal conservations and cannot be recovered from a per-guide table at all**; it
needs a per-genome presence matrix. As a worked case study the pipeline audits the
clinical-stage EBT-104 guide pair of Amrani et al. (2024) on real data, finds a
recurrent single-base variant that removes one of its two target sites from 25 of 183
sequenced isolates, and shows that changing one guide improves every axis that can be
computed here. The audit is a demonstration of the two methodological points on
clinical-stage guides, not an indictment: the published workflow could not have seen the
gap, and the repository withdraws two of its own recommendations for the same reason.

Everything is pure Python plus Biopython, pandas, numpy and matplotlib. No aligner, no
conda, no WSL, no compiled tooling beyond pip wheels. Conservation is defined as exact
substring presence rather than alignment identity, so there is no gap penalty, no
substitution matrix, no heuristic and no random seed anywhere in the measurement.

---

## What it found

All numbers below are from the pinned run of **2026-09-07**, against **183 complete
HSV-1 genomes**. Every one of them is asserted by [`verify.py`](verify.py) against the
version-controlled result files, and the section each belongs to in the manuscript is
given so you can find the argument behind it.

| | value | where |
|---|---|---|
| Complete HSV-1 genomes in GenBank at the retrieval date | **183** (147,898–159,092 bp) | Methods 5.1 |
| SpCas9 `NGG` candidate sites across 7 essential genes | **4,777** | Results 2.1 |
| ... present in **all** 183 genomes | **833** (17.4 %) | Results 2.1 |
| ... and also passing poly-T, homopolymer and GC filters | **644** | Results 2.1 |
| SaCas9 `NNGRRT` candidate sites, same 7 genes | **448** (67 perfect, 15.0 %) | Figure 3 |
| Clopper–Pearson 95 % upper limit on 0 absences in 183 | **0.0162** | Results 2.1 |
| ... so one site cannot be certified below | **~1.6 × 10⁻²** | Results 2.1 |
| ... and a pair cannot be certified below | **~2.6 × 10⁻⁴** | Results 2.1 |
| Certifiable minimum guide count at a 10⁻⁶ threshold | **4**, not the 3 a point estimate implies | Table 1 |
| EBT-104 lead pair marginal conservations | 0.847 (ICP0g2) and 0.978 (ICP27g1) | Table 2 |
| ... independence would predict joint intactness | 0.8285 | Results 2.2 |
| ... **measured** joint intactness | **0.8251** (151/183) — *below either marginal* | Results 2.2 |
| Isolates lacking the ICP0g2 target site | **28**; 25 carry an identical G→A at spacer position 9, 0 are N-gaps, 3 unresolved | Results 2.3, Figure 5 |
| ... collapsing to independent lineages at 99.9 % identity | **13 clusters**, sampled 1967–2020, ≥ 2 continents | Results 2.3 |
| ICP0 `NNGRRT` sites better conserved than ICP0g2 | **30 of 46**; 12 also pass every filter; 4 are perfect | Results 2.3 |
| P(escape) of the published pair, per exposed genome | **9.720 × 10⁻⁴**; **96.2 %** of it from the 32 isolates that already lost a site | Table 3 |
| Same architecture, one guide swapped (RL2_5335+ + ICP27g1) | joint **0.978**, P(escape) **2.515 × 10⁻⁴**, fewer off-targets on every measure | Table 4 |
| Guide-set rankings stable across the sensitivity sweep | **27 of 37** settings, spanning 3.1 orders of magnitude in absolute P(escape) | Results 2.7 |

Two of this repository's own recommendations were withdrawn during the work, one of them
refuted by its own subsequent off-target screening. Both withdrawals are documented in
[`results/recommendation.md`](results/recommendation.md) rather than quietly edited out.

---

## Reproduce

### 0. Requirements

| | |
|---|---|
| Python | **3.14** (developed and pinned on CPython 3.14.5, 64-bit) |
| OS | developed on Windows 11; the code is pure Python and platform-independent |
| Disk, stages 1–7 | **~250 MB** (183 genome FASTAs ~28 MB, gene-level corpus ~80 MB, results ~8 MB, plus caches) |
| Disk, stage 8 | **+ ~4.1 GB** (see the warning below) |
| Network | needed only by stage 1 and by the gene-corpus half of stage 5. Everything else runs offline from the caches those two populate. |

```bash
git clone https://github.com/ivanjalid1/hsv1-crispr-escape.git
cd hsv1-crispr-escape

python -m venv .venv
# Windows
.venv\Scripts\python.exe -m pip install --only-binary=:all: -r requirements.txt
# macOS / Linux
.venv/bin/python -m pip install --only-binary=:all: -r requirements.txt
```

`--only-binary=:all:` is recommended: it makes pip fail loudly if a cp314 wheel is ever
missing, rather than silently attempting a source build that would need a C toolchain.
All four scientific packages ship prebuilt `cp314` wheels at the pinned versions
(`numpy 2.5.3`, `pandas 3.0.5`, `biopython 1.88`, `matplotlib 3.11.1`). Note that pandas
3.0 is a major release with behavioural changes from 2.x; this code is written against
3.0 and was run against it.

### 1. Set the required environment variables

```powershell
# PowerShell
$env:NCBI_EMAIL = "you@example.org"      # REQUIRED
$env:NCBI_API_KEY = "<your key>"         # optional
```

```bash
# bash / zsh
export NCBI_EMAIL="you@example.org"      # REQUIRED
export NCBI_API_KEY="<your key>"         # optional
```

| variable | required | why |
|---|---|---|
| `NCBI_EMAIL` | **yes**, for any stage that contacts NCBI | NCBI Entrez requires a contact address so it can reach you if a script misbehaves; it is [their stated condition of use](https://www.ncbi.nlm.nih.gov/books/NBK25497/). If it is unset the pipeline exits immediately with instructions and does **not** contact NCBI. |
| `NCBI_API_KEY` | no | Raises the NCBI rate limit from 3 to 10 requests/second; the pipeline adjusts its own throttle accordingly (0.40 s → 0.12 s between requests). Get one at <https://www.ncbi.nlm.nih.gov/account/settings/>. Stage 1 takes roughly a minute either way. |

**No email address is hardcoded anywhere in this repository**, and none is written into
any output file. Both variables are read from the environment only. `.env` files are
gitignored.

### 2. Run the pipeline

Commands are given for Windows; substitute `.venv/bin/python` elsewhere. Runtimes are
wall-clock on a 2024-era laptop and are dominated by I/O, not CPU, except stage 8.

```powershell
# --- fast smoke test: 5 genomes end to end, proves the install works -----------
.venv\Scripts\python.exe run_pipeline.py --limit 5                 # ~10 s, ~1 MB

# --- THE MAIN RUN: stages 1-4, SpCas9, all 183 genomes -------------------------
.venv\Scripts\python.exe run_pipeline.py                           # ~80 s cold, ~5 s warm
#   downloads ~28 MB, writes ~3 MB of tables
#   -> results/summary.json, results/guides_ranked.tsv  (the 4,777 / 833 / 644 numbers)

# --- stage 5: denominator robustness -------------------------------------------
.venv\Scripts\python.exe run_pipeline.py --robustness              # ~6 min, downloads ~80 MB once
.venv\Scripts\python.exe run_pipeline.py --robustness --skip-gene-corpus   # parts A+B only, offline, ~90 s

# --- stage 6: SaCas9 head-to-head vs Amrani et al. 2024 ------------------------
.venv\Scripts\python.exe run_pipeline.py --skip-fetch --benchmark-sacas9   # ~50 s, offline

# --- stage 7: multiplex escape-probability model -------------------------------
.venv\Scripts\python.exe run_pipeline.py --skip-fetch --escape     # ~3 min, offline

# --- reconcile stages 6-8 into one recommendation ------------------------------
.venv\Scripts\python.exe src\reconcile.py                          # ~2 s, offline

# --- figures 1-5 ---------------------------------------------------------------
.venv\Scripts\python.exe src\figures.py                            # ~30 s, offline
#   -> figures/fig1..fig5 .png and .pdf
```

> ### ⚠ Stage 8 needs a 4 GB human genome download. Read this before starting it.
>
> The human off-target screen requires the Ensembl release-116 GRCh38 primary assembly
> and its GTF. That is **~1.0 GB downloaded** and **~4.1 GB occupied on disk** once
> decompressed into the scanning index, and the scan itself takes **~6 minutes** (319 s
> measured for 50 guides × 2 strands × 3.10 Gb). None of it is needed for any
> conservation, joint-coverage or escape result in the paper.
>
> The download **never happens implicitly**. If the cache is absent, stage 8 logs why it
> is skipping and the pipeline continues. To build it you must ask, explicitly and once:
>
> ```powershell
> .venv\Scripts\python.exe src\offtarget.py --stage fetch    # ~1.0 GB down, ~4.1 GB on disk
> .venv\Scripts\python.exe src\offtarget.py                  # ~6 min
> ```
>
> The download is verified three independent ways — Ensembl's published BSD `sum`
> checksum, the gzip CRC-32/ISIZE trailer, and a locally computed SHA-256 recorded in
> `data/genome/genome_manifest.json` — and the release is **pinned, not "current"**, so a
> floating release cannot silently change the numbers.

### 3. Run the tests

**123 tests across eight files, all offline, no network and no genome download.**

```powershell
.venv\Scripts\python.exe tests\test_core.py        # 7   conservation scanner
.venv\Scripts\python.exe tests\test_nuclease.py    # 11  PAM model; scanner == brute force, every PAM
.venv\Scripts\python.exe tests\test_robustness.py  # 23  stage 5
.venv\Scripts\python.exe tests\test_benchmark.py   # 6   stage 6
.venv\Scripts\python.exe tests\test_escape.py      # 31  stage 7
.venv\Scripts\python.exe tests\test_offtarget.py   # 17  stage 8 (~135 s; no genome needed)
.venv\Scripts\python.exe tests\test_reconcile.py   # 8   the selection rule, pinned
.venv\Scripts\python.exe tests\test_figures.py     # 20  figures render byte-identically
```

The two performance-critical scanners are not merely tested, they are *proved*:
`tests/test_nuclease.py` asserts set-identity between the jump scanner and a naive
both-strand substring search for every supported PAM on randomised sequences, GC-rich
sequences and the real HSV-1 reference, and `tests/test_offtarget.py` asserts set
identity between the pigeonhole scanner and a structurally independent brute-force scan
over the whole of chromosome 21. The optimisations are not approximations.

---

## Verify the headline numbers

[`verify.py`](verify.py) checks every number in the tables above against the
version-controlled result files, and recomputes the Clopper–Pearson resolution floor
from first principles. **The fast path needs no network, no downloads and under a
second.**

```powershell
# FAST PATH -- offline, < 1 s, nothing to download. Start here.
.venv\Scripts\python.exe verify.py
#   82 assertions against the pinned, version-controlled results.
#   Answers: does the archived record actually say what the paper says?
```

```powershell
# SLOWER -- recomputes stages 2-4 and checks byte-identity, ~30 s, offline
# but needs data/raw/ populated (run `run_pipeline.py` once first).
.venv\Scripts\python.exe verify.py --rerun

# SLOWEST -- also recomputes stages 6 and 7, ~5 min, offline.
.venv\Scripts\python.exe verify.py --rerun --full
```

`--rerun` re-runs the pipeline with `--force` and compares every regenerated file
against the SHA-256 recorded in [`results/CHECKSUMS.sha256`](results/CHECKSUMS.sha256)
at pin time. That file also pins the large outputs that are *not* version-controlled, so
byte-identity can be checked for them without shipping 4 MB of TSV.

**Stage 8 is deliberately outside every `verify.py` path**, because it costs 4 GB and
six minutes and no conservation, joint-coverage or escape claim depends on it. Its
pinned outputs are still checked on the fast path, from the version-controlled
`results/offtarget_summary.tsv`. To recompute it yourself, see the warning box above.

---

## Data provenance, and why the retrieval date is the pin

**This analysis is pinned to a retrieval date, not to a query.** That distinction is the
single most important thing to understand before re-running anything.

The complete-genome corpus was retrieved from NCBI Nucleotide on
**2026-09-07T21:56:25Z** with the exact query:

```
txid10298[Organism:exp] AND 145000:160000[SLEN] AND biomol_genomic[PROP]
  AND "complete genome"[Title] NOT patent[PROP]
```

NCBI reported and returned **183** records (147,898–159,092 bp). A second, gene-level
corpus of 4,945 further HSV-1 records was built with two further queries, given in
Methods 5.1 and written verbatim into the stage-5 report.

**Re-running that query today will return MORE than 183 genomes**, because GenBank
grows. Every number in the paper and in the tables above is a function of the corpus
that existed at that timestamp. A reviewer who re-runs the pipeline in a year and gets
different counts has not found an irreproducibility; they have measured a different,
larger corpus. That is a property of the field's evidence base — indeed it is the point
of the resolution-floor argument, which only improves as the corpus grows — and not a
defect of the code.

The pinned record is therefore two files, both version-controlled:

| file | what it pins |
|---|---|
| [`data/manifest.tsv`](data/manifest.tsv) | one row per genome: accession, length, ambiguous-base count, strain, isolate, country, collection date, host, isolation source, definition, completeness label, **SHA-256 of the downloaded FASTA**, plus the exact Entrez query and the UTC retrieval timestamp in every row |
| [`results/run_log.json`](results/run_log.json) | the full **sorted accession list actually used**, the exact query, the retrieval timestamp, the command line, Python and package versions, which stages executed, and the summary statistics |

To reproduce the paper's numbers exactly rather than to measure today's GenBank, fetch
the accessions listed in `results/run_log.json` and verify each FASTA against the
SHA-256 in `data/manifest.tsv`. To measure today's GenBank instead, just run
`run_pipeline.py` and expect the counts to move.

Two further pins, for the same reason:
[`data/gene_corpus_manifest.tsv`](data/gene_corpus_manifest.tsv) records the 4,945
gene-level records per-record, and the human assembly is pinned to **Ensembl release
116** with a recorded SHA-256 (`d8c3af00…c187b92`, 881,964,081 bytes) rather than to
"current".

Nothing is ever padded or imputed: if NCBI returns fewer records than its own count
reports, the shortfall is logged as a warning and the manifest contains exactly what was
retrieved.

---

## What is in `results/`, and what is not

`results/` is regenerable, so the default is to leave it out of version control. The
exceptions are tracked deliberately, and the rule is: **a reviewer must be able to check
any number in the paper without running anything, and must be able to regenerate
everything.**

**Tracked** (~340 KB total) — small, decision-bearing, or provenance:

* `recommendation.md` and `recommendation_table.tsv` — the reconciliation of stages 6, 7
  and 8, which proposed three different guides on three different axes. The table is
  computed; the selection rule and the withdrawals are an argument, and an argument that
  lives only in an ignored directory is an argument nobody can review.
* `summary.json`, `run_log.json`, `CHECKSUMS.sha256` — the provenance record of the
  pinned run, and the SHA-256 of every regenerable file it produced.
* The four stage reports — `robustness_report.md`, `sacas9_benchmark_report.md`,
  `escape_model_report.md`, `offtarget_report.md` — each of which carries claims the
  paper cites.
* The four machine-readable stage summaries (`*_summary.json`) — every headline number
  in a form `verify.py` can assert against.
* Small decision-bearing tables, all under 30 KB: `escape_k_curve*.tsv`,
  `escape_guide_sets.tsv`, `escape_sensitivity*.tsv`, `escape_tolerance_qc.tsv`,
  `offtarget_summary.tsv`, `offtarget_annotated.tsv` (the coding-exon context that
  withdrew one of our own recommendations), `sacas9_benchmark_pool.tsv`,
  `robustness_redundancy.tsv` and the four tiny robustness tables.

**Not tracked, regenerate them** — bulk per-guide and per-site enumerations totalling
~4 MB, plus the cached FASTAs and the 4 GB human genome index:

`guides_candidates.tsv`, `conservation.tsv`, `guides_ranked.tsv`,
`escape_site_parameters.tsv`, `escape_codon_tolerance.tsv`, `offtarget_sites.tsv`,
`robustness_conservation_modes.tsv`, `robustness_gene_level.tsv`,
`sacas9_benchmark_pairs.tsv`, the alternative-nuclease re-runs under `results/sacas9*/`,
and everything under `data/raw/` and `data/genome/`.

Every one of those is pinned by SHA-256 in `results/CHECKSUMS.sha256`, so
`python verify.py --rerun` can confirm that a regenerated copy is byte-identical to the
one the paper was written from — without the repository having to carry it.

One file is *deliberately* left untracked even though it is small:
`data/fetch_summary.json` describes whichever fetch ran last, including a smoke test
against a different taxon, and a stale copy of it in version control would be actively
misleading. The provenance records are `data/manifest.tsv` and `results/run_log.json`.

Re-running the pipeline will overwrite the tracked files in `results/`. That is intended:
`git diff` then shows you exactly what moved relative to the pinned run.

---

## Pipeline stages

Stages 1–4 always run. Stages 5–8 are opt-in, and each skips with an explanatory message
if its inputs are absent rather than failing.

| Stage | Module | Flag | Cost | Network | Output |
|---|---|---|---|---|---|
| 1 | `src/fetch_genomes.py` | (always) | ~60 s, ~28 MB | **yes** | `data/raw/*.fasta`, `data/manifest.tsv` |
| 2 | `src/extract_guides.py` | (always) | ~5 s | no | `results/guides_candidates.tsv` |
| 3 | `src/conservation.py` | (always) | ~3 s | no | `results/conservation.tsv` |
| 4 | `src/report.py` | (always) | ~2 s | no | `results/guides_ranked.tsv`, `results/summary.json` |
| 5 | `src/robustness.py` | `--robustness` | ~6 min, ~80 MB | **yes** (gene corpus only; `--skip-gene-corpus` for parts A+B offline) | `results/robustness_report.md` + tables |
| 6 | `src/benchmark_sacas9.py` | `--benchmark-sacas9` | ~50 s | no | `results/sacas9_benchmark_report.md` + tables |
| 7 | `src/escape.py` | `--escape` | ~3 min | no | `results/escape_model_report.md` + tables |
| 8 | `src/offtarget.py` | `--offtarget` | **~6 min, ~4.1 GB disk** | **yes, once, explicitly** | `results/offtarget_report.md` + tables |
| — | `src/reconcile.py` | (standalone) | ~2 s | no | `results/recommendation_table.tsv` |
| — | `src/figures.py` | (standalone) | ~30 s | no | `figures/fig1–fig5 .png/.pdf` |
| — | `run_pipeline.py` | — | — | — | `results/run_log.json` |

Outputs are **namespaced by nuclease**. SpCas9 keeps the historical paths
(`results/guides_ranked.tsv`); any other nuclease writes to `results/<tag>/`
(`results/sacas9/`, `results/sacas9-nngrrn/`, `results/sacas9-20nt/`), so runs under
different nucleases never overwrite each other. Stages 6, 7 and 8 are nuclease-fixed by
construction — they exist to compare against Amrani et al.'s SaCas9 guides in that
paper's own grammar — and always write to `results/`.

Useful flags: `--force` recomputes every stage (needed if you change a parameter the
staleness check cannot see), `--skip-fetch` reuses an existing manifest without
contacting NCBI, `--refresh` bypasses the download cache, `--limit N` takes the first N
accessions after sorting (so it is deterministic), `-v` enables debug logging. Every
stage is also runnable on its own: `python src/fetch_genomes.py --help`, and so on.
Per-stage narrative documentation is under
[Pipeline stages in detail](#pipeline-stages-in-detail).

Other supported configurations:

```powershell
.venv\Scripts\python.exe run_pipeline.py --include-partial                      # admit near-full-length "partial genome" isolates
.venv\Scripts\python.exe run_pipeline.py --max-ambiguous-fraction 0.001         # drop poorly resolved assemblies from the denominator
.venv\Scripts\python.exe run_pipeline.py --taxid 10310 --reference NC_001798    # HSV-2 instead of HSV-1
.venv\Scripts\python.exe run_pipeline.py --genes UL30,UL29,UL54                 # a different gene set
.venv\Scripts\python.exe run_pipeline.py --nuclease sacas9 --skip-fetch         # SaCas9, 21 nt + NNGRRT
.venv\Scripts\python.exe run_pipeline.py --nuclease sacas9 --pam NNGRRN --skip-fetch
.venv\Scripts\python.exe run_pipeline.py --nuclease sacas9 --spacer-length 20 --skip-fetch
```

---

## How to cite

Please cite **both** the preprint and the archived software release.

* **Preprint** — Ivan Heredia Jalid. *What a genome corpus can and cannot certify about CRISPR
  antiviral escape: resolution floors, joint coverage, and a worked audit of an HSV-1
  guide pair.* bioRxiv (2026). doi:`[BIORXIV-DOI]`
* **Software** — Ivan Heredia Jalid. *Conserved CRISPR-Cas9 target sites in HSV-1* (version
  1.0.0). Zenodo (2026). doi:`10.5281/zenodo.22664837` —
  <https://doi.org/10.5281/zenodo.22664837>. This is the **concept DOI**: it always
  resolves to the most recent archived version. Zenodo also mints a version-specific
  DOI for each release; cite the concept DOI unless you need to pin a single version.

GitHub renders [`CITATION.cff`](CITATION.cff) as a "Cite this repository" button that
emits both, and Zenodo reads [`CITATION.cff`](CITATION.cff) and
[`.zenodo.json`](.zenodo.json) when minting the DOI from a release tag. The Zenodo DOI
is now filled in everywhere; the only placeholder left is the bioRxiv DOI above, which
is filled once the preprint is posted. The checklist is [`PUBLISH.md`](PUBLISH.md).

If you use only the escape model or only the off-target scanner, please still cite the
preprint — the arguments they implement are what make the numbers interpretable.

---

## Licence

**MIT** — see [`LICENSE`](LICENSE). Permissive, OSI-approved, and compatible with every
dependency (numpy, pandas, matplotlib, contourpy, cycler, kiwisolver and pyparsing are
BSD-family; biopython uses the BSD-style Biopython License; six, fonttools and pyparsing
are MIT; pillow is MIT-CMU/HPND; python-dateutil, tzdata and packaging are Apache-2.0 or
dual-licensed). No dependency carries a copyleft obligation, so nothing forces a
stronger licence.

`refs/` contains third-party material redistributed under **its own** licences — the
verbatim CC BY-NC-ND 4.0 full text of Amrani et al. (2024), and a CC BY 4.0 derived
table from the source data of Ramadoss et al. (2025). The MIT grant does not extend to
them; the full notice is in [`NOTICE.md`](NOTICE.md). It is a separate file because
GitHub's licence detector only reports **MIT** when `LICENSE` is the unmodified MIT
text, so the third-party notice cannot live inside it. Sequence data is retrieved at
run time from NCBI and Ensembl and is not redistributed here.

---

## Repository layout

```
run_pipeline.py            single entrypoint: stages 1-4 plus opt-in 5-8, and the run log
verify.py                  check the headline numbers reproduce (fast path, offline, < 1 s)
requirements.txt           pinned, installed and tested on CPython 3.14.5
LICENSE                    MIT, unmodified so GitHub detects it
NOTICE.md                  third-party licences covering refs/ (CC BY-NC-ND 4.0, CC BY 4.0)
CITATION.cff               "Cite this repository"; read by GitHub and by Zenodo
.zenodo.json               deposition metadata for the archived DOI
PUBLISH.md                 the author's checklist for the GitHub push and the Zenodo mint

src/                       the pipeline (see "Repository layout" detail below)
tests/                     123 offline tests across eight files
manuscript/                manuscript.md, references.md
figures/                   fig1-fig5, .png and .pdf
refs/                      third-party source material, with extraction notes
results/                   pinned outputs (tracked by exception; see the policy above)
data/                      manifests (tracked); raw downloads and the GRCh38 cache (not)
```

---

# Reference documentation

*Everything below is the working documentation of the analysis: how each measurement is
defined, what was actually observed at each stage, and where the method fails. It was
written as the work proceeded and is kept as the detailed record behind the summary
above.*

## Method

**Conservation is defined as exact substring presence, not alignment identity.**

For every candidate guide we take the target site `protospacer + PAM` -- a 23-mer for
SpCas9 (`20 + NGG`), a 27-mer for SaCas9 (`21 + NNGRRT`) -- and ask, for each
downloaded genome, a single yes/no question: *does this exact site occur verbatim in
that genome, on either strand?*

```
conservation_fraction = n_strains_present / n_strains_total
```

This is the right question for the biology. Cas9 cleavage requires an intact PAM and
near-perfect protospacer complementarity, so a site that differs by even one base is,
for guide-design purposes, a different site. Framing conservation this way:

* removes every external binary dependency and every alignment parameter (gap
  penalties, substitution matrices) from the analysis;
* is exactly reproducible — there is no heuristic and no random seed;
* is conservative. Exact matching can only ever *under*-state conservation, never
  over-state it, so a guide reported at 100% really is present in all genomes.

The trade-off, stated plainly: this method cannot distinguish "the site is genuinely
absent/variant in that isolate" from "that region of that assembly is unresolved (N)
or truncated". Stage 5 (`--robustness`) removes that limitation by locating the
homologous region with k-mer anchors and reporting the two cases separately; the
measured size of the effect is under
[Robustness](#robustness-what-the-denominator-is-worth) and in
[Limitations](#limitations).

### Scanning algorithm

The naive approach — search each of G guides in each of S genomes of length L — is
`O(G x S x L)`; with ~4,800 guides and 183 genomes of 152 kb that is far too slow.
Storing every 23-mer of every genome in a set is `O(S x L)` entries — tens of
millions of strings, gigabytes of RAM.

Instead `src/conservation.py` inverts the problem and scans each genome once:

1. Build two small dictionaries keyed on the guide 23-mers themselves —
   `fwd[t23] -> guide_id` and `rev[revcomp(t23)] -> guide_id`. Together these cover
   both strands for every guide. Memory is `O(G)`: a few thousand entries.
2. Every target site carries the *literal* positions of its PAM pattern at fixed
   offsets: an SpCas9 site always has `GG` at offset 21, a SaCas9 `NNGRRT` site
   always has `G` at offset 23 and `T` at offset 26. Rather than testing all L
   positions we jump between occurrences of one such literal anchor using
   `str.find`, which runs in optimised C. Only those candidate positions produce a
   slice and a dict lookup. The anchors are derived from the PAM pattern in
   `src/nuclease.py`, so no PAM literal is hardcoded in the scanner; where several
   anchors are equally long the scanner picks, per sequence, the one that occurs
   least often (in a 68% GC genome the `T` of `NNGRRT` is ~3x rarer than the `G`).
   A PAM with no literal position at all falls back to an exhaustive positional
   scan -- slower, still exact.

Work is `O(S x L)` with a small constant, memory is `O(G)`. Measured: **~95x faster
than brute force**, scoring 4,777 guides against 183 genomes in about 3 seconds.
`tests/test_core.py::test_scan_matches_bruteforce` proves on randomised sequences
that the jump scan returns *exactly* the same result set as a naive both-strand
search — the optimisation is not an approximation. Because the anchors are now
PAM-derived rather than hardcoded, that proof is repeated for **every supported PAM**
in `tests/test_nuclease.py`: `NGG`, `NAG`, `NNGRRT`, `NNGRRN`, a literal-free `NNN`
(which exercises the exhaustive fallback) and two synthetic grammars, each against a
naive both-strand substring search on randomised sequences, on GC-rich sequences and
on the real HSV-1 reference. `scan_genome_naive` -- that reference implementation --
ships in `src/conservation.py` rather than only in the tests.

### The nuclease model (`src/nuclease.py`)

A `Nuclease` is `spacer_length` + an IUPAC `pam` pattern 3' of the protospacer, and
everything else is derived from it: footprint length, PAM matching on both strands
(the minus-strand PAM is matched as the reverse-complement *pattern* at the low end
of the plus-strand window), the scanner's literal anchors, the predicted cut site,
and the output namespace tag.

| nuclease | spacer | PAM | footprint | selector |
|---|---|---|---|---|
| SpCas9 | 20 nt | `NGG` | 23 nt | `--nuclease spcas9` (default) |
| SaCas9 | 21 nt | `NNGRRT` | 27 nt | `--nuclease sacas9` |
| SaCas9, permissive | 21 nt | `NNGRRN` | 27 nt | `--nuclease sacas9 --pam NNGRRN` |
| SaCas9, Amrani et al. spacer length | 20 nt | `NNGRRT` | 26 nt | `--nuclease sacas9 --spacer-length 20` |

`--pam` accepts any IUPAC string. 5'-PAM nucleases (Cas12a) are **rejected with an
explicit error** rather than half-supported: their cut geometry is staggered and
PAM-distal, and faking it would produce plausible-looking wrong coordinates.

Defaults are unchanged, and this was verified rather than assumed: after the
refactor, stage 2 re-run under SpCas9 produces the same 4,777 guides with every
original column byte-identical, and stages 3-4 reproduce `conservation.tsv`,
`guides_ranked.tsv` and `summary.json` byte-for-byte.

---

## Pipeline stages in detail

The summary table, with runtimes and disk cost, is under
[Pipeline stages](#pipeline-stages) above. This section is the per-stage narrative:
what each one actually does and which decisions inside it are load-bearing.

**A note on stage numbering.** Stages 7 and 8 were developed concurrently and both
were called "stage 7" in their first drafts. This repository fixes **escape = stage
7, off-target = stage 8**, following the order in which they were committed. Every
module docstring, report header, CLI argument group and test file now uses that
numbering; if you find a "stage 7" attached to off-target screening anywhere, it is
stale.

**Stages 7 and 8 are opt-in and are skipped when their inputs are absent.** Stage 7
needs the stage-1 manifest and the cached reference GenBank record; stage 8 needs the
stage-6 candidate pool and the GRCh38 cache. If a prerequisite is missing the stage
logs *why* it is skipping and the pipeline continues. The GRCh38 cache in particular
is **never built implicitly**: it is a ~1.0 GB download occupying ~4.1 GB on disk, so
it must be requested with `python src/offtarget.py --stage fetch` or
`--offtarget-download`. Stages 1-4 are unaffected by either flag and reproduce
byte-identically with or without them.

**1. Fetch.** Queries NCBI Nucleotide, sorts accessions before applying `--limit` (so
`--limit 5` is deterministic), downloads FASTA per accession with caching, batching,
rate limiting and exponential-backoff retry, and builds a manifest. Metadata (strain,
country, collection date, host, isolation source) is pulled cheaply by requesting a
1 bp GenBank slice — that returns the full header and `source` feature for ~2 kB
instead of the ~400 kB of a full annotated flat file.

The default query, built from `--taxid`, `--min-length`, `--max-length`:

```
txid10298[Organism:exp] AND 145000:160000[SLEN] AND biomol_genomic[PROP]
  AND "complete genome"[Title] NOT patent[PROP]
```

It is written verbatim into `data/manifest.tsv` (column `entrez_query`) along with
the UTC retrieval timestamp, and repeated in `results/run_log.json` together with the
full sorted list of accessions used and the exact package versions. `--query`
overrides it entirely.

**2. Extract.** Reads the annotated GenBank record for the reference (default
`NC_001806`, HSV-1 strain 17) and enumerates every site of the selected nuclease on
both strands within the target genes. **Gene coordinates are parsed from the annotation — no coordinates
are hardcoded.** Default targets, all essential or high-value antiviral targets:

| Gene | Product |
|------|---------|
| UL30 | DNA polymerase catalytic subunit |
| UL19 | major capsid protein VP5 |
| UL5  | helicase-primase helicase subunit |
| UL52 | helicase-primase primase subunit |
| UL29 | ICP8 single-stranded DNA binding protein |
| RL2  | ICP0 ubiquitin E3 ligase |
| UL54 | ICP27 multifunctional expression regulator |

Two annotation details are handled explicitly. RL2/ICP0 is **spliced**, so each exon
is treated as a separate contiguous genomic segment — every reported 23-mer is a
contiguous stretch of genomic DNA (what Cas9 actually cuts) that also lies inside
coding sequence. RL2 is also **present twice**, in the TRL and IRL repeats; the two
copies yield identical 23-mers, which are collapsed to one guide with every reference
position kept in `ref_all_positions` and counted in `n_reference_copies`.

**3. Conserve.** As described under [Method](#method).

**4. Report.** Joins, flags, ranks. Ranking is deterministic — identical inputs give
byte-identical output.

**5. Robustness (opt-in, `--robustness`).** Asks what the stage-3 denominator is
actually worth. Three independent attacks on it, all deterministic:

* **Redundancy.** A mash-style bottom-*s* MinHash sketch of canonical *k*-mers per
  genome (pure Python + numpy; splitmix64 mixing, *not* Python's salted `hash()`, so
  sketches are reproducible across processes), all-pairs Jaccard, an ANI-like
  identity, single-linkage clustering, and conservation recomputed on one
  representative per cluster. The cluster count is an effective sample size.
* **Ambiguity.** A sweep over `--max-ambiguous-fraction`, plus an *N-tolerant* match
  mode in which a site whose homologous region in a genome is unresolved is scored
  UNKNOWN and dropped from that guide's denominator instead of counted as a mismatch.
* **Denominator.** The same guides re-measured against every HSV-1 nucleotide record
  in GenBank that is shorter than a genome, plus the near-full-length "partial
  genome" records — 4,945 records in this run, cached under `data/raw/genes/`.

Homology for the last two is established without an aligner, by exact *k*-mer anchors
and greedy colinear chaining (the seed-and-chain half of BLAST/minimap2). **A record
counts toward a guide only if one chain has an anchor entirely left of the guide
footprint and another entirely right of it** — those two anchors bracket the
homologous stretch exactly, with no padding constant. Everything between two anchors
of one chain counts as covered *even where it is divergent*, so a guide is never
dropped merely for sitting in a variable region; that bias would have manufactured a
flattering answer. A non-matching pair is then ABSENT (bracket short and fully
resolved — real variation), or UNKNOWN (bracket contains an `N`, or is implausibly
wide, or the record does not reach both sides). Only PRESENT and ABSENT enter the
corrected denominator.

**6. SaCas9 head-to-head (opt-in, `--benchmark-sacas9`).** Re-enumerates the ICP0/RL2
and ICP27/UL54 site space in Amrani et al.'s own grammar (20 nt spacer + `NNGRRT`),
recovers all four of their published guides from it — asserted in
`tests/test_benchmark.py`, because if they were not recovered every rank would be
meaningless — and ranks them against the alternatives on both conservation
denominators, including the *joint* conservation of every ICP0 × ICP27 pair.
[Details](#sacas9-head-to-head-with-amrani-et-al-2024).

**7. Escape model (opt-in, `--escape`).** Computes P(escape) for arbitrary guide sets
over the empirical 183-genome presence matrix, with a per-site NHEJ repair term scored
against a per-codon tolerance profile *measured* from cross-strain amino-acid
variation. Deterministic closed form plus a seeded Monte Carlo cross-check.
[Details](#escape-probability-stage-7).

**8. Human off-target screening (opt-in, `--offtarget`).** Exhaustive substitution
search of the guide set against GRCh38, both strands, ≤ 4 mismatches, `NNGRRN` with
`NNGRRT` flagged as its subset, coding-exon annotation from the Ensembl GTF. Runs only
if the genome cache is already present. [Details](#human-off-target-screening-stage-8).

---

## Output columns (`results/guides_ranked.tsv`)

| Column | Meaning |
|--------|---------|
| `rank` | 1 = best. Sorted by conservation, then filter pass, then poly-T, homopolymer run, distance of GC from the centre of the allowed band, then gene and coordinate. |
| `guide_id` | Stable id: `<gene>_<ref_start><strand>` |
| `gene`, `gene_product` | From the reference annotation. `A\|B` if one 23-mer is shared by two genes. |
| `strand` | Strand the protospacer lies on (`+`/`-`) |
| `position_in_reference`, `ref_start`, `ref_end` | 1-based inclusive plus-strand coordinates spanning the whole 23-mer |
| `cut_site_ref` | Predicted blunt cut position, reported as `ref_end - pam_length - 2` (plus strand) / `ref_start + pam_length + 1` (minus). For SpCas9 that is `e-5`/`s+4`, exactly as in the original implementation. Note this sits one base inside the canonical blunt cut (which is between protospacer positions 17 and 18, i.e. `e-6`/`e-5`); the original convention was kept deliberately so that SpCas9 output stays reproducible. The column is informational and is not used in ranking, filtering or any analysis in this repository |
| `protospacer`, `pam`, `target_23mer` | spacer / PAM / the concatenation actually searched for. `target_23mer` keeps its historical name for back-compatibility but its **length is nuclease-dependent**: 23 nt for SpCas9, 26-27 nt for SaCas9 |
| `nuclease`, `pam_pattern` | provenance columns written by stage 2 (`spcas9`/`NGG`, `sacas9`/`NNGRRT`, ...). Stage 3 reads them back, so conservation can never be scored under a different grammar than the one that enumerated the sites |
| `gc_content`, `gc_in_range` | GC fraction of the protospacer; within `--gc-min`/`--gc-max` (default 0.35–0.75) |
| `max_homopolymer_run` | Longest single-base run in the protospacer |
| `has_polyT` | **`TTTT` or longer in the protospacer.** RNA polymerase III terminates at a run of >=4 T, so a U6/H1-driven sgRNA containing one is transcribed truncated and is non-functional. A hard exclusion, not a soft penalty. |
| `conservation_fraction`, `n_strains_present`, `n_strains_total` | The conservation result |
| `passes_filters` | AND of: conservation >= `--min-conservation` (default 1.0), not `has_polyT`, `max_homopolymer_run` <= `--max-homopolymer`, `gc_in_range` |
| `n_reference_copies`, `ref_all_positions` | Repeat-region multiplicity |
| `absent_in` | Up to 25 accessions lacking the site — makes any non-perfect score auditable |

---

## Results actually observed

Full HSV-1 run under the SpCas9 defaults, 2026-09-07, query as above. **These are the
real numbers produced by the run; nothing is padded or extrapolated.** The same run
under SaCas9 is [below](#the-same-run-under-sacas9).

* NCBI reported and returned **183** complete HSV-1 genomes (147,898–159,092 bp).
* **4,777** unique candidate 23-mers across the 7 target genes (RL2 713, UL19 959,
  UL29 783, UL30 794, UL5 483, UL52 673, UL54 372).
* **833 / 4,777 (17.4%)** are present in **all 183** genomes.
* **644** are perfectly conserved *and* pass the poly-T, homopolymer and GC filters.

Cumulative conservation distribution:

| conservation_fraction | guides |
|---|---|
| = 1.00 | 833 |
| >= 0.99 | 1,465 |
| >= 0.95 | 3,431 |
| >= 0.90 | 3,892 |
| >= 0.50 | 4,442 |
| = 0.00 | 0 |

Perfect conservation by gene — the ordering is biologically sensible, with the
essential structural gene most conserved and the GC-rich immediate-early repeat gene
least:

| Gene | guides | perfectly conserved | % |
|------|--------|---------------------|---|
| UL19 (VP5) | 959 | 247 | 25.8 |
| UL5 | 483 | 94 | 19.5 |
| UL52 | 673 | 121 | 18.0 |
| UL29 (ICP8) | 783 | 124 | 15.8 |
| UL54 (ICP27) | 372 | 52 | 14.0 |
| UL30 (Pol) | 794 | 104 | 13.1 |
| RL2 (ICP0) | 713 | 91 | 12.8 |

Also observed: 119 guides carry a poly-T terminator and 1,366 fall outside the GC
band, which is why the perfectly-conserved count drops from 833 to 644 after
filtering.

### The same run under SaCas9

`python run_pipeline.py --nuclease sacas9 --skip-fetch`, same 183 genomes, same 7
genes, outputs under `results/sacas9/`:

| nuclease / PAM | candidates | perfectly conserved | % | pass all filters |
|---|---|---|---|---|
| SpCas9 `NGG` (20 nt) | 4,777 | 833 | 17.4 | 644 |
| SaCas9 `NNGRRT` (21 nt) | 448 | 67 | 15.0 | 45 |
| SaCas9 `NNGRRT` (20 nt) | 449 | 69 | 15.4 | 58 |
| SaCas9 `NNGRRN` (21 nt) | 3,526 | 496 | 14.1 | 338 |

**The headline is the first column, not the third.** Requiring `NNGRRT` costs an
order of magnitude in candidate sites (448 vs 4,777) while the *fraction* that is
perfectly conserved barely moves. Anyone designing a SaCas9 therapeutic against
HSV-1 is choosing from roughly a tenth of the options an SpCas9 design has, which is
the structural reason a SaCas9 paper can end up with a mediocre guide without doing
anything wrong procedurally.

Perfect conservation by gene under SaCas9 `NNGRRT` (21 nt), for comparison with the
SpCas9 table above:

| Gene | guides | perfectly conserved | % |
|------|--------|---------------------|---|
| UL19 (VP5) | 77 | 16 | 20.8 |
| UL52 | 77 | 15 | 19.5 |
| UL5 | 57 | 11 | 19.3 |
| RL2 (ICP0) | 46 | 6 | 13.0 |
| UL29 (ICP8) | 63 | 7 | 11.1 |
| UL30 (Pol) | 93 | 10 | 10.8 |
| UL54 (ICP27) | 35 | 2 | 5.7 |

---

## Robustness: what the denominator is worth

Stage 5, run 2026-09-07. Full write-up with every table in
`results/robustness_report.md`; the numbers below are the headline findings, and they
include the ones that do not flatter the result.

### A — effective sample size

The 183 genomes behave like **132** independent ones at 99.9% estimated identity
(24 multi-member clusters; the largest, 9 members, is the strain-17 family:
`NC_001806`, `X14112`, `BK012101`, `JN555585`, three `MN1593xx` and two `OZ34xxxx`).
Recomputing on one representative per cluster gives **857 perfectly conserved guides
over 132 genomes**, versus 833 over 183.

Note the sign. De-duplication *raises* the count, because a duplicate contributes no
new sequence but can contribute new assembly noise. That is a caution about
interpreting the number, not a licence to prefer it: pairwise identity across the set
runs 98.1–100%, and single linkage at 99.5% would collapse 183 genomes to 20 and
"raise" conservation to 1,948 guides. The threshold matters more than the result, so
the whole sweep is reported rather than one chosen value.

### B — how much of "not conserved" is "not sequenced"

Ambiguity sweep (both columns must always be quoted together):

| `--max-ambiguous-fraction` | genomes | perfectly conserved |
|---|---|---|
| none | 183 | 833 |
| 0.01 | 169 | 892 |
| 0.001 | 157 | 945 |
| 0.0001 | 141 | 970 |
| 0 | 83 | 1,294 |

Of the 874,191 (guide, genome) pairs, 787,938 are PRESENT, 84,007 are ABSENT with the
homologous region demonstrably present and fully resolved, and only 2,246 (0.26%) are
genuinely UNKNOWN. Scoring those as UNKNOWN rather than MISMATCH moves perfect
conservation from **833 to 931**; 98 guides are perfect only under that reading.

So the honest statement is that assembly ambiguity accounts for a **~12% relative**
change in the perfectly-conserved count, not for the bulk of the non-conservation.
84,007 ABSENT calls are real sequence variation.

### C — 183 complete genomes vs a per-gene corpus

4,945 further HSV-1 records were downloaded (4,552 sub-genomic + 393 near-full-length
"partial genome"), 4,313 of which anchor-map to the reference.

**The most important finding is about the corpus, not the guides.** Only **1,021 of
the 4,552** sub-genomic records touch any of the seven target genes at all — the rest
are thymidine-kinase, glycoprotein-G and glycoprotein-B typing fragments from
unrelated surveys. For six of the seven genes the denominator is therefore dominated
by the 393 partial genomes, which are the *same kind of evidence* as the 183 complete
genomes, separated from them only by a title convention. Only **UL30** has a real
clinical-amplicon corpus (574 sub-genomic records), because UL30 is what gets
sequenced for aciclovir-resistance genotyping.

Fate of the 833 perfectly-conserved guides, by corpus tier:

| tier | guides with >=10 records | of the 833: <0.99 | <0.95 | <0.70 | still 1.00 | Pearson r |
|---|---|---|---|---|---|---|
| whole corpus | 4,777 | 119 | 0 | 0 | 183 | 0.992 |
| sub-genomic records only | 4,769 | 168 | **15** | 0 | 586 | 0.929 |
| partial genomes only | 4,777 | 145 | 0 | 0 | 268 | 0.992 |
| no coverage correction | 4,777 | 833 | 823 | 0 | 0 | 0.985 |

**The answer to the headline question: 15 of 833 (1.8%).** On the only tier that is
genuinely independent of the complete-genome set, 15 guides that look perfect at
n=183 fall below 95%; none falls below 70%. They are 13 in UL19, 1 in UL29, 1 in
UL52. Minimum over the perfect set is 0.889.

**The largest single effect in the entire analysis is the coverage correction.**
Without it — dividing by every record considered rather than by the records that
actually span the site, which is how the competing paper describes its calculation —
823 of the 833 fall below 95% and mean conservation drops from 0.90 to 0.82. That is
almost entirely artefact: a 600 bp amplicon does not contain a guide 40 kb away.
Which cuts against us as much as for us: it means published conservation percentages
computed that way are probably *understated*, not that ours are better.

Statistical power is the other caveat. Median sub-genomic records per guide is ~40
outside UL30 (326) and UL5 (107). The exact binomial 95% lower bound on 30/30 is
about 0.88, so "100% over the sub-genomic tier" is a much weaker claim for UL19 or
RL2 than for UL30.

### D — the published competitor guides, scored by our method

The four SaCas9 guides of Amrani et al. 2024 (Table 1), scored with our own exact
26-mer (20 nt spacer + NNGRRT PAM) presence test against our own corpora. This is not
a quotation of their numbers; it is the same measurement applied to their sequences.

| guide | complete genomes (n=183) | gene-level corpus |
|---|---|---|
| ICP0g1 | 0.984 | 0.937 |
| ICP0g2 | 0.847 | 0.807 |
| ICP27g1 | 0.978 | 0.973 |
| ICP27g2 | 0.973 | 0.937 |

### Verdict

"Complete genomes are the better denominator" is **defensible but must be
qualified**, and the qualification is not the one we expected. The 183-genome
estimate is not meaningfully inflated by redundancy or ambiguity, and it survives a
widened denominator with 1.8% attrition. But for six of seven target genes GenBank
simply does not hold an independent per-gene corpus to test against, so for those
genes "we validated against thousands of sequences" would be an overstatement — the
extra sequences are mostly more whole genomes. The reportable claim is a
two-denominator one: perfectly conserved over 183 complete genomes **and** >=95% over
the gene-level corpus, with the per-guide denominator stated.

---

## SaCas9 head-to-head with Amrani et al. 2024

Stage 6 (`--benchmark-sacas9`, full write-up in
`results/sacas9_benchmark_report.md`). Amrani et al. 2024
(*Mol Ther Methods Clin Dev* 32:101303, EBT-104) selected four SaCas9 guides -- two
in ICP0/RL2, two in ICP27/UL54 -- on a `>70%` conservation criterion against a
Feb-2022 ViPR CDS snapshot plus an off-target count, and took **ICP0g2 + ICP27g1**
forward as the clinical pair. Stage 5 scored those four guides with our method and
found ICP0g2 at 0.847. Stage 6 asks the question that number provokes: *within the
SaCas9 site space they were actually working in, what else was available?*

The candidate pool is enumerated in their exact grammar -- 20 nt spacer + `NNGRRT`,
the site definition printed in their Table 1 -- inside the annotated CDS of the same
two genes. All four published guides are recovered from that pool, which is asserted
in `tests/test_benchmark.py`; if they were not, every rank below would be
meaningless.

### Where their guides rank

| guide | gene | rank in gene pool | conservation (n=183) | gene-level corpus | GC | passes filters |
|---|---|---|---|---|---|---|
| ICP0g1 | ICP0/RL2 | 12 / 46 | 0.984 | 0.983 | 0.80 | no (GC 0.80, homopolymer 5) |
| **ICP0g2** (lead) | ICP0/RL2 | **31 / 46** | **0.847** | 0.898 | 0.60 | yes |
| **ICP27g1** (lead) | ICP27/UL54 | **12 / 35** | **0.978** | 0.991 | 0.60 | yes |
| ICP27g2 | ICP27/UL54 | 19 / 35 | 0.973 | 0.995 | 0.65 | yes |

### The central result

* **30 of the 46 `NNGRRT` sites in ICP0 are better conserved than ICP0g2.**
  **12** of those also pass every standard filter (no poly-T terminator,
  homopolymer run <= 4, GC 0.35-0.75), and **4** are present in all 183 genomes.
  All 12 lie in both ICP0 repeat copies, so they preserve the three-DSB design that
  motivates the paper.
* Their `>70%` threshold admits **69 of the 81** sites across both genes. The
  criticism is not that they broke their own rule -- all four guides clear it -- but
  that the rule barely discriminates.
* **Pairs.** Their lead pair is intact at both sites in **151 / 183 genomes
  (0.825)** -- *below either guide's own conservation*, because the two guides fail
  in different isolates. 258 filter-passing pairs beat it and 8 are intact in all
  183. Keeping ICP27g1 exactly as published and swapping only the ICP0 guide for
  `RL2_3441+` — perfectly conserved on *both* denominators, GC 0.50, local GC 0.635,
  present in both ICP0 copies — takes the pair from 0.825 to **0.978** (28 more
  isolates covered).
  This joint number is invisible to a per-guide conservation table, which is the
  form in which their method reports its selection.

  > **`RL2_3441+` is WITHDRAWN as the named recommendation.** The sentence above is
  > kept because the claim was made here, and its conservation arithmetic is still
  > correct. Stage 8 screened it against GRCh38 and found it *dirtier than ICP0g2*:
  > 19 NNGRRT sites at ≤ 4 mismatches against ICP0g2's 5, and 6 at ≤ 3 mismatches
  > against ICP0g2's 0, one of them inside *ABL1*. Several ICP0 sites tie exactly on
  > the conservation figure quoted, so naming this one was always a tie-break, not a
  > finding. See [the reconciled recommendation](#the-reconciled-recommendation) and
  > `results/recommendation.md`.

### The feasibility check, which does not go the way you would expect

ICP0 is GC-rich and repeat-associated, so the obvious objection is that the
better-conserved sites are unusable GC-extremes. They are not. The 12 filter-passing
alternatives span GC 0.50-0.75 with a median local (200 bp) GC of 0.69, against
**0.84 local GC around ICP0g2 itself** -- a figure that independently corroborates
Amrani et al.'s own report that the ICP0g2 site could not be amplified or sequenced
because of ~85% GC. On the one feasibility axis visible in the data, and the one they
themselves flagged as a problem, the alternatives are *better*, not worse.

### What this does not show

Stage 6 ranks on cross-isolate conservation and on nothing else. Two of the axes it
cannot see have since been measured, and one of them changed the answer:

* **Off-target burden is now measured** — stage 8 (`src/offtarget.py`), against the
  Ensembl 116 GRCh38 primary assembly, ≤ 4 mismatches, both PAM variants, both
  strands, with coding-exon annotation. It **overturned this stage's named guide**
  (`RL2_3441+`) and vindicated their ICP27 choice. It does not model DNA/RNA bulges,
  which Cas-OFFinder does, and it includes no cell-based validation, whereas Amrani
  et al. ran GUIDE-seq.
* **Escape probability is now modelled** — stage 7 (`src/escape.py`). It turns out to
  carry almost no information about *which* ICP0 guide to choose, for a reason worth
  knowing; see [the reconciled recommendation](#the-reconciled-recommendation).
* **Activity is still not predicted.** Conservation says a site exists in an isolate,
  not that SaCas9 cuts it efficiently. Amrani et al. screened six pairwise
  combinations in Vero cells and chose on measured antiviral activity. AAV packaging,
  synthesis feasibility and unpublished screening failures remain invisible here.

The defensible claim from this stage alone is therefore narrow and specific: **on
cross-isolate conservation -- the axis they themselves selected on -- their lead ICP0
guide is beaten by a quarter of its own gene's SaCas9 sites, and their lead pair by
258 filter-passing pairs.** Whether any of those alternatives is a better *drug* is
not settled by this analysis, and the report says so at length. What *is* now settled
is that conservation alone picks the wrong one of them.

---

## Escape probability (stage 7)

Full write-up in `results/escape_model_report.md`; `src/escape.py`, opt-in with
`--escape`. Amrani et al. justify a two-guide design by citing HIV CRISPR-escape
literature and asserting that targeting "two or more" sites prevents escape. They
never quantify it. This stage does.

A viral genome escapes multiplex editing if it ends up uncleavable at *every*
targeted site while still encoding functional essential proteins. Per site that
decomposes into a **measured** term — the per-genome presence matrix across the 183
complete genomes, so correlated failure across sites is reproduced rather than
assumed away — and a **repair** term `q_i`, the probability that NHEJ produces an
indel that destroys recognition and leaves the protein working. Frameshifts in an
essential gene are not escape. In-frame indels are scored against a per-codon
tolerance profile *measured from cross-strain amino-acid variation in the same 183
genomes*: a codon that varies among viable clinical isolates is demonstrably
tolerant. The DEFAULT NHEJ indel-length spectrum is an explicitly labelled
ASSUMPTION with an override and a sensitivity sweep; no citation is invented for it.
The sweep also includes one spectrum that is measured rather than posited --
`neuronal_nhej`, the pooled CRISPResso2 net-length histogram deposited with Ramadoss
et al. 2025 Fig. 1d, from post-mitotic human iPSC-derived neurons (in-frame fraction
0.090 against 0.199 for the default). It lowers escape, leaves every selected guide
set unchanged, and stays a SCENARIO rather than becoming the default, because it was
measured with SpCas9 at one human locus in culture and not with SaCas9 in a
trigeminal ganglion. See `refs/ramadoss2025_notes.md`.

| set | k | joint conservation | P(escape) per exposed genome |
|---|---|---|---|
| Amrani lead pair (ICP0g2 + ICP27g1), as published | 2 | 0.825 | 9.720e-04 |
| best k=2 from the same SaCas9 site space | 2 | 1.000 | 2.260e-05 |
| best k=3 | 3 | 1.000 | 1.075e-07 |

**The dominant term at k=2 is guide *choice*, not guide *count*.** 96% of the
published pair's escape probability comes from the 32 isolates that have already lost
ICP0g2; the same two-guide architecture with a better-conserved ICP0 guide clears the
same thresholds. Two findings cut against the tidy version of that story:

* **The sampling resolution of the corpus, not the NHEJ model, sets the floor.** A
  site absent in 0 of 183 isolates is still consistent with a population absence
  frequency up to 0.0162 (Clopper–Pearson 95%, the rule of three). At that bound a
  single site cannot be certified below ~1.6e-02 and a pair below ~2.6e-04. Any claim
  that a two-guide design achieves escape below about 3e-04 is an extrapolation
  beyond the evidence — **including claims made by this model**.
* **`theta`, the fraction of viral genomes never exposed to an active nuclease, is
  reported separately and never folded into P(escape)**, because guide count cannot
  touch it. The only *in vivo* on-target measurement in Amrani et al. is ~1% indels
  in trigeminal ganglia.

---

## Human off-target screening (stage 8)

Full write-up in `results/offtarget_report.md`; `src/offtarget.py`, opt-in with
`--offtarget`. Ensembl release-116 GRCh38 primary assembly (unmasked, not
soft/hard-masked; primary, not toplevel), 3,099,750,718 bases over 194 sequences,
both strands, up to 4 mismatches, PAM `NNGRRN` with `NNGRRT` flagged as its subset,
coding-exon annotation from the Ensembl 116 GTF for sites at ≤ 3 mismatches.

The search is exact, not heuristic. PAM sparsity discards 15/16 of the genome with
three vectorised slice comparisons; the pigeonhole principle (m+1 chunks, at least
one of which must match exactly) prunes ~98% of what is left; survivors are re-read
from the genome and counted exactly. `tests/test_offtarget.py` proves it the way
`tests/test_nuclease.py` proves the stage-3 scanner — a structurally independent
brute-force scan over the whole of chromosome 21, asserting set identity of the hits.
**Bulges are not modelled**, which Cas-OFFinder does model; that is stated in the
report in those words.

It was run because it was capable of invalidating this repository's own
recommendation, and it partly did:

| | ICP0g2 (theirs) | RL2_3441+ (stage 6's pick) | RL2_5335+ |
|---|---|---|---|
| NNGRRT sites ≤ 4 mm | 5 | **19** | **3** |
| NNGRRT sites ≤ 3 mm | 0 | **6** | 0 |
| NNGRRN sites ≤ 4 mm | 56 | 78 | 41 |
| NNGRRN sites ≤ 3 mm | 5 | 9 | 2 |
| coding-exon hits ≤ 3 mm | 0 | 0 | 0 |
| HSV-1 conservation (n=183) | 0.847 | 1.000 | 1.000 |

One of `RL2_3441+`'s 3-mismatch NNGRRT sites is inside *ABL1*. Independently, the
screen corroborates Amrani et al.: they report no ≤ 3-mismatch sites for either lead
guide, and under the canonical NNGRRT PAM this screen finds exactly zero for both,
computed from a different assembly download by a different algorithm. It also
vindicates their ICP27 choice — ICP27g1 has the cleanest profile of every ICP27
candidate screened, and none of the 9 better-conserved filter-passing alternatives
matches it on all four measures.

---

## The reconciled recommendation

Stages 6, 7 and 8 each proposed a *different* ICP0 replacement, because each
optimised a different axis. Reconciling them is
[`results/recommendation.md`](results/recommendation.md), which carries the joint
table, the selection rule stated before it is applied, and an explicit
"recommendations withdrawn and why" section.

The single recommendation is **`RL2_5335+` + `ICP27g1` (`UL54_115156+`)** — one guide
changed from the published pair:

| | ICP0g2 + ICP27g1 (published) | **RL2_5335+ + ICP27g1** |
|---|---|---|
| conservation, 183 complete genomes | 0.847 | **1.000** |
| conservation, gene-level corpus | 0.898 | **1.000** |
| conservation, sub-genomic tier | 0.952 (21 records) | **1.000 (36 records)** |
| joint conservation of the pair | 0.825 (151/183) | **0.978 (179/183)** |
| P(escape) per exposed genome | 9.720e-04 | **2.515e-04** |
| NNGRRT off-targets ≤ 4 mm / ≤ 3 mm | 5 / 0 | **3 / 0** |
| NNGRRN off-targets ≤ 4 mm / ≤ 3 mm | 56 / 5 | **41 / 2** |
| coding-exon hits ≤ 3 mm | 0 | 0 |
| local GC, 200 bp | 0.840 | **0.725** |

Two results from that reconciliation are worth stating here because they are easy to
get wrong:

* **Escape probability contributes no information at all to this choice.** All four
  perfectly conserved, filter-passing ICP0 candidates have identical `q` to ten
  significant figures and identical joint conservation with ICP27g1, so they are
  *exactly tied* at P(escape) = 2.515e-04. Stage 7's apparent preference for
  `RL2_5080+` was an `argmin` over four values that differ only in the 12th
  significant digit — floating-point summation order, not biology.
* **There is no candidate that dominates on every axis.** `RL2_5335+` carries two
  NNGRRN sites at ≤ 3 mismatches where `RL2_3364+` and `RL2_5080+` carry one, so the
  Pareto set has three members and the choice needs one stated priority: that the
  canonical `NNGRRT` tier, which is what SaCas9 actually cleaves efficiently and what
  all four published on-target sites use, outweighs the deliberately permissive
  `NNGRRN` search tier. That priority is a biological judgement, it is stated rather
  than smuggled in, and `results/recommendation.md` reports what the answer would
  have been under the opposite one.

---

## Limitations

Read these before using any guide from this table.

1. **Off-target screening is implemented, but only for the SaCas9 benchmark guide
   set.** `src/offtarget.py` (stage 8) is no longer a stub: it screens guides against
   the Ensembl release-116 GRCh38 primary assembly, both strands, genome-wide, at up
   to 4 mismatches under `NNGRRN` and `NNGRRT`, and annotates sites at ≤ 3 mismatches
   against the Ensembl GTF. What it has actually been *run* on is the 25-guide
   stage-6 benchmark set in RL2 and UL54 (plus their 21-nt variants), because that is
   the set the head-to-head is about. **None of the 4,777 SpCas9 candidates in
   `results/guides_ranked.tsv` has been screened**, so for those guides the original
   caveat stands unchanged: *"Guides were not screened for off-target activity
   against the human genome; conservation ranking reflects on-target coverage across
   HSV-1 isolates only."* For the guides that *were* screened the correct sentence is
   the one `src/offtarget.py`'s `caveat()` emits, and it names the two things the
   screen does not do: **DNA/RNA bulges are not modelled** (Cas-OFFinder models them)
   and **no cell-based validation** (e.g. GUIDE-seq) was performed. The run log
   records which of the two sentences applies to that run.

2. **Assembly ambiguity depresses conservation — now quantified.** 100 of the 183
   genomes contain at least one non-ACGT base; 26 have more than 0.1%. A site
   overlapping an `N` cannot match exactly and is scored absent, so some "not
   conserved" calls are assembly artefacts. Re-running with
   `--max-ambiguous-fraction 0.001` drops 26 genomes and raises perfect conservation
   from **833/183 genomes to 945/157 genomes**. Stage 5 measures the effect directly
   rather than by proxy: only **2,246 of 874,191** (guide, genome) pairs are actually
   unresolvable, and N-tolerant scoring moves the count 833 -> 931. The remaining
   84,007 non-matches are real variation. Per-genome `n_ambiguous` and
   `ambiguous_fraction` are in the manifest.

3. **The genome set is not 183 independent clinical isolates — now quantified.** It
   mixes clinical isolates with long-passaged laboratory strains, and the strain-17
   family alone accounts for 9 records. Stage 5 puts the effective sample size at
   **132** at 99.9% identity. No two records are byte-identical, so nothing is
   silently de-duplicated. Contrary to the expectation recorded here originally,
   de-duplication *raises* the perfectly-conserved count (833 -> 857) rather than
   lowering it, so redundancy is not the source of optimism it was assumed to be.
   The clustering threshold, however, dominates the answer (see the sweep), so the
   effective-n figure should always be quoted with its threshold.

4. **"Complete genome" is a title convention, not a guarantee.** The default query
   trusts the submitter's title. A further **393** HSV-1 records of 145–160 kb are
   deposited as *"partial genome"* — mostly genuine near-full-length clinical
   isolates whose terminal/internal repeats were not resolved. They are excluded by
   default and admitted by `--include-partial`. Stage 5 scores them explicitly as
   their own tier: they change nothing (0 of 833 perfect guides fall below 95% on
   them), which is itself the finding — they are the same kind of evidence as the
   183, so admitting them widens the denominator without widening the evidence.

5. **The reference defines the candidate space.** Guides are enumerated only from
   `NC_001806`. A site conserved across all other isolates but absent from strain 17
   is never considered.

6. **Metadata is sparse.** Of 183 records, 144 report a strain, 108 a country and 82
   a collection date. Geographic or temporal stratification of conservation is not
   currently well supported by the available metadata.

7. **Cutting is not the whole story.** Perfect conservation and clean expression
   flags say nothing about chromatin accessibility on the incoming or latent viral
   genome or about editing efficiency. The last item on this list — the rate at which
   NHEJ repair generates cut-resistant escape variants, which is the actual
   motivation for multiplex design — is now modelled by stage 7, but *modelled* is
   the operative word: its DEFAULT NHEJ indel-length spectrum is a labelled
   ASSUMPTION rather than a measurement, its absolute escape probabilities move over
   ~3 orders of magnitude across the plausible parameter space, and with 183 genomes
   a two-guide set cannot be certified below ~2.6e-04 whatever the repair term says.
   The sweep now includes the one spectrum measured in the right cell class
   (`neuronal_nhej`, post-mitotic human neurons), which moves escape DOWN and the
   selected sets not at all. What is stable across the sweep is the ranking of guide
   sets, not the numbers.

8. **The per-gene corpus is thinner than GenBank's record counts suggest.** Only
   1,021 of 4,552 sub-genomic HSV-1 records cover any of the seven target genes, and
   outside UL30 the median guide has ~40 independent records behind it. Statements of
   the form "conserved across thousands of sequences" are not supportable for six of
   the seven genes; see the stage-5 report for the per-gene denominators.

9. **A stage-6 rank is still not a guide-selection decision — and acting as though
   it were is a mistake this repository actually made.** Stage 6 ranks on
   conservation alone and, on that basis, named `RL2_3441+`. Stage 8 then screened it
   and found it dirtier against the human genome than the guide it was proposed to
   replace, so that recommendation is withdrawn (see
   [the reconciled recommendation](#the-reconciled-recommendation) and
   `results/recommendation.md`). Two of the three axes Amrani et al. selected on are
   now visible here — conservation and off-target burden — and stage 7 adds a fourth
   they did not quantify at all. The one that remains invisible is **measured
   antiviral activity in cells**, which is what they actually chose their lead pair
   on. A rank from any single stage of this pipeline is an argument about a selection
   *metric*, not a claim that a better guide exists all things considered.

10. **Coverage determination is anchor-based, not alignment-based.** A record is
   admitted to a guide's denominator only when exact 25-mer anchors bracket the guide
   footprint within one colinear chain. This is deliberately conservative and is
   validated by `tests/test_robustness.py` on synthetic truncations, SNPs, N blocks
   and repeat duplications, but it is a heuristic: a record whose homologous region is
   real yet too divergent to anchor on either side is scored UNKNOWN rather than
   ABSENT, which biases the corrected conservation slightly upward. The
   `no coverage correction` tier in the stage-5 report brackets the effect from the
   other side.

---

## Reproducibility: the full output inventory

* Deterministic throughout: accessions sorted before `--limit`, stable sorts in every
  table, no random seeds anywhere in the pipeline.
* `results/run_log.json` records the full sorted accession list, the exact Entrez
  query, the retrieval timestamp, the command line, Python and package versions,
  which stages executed, and the summary statistics.
* `data/manifest.tsv` records per genome: accession, length, ambiguous base count,
  strain, country, collection date, host, isolation source, definition, completeness
  label, SHA-256 of the downloaded FASTA, plus the query and retrieval date.
* Downloads are cached in `data/raw/`, so re-runs never re-download and are offline.
* If NCBI returns fewer records than its own count reports, the shortfall is logged as
  a warning and the manifest contains exactly what was actually retrieved. Nothing is
  padded, imputed or invented.

Stage 5 adds to the record: `results/robustness_report.md` (human-readable),
`results/robustness_summary.json` (machine-readable, and embedded in
`results/run_log.json` when the stage runs inside the pipeline), and the supporting
tables `robustness_redundancy.tsv`, `robustness_effective_n.tsv`,
`robustness_ambiguity_sweep.tsv`, `robustness_conservation_modes.tsv`,
`robustness_gene_level.tsv`, `robustness_gene_level_by_gene.tsv`,
`robustness_gene_records.tsv`, `robustness_benchmark_guides.tsv`, plus the per-record
`data/gene_corpus_manifest.tsv`. Both Entrez queries used to build the gene-level
corpus are written into the report and the JSON summary.

Stage 7 adds `results/escape_model_report.md`, `results/escape_summary.json` and the
supporting tables `escape_site_parameters.tsv` (per-site q and cut codon),
`escape_k_curve.tsv`, `escape_k_curve_one_per_gene.tsv`, `escape_guide_sets.tsv`,
`escape_codon_tolerance.tsv` (the measured per-codon cross-strain variation),
`escape_sensitivity.tsv`, `escape_sensitivity_sets.tsv` and `escape_tolerance_qc.tsv`.
It is deterministic — closed form throughout, with a seeded Monte Carlo cross-check
that must agree with it — and runs offline from the caches populated by stages 1
and 5.

Stage 8 adds `results/offtarget_report.md`, `results/offtarget_summary.json`,
`results/offtarget_summary.tsv` (per-guide counts at each mismatch level under both
PAM variants), `results/offtarget_sites.tsv` (every site found) and
`results/offtarget_annotated.tsv` (genomic context for sites at ≤ 3 mismatches). The
GRCh38 download is verified three independent ways — Ensembl's published BSD `sum`,
the gzip CRC-32/ISIZE trailer, and a SHA-256 recorded in
`data/genome/genome_manifest.json` — and the assembly release is pinned, not
"current", so a floating release cannot silently change the numbers.

`results/recommendation.md` reconciles stages 6, 7 and 8 into one recommendation. It
is version-controlled because it is a written argument across stages rather than a
regenerable table; the tracking policy for the rest of `results/` is under
[What is in `results/`, and what is not](#what-is-in-results-and-what-is-not). Its
numbers are not typed by hand: `python src/reconcile.py` joins the three stages' outputs on `guide_id`,
computes the one quantity none of them produced — P(escape) for each ICP0 candidate
paired with ICP27g1 — evaluates Pareto dominance mechanically, applies the document's
selection rule from constants in that module, and writes
`results/recommendation_table.tsv`. It **refuses to report anything unless its rebuilt
escape model reproduces stage 7's own Amrani lead-pair P(escape) to within 1e-12**.

Stage 6 adds `results/sacas9_benchmark_report.md`,
`results/sacas9_benchmark_summary.json` (also embedded in `results/run_log.json` when
the stage runs inside the pipeline), `results/sacas9_benchmark_pool.tsv` (every
SaCas9 site in RL2 and UL54 with both conservation measures, GC context, filter flags
and rank) and `results/sacas9_benchmark_pairs.tsv` (all 1,610 ICP0 x ICP27 pairs with
joint conservation). It runs entirely offline from the caches populated by stages 1
and 5.

`data/raw/` and `data/genome/` are gitignored, as are the bulk per-guide and per-site
tables under `results/`: they are regenerable, and `results/CHECKSUMS.sha256` pins
every one of them by SHA-256 so a regenerated copy can be checked byte-for-byte. The
provenance records `data/manifest.tsv`, `data/gene_corpus_manifest.tsv` and
`results/run_log.json` **are** version-controlled — those files plus this repository
fully determine the results. The full policy, and the reasoning behind each side of
it, is under [What is in `results/`, and what is not](#what-is-in-results-and-what-is-not).

One caveat on `results/run_log.json`: it is the log of the invocation that produced
the pinned stage-1–4 and stage-7 tables (`--skip-fetch --force --escape`), so its
`robustness`, `sacas9_benchmark` and `offtarget` keys are `null`. Those stages were
run separately and their own records are `results/robustness_summary.json`,
`results/sacas9_benchmark_summary.json` and `results/offtarget_summary.json`, each
with its own generation timestamp. The corpus pin — query, retrieval timestamp and
the full sorted accession list — is in `run_log.json` and is the same for all of
them.

---

## Source layout in detail

```
run_pipeline.py          single entrypoint, stages 1-4 plus opt-in 5-8, run log
verify.py                headline-number verification; fast path is offline and < 1 s
src/common.py            paths, logging, Entrez config, rate limit, retry, seq utils
src/fetch_genomes.py     stage 1 - NCBI retrieval + manifest
src/extract_guides.py    stage 2 - SpCas9 site enumeration from GenBank annotation
src/conservation.py      stage 3 - alignment-free exact-match conservation scoring
src/report.py            stage 4 - flags, ranking, guides_ranked.tsv
src/robustness.py        stage 5 - redundancy, N-sensitivity, denominator sensitivity
src/benchmark_sacas9.py  stage 6 - SaCas9 head-to-head vs Amrani et al. 2024
src/escape.py            stage 7 - multiplex escape-probability model
src/offtarget.py         stage 8 - human GRCh38 off-target screening (implemented)
src/reconcile.py         joins stages 6-8 into one table; the arithmetic behind
                         results/recommendation.md (a reader, not a stage)
src/nuclease.py          nuclease/PAM model (IUPAC, both strands, scanner anchors)
tests/test_core.py       offline unit tests (no network)
tests/test_robustness.py offline unit tests for stage 5 (no network)
tests/test_nuclease.py   PAM model + scanner-equivalence proof for every PAM
tests/test_benchmark.py  offline unit tests for stage 6 (no network)
tests/test_escape.py     offline unit tests for stage 7 (no network)
tests/test_offtarget.py  offline unit tests for stage 8 (no network, no genome needed)
tests/test_reconcile.py  the selection rule of results/recommendation.md, pinned
tests/test_figures.py    figures 1-5 render byte-identically from the pinned tables
src/figures.py           figures 1-5 (.png and .pdf), rendered from results/
results/recommendation.md  the reconciled single recommendation across stages 6-8
results/CHECKSUMS.sha256   SHA-256 of every deterministic output of the pinned run
requirements.txt         pinned, installed and tested on CPython 3.14.5
LICENSE                  MIT, unmodified so GitHub detects it
NOTICE.md                third-party licences covering refs/
CITATION.cff             "Cite this repository"; read by GitHub and by Zenodo
.zenodo.json             deposition metadata for the archived DOI
PUBLISH.md               the author's checklist for the GitHub push and Zenodo mint
```
