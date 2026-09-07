# Conserved CRISPR-Cas9 target sites in HSV-1

Phase 1 of a computational study on escape-resistant multiplex guide RNA design for
herpes simplex virus. This pipeline identifies SpCas9 target sites in essential
HSV-1 genes that are **perfectly conserved across every publicly available complete
HSV-1 genome**, using an alignment-free exact-match method.

It is deliberately dependency-light: pure Python plus Biopython, pandas and numpy.
No conda, no WSL, no MAFFT/MUSCLE/Clustal, no compiled tooling beyond what pip
wheels provide.

---

## Method

**Conservation is defined as exact substring presence, not alignment identity.**

For every candidate guide we take the 23-mer `protospacer(20) + PAM(NGG, 3)` and ask,
for each downloaded genome, a single yes/no question: *does this exact 23-mer occur
verbatim in that genome, on either strand?*

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
2. Every SpCas9 target ends in `NGG`, so on the plus strand it must contain `GG` at
   offset 21, and its reverse complement must start with `CC`. Rather than testing
   all L positions we jump between occurrences of `GG` and `CC` using `str.find`,
   which runs in optimised C. Only those candidate positions produce a slice and a
   dict lookup.

Work is `O(S x L)` with a small constant, memory is `O(G)`. Measured: **~95x faster
than brute force**, scoring 4,777 guides against 183 genomes in about 3 seconds.
`tests/test_core.py::test_scan_matches_bruteforce` proves on randomised sequences
that the jump scan returns *exactly* the same result set as a naive both-strand
search — the optimisation is not an approximation.

---

## Install

Verified on **Windows 11, CPython 3.14.5 (64-bit), pip 26.1.1**.

```powershell
python -m venv .venv
.venv\Scripts\python.exe -m pip install --only-binary=:all: -r requirements.txt
```

`--only-binary=:all:` is recommended: it makes pip fail loudly if a cp314 wheel is
ever missing rather than silently attempting a source build that would need a C
toolchain.

Python 3.14 wheel availability was tested, not assumed. All three scientific
packages ship prebuilt `cp314-win_amd64` wheels at the pinned versions:

| Package   | Version | Wheel                                     |
|-----------|---------|-------------------------------------------|
| numpy     | 2.5.3   | `numpy-2.5.3-cp314-cp314-win_amd64.whl`   |
| pandas    | 3.0.5   | `pandas-3.0.5-cp314-cp314-win_amd64.whl`  |
| biopython | 1.88    | `biopython-1.88-cp314-cp314-win_amd64.whl`|

Note that pandas 3.0 is a major release with behavioural changes from 2.x; the code
here is written against 3.0 and was run against it.

---

## NCBI credentials (required)

NCBI Entrez requires a contact email address. **No email address is hardcoded
anywhere in this repository.** It is read only from the environment:

```powershell
$env:NCBI_EMAIL = "you@example.org"      # required
$env:NCBI_API_KEY = "<your key>"         # optional
```

```bash
export NCBI_EMAIL="you@example.org"
export NCBI_API_KEY="<your key>"
```

If `NCBI_EMAIL` is unset the pipeline exits immediately with instructions and does
not contact NCBI. `NCBI_API_KEY` is optional; supplying one raises the NCBI rate
limit from 3 to 10 requests/second and the pipeline adjusts its throttle
accordingly. Get a key at <https://www.ncbi.nlm.nih.gov/account/settings/>.

Never commit these values. `.env` files are gitignored.

---

## Usage

```bash
# fast smoke test (5 genomes, ~10 s)
python run_pipeline.py --limit 5

# full HSV-1 run (~80 s cold, ~5 s warm from cache)
python run_pipeline.py

# admit near-full-length clinical isolates deposited as "partial genome"
python run_pipeline.py --include-partial

# exclude poorly resolved assemblies from the denominator
python run_pipeline.py --max-ambiguous-fraction 0.001

# HSV-2 instead of HSV-1
python run_pipeline.py --taxid 10310 --reference NC_001798

# a different gene set
python run_pipeline.py --genes UL30,UL29,UL54

# stage 5: robustness / denominator-sensitivity analysis (downloads ~80 MB once)
python run_pipeline.py --robustness

# ... parts A and B only, no further network access
python run_pipeline.py --robustness --skip-gene-corpus

# offline unit tests (no network)
python tests/test_core.py
python tests/test_robustness.py
```

Each stage is also runnable on its own: `python src/fetch_genomes.py --help`, etc.

Useful flags: `--force` recomputes every stage (needed if you change parameters that
the staleness check cannot see), `--skip-fetch` reuses an existing manifest without
contacting NCBI, `--refresh` bypasses the download cache, `-v` enables debug logging.

---

## Pipeline stages

| Stage | Module | Output |
|-------|--------|--------|
| 1 | `src/fetch_genomes.py` | `data/raw/*.fasta`, `data/manifest.tsv` |
| 2 | `src/extract_guides.py` | `results/guides_candidates.tsv` |
| 3 | `src/conservation.py` | `results/conservation.tsv` |
| 4 | `src/report.py` | `results/guides_ranked.tsv`, `results/summary.json` |
| 5 | `src/robustness.py` | `results/robustness_report.md` + supporting TSVs (opt-in, `--robustness`) |
| — | `run_pipeline.py` | `results/run_log.json` |

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
`NC_001806`, HSV-1 strain 17) and enumerates every SpCas9 site on both strands within
the target genes. **Gene coordinates are parsed from the annotation — no coordinates
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

---

## Output columns (`results/guides_ranked.tsv`)

| Column | Meaning |
|--------|---------|
| `rank` | 1 = best. Sorted by conservation, then filter pass, then poly-T, homopolymer run, distance of GC from the centre of the allowed band, then gene and coordinate. |
| `guide_id` | Stable id: `<gene>_<ref_start><strand>` |
| `gene`, `gene_product` | From the reference annotation. `A\|B` if one 23-mer is shared by two genes. |
| `strand` | Strand the protospacer lies on (`+`/`-`) |
| `position_in_reference`, `ref_start`, `ref_end` | 1-based inclusive plus-strand coordinates spanning the whole 23-mer |
| `cut_site_ref` | Predicted blunt cut, 3 bp 5' of the PAM |
| `protospacer`, `pam`, `target_23mer` | 20 nt / 3 nt / the concatenation actually searched for |
| `gc_content`, `gc_in_range` | GC fraction of the protospacer; within `--gc-min`/`--gc-max` (default 0.35–0.75) |
| `max_homopolymer_run` | Longest single-base run in the protospacer |
| `has_polyT` | **`TTTT` or longer in the protospacer.** RNA polymerase III terminates at a run of >=4 T, so a U6/H1-driven sgRNA containing one is transcribed truncated and is non-functional. A hard exclusion, not a soft penalty. |
| `conservation_fraction`, `n_strains_present`, `n_strains_total` | The conservation result |
| `passes_filters` | AND of: conservation >= `--min-conservation` (default 1.0), not `has_polyT`, `max_homopolymer_run` <= `--max-homopolymer`, `gc_in_range` |
| `n_reference_copies`, `ref_all_positions` | Repeat-region multiplicity |
| `absent_in` | Up to 25 accessions lacking the site — makes any non-perfect score auditable |

---

## Results actually observed

Full HSV-1 run, 2026-09-07, query as above. **These are the real numbers produced by
the run; nothing is padded or extrapolated.**

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

## Limitations

Read these before using any guide from this table.

1. **No off-target screening.** `src/offtarget.py` is a clearly marked stub. No guide
   here has been checked against the human genome. Any manuscript must state:
   *"Guides were not screened for off-target activity against the human genome;
   conservation ranking reflects on-target coverage across HSV-1 isolates only."*
   The module documents the phase-2 plan and the fact that a 3.1 Gb mismatch-tolerant
   search cannot honour the pure-Python constraint — that conflict is an explicit
   phase-2 decision, not an oversight.

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
   genome, editing efficiency, or the rate at which NHEJ repair generates
   cut-resistant escape variants — the last being the actual motivation for multiplex
   design.

8. **The per-gene corpus is thinner than GenBank's record counts suggest.** Only
   1,021 of 4,552 sub-genomic HSV-1 records cover any of the seven target genes, and
   outside UL30 the median guide has ~40 independent records behind it. Statements of
   the form "conserved across thousands of sequences" are not supportable for six of
   the seven genes; see the stage-5 report for the per-gene denominators.

9. **Coverage determination is anchor-based, not alignment-based.** A record is
   admitted to a guide's denominator only when exact 25-mer anchors bracket the guide
   footprint within one colinear chain. This is deliberately conservative and is
   validated by `tests/test_robustness.py` on synthetic truncations, SNPs, N blocks
   and repeat duplications, but it is a heuristic: a record whose homologous region is
   real yet too divergent to anchor on either side is scored UNKNOWN rather than
   ABSENT, which biases the corrected conservation slightly upward. The
   `no coverage correction` tier in the stage-5 report brackets the effect from the
   other side.

---

## Reproducibility

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

`data/raw/`, `results/` and the run-specific manifests are gitignored: they are
regenerable, and a stale or smoke-test copy in version control would be misleading.
Archive `data/manifest.tsv`, `data/gene_corpus_manifest.tsv` and
`results/run_log.json` alongside a manuscript — those files plus this repository fully
determine the results.

---

## Layout

```
run_pipeline.py          single entrypoint, all four stages, run log
src/common.py            paths, logging, Entrez config, rate limit, retry, seq utils
src/fetch_genomes.py     stage 1 - NCBI retrieval + manifest
src/extract_guides.py    stage 2 - SpCas9 site enumeration from GenBank annotation
src/conservation.py      stage 3 - alignment-free exact-match conservation scoring
src/report.py            stage 4 - flags, ranking, guides_ranked.tsv
src/robustness.py        stage 5 - redundancy, N-sensitivity, denominator sensitivity
src/offtarget.py         STUB - human off-target screening, phase 2
tests/test_core.py       offline unit tests (no network)
tests/test_robustness.py offline unit tests for stage 5 (no network)
requirements.txt         pinned, installed and tested on CPython 3.14.5
```
