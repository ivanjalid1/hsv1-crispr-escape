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
or truncated". See [Limitations](#limitations).

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

# offline unit tests (no network)
python tests/test_core.py
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

## Limitations

Read these before using any guide from this table.

1. **No off-target screening.** `src/offtarget.py` is a clearly marked stub. No guide
   here has been checked against the human genome. Any manuscript must state:
   *"Guides were not screened for off-target activity against the human genome;
   conservation ranking reflects on-target coverage across HSV-1 isolates only."*
   The module documents the phase-2 plan and the fact that a 3.1 Gb mismatch-tolerant
   search cannot honour the pure-Python constraint — that conflict is an explicit
   phase-2 decision, not an oversight.

2. **Assembly ambiguity depresses conservation.** 100 of the 183 genomes contain at
   least one non-ACGT base; 26 have more than 0.1%. A site overlapping an `N` cannot
   match exactly and is scored absent, so some "not conserved" calls are assembly
   artefacts rather than real sequence variation. Re-running with
   `--max-ambiguous-fraction 0.001` drops 26 genomes and raises perfect conservation
   from **833/183 genomes to 945/157 genomes** — a real effect worth reporting in the
   methods. Per-genome `n_ambiguous` and `ambiguous_fraction` are in the manifest.

3. **The genome set is not 183 independent clinical isolates.** It mixes clinical
   isolates with long-passaged laboratory strains, and strain 17 alone is deposited
   four times (`NC_001806`, `X14112`, `BK012101`, `JN555585`) with near-identical
   sequence. No two records are byte-identical, so nothing is silently de-duplicated,
   but the effective sample size is smaller than 183 and conservation fractions are
   correspondingly optimistic. A curated, de-duplicated subset is a sensible phase-2
   refinement.

4. **"Complete genome" is a title convention, not a guarantee.** The default query
   trusts the submitter's title. A further ~393 HSV-1 records of 145–160 kb are
   deposited as *"partial genome"* — these are mostly genuine near-full-length
   clinical isolates whose terminal/internal repeats were not resolved. They are
   excluded by default (that is what makes the "complete genomes" claim defensible)
   and admitted by `--include-partial`, which roughly triples the dataset. Every
   record's own label is preserved in the `completeness_label` manifest column.

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

`data/raw/`, `results/` and the run-specific manifest are gitignored: they are
regenerable, and a stale or smoke-test copy in version control would be misleading.
Archive `data/manifest.tsv` and `results/run_log.json` alongside a manuscript — those
two files plus this repository fully determine the results.

---

## Layout

```
run_pipeline.py          single entrypoint, all four stages, run log
src/common.py            paths, logging, Entrez config, rate limit, retry, seq utils
src/fetch_genomes.py     stage 1 - NCBI retrieval + manifest
src/extract_guides.py    stage 2 - SpCas9 site enumeration from GenBank annotation
src/conservation.py      stage 3 - alignment-free exact-match conservation scoring
src/report.py            stage 4 - flags, ranking, guides_ranked.tsv
src/offtarget.py         STUB - human off-target screening, phase 2
tests/test_core.py       offline unit tests (no network)
requirements.txt         pinned, installed and tested on CPython 3.14.5
```
