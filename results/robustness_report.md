# Robustness and sensitivity analysis of HSV-1 guide conservation

Generated 2026-09-07T23:34:47Z by `src/robustness.py`. Reference NC_001806.2. 4777 candidate guides, 183 complete genomes.

Every number below was produced by the run that wrote this file. Nothing is extrapolated, and no threshold was chosen after seeing its effect: the identity and ambiguity grids are module-level constants.

## Headline

* Baseline (stage 3): **833 / 4777** guides are present in all 183 complete genomes.
* De-duplicating near-identical genomes at 99.90% estimated identity leaves an effective sample size of **132** independent genomes, and the perfectly-conserved count moves to **857** (+24).
* Scoring N as UNKNOWN rather than MISMATCH raises the perfectly-conserved count to **931** (+98); 98 guides are perfect only because unresolved positions were removed from their denominator.
* Against the gene-level corpus (4945 records, 4313 of which map to the reference), of the **833** guides that are 100% conserved at n=183 and have at least 10 covering records, **0** fall below 95% and **0** fall below the 70% threshold used by Amrani et al. 2024.
* **Caveat that changes the reading:** most of that corpus is 393 near-full-length `partial genome` assemblies, not independent per-gene records. Restricted to genuinely sub-genomic records, only **4769 of 4777** guides have a usable denominator at all; of the **833** perfect ones among them, **15** fall below 95%.

## A. Redundancy and effective sample size

Method: bottom-5000 MinHash sketch of canonical 21-mers per genome (splitmix64 mixing, no salted hashing), all-pairs Jaccard by the mash estimator, Jaccard converted to an ANI-like identity, single-linkage clustering. No alignment, no all-pairs full comparison.

Pairwise identity across the set: min 0.98134, median 0.99369, max 1.00000.

### Effective sample size vs clustering threshold

| ani_threshold | n_clusters | n_genomes_collapsed | largest_cluster | n_perfectly_conserved | pct_perfectly_conserved |
|---|---|---|---|---|---|
| 0.9999 | 153 | 30 | 5 | 833 | 17.44 |
| 0.999 | 132 | 51 | 9 | 857 | 17.94 |
| 0.998 | 107 | 76 | 18 | 906 | 18.97 |
| 0.995 | 20 | 163 | 151 | 1948 | 40.78 |
| 0.99 | 5 | 178 | 179 | 2922 | 61.17 |
| 0.98 | 1 | 182 | 183 | 4356 | 91.19 |

At the primary threshold (0.999), 75 genomes fall into 24 multi-member clusters; the rest are singletons.

Largest clusters (deposited near-duplicates):

| cluster_id | cluster_size | accession | strain | nearest_accession | nearest_identity |
|---|---|---|---|---|---|
| 1 | 9 | BK012101.1 | 17 | JN555585.1 | 1 |
| 1 | 9 | JN555585.1 | 17 | BK012101.1 | 1 |
| 1 | 9 | MN159377.1 | nan | BK012101.1 | 0.9998 |
| 1 | 9 | MN159378.1 | nan | BK012101.1 | 0.9998 |
| 1 | 9 | MN159379.1 | nan | BK012101.1 | 0.9998 |
| 1 | 9 | NC_001806.2 | 17 | BK012101.1 | 1 |
| 1 | 9 | OZ348666.1 | nan | OZ348992.1 | 1 |
| 1 | 9 | OZ348992.1 | nan | OZ348666.1 | 1 |
| 1 | 9 | X14112.1 | 17 | BK012101.1 | 0.9996 |
| 26 | 6 | MH999846.1 | HSV-H1215 | ON007142.1 | 0.9992 |
| 26 | 6 | ON007136.1 | v48_d101_oral | ON007142.1 | 1 |
| 26 | 6 | ON007142.1 | v49_d257_gen_les | ON007136.1 | 1 |
| 26 | 6 | ON007149.1 | v49_d349_gen_les | ON007136.1 | 1 |
| 26 | 6 | OP297875.1 | v49_d257_cu_gen_les | ON007142.1 | 0.9999 |
| 26 | 6 | OR771672.1 | nan | MH999846.1 | 0.999 |
| 4 | 5 | JQ673480.1 | KOS | KT899744.1 | 0.9999 |
| 4 | 5 | JQ780693.1 | KOS | JQ673480.1 | 0.9998 |
| 4 | 5 | KT887224.1 | nan | KT899744.1 | 0.9999 |
| 4 | 5 | KT899744.1 | nan | KT887224.1 | 0.9999 |
| 4 | 5 | MF156584.1 | nan | JQ673480.1 | 0.9992 |
| 37 | 5 | ON007143.1 | v46_d349_gen_les2 | ON007151.1 | 1 |
| 37 | 5 | ON007144.1 | v46_d345_gen | ON007151.1 | 1 |
| 37 | 5 | ON007151.1 | v46_d349_gen_les3 | ON007143.1 | 1 |
| 37 | 5 | ON007158.1 | v46_d349_gen_les1 | ON007151.1 | 1 |
| 37 | 5 | OP297868.1 | HSV1-v46_d349_cu_gen_les | ON007151.1 | 1 |
| 123 | 5 | PZ169677.1 | A_d6_skin_culture | PZ169679.1 | 1 |
| 123 | 5 | PZ169678.1 | A_d6_skin | PZ169677.1 | 1 |
| 123 | 5 | PZ169679.1 | A_d6_swab_culture | PZ169677.1 | 1 |
| 123 | 5 | PZ169681.1 | A_d7_eye_culture | PZ169677.1 | 0.9999 |
| 123 | 5 | PZ169683.1 | A_d7_plasma | PZ169679.1 | 0.9999 |
| 33 | 4 | MN136524.1 | McKrae | PX969520.1 | 0.9998 |
| 33 | 4 | OR723971.1 | nan | MN136524.1 | 0.9998 |
| 33 | 4 | PX763612.1 | nan | MN136524.1 | 0.9996 |
| 33 | 4 | PX969520.1 | nan | MN136524.1 | 0.9998 |

**Effect on the result.** Recomputing conservation on one representative per cluster at 0.999: 857 perfectly conserved guides over 132 genomes, versus 833 over 183. Interpretation is in the closing section.

Per-genome cluster assignments: `results/robustness_redundancy.tsv`.

## B. Ambiguity (N) sensitivity

### B1. Sweeping --max-ambiguous-fraction

Dropping poorly resolved assemblies raises conservation and shrinks the denominator at the same time. Both columns must be quoted together; either one alone is misleading.

| max_ambiguous_fraction | n_genomes | n_genomes_dropped | n_perfectly_conserved | pct_perfectly_conserved |
|---|---|---|---|---|
| none | 183 | 0 | 833 | 17.44 |
| 0.05 | 183 | 0 | 833 | 17.44 |
| 0.02 | 175 | 8 | 860 | 18 |
| 0.01 | 169 | 14 | 892 | 18.67 |
| 0.005 | 166 | 17 | 903 | 18.9 |
| 0.002 | 163 | 20 | 921 | 19.28 |
| 0.001 | 157 | 26 | 945 | 19.78 |
| 0.0005 | 153 | 30 | 950 | 19.89 |
| 0.0002 | 144 | 39 | 965 | 20.2 |
| 0.0001 | 141 | 42 | 970 | 20.31 |
| 0 | 83 | 100 | 1294 | 27.09 |

### B2. N-tolerant matching

For every (guide, genome) pair where the exact 23-mer is absent, the homologous region is located by anchor chaining and inspected. Only pairs whose homologous region is present and fully resolved are counted as real mismatches.

| classification | pairs | pct_of_all_pairs |
|---|---|---|
| present | 787938 | 90.13 |
| absent | 84007 | 9.61 |
| unknown_ambiguous | 174 | 0.02 |
| unknown_structural | 2072 | 0.237 |
| not_covered | 0 | 0 |

Strict (stage 3) perfectly conserved: **833**. N-tolerant perfectly conserved: **931**. Guides that become perfect only under N-tolerance: **98**.

Per-guide breakdown: `results/robustness_conservation_modes.tsv`.

## C. The denominator question: 183 complete genomes vs the gene-level corpus

### Corpus construction

* `subgenomic`: `txid10298[Organism:exp] AND biomol_genomic[PROP] NOT patent[PROP] AND 100:145000[SLEN]` -> NCBI reported 4552, 4552 used.
* `partial_genome`: `txid10298[Organism:exp] AND biomol_genomic[PROP] NOT patent[PROP] AND 145000:160000[SLEN] NOT "complete genome"[Title]` -> NCBI reported 393, 393 used.

4945 records are on disk under `data/raw/genes/`; 4313 of them anchor-map to the reference genome at all. Records that do not map contribute nothing and are not counted anywhere.

| record_class | records_downloaded | records_mapping_to_reference | records_covering_a_target_gene | median_bp |
|---|---|---|---|---|
| subgenomic | 4552 | 3920 | 1021 | 1131 |
| partial_genome | 393 | 393 | 393 | 151666 |

**Read that table before reading anything below it.** Most HSV-1 sub-genomic records in GenBank are not about our genes at all -- they are thymidine kinase, glycoprotein G and glycoprotein B typing fragments. The sub-genomic corpus that actually touches the seven target genes is far smaller than the raw record count suggests, and for six of the seven genes the near-full-length `partial genome` records dominate the denominator. Those are the same KIND of evidence as the 183 complete genomes -- whole assemblies separated from them only by a title convention -- so pooling everything would make the widened denominator look far more independent than it is. Every result below is therefore also reported for the sub-genomic records alone.

Gene keyword queries were deliberately NOT used: `RL2[All Fields]` returns 5 HSV-1 records and `ICP27`-style terms only 10, so a keyword-built corpus would have missed almost everything. Gene assignment comes from anchor mapping instead, so a record's own labelling is never trusted.

### Coverage: distinguishing ABSENT from NOT COVERED

A record enters a guide's denominator only if one of its anchor chains contains an anchor entirely left of the guide footprint and another entirely right of it. Those two anchors bracket the homologous stretch exactly. Everything between two anchors of the same chain counts as covered even where it is divergent, so variable sites are not silently discarded. A bracket wider than 2000 bp, or a bracket containing a non-ACGT base, is scored UNKNOWN, not ABSENT.

### Records per gene

| gene | records_covering_gene | subgenomic | partial_genome | median_record_bp |
|---|---|---|---|---|
| UL30 | 967 | 574 | 393 | 3708 |
| UL19 | 432 | 39 | 393 | 151532 |
| UL5 | 706 | 313 | 393 | 150308 |
| UL52 | 494 | 101 | 393 | 151313 |
| UL29 | 455 | 62 | 393 | 151463 |
| RL2 | 445 | 52 | 393 | 151485 |
| UL54 | 444 | 51 | 393 | 151490 |

### Conservation over the gene-level corpus

4777 of 4777 guides have at least 10 covering records (median 436, max 930).

* Mean conservation over complete genomes: 0.9013
* Mean conservation over the gene-level corpus: 0.9002
* Pearson r = 0.9919, Spearman rho = 0.8345

**Fate of the perfectly conserved guides.**

| criterion | guides |
|---|---|
| 100% at n=183 (with >= 10 covering records) | 833 |
| ...of which < 0.99 at gene level | 119 |
| ...of which < 0.95 at gene level | 0 |
| ...of which < 0.90 at gene level | 0 |
| ...of which < 0.70 at gene level | 0 |
| ...of which still exactly 100% at gene level | 183 |

Median gene-level conservation of the perfect set: 0.9958; minimum 0.96.

### The same question, one corpus tier at a time

`subgenomic_only` is the tier that matters most: it is the only one built from evidence of a different kind than the 183 complete genomes. `no_coverage_correction` is not an estimate at all -- it divides by every record considered for the guide, spanning the site or not, which is how Amrani et al. describe their own calculation. The gap between it and `whole_corpus` is the size of the coverage correction, i.e. how much apparent 'variation' is really just records that stop short of the site.

| tier | guides_with_>=N_records | perfect_at_183 | below_0.99 | below_0.95 | below_0.70 | still_100pct | pearson_r |
|---|---|---|---|---|---|---|---|
| whole_corpus | 4777 | 833 | 119 | 0 | 0 | 183 | 0.9919 |
| subgenomic_only | 4769 | 833 | 168 | 15 | 0 | 586 | 0.9285 |
| partial_genomes_only | 4777 | 833 | 145 | 0 | 0 | 268 | 0.9919 |
| no_coverage_correction | 4777 | 833 | 833 | 823 | 0 | 0 | 0.985 |

Per gene, on the sub-genomic tier alone -- this is where the denominators genuinely differ, and a pooled figure hides it:

| gene | guides | median_subgenomic_records | max_subgenomic_records | perfect_at_183_with_usable_n | of_those_below_0.95_subgenomic | min_subgenomic_conservation_of_perfect |
|---|---|---|---|---|---|---|
| RL2 | 713 | 27 | 38 | 91 | 0 | 0.9667 |
| UL19 | 959 | 34 | 36 | 247 | 13 | 0.8889 |
| UL29 | 783 | 43 | 56 | 124 | 1 | 0.9302 |
| UL30 | 794 | 326 | 538 | 104 | 0 | 0.9601 |
| UL5 | 483 | 107 | 311 | 94 | 0 | 0.9533 |
| UL52 | 673 | 83 | 101 | 121 | 1 | 0.9167 |
| UL54 | 372 | 42 | 44 | 52 | 0 | 0.9762 |

Note the denominators: the median guide has 43 sub-genomic records behind it (range of per-gene medians 27 to 326). A guide scored 100% over ~30 records is not the same claim as one scored 100% over 500; the exact binomial 95% lower bound on 30/30 is about 0.88. Only UL30 (and to a lesser extent UL5) has enough independent records for a strong per-gene statement.

Per-guide gene-level numbers: `results/robustness_gene_level.tsv`. Per-record corpus manifest: `data/gene_corpus_manifest.tsv`.

## D. Head-to-head with the published Excision BioTherapeutics guides

The four SaCas9 guides of Amrani et al. 2024 (Table 1) scored with OUR method against OUR corpora -- exact 26-mer (20 nt spacer + NNGRRT PAM) presence on either strand. This is not a quotation of their reported numbers; it is the same measurement applied to their sequences.

| guide | spacer | pam | present_in_complete_genomes_n183 | conservation_complete_genomes | gene_level_records_covering | gene_level_present | conservation_gene_level |
|---|---|---|---|---|---|---|---|
| Amrani2024_ICP0g1 | GTACCCGACGGCCCCCGCGT | CGGAGT | 180 | 0.9836 | 445 | 417 | 0.9371 |
| Amrani2024_ICP0g2 | CTCAGGCCGCGAACCAAGAA | CAGAGT | 155 | 0.847 | 445 | 359 | 0.8067 |
| Amrani2024_ICP27g1 | AATCCTAGACACGCACCGCC | AGGAGT | 179 | 0.9781 | 444 | 432 | 0.973 |
| Amrani2024_ICP27g2 | TCGCCAGCGTCATTAGCGGG | GGGGGT | 178 | 0.9727 | 444 | 416 | 0.9369 |

## What this means for guide selection

**Redundancy is real but small.** 183 deposited genomes behave like 132 at 99.90% identity. Collapsing them changes the perfectly-conserved count by +24 (2.9% of the baseline). The direction is the informative part: removing near-duplicates usually *raises* the count, because a duplicate adds no new sequence but can add new assembly noise. A count that rises on de-duplication is not evidence of a more conserved virus, only of a less redundant denominator.

**A large part of 'not conserved' is 'not sequenced'.** Of 874191 (guide, genome) pairs, 2246 could not be resolved into present/absent by exact matching plus anchor mapping. Treating those as UNKNOWN instead of MISMATCH moves the perfect count from 833 to 931. Both numbers are defensible; what is not defensible is quoting one without the other. The strict number is a lower bound and is the one to publish; the N-tolerant number is the right one for ranking guides against each other, because it does not punish a guide for happening to sit under someone else's assembly gap.

**The widened denominator is smaller than it looks, and this is the most important caveat in this report.** 4552 sub-genomic HSV-1 records were downloaded, but only 1021 of them cover any of the seven target genes; the rest are thymidine-kinase, glycoprotein-G and glycoprotein-B typing fragments from unrelated surveys. The remaining 393 covering records are near-full-length `partial genome` assemblies, i.e. the same kind of evidence as the 183 complete genomes, excluded from stage 1 only by a title convention. So for six of the seven genes this is not really a test against a larger, messier, independent corpus: it is mostly a test against 393 more whole genomes. UL30 is the exception -- it has a genuine clinical-amplicon corpus, because UL30 is the gene sequenced for aciclovir-resistance genotyping.

**The headline number.** Of the 833 guides that are 100% conserved across 183 complete genomes and have at least 10 covering records in the gene-level corpus, 0 (0.0%) fall below 95% when the denominator is widened, and 0 fall below the 70% bar used by Amrani et al. 183 remain at exactly 100%. Correlation between the two conservation measures is r = 0.9919 (Spearman rho = 0.8345).

**The same number on sub-genomic records only** -- the tier that is actually independent of the complete-genome set. 4769 guides clear the 10-record bar there (against 4777 on the whole corpus), of which 833 are perfect at n=183: 15 (1.8%) drop below 95% and 0 below 70%. Where the two tiers disagree, this is the one to believe, and the shrunken guide count is itself the finding: for most of our genes GenBank simply does not hold an independent per-gene corpus to test against.

**The coverage correction is the single biggest effect in this whole analysis, larger than redundancy and ambiguity combined.** Without it -- dividing by every record considered rather than by the records that actually span the site -- 823 of the 833 perfect guides fall below 95% and the mean conservation drops to 0.8248. Almost all of that is an artefact: a 600 bp amplicon does not contain a guide 40 kb away, and a partial genome missing its repeats does not contain an RL2 guide. Any per-gene conservation number computed from a heterogeneous record set without a coverage check -- including a >70%-style threshold applied to a ViPR-scale corpus -- is measuring record length as much as sequence variation. This cuts both ways: it means the competing published numbers are probably pessimistic, not that ours are better.

*(The verdict below is taken from the sub-genomic tier.)*

**The position survives but needs qualification.** Most perfectly conserved guides stay high at gene-level n, but a non-trivial minority does not, and those are exactly the ones that would fail in the clinic. The defensible framing is a two-denominator one: report conservation over complete genomes AND over the gene-level corpus, and select only guides that are high in both.

**Practical rule for guide selection.** Rank on the intersection: perfectly conserved across complete genomes, still >= 95% over the gene-level corpus with a real denominator behind it (not 3 records), no poly-T, GC in band. Guides whose gene-level denominator is small should be labelled 'insufficient evidence', not silently promoted -- a guide with 2/2 records is not more conserved than one with 380/400.

## Parameters

| parameter | value |
|---|---|
| minhash_k | 21 |
| minhash_sketch | 5000 |
| dedup_identity | 0.999 |
| anchor_k | 25 |
| anchor_step | 4 |
| anchor_max_gap | 3000 |
| anchor_max_drift | 60 |
| anchor_max_occurrences | 4 |
| max_bracket_span | 2000 |
| min_records_per_guide | 10 |
| gene_corpus_min_length | 100 |
| gene_corpus_max_length | 145000 |
| partial_genomes_included | True |
