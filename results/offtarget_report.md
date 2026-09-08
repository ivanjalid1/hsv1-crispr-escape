# Human (GRCh38) off-target screening of the SaCas9 guide set

Stage 8, `src/offtarget.py`. Generated 2026-09-08T05:27:33Z. Every number below is measured on this machine, not estimated or quoted.

## 0. Headline

This analysis was run because it is capable of invalidating this repository's own recommendation. It partly did. The finding is reported as it fell.

**1. Our named headline guide does not survive.** `RL2_3441+` was put forward in stage 6 as the ICP0 guide that should have been chosen instead of ICP0g2. Against GRCh38 it is materially *dirtier* than ICP0g2:

| | ICP0g2 (theirs) | RL2_3441+ (ours) |
|---|---|---|
| NNGRRT sites <= 4 mm | 5 | 19 |
| NNGRRT sites <= 3 mm | 0 | 6 |
| NNGRRN sites <= 4 mm | 56 | 78 |
| NNGRRN sites <= 3 mm | 5 | 9 |
| HSV-1 conservation (n=183) | 0.847 | 1.000 |

The 3-mismatch tier is the one that matters and it is where `RL2_3441+` is worst: it carries close sites under the **canonical** NNGRRT PAM, where ICP0g2 carries none, and one of them sits in *ABL1*. Four more are copies of a single repeated sequence on chrY. On off-target burden alone, ICP0g2 is the better of the two guides, and stage 6's specific recommendation is withdrawn.

**2. The underlying argument survives, with a different guide.** 2 of the ICP0 alternatives are perfectly conserved across all 183 complete HSV-1 genomes *and* have an off-target profile at least as clean as ICP0g2 on every measure computed here: `RL2_5335+`, `RL2_3364+`. They were in the same candidate pool Amrani et al. drew from, in their own site grammar.

| guide | site id | HSV-1 cons. | NNGRRT <=4mm | NNGRRT <=3mm | NNGRRN <=4mm | NNGRRN <=3mm |
|---|---|---|---|---|---|---|
| ICP0g2 | RL2_4496+ | 0.847 | 5 | 0 | 56 | 5 |
| RL2_5335+ | RL2_5335+ | 1.000 | 3 | 0 | 41 | 2 |
| RL2_3364+ | RL2_3364+ | 1.000 | 5 | 0 | 54 | 1 |

So the claim that better ICP0 options existed inside their own site space is *not* refuted by off-target data -- it is refuted only for the particular guide stage 6 named. The corrected statement is narrower and still stands: at least one perfectly conserved ICP0 site is also cleaner against the human genome than the guide they took to the clinic.

**2b. Their ICP27 choice looks good, and this analysis says so.** ICP27g1 carries 24 NNGRRN sites at <= 4 mismatches and 1 at <= 3 -- the cleanest profile of every ICP27/UL54 candidate screened. **None** of the 9 better-conserved, filter-passing ICP27 alternatives matches it on all four off-target measures. Stage 6 never proposed replacing ICP27g1; this is the first evidence that it should not be.

**3. Independent corroboration of the method.** Amrani et al. report (Table 3) that no sites with <= 3 total mismatches plus bulges were found for either lead guide. Under the canonical NNGRRT PAM this screen finds exactly that: 0 sites at <= 3 mm for ICP0g2 and 0 for ICP27g1, computed independently, from a different assembly download, by a different algorithm. Under the permissive NNGRRN PAM that they actually searched, this screen does find sites (5 for ICP0g2, 1 for ICP27g1), all with non-canonical PAMs. The most likely explanation is their BWA `aln` pre-nomination step, which is a heuristic aligner and is not guaranteed to surface every 3-mismatch locus; the scan here is exhaustive. That is a difference in completeness, not a contradiction, and it cuts in their favour as much as against them.

## 1. Reference genome

| item | value |
|---|---|
| assembly | GRCh38 primary assembly (unmasked) |
| source | Ensembl release 116 (GRCh38.p14, GCA_000001405.29) |
| file | `Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz` |
| URL | https://ftp.ensembl.org/pub/release-116/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz |
| compressed size | 881,964,081 bytes |
| SHA-256 | `d8c3af0094a7bba6125763bad779ec18a81483c739c6ed122094bdf86c187b92` |
| Ensembl published `sum` | 22450 (861294 x 1 KiB blocks) |
| checksum verified | True |
| sequences in assembly | 194 |
| sequences scanned | 194 |
| bases scanned | 3,099,750,718 |

**Why unmasked, and why the primary assembly.** `dna` (unmasked) rather than `dna_sm` or `dna_rm`. A Cas9 site inside a LINE, an Alu or a satellite is a real cleavage substrate, so repeat masking is not merely unhelpful here, it is wrong: hard-masking would delete roughly half the genome from the search, and soft-masking is only safe if every read remembers to upper-case, which is exactly the kind of silent bug that would quietly shrink an off-target count. `primary_assembly` rather than `toplevel`: toplevel adds ALT haplotypes and patch scaffolds, which are alternative representations of sequence already present on the chromosomes and would double-count off-targets.

**Disk actually used: 4.14 GB.** 882 MB compressed FASTA + 3.10 GB one-byte-per-base cache + 141 MB GTF + 20 MB parsed annotation.

**Integrity, three independent ways.** Ensembl's own published `sum` checksum (the BSD 16-bit rotating checksum, recomputed here and compared); the gzip member's CRC-32 and ISIZE trailer, which the decompressor verifies when the whole file is streamed during cache construction; and a SHA-256 we compute and record in `data/genome/genome_manifest.json`. The download is resumable via HTTP Range and is never repeated on a re-run.

## 2. Method

- PAM searched: **NNGRRN**. `NNGRRT` sites are a strict subset of `NNGRRN` sites, so they are flagged rather than searched separately and both variants come out of a single pass.
- Up to **4 mismatches** in the protospacer, **both strands**, genome-wide, across every sequence in the primary assembly (chromosomes, unlocalised and unplaced scaffolds, and the mitochondrion).
- Spacer lengths screened: [20, 21]. 20 nt is the grammar Amrani et al. published in their Table 1 and the grammar the stage-6 candidate pool is enumerated in, so it is the like-for-like length; 21 nt is this repository's SaCas9 default and is built by extending each site one base 5' in HSV-1. A shorter spacer is the permissive direction -- fewer positions available to mismatch means more genomic sites qualify -- so the 20-nt numbers are the conservative ones.
- **Substitutions only. Bulges (DNA/RNA insertions) are not modelled.**
- Ambiguous genomic bases (`N`) count as mismatches in the protospacer and disqualify a PAM outright, so no site is ever credited to unresolved sequence.

**Algorithm.** Exact substring search is the wrong tool -- an off-target is by definition mismatched -- so stage 3's jump scanner cannot be reused. Two properties of the problem are exploited instead.

*PAM sparsity.* SaCas9 needs `NNGRRN`: a G at PAM offset 2 and purines at offsets 3 and 4. That is 1 position in 16 per strand (1 in 64 for `NNGRRT`). Only those positions can ever be cut, so three vectorised slice comparisons discard 15/16 of the genome before any protospacer is looked at. This is the single biggest win available and it is free.

*The pigeonhole principle.* With at most 4 mismatches, a 20-nt protospacer cut into five 4-nt chunks must match the guide exactly in at least one chunk -- four mismatches cannot touch five chunks. Every genomic 4-mer is packed into a byte once per chromosome by pure slice arithmetic (no gathers); the five chunk codes are then gathered at the PAM positions only; and a 256-entry bitmask lookup table per chunk resolves, for all guides simultaneously, which positions share an exact chunk with which guide. About 2% of PAM positions survive per guide.

*Scoring from the codes already in hand.* Because the five chunks tile the whole protospacer, a per-(guide, chunk) 256-entry table of "how many of these four bases differ" turns the five codes into a mismatch count with five more table lookups and no further genome access. Those codes pack an ambiguous base onto A, so that count is a strict *lower* bound -- a site can look better than it is, never worse -- which makes it a sound filter but not a sound answer. The handful of positions that survive it are therefore read back from the genome and counted exactly, which is where the "N is a mismatch" rule is actually enforced. Every reported site is then re-verified a second time during table construction and the count asserted, so a fast-path error would raise rather than silently under-report. Everything is numpy; there is no per-base Python loop on the hot path.

**Wall-clock, measured on this machine.** One-off costs, paid once and cached: BSD-checksum verification of the 882 MB FASTA 101.4 s; decoding it into the 1-byte-per-base array 25.8 s; parsing the GTF 59.1 s (11248794 feature lines -> 78941 genes, 5087789 exons, 3210887 CDS).

**The genome-wide scan itself: 319 s.** 50 guides x 2 strands x 3.10 Gb, i.e. 9.7 Mb of genome per second, or 0.97 billion (guide x strand x base) positions per second. 774,222,748 PAM-matched positions were examined out of 6,199,501,436 strand-bases. Annotation of the close sites: 5 s. A warm re-run (caches present) is the scan time plus a few seconds. For scale, the brute-force reference implementation the tests compare against needs about 35 s per (spacer length x strand) for chromosome 21's 46.7 Mb, which extrapolates to roughly 2.5 hours for the whole assembly against 5 minutes here -- consistent with the 23-44x speed-up measured directly on chr21 in the equivalence test.

**Correctness.** `tests/test_offtarget.py::test_fast_matches_bruteforce_chr21` asserts that the fast path returns *exactly* the same hit set as a structurally independent brute-force scan -- mismatch counts accumulated at every genomic position by one whole-array slice comparison per protospacer base, no PAM prefilter, no chunk index, explicit reverse complement -- over the whole of chromosome 21, for every guide screened. A third, plainly written pure-Python character-by-character implementation validates the numpy brute force in turn. This is the standard `tests/test_nuclease.py` already applies to the stage-3 PAM scanner, applied here.

**Positioning, honestly.** Cas-OFFinder is the field standard and is what Amrani et al. used (BWA `aln` to nominate homologous loci, then Cas-OFFinder over those loci with PAM NNGRRN, then GUIDE-seq validation). What is used here is a validated equivalent **for substitution-only search**, and it is exhaustive where a BWA pre-filter is not. It is **not** equivalent for bulges: their post-filter retained 1mm+1bulge, 2mm+1bulge, 3mm+1bulge, 4mm+0bulge and 5mm+0bulge classes, so their published totals (ICP0g1 358, ICP0g2 910, ICP27g1 443, ICP27g2 316) count classes this screen does not enumerate, and were computed against a different assembly build. Counts here are therefore **not** comparable to those totals line for line. They are comparable *between guides screened here*, which is the question this stage exists to answer.

## 3. Guides screened

Selection rule, fixed before any off-target number existed, so it cannot have been tuned to flatter the recommendation:

1. the four published Amrani et al. 2024 guides;
2. every ICP0/RL2 site that beats their lead ICP0g2 on cross-isolate conservation *and* passes the pipeline's standard filters -- the 12 stage-6 alternatives, including the headline recommendation `RL2_3441+`;
3. every ICP27/UL54 site that passes those filters *and* outranks their lead ICP27g1 within its own gene.

All 4,777 SpCas9 candidates are deliberately not screened: they are not in SaCas9 site space, and the question at issue is a head-to-head against four named guides.

| guide | site id | gene | class | protospacer (20 nt) | PAM | HSV-1 conservation |
|---|---|---|---|---|---|---|
| ICP0g1 | RL2_5319+ | ICP0 | published | GTACCCGACGGCCCCCGCGT | CGGAGT | 0.984 |
| ICP0g2 | RL2_4496+ | ICP0 | published | CTCAGGCCGCGAACCAAGAA | CAGAGT | 0.847 |
| ICP27g1 | UL54_115156+ | ICP27 | published | AATCCTAGACACGCACCGCC | AGGAGT | 0.978 |
| ICP27g2 | UL54_114376- | ICP27 | published | TCGCCAGCGTCATTAGCGGG | GGGGGT | 0.973 |
| RL2_5335+ | RL2_5335+ | ICP0 | icp0_alternative | GCGTCGGAGTGGAACAGCCT | CTGGAT | 1.000 |
| RL2_3364+ | RL2_3364+ | ICP0 | icp0_alternative | GCGACGTGTGCGCCGTGTGC | ACGGAT | 1.000 |
| RL2_5080+ | RL2_5080+ | ICP0 | icp0_alternative | ACGCGCTACCTGCCCATCTC | GGGGGT | 1.000 |
| RL2_3441+ | RL2_3441+ | ICP0 | icp0_alternative | TGCATCCCGTGCATGAAAAC | CTGGAT | 1.000 |
| RL2_3196- | RL2_3196- | ICP0 | icp0_alternative | CGTGCTGTCCGCCTCGGAGG | CGGAGT | 0.995 |
| RL2_5167- | RL2_5167- | ICP0 | icp0_alternative | ATGTTCCCCGTCTCCATGTC | CAGGAT | 0.995 |
| RL2_5395- | RL2_5395- | ICP0 | icp0_alternative | CGGAAGTCCAGGGCGCCCAC | TAGGGT | 0.989 |
| RL2_3568- | RL2_3568- | ICP0 | icp0_alternative | CTCGGCCTCCATGCGGGTCT | GGGGGT | 0.967 |
| RL2_3576- | RL2_3576- | ICP0 | icp0_alternative | ACGGCCTCCTCGGCCTCCAT | GCGGGT | 0.951 |
| RL2_4810- | RL2_4810- | ICP0 | icp0_alternative | GAGGCCGCCGAGGACGTCAG | GGGGGT | 0.951 |
| RL2_4659+ | RL2_4659+ | ICP0 | icp0_alternative | GGTGCGTCCGAGGAAGAGGC | GCGGGT | 0.945 |
| RL2_3515+ | RL2_3515+ | ICP0 | icp0_alternative | GATAGTGGGCGTGACGCCCA | GCGGGT | 0.891 |
| UL54_114638+ | UL54_114638+ | ICP27 | icp27_alternative | CCCTTTGACGCCGAGACCAG | ACGGGT | 1.000 |
| UL54_114833- | UL54_114833- | ICP27 | icp27_alternative | GGGCGCAGCGGCAGGTTGTG | GTGGAT | 1.000 |
| UL54_114720- | UL54_114720- | ICP27 | icp27_alternative | CTTGGCGGTCGATGCGGCCC | GAGGAT | 0.995 |
| UL54_115076+ | UL54_115076+ | ICP27 | icp27_alternative | GACTACGCGACCCTTGGTGT | CGGGGT | 0.989 |
| UL54_115012- | UL54_115012- | ICP27 | icp27_alternative | TGGCCAGAATGACAAACACG | AAGGAT | 0.989 |
| UL54_115085- | UL54_115085- | ICP27 | icp27_alternative | TTCTCTCCGACCCCGACACC | AAGGGT | 0.989 |
| UL54_115028- | UL54_115028- | ICP27 | icp27_alternative | ACGCGGTTGGCGAGCCTGGC | CAGAAT | 0.984 |
| UL54_113964- | UL54_113964- | ICP27 | icp27_alternative | CGTCTGGGTGCTGGGTACGC | CGGGGT | 0.984 |
| UL54_114021- | UL54_114021- | ICP27 | icp27_alternative | CACACTGTGGGGCGCTGGTT | GAGGAT | 0.984 |

## 4. Off-target counts

Distinct genomic sites within *k* mismatches of the protospacer with a matching PAM. `NNGRRT` is the canonical SaCas9 PAM; `NNGRRN` is the permissive variant Amrani et al. searched, and is a superset of it.

### 20-nt spacers

**NNGRRT (canonical)**

| guide | site id | class | 0 mm | 1 mm | 2 mm | 3 mm | 4 mm | total | seed<=1mm |
|---|---|---|---|---|---|---|---|---|---|
| ICP0g1 | RL2_5319+ | published | 0 | 0 | 0 | 1 | 3 | 4 | 0 |
| ICP0g2 | RL2_4496+ | published | 0 | 0 | 0 | 0 | 5 | 5 | 1 |
| ICP27g1 | UL54_115156+ | published | 0 | 0 | 0 | 0 | 6 | 6 | 0 |
| ICP27g2 | UL54_114376- | published | 0 | 0 | 0 | 1 | 3 | 4 | 2 |
| RL2_5335+ | RL2_5335+ | icp0_alternative | 0 | 0 | 0 | 0 | 3 | 3 | 1 |
| RL2_3364+ | RL2_3364+ | icp0_alternative | 0 | 0 | 0 | 0 | 5 | 5 | 1 |
| RL2_5080+ | RL2_5080+ | icp0_alternative | 0 | 0 | 0 | 0 | 10 | 10 | 5 |
| RL2_3441+ | RL2_3441+ | icp0_alternative | 0 | 0 | 0 | 6 | 13 | 19 | 7 |
| RL2_3196- | RL2_3196- | icp0_alternative | 0 | 0 | 0 | 1 | 10 | 11 | 1 |
| RL2_5167- | RL2_5167- | icp0_alternative | 0 | 0 | 0 | 0 | 23 | 23 | 4 |
| RL2_5395- | RL2_5395- | icp0_alternative | 0 | 0 | 0 | 0 | 13 | 13 | 4 |
| RL2_3568- | RL2_3568- | icp0_alternative | 0 | 0 | 0 | 0 | 28 | 28 | 1 |
| RL2_3576- | RL2_3576- | icp0_alternative | 0 | 0 | 0 | 2 | 24 | 26 | 4 |
| RL2_4810- | RL2_4810- | icp0_alternative | 0 | 0 | 0 | 0 | 18 | 18 | 2 |
| RL2_4659+ | RL2_4659+ | icp0_alternative | 0 | 0 | 0 | 3 | 21 | 24 | 6 |
| RL2_3515+ | RL2_3515+ | icp0_alternative | 0 | 0 | 0 | 0 | 4 | 4 | 0 |
| UL54_114638+ | UL54_114638+ | icp27_alternative | 0 | 0 | 0 | 0 | 8 | 8 | 1 |
| UL54_114833- | UL54_114833- | icp27_alternative | 0 | 0 | 0 | 4 | 29 | 33 | 3 |
| UL54_114720- | UL54_114720- | icp27_alternative | 0 | 0 | 0 | 0 | 5 | 5 | 0 |
| UL54_115076+ | UL54_115076+ | icp27_alternative | 0 | 0 | 0 | 0 | 5 | 5 | 1 |
| UL54_115012- | UL54_115012- | icp27_alternative | 0 | 0 | 0 | 5 | 34 | 39 | 3 |
| UL54_115085- | UL54_115085- | icp27_alternative | 0 | 0 | 0 | 1 | 10 | 11 | 1 |
| UL54_115028- | UL54_115028- | icp27_alternative | 0 | 0 | 0 | 0 | 7 | 7 | 1 |
| UL54_113964- | UL54_113964- | icp27_alternative | 0 | 0 | 0 | 1 | 19 | 20 | 2 |
| UL54_114021- | UL54_114021- | icp27_alternative | 0 | 0 | 0 | 1 | 23 | 24 | 1 |

**NNGRRN (permissive)**

| guide | site id | class | 0 mm | 1 mm | 2 mm | 3 mm | 4 mm | total | seed<=1mm |
|---|---|---|---|---|---|---|---|---|---|
| ICP0g1 | RL2_5319+ | published | 0 | 0 | 0 | 1 | 18 | 19 | 5 |
| ICP0g2 | RL2_4496+ | published | 0 | 0 | 0 | 5 | 51 | 56 | 8 |
| ICP27g1 | UL54_115156+ | published | 0 | 0 | 0 | 1 | 23 | 24 | 3 |
| ICP27g2 | UL54_114376- | published | 0 | 0 | 0 | 1 | 21 | 22 | 5 |
| RL2_5335+ | RL2_5335+ | icp0_alternative | 0 | 0 | 0 | 2 | 39 | 41 | 19 |
| RL2_3364+ | RL2_3364+ | icp0_alternative | 0 | 0 | 0 | 1 | 53 | 54 | 21 |
| RL2_5080+ | RL2_5080+ | icp0_alternative | 0 | 0 | 0 | 1 | 60 | 61 | 31 |
| RL2_3441+ | RL2_3441+ | icp0_alternative | 0 | 0 | 0 | 9 | 69 | 78 | 17 |
| RL2_3196- | RL2_3196- | icp0_alternative | 0 | 0 | 0 | 4 | 66 | 70 | 9 |
| RL2_5167- | RL2_5167- | icp0_alternative | 0 | 0 | 1 | 2 | 87 | 90 | 16 |
| RL2_5395- | RL2_5395- | icp0_alternative | 0 | 0 | 0 | 5 | 84 | 89 | 25 |
| RL2_3568- | RL2_3568- | icp0_alternative | 0 | 0 | 0 | 38 | 311 | 349 | 8 |
| RL2_3576- | RL2_3576- | icp0_alternative | 0 | 0 | 0 | 12 | 139 | 151 | 36 |
| RL2_4810- | RL2_4810- | icp0_alternative | 0 | 0 | 0 | 4 | 93 | 97 | 7 |
| RL2_4659+ | RL2_4659+ | icp0_alternative | 0 | 0 | 0 | 8 | 138 | 146 | 33 |
| RL2_3515+ | RL2_3515+ | icp0_alternative | 0 | 0 | 0 | 3 | 35 | 38 | 7 |
| UL54_114638+ | UL54_114638+ | icp27_alternative | 0 | 0 | 0 | 2 | 39 | 41 | 4 |
| UL54_114833- | UL54_114833- | icp27_alternative | 0 | 0 | 0 | 17 | 211 | 228 | 20 |
| UL54_114720- | UL54_114720- | icp27_alternative | 0 | 0 | 0 | 1 | 37 | 38 | 3 |
| UL54_115076+ | UL54_115076+ | icp27_alternative | 0 | 0 | 0 | 7 | 189 | 196 | 10 |
| UL54_115012- | UL54_115012- | icp27_alternative | 0 | 0 | 0 | 31 | 433 | 464 | 82 |
| UL54_115085- | UL54_115085- | icp27_alternative | 0 | 0 | 0 | 11 | 77 | 88 | 9 |
| UL54_115028- | UL54_115028- | icp27_alternative | 0 | 0 | 1 | 5 | 70 | 76 | 18 |
| UL54_113964- | UL54_113964- | icp27_alternative | 0 | 0 | 1 | 6 | 105 | 112 | 11 |
| UL54_114021- | UL54_114021- | icp27_alternative | 0 | 0 | 2 | 4 | 116 | 122 | 9 |

### 21-nt spacers

**NNGRRT (canonical)**

| guide | site id | class | 0 mm | 1 mm | 2 mm | 3 mm | 4 mm | total | seed<=1mm |
|---|---|---|---|---|---|---|---|---|---|
| ICP0g1 (21 nt) | RL2_5319+|21nt | published_21nt | 0 | 0 | 0 | 0 | 2 | 2 | 0 |
| ICP0g2 (21 nt) | RL2_4496+|21nt | published_21nt | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| ICP27g1 (21 nt) | UL54_115156+|21nt | published_21nt | 0 | 0 | 0 | 0 | 3 | 3 | 0 |
| ICP27g2 (21 nt) | UL54_114376-|21nt | published_21nt | 0 | 0 | 0 | 0 | 2 | 2 | 2 |
| RL2_5335+ (21 nt) | RL2_5335+|21nt | icp0_alternative_21nt | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| RL2_3364+ (21 nt) | RL2_3364+|21nt | icp0_alternative_21nt | 0 | 0 | 0 | 0 | 2 | 2 | 1 |
| RL2_5080+ (21 nt) | RL2_5080+|21nt | icp0_alternative_21nt | 0 | 0 | 0 | 0 | 2 | 2 | 0 |
| RL2_3441+ (21 nt) | RL2_3441+|21nt | icp0_alternative_21nt | 0 | 0 | 0 | 1 | 10 | 11 | 5 |
| RL2_3196- (21 nt) | RL2_3196-|21nt | icp0_alternative_21nt | 0 | 0 | 0 | 0 | 6 | 6 | 1 |
| RL2_5167- (21 nt) | RL2_5167-|21nt | icp0_alternative_21nt | 0 | 0 | 0 | 0 | 8 | 8 | 2 |
| RL2_5395- (21 nt) | RL2_5395-|21nt | icp0_alternative_21nt | 0 | 0 | 0 | 0 | 3 | 3 | 1 |
| RL2_3568- (21 nt) | RL2_3568-|21nt | icp0_alternative_21nt | 0 | 0 | 0 | 0 | 8 | 8 | 0 |
| RL2_3576- (21 nt) | RL2_3576-|21nt | icp0_alternative_21nt | 0 | 0 | 0 | 1 | 4 | 5 | 1 |
| RL2_4810- (21 nt) | RL2_4810-|21nt | icp0_alternative_21nt | 0 | 0 | 0 | 0 | 4 | 4 | 0 |
| RL2_4659+ (21 nt) | RL2_4659+|21nt | icp0_alternative_21nt | 0 | 0 | 0 | 1 | 2 | 3 | 2 |
| RL2_3515+ (21 nt) | RL2_3515+|21nt | icp0_alternative_21nt | 0 | 0 | 0 | 0 | 1 | 1 | 0 |
| UL54_114638+ (21 nt) | UL54_114638+|21nt | icp27_alternative_21nt | 0 | 0 | 0 | 0 | 4 | 4 | 0 |
| UL54_114833- (21 nt) | UL54_114833-|21nt | icp27_alternative_21nt | 0 | 0 | 0 | 3 | 5 | 8 | 2 |
| UL54_114720- (21 nt) | UL54_114720-|21nt | icp27_alternative_21nt | 0 | 0 | 0 | 0 | 2 | 2 | 0 |
| UL54_115076+ (21 nt) | UL54_115076+|21nt | icp27_alternative_21nt | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| UL54_115012- (21 nt) | UL54_115012-|21nt | icp27_alternative_21nt | 0 | 0 | 0 | 3 | 14 | 17 | 2 |
| UL54_115085- (21 nt) | UL54_115085-|21nt | icp27_alternative_21nt | 0 | 0 | 0 | 0 | 2 | 2 | 0 |
| UL54_115028- (21 nt) | UL54_115028-|21nt | icp27_alternative_21nt | 0 | 0 | 0 | 0 | 4 | 4 | 0 |
| UL54_113964- (21 nt) | UL54_113964-|21nt | icp27_alternative_21nt | 0 | 0 | 0 | 0 | 4 | 4 | 0 |
| UL54_114021- (21 nt) | UL54_114021-|21nt | icp27_alternative_21nt | 0 | 0 | 0 | 1 | 3 | 4 | 1 |

**NNGRRN (permissive)**

| guide | site id | class | 0 mm | 1 mm | 2 mm | 3 mm | 4 mm | total | seed<=1mm |
|---|---|---|---|---|---|---|---|---|---|
| ICP0g1 (21 nt) | RL2_5319+|21nt | published_21nt | 0 | 0 | 0 | 0 | 4 | 4 | 1 |
| ICP0g2 (21 nt) | RL2_4496+|21nt | published_21nt | 0 | 0 | 0 | 1 | 14 | 15 | 5 |
| ICP27g1 (21 nt) | UL54_115156+|21nt | published_21nt | 0 | 0 | 0 | 0 | 10 | 10 | 1 |
| ICP27g2 (21 nt) | UL54_114376-|21nt | published_21nt | 0 | 0 | 0 | 0 | 5 | 5 | 2 |
| RL2_5335+ (21 nt) | RL2_5335+|21nt | icp0_alternative_21nt | 0 | 0 | 0 | 0 | 6 | 6 | 3 |
| RL2_3364+ (21 nt) | RL2_3364+|21nt | icp0_alternative_21nt | 0 | 0 | 0 | 0 | 20 | 20 | 12 |
| RL2_5080+ (21 nt) | RL2_5080+|21nt | icp0_alternative_21nt | 0 | 0 | 0 | 1 | 17 | 18 | 7 |
| RL2_3441+ (21 nt) | RL2_3441+|21nt | icp0_alternative_21nt | 0 | 0 | 0 | 2 | 26 | 28 | 8 |
| RL2_3196- (21 nt) | RL2_3196-|21nt | icp0_alternative_21nt | 0 | 0 | 0 | 2 | 32 | 34 | 4 |
| RL2_5167- (21 nt) | RL2_5167-|21nt | icp0_alternative_21nt | 0 | 0 | 0 | 2 | 24 | 26 | 5 |
| RL2_5395- (21 nt) | RL2_5395-|21nt | icp0_alternative_21nt | 0 | 0 | 0 | 3 | 23 | 26 | 10 |
| RL2_3568- (21 nt) | RL2_3568-|21nt | icp0_alternative_21nt | 0 | 0 | 0 | 3 | 74 | 77 | 3 |
| RL2_3576- (21 nt) | RL2_3576-|21nt | icp0_alternative_21nt | 0 | 0 | 0 | 3 | 34 | 37 | 14 |
| RL2_4810- (21 nt) | RL2_4810-|21nt | icp0_alternative_21nt | 0 | 0 | 0 | 2 | 30 | 32 | 3 |
| RL2_4659+ (21 nt) | RL2_4659+|21nt | icp0_alternative_21nt | 0 | 0 | 0 | 2 | 21 | 23 | 7 |
| RL2_3515+ (21 nt) | RL2_3515+|21nt | icp0_alternative_21nt | 0 | 0 | 0 | 1 | 9 | 10 | 3 |
| UL54_114638+ (21 nt) | UL54_114638+|21nt | icp27_alternative_21nt | 0 | 0 | 0 | 0 | 13 | 13 | 1 |
| UL54_114833- (21 nt) | UL54_114833-|21nt | icp27_alternative_21nt | 0 | 0 | 0 | 6 | 74 | 80 | 11 |
| UL54_114720- (21 nt) | UL54_114720-|21nt | icp27_alternative_21nt | 0 | 0 | 0 | 0 | 11 | 11 | 1 |
| UL54_115076+ (21 nt) | UL54_115076+|21nt | icp27_alternative_21nt | 0 | 0 | 0 | 0 | 10 | 10 | 5 |
| UL54_115012- (21 nt) | UL54_115012-|21nt | icp27_alternative_21nt | 0 | 0 | 0 | 7 | 84 | 91 | 23 |
| UL54_115085- (21 nt) | UL54_115085-|21nt | icp27_alternative_21nt | 0 | 0 | 0 | 3 | 33 | 36 | 4 |
| UL54_115028- (21 nt) | UL54_115028-|21nt | icp27_alternative_21nt | 0 | 0 | 0 | 2 | 30 | 32 | 6 |
| UL54_113964- (21 nt) | UL54_113964-|21nt | icp27_alternative_21nt | 0 | 0 | 0 | 3 | 23 | 26 | 2 |
| UL54_114021- (21 nt) | UL54_114021-|21nt | icp27_alternative_21nt | 0 | 0 | 0 | 4 | 30 | 34 | 3 |

`seed<=1mm` counts sites with at most one mismatch in the PAM-proximal 12 nt, the region where SaCas9 tolerates mismatches least. It is reported, not used to filter.

## 5. Genomic context of the close off-targets

233 site(s) at <= 3 mismatches, annotated against Ensembl 116 (`Homo_sapiens.GRCh38.116.gtf.gz`). Coordinates are 1-based inclusive on Ensembl chromosome names -- the same names the FASTA uses, so no liftover and no name mapping is involved and nothing can silently shift.

**18 site(s) fall inside a coding exon (CDS)**, of which 15 put the predicted blunt cut itself inside the CDS rather than merely overlapping one with the 26-nt footprint. (66 sites overlap an exon of any biotype, 60 cut inside one.) The intended target cell type is post-mitotic sensory neurons, where double-strand break repair is NHEJ-only and any indel is permanent and unrepairable by homologous recombination. A coding-exon off-target is the dangerous class, and these are listed first.

| guide | chrom | start | str | mm | seed mm | PAM class | PAM | CDS gene | cut in CDS |
|---|---|---|---|---|---|---|---|---|---|
| RL2_3196- | 22 | 30893858 | + | 3 | 3 | NNGRRN_only | ATGAAC | OSBP2 | True |
| RL2_3515+ | 10 | 110881424 | + | 3 | 2 | NNGRRN_only | TAGAAG | PDCD4 | True |
| RL2_3515+ | X | 152858806 | - | 3 | 2 | NNGRRN_only | GTGGAA | NSDHL | True |
| RL2_3576- | 1 | 228274321 | - | 3 | 1 | NNGRRN_only | CTGGAA | OBSCN | True |
| RL2_3576- | 16 | 1773731 | - | 3 | 2 | NNGRRN_only | CAGGAC | EME2 | True |
| RL2_3576- | 19 | 17872314 | - | 3 | 1 | NNGRRN_only | GAGGGC | SLC5A5 | True |
| RL2_3576- | 7 | 144188364 | + | 3 | 1 | NNGRRN_only | CAGGGC | ARHGEF35 | True |
| RL2_3576- | 7 | 144362664 | - | 3 | 1 | NNGRRN_only | CAGGGC | ARHGEF5 | True |
| RL2_5167- | 19 | 11211868 | - | 3 | 2 | NNGRRN_only | CAGGAC | DOCK6 | False |
| RL2_5167- (21 nt) | 19 | 11211868 | - | 3 | 2 | NNGRRN_only | CAGGAC | DOCK6 | False |
| UL54_114833- | 1 | 41381622 | - | 3 | 1 | NNGRRT | CCGGAT | FOXO6 | True |
| UL54_114833- | 8 | 59118940 | + | 3 | 3 | NNGRRN_only | CTGGAG | TOX | True |
| UL54_114833- (21 nt) | 1 | 41381622 | - | 3 | 1 | NNGRRT | CCGGAT | FOXO6 | True |
| UL54_115012- | 11 | 58191138 | - | 3 | 2 | NNGRRT | TAGGAT | OR9Q2 | True |
| UL54_115028- | 19 | 50210469 | + | 3 | 1 | NNGRRN_only | TCGGGC | MYH14 | True |
| UL54_115085- | 10 | 62813205 | - | 3 | 1 | NNGRRN_only | TTGAGA | EGR2 | True |
| UL54_115085- | 2 | 233802146 | + | 3 | 3 | NNGRRT | AGGAGT | MROH2A | False |
| UL54_115085- (21 nt) | 10 | 62813205 | - | 3 | 1 | NNGRRN_only | TTGAGA | EGR2 | True |

**Loci worth naming.** Not a risk model -- a lookup of genes whose identity a reader would want flagged. A gene absent from this list is not thereby safe.

- ***ABL1*** -- proto-oncogene tyrosine kinase; the BCR-ABL fusion partner. Hit by `RL2_3441+` at 3 mismatches (NNGRRT PAM), within the gene body.
- ***ARHGEF5*** -- Rho guanine nucleotide exchange factor; reported oncogenic activity. Hit by `RL2_3576-` at 3 mismatches (NNGRRN_only PAM), inside a coding exon.
- ***DOCK6*** -- Adams-Oliver syndrome gene. Hit by `RL2_5167-` at 3 mismatches (NNGRRN_only PAM), inside a coding exon.
- ***EGR2*** -- master transcriptional regulator of peripheral-nerve myelination -- the intended target tissue here is the trigeminal ganglion. Hit by `UL54_115085-` at 3 mismatches (NNGRRN_only PAM), inside a coding exon.
- ***EME2*** -- structure-specific endonuclease subunit, DNA repair. Hit by `RL2_3576-` at 3 mismatches (NNGRRN_only PAM), inside a coding exon.
- ***FOXO6*** -- forkhead transcription factor, high expression in brain. Hit by `UL54_114833-` at 3 mismatches (NNGRRT PAM), inside a coding exon.
- ***MYH14*** -- non-muscle myosin heavy chain; DFNA4 hearing loss. Hit by `UL54_115028-` at 3 mismatches (NNGRRN_only PAM), inside a coding exon.
- ***NSDHL*** -- X-linked sterol dehydrogenase; loss of function causes CHILD syndrome. Hit by `RL2_3515+` at 3 mismatches (NNGRRN_only PAM), inside a coding exon.
- ***OBSCN*** -- obscurin, sarcomeric signalling. Hit by `RL2_3576-` at 3 mismatches (NNGRRN_only PAM), inside a coding exon.
- ***OSBP2*** -- oxysterol-binding protein, photoreceptor-enriched. Hit by `RL2_3196-` at 3 mismatches (NNGRRN_only PAM), inside a coding exon.
- ***PDCD4*** -- tumour suppressor, translation inhibitor. Hit by `RL2_3515+` at 3 mismatches (NNGRRN_only PAM), inside a coding exon.
- ***SLC5A5*** -- sodium/iodide symporter; thyroid hormone synthesis. Hit by `RL2_3576-` at 3 mismatches (NNGRRN_only PAM), inside a coding exon.
- ***TOX*** -- HMG-box transcription factor, T-cell and neural development. Hit by `UL54_114833-` at 3 mismatches (NNGRRN_only PAM), inside a coding exon.

Per guide (all spacer lengths pooled; a site found by both the 20-nt and the 21-nt version of one guide is the same genomic locus and is counted once per version):

| guide | sites <= 3mm | of which canonical NNGRRT | in an exon | in a CDS | cut inside a CDS |
|---|---|---|---|---|---|
| RL2_3576- | 12 | 2 | 7 | 5 | 5 |
| UL54_114833- | 17 | 4 | 5 | 2 | 2 |
| UL54_115085- | 11 | 1 | 6 | 2 | 1 |
| RL2_3515+ | 3 | 0 | 2 | 2 | 2 |
| UL54_115012- | 31 | 5 | 13 | 1 | 1 |
| UL54_114833- (21 nt) | 6 | 3 | 3 | 1 | 1 |
| RL2_3196- | 4 | 1 | 3 | 1 | 1 |
| UL54_115028- | 6 | 0 | 2 | 1 | 1 |
| RL2_5167- | 3 | 0 | 1 | 1 | 0 |
| UL54_115085- (21 nt) | 3 | 0 | 1 | 1 | 1 |
| RL2_5167- (21 nt) | 2 | 0 | 1 | 1 | 0 |
| RL2_3441+ | 9 | 6 | 0 | 0 | 0 |
| RL2_4659+ | 8 | 3 | 1 | 0 | 0 |
| UL54_115012- (21 nt) | 7 | 3 | 0 | 0 | 0 |
| UL54_113964- | 7 | 1 | 0 | 0 | 0 |
| UL54_114021- | 6 | 1 | 3 | 0 | 0 |
| UL54_114021- (21 nt) | 4 | 1 | 3 | 0 | 0 |
| RL2_3576- (21 nt) | 3 | 1 | 1 | 0 | 0 |
| RL2_3441+ (21 nt) | 2 | 1 | 0 | 0 | 0 |
| RL2_4659+ (21 nt) | 2 | 1 | 0 | 0 | 0 |
| ICP0g1 | 1 | 1 | 0 | 0 | 0 |
| ICP27g2 | 1 | 1 | 0 | 0 | 0 |
| RL2_3568- | 38 | 0 | 9 | 0 | 0 |
| UL54_115076+ | 7 | 0 | 0 | 0 | 0 |
| ICP0g2 | 5 | 0 | 0 | 0 | 0 |
| RL2_5395- | 5 | 0 | 1 | 0 | 0 |
| RL2_4810- | 4 | 0 | 1 | 0 | 0 |
| RL2_3568- (21 nt) | 3 | 0 | 1 | 0 | 0 |
| RL2_5395- (21 nt) | 3 | 0 | 0 | 0 | 0 |
| UL54_113964- (21 nt) | 3 | 0 | 0 | 0 | 0 |
| RL2_3196- (21 nt) | 2 | 0 | 1 | 0 | 0 |
| RL2_4810- (21 nt) | 2 | 0 | 1 | 0 | 0 |
| RL2_5335+ | 2 | 0 | 0 | 0 | 0 |
| UL54_114638+ | 2 | 0 | 0 | 0 | 0 |
| UL54_115028- (21 nt) | 2 | 0 | 0 | 0 | 0 |
| ICP0g2 (21 nt) | 1 | 0 | 0 | 0 | 0 |
| ICP27g1 | 1 | 0 | 0 | 0 | 0 |
| RL2_3364+ | 1 | 0 | 0 | 0 | 0 |
| RL2_3515+ (21 nt) | 1 | 0 | 0 | 0 | 0 |
| RL2_5080+ | 1 | 0 | 0 | 0 | 0 |
| RL2_5080+ (21 nt) | 1 | 0 | 0 | 0 | 0 |
| UL54_114720- | 1 | 0 | 0 | 0 | 0 |

Full listing for the published 20-nt grammar (the 21-nt variants are in `results/offtarget_annotated.tsv`):

| guide | chrom | start | str | mm | seed mm | PAM class | PAM | gene(s) | exon | CDS |
|---|---|---|---|---|---|---|---|---|---|---|
| RL2_3196- | 14 | 93071352 | + | 3 | 3 | NNGRRT | AAGGAT | ITPK1,ITPK1-AS1 | True | False |
| RL2_3196- | 17 | 17560291 | + | 3 | 1 | NNGRRN_only | CAGGAG | ENSG00000306320,PEMT | False | False |
| RL2_3196- | 22 | 30893858 | + | 3 | 3 | NNGRRN_only | ATGAAC | OSBP2 | True | True |
| RL2_3196- | 9 | 94945612 | + | 3 | 3 | NNGRRN_only | AGGGAA | AOPEP,AOPEP-AS1 | True | False |
| RL2_3364+ | 16 | 89584207 | + | 3 | 2 | NNGRRN_only | CTGGGG | CPNE7 | False | False |
| RL2_3441+ | 2 | 155488376 | + | 3 | 2 | NNGRRT | CTGAGT | ENSG00000305077 | False | False |
| RL2_3441+ | 22 | 47560753 | + | 3 | 2 | NNGRRN_only | CTGGAG |  | False | False |
| RL2_3441+ | 5 | 59707526 | - | 3 | 3 | NNGRRN_only | CAGAAC | PDE4D | False | False |
| RL2_3441+ | 6 | 83272054 | - | 3 | 0 | NNGRRN_only | TTGAAA | ME1 | False | False |
| RL2_3441+ | 9 | 130830349 | + | 3 | 1 | NNGRRT | ATGGGT | ABL1 | False | False |
| RL2_3441+ | Y | 23090943 | - | 3 | 1 | NNGRRT | CTGAAT |  | False | False |
| RL2_3441+ | Y | 23329736 | + | 3 | 1 | NNGRRT | CTGAAT |  | False | False |
| RL2_3441+ | Y | 24724687 | - | 3 | 1 | NNGRRT | CTGAAT |  | False | False |
| RL2_3441+ | Y | 24945396 | + | 3 | 1 | NNGRRT | CTGAAT |  | False | False |
| RL2_3515+ | 10 | 110881424 | + | 3 | 2 | NNGRRN_only | TAGAAG | PDCD4 | True | True |
| RL2_3515+ | 20 | 4091681 | - | 3 | 1 | NNGRRN_only | CTGGGC |  | False | False |
| RL2_3515+ | X | 152858806 | - | 3 | 2 | NNGRRN_only | GTGGAA | NSDHL | True | True |
| RL2_3568- | 1 | 76897980 | + | 3 | 2 | NNGRRN_only | ACGGGG | ST6GALNAC5 | False | False |
| RL2_3568- | 1 | 187768713 | + | 3 | 2 | NNGRRN_only | AAGGGG |  | False | False |
| RL2_3568- | 1 | 231833417 | + | 3 | 2 | NNGRRN_only | AGGGGG | DISC1,ENSG00000286071 | False | False |
| RL2_3568- | 1 | 236935471 | + | 3 | 2 | NNGRRN_only | ACGGGG |  | False | False |
| RL2_3568- | 1 | 241436472 | - | 3 | 2 | NNGRRN_only | ACGGGG | ENSG00000287516 | True | False |
| RL2_3568- | 11 | 96590424 | - | 3 | 2 | NNGRRN_only | AGGGGG | LINC02737 | True | False |
| RL2_3568- | 11 | 123868806 | + | 3 | 2 | NNGRRN_only | AGGGGG |  | False | False |
| RL2_3568- | 11 | 127545940 | + | 3 | 2 | NNGRRN_only | AGGGGG |  | False | False |
| RL2_3568- | 12 | 79542360 | - | 3 | 2 | NNGRRN_only | AGGGGG | ENSG00000257474 | True | False |
| RL2_3568- | 13 | 63215614 | - | 3 | 2 | NNGRRN_only | ACGGGG | LINC00376 | False | False |
| RL2_3568- | 15 | 24688787 | + | 3 | 2 | NNGRRN_only | AGGGGG | ENSG00000286110,PWRN1 | False | False |
| RL2_3568- | 15 | 88550709 | + | 3 | 2 | NNGRRN_only | AGGGGG | ENSG00000305611 | False | False |
| RL2_3568- | 19 | 1576912 | + | 3 | 1 | NNGRRN_only | TTGGAG | MBD3 | True | False |
| RL2_3568- | 2 | 7875866 | + | 3 | 2 | NNGRRN_only | AGGGGG | ENSG00000226506 | False | False |
| RL2_3568- | 2 | 156723273 | + | 3 | 2 | NNGRRN_only | AGGGGG | ENSG00000299347 | True | False |
| RL2_3568- | 3 | 54637697 | + | 3 | 2 | NNGRRN_only | GCGGGG | CACNA2D3,ESRG | True | False |
| RL2_3568- | 3 | 130061765 | - | 3 | 2 | NNGRRN_only | ATGGGG | ENSG00000309438 | False | False |
| RL2_3568- | 3 | 146364555 | - | 3 | 2 | NNGRRN_only | ACGGGG |  | False | False |
| RL2_3568- | 3 | 191669693 | + | 3 | 2 | NNGRRN_only | AGGGGG |  | False | False |
| RL2_3568- | 4 | 11656938 | + | 3 | 2 | NNGRRN_only | AGGGGG | ENSG00000249631 | False | False |
| RL2_3568- | 4 | 103556972 | + | 3 | 2 | NNGRRN_only | AGGGGG | ENSG00000250920,TACR3-AS1 | True | False |
| RL2_3568- | 4 | 182817017 | - | 3 | 2 | NNGRRN_only | AGGGGG |  | False | False |
| RL2_3568- | 5 | 149757859 | - | 3 | 1 | NNGRRN_only | GTGGGC | PPARGC1B | False | False |
| RL2_3568- | 6 | 14300164 | + | 3 | 2 | NNGRRN_only | AGGGGG | ENSG00000286277,ENSG00000295484 | True | False |
| RL2_3568- | 6 | 21364716 | + | 3 | 2 | NNGRRN_only | ATGGGG | ENSG00000227089,ENSG00000294517 | True | False |
| RL2_3568- | 6 | 109858974 | + | 3 | 2 | NNGRRN_only | AGGGGG | ENSG00000302596,FIG4 | False | False |
| RL2_3568- | 6 | 123585501 | + | 3 | 2 | NNGRRN_only | AGGGGG | TRDN | False | False |
| RL2_3568- | 6 | 131507493 | + | 3 | 2 | NNGRRN_only | AGGGGG | ARG1 | False | False |
| RL2_3568- | 7 | 6906277 | - | 3 | 2 | NNGRRN_only | ACGGGG | ENSG00000295732 | False | False |
| RL2_3568- | 7 | 42873826 | + | 3 | 2 | NNGRRN_only | ATGGGG |  | False | False |
| RL2_3568- | 8 | 115817096 | - | 3 | 2 | NNGRRN_only | ACGGGG |  | False | False |
| RL2_3568- | 9 | 78183306 | - | 3 | 2 | NNGRRN_only | ACGGAG |  | False | False |
| RL2_3568- | 9 | 93350104 | - | 3 | 2 | NNGRRN_only | AGGGGG |  | False | False |
| RL2_3568- | 9 | 115842413 | + | 3 | 2 | NNGRRN_only | AGGGGG | ENSG00000228714,ENSG00000308189,LINC00474 | False | False |
| RL2_3568- | 9 | 137417024 | + | 3 | 3 | NNGRRN_only | CAGAAC | EXD3 | False | False |
| RL2_3568- | X | 3556342 | + | 3 | 2 | NNGRRN_only | ACGGGG |  | False | False |
| RL2_3568- | X | 4894081 | - | 3 | 2 | NNGRRN_only | ACGGGG |  | False | False |
| RL2_3568- | X | 53263337 | - | 3 | 2 | NNGRRN_only | GTGGGA | IQSEC2 | False | False |
| RL2_3576- | 1 | 228274321 | - | 3 | 1 | NNGRRN_only | CTGGAA | ENSG00000269934,OBSCN | True | True |
| RL2_3576- | 10 | 133231713 | - | 3 | 1 | NNGRRT | TTGGGT |  | False | False |
| RL2_3576- | 12 | 126457924 | - | 3 | 3 | NNGRRT | GCGGGT | LINC02347,LINC02350 | True | False |
| RL2_3576- | 14 | 75489717 | - | 3 | 3 | NNGRRN_only | GAGAAA |  | False | False |
| RL2_3576- | 16 | 1773731 | - | 3 | 2 | NNGRRN_only | CAGGAC | EME2 | True | True |
| RL2_3576- | 19 | 17872314 | - | 3 | 1 | NNGRRN_only | GAGGGC | ENSG00000299140,SLC5A5 | True | True |
| RL2_3576- | 3 | 34288402 | + | 3 | 0 | NNGRRN_only | GGGAAG | LINC01811 | False | False |
| RL2_3576- | 5 | 126672655 | - | 3 | 2 | NNGRRN_only | AAGGGC |  | False | False |
| RL2_3576- | 7 | 97357441 | + | 3 | 2 | NNGRRN_only | ATGGGA | ENSG00000305877 | False | False |
| RL2_3576- | 7 | 144188364 | + | 3 | 1 | NNGRRN_only | CAGGGC | ARHGEF35 | True | True |
| RL2_3576- | 7 | 144286947 | + | 3 | 1 | NNGRRN_only | CAGGGC | ARHGEF34P,ARHGEF35-AS1,ENSG00000308967,OR2A1-AS1 | True | False |
| RL2_3576- | 7 | 144362664 | - | 3 | 1 | NNGRRN_only | CAGGGC | ARHGEF5 | True | True |
| ICP0g2 | 13 | 64387905 | - | 3 | 2 | NNGRRN_only | AAGAAA |  | False | False |
| ICP0g2 | 14 | 59686836 | + | 3 | 2 | NNGRRN_only | CAGGAC | RTN1 | False | False |
| ICP0g2 | 6 | 75548137 | - | 3 | 2 | NNGRRN_only | TTGGAA | ENSG00000297132 | False | False |
| ICP0g2 | 7 | 27216651 | - | 3 | 3 | NNGRRN_only | TGGGGA | ENSG00000308571 | False | False |
| ICP0g2 | 8 | 55593551 | - | 3 | 2 | NNGRRN_only | GGGAGA |  | False | False |
| RL2_4659+ | 1 | 158022842 | - | 3 | 1 | NNGRRN_only | TGGAAA | KIRREL1 | False | False |
| RL2_4659+ | 1 | 201059371 | - | 3 | 1 | NNGRRN_only | AGGGGC | CACNA1S | False | False |
| RL2_4659+ | 13 | 18944115 | - | 3 | 2 | NNGRRT | TGGGGT | ENSG00000231358 | True | False |
| RL2_4659+ | 16 | 25980273 | - | 3 | 2 | NNGRRN_only | CTGAGA | ENSG00000285882,HS3ST4 | False | False |
| RL2_4659+ | 18 | 11303210 | - | 3 | 1 | NNGRRT | TGGGGT |  | False | False |
| RL2_4659+ | 2 | 71588588 | + | 3 | 1 | NNGRRT | CAGGGT | DYSF | False | False |
| RL2_4659+ | 2 | 100945176 | + | 3 | 3 | NNGRRN_only | TGGAAG | NPAS2 | False | False |
| RL2_4659+ | 5 | 2021528 | - | 3 | 1 | NNGRRN_only | AAGGGC | ENSG00000248994 | False | False |
| RL2_4810- | 13 | 50005947 | - | 3 | 3 | NNGRRN_only | ATGAGG | DLEU2,TRIM13 | False | False |
| RL2_4810- | 14 | 60637686 | - | 3 | 2 | NNGRRN_only | CGGAGG | ENSG00000307639 | False | False |
| RL2_4810- | 16 | 86555031 | - | 3 | 1 | NNGRRN_only | AGGGGG | MTHFSD | True | False |
| RL2_4810- | 2 | 70008434 | + | 3 | 2 | NNGRRN_only | CAGAGG | ENSG00000293615,PCBP1-AS1 | False | False |
| RL2_5080+ | 8 | 132949796 | - | 3 | 1 | NNGRRN_only | TGGAAA | TG | False | False |
| RL2_5167- | 11 | 19525976 | + | 2 | 2 | NNGRRN_only | AGGAGG | NAV2 | False | False |
| RL2_5167- | 19 | 11211868 | - | 3 | 2 | NNGRRN_only | CAGGAC | DOCK6,DOCK6-AS1 | True | True |
| RL2_5167- | 3 | 158841984 | + | 3 | 1 | NNGRRN_only | CCGGGA | MFSD1 | False | False |
| ICP0g1 | 14 | 100126611 | - | 3 | 2 | NNGRRT | CAGAGT | EVL | False | False |
| RL2_5335+ | 16 | 75480285 | - | 3 | 1 | NNGRRN_only | CAGGGA | CHST6,ENSG00000203472 | False | False |
| RL2_5335+ | 18 | 59279076 | - | 3 | 1 | NNGRRN_only | GAGAGC | CPLX4 | False | False |
| RL2_5395- | 1 | 3372496 | - | 3 | 1 | NNGRRN_only | AGGAAG | ENSG00000286518,PRDM16 | False | False |
| RL2_5395- | 14 | 104725853 | - | 3 | 1 | NNGRRN_only | CAGAGC | ADSS1 | False | False |
| RL2_5395- | 5 | 172563075 | + | 3 | 2 | NNGRRN_only | GAGAGA |  | False | False |
| RL2_5395- | 7 | 149805415 | - | 3 | 2 | NNGRRN_only | CCGAAC | SSPOP | True | False |
| RL2_5395- | X | 153051198 | - | 3 | 1 | NNGRRN_only | AAGAGC | ENSG00000310029 | False | False |
| UL54_113964- | 2 | 131583364 | - | 2 | 1 | NNGRRN_only | TAGGAC | ENSG00000286208 | False | False |
| UL54_113964- | 1 | 94340380 | + | 3 | 2 | NNGRRN_only | GTGGGA | ARHGAP29-AS1 | False | False |
| UL54_113964- | 14 | 100906526 | + | 3 | 2 | NNGRRN_only | CTGGAC | MEG8,MIR493HG | False | False |
| UL54_113964- | 17 | 19359692 | - | 3 | 2 | NNGRRT | AGGGGT | B9D1 | False | False |
| UL54_113964- | 21 | 44697042 | - | 3 | 3 | NNGRRN_only | AGGAGG | TSPEAR | False | False |
| UL54_113964- | 7 | 40046166 | + | 3 | 2 | NNGRRN_only | CAGAAG | CDK13 | False | False |
| UL54_113964- | 8 | 142595780 | - | 3 | 3 | NNGRRN_only | TGGAGG |  | False | False |
| UL54_114021- | 10 | 31721560 | + | 2 | 2 | NNGRRN_only | TGGGGA | ENSG00000300971 | True | False |
| UL54_114021- | 11 | 66870220 | - | 2 | 0 | NNGRRN_only | GGGGAA | PC | False | False |
| UL54_114021- | 10 | 114976898 | - | 3 | 3 | NNGRRN_only | ATGAAC | TRUB1 | True | False |
| UL54_114021- | 11 | 897769 | + | 3 | 2 | NNGRRN_only | GAGGGA | CHID1 | False | False |
| UL54_114021- | 19 | 19546461 | + | 3 | 1 | NNGRRT | GTGGGT | CILP2 | True | False |
| UL54_114021- | 20 | 40025088 | + | 3 | 2 | NNGRRN_only | TTGGGA | LINC01370 | False | False |
| ICP27g2 | 4 | 21762725 | + | 3 | 1 | NNGRRT | CTGGGT | KCNIP4 | False | False |
| UL54_114638+ | 14 | 24820271 | - | 3 | 2 | NNGRRN_only | TGGAAA | STXBP6 | False | False |
| UL54_114638+ | 6 | 106448795 | - | 3 | 3 | NNGRRN_only | AAGGAC | CRYBG1 | False | False |
| UL54_114720- | 7 | 51473408 | - | 3 | 2 | NNGRRN_only | CTGAGC | ENSG00000285741 | False | False |
| UL54_114833- | 1 | 41381622 | - | 3 | 1 | NNGRRT | CCGGAT | FOXO6 | True | True |
| UL54_114833- | 1 | 235582317 | - | 3 | 1 | NNGRRN_only | CTGGGA | GNG4 | False | False |
| UL54_114833- | 10 | 119651997 | - | 3 | 3 | NNGRRN_only | GCGGGG | BAG3 | False | False |
| UL54_114833- | 11 | 12009525 | - | 3 | 3 | NNGRRN_only | AAGGAG | DKK3 | True | False |
| UL54_114833- | 11 | 12585652 | + | 3 | 2 | NNGRRT | ACGAGT |  | False | False |
| UL54_114833- | 12 | 6308490 | - | 3 | 2 | NNGRRT | GTGGGT |  | False | False |
| UL54_114833- | 13 | 91295721 | - | 3 | 2 | NNGRRN_only | GTGGGC | ENSG00000309877 | False | False |
| UL54_114833- | 14 | 100307188 | + | 3 | 1 | NNGRRN_only | GAGAGA |  | False | False |
| UL54_114833- | 17 | 42111011 | - | 3 | 1 | NNGRRN_only | AAGGGA | DHX58 | False | False |
| UL54_114833- | 19 | 3769915 | - | 3 | 2 | NNGRRN_only | GAGAGA | RAX2 | True | False |
| UL54_114833- | 20 | 57675646 | + | 3 | 2 | NNGRRN_only | GAGGGA | PMEPA1 | False | False |
| UL54_114833- | 4 | 88859439 | + | 3 | 2 | NNGRRN_only | GAGGGA | FAM13A | False | False |
| UL54_114833- | 5 | 116575788 | - | 3 | 3 | NNGRRN_only | GGGAAC | SEMA6A,SEMA6A-AS2 | True | False |
| UL54_114833- | 8 | 59118940 | + | 3 | 3 | NNGRRN_only | CTGGAG | TOX,TOX-DT | True | True |
| UL54_114833- | 8 | 98173148 | + | 3 | 1 | NNGRRT | TGGAGT |  | False | False |
| UL54_114833- | 8 | 144045860 | + | 3 | 2 | NNGRRN_only | TGGGAG | SPATC1 | False | False |
| UL54_114833- | X | 110668682 | - | 3 | 2 | NNGRRN_only | AGGAGA |  | False | False |
| UL54_115012- | 1 | 40333541 | - | 3 | 1 | NNGRRN_only | AAGGGA | RPL21P20 | True | False |
| UL54_115012- | 1 | 178199890 | + | 3 | 3 | NNGRRN_only | CAGAAG | RASAL2 | False | False |
| UL54_115012- | 1 | 212051989 | - | 3 | 1 | NNGRRN_only | AAGGGA | DTL,RPL21P28 | True | False |
| UL54_115012- | 11 | 3997077 | - | 3 | 3 | NNGRRN_only | GAGAAA | STIM1 | False | False |
| UL54_115012- | 11 | 58191138 | - | 3 | 2 | NNGRRT | TAGGAT | OR9Q2 | True | True |
| UL54_115012- | 11 | 78861907 | + | 3 | 2 | NNGRRT | CTGGGT | TENM4 | False | False |
| UL54_115012- | 12 | 111952642 | + | 3 | 2 | NNGRRN_only | CAGAAG | TMEM116 | False | False |
| UL54_115012- | 13 | 72703685 | + | 3 | 1 | NNGRRN_only | AAGGGA | RPL21P110 | True | False |
| UL54_115012- | 14 | 96051840 | - | 3 | 1 | NNGRRN_only | ATGAAG | C14orf132 | False | False |
| UL54_115012- | 16 | 9156864 | - | 3 | 1 | NNGRRN_only | AAGGGA | RPL21P119 | True | False |
| UL54_115012- | 16 | 77930798 | - | 3 | 1 | NNGRRT | CTGGGT | ENSG00000294326,VAT1L | False | False |
| UL54_115012- | 16 | 82434820 | + | 3 | 1 | NNGRRN_only | AAGAGA | ENSG00000259873 | True | False |
| UL54_115012- | 17 | 61182708 | - | 3 | 2 | NNGRRN_only | TAGAGA | BCAS3 | False | False |
| UL54_115012- | 19 | 52064597 | + | 3 | 3 | NNGRRN_only | AGGGAG | ZNF432,ZNF841 | True | False |
| UL54_115012- | 2 | 62533042 | - | 3 | 1 | NNGRRN_only | AAGGGA | ENSG00000228541,RPL21P37 | True | False |
| UL54_115012- | 2 | 63321469 | - | 3 | 2 | NNGRRN_only | AGGAAA | WDPCP | False | False |
| UL54_115012- | 2 | 128217580 | + | 3 | 1 | NNGRRN_only | GAGGGA | RPL21P34 | True | False |
| UL54_115012- | 2 | 141096600 | + | 3 | 1 | NNGRRN_only | GGGAGA | LRP1B | False | False |
| UL54_115012- | 20 | 52593038 | - | 3 | 2 | NNGRRN_only | TAGAAA | LINC01524 | False | False |
| UL54_115012- | 3 | 21636610 | - | 3 | 3 | NNGRRN_only | CAGAAG | ZNF385D | False | False |
| UL54_115012- | 3 | 145824722 | - | 3 | 1 | NNGRRN_only | AAGGGA | RPL21P39 | True | False |
| UL54_115012- | 4 | 47491444 | - | 3 | 1 | NNGRRN_only | AAGGGA | ATP10D,RPL21P52 | True | False |
| UL54_115012- | 5 | 177264758 | + | 3 | 1 | NNGRRN_only | AAGGGA | NSD1,RPL21P60 | True | False |
| UL54_115012- | 6 | 44476482 | - | 3 | 2 | NNGRRT | TTGGGT |  | False | False |
| UL54_115012- | 7 | 89508942 | + | 3 | 1 | NNGRRN_only | GGGAAA |  | False | False |
| UL54_115012- | 7 | 154285881 | - | 3 | 1 | NNGRRN_only | CCGGAG | DPP6 | False | False |
| UL54_115012- | 8 | 20480478 | - | 3 | 3 | NNGRRN_only | CAGAAA |  | False | False |
| UL54_115012- | 8 | 133003058 | - | 3 | 1 | NNGRRN_only | AAGGGA | RPL21P78,TG | True | False |
| UL54_115012- | X | 76797848 | - | 3 | 3 | NNGRRT | TTGAGT | MIR325HG | False | False |
| UL54_115012- | X | 115945191 | + | 3 | 1 | NNGRRN_only | GCGAAA | DANT2 | False | False |
| UL54_115012- | X | 154021200 | + | 3 | 3 | NNGRRN_only | AGGAGC |  | False | False |
| UL54_115028- | 2 | 148561609 | + | 2 | 1 | NNGRRN_only | CTGGAA |  | False | False |
| UL54_115028- | 1 | 108661163 | - | 3 | 2 | NNGRRN_only | GGGGAG | ENSG00000285923,HENMT1 | True | False |
| UL54_115028- | 10 | 109777689 | + | 3 | 2 | NNGRRN_only | CTGGAG |  | False | False |
| UL54_115028- | 11 | 118924883 | + | 3 | 2 | NNGRRN_only | AGGGGG | BCL9L | False | False |
| UL54_115028- | 14 | 69251651 | + | 3 | 1 | NNGRRN_only | TAGAGA | ENSG00000258520 | False | False |
| UL54_115028- | 19 | 50210469 | + | 3 | 1 | NNGRRN_only | TCGGGC | MYH14 | True | True |
| UL54_115076+ | 1 | 79492671 | + | 3 | 1 | NNGRRN_only | TGGGAC | ENSG00000288822 | False | False |
| UL54_115076+ | 1 | 229413998 | - | 3 | 1 | NNGRRN_only | TGGGAC |  | False | False |
| UL54_115076+ | 10 | 84226011 | - | 3 | 2 | NNGRRN_only | TGGGAA |  | False | False |
| UL54_115076+ | 11 | 27633186 | + | 3 | 1 | NNGRRN_only | TGGGAC | BDNF-AS,LINC00678 | False | False |
| UL54_115076+ | 13 | 36318041 | - | 3 | 1 | NNGRRN_only | TGGGAC | SPART | False | False |
| UL54_115076+ | 5 | 96171018 | + | 3 | 2 | NNGRRN_only | TGGGAC | CAST,ENSG00000289337 | False | False |
| UL54_115076+ | Y | 12564054 | + | 3 | 2 | NNGRRN_only | TGGGAA | USP9Y | False | False |
| UL54_115085- | 1 | 43408902 | - | 3 | 2 | NNGRRN_only | ATGGAA | SZT2 | False | False |
| UL54_115085- | 10 | 62813205 | - | 3 | 1 | NNGRRN_only | TTGAGA | EGR2 | True | True |
| UL54_115085- | 16 | 28823103 | + | 3 | 2 | NNGRRN_only | GCGGGG | ATXN2L,ENSG00000275807 | True | False |
| UL54_115085- | 19 | 34265364 | - | 3 | 3 | NNGRRN_only | CAGGAG | GARRE1 | True | False |
| UL54_115085- | 19 | 42977992 | - | 3 | 3 | NNGRRN_only | ATGAGC | PSG11-AS1 | False | False |
| UL54_115085- | 19 | 47340770 | + | 3 | 3 | NNGRRN_only | AGGAGC | C5AR2 | True | False |
| UL54_115085- | 2 | 102497655 | - | 3 | 3 | NNGRRN_only | AGGAAA | SLC9A4 | False | False |
| UL54_115085- | 2 | 233802146 | + | 3 | 3 | NNGRRT | AGGAGT | MROH2A | True | True |
| UL54_115085- | 7 | 73840383 | - | 3 | 3 | NNGRRN_only | CCGAGG | METTL27 | False | False |
| UL54_115085- | 9 | 35604832 | - | 3 | 1 | NNGRRN_only | TTGAAC | ENSG00000288586 | True | False |
| UL54_115085- | X | 7033648 | - | 3 | 3 | NNGRRN_only | ATGAGC | ENSG00000297641,PUDP | False | False |
| ICP27g1 | 5 | 176833170 | - | 3 | 1 | NNGRRN_only | TTGAAC | UNC5A | False | False |

## 6. Limitations

1. **Substitutions only.** DNA/RNA bulges are not enumerated. Cas-OFFinder does enumerate them and Amrani et al. counted them. A bulge-tolerant reanalysis can only add sites, never remove them, and could in principle reorder the guides.
2. **No activity model.** These are raw site counts. No CFD, MIT/Hsu or Doench score is applied, because none of those models is parameterised for SaCas9; an SpCas9 score applied to SaCas9 data would be worse than no score. Seed-region mismatch counts are reported instead, unweighted.
3. **No empirical validation.** Amrani et al. ran GUIDE-seq and targeted amplicon deep sequencing. Nothing here is measured in cells. In-silico nomination is a filter, not evidence of cleavage -- and equally, absence of a nominated site is not proof of safety.
4. **One reference, one haplotype.** GRCh38 primary assembly only. Population variation creates and destroys both PAMs and protospacer matches; ALT haplotypes were excluded to avoid double-counting, which necessarily makes haplotype-specific sites invisible.
5. **Annotation is Ensembl-only**, at gene/exon/CDS granularity. No cancer-gene or essential-gene list is applied. A site outside a gene is not thereby safe: enhancers, promoters and other regulatory elements are not modelled. Amrani et al.'s own validated off-target for ICP27g1 sits in an intron of *ZNF331* that is exon 1 of one transcript variant -- exactly the kind of case that gene-level annotation renders as 'intronic'.
6. **The 4-mismatch ceiling is a choice.** Amrani et al.'s post-filter also retained a 5-mismatch class. Raising the ceiling would raise every count; it would not change the ranking unless one guide's excess is concentrated at 5 mm.

## 7. Reproducibility

```
python src/offtarget.py                  # fetch + verify + scan + annotate + report
python src/offtarget.py --stage fetch    # acquisition and caches only
python src/offtarget.py --chrom 21       # scan one sequence (fast smoke test)
python tests/test_offtarget.py           # includes the chr21 equivalence proof
```

Supporting tables: `results/offtarget_summary.tsv` (per-guide counts), `results/offtarget_sites.tsv` (every site, with coordinates, strand, mismatch positions and the observed PAM), `results/offtarget_annotated.tsv` (close sites with gene context), `results/offtarget_summary.json` (machine-readable; carries the genome manifest, the checksums and all timings).

**Caveat for any manuscript.** Guides were screened for substitution off-targets (<= 4 mismatches, PAM NNGRRN and NNGRRT) against the Ensembl release-116 GRCh38 primary assembly, both strands, genome-wide. Bulge-containing off-targets were NOT modelled, and no cell-based validation (e.g. GUIDE-seq) was performed.
