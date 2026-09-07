# Amrani et al. 2024 — EBT-104 (Excision BioTherapeutics) — methodological notes

**Citation.** Amrani N, Luk K, Singh P, Shipley M, Isik M, Donadoni M, Bellizzi A, Khalili K,
Sariyer IK, Neumann D, Gordon J, Ruan G-X. "CRISPR-Cas9-mediated genome editing delivered by a
single AAV9 vector inhibits HSV-1 reactivation in a latent rabbit keratitis model."
*Mol Ther Methods Clin Dev* 2024;32(3):101303. doi:10.1016/j.omtm.2024.101303.
Epub 14 Aug 2024; collection 12 Sep 2024. PMID 39610766 · PMCID PMC11602521.

**URLs.**
- Full text (open access): https://pmc.ncbi.nlm.nih.gov/articles/PMC11602521/
- Europe PMC: https://europepmc.org/articles/PMC11602521
- Publisher: https://doi.org/10.1016/j.omtm.2024.101303
- XML archived locally: `./amrani2024_full.xml` ; plain text: `./amrani2024_full.txt`

**Affiliations.** (1) Excision BioTherapeutics, Watertown MA; (2) Univ. Wisconsin-Madison,
Ophthalmology & Visual Sciences (Neumann); (3) Temple Univ. Lewis Katz School of Medicine,
Center for Neurovirology and Gene Editing (Khalili, Sariyer).
Corresponding authors: D. Neumann (dneumann3@wisc.edu), G-X. Ruan (eruan@excisionbio.com).

---

## 1. Nuclease — RESOLVED

**SaCas9** (*Staphylococcus aureus* Cas9). NOT SpCas9, NOT CasX.

- PAM used in their off-target search: **`NNGRRN`** (explicit Cas-OFFinder run parameter in Methods).
- All four Table 1 target sites end in an `NNGRRT`-class PAM: ICP0g1 `CGGAGT`, ICP0g2 `CAGAGT`,
  ICP27g1 `AGGAGT`, ICP27g2 `GGGGGT`.
- Spacers are **20 nt** (shorter than the 21–22 nt often used for SaCas9).
- SpCas9 appears only as an incidental reagent: SpCas9 RNP + AAVS1 gRNA used to knock the reporter
  construct into HEK293FT cells. Nothing to do with the therapeutic guides. The "CasX" claim from
  the earlier second-hand check is **wrong** — no CasX anywhere in this paper.
- **Implication for us:** their guide sites live in `NNGRRT`-PAM space. An SpCas9 (`NGG`) guide set
  occupies a *different* site space. Site-for-site benchmarking is not apples-to-apples. Either
  (a) re-run our pipeline under SaCas9 PAM constraints to produce a like-for-like comparison, or
  (b) frame the comparison explicitly as cross-nuclease and compare *methodology + escape
  probability*, not raw site identity.

## 2. Genes targeted — CONFIRMED

**ICP0 (RL2)** and **ICP27 (UL54)**, both immediate-early genes. Confirmed by the literal regex in
Methods: `"ICP0$|ICP0.|RL2$|RL2."` and `"ICP27$|ICP27.|UL54$|UL54."`.

Rationale as stated:
1. Both encode critical IE proteins required for productive infection.
2. **ICP0 overlaps the LAT locus** and ICP27 sits in an adjacent region — argued that chromatin
   accessibility at these loci during latency matters for editing latent genomes.
3. **ICP0 is duplicated** in the HSV-1 genome, so a 2-guide pair introduces **three DSBs**, favouring
   large deletions/excision rather than a single repairable indel.

## 3. Computational guide-selection method — VERBATIM DETAILS (the audit target)

Database and denominator:
- Source: **ViPR / Virus Pathogen Resource**, cited as `https://www.bv-brc.org/view/Virus/10239`,
  **Version Feb 2022**.
- Raw pull: **46,027 HSV-1 strains / 62,554 GenBank IDs**. They downloaded **strain info + CDS FASTA
  sequences** — i.e. **CDS-level, not complete genomes**.
- **Filter step:** the GenBank IDs and CDS FASTAs were filtered using **9,409 GenBank IDs** retrieved
  from the **NCBI Nucleotide** database with the keyword **"human alphaherpesvirus 1."**
- Gene assignment: pure **annotation-string regex** on the labels (see §2). No HMM, no ortholog
  calling, no sequence-level verification described.
- **Conservation metric, exact wording:** "The percentage of guide conservation was calculated by
  dividing the strain counts for each guide RNA by the total strain counts of the ICP0 or ICP27
  in the database."
  → The denominator is the **per-gene strain count among annotated entries surviving the
  9,409-ID filter** — NOT 46,027, NOT 62,554, and NOT a count of complete genomes.
- **Threshold: `>70%`**, strictly greater than 70. VERIFIED from the XML source markup (`&gt;70%`).
  Full criterion: guides "selected based on sequence match to HSV-1 strains in our curated database
  (>70%) and number of in silico nominated off-target sites (Table 2)."
- Consensus construction: all ICP0 / ICP27 CDS FASTAs → **Clustal Omega** MSA (BioPython
  `ClustalOmegaCommandLine`) → consensus via BioPython `AlignIO` with **threshold = 0.5**, i.e. emit
  the most common base if >50%, otherwise emit **`N`**. These consensus CDSs were synthesized
  (GenScript) into the HEK293FT AAVS1 reporter used for all cell-line screening.

### Methodological weak points visible in their own description
- **Per-gene denominators are never reported.** How many strains actually carry an annotated ICP0
  vs ICP27 after the 9,409-ID filter is not stated anywhere. The `>70%` figure is therefore
  uninterpretable as published, and ICP0 and ICP27 almost certainly have *different* denominators.
- **No dereplication.** GenBank HSV-1 entries are heavily redundant (partial CDS records, duplicate
  submissions, lab-passaged strains, single sequencing projects contributing many records). No
  clustering, no de-duplication, no clinical-vs-lab-isolate stratification is described.
- **Annotation-string gene assignment silently drops strains.** Any genome lacking an `ICP0`/`RL2`
  or `ICP27`/`UL54` annotation string disappears from *both* numerator and denominator, which
  inflates apparent conservation.
- **The 0.5 consensus threshold** means the reporter sequence they screened against is a synthetic
  majority-rule construct that may match no real isolate, with ambiguous positions collapsed to `N`.
- **Tiny candidate set.** Only **4 guides total** (2 per gene) were ever evaluated. No candidate-pool
  size, no ranking function, no scoring model beyond ">70% conservation + few nominated off-targets."
- **Geographic/clade coverage never assessed.** No statement about whether the surviving strain set
  is representative of circulating diversity.

### Off-target pipeline (for completeness)
In-house **Nextflow** pipeline: `BWA aln` to find hg38 sites homologous to the spacer →
**Cas-OFFinder** for all guide/target alignments at those sites, PAM `NNGRRN` →
post-filter retaining sites with (1) 1mm+1bulge, (2) 2mm+1bulge, (3) 3mm+1bulge, (4) 4mm+0bulge,
(5) 5mm+0bulge. Empirical validation by **GUIDE-seq** in the reporter line; on/off-target indels by
targeted amplicon deep sequencing (MiSeq v2) analyzed with **CRISPResso2**.

Table 2 (in-silico nominated off-target counts across the mismatch/bulge classes, totals last):
ICP0g1 = 358; ICP0g2 = 910; ICP27g1 = 443; ICP27g2 = 316.
No sites with ≤3 total mismatches+bulges were found for either lead guide (Table 3).

## 4. Guide RNA sequences — PUBLISHED IN TABLE 1

| Guide | Spacer (20 nt, as RNA) | Spacer (DNA) | PAM |
|---|---|---|---|
| ICP0g1  | GUACCCGACGGCCCCCGCGU | GTACCCGACGGCCCCCGCGT | CGGAGT |
| ICP0g2  | CUCAGGCCGCGAACCAAGAA | CTCAGGCCGCGAACCAAGAA | CAGAGT |
| ICP27g1 | AAUCCUAGACACGCACCGCC | AATCCTAGACACGCACCGCC | AGGAGT |
| ICP27g2 | UCGCCAGCGUCAUUAGCGGG | TCGCCAGCGTCATTAGCGGG | GGGGGT |

Table 1 also gives 5-nt lower-case flanks; e.g. ICP0g1 target reads `cccga` + spacer + PAM + `ggaac`;
ICP0g2 `cctgg` + spacer + PAM + `ctgtg`; ICP27g1 `atcga` + spacer + PAM + `gttcg`;
ICP27g2 `ggcaa` + spacer + PAM + `gcttg`.

**Lead clinical pair = ICP0g2 + ICP27g1.** The construct is named `SaCas9-g2g1`
(AAV2-SaCas9-g2g1 in vitro; AAV8(Y733F)- and AAV9-SaCas9-g2g1 in vivo). Six pairwise combinations
of the four guides were screened in Vero cells before g2g1 was selected.

Genomic locations of the four guides are in Figure S2A (supplementary, not in the XML).
**Per-guide conservation percentages are NOT given** in the main text or Table 1 — only the `>70%`
selection criterion. If individual values exist they are in the supplementary (Figures S1A/S1B and
the Table S-series), which is a separate file.

## 5. Escape modelling — THEY MOTIVATE IT BUT DO NOT MODEL, MEASURE, OR QUANTIFY IT

- **Introduction (the entire escape argument):** previous studies indicated that using a single gRNA
  to eliminate HIV-1 led to the creation of escape mutants (refs 15–18); use of a single gRNA against
  HSV-1 may lead to a similar phenomenon after DNA repair; therefore simultaneously targeting two or
  more viral sites, which can induce large DNA deletions, can more effectively prevent viral
  replication and avoid virus survival.
  Refs 15–18 = Wang/Berkhout/Das *Mol Ther* 2016; Yoder & Bundschuh *Sci Rep* 2016;
  Liu/Berkhout/Das *Viruses* 2021; Wang/Liang (CRISPR/Cas9-derived HIV escape).
- **Design rationale:** 2 guides × duplicated ICP0 = **3 DSBs** → excision/large deletion,
  "effectively preventing viral replication and escape following DNA repair."
- **No escape-variant deep sequencing of survivors.** The only sequencing of residual virus is
  on-target amplicon-seq at the **ICP27g1 site only**: mean indel **0.78%** (low dose) and **1.65%**
  (high dose) in trigeminal ganglia at day 40 (Figure S6); discussion rounds this to "~1%".
- **ICP0g2 site could not be amplified or sequenced at all** — GC content ~85% around the target.
  So half the on-target outcome data is simply missing, and the ICP0 arm is unassessed.
- They acknowledge their assay **cannot distinguish escape from clearance**: residual DNA "may have
  undergone repair and possesses gene editing outcomes (InDel or excision)"; and they speculate that
  episomes receiving 2–3 DSBs were degraded/not recircularized and therefore invisible to
  amplicon-seq. The DNA they sequenced "was likely either unedited DNA or DNA that was repaired."
- **No justification for n=2 guides specifically.** "Two or more" is asserted; two is what was built.
  No comparison of 1 vs 2 vs 3 vs 4 guides, no escape-rate estimate, no probabilistic model, no
  escape kinetics, no serial-passage/escape-selection experiment in vitro.
- Timepoint caveat is theirs: indel analysis at 40 days post-AAV / 12 days post-reactivation stimulus;
  "Analyzing at an earlier time point might reveal higher InDel rates."

## 6. Stated limitations (Discussion)

1. Only two delivery routes tested (corneal scarification, IV); alternative/combined routes = future work.
2. Non-viral delivery to the peripheral nervous system unexplored — "much remains unknown."
3. Possible **AAV integration into CRISPR-induced HSV-1 DNA breaks** is unknown; argued low-risk in
   theory but "more in-depth analysis is warranted."
4. Residual off-target editing at one site per guide: **ICP0g2 OT1 = 0.06%** (intergenic, chr8);
   **ICP27g1 OT1 = 0.54%** (intron of *ZNF331*, chr19, but maps to exon 1 of one transcript variant;
   ZNF331 is a transcriptional repressor frequently methylated in oesophageal/gastric/colorectal cancer).
5. Indel analysis at a single late timepoint only.
6. ICP0g2 amplicon not sequenceable (GC ~85%).
7. Unknown what fraction of HSV-1-infected TG neurons are AAV9-transduced; unknown viral-load
   reduction threshold required for clinical benefit. Explicitly flagged as future work.
8. Aspiration stated as future work: identify guides conserved in **both HSV-1 and HSV-2** for a
   "one drug fits all" herpes product.

## 7. Related / companion work

- Program name: **EBT-104** (Excision BioTherapeutics). Company press release announcing this paper:
  14 Aug 2024 (globenewswire; excision.bio press release #45).
- **ASGCT 2024** (Apr 2024): EBT-104 HSV-1 keratitis data presented — precursor to this paper.
- **ASGCT 2025** (13 May 2025, New Orleans): additional EBT-104 HSV data presented alongside their
  HBV program. Conference abstract only; no peer-reviewed follow-up located.
- Khalili / Sariyer (Temple) are co-authors; the same lab's HIV CRISPR work (EBT-101) is the
  methodological lineage for the multiplex-guide/excision strategy.
- **Adjacent competitor line, different group** (Cai / Hong et al.), which Amrani et al. critique:
  - *Nat Biotechnol* 2021;39:567–577 — "Targeting herpes simplex virus with CRISPR-Cas9 cures
    herpetic stromal keratitis in mice."
  - *Mol Ther* 2023;31:3163–3175 — "In vivo CRISPR gene editing in patients with herpetic stromal
    keratitis" (3 patients, lentiviral delivery, concurrent corneal transplantation). Amrani et al.
    criticise both the random-integration risk of lentivirus and the confounding of the concurrent
    transplant.

---

## Bottom line for positioning

Their computational method is thin and fully disclosed: one ViPR snapshot (Feb 2022), an
annotation-string filter down to 9,409 NCBI Nucleotide IDs, a per-gene strain-count ratio with an
**unreported denominator**, a `>70%` cut, and four candidate guides. The escape argument is a
**citation to HIV literature plus a qualitative "more cuts is better"** — no quantitative escape
model, no escape sequencing of survivors, no justification of n=2, and no measurement at all at the
ICP0 site. Both the denominator/conservation-methodology audit and quantitative escape-probability
optimisation are unoccupied ground.

Caveat for benchmarking: they are **SaCas9 (`NNGRRT`)**, so a like-for-like site comparison requires
running our method under SaCas9 PAM constraints.
