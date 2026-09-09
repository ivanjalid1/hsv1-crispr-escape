# What a genome corpus can and cannot certify about CRISPR antiviral escape: resolution floors, joint coverage, and a worked audit of an HSV-1 guide pair

**Ivan Heredia Jalid**

Universidad Católica de Córdoba, Córdoba, Argentina

ORCID: 0009-0003-6702-1295

Correspondence: Ivan Heredia Jalid <ivanjalid@gmail.com>

Preprint. Not peer reviewed.

---

## Abstract

Multiplex CRISPR antivirals rest on an escape argument: cut a viral genome at several conserved sites at once, and no single repair event can restore an uncuttable, viable genome. That argument is almost always supported by per-guide conservation percentages computed against a public sequence corpus. We show that this evidentiary chain has two structural weaknesses, both properties of the method rather than of any design. First, the corpus imposes a resolution floor: a site absent in 0 of 183 complete HSV-1 genomes is still consistent with a population absence frequency of 0.0162 (Clopper–Pearson 95% upper limit). Sampling resolution, not the repair model, sets how low an escape probability can be *demonstrated* — a single site cannot be certified below ~1.6 × 10⁻², a pair below ~2.6 × 10⁻⁴, and the certifiable minimum guide count for a 10⁻⁶ per-genome threshold is four, not the three a point estimate implies. Second, joint intactness is not the product of marginal conservations and cannot be recovered from a per-guide table; it requires a per-genome presence matrix. A published clinical lead pair with marginals 0.847 and 0.978 has joint intactness 0.825 — below either marginal — because its guides fail in different isolates. As a case study we audit the EBT-104 guide set (Amrani et al., 2024): their ICP0g2 site carries a recurrent G→A substitution at spacer position 9 in 25 of 183 genomes, and 96% of the pair's modelled escape traces to isolates that have already lost it. Changing one guide raises joint intactness to 0.978 and improves every off-target measure. Absolute escape probabilities span 3.1 orders of magnitude while guide-set rankings hold in 27 of 37 swept settings — including a spectrum measured in post-mitotic human neurons, which lowers escape without changing any selected set: the model ranks, it does not measure.

---

## 1. Introduction

CRISPR nucleases are being developed as antivirals against persistent DNA viruses, where a latent reservoir is not reachable by conventional antivirals. Herpes simplex virus type 1 (HSV-1) is a leading target: it establishes lifelong latency in sensory ganglia, and no licensed therapy eliminates the latent genome (James et al., 2020; Cohen, 2020). Several groups have taken CRISPR anti-HSV programs into animals and, in one case, into patients (Yin et al., 2021; Wei et al., 2023).

Every such program must answer a question that is quantitative in form but is usually answered qualitatively: how many guides, and which ones, are enough to prevent viral escape? The standard argument runs by analogy to HIV, where single-guide CRISPR pressure has been shown to *accelerate* the emergence of resistant proviruses, because imprecise non-homologous end joining (NHEJ) at the cut site generates exactly the mutations that destroy guide recognition while leaving the genome viable (Wang et al., 2016). The inference drawn is that two or more simultaneous cuts, ideally producing excision of the intervening segment, prevent escape. That inference is plausible. It is also, in the published HSV-1 literature, unquantified.

The evidence offered in support is normally a table of per-guide conservation percentages: the fraction of sequences in some public corpus that contain the target site. The claim "conserved in >95% of strains" is doing two distinct jobs in such a table. It is being read as a statement about coverage — the design will work in almost every isolate — and as a statement about escape — an absent site is one the virus already tolerates losing. Both readings are load-bearing, and both are more fragile than they look.

This paper makes three points, in decreasing order of generality.

**First, the corpus sets a resolution floor.** Conservation is estimated from a finite sample. A site never observed absent has an estimated absence frequency of exactly zero, and a model that takes that estimate literally will report escape probabilities that the data cannot support. The correct statement is an interval, and its upper limit propagates through any multiplex calculation. We compute that propagation and find that it, not the repair term, is what binds. This is a limitation of the field's methodology, and it applies to every design evaluated this way, including our own.

**Second, joint coverage is not the product of marginals, and per-guide tables cannot see it.** Two guides at 0.90 conservation each are jointly intact in between 0.80 and 0.90 of isolates depending entirely on whether their failures coincide. Only a per-genome presence matrix distinguishes the cases. This matters for both the coverage reading and the escape reading, and — as we show — it matters with *opposite signs* for the two, which is easy to get backwards.

**Third, as a worked case study**, we audit a specific published program: the SaCas9-based EBT-104 design of Amrani et al. (2024), which targets ICP0 (RL2) and ICP27 (UL54) and took the guide pair ICP0g2 + ICP27g1 into a rabbit latency model. This audit is a demonstration of the first two points on real clinical-stage guides, not an indictment. Amrani et al. did not err procedurally: their stated >70% conservation criterion was met by all four of their candidates; their lead pair was selected on measured antiviral activity in cells, which is a legitimate and arguably superior axis to any in-silico ranking; and their efficacy screen used a single laboratory strain, which by construction contains all four target sites and therefore *could not* have detected a cross-strain coverage gap. That is a structural blind spot in a standard workflow, not negligence.

We also report two of our own recommendations that we have withdrawn, one of them refuted by our own subsequent screening. That disclosure is not decoration: it is the strongest available evidence for this paper's central methodological claim, which is that single-axis guide selection is fragile regardless of who performs it.

Throughout, we use an alignment-free exact-match definition of conservation. A target site is present in a genome if the full protospacer-plus-PAM sequence occurs verbatim on either strand. This is the right question for the biology — Cas9 cleavage requires an intact PAM and near-perfect protospacer complementarity, so a site differing by one base is a different site — and it removes every alignment parameter and heuristic from the measurement. It is conservative: exact matching can under-state conservation but never over-state it.

---

## 2. Results

### 2.1 The corpus resolution floor

We retrieved every complete HSV-1 genome in GenBank (183 records, 147,898–159,092 bp; Methods) and enumerated all SpCas9 (`NGG`) and SaCas9 (`NNGRRT`) target sites within seven essential or high-value genes. Under SpCas9, 833 of 4,777 candidate 23-mers (17.4%) are present in all 183 genomes, of which 644 also pass standard poly-T, homopolymer and GC filters; under the SaCas9 `NNGRRT` grammar the same corpus yields 67 of 448 (Figure 3).

A guide present in 183 of 183 genomes has a *measured* absence frequency of zero. It does not have a true absence frequency of zero. The exact binomial (Clopper–Pearson) 95% upper limit on 0 events in 183 trials is **0.0162** (Clopper and Pearson, 1934) — essentially the familiar rule of three, 3/*n* = 0.0164 (Hanley and Lippman-Hand, 1983). Every perfectly conserved site in this corpus is therefore consistent with being absent in up to 1.6% of circulating virus.

We propagated that bound through a multiplex escape model (Section 2.7, Methods), which computes, for a viral genome drawn from the empirical isolate population and conditional on exposure to an active nuclease, the probability that it ends up uncuttable at *every* targeted site while still encoding functional essential proteins. Per site this decomposes into a measured term — presence or absence in that genome, from the per-genome matrix — and a repair term *q*, the probability that NHEJ produces an indel which destroys recognition and leaves the protein working. Re-running the model with every site's absence frequency replaced by its 95% upper bound gives the certifiable, as opposed to point-estimate, escape probability (Table 1, Figure 1).

**Table 1.** Point estimate versus certifiable bound, best guide set at each *k*.

| *k* | P(escape), measured absences | P(escape), 95% upper bound on unobserved absence | ratio |
|---|---|---|---|
| 1 | 4.754 × 10⁻³ | 2.091 × 10⁻² | 4.4 |
| 2 | 2.260 × 10⁻⁵ | 4.374 × 10⁻⁴ | 19.4 |
| 3 | 1.075 × 10⁻⁷ | 9.148 × 10⁻⁶ | 85.1 |
| 4 | 5.109 × 10⁻¹⁰ | 1.913 × 10⁻⁷ | 374.5 |
| 5 | 2.429 × 10⁻¹² | 4.001 × 10⁻⁹ | 1,647 |
| — published lead pair | 9.720 × 10⁻⁴ | 1.137 × 10⁻² | 11.7 |

Two consequences follow.

**The floor is set by the corpus, not by repair biology.** With 183 genomes, a single site cannot be certified below ~1.6 × 10⁻² and a pair below ~2.6 × 10⁻⁴ (0.0162², the pure presence-term product) from cross-isolate evidence alone. Any claim that a two-guide design achieves an escape probability below about 3 × 10⁻⁴ is an extrapolation beyond the evidence, no matter how the repair term is parameterised — **including claims made by our own model**. Note that this is a statement about *demonstrability*, not about the truth: the design may well be better than that. The evidence cannot show it.

**The certifiable minimum guide count differs from the point estimate.** For a per-genome escape threshold of 10⁻³, both the point estimate and the certifiable bound require *k* = 2. For 10⁻⁶ the point estimate requires *k* = 3 while the bound requires ***k* = 4** — the divergence visible in Figure 1. The same one-guide gap appears at a 10⁻⁷ threshold (4 versus 5). A design justified on a point estimate at the 10⁻⁶ level is therefore one guide short of what its own evidence base can certify.

The bound is deliberately pessimistic in one respect: it assumes unobserved absences are independent across sites, the direction that maximises escape. It does not shrink with better modelling — only with more genomes — and Limitation 9 shows that widening the corpus is harder than record counts suggest.

### 2.2 Joint coverage is not the product of marginals

Multiplex guide sets are conventionally reported as a list of per-guide conservation values. That representation cannot answer the question a multiplex design actually poses: in what fraction of isolates are *all* the targeted sites simultaneously intact? For two guides, joint intactness is bounded above by the smaller marginal and below by (*a* + *b* − 1), and which point in that interval is realised depends entirely on whether the two guides fail in the same isolates or in different ones. Only a per-genome presence matrix — a binary genome × site array, never collapsed to column means — can distinguish them.

The published EBT-104 lead pair is a clean demonstration. Its two guides have marginal conservations of **0.847** (ICP0g2, present in 155/183 genomes) and **0.978** (ICP27g1, 179/183). If presence were independent, joint intactness would be 0.8285. The measured joint intactness is **0.8251** — 151 of 183 genomes — which is *below either marginal* and sits essentially exactly at the theoretical lower bound of 0.825 (Figure 2). The two guides' failures are almost perfectly non-overlapping: the union of failures is very nearly the sum.

This has two consequences with opposite signs, and getting the sign right requires the joint matrix:

- **Coverage is worse than the marginals suggest.** 32 of 183 isolates lose at least one of the two sites. A reader of the per-guide table would infer that the worse guide, at 0.847, bounds the design; in fact the pair is worse than that guide alone.
- **Escape is *better* than an independence calculation would say.** Escape requires *both* sites to fail in the *same* genome; spreading failures across different isolates makes that coincidence rarer. Assuming independence gives P(escape) = 4.269 × 10⁻³ against the empirical 9.720 × 10⁻⁴ — a relative error of **3.39**, i.e. a 4.4-fold over-statement.

Neither effect is recoverable from marginal percentages. A per-guide conservation table is not a reduced form of the joint matrix; it is a strictly weaker object. We recommend that multiplex designs report joint intactness over a per-genome presence matrix as the primary conservation statistic, with marginals as supporting detail. The underlying computational object is not novel — per-sequence presence matrices optimised jointly rather than marginally are established practice in viral diagnostic design, where maximising the fraction of sequences covered by a designed set is the explicit objective (Metsky et al., 2022), and breadth of coverage across circulating variants has likewise been treated as a design target for CRISPR antiviral guides (Bagchi et al., 2022). The antiviral guide-selection literature has largely not adopted it.

### 2.3 Case study: the EBT-104 ICP0 guide

Amrani et al. (2024) selected four SaCas9 guides — two in ICP0/RL2, two in ICP27/UL54 — against a February 2022 ViPR CDS snapshot filtered to 9,409 NCBI Nucleotide identifiers, using a stated criterion of >70% conservation plus a count of in-silico nominated off-target sites, and took **ICP0g2 + ICP27g1** forward as the clinical pair (construct SaCas9-g2g1).

To make a like-for-like comparison we re-enumerated the ICP0 and ICP27 site space in their own grammar — 20 nt spacer plus `NNGRRT` PAM, as printed in their Table 1 — inside the annotated coding sequence of the same two genes in the reference strain. All four published guides are recovered from that pool, which is asserted in the test suite; if they were not, every rank below would be meaningless. The pool contains **46 ICP0 sites and 35 ICP27 sites**. It is small: in a 68%-GC genome an `NNGRRT` PAM occurs roughly once per 70 bp per strand against once per ~9 bp for SpCas9's `NGG`, so requiring `NNGRRT` costs an order of magnitude in candidate sites (448 versus 4,777 across the seven-gene panel) while the *fraction* perfectly conserved barely moves — 15.0% versus 17.4% (Figure 3). This search-space contraction is the structural reason a SaCas9 program can arrive at a mediocre guide without any procedural error.

**Where their guides fall** (Table 2), scored by our exact-match method against the 183 complete genomes.

**Table 2.** The four published guides, measured by our method.

| guide | site | rank in gene pool | conservation (*n* = 183) | gene-level corpus | GC | local GC (200 bp) | passes filters |
|---|---|---|---|---|---|---|---|
| ICP0g1 | RL2_5319+ | 12 / 46 | 0.9836 | 0.9835 | 0.80 | 0.735 | no (GC, homopolymer 5) |
| **ICP0g2** (lead) | RL2_4496+ | **31 / 46** | **0.8470** | 0.8975 | 0.60 | **0.840** | yes |
| **ICP27g1** (lead) | UL54_115156+ | **12 / 35** | **0.9781** | 0.9908 | 0.60 | 0.590 | yes |
| ICP27g2 | UL54_114376− | 19 / 35 | 0.9727 | 0.9952 | 0.65 | 0.795 | yes |

Thirty of the 46 `NNGRRT` sites in ICP0 are better conserved than ICP0g2. Twelve of those also pass every standard filter, and four are present in all 183 genomes. All 12 lie in both ICP0 repeat copies, so they preserve the three-double-strand-break architecture that motivates the design.

Their >70% criterion is met by all four guides — and by 69 of the 81 sites across both genes. The substantive observation is not that they broke their own rule, but that the rule barely discriminates. A threshold that admits 85% of the candidate space is not performing selection.

**A recurrent, geographically and temporally dispersed variant.** The 28 isolates in which the ICP0g2 site is absent are not assembly noise. We located the homologous region in each by exact 25-mer flanking anchors (Methods) and read out the variant directly. Of the 28: **25 carry a single G→A substitution at spacer position 9** — the *identical* substitution in every case, converting `CTCAGGCC` **G** `CGAACCAAGAA` to `CTCAGGCC` **A** `CGAACCAAGAA` with the `CAGAGT` PAM unchanged; **0 have an ambiguity (N) gap** in the bracketed region; and **3 could not be bracketed** by flanking anchors and are recorded as unresolved rather than absent (Figure 5). There is therefore no evidence that this gap is an artefact of assembly quality; it is a real, recurrent polymorphism.

The carriers are dispersed, and they are not 25 independent observations. Country of origin is recorded for 8 of the 25 (7 USA, 1 China) and collection date for the same 8, spanning **1967 to 2020**; the remaining 17 carry neither in the GenBank source record, which reflects the sparsity of HSV-1 metadata generally (only 108 of 183 genomes report a country, 82 a collection date). Applying the 99.9%-identity clustering described in Methods 5.5, the 25 records collapse to **13 distinct clusters** — several are longitudinal or multi-site samples from one subject, including five records from one subject and three from another in a published serial-sampling series, and two anatomical sites from a single 1982 isolate. The defensible statement is therefore: *the variant is carried by at least 13 independent HSV-1 lineages sampled across at least five decades and at least two continents.* It is not a rare private polymorphism.

Position 9 sits in the PAM-distal half of the spacer, where SaCas9 tolerates mismatches better than in the seed. We make no claim that the substitution abolishes cleavage — that is a wet-lab question. Our claim is definitional and narrower: under the exact-match criterion, and under any criterion treating a one-base spacer mismatch as a different site, 15.3% of sequenced isolates do not carry the intended ICP0g2 target.

**Why a standard workflow could not have seen this.** Amrani et al. screened six pairwise guide combinations in Vero cells against HSV-1 strain 17Syn+ and selected g2g1 on measured antiviral activity. 17Syn+ is the strain from which the reference genome derives and by construction contains all four target sites, so an efficacy screen in a single laboratory strain *cannot* detect a cross-strain coverage gap however well it is executed. Their in-vivo on-target amplicon sequencing could not close the gap either: the ICP0g2 amplicon could not be sequenced at all because of the region's very high GC content — our independent measurement gives 0.840 GC across the 200 bp around the site, corroborating their reported ~85% — so the ICP0 arm has no on-target outcome data. This is the blind spot, and it is a property of the workflow, not of the investigators.

### 2.4 Where the escape probability of the published pair comes from

We scored the lead pair under the escape model. It gives P(escape) = **9.720 × 10⁻⁴** per exposed viral genome under baseline parameters. The decomposition is exact, not fitted:

**Table 3.** Contributions to the published pair's escape probability.

| isolate class | isolates | contribution | share |
|---|---|---|---|
| both sites already absent (escape is free) | 0 | 0 | 0% |
| one site already absent (one repair escape needed) | 32 | 9.348 × 10⁻⁴ | **96.2%** |
| both sites intact (two repair escapes needed) | 151 | 3.721 × 10⁻⁵ | 3.8% |

**96.2% of the pair's modelled escape comes from the 32 isolates that have already lost one of the two sites** — almost all of them the ICP0g2 carriers of Section 2.3. In those isolates the design is effectively a single-guide design, and only one NHEJ escape event is required.

The design conclusion this supports is encouraging rather than damning: **the weakness of the published pair is not that *k* = 2 is too few. It is guide choice.** The same architecture with a better-conserved ICP0 guide clears the same thresholds; adding guides is not the only remedy and is not the cheapest one. More generally, at every site with even moderate cross-isolate conservation, pre-existing absence dominates repair escape by orders of magnitude — the median per-site repair term across the pool is *q* = 4.75 × 10⁻³, because a frameshift in an essential gene is not escape and the in-frame fraction of the assumed indel spectrum is only 0.199. Multiplexing buys protection mainly by covering isolates a single guide misses, not by making NHEJ escape combinatorially unlikely.

**One conclusion runs the other way, and we state it.** Amrani et al. argue that because ICP0 is duplicated in the TRL and IRL repeats, a two-guide pair delivers three double-strand breaks, favouring excision. Our model's `all-copies` repeat regime — in which a guide cutting both ICP0 copies requires an independent viable escape at *each* — reproduces exactly that argument quantitatively: under it, a **single** ICP0 guide reaches P(escape) = 2.260 × 10⁻⁵ at *k* = 1, the value two guides elsewhere require. In that regime their argument is correct. Our default does not assume it, because a redundant gene can absorb a broken copy. Which regime holds is the single highest-value experiment this analysis points to.

### 2.5 A one-guide substitution improves every measured axis

Keeping ICP27g1 exactly as published and replacing only the ICP0 guide with **RL2_5335+** (`GCGTCGGAGTGGAACAGCCT` + `CTGGAT`; reference positions 5,335–5,360; present in both ICP0 repeat copies) improves every measure we can compute (Table 4, Figure 4).

**Table 4.** The published pair versus the single-substitution pair. All values measured here.

| | ICP0g2 + ICP27g1 (published) | **RL2_5335+ + ICP27g1** |
|---|---|---|
| conservation, 183 complete genomes | 0.847 | **1.000** |
| conservation, gene-level corpus | 0.8975 (400 records) | **1.000 (429 records)** |
| conservation, sub-genomic tier only | 0.9524 (21 records) | **1.000 (36 records)** |
| **joint** intactness of the pair | 0.8251 (151/183) | **0.9781 (179/183)** |
| P(escape) per exposed genome | 9.720 × 10⁻⁴ | **2.5145 × 10⁻⁴** (3.87× lower) |
| GRCh38 `NNGRRT` sites ≤ 4 mm / ≤ 3 mm | 5 / 0 | **3 / 0** |
| GRCh38 `NNGRRN` sites ≤ 4 mm / ≤ 3 mm | 56 / 5 | **41 / 2** |
| coding-exon hits ≤ 3 mm (either PAM) | 0 | 0 |
| 21-nt variant, `NNGRRT` / `NNGRRN` ≤ 4 mm | 0 / 15 | **0 / 6** |
| local GC, 200 bp | 0.840 | **0.725** |

Selection used hard gates stated before they were applied — perfect conservation across all 183 genomes; ≥ 0.95 on the gene-level corpus; passes standard filters; present in both ICP0 repeat copies; zero GRCh38 off-targets at ≤ 3 mismatches under the canonical `NNGRRT` PAM; zero coding-exon hits at ≤ 3 mismatches — followed by a single primary discriminator, total `NNGRRT` off-target burden at ≤ 4 mismatches. Three of 14 screened candidates survive the gates, and `RL2_5335+` wins the discriminator outright (3 versus 5 and 10 sites).

Two features of this answer are not tidy, and we report them rather than smoothing them.

**The escape model contributes no information to the choice among the alternatives.** All four perfectly conserved, filter-passing ICP0 candidates have identical presence vectors and identical local codon-tolerance vectors, and therefore identical *q* to ten significant figures (spread 5.9 × 10⁻¹² absolute, 6.3 × 10⁻¹⁰ relative — floating-point summation order, not biology). They are *exactly tied* at P(escape) = 2.5145 × 10⁻⁴. The escape axis separates ICP0g2 from any perfectly conserved candidate by 3.87×, and is dead thereafter. Moreover, the recommended pair sits essentially *at* the resolution floor of Section 2.1 (2.5145 × 10⁻⁴ against a certifiable floor of ~2.6 × 10⁻⁴). Any further improvement in escape probability from guide choice alone would be extrapolation past what this corpus supports — which is itself an argument for choosing on the off-target axis, where the measurements are not at their floor.

**No candidate dominates on every axis.** Testing Pareto dominance mechanically over nine conservation, escape and off-target axes leaves a three-member front. The single measure that blocks domination is `NNGRRN` sites at ≤ 3 mismatches, where `RL2_5335+` has 2 and the other two survivors have 1 each; both of `RL2_5335+`'s are intronic, non-coding, with non-canonical PAMs. The choice therefore rests on one explicitly stated biological priority — that the canonical `NNGRRT` tier, which is what SaCas9 cleaves efficiently and which all four published on-target sites use, outweighs the deliberately permissive `NNGRRN` search tier. We state the supporting evidence exactly: SaCas9's original characterisation reports `NNGRRT` as the PAM at which it cleaves genomic targets most efficiently, but also reports degeneracy at the sixth position, concluding that all `NNGRR` PAMs can be cleaved and should be considered as potential targets (Ran et al., 2015). Our priority is therefore a weighting of the canonical tier over the permissive one, not a claim that `NNGRRN` off-targets are inert. Under the opposite priority the answer changes and does not even resolve to a single guide without a further arbitrary tie-break.

**Their ICP27 choice is vindicated.** We examined the assumption that ICP27g1 should be held fixed, and the off-target screen supports it independently. ICP27g1 has the cleanest profile of every ICP27/UL54 candidate screened: 24 `NNGRRN` sites at ≤ 4 mismatches, 1 at ≤ 3, and zero coding-exon hits. **None** of the nine better-conserved, filter-passing ICP27 alternatives matches it on all four off-target measures. The obvious "better" swap on conservation grounds, `UL54_114833−` (conservation 1.000, and it would take the pair's joint intactness to 1.000), carries 33 `NNGRRT` sites at ≤ 4 mismatches, 4 at ≤ 3, and 2 coding-exon hits — one in *FOXO6*, one in *TOX*. Changing both guides buys 4 more isolates at the cost of a materially dirtier ICP27 guide. Independently, our screen corroborates their published result: they report no sites at ≤ 3 total mismatches plus bulges for either lead guide, and under the canonical `NNGRRT` PAM our exhaustive scan finds exactly zero for both, computed from a different assembly download by a different algorithm.

### 2.6 Two of our own recommendations, withdrawn

We report two withdrawn internal recommendations, because they are the clearest available evidence for this paper's thesis.

**`RL2_3441+` — withdrawn, refuted.** Ranking the ICP0 pool on cross-isolate conservation alone identified `RL2_3441+` as the guide that should replace ICP0g2. Nothing measured in support of that was wrong: it is perfectly conserved on both denominators, passes every filter, sits in both ICP0 repeat copies, and takes the pair's joint intactness from 0.825 to 0.978. What was wrong was treating a conservation rank as a guide selection. Subsequent GRCh38 screening found it materially *dirtier* than the guide it was proposed to replace: **19** `NNGRRT` sites at ≤ 4 mismatches against ICP0g2's 5, and **6** at ≤ 3 mismatches against ICP0g2's **0** — one inside *ABL1* (chr9:130,830,349, 3 mismatches, 1 in the seed, PAM `ATGGGT`), four more copies of a single repeated sequence on chrY. Proposing a guide worse than the incumbent on an axis the incumbent was partly selected on is a straightforward failure, and the recommendation is withdrawn.

**`RL2_5080+` — withdrawn, superseded.** Ranking on escape probability alone named `RL2_5080+`. As Section 2.5 shows, this was `argmin` over four probabilities that agree to ten significant figures and differ only in floating-point summation order. It was never a real preference. `RL2_5080+` remains a respectable guide — it clears every hard gate and sits on the Pareto front — but it carries 10 `NNGRRT` sites at ≤ 4 mismatches against `RL2_5335+`'s 3, with 5 carrying ≤ 1 seed mismatch against 1.

The generalisable lesson is that **neither guide was ever uniquely indicated even within the stage that named it.** Four ICP0 sites tie exactly on the conservation figure of merit; `RL2_3441+` was named because it sorts first by identifier. A recommendation produced by a tie-break should not be reported as a finding without saying that it was one. This is the same failure mode, on a smaller scale, as selecting on a threshold that admits 85% of the candidate space.

### 2.7 The model ranks; its absolute probabilities are not measurements

The escape model has two measured inputs — the per-genome presence matrix, and a per-codon amino-acid tolerance profile derived from cross-strain variation among the same 183 isolates — and several parameters that are explicitly labelled assumptions: the NHEJ indel-length spectrum, two tolerance constants, the frameshift-viability assumption, and the repeat-copy regime. No default was chosen after seeing its effect on any conclusion. Every assumption is swept, and the sweep now includes an indel spectrum measured in post-mitotic human neurons rather than posited (point 4 below).

The result is a genuinely mixed one and we report all four parts of it together.

1. **Absolute probabilities are not usable as measurements.** Across the swept parameter space the escape probability of a fixed guide set moves by up to **3.1 orders of magnitude**, almost all of it driven by `tol_invariant` — the probability that an in-frame indel is tolerated at a codon showing no amino-acid variation across 183 isolates. That quantity is precisely what the data cannot pin down: absence of observed variation bounds the frequency of tolerated variants but does not measure the tolerance. **No number in this paper should be quoted as an escape rate.**
2. **The chosen guide *sets* are stable.** **27 of 37** parameter settings return exactly the same best *k* = 2, *k* = 3 and *k* = 4 sets as the baseline.
3. **The per-site ordering is only partly stable, and the exceptions are informative.** The rank correlation of per-site *q* falls as low as ρ = 0.441, with top-20 overlap as low as 0%, worst under `repeat_model`. Two mechanisms cause this and neither is a defect: most sites sit at exactly the same *q*, so ordering is dominated by ties and any reweighting reshuffles them; and the `repeat_model` and ICP0-frameshift-viability settings genuinely reorder ICP0 against everything else, because they are assumptions *about ICP0 specifically*. The correct reading is that **the relative merit of an ICP0 guide is not determined by the available data** — it depends on whether cutting both repeat copies must be escaped twice, and on whether an ICP0 frameshift is lethal. Both are experiments, not parameters.

4. **The one swept spectrum that is measured rather than posited moves the absolute scale and leaves the selection untouched.** The therapy targets post-mitotic trigeminal-ganglion neurons, while our default spectrum encodes only qualitative features that are uncontroversial in the field — a dominant +1 insertion, short deletions with a decaying tail — all of which come from work in dividing cells. A reviewer should press on that, so we sweep the measurement. Ramadoss et al. (2025) deposit raw CRISPResso2 indel histograms with Figure 1d of their paper — six replicates of SpCas9 RNP editing in human iPSC-derived neurons beside the genetically identical dividing iPSCs, 118,722 and 75,183 reads respectively — and we take the neuronal arm verbatim as the scenario `neuronal_nhej` (Methods 5.7). **Its in-frame fraction is 0.090, against 0.199 for our default and 0.289 for the isogenic dividing arm — and below the 0.10–0.50 band we had been sweeping.** The measured post-mitotic spectrum lay outside the range this model previously explored.

    The effect runs in the direction that helps, and it is easy to get backwards. "A narrower distribution of smaller indels" reads like a gentler outcome, but ±1 and ±2 nt account for **77%** of measured neuronal indels and every one of them frameshifts, which in an essential gene is lethal to the virus and therefore is not escape. Escape consequently *falls*, and falls harder at larger *k* because *q* multiplies across sites: the published pair goes 9.720 × 10⁻⁴ → 3.642 × 10⁻⁴ (2.7× lower), the best *k* = 2 set 6.8× lower, *k* = 3 17.7× lower, *k* = 4 46× lower. Our reconciled pair (Section 2.5) goes 2.5145 × 10⁻⁴ → 8.616 × 10⁻⁵, so its 3.9× advantage over the published pair widens to 4.2×. **Under the one spectrum measured in post-mitotic human neurons, our default is the conservative choice** — it over-states escape — which is why we leave it as the default rather than adopting the neuronal one.

    Selection does not move at all. Rank correlation of per-site *q* against baseline is ρ = 0.999996, top-20 site overlap 100%, the best sets at every *k* from 1 to 6 are identical to baseline, and the minimum-*k* answers at both the 10⁻³ and 10⁻⁶ thresholds are unchanged (*k* = 2 and *k* = 3). No conclusion in this paper changes; three of them are strengthened.

The rule for using this model is therefore: use it to choose between guide sets; do not quote its absolute escape probabilities as if they were measurements. The minimum-*k* answers inherit the absolute-scale uncertainty and must be read together with the resolution floor of Section 2.1. The closed-form calculation was cross-checked against a seeded 2,000,000-draw Monte Carlo simulation of the same generative process; the two agree within 1.2 standard errors for every set above the simulation's resolution.

**A corroborating measurement fell out of the tolerance extraction.** The coordinate map detects net in-coding-sequence length polymorphism directly. Across all seven target genes, 178 isolate–gene extractions carry a net length change, 62 of them a multiple of three — and the distribution is extremely uneven, which is the finding: **156 of the 175 isolates whose ICP0/RL2 coding sequence could be read carry a net length change in it**, against single figures for every other target gene. In-frame length variation in ICP0 is not a theoretical construct; it is the normal state of circulating, viable HSV-1, observed with no CRISPR involved. That is exactly the mechanism the escape model is about, and it lands on the gene the published design leans on hardest.

### 2.8 Delivery, not escape, dominates persistence at reported in-vivo editing rates

Escape probability is a per-genome quantity conditional on exposure to an active nuclease. It says nothing about the fraction of viral genomes that never meet one — untransduced neurons, silenced episomes, inaccessible chromatin. We denote that fraction θ and report it separately, because guide count cannot touch it: it is *k*-independent by construction and adds a floor no amount of multiplexing lowers.

**Table 5.** Persistence as a function of θ.

| θ | P(persist), published pair | P(persist), best *k* = 2 | P(persist), best *k* = 8 |
|---|---|---|---|
| 0 | 9.720 × 10⁻⁴ | 2.260 × 10⁻⁵ | 2.610 × 10⁻¹⁹ |
| 0.5 | 5.005 × 10⁻¹ | 5.000 × 10⁻¹ | 5.000 × 10⁻¹ |
| 0.9 | 9.001 × 10⁻¹ | 9.000 × 10⁻¹ | 9.000 × 10⁻¹ |
| 0.99 | 9.900 × 10⁻¹ | 9.900 × 10⁻¹ | 9.900 × 10⁻¹ |

The only in-vivo on-target measurement reported by Amrani et al. is a mean indel frequency of 0.78–1.65% in trigeminal ganglia at day 40, which they round to ~1% in discussion. Read literally as an editing rate, that corresponds to θ ≈ 0.99. At that θ, Table 5 is flat: **every guide set, at every *k*, gives essentially the same persistence, because persistence is delivery-limited and not escape-limited.**

This is not an argument against multiplexing. It is an argument that the two terms must be reported separately, which no published HSV-1 CRISPR work currently does. Escape probability governs whether the edited fraction can regenerate a resistant population under selection; θ governs how large the unedited fraction is. Conflating them lets a delivery problem read as a design problem, and makes guide-count arguments look more consequential than they are.

---

## 3. Discussion

The three findings compose into a single methodological recommendation, which we state as a checklist because that is the form in which it is usable.

**Report the interval, not the point.** Conservation measured against a finite corpus is an estimate with a one-sided lower bound that matters. A guide reported as "100% conserved" over *n* genomes should be reported as "absent in 0/*n*; 95% upper bound on population absence frequency = 1 − 0.05^(1/*n*)". For *n* = 183 that is 0.0162, and it propagates into any multiplex escape claim. The practical consequence is that escape-probability targets below about 10⁻⁴ for a two-guide design are not certifiable against any currently available HSV-1 corpus, and that a design targeting 10⁻⁶ needs four guides rather than three to be certifiable at that threshold. Because this floor scales with corpus size, it is also a concrete argument for corpus investment: doubling the number of complete genomes roughly halves the bound.

**Report joint intactness from a per-genome matrix.** A per-guide conservation table is a strictly weaker object than the matrix it summarises, and the direction of the error is not predictable a priori. In the case examined here, marginals of 0.847 and 0.978 give a joint of 0.825, at the theoretical lower bound. Computing this costs nothing beyond keeping the matrix rather than collapsing it.

**Do not select on one axis.** This is the lesson our own withdrawn recommendations teach most sharply. Conservation ranking named a guide that human off-target screening then found dirtier than the incumbent, with six near-cognate sites under the canonical PAM where the incumbent had none, one in a proto-oncogene. Escape ranking named a guide by floating-point noise across a four-way exact tie. Neither error was detectable within the stage that made it. Multi-axis selection with hard gates stated in advance, and with the discriminating priority written down rather than smuggled in, is the minimum defensible procedure.

**Separate delivery from escape.** At ~1% in-vivo editing, no guide count changes persistence. Reporting the two terms together lets a delivery limitation be read as a design achievement, or vice versa.

We want to be exact about what this analysis does and does not say about the program we audited. It says that within their own SaCas9 site space, one guide substitution improves cross-isolate conservation, joint pair intactness, modelled escape probability and human off-target burden simultaneously, and that their lead ICP0 guide is absent in 15.3% of sequenced isolates because of a single recurrent substitution their workflow was structurally unable to detect. It does **not** say that `RL2_5335+` is a better drug. Amrani et al. chose their pair on measured antiviral activity in cells — evidence of a kind this analysis contains none of, and the single most likely reason a well-founded computational recommendation could still be wrong. Their ICP27 choice is independently vindicated here, and their ICP0-duplication argument is correct under a repeat regime that our own model reproduces. The failure mode we describe is a property of the standard workflow: a conservation threshold that barely discriminates, a per-guide table that cannot see joint coverage, and an efficacy screen in a single laboratory strain that by construction contains every target site.

The obvious extensions are wet-lab ones. Cleavage of the position-9 G→A variant by ICP0g2 is directly testable and would convert a definitional claim about site identity into a functional one. Whether an ICP0 frameshift is lethal, and whether cutting both repeat copies requires two independent escapes, are the two parameters that dominate our sensitivity analysis and are experiments rather than modelling choices.

---

## 4. Limitations

We state these at length because several of them bear directly on the numbers above.

1. **Substitutions only; no bulges modelled.** The off-target screen enumerates substitution mismatches (≤ 4) under `NNGRRN` and `NNGRRT`. DNA and RNA bulges are not enumerated, and Cas-OFFinder — the field standard, and the tool Amrani et al. used (Bae et al., 2014) — does enumerate them in its current implementation, although bulge search is a later addition to the tool and is not described in the 2014 publication; off-target cleavage at sites carrying a DNA or RNA bulge is itself experimentally documented (Lin et al., 2014). A bulge-tolerant reanalysis can only add sites, never remove them, and could in principle reorder the guides. Consequently, our off-target counts are **not** comparable line-for-line to their published totals (ICP0g1 358, ICP0g2 910, ICP27g1 443, ICP27g2 316), which count mismatch-plus-bulge classes we do not enumerate and were computed against a different assembly build. They are comparable *between guides screened here*, which is the question the screen exists to answer.

2. **One human haplotype.** GRCh38 primary assembly only; ALT haplotypes were excluded to avoid double-counting sequence already present on the chromosomes. Population variation creates and destroys both PAMs and protospacer matches, so haplotype-specific off-target sites are invisible here.

3. **No activity prediction.** Conservation says a site exists in an isolate, not that SaCas9 cuts it efficiently. No CFD, MIT/Hsu or Doench-style score is applied, because none of those models is parameterised for SaCas9 and an SpCas9 score applied to SaCas9 data would be worse than no score: the mismatch-position weights underlying MIT/Hsu-style scores were fitted to SpCas9 activity in human cells (Hsu et al., 2013), and the CFD and Rule Set 2 models likewise to SpCas9 datasets (Doench et al., 2016). Seed-region mismatch counts are reported instead, unweighted and not used to filter.

4. **No wet-lab validation of anything.** No cleavage assay, no GUIDE-seq, no CIRCLE-seq, no cell or animal work (contrast the cell-based, genome-wide off-target profiling of GUIDE-seq, Tsai et al., 2015, which Amrani et al. performed). In-silico nomination is a filter, not evidence of cleavage; equally, absence of a nominated site is not proof of safety.

5. **The default indel spectrum remains a labelled assumption, and the one measured spectrum we can sweep is not measured in our cells.** The default NHEJ indel-length distribution used in the repair term is an explicitly labelled assumption with an override and a sensitivity sweep, not an attribution to an invented citation. Published empirical spectra exist — repair outcomes at Cas9-induced breaks have been profiled at scale and are reproducible but strongly target-sequence dependent rather than universal (van Overbeek et al., 2016; Chen et al., 2019) — so there is no single correct default to adopt, and substituting one published spectrum would trade a labelled assumption for an unlabelled one. It acts almost entirely through its in-frame fraction (0.199), swept over 0.10–0.50. Microhomology-mediated end joining produces a different, often larger and more predictable deletion spectrum with a different in-frame fraction, and its usage varies by cell type (Sfeir and Symington, 2015). That is not hypothetical for the target cell here, and the relevant measurement is now inside the analysis rather than outside it. Ramadoss et al. (2025) show that resting human neurons resolve Cas9-induced breaks over a substantially longer window than isogenic dividing cells, upregulate non-canonical repair factors, and favour NHEJ over MMEJ, yielding a narrower distribution of smaller indels; **their deposited Figure 1d source data is now swept as the scenario `neuronal_nhej`** (Sections 2.7, 5.7), and it lowers escape rather than raising it, without changing any selected guide set. Four caveats keep this a limitation rather than a resolution. (i) **A mouse or rabbit trigeminal ganglion is not a cultured human neuron.** Their cells are iPSC-derived neurons in culture; the therapy acts on latently infected sensory neurons in vivo in animal models, and eventually in humans. (ii) **Their nuclease is not ours.** The whole dataset is SpCas9; the design under audit uses SaCas9. Both are blunt cutters at an analogous position relative to the PAM, but we are aware of no SaCas9 spectrum in post-mitotic neurons against which to check that. (iii) **It is one guide at one locus.** The pooled histogram is sgRNA B2Mg1 at human *B2M*. The paper's own central finding is that outcome distributions are strongly guide-dependent — the perturbations that shift outcomes for three of their sgRNAs fail for a fourth, insertion-biased one — so a single-locus in-frame fraction is one draw from an unmeasured spread, not a constant. (iv) **The model is still not pathway-aware.** It applies one spectrum per run to every site identically, with no dependence on local microhomology content, which is the very mechanism generating the guide-to-guide differences. The spectrum therefore remains an assumption of this model; what has changed is that its plausible range is now anchored at one end by a measurement in the right cell class rather than by parametric guesswork, and that measurement runs in the direction that makes our default conservative.

6. **Effective sample size is 132–153, not 183, and longitudinal sampling inflates record counts.** The 183 deposited genomes behave like **132** independent ones at 99.9% estimated identity (24 multi-member clusters; the largest, 9 members, is the strain-17 family) and **153** at 99.99%; pairwise identity across the set runs 0.98134–1.00000. Note the sign: de-duplication *raises* the perfectly conserved count (833 → 857 at 99.9%), because a duplicate contributes no new sequence but can contribute new assembly noise. The threshold dominates the answer — single linkage at 99.5% would collapse 183 genomes to 20 and "raise" conservation to 1,948 guides — so effective *n* must always be quoted with its threshold. All headline numbers here use the full *n* = 183 denominator, the conservative choice. The same effect applies within the variant analysis: the 25 records carrying the ICP0g2 variant represent **13** independent clusters, not 25 individuals. Any count of records in a public corpus is an upper bound on the independent observations behind it.

7. **The genome set mixes clinical isolates with laboratory strains,** and metadata is sparse: of 183 records, 144 report a strain, 108 a country and 82 a collection date. Geographic and temporal stratification of conservation is therefore poorly supported, which is why Section 2.3 reports carrier geography for only the 8 records that have it.

8. **Assembly ambiguity depresses conservation; we quantify it rather than assume it away.** 100 of 183 genomes contain at least one non-ACGT base and 26 exceed 0.1%. Of 874,191 (guide, genome) pairs, 787,938 are present, 84,007 are absent with the homologous region demonstrably present and fully resolved, and only **2,246 (0.26%)** are genuinely unresolvable. Scoring those as unknown rather than mismatch moves the perfectly conserved count from 833 to 931. Ambiguity therefore accounts for a ~12% relative change, not for the bulk of non-conservation; for the specific ICP0g2 case, zero of the 28 absences are ambiguity gaps.

9. **The widened denominator is smaller than record counts suggest.** We re-measured every guide against 4,945 additional HSV-1 records (4,552 sub-genomic plus 393 near-full-length "partial genome" assemblies), 4,313 of which anchor-map to the reference. Only **1,021 of the 4,552** sub-genomic records touch any of the seven target genes; the rest are thymidine-kinase, glycoprotein-G and glycoprotein-B typing fragments from unrelated surveys. For six of seven genes the widened denominator is therefore dominated by the 393 partial genomes — the same *kind* of evidence as the 183 complete genomes, separated only by a title convention. Only UL30 has a genuine clinical-amplicon corpus, because UL30 is what is sequenced for aciclovir-resistance genotyping. On the sub-genomic tier alone — the only tier independent of the complete-genome set — 15 of the 833 perfect guides (1.8%) fall below 0.95, none below 0.70, minimum 0.889. "Conserved across thousands of sequences" is not supportable for six of the seven genes.

10. **Coverage determination is anchor-based, not alignment-based,** and it is the largest single effect in the analysis. A record enters a guide's denominator only when exact 25-mer anchors bracket the guide footprint within one colinear chain. Without this correction — dividing by every record considered rather than by records that actually span the site, which is how the competing paper describes its own calculation — 823 of the 833 perfect guides fall below 0.95 and mean conservation drops from 0.90 to 0.82. That is almost entirely artefact: a 600 bp amplicon does not contain a guide 40 kb away. **This cuts against us as much as for us: it means published conservation percentages computed without a coverage check are probably understated, not that ours are better.** The correction is itself conservative in one direction: a record whose homologous region is real but too divergent to anchor on either side is scored unknown rather than absent, which biases corrected conservation slightly upward.

11. **Two of our own analyses give different gene-level conservation for the same guide, and we resolve it explicitly.** The denominator-robustness analysis scores ICP0g2 at 0.807 over 445 records; the head-to-head analysis scores it at **0.8975 over 400 records**. The difference is the admission criterion: the former admits a record to a *published* guide's denominator when the record's anchor chains overlap the **gene**, the latter requires anchors to bracket the guide's own **site footprint**. The looser rule admits 45 records that reach RL2 but not the ICP0g2 footprint and scores them as absences. **We adopt the stricter per-site criterion throughout this paper**, applied uniformly to the published guides and to every alternative; it is the correct denominator, and the looser figure is reported here only so the two are not left contradicting each other.

12. **One reference genome defines the candidate space.** Sites are enumerated only from NC_001806.2 and only inside annotated coding sequence. A site conserved across every other isolate but absent from strain 17 is never considered, and a spliced-exon junction is never a candidate because the enumerator works on contiguous genomic segments.

13. **The escape model cannot see selection, fitness, kinetics, multi-cut excision or within-host structure.** Escape here is a per-genome probability of producing a viable uncuttable genome, not a rate of outgrowth; a crippled escape variant counts the same as a fit one. Two simultaneous double-strand breaks frequently excise the intervening fragment, which for guides in different genes is almost certainly lethal to the virus and would push escape *below* our estimates, so excision is an unmodelled escape-suppressing mechanism. Finally, the between-host isolate population is used as a proxy for the diversity a single patient's virus presents; within-host diversity is lower, which makes the presence term pessimistic for an individual and appropriate for a population product.

---

## 5. Methods

All analyses are deterministic: accessions are sorted before any limit is applied, all sorts are stable, and the only randomness anywhere is a seeded Monte Carlo cross-check. The implementation is pure Python with Biopython 1.88, pandas 3.0.5 and numpy 2.5.3 on CPython 3.14.5 (Windows 11). No external aligner, no compiled tooling.

### 5.1 Genome corpus

Complete HSV-1 genomes were retrieved from NCBI Nucleotide with the exact query:

```
txid10298[Organism:exp] AND 145000:160000[SLEN] AND biomol_genomic[PROP]
  AND "complete genome"[Title] NOT patent[PROP]
```

Retrieval date **2026-09-07T21:56:25Z**. NCBI reported and returned **183** records (147,898–159,092 bp). The query, the retrieval timestamp, the full sorted accession list and a SHA-256 per downloaded FASTA are recorded in the run manifest. Per-genome metadata (strain, isolate, country, collection date, host, isolation source, ambiguous-base count) was obtained by requesting a 1 bp GenBank slice per accession, which returns the full header and `source` feature.

A second, gene-level corpus was built with two further queries:

```
txid10298[Organism:exp] AND biomol_genomic[PROP] NOT patent[PROP] AND 100:145000[SLEN]
txid10298[Organism:exp] AND biomol_genomic[PROP] NOT patent[PROP]
  AND 145000:160000[SLEN] NOT "complete genome"[Title]
```

returning 4,552 sub-genomic and 393 near-full-length "partial genome" records respectively (4,945 total, 4,313 anchor-mapping to the reference).

### 5.2 Site enumeration

The annotated GenBank record for **NC_001806.2** (HSV-1 strain 17) defines the candidate space. Gene coordinates are parsed from the annotation; none are hardcoded. Target genes: UL30 (DNA polymerase), UL19 (VP5), UL5 and UL52 (helicase–primase), UL29 (ICP8), RL2 (ICP0) and UL54 (ICP27). RL2 is spliced, so each exon is treated as a separate contiguous genomic segment, guaranteeing that every reported site is a contiguous stretch of genomic DNA lying inside coding sequence. RL2 is also present twice, in the TRL and IRL repeats; the two copies yield identical sites, which are collapsed to one guide with every reference position retained and the copy number recorded.

A nuclease is defined as a spacer length plus an IUPAC PAM pattern 3′ of the protospacer; footprint length, both-strand PAM matching, scanner anchors and cut-site prediction are all derived from it. Configurations used: SpCas9 (20 nt, `NGG`, 23 nt footprint), SaCas9 (21 nt, `NNGRRT`, 27 nt), SaCas9 in the published grammar of Amrani et al. (20 nt, `NNGRRT`, 26 nt), and SaCas9 permissive (21 nt, `NNGRRN`). 5′-PAM nucleases are rejected with an explicit error rather than approximated.

### 5.3 Conservation

Conservation is exact substring presence, not alignment identity: for each candidate target site (protospacer + PAM) and each genome, does that exact sequence occur verbatim on either strand? The scan inverts the naive loop — two dictionaries keyed on the guide *k*-mers themselves (forward and reverse complement) are built once, and each genome is scanned once by jumping between occurrences of a literal anchor derived from the PAM pattern (for `NNGRRT`, the terminal `T`, chosen per sequence as the rarest available literal). Work is O(*S* × *L*) with memory O(*G*). A PAM with no literal position falls back to an exhaustive positional scan. Equivalence to a naive both-strand substring search is asserted for every supported PAM on randomised sequences, GC-rich sequences and the real reference.

Standard filters, applied identically everywhere: no `TTTT` or longer run in the protospacer (RNA polymerase III terminates within a run of thymidines; in the type-3 U6, H1 and 7SK promoters used to express sgRNAs, T₄ is the minimal termination signal and termination becomes efficient only at T₆ or longer (Gao et al., 2018), so a U6/H1-driven sgRNA whose spacer contains `TTTT` is partly transcribed truncated — we treat this as a hard exclusion, not a penalty); longest homopolymer run ≤ 4; GC fraction within 0.35–0.75.

### 5.4 Coverage-corrected conservation on the gene-level corpus

Homology to the reference is established without an aligner, by exact 25-mer anchors (step 4, maximum gap 3,000 bp, maximum drift 60, maximum 4 occurrences per anchor) and greedy colinear chaining — the seed-and-chain half of a standard alignment pipeline. **A record counts toward a guide only if one chain has an anchor entirely left of the guide footprint and another entirely right of it**, so the two anchors bracket the homologous stretch exactly with no padding constant. Everything between two anchors of one chain counts as covered even where divergent. A non-matching bracketed pair is scored ABSENT if the bracket is short (≤ 2,000 bp) and fully resolved, and UNKNOWN if it contains a non-ACGT base, is implausibly wide, or the record does not reach both sides. Only PRESENT and ABSENT enter the corrected denominator; a minimum of 10 covering records is required for a guide to be reported on this corpus.

### 5.5 Effective sample size and ambiguity sensitivity

Redundancy: a bottom-5000 MinHash sketch of canonical 21-mers per genome, using splitmix64 mixing rather than Python's salted `hash()` so sketches are reproducible across processes; all-pairs Jaccard by the mash estimator; Jaccard converted to an ANI-like identity; single-linkage clustering (Ondov et al., 2016). Conservation is then recomputed on one representative per cluster. The identity grid (0.98, 0.99, 0.995, 0.998, 0.999, 0.9999) is a module-level constant fixed before the analysis. Primary threshold 0.999.

Ambiguity: a sweep over a maximum ambiguous-base fraction per genome, plus an N-tolerant match mode in which a site whose homologous region in a genome is unresolved is scored UNKNOWN and dropped from that guide's denominator rather than counted as a mismatch.

### 5.6 The ICP0g2 variant analysis

For the ICP0g2 target site (26-mer `CTCAGGCCGCGAACCAAGAA` + `CAGAGT`, reference positions 4,496–4,521), the 30 bp immediately flanking the footprint in the reference were used as left and right anchors. In each of the 28 genomes lacking an exact match, both strands were searched for both anchors; where both were found with an intervening span under 200 bp, the bracketed sequence was extracted and compared base by base to the reference site. Genomes in which both anchors could not be located were recorded as unresolved, not absent. Carrier independence was assessed by mapping each carrier accession to its 99.9%-identity cluster from the redundancy analysis.

### 5.7 Escape model

For one viral genome drawn from the empirical isolate population, conditional on exposure to an active nuclease:

P(escape | S) = (1/N) Σ_genomes Π_{sites i ∈ S} f(g, i), where f(g, i) = 1 if site *i* is absent in genome *g* (it can never be cut, and the isolate is viable by construction) and f(g, i) = q_i if present.

q_i = p_disrupt × Σ_{L ≠ 0} [P(L)/(1 − P(0))] × viability(L, i), with viability equal to the gene's frameshift viability when L mod 3 ≠ 0 and, when L mod 3 = 0, the mean over placements of the product of per-codon tolerances over the |L|/3 codons deleted or inserted around the cut codon. A codon is scored tolerant (`tol_variable` = 0.5) if at least one of the 183 isolates carries a different amino acid there, and intolerant (`tol_invariant` = 0.05) otherwise. Conditioning on L ≠ 0 is not a convenience: a persistently expressed nuclease re-cuts a perfectly repaired site, so the terminal state of a genome that keeps its target intact is continued exposure, not escape.

For a gene present in *c* copies, the default `redundant` regime gives q = 1 − (1 − q_single)^c; the `single` and `all-copies` regimes are also reported. Baseline frameshift viability is 0.0 for every gene — the assumption most favourable to a two-guide design. Minimum separation between two same-gene sites is 150 bp, so that two sites close enough for one repair event to destroy both are never scored as independent.

The per-codon tolerance profile is measured: every codon of every isolate is read through an alignment-free reference→isolate coordinate map (the same 25-mer anchor and chaining machinery, with maximal constant-offset runs) and translated. An isolate contributes nothing to a gene unless ≥ 95% of that gene's codons map and the translated protein contains no internal stop, so a broken coordinate map cannot manufacture apparent tolerance. The inference is one-directional and should be read as such: observed variation proves tolerance; absence of observed variation over 183 isolates bounds the frequency of tolerated variants but does not prove intolerance.

Guide-set search is exhaustive for *k* ≤ 2 and, for larger *k*, multi-start greedy plus exhaustive local swap with fixed scan order and index-based tie-breaks; at *k* = 3 the heuristic was checked against brute force over the best 60 sites (34,220 sets) and matched exactly. The closed form is cross-checked by a Monte Carlo drawing a genome uniformly from the empirical population and then a Bernoulli repair outcome per present site, `numpy.random.default_rng(20240814)`, 2,000,000 draws.

Certifiable bounds replace each site's absence count *a* out of *n* by the Clopper–Pearson 95% upper limit, which for *a* = 0 is 1 − 0.05^(1/n) = 0.0162 at *n* = 183.

Sensitivity: every assumption is swept one at a time holding the rest at baseline — `tol_variable` over 0.10–1.00, `tol_invariant` over 0.001–0.50, minimum variant strains 1–5, six indel spectra, in-frame fraction 0.10–0.50, three repeat regimes, ICP0 frameshift viability 0–1, and `p_disrupt` 0.8–1.0 — for 37 settings in total.

Five of the six swept spectra are parametric assumptions. The sixth, `neuronal_nhej`, is measured: it is the pooled CRISPResso2 net-length histogram (118,722 reads over six replicates) deposited as source data with Figure 1d of Ramadoss et al. (2025), from SpCas9 RNP editing of human iPSC-derived post-mitotic neurons, archived in this repository as `refs/ramadoss2025_fig1d_indel_histogram.tsv` with the extraction method and its limits in `refs/ramadoss2025_notes.md`. It enters as a scenario and not as the default, for the reasons given in Limitation 5; the default spectrum and every default-parameter output are unchanged by its addition.

### 5.8 Human off-target screening

Reference: **Ensembl release 116**, GRCh38 primary assembly, **unmasked** (`Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz`; 881,964,081 bytes; SHA-256 `d8c3af0094a7bba6125763bad779ec18a81483c739c6ed122094bdf86c187b92`), 194 sequences, **3,099,750,718 bases scanned** (Yates et al., 2026). Unmasked rather than soft- or hard-masked, because a Cas9 site inside a LINE, Alu or satellite is a real cleavage substrate; primary assembly rather than toplevel, because toplevel adds ALT haplotypes and patch scaffolds that would double-count off-targets. Integrity was verified three independent ways: Ensembl's published BSD `sum` checksum, the gzip CRC-32/ISIZE trailer, and a locally computed SHA-256. The release is pinned, not "current".

Search: PAM `NNGRRN`, with `NNGRRT` sites flagged as its strict subset so both come out of one pass; up to **4 mismatches** in the protospacer; both strands; spacer lengths 20 and 21 nt. Substitutions only. Ambiguous genomic bases count as mismatches and disqualify a PAM outright, so no site is ever credited to unresolved sequence.

Algorithm: three vectorised slice comparisons discard 15/16 of the genome using PAM sparsity; the pigeonhole principle (a 20-nt protospacer cut into five 4-nt chunks must match exactly in at least one chunk under ≤ 4 mismatches) prunes ~98% of what remains, resolved for all guides simultaneously by 256-entry bitmask lookup tables; survivors are re-read from the genome and counted exactly. Correctness is asserted against a structurally independent brute-force scan over the whole of chromosome 21 for every guide screened, itself validated against a plainly written character-by-character implementation. Measured scan time 319 s for 50 guides × 2 strands × 3.10 Gb.

Sites at ≤ 3 mismatches were annotated against the Ensembl 116 GTF at gene, exon and CDS granularity, on the same chromosome names the FASTA uses, so no liftover is involved.

Guides screened, with the selection rule fixed before any off-target number existed: the four published guides; every ICP0 site beating ICP0g2 on conservation *and* passing the standard filters (12 sites); and every ICP27 site passing those filters and outranking ICP27g1 within its own gene (9 sites).

### 5.9 Reconciliation and selection rule

The three axes were joined on guide identifier and the selection rule of Section 2.5 applied mechanically from constants in the reconciliation module, so that the written rule and the arithmetic cannot drift apart. The reconciliation refuses to report anything unless its rebuilt escape model reproduces the canonical model's published lead-pair P(escape) to within 10⁻¹²; observed 0.00097198790268153488 against 0.0009719879026815389. Pareto dominance over the nine-axis set is computed, not asserted.

---

## 6. Figures

**Figure 1. Escape probability versus guide count, point estimate and certifiable bound.**
P(escape) per exposed viral genome for the best guide set at each *k* from 1 to 8, on a logarithmic axis. The lower trace is the point estimate from measured site absences in the 183-genome presence matrix; the upper trace re-runs the identical model with every site's absence frequency replaced by its Clopper–Pearson 95% upper limit (0.0162 for a site absent in 0/183). Horizontal reference lines mark per-genome thresholds of 10⁻³ and 10⁻⁶. The traces diverge with *k*, from 4.4-fold at *k* = 1 to 374-fold at *k* = 4, because each additional site multiplies in another bounded absence term. The key feature is the crossing of the 10⁻⁶ line: the point estimate clears it at *k* = 3 (1.075 × 10⁻⁷) while the certifiable bound does not until ***k* = 4** (1.913 × 10⁻⁷, against 9.148 × 10⁻⁶ at *k* = 3). The published lead pair is marked separately at 9.720 × 10⁻⁴ (bound 1.137 × 10⁻²). The limit illustrated is a property of corpus size, not of the repair model: it moves only when more genomes are sequenced.

**Figure 2. Per-genome presence matrix for the published lead pair, and marginals versus joint.**
(a) A 183-row binary matrix for the two guides of the EBT-104 lead pair, ICP0g2 (RL2_4496+) and ICP27g1 (UL54_115156+), one row per complete HSV-1 genome, ordered so failures are visible; a filled cell denotes exact presence of the full protospacer-plus-PAM on either strand. Twenty-eight genomes lack ICP0g2 and four lack ICP27g1, and **the two failure sets do not overlap**: no genome lacks both. (b) The consequence, as a bar chart: marginals 0.8470 (155/183) and 0.9781 (179/183); the product of marginals, 0.8285, which independence would predict; and the measured joint intactness, **0.8251** (151/183) — below either marginal and essentially exactly at the theoretical lower bound *a* + *b* − 1. The gap between product and measured joint is the quantity no per-guide table can express; assuming independence would over-state this pair's escape probability 4.4-fold.

**Figure 3. Conservation distribution across the SpCas9 and SaCas9 candidate spaces.**
Cumulative distribution of conservation fraction across the 183 complete genomes for every candidate target site in the seven-gene panel, for SpCas9 `NGG` (20 nt spacer, *n* = 4,777 sites) and SaCas9 `NNGRRT` (21 nt spacer, *n* = 448 sites). The curves have nearly the same shape — 17.4% of SpCas9 and 15.0% of SaCas9 sites are present in all 183 genomes — but the SaCas9 curve is drawn from an order of magnitude fewer sites. Inset: absolute counts at each tier (perfectly conserved / ≥ 0.99 / ≥ 0.95: SpCas9 833 / 1,465 / 3,431; SaCas9 67 / 104 / 290) and counts surviving the poly-T, homopolymer and GC filters (644 versus 45). The point is the search-space contraction, not the fraction: requiring `NNGRRT` leaves roughly a tenth of the options an SpCas9 design has.

**Figure 4. Multi-axis comparison of ICP0 candidates, with Pareto front.**
The 14 ICP0/RL2 candidates screened on all three axes: conservation across the 183 complete genomes on the *x* axis, canonical `NNGRRT` off-target burden at ≤ 4 mismatches (lower is better) on the *y* axis, marker area proportional to local 200 bp GC. Points failing any hard gate are open; gate survivors filled. Three guides are labelled: **ICP0g2** (RL2_4496+; conservation 0.847, 5 sites, local GC 0.840 — the published clinical guide), **RL2_5335+** (1.000, 3 sites, local GC 0.725 — the recommendation here), and **RL2_3441+** (1.000, 19 sites of which 6 at ≤ 3 mismatches including one in *ABL1*, local GC 0.635 — named on conservation alone and subsequently withdrawn). The Pareto front over the nine-axis set {conservation on three denominators, P(escape), `NNGRRT` ≤ 3 and ≤ 4, `NNGRRN` ≤ 3 and ≤ 4, coding-exon hits} has three members — RL2_5335+, RL2_3364+, RL2_5080+ — and is outlined; RL2_5335+ Pareto-dominates both ICP0g2 and RL2_3441+ outright. No candidate dominates on every axis, which is why an explicit priority between the two PAM tiers is required and is stated.

**Figure 5. The recurrent ICP0g2 variant: position in the spacer, and distribution of carriers.**
(a) The 26-nt ICP0g2 target site (20 nt spacer + `CAGAGT` PAM) drawn base by base, with the recurrent substitution marked at **spacer position 9 (G→A)**. All 25 carrier genomes carry the identical substitution; the PAM is unchanged in every case. (b) Disposition of the 28 genomes lacking the exact site, from the flanking-anchor readout: **25 real single-substitution variants, 0 ambiguity (N) gaps, 3 unresolved** (anchors not locatable). There is no evidence the gap is an assembly artefact. (c) Carrier metadata where available: collection year for the 8 carriers with a recorded date, spanning **1967–2020**, and country for the same 8 (7 USA, 1 China); 17 report neither. (d) Independence: the 25 carrier records collapse to **13 distinct clusters** at 99.9% estimated identity, several carrier sets being longitudinal or multi-site samples from one subject.

---

## 7. Data and Code Availability

All analyses are computational and use only public data.

**Input data.** 183 complete HSV-1 genomes from NCBI Nucleotide, retrieved 2026-09-07T21:56:25Z with the exact query in Section 5.1; 4,945 further HSV-1 records for the gene-level corpus, with the two queries given there; the annotated GenBank record NC_001806.2; the Ensembl release-116 GRCh38 primary assembly (unmasked) and the matching Ensembl 116 GTF, with the URL and SHA-256 in Section 5.8.

**Provenance records sufficient to reproduce every number.** The per-genome manifest (accession, length, ambiguous-base count, strain, country, collection date, host, isolation source, definition, completeness label, SHA-256 of the downloaded FASTA, plus the query and retrieval timestamp); the per-record gene-corpus manifest; and the run log (full sorted accession list, exact Entrez query, retrieval timestamp, command line, Python and package versions, stages executed, summary statistics, and the genome manifest and checksums for the human assembly). These files plus the code fully determine the results.

**Code.** All analysis code is public at `https://github.com/ivanjalid1/hsv1-crispr-escape` and permanently archived at Zenodo (Heredia Jalid, 2026), doi:`10.5281/zenodo.22664837`, which resolves at `https://doi.org/10.5281/zenodo.22664837`. That is the concept DOI and always resolves to the most recent archived version; the release reported here is v1.0.0. The pipeline is dependency-light — Python plus Biopython, pandas and numpy — with no external aligner and no compiled tooling. Every stage is separately runnable and offline-reproducible from cached downloads. Correctness of the two performance-critical scanners is proved by set-identity assertions against structurally independent brute-force implementations (all supported PAMs for the conservation scanner; the whole of chromosome 21 for the off-target scanner).

**Derived tables.** Per-guide conservation and ranking; redundancy, ambiguity and gene-level corpus tables; the SaCas9 candidate pool and all pairwise joint conservations; per-site escape parameters, the *k*-curve, the guide-set table, the measured per-codon tolerance profile and the full sensitivity sweep; per-guide off-target counts at each mismatch level under both PAM variants, every site found, and genomic context for sites at ≤ 3 mismatches; and the joint reconciliation table. These are regenerated by the pipeline rather than version-controlled, so that a stale copy cannot be mistaken for a current one.

---

## 8. Author Contributions

**Ivan Heredia Jalid** (Universidad Católica de Córdoba, Córdoba, Argentina; ORCID 0009-0003-6702-1295) is the sole author and conceived the study, wrote all code, performed all analyses, and wrote the manuscript.

## 9. Competing Interests

The author declares no competing financial or non-financial interests. This work was not funded by, commissioned by, or performed in consultation with any commercial entity, including any party with an interest in the program analysed as a case study.

## 10. Acknowledgements

We thank the authors of Amrani et al. (2024) for publishing their guide sequences, selection criterion and off-target pipeline in enough detail to be independently re-measured. An audit of this kind is only possible against work that is fully disclosed, and the disclosure is to their credit.
