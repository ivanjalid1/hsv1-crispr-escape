# A quantitative escape-probability model for multiplex CRISPR editing of HSV-1

Generated 2026-09-08T16:12:47Z · nuclease sacas9 (20 nt spacer + NNGRRT PAM) · 183 complete genomes · 449 candidate sites in 7 genes.

## 0. Headline

- **The Amrani et al. lead pair (ICP0g2 + ICP27g1) scores P(escape) = 9.720e-04 per exposed viral genome** under the baseline parameters, driven almost entirely by pre-existing site absence: the pair is intact in only 0.825 of the 183 sequenced isolates.
- The best two-guide set this pipeline can propose from the same SaCas9 site space scores 2.260e-05 -- 43x lower -- **without changing the number of guides**. The dominant term at k=2 is guide *choice*, not guide *count*.
- Minimum k on the point estimate: **k = 2** for P(escape) < 1e-3, **k = 3** for < 1e-6. But 183 genomes cannot CERTIFY those numbers: a site absent in 0/183 isolates is still consistent with a population absence frequency up to 0.0162, and at that bound the same thresholds need **k = 2** and **k = 4**. **The sampling resolution of the genome corpus, not the NHEJ repair model, is what limits how low an escape probability can be demonstrated.** Section 7.1.
- **Is the published two-guide design adequate?** On this model it is adequate at a 1e-3 per-genome threshold and not adequate at 1e-6, and the certifiable bound (1.137e-02) is weaker still. The failure is specific and fixable: 96% of its escape probability comes from the 32 isolates that have already lost ICP0g2, so the same k=2 architecture with a better-conserved ICP0 guide clears the same thresholds. Adding guides is not the only remedy and is not the cheapest one.
- Escape probability is NOT the product of marginal conservations. The empirical joint matrix is used throughout; section 7.3 gives the exact size and SIGN of the error an independence assumption would have made -- for this pair, independence would have over-stated escape 3.4-fold.
- Absolute escape probabilities move over 3.1 orders of magnitude across the plausible parameter space, while the chosen guide SETS are stable across 27 of 37 parameter settings. That combination -- stable ranking, unstable absolute scale -- is the honest and still-useful result, and section 8.1 says so at length.
- **The indel spectrum measured in post-mitotic human neurons is now one of those settings, and it cuts our way.** Ramadoss et al. 2025's deposited Fig. 1d source data, swept here as `spectrum=neuronal_nhej`, has an in-frame fraction of 0.090 -- below the 0.10-0.50 band swept previously, because +-1 and +-2 nt account for 77% of neuronal indels and every one of them frameshifts. Escape FALLS (2.7x for the published pair), not rises, and the selected sets are unchanged at every k. The default spectrum is deliberately left where it was: it is the more pessimistic of the two, and the neuronal measurement is in cultured human neurons cut by SpCas9, not in a trigeminal ganglion cut by SaCas9.
- One conclusion runs the competitor's way. Under the repeat model in which both ICP0 copies must independently produce a viable escape allele, a single ICP0 guide is worth two guides elsewhere -- which is exactly the 'ICP0 is duplicated, so two guides give three DSBs' argument, quantified. Our default does not assume that regime, and which regime holds is the highest-value experiment this analysis points to.

## 1. The model


For one viral genome drawn from the empirical isolate population, and conditional on
that genome being exposed to an active nuclease:

    P_escape(S) = (1/N) * SUM over genomes g of  PROD over sites i in S of f(g, i)

    f(g, i) = 1     if site i is ABSENT in genome g   -- it can never be cut,
                                                         and the isolate is viable
                                                         by construction
            = q_i   if site i is PRESENT in genome g  -- it must acquire a repair
                                                         product that is both
                                                         uncuttable and viable

    q_i = p_disrupt * SUM over L != 0 of  P(L)/(1 - P(0)) * viability(L, i)

    viability(L, i) = frameshift_viability[gene(i)]        if L mod 3 != 0
                    = mean over placements of PROD t_j     if L mod 3 == 0
                      over the |L|/3 codons deleted/inserted around the cut codon

    t_j = tol_variable    if codon j varies in amino acid across the isolates
        = tol_invariant   otherwise

For a gene present in `c` copies (RL2/ICP0 sits in both the TRL and IRL repeats), the
site fails only if EVERY copy has become uncuttable, and the virus survives only if at
least one copy is still functional; under the default `redundant` model that gives
q = 1 - (1 - q_single)^c.

Conditioning on L != 0 is not a convenience: a persistently expressed nuclease re-cuts
a perfectly repaired site, so the terminal state of a genome that keeps its target
intact is not escape but continued exposure. It also means `p_disrupt`, being constant
in L, cancels out of the conditional distribution -- it is kept explicit so that the
assumption is visible.

Two limiting cases are asserted in `tests/test_escape.py` because they are what makes
the formula checkable rather than merely plausible:
  * k = 1 with q = 0 gives exactly 1 - conservation_fraction(i);
  * sites present in every genome give exactly PROD q_i, so perfectly conserved
    independent sites multiply.

## 2. Parameters, their provenance and the direction of their effect

| parameter | value | source | effect if wrong |
|---|---|---|---|
| site presence matrix | 183 genomes x 449 sites | MEASURED (stage 3 semantics) | dominates P_escape at small k; under-states conservation if assemblies are unresolved, which inflates escape |
| per-codon amino-acid variation | 183-isolate profile | MEASURED (this module) | sets which in-frame indels are survivable; a wrong frame would randomise it, so the frame is asserted by translation |
| indel length spectrum | parametric-default | ASSUMPTION -- no citation invented | acts almost entirely through its in-frame fraction (0.199); swept 0.10-0.50, plus the measured post-mitotic-neuron scenario neuronal_nhej (0.090), which sits BELOW that range |
| tol_variable | 0.5 | ASSUMPTION | scales q at variable codons; swept 0.10-1.00 |
| tol_invariant | 0.05 | ASSUMPTION | scales q at invariant codons; swept 0.001-0.50 |
| min_variant_strains | 1 | ASSUMPTION | moves codons between the two tolerance classes; swept 1-5 |
| frameshift viability | 0.0 for every gene | ASSUMPTION (baseline is the one most favourable to a 2-guide design) | RL2/ICP0 is not strictly essential; setting RL2=1.0 raises escape at every ICP0 site by ~1/(in-frame fraction) |
| repeat model | redundant | ASSUMPTION (all three variants reported) | changes ICP0 site q by up to two orders of magnitude |
| p_disrupt | 1 | ASSUMPTION | cancels in the conditioning while constant in length; only a length-dependent version would matter |
| min site separation | 150 bp | MODELLING CHOICE | prevents two same-gene sites close enough for one repair event to destroy both being scored as independent |

The only parameters that are MEASURED are the ones taken from this repository's own genome corpus. Everything labelled ASSUMPTION is swept in section 8 and no default was chosen after seeing its effect on the conclusion.

## 3. Measured input A -- the per-genome site presence matrix


The matrix is 183 genomes x 449 sites of exact protospacer+PAM presence on
either strand -- identical semantics to stage 3, so the column means ARE this
repository's published `conservation_fraction`. It is kept per genome, never collapsed
to marginals, because that is the only way correlated failure across sites is
representable. MEASURED.

| statistic | value |
|---|---|
| genomes | 183 |
| sites | 449 |
| sites present in all genomes | 69 |
| median site conservation | 0.9727 |
| mean site conservation | 0.8838 |
| distinct presence patterns (rows) | 141 |

## 4. Measured input B -- per-codon tolerance from cross-strain variation


For each gene, every codon of every one of the 183 isolates is read through an
alignment-free reference->isolate coordinate map (exact 25-mer anchors, greedy
colinear chaining, maximal constant-offset runs -- the stage-5 machinery, reused
verbatim) and translated. A codon counts as TOLERANT if at least
1 isolate carries a different amino acid there: a viable
clinical isolate carrying a substitution is a direct demonstration that the position
admits change. An isolate contributes nothing to a gene unless >=95% of that gene's
codons map and the translated protein has no internal stop, so a broken coordinate map
cannot manufacture apparent tolerance. MEASURED.

The inference actually made is one-directional and should be read as such: observed
variation proves tolerance; absence of observed variation over 183 isolates bounds
the frequency of tolerated variants but does not prove intolerance. `tol_invariant`
is exactly that unknown, and it is swept over three orders of magnitude in section 8.

| gene | essential | codons | isolates used | isolates dropped (QC) | variable codons | variable fraction | median isolates resolved/codon | isolates with in-CDS length polymorphism | ... of which in-frame |
|---|---|---|---|---|---|---|---|---|---|
| RL2 | False | 776 | 175 | 8 | 84 | 0.1082 | 175 | 156 | 45 |
| UL19 | True | 1375 | 182 | 1 | 89 | 0.06473 | 182 | 4 | 3 |
| UL29 | True | 1197 | 181 | 2 | 68 | 0.05681 | 181 | 1 | 1 |
| UL30 | True | 1236 | 183 | 0 | 127 | 0.1028 | 183 | 6 | 5 |
| UL5 | True | 883 | 182 | 1 | 53 | 0.06002 | 182 | 4 | 2 |
| UL52 | True | 1059 | 183 | 0 | 75 | 0.07082 | 183 | 6 | 6 |
| UL54 | True | 513 | 179 | 4 | 47 | 0.09162 | 179 | 1 | 0 |


**An independent corroboration falls out of the same extraction.** Length polymorphism
inside the coding sequence is detected directly by the coordinate map (a change of
offset between the start and the end of a CDS is an indel in that isolate).
Summed over all seven genes, 178 isolate-gene extractions carry a net
in-CDS length change and 62 of those are a multiple of three. The
distribution is extremely uneven and the uneven part is the finding: see the
"isolates with in-CDS length polymorphism" column above. **156 of the
175 isolates whose ICP0/RL2 coding sequence could be read carry a net length
change in it**, against single figures for every other target gene. In-frame length
variation in ICP0 is therefore not a theoretical construct -- it is the normal state
of circulating, viable HSV-1, observed with no CRISPR involved.

That is exactly the mechanism the escape model is about, and it lands on the gene the
published design leans on hardest.

Note the direction of the resulting bias. Isolates whose gene is length-polymorphic
are the *most* indel-tolerant ones, and they are also the ones most likely to be
dropped by the coordinate-map QC. The tolerance profile is therefore, if anything, an
under-estimate, which makes the escape estimates conservative in the direction that
flatters a small-k design.

## 5. Per-site escape probabilities

| gene | sites | median_q | min_q | max_q | median_local_tolerance | copies | essential |
|---|---|---|---|---|---|---|---|
| RL2 | 46 | 0.009486 | 0.009486 | 0.1148 | 0.09091 | 2 | False |
| UL19 | 77 | 0.004754 | 0.004754 | 0.05422 | 0.05 | 1 | True |
| UL29 | 63 | 0.004754 | 0.004754 | 0.06061 | 0.05 | 1 | True |
| UL30 | 93 | 0.004754 | 0.004754 | 0.05912 | 0.05 | 1 | True |
| UL5 | 58 | 0.004754 | 0.004754 | 0.04754 | 0.05 | 1 | True |
| UL52 | 77 | 0.004754 | 0.004754 | 0.05289 | 0.09091 | 1 | True |
| UL54 | 35 | 0.004754 | 0.004754 | 0.05389 | 0.09091 | 1 | True |


`q_site` is the per-site repair-escape probability. It is small everywhere -- the
median across the pool is 0.00475 -- because a
frameshift in an essential gene is not escape and because the in-frame fraction of the
NHEJ spectrum is only 0.199. The consequence is structural
and worth stating plainly: **at every site with even moderate cross-isolate
conservation, pre-existing absence dominates repair escape by orders of magnitude.**
Multiplexing buys protection mainly by covering isolates that a single guide misses,
not by making NHEJ escape combinatorially unlikely.

RL2/ICP0 is the exception in two directions at once. It is duplicated, so under the
default redundant-copy model a broken copy can be absorbed; and it is the one target
gene that is not strictly essential, so the baseline assumption that a frameshift
kills the virus is least secure exactly there. Both effects raise escape at ICP0
sites, and both are swept in section 8.

## 6. The Amrani et al. lead pair, scored

| amrani guide | our site id | gene | conservation (n=183) | cut codon | codon varies | q_site | P(escape) alone |
|---|---|---|---|---|---|---|---|
| ICP0g1 | RL2_5319+ | RL2 | 0.9836 | 724 | False | 0.009506 | 0.02574 |
| ICP0g2 | RL2_4496+ | RL2 | 0.847 | 450 | False | 0.009486 | 0.161 |
| ICP27g1 | UL54_115156+ | UL54 | 0.9781 | 479 | False | 0.004754 | 0.02651 |
| ICP27g2 | UL54_114376- | UL54 | 0.9727 | 216 | False | 0.004754 | 0.03195 |
| LEAD PAIR ICP0g2+ICP27g1 | RL2_4496++UL54_115156+ | RL2+UL54 | 0.8251 | - | - | - | 0.000972 |

**Where that number comes from.** The three isolate classes contribute separately and the split is exact, not a fit:

| isolate class | isolates | contribution to P(escape) | share |
|---|---|---|---|
| both sites already absent (escape is free) | 0 | 0 | 0 |
| one site already absent (one repair escape needed) | 32 | 9.348e-04 | 0.9617 |
| both sites intact (two repair escapes needed) | 151 | 3.721e-05 | 0.03829 |


The pair's escape probability is 9.720e-04 per exposed viral genome, and the
decomposition above is the whole argument in one table. 96.2%
of it comes from the 32 isolates in which ONE of the two sites is already
absent before any editing -- in those isolates the design is effectively a one-guide
design, and only a single NHEJ escape event is required. A further
3.8% comes from the 151 isolates where both
sites are intact and two independent escapes are needed; 0.0%
comes from isolates where neither site is present at all.

That is the number their design argument needed and never computed, and it says
something specific: **the weakness of the published pair is not that k=2 is too few,
it is that one of the two guides is absent in 15.3%
of sequenced isolates.** Escape at k=2 here is dominated by pre-existing sequence
variation, not by NHEJ.

The comparison that matters is therefore not against k=3 or k=4 but against the best
pair available in the same site space and the same nuclease grammar, which section 6.1
gives. If a better-chosen pair reaches the escape probability that the published pair
would need three or four guides to reach, then the honest conclusion is that the
two-guide architecture was never the problem.

### 6.1 Alternatives at the same k

| set | k | guides | genes | joint conservation | P(escape) | fold better than Amrani pair |
|---|---|---|---|---|---|---|
| Amrani lead pair (as published) | 2 | RL2_4496+;UL54_115156+ | RL2;UL54 | 0.8251 | 0.000972 | 1 |
| best k=2 (this model) | 2 | UL19_36407+;UL19_37448- | UL19 | 1 | 2.26e-05 | 43 |
| best k=3 (this model) | 3 | UL19_36407+;UL19_37448-;UL19_37908+ | UL19 | 1 | 1.075e-07 | 9045 |
| best k=4 (this model) | 4 | UL19_36407+;UL19_37448-;UL19_37908+;UL19_38322+ | UL19 | 1 | 5.109e-10 | 1.903e+06 |
| best k=2 restricted to ICP0+ICP27 (their gene choice) | 2 | UL54_114638+;UL54_114833- | UL54 | 1 | 2.26e-05 | 43 |
| best k=2, at most 1 guide(s) per gene | 2 | UL19_36407+;UL29_58709- | UL19;UL29 | 1 | 2.26e-05 | 43 |
| best k=3, at most 1 guide(s) per gene | 3 | UL19_36407+;UL29_58709-;UL30_64524+ | UL19;UL29;UL30 | 1 | 1.075e-07 | 9045 |
| best k=4, at most 1 guide(s) per gene | 4 | UL19_36407+;UL29_58709-;UL30_64524+;UL5_12570+ | UL19;UL29;UL30;UL5 | 1 | 5.109e-10 | 1.903e+06 |
| best pair keeping ICP27g1 (partner from any gene) | 2 | UL19_36407+;UL54_115156+ | UL19;UL54 | 0.9781 | 0.000126 | 7.713 |
| best pair keeping ICP27g1, ICP0 guide swapped only | 2 | RL2_5080+;UL54_115156+ | RL2;UL54 | 0.9781 | 0.0002515 | 3.865 |

## 7. Escape probability versus number of guides

**Unconstrained best set at each k** -- the model is free to put every guide in the same gene if that is optimal, and at small k it does.

| k | P(escape), best k-set | P(escape) if presence were independent | joint conservation of the set | guides | genes | below_0.001 | below_1e-06 |
|---|---|---|---|---|---|---|---|
| 1 | 4.754e-03 | 4.754e-03 | 1 | UL19_36407+ | UL19 | False | False |
| 2 | 2.260e-05 | 2.260e-05 | 1 | UL19_36407+;UL19_37448- | UL19;UL19 | True | False |
| 3 | 1.075e-07 | 1.075e-07 | 1 | UL19_36407+;UL19_37448-;UL19_37908+ | UL19;UL19;UL19 | True | True |
| 4 | 5.109e-10 | 5.109e-10 | 1 | UL19_36407+;UL19_37448-;UL19_37908+;UL19_38322+ | UL19;UL19;UL19;UL19 | True | True |
| 5 | 2.429e-12 | 2.429e-12 | 1 | UL19_36407+;UL19_37448-;UL19_37908+;UL19_38322+;UL19_39276+ | UL19;UL19;UL19;UL19;UL19 | True | True |
| 6 | 1.155e-14 | 1.155e-14 | 1 | UL19_36407+;UL19_37448-;UL19_37908+;UL19_38322+;UL19_39276+;UL29_58709- | UL19;UL19;UL19;UL19;UL19;UL29 | True | True |
| 7 | 5.490e-17 | 5.490e-17 | 1 | UL19_36407+;UL19_37448-;UL19_37908+;UL19_38322+;UL19_39276+;UL29_58709-;UL29_61292+ | UL19;UL19;UL19;UL19;UL19;UL29;UL29 | True | True |
| 8 | 2.610e-19 | 2.610e-19 | 1 | UL19_36407+;UL19_37448-;UL19_37908+;UL19_38322+;UL19_39276+;UL29_58709-;UL29_61292+;UL29_61578+ | UL19;UL19;UL19;UL19;UL19;UL29;UL29;UL29 | True | True |

**Best set at each k under a design policy of at most 1 guide(s) per gene**, which is what a real multiplex construct would be built to.

| k | P(escape), best k-set | P(escape) if presence were independent | joint conservation of the set | guides | genes | below_0.001 | below_1e-06 |
|---|---|---|---|---|---|---|---|
| 1 | 4.754e-03 | 4.754e-03 | 1 | UL19_36407+ | UL19 | False | False |
| 2 | 2.260e-05 | 2.260e-05 | 1 | UL19_36407+;UL29_58709- | UL19;UL29 | True | False |
| 3 | 1.075e-07 | 1.075e-07 | 1 | UL19_36407+;UL29_58709-;UL30_64524+ | UL19;UL29;UL30 | True | True |
| 4 | 5.109e-10 | 5.109e-10 | 1 | UL19_36407+;UL29_58709-;UL30_64524+;UL5_12570+ | UL19;UL29;UL30;UL5 | True | True |
| 5 | 2.429e-12 | 2.429e-12 | 1 | UL19_36407+;UL29_58709-;UL30_64524+;UL5_12570+;UL52_109245+ | UL19;UL29;UL30;UL5;UL52 | True | True |
| 6 | 1.155e-14 | 1.155e-14 | 1 | UL19_36407+;UL29_58709-;UL30_64524+;UL5_12570+;UL52_109245+;UL54_114833- | UL19;UL29;UL30;UL5;UL52;UL54 | True | True |
| 7 | 1.095e-16 | 1.095e-16 | 1 | RL2_5080+;UL19_36407+;UL29_58709-;UL30_64524+;UL5_12570+;UL52_109245+;UL54_114833- | RL2;UL19;UL29;UL30;UL5;UL52;UL54 | True | True |


Read the curve, not any single number. Escape falls steeply from k=1 to k=2, then the
returns collapse once the set covers every sequenced isolate: beyond that point each
additional guide multiplies in only its own q_i, which is small but no longer
improving anything about coverage.

The curve is computed on the same empirical genome population throughout, so the
flattening is a real property of the isolate set, not a numerical artefact. It is also
the reason a "how many guides do you need" answer cannot be given as a single integer
without stating the threshold, which is what section 7.2 does.

### 7.0 Is the search actually finding the optimum?

| check | method | result |
|---|---|---|
| k=1, k=2 search | exhaustive over the whole allowed pool | exact by construction |
| k=3 search | brute force over the best 60 sites (60C3 = 34,220 sets) | brute force 1.075e-07 vs heuristic 1.075e-07 -- heuristic matches |

### 7.1 The resolution floor -- what this many genomes can certify

| set | k | P(escape), measured absences | P(escape), 95% upper bound on unobserved absence | ratio |
|---|---|---|---|---|
| Amrani lead pair | 2 | 9.720e-04 | 1.137e-02 | 11.7 |
| best k=1 | 1 | 4.754e-03 | 2.091e-02 | 4.399 |
| best k=2 | 2 | 2.260e-05 | 4.374e-04 | 19.35 |
| best k=3 | 3 | 1.075e-07 | 9.148e-06 | 85.13 |
| best k=4 | 4 | 5.109e-10 | 1.913e-07 | 374.5 |
| best k=5 | 5 | 2.429e-12 | 4.001e-09 | 1647 |
| best k=6 | 6 | 1.155e-14 | 8.368e-11 | 7246 |
| best k=7 | 7 | 5.490e-17 | 1.750e-12 | 3.188e+04 |
| best k=8 | 8 | 2.610e-19 | 3.660e-14 | 1.402e+05 |


This is the limit of what the data can support, and it is easy to miss.

A site absent in 0 of 183 sequenced isolates has a MEASURED absence frequency of
exactly zero, so the model assigns it no pre-existing-absence escape at all. But
0/183 is consistent with a true population absence frequency of anything up to
0.0162 (Clopper-Pearson 95% upper limit; the familiar rule of three). The
right-hand column re-runs the model with every site's absence replaced by that bound.

The consequence is structural rather than parametric: for well-chosen sets it is the
**sampling resolution of the genome corpus, not the NHEJ model, that sets the floor on
demonstrable escape probability**. With 183 genomes a single site cannot be
certified below ~1.6e-02, a pair below ~2.6e-04, and so
on. Any claim that a two-guide design achieves an escape probability below about
3e-04 is therefore an extrapolation beyond the evidence, no matter
how the repair term is parameterised -- and that includes claims made by this model.
This bound is one-sided and pessimistic on purpose: it assumes the unobserved absences
are independent across sites, which is the direction that maximises escape.

### 7.2 Minimum k for stated thresholds

| threshold on P(escape) per genome | minimum k | achieved P(escape) |
|---|---|---|
| 1.000e-03 | 2 | 2.260e-05 |
| 1.000e-06 | 3 | 1.075e-07 |
| 1.000e-03  (<1 escape genome in a population of 1000) | 2 | 2.260e-05 |
| 1.000e-05  (<1 escape genome in a population of 100000) | 3 | 1.075e-07 |
| 1.000e-07  (<1 escape genome in a population of 1e+07) | 4 | 5.109e-10 |

The same question asked of the evidence rather than of the point estimate. The right-hand column is the smallest k that clears the threshold even at the 95% upper bound on unobserved site absence, and it is the number a manuscript should quote.

| threshold on P(escape) per genome | minimum k, point estimate | minimum k the corpus can certify | gap |
|---|---|---|---|
| 1.000e-03 | 2 | 2 | 0 |
| 1.000e-06 | 3 | 4 | 1 |
| 1.000e-05 | 3 | 3 | 0 |
| 1.000e-07 | 4 | 5 | 1 |


The population-size rows are the ones with a clinical reading. A per-genome escape
probability only becomes an escape *event* when multiplied by the number of viral
genomes in the compartment. The latent HSV-1 genome load per host is an ASSUMPTION
here -- it is exposed as `--viral-population-size` and tabulated over four orders of
magnitude rather than being asserted -- but the structure is unambiguous: the required
per-genome threshold is 1/N, so every tenfold increase in viral load costs roughly one
extra order of magnitude of escape suppression, and the k needed grows with it.

### 7.3 Correlated failure, and what independence would have got wrong

| set | k | joint conservation (measured) | product of marginals | measured - product | P(escape) empirical | P(escape) independent | relative error of independence |
|---|---|---|---|---|---|---|---|
| Amrani lead pair | 2 | 0.8251 | 0.8285 | -0.003344 | 9.720e-04 | 4.269e-03 | 3.392 |
| best k=1 | 1 | 1 | 1 | 0 | 4.754e-03 | 4.754e-03 | 0 |
| best k=2 | 2 | 1 | 1 | 0 | 2.260e-05 | 2.260e-05 | 1.499e-16 |
| best k=3 | 3 | 1 | 1 | 0 | 1.075e-07 | 1.075e-07 | 0 |
| best k=4 | 4 | 1 | 1 | 0 | 5.109e-10 | 5.109e-10 | 0 |
| best k=5 | 5 | 1 | 1 | 0 | 2.429e-12 | 2.429e-12 | 1.663e-16 |
| best k=6 | 6 | 1 | 1 | 0 | 1.155e-14 | 1.155e-14 | -1.366e-16 |
| best k=7 | 7 | 1 | 1 | 0 | 5.490e-17 | 5.490e-17 | 3.368e-16 |
| best k=8 | 8 | 1 | 1 | 0 | 2.610e-19 | 2.610e-19 | 0 |


This is the check the model exists to make, and it has a sign that is easy to get
backwards, so it is worth spelling out.

`joint conservation (measured)` versus `product of marginals` is the exact, empirical
statement of how far site presence is from independent. When the measured joint is
BELOW the product -- which is the case for the Amrani lead pair, reproducing the
stage-6 finding -- the two guides fail in DIFFERENT isolates. That has two opposite
consequences and both are real:

* **Coverage is worse than the marginals suggest.** More isolates lose at least one
  of the two sites (0.825 of isolates keep both, against marginals of
  0.847 and 0.978).
  This is the efficacy cost, and it is invisible to a per-guide conservation
  table -- which is exactly the form in which their selection method reports itself.
* **Escape is *better* than an independence calculation would say.** Escape needs
  BOTH sites to fail in the same genome, and spreading the failures across different
  isolates makes that coincidence rarer. The `relative error of independence` column
  is positive for the lead pair: assuming independence would have OVER-stated its
  escape probability.

Neither effect can be recovered from marginal conservation percentages. Getting the
sign right requires the joint matrix, which is why the model never collapses it.

## 8. Sensitivity analysis


Every uncertain parameter is swept one at a time, holding the rest at baseline. Two
different things are reported and they behave very differently: the ABSOLUTE escape
probability, and the RANK ORDER of candidate sites.

| group | label | inframe_fraction | median_q_site | p_escape_amrani_pair | p_escape_best_k1 | p_escape_best_k2 | p_escape_best_k3 | p_escape_best_k4 | min_k_below_0.001 | min_k_below_1e-06 | spearman_q_vs_baseline |
|---|---|---|---|---|---|---|---|---|---|---|---|
| baseline | baseline | 0.1986 | 4.754e-03 | 9.720e-04 | 4.754e-03 | 2.260e-05 | 1.075e-07 | 5.109e-10 | 2 | 3 | 1 |
| tol_variable | tol_variable=0.1 | 0.1986 | 4.754e-03 | 9.720e-04 | 4.754e-03 | 2.260e-05 | 1.075e-07 | 5.109e-10 | 2 | 3 | 0.9964 |
| tol_variable | tol_variable=0.25 | 0.1986 | 4.754e-03 | 9.720e-04 | 4.754e-03 | 2.260e-05 | 1.075e-07 | 5.109e-10 | 2 | 3 | 0.9965 |
| tol_variable | tol_variable=0.5 | 0.1986 | 4.754e-03 | 9.720e-04 | 4.754e-03 | 2.260e-05 | 1.075e-07 | 5.109e-10 | 2 | 3 | 1 |
| tol_variable | tol_variable=0.75 | 0.1986 | 4.754e-03 | 9.720e-04 | 4.754e-03 | 2.260e-05 | 1.075e-07 | 5.109e-10 | 2 | 3 | 1 |
| tol_variable | tol_variable=1 | 0.1986 | 4.754e-03 | 9.720e-04 | 4.754e-03 | 2.260e-05 | 1.075e-07 | 5.109e-10 | 2 | 3 | 1 |
| tol_invariant | tol_invariant=0.001 | 0.1986 | 9.278e-05 | 1.827e-05 | 9.278e-05 | 8.608e-09 | 7.986e-13 | 7.410e-17 | 1 | 2 | 0.9543 |
| tol_invariant | tol_invariant=0.01 | 0.1986 | 9.319e-04 | 1.847e-04 | 9.319e-04 | 8.685e-07 | 8.094e-10 | 7.543e-13 | 1 | 2 | 0.9824 |
| tol_invariant | tol_invariant=0.05 | 0.1986 | 4.754e-03 | 9.720e-04 | 4.754e-03 | 2.260e-05 | 1.075e-07 | 5.109e-10 | 2 | 3 | 1 |
| tol_invariant | tol_invariant=0.15 | 0.1986 | 1.503e-02 | 3.323e-03 | 1.503e-02 | 2.260e-04 | 3.397e-06 | 5.107e-08 | 2 | 4 | 0.994 |
| tol_invariant | tol_invariant=0.3 | 0.1986 | 3.277e-02 | 8.168e-03 | 3.277e-02 | 1.074e-03 | 3.518e-05 | 1.153e-06 | 3 | 5 | 0.9839 |
| tol_invariant | tol_invariant=0.5 | 0.1986 | 6.227e-02 | 1.836e-02 | 6.227e-02 | 3.877e-03 | 2.414e-04 | 1.503e-05 | 3 | 5 | 0.4406 |
| min_variant_strains | min_variant_strains=1 | 0.1986 | 4.754e-03 | 9.720e-04 | 4.754e-03 | 2.260e-05 | 1.075e-07 | 5.109e-10 | 2 | 3 | 1 |
| min_variant_strains | min_variant_strains=2 | 0.1986 | 4.754e-03 | 9.720e-04 | 4.754e-03 | 2.260e-05 | 1.075e-07 | 5.109e-10 | 2 | 3 | 0.7124 |
| min_variant_strains | min_variant_strains=3 | 0.1986 | 4.754e-03 | 9.720e-04 | 4.754e-03 | 2.260e-05 | 1.075e-07 | 5.109e-10 | 2 | 3 | 0.6383 |
| min_variant_strains | min_variant_strains=5 | 0.1986 | 4.754e-03 | 9.720e-04 | 4.754e-03 | 2.260e-05 | 1.075e-07 | 5.109e-10 | 2 | 3 | 0.5516 |
| indel_spectrum | spectrum=default | 0.1986 | 4.754e-03 | 9.720e-04 | 4.754e-03 | 2.260e-05 | 1.075e-07 | 5.109e-10 | 2 | 3 | 1 |
| indel_spectrum | spectrum=deletion-heavy | 0.2676 | 3.489e-03 | 7.062e-04 | 3.489e-03 | 1.218e-05 | 4.249e-08 | 1.482e-10 | 2 | 3 | 0.9976 |
| indel_spectrum | spectrum=insertion-heavy | 0.1147 | 2.903e-03 | 5.849e-04 | 2.903e-03 | 8.430e-06 | 2.448e-08 | 7.107e-11 | 2 | 3 | 1 |
| indel_spectrum | spectrum=short-indels | 0.1397 | 5.640e-03 | 1.161e-03 | 5.640e-03 | 3.181e-05 | 1.794e-07 | 1.012e-09 | 2 | 3 | 0.9543 |
| indel_spectrum | spectrum=uniform-1to20 | 0.3 | 2.632e-03 | 5.290e-04 | 2.632e-03 | 6.925e-06 | 1.822e-08 | 4.796e-11 | 2 | 3 | 0.9543 |
| indel_spectrum | spectrum=neuronal_nhej | 0.09034 | 1.824e-03 | 3.642e-04 | 1.824e-03 | 3.326e-06 | 6.066e-09 | 1.106e-11 | 2 | 3 | 1 |
| inframe_fraction | inframe_fraction=0.1 | 0.1 | 2.394e-03 | 4.803e-04 | 2.394e-03 | 5.732e-06 | 1.372e-08 | 3.286e-11 | 2 | 3 | 1 |
| inframe_fraction | inframe_fraction=0.2 | 0.2 | 4.788e-03 | 9.792e-04 | 4.788e-03 | 2.293e-05 | 1.098e-07 | 5.257e-10 | 2 | 3 | 1 |
| inframe_fraction | inframe_fraction=0.25 | 0.25 | 5.986e-03 | 1.236e-03 | 5.986e-03 | 3.583e-05 | 2.144e-07 | 1.284e-09 | 2 | 3 | 1 |
| inframe_fraction | inframe_fraction=0.33 | 0.33 | 7.901e-03 | 1.656e-03 | 7.901e-03 | 6.242e-05 | 4.932e-07 | 3.897e-09 | 2 | 3 | 1 |
| inframe_fraction | inframe_fraction=0.5 | 0.5 | 1.197e-02 | 2.587e-03 | 1.197e-02 | 1.433e-04 | 1.716e-06 | 2.054e-08 | 2 | 4 | 1 |
| repeat_model | repeat_model=redundant | 0.1986 | 4.754e-03 | 9.720e-04 | 4.754e-03 | 2.260e-05 | 1.075e-07 | 5.109e-10 | 2 | 3 | 1 |
| repeat_model | repeat_model=single | 0.1986 | 4.754e-03 | 8.500e-04 | 4.754e-03 | 2.260e-05 | 1.075e-07 | 5.109e-10 | 2 | 3 | 0.9209 |
| repeat_model | repeat_model=all-copies | 0.1986 | 4.754e-03 | 7.280e-04 | 2.260e-05 | 5.109e-10 | 1.155e-14 | 5.490e-17 | 1 | 2 | 0.525 |
| frameshift_viability_RL2 | RL2_frameshift_viability=0 | 0.1986 | 4.754e-03 | 9.720e-04 | 4.754e-03 | 2.260e-05 | 1.075e-07 | 5.109e-10 | 2 | 3 | 1 |
| frameshift_viability_RL2 | RL2_frameshift_viability=0.05 | 0.1986 | 4.754e-03 | 2.987e-03 | 4.754e-03 | 2.260e-05 | 1.075e-07 | 5.109e-10 | 2 | 3 | 0.9917 |
| frameshift_viability_RL2 | RL2_frameshift_viability=0.25 | 0.1986 | 4.754e-03 | 1.022e-02 | 4.754e-03 | 2.260e-05 | 1.075e-07 | 5.109e-10 | 2 | 3 | 0.9917 |
| frameshift_viability_RL2 | RL2_frameshift_viability=1 | 0.1986 | 4.754e-03 | 2.554e-02 | 4.754e-03 | 2.260e-05 | 1.075e-07 | 5.109e-10 | 2 | 3 | 0.9917 |
| p_disrupt | p_disrupt=0.8 | 0.1986 | 3.803e-03 | 7.717e-04 | 3.803e-03 | 1.447e-05 | 5.502e-08 | 2.093e-10 | 2 | 3 | 1 |
| p_disrupt | p_disrupt=0.9 | 0.1986 | 4.279e-03 | 8.715e-04 | 4.279e-03 | 1.831e-05 | 7.834e-08 | 3.352e-10 | 2 | 3 | 1 |
| p_disrupt | p_disrupt=1 | 0.1986 | 4.754e-03 | 9.720e-04 | 4.754e-03 | 2.260e-05 | 1.075e-07 | 5.109e-10 | 2 | 3 | 1 |

### 8.1 Does the ranking survive?

| group | settings | min P(escape) Amrani pair | max P(escape) Amrani pair | orders of magnitude spanned | min Spearman rho of site q | min top-20 site overlap | settings whose best k=2/3/4 sets all match baseline |
|---|---|---|---|---|---|---|---|
| baseline | 1 | 9.720e-04 | 9.720e-04 | 0 | 1 | 1 | 1 |
| frameshift_viability_RL2 | 4 | 9.720e-04 | 2.554e-02 | 1.42 | 0.9917 | 1 | 4 |
| indel_spectrum | 6 | 3.642e-04 | 1.161e-03 | 0.5036 | 0.9543 | 0.7 | 4 |
| inframe_fraction | 5 | 4.803e-04 | 2.587e-03 | 0.7313 | 1 | 1 | 5 |
| min_variant_strains | 4 | 9.720e-04 | 9.720e-04 | 9.643e-17 | 0.5516 | 0.45 | 1 |
| p_disrupt | 3 | 7.717e-04 | 9.720e-04 | 0.1002 | 1 | 1 | 3 |
| repeat_model | 3 | 7.280e-04 | 9.720e-04 | 0.1255 | 0.525 | 0 | 1 |
| tol_invariant | 6 | 1.827e-05 | 1.836e-02 | 3.002 | 0.4406 | 0.35 | 3 |
| tol_variable | 5 | 9.720e-04 | 9.720e-04 | 9.876e-06 | 0.9964 | 0.9 | 5 |


**This is the honest headline of the sensitivity analysis, and it is not a clean
"the ranking is robust" result.** It has four parts and they should be read together.

**1. Absolute probabilities are not usable as measurements.** Across the swept space
the escape probability of a fixed guide set moves by up to 3.1 orders of
magnitude -- almost all of it driven by `tol_invariant`, the one quantity the data
genuinely cannot pin down (absence of observed variation over 183 isolates bounds
the tolerated-variant frequency but does not measure the tolerance). No number in this
report should be quoted as an escape rate.

**2. The chosen guide SETS are stable.** 27 of 37 parameter
settings return exactly the same best k=2, k=3 and k=4 sets as the baseline. Where
the design recommendation is concerned, the model is doing something reproducible.

**3. The per-site ORDERING is only partly stable, and the exceptions are informative
rather than noise.** The rank correlation of per-site q falls as low as
0.441 and the top-20 overlap as low as 0%, worst in
`repeat_model`. Two mechanisms cause this and neither is a defect:

* Most sites sit at exactly the same q (the invariant-codon value), so the ordering is
  dominated by ties and any reweighting reshuffles them. That is why the top-20
  overlap statistic is reported next to Spearman -- a large rank movement among tied
  sites is not a change of recommendation.
* The `repeat_model` and `frameshift_viability_RL2` settings genuinely reorder ICP0
  against everything else, because they are assumptions ABOUT ICP0 specifically. The
  correct reading is not that the model is unstable but that **the relative merit of
  an ICP0 guide is not determined by the available data** -- it depends on whether
  cutting both repeat copies must be escaped twice, and on whether an ICP0 frameshift
  is lethal. Both are experiments, not parameters.

**4. The one swept spectrum that is not an assumption moves the absolute scale and
leaves the selection untouched.** `spectrum=neuronal_nhej` is the pooled CRISPResso2
net-length histogram deposited as source data with Ramadoss et al. 2025 Fig. 1d --
118,722 reads over six replicates of SpCas9 RNP editing in post-mitotic human
iPSC-derived neurons, next to the genetically identical dividing iPSCs. Its shape is
MEASURED; its applicability to SaCas9 in a latently infected trigeminal ganglion is
not, which is why it is a scenario here and not the default. Three things follow.

* Its in-frame fraction is **0.090**, against
  0.199 for the default -- and BELOW the 0.10-0.50 band
  the `inframe_fraction` group sweeps. The measured post-mitotic spectrum was outside
  the range this model previously explored.
* Escape therefore FALLS, and it falls harder the more guides you use, because q
  multiplies across sites: the published pair goes
  9.720e-04 -> 3.642e-04
  (2.67x lower),
  and the best sets fall 3x at k=1, 7x at k=2, 18x at k=3, 46x at k=4. The direction is worth being explicit about
  because it is easy to get backwards: "a narrower distribution of smaller indels"
  reads like a gentler outcome, but +-1 and +-2 nt are 77% of the measured
  neuronal indels and every one of them is a frameshift. A smaller indel in an
  essential gene is a MORE lethal indel. If neurons repair the way Ramadoss et al.
  measured, this model's default is conservative -- it over-states escape.
* The selection does not move. Spearman rho of per-site q against baseline is
  0.999996, top-20 overlap
  100%, the best sets at every k
  match the baseline exactly,
  and the minimum-k answers at both thresholds are
  unchanged. This is the
  same pattern as the rest of the sweep, now demonstrated against measured rather than
  posited numbers: the absolute scale is not a measurement, and the ranking is what the
  model is for.

The isogenic dividing arm of the same experiment is available as `ipsc_dividing` and is
deliberately NOT swept -- it is the control that shows the shift is a cell-type effect
and not a batch effect (in-frame fraction 0.289 in the iPSCs against
0.090 in the neurons -- same experiment, same guide,
same Cas9 dose, cells differing only in whether they had been differentiated).
Sweeping both arms would double-count one experiment as two independent constraints on
the same parameter.

**The one result that flips a design conclusion is worth stating on its own.** Under
`repeat_model=all-copies` -- the regime in which a guide cutting both ICP0 repeat
copies requires an independent viable escape at each -- a SINGLE ICP0 guide reaches
the escape probability that two guides elsewhere reach (2.260e-05 at k=1). That
is Amrani et al.'s "ICP0 is duplicated, so two guides give three DSBs" argument,
quantified, and in that regime it is correct. It is also the regime our default does
NOT assume, because a redundant gene can absorb a broken copy. Which of the two holds
is the single highest-value experiment this analysis points to.

Stated as a rule for using this model: **use it to choose between guide sets, and do
not quote its absolute escape probabilities as if they were measurements.** The
minimum-k answers inherit the absolute-scale uncertainty and should be read together
with the resolution floor in section 7.1.

## 9. Closed form versus simulation

| set | k | closed form | Monte Carlo | MC standard error | abs(difference) / SE | resolvable at this many draws |
|---|---|---|---|---|---|---|
| Amrani lead pair | 2 | 9.720e-04 | 9.975e-04 | 2.232e-05 | 1.143 | yes |
| best k=1 | 1 | 4.754e-03 | 4.686e-03 | 4.829e-05 | 1.423 | yes |
| best k=2 | 2 | 2.260e-05 | 1.950e-05 | 3.122e-06 | 0.9937 | yes |
| best k=3 | 3 | 1.075e-07 | 5.000e-07 | 5.000e-07 | 0.7851 | no (below MC resolution) |
| best k=4 | 4 | 5.109e-10 | 0 | 7.071e-10 | 0.7225 | no (below MC resolution) |


The Monte Carlo draws a genome uniformly from the empirical population and then a
Bernoulli repair outcome per present site, with seed 20240814 and
2,000,000 draws. It is a simulation of the same generative process the closed
form integrates, so agreement is a correctness check on the implementation, not
evidence about biology. Sets whose closed-form probability is far below the Monte
Carlo resolution (~2.5e-06) are marked as such rather than being
reported as spuriously equal to zero.

## 10. The delivery term, which no number of guides can fix

| theta (genomes never exposed to nuclease) | P(persist), Amrani pair | P(persist), best k=2 | P(persist), best k=8 |
|---|---|---|---|
| 0 | 9.720e-04 | 2.260e-05 | 2.610e-19 |
| 0.5 | 5.005e-01 | 5.000e-01 | 5.000e-01 |
| 0.9 | 9.001e-01 | 9.000e-01 | 9.000e-01 |
| 0.99 | 9.900e-01 | 9.900e-01 | 9.900e-01 |


theta is the fraction of viral genomes that never meet an active nuclease at all --
untransduced neurons, silenced AAV episomes, genomes in an inaccessible chromatin
state. It is k-independent by construction, so it adds a floor that no amount of
multiplexing can lower.

This matters because the only in vivo on-target measurement in Amrani et al. is a mean
indel frequency of ~0.78-1.65% in trigeminal ganglia at day 40, which -- read
literally as an editing rate -- corresponds to theta near 0.99. At that theta the
table above is flat: every guide set, at every k, gives essentially the same
persistence, because persistence is delivery-limited and not escape-limited.

That is not an argument against multiplexing. It is an argument that the two questions
must be reported separately, which no published HSV-1 CRISPR work currently does.
Escape probability governs whether the edited fraction can regenerate a resistant
population under selection; theta governs how large the unedited fraction is. This
model is about the first, and is silent about the second.

## 11. What this model cannot tell you

1. **It cannot tell you an absolute escape rate.** The NHEJ indel spectrum and the two tolerance constants are assumptions, and section 8 shows the absolute answer moves by orders of magnitude across their plausible range. Every absolute number in this report is conditional on the parameter table in section 2.

2. **It cannot substitute cross-strain substitution tolerance for indel tolerance.** A codon that tolerates a substitution among clinical isolates may still not tolerate deletion. The profile is an upper bound on tolerance, so it over-states escape -- the safe direction for a guide-count argument, but a real limitation on any claim about a specific site.

3. **It cannot see selection, fitness or kinetics.** Escape here is a per-genome probability of producing a viable uncuttable genome, not a rate of outgrowth. A viable-but-crippled escape variant counts the same as a fully fit one. Real escape dynamics depend on replication rate, bottlenecks and immune pressure, none of which are modelled.

4. **It cannot model multi-cut excision.** Two simultaneous DSBs frequently excise the intervening fragment. For guides in different genes that is certainly lethal to the virus and would push escape BELOW our estimates; for the two ICP0 repeat copies it would excise the entire unique long region. Excision is therefore an unmodelled escape-suppressing mechanism, and the `all-copies` repeat model in section 8 is only a crude proxy for it.

5. **It cannot resolve whether ICP0 frameshifts are lethal.** ICP0 is a virulence and reactivation factor, not a strictly essential gene. The baseline treats an ICP0 frameshift as lethal because that is the assumption most favourable to the published two-guide design; the sweep shows what happens if it is not, and the honest answer is that this is a question for an experiment, not for a model.

6. **It cannot see off-target activity, delivery, packaging or expression.** Nothing in THIS module looks at the human genome. Stage 8 (`src/offtarget.py`) does, but only for the stage-6 SaCas9 benchmark guide set in RL2 and UL54 -- not for the UL19/UL29/UL30/UL5/UL52 sites this model puts at the top of its k-curve, none of which has been screened. A set this model calls optimal may still be undeliverable, unsynthesisable or unsafe, and results/recommendation.md is where the three axes are reconciled for the guides that HAVE been screened on all of them.

7. **The denominator is 183 complete genomes, not the circulating population.** Stage 5 puts the effective sample size nearer 132 at 99.9% identity, the set mixes clinical isolates with laboratory strains, and rare variants below ~1/183 frequency are invisible. A site called perfectly conserved here has an exact binomial 95% lower bound of 0.984, not 1.0, so P(escape) values below roughly 1.6e-02 from the presence term alone are extrapolations beyond the resolution of the data. Section 7.1 quantifies this properly; it is the single most binding limitation in the whole analysis.

8. **It cannot model within-host quasispecies structure.** The isolate population is used as a proxy for the diversity a single patient's virus presents. Within-host diversity is lower; between-host diversity is what is measured here. That makes the presence term pessimistic for a single patient and appropriate for a product intended for a population.

9. **It models repair pathway choice only as a swept scenario, not as a mechanism.** Microhomology-mediated end joining produces predictable, often larger deletions with a different in-frame fraction, and its use varies by cell type. Section 8 now sweeps a spectrum measured in post-mitotic human neurons (`neuronal_nhej`) alongside the assumed default, which bounds the consequence of that pathway shift for THIS model. It does not make the model pathway-aware: there is still one spectrum per run, applied to every site identically, with no dependence on the local microhomology content of the target -- which is the very thing the source paper shows drives the guide-to-guide differences.

## 12. Reproducibility


Deterministic given the seed. The closed-form path uses no randomness at all; the
Monte Carlo uses `numpy.random.default_rng(20240814)` and is reported alongside the
closed form so that any disagreement is visible. Guide-set search is deterministic:
exhaustive for k<=2, and for larger k a multi-start greedy plus exhaustive local swap
with fixed scan order and index-based tie-breaks.

Command line: `python src/escape.py --skip-fetch --force --escape`
Reference: NC_001806.2. Genomes: 183. Nuclease: sacas9 (20 nt spacer + NNGRRT PAM).
Indel spectrum: parametric-default (ASSUMPTION), in-frame fraction
0.1986, mean deletion length
6.00 nt.
Anchor parameters for the tolerance profile: k=25, step=6,
max unanchored bridge=60 nt, minimum resolved fraction 0.95.

