# Ramadoss et al. 2025 — CRISPR repair outcomes in nondividing human cells — extraction notes

**Citation.** Ramadoss GN, Namaganda SJ, Kumar MM, Hamilton JR, Sharma R, Chow KG, Workley LA,
Macklin BL, Sun M, Ha AS, Liu JC, Fellmann C, Watry HL, Dierks PH, Bose RS, Jin J, Perez BS,
Sandoval Espinoza CR, Matia MP, Lu SH, Judge LM, Shy BR, Nussenzweig A, Adamson B, Murthy N,
Doudna JA, Kampmann M, Conklin BR. "Characterizing and controlling CRISPR repair outcomes in
nondividing human cells." *Nature Communications* 2025;16(1):9883.
doi:10.1038/s41467-025-66058-3. PMID 41249169 · PMCID PMC12623481.

**URLs.**
- Full text (open access): https://europepmc.org/articles/PMC12623481
- Machine-readable full text used here:
  https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12623481/fullTextXML
- Source data (FigShare, cited in the paper's own Data Availability statement):
  FigShare article 30366298 (doi:10.6084/m9.figshare.30366298.v1, published 2025-10-15),
  file 58769638, `Source Data.zip`, 2,273,223,940 bytes.
  https://figshare.com/articles/dataset/Source_Data_for_Characterizing_and_controlling_CRISPR_repair_outcomes_in_nondividing_human_cells_Ramadoss_et_al_Nature_Communications/30366298

---

## 1. What the PAPER TEXT reports quantitatively about indel length — very little

This matters, so it is stated first. The narrative text of the paper makes the
indel-length claim **qualitatively**, not numerically:

- "while B2Mg1-edited iPSCs displayed a broad range of indels, neurons exhibited a much
  narrower distribution of outcomes (Fig. 1d)"
- "In iPSCs, the most prevalent indel outcomes were larger deletions typically associated
  with MMEJ ... In neurons, the most prevalent outcomes were those usually attributed to
  NHEJ: small indels associated with NHEJ processing, and unedited outcomes"
- "for every sgRNA we tested, the ratio of insertions to deletions was significantly higher
  in neurons than iPSCs (Supplementary Fig. 6)"

The only length-resolved numbers stated in running text are about **perturbations**, not
about the untreated neuronal baseline: RNR (RRM2) inhibition with 3AP gives "a ~50% increase
in total indels" and "a ~3-fold increase in single-base deletions"; siRRM1 increases the
frequency of single-base deletions "by ~75%". None of those is a spectrum.

**The indel-length distribution itself is presented only as figure panels (Fig. 1d,
Supplementary Fig. 6).** It cannot be read numerically off the figures. Had that been the
end of it, no measured spectrum would have been extractable, and the scenario below would
have had to be an explicitly *derived-from-qualitative-description* construction rather than
a measured one.

## 2. What WAS extracted, and how — the deposited source data

The paper deposits per-figure raw CRISPResso2 output. `Source Data/Figure 1d.zip` inside the
FigShare archive contains twelve CRISPResso2 run directories:

    Figure 1d/Neurons/CRISPResso_on_12194dedvr{1..6}/Indel_histogram.txt
    Figure 1d/iPSCs/CRISPResso_on_1122-ipscs-04d-edv-r{1..6}/Indel_histogram.txt

`Indel_histogram.txt` is CRISPResso2's read count per **net length change** of the aligned
read relative to reference (negative = net deletion, positive = net insertion, 0 = no net
length change). That is exactly the quantity this repository's `IndelSpectrum` is defined
over, so no re-interpretation was needed. One consequence of the "net" definition is worth
stating: a read carrying a compensating insertion and deletion of equal size lands in the
0 bin, so the reported no-net-change mass is a slight over-count of true indel-free repair.
The model conditions the 0 bin away, so this affects reporting only.

Retrieval was by HTTP range request against the FigShare S3 object: the ZIP end-of-central-
directory and central directory were read from the tail, the two member entries located by
offset, and only those bytes fetched. The 2.27 GB archive was not downloaded in full.

Pooled counts (6 replicates per arm, summed) are archived verbatim in
`refs/ramadoss2025_fig1d_indel_histogram.tsv` (SHA-256 `d263b5698e9592e4584d072fb55c286ad655933aa82fe507b0032c46bcf0e1a7`),
columns `length`, `neuron_reads`, `ipsc_reads`, covering net lengths -71 .. +112.
Totals: 118,722 neuron reads; 75,183 iPSC reads.

## 3. Summary statistics computed from those counts

Both arms are the SAME experiment: equal doses of VLP-delivered SpCas9 RNP with sgRNA
**B2Mg1** targeting the human **B2M** locus, in human iPSC-derived neurons and in the
**genetically identical** parental iPSCs, harvested 4 days post-transduction. The pairing is
what makes the contrast a cell-type effect rather than a batch effect.

| statistic (conditioned on net length != 0 unless noted) | neurons | iPSCs |
|---|---|---|
| pooled reads | 118,722 | 75,183 |
| reads with an indel (net length != 0) | 54,881 | 67,537 |
| P(no net length change) — reported only, conditioned away by the model | 0.5377 | 0.1017 |
| **in-frame fraction (net length divisible by 3)** | **0.0903** | **0.2893** |
| insertions as a share of indels | 0.4030 | 0.1606 |
| insertion : deletion ratio | 0.675 | 0.191 |
| +1 nt as a share of all indels | 0.2940 | 0.0533 |
| +1 nt as a share of insertions | 0.7296 | 0.3318 |
| mean deletion length (nt) | 6.60 | 12.72 |
| median absolute indel length (nt) | 1 | 10 |
| mean absolute indel length (nt) | 4.63 | 11.85 |
| share of indels with absolute length <= 3 nt | 0.807 | 0.231 |

Per-replicate spread (n = 6 each), to show the pooled numbers are not driven by one well:

- neuron in-frame fraction: 0.0795, 0.0876, 0.0888, 0.0904, 0.0918, 0.1109 (mean 0.0915, SD 0.0104)
- iPSC in-frame fraction: 0.2819, 0.2861, 0.2882, 0.2898, 0.2929, 0.2965 (mean 0.2892, SD 0.0051)
- neuron insertion:deletion: 0.616, 0.625, 0.666, 0.683, 0.697, 0.749
- iPSC insertion:deletion: 0.180, 0.184, 0.189, 0.195, 0.197, 0.206

Every qualitative claim in the paper text is reproduced by these numbers: the neuronal
distribution is narrower (median 1 nt vs 10 nt; 81% of indels <= 3 nt vs 23%), the
insertion:deletion ratio is 3.5x higher in neurons, and the large MMEJ-like deletions that
dominate the iPSC arm are a minor component in neurons.

The four most probable neuronal outcomes are +1 (0.294), -1 (0.280), +2 (0.099) and -2
(0.099) — together 77% of all indels, and **all four are frameshifts**. That is the mechanism
behind the single summary statistic the escape model is sensitive to: a narrower spectrum of
smaller indels is a MORE frameshifting spectrum, not a less disruptive one.

## 4. What could NOT be extracted

- **No spectrum for SaCas9.** The whole dataset is SpCas9. The therapy under audit uses
  SaCas9. Both are blunt cutters at an analogous position relative to the PAM, so a similar
  outcome distribution is plausible, but that is an assumption and not a measurement.
- **No spectrum from trigeminal-ganglion neurons, and none from any animal.** These are
  cultured human iPSC-derived neurons. The therapy targets mouse and rabbit — and eventually
  human — trigeminal ganglia in vivo, in latently infected sensory neurons.
- **No locus generality.** The pooled histogram above is one sgRNA at one locus (B2Mg1 at
  B2M). The paper's own central caveat is that "each sgRNA had a different intrinsic
  distribution of available indel types", and that RNR inhibition helps for three sgRNAs but
  not for the insertion-biased HSPB1g2. Supplementary Fig. 6 covers the additional sgRNAs;
  its source data is a separate 123 MB member that was not extracted, so the sgRNA-to-sgRNA
  spread of the neuronal in-frame fraction is NOT characterised here. A single-locus estimate
  is used, and it should be read as one draw from that unmeasured spread.
- **No MMEJ/NHEJ pathway fractions as numbers.** The paper assigns outcomes to pathways by
  qualitative signature (large microhomology-flanked deletions = MMEJ-like; small indels and
  indel-free repair = NHEJ-like). It reports no numeric MMEJ:NHEJ split, and none is invented
  here.
- **Timing is described but not as a single rate constant.** Neurons: indels keep rising for
  at least 16 days post-RNP and up to ~2 weeks post-transduction; DSB repair signal at the
  target decays by only ~50% one week after RNP delivery; Cas9 protein remains detectable up
  to 30 days. iPSCs plateau within a few days; the dividing-cell repair half-life they quote
  from the literature is 1-10 h. The 4-day harvest used for Fig. 1d is therefore taken BEFORE
  the neuronal plateau, which is a caveat on the neuronal no-net-change mass (0.538) — it is
  an upper bound on the eventual unedited fraction. It does not bias the *conditional*
  spectrum in any direction the data can resolve, and the model conditions on an indel having
  occurred.
- **Repair factors upregulated.** Named and validated: RRM2/RRM1 (ribonucleotide reductase —
  canonically inactive in nondividing cells, yet inhibiting it shifts neuronal outcomes from
  insertions toward deletions and triples 1 bp deletions), POLL (Pol lambda), XRCC5 (Ku80).
  These are perturbations, not baseline parameters, and none is used in the model.

## 5. How this is used in this repository

`src/escape.py` exposes the pooled neuronal column as the built-in spectrum
**`neuronal_nhej`**, added to the sensitivity sweep only. **The default spectrum is
unchanged and the default pipeline outputs are byte-identical.** The scenario's provenance
string is MEASURED with respect to Ramadoss et al.'s neurons and ASSUMPTION with respect to
this therapy's cells and nuclease; both halves are printed in the report's parameter table.
