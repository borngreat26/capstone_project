# Weekly Tasks: Results and Takeaways

Findings from DSC 180A Weeks 1–6. The code and plots are in the weekly notebook
([public copy](notebooks/weekly_tasks_public.ipynb)). All associations come from additive models with no covariates,
in people of European ancestry and one cell type. They are associations, not proof of a causal effect.

## Week 1: Reading the Data

The chromosome 1 reference files contain **98,642 SNPs and 489 people**. `.bed` stores genotypes, `.bim` describes the
variants, and `.fam` describes the people.

## Week 2: Population Genetics

| Scope | Eligible SNPs | Rare (< 1%) | Low-frequency (1% < MAF < 5%) | Common (> 5%) |
|---|---:|---:|---:|---:|
| Chromosome 1 | 98,642 | 1,235 (1.2520%) | 9,300 (9.4280%) | 88,107 (89.3200%) |
| Autosomes 1–22 | 1,190,321 | 13,385 (1.1245%) | 106,159 (8.9185%) | 1,070,777 (89.9570%) |

Percentages describe SNP loci, not allele copies. Both scopes have zero SNPs exactly at 0.01 or 0.05, zero excluded or
monomorphic loci, and zero missing genotype calls.

There are **489 genotyped people, 462 expression samples, and 344 people in both**, with 145 genotype-only and 118
expression-only people. There are no duplicate IDs, and sample order is consistent across chromosomes.

**Takeaways**

- Most SNPs in this panel are common, but the panel's own SNP selection affects that share. These are EUR sample
  frequencies, not universal ones.
- An eQTL regression must compare each person's genotype with that same person's expression, so the files are aligned by
  ID, never by position. `data/week2_results/matched_sample_ids.tsv` records the pairing.
- Rare alleles have fewer carriers, so tests of them have less power at the same sample size and effect size.

## Week 3: cis-eQTL Mapping for PARK7

| Quantity | Result |
|---|---|
| Gene | PARK7, ENSG00000116288.8 |
| Gene body | GRCh37 chr1:8,014,351–8,045,565; plus strand |
| Cis-window | GRCh37 chr1:7,514,351–8,545,565; inclusive |
| Matched people | 344 |
| Candidate / retained SNPs | 437 / 437 |
| Missing expression / genotype calls | 0 / 0 |
| Constant SNPs excluded | 0 |
| Nominal p < 0.05 | 44 SNPs |
| BH FDR < 0.05 / Bonferroni < 0.05 | 0 / 0 SNPs |
| Leading SNP | rs11805932, chr1:7,521,479; counted allele G |
| Leading beta / standard error | −19.607803 / 5.905207 |
| Leading p / BH-adjusted p | 0.000995760 / 0.276679 |
| Leading Bonferroni-adjusted p | 0.435147 |
| Leading genotype group sizes | 0 G: 1; 1 G: 41; 2 G: 302 |

**Takeaways**

- **No SNP passes either multiple-testing correction,** so this analysis does not establish a significant cis-eQTL for
  PARK7. It also does not prove PARK7 has no genetic regulation: sample size, tissue, SNP panel, and the unadjusted model
  all limit what it can show.
- The strongest nominal SNP predicts about **19.61 fewer expression units per additional G allele**. Only one person has
  zero G alleles, so that slope needs caution.
- A positive beta means higher expression per copy of the named counted allele. Counting the other allele flips the sign
  and leaves the p-value unchanged. A p-value is not the probability that a SNP is causal.
- The local annotation lists PARK7 at 7,954,290–7,985,505, which is GRCh38. Our genotypes are GRCh37, so the window uses
  the GENCODE v12 GRCh37 record instead. PARK7 was chosen before looking at its association results.
- [GTEx](https://pmc.ncbi.nlm.nih.gov/articles/PMC7737656/) studies these genetic effects across human tissues, where
  they can be shared or differ. This one-gene exercise shows the same ideas on a smaller scale.

## Week 4: Chromosome 1 eQTLs and Genetic Prediction

### Mapping

All 2,479 expressed chromosome 1 genes have exact-version GRCh37 annotation. We analyzed **2,443 genes** with
**1,004,372 SNP–gene tests**. The other 36 genes had no eligible cis-SNPs.

| Mapping result | SNP–gene pairs | Genes with at least one passing pair |
|---|---:|---:|
| Nominal p < 0.05 | 86,234 | 2,368 |
| BH across all chromosome 1 pairs, < 0.05 | 10,621 | 746 |
| Bonferroni across all chromosome 1 pairs, < 0.05 | 2,446 | 183 |

8.6% of tests are nominally significant, but only 0.24% pass Bonferroni. Most genes have a nominal hit, so one alone says
little about a gene. The median number of SNPs tested per gene is **422** (IQR **335–497**), with a mean of **411.12** and
a range of **1–1,394**. These count the reference panel's markers, not all variants in the region.

### Prediction

**Gene selected: FLVCR1-AS1 (ENSG00000198468.2)**, a long non-coding RNA at GRCh37 chr1:213,025,450–213,031,430. A scan in
the 275 training people alone (1,004,175 tests; 143 genes with a Bonferroni-significant SNP) gave it the smallest p-value:
**p = 2.92 × 10⁻⁸⁶**, Bonferroni-adjusted **2.93 × 10⁻⁸⁰**. It stands well above every other gene, and no tie-break was
needed. The gap near 122–144 Mb in the gene-selection figure is the centromere region, where no analyzed genes begin.
Checked only afterward, it is also among the 183 full-sample Bonferroni genes.

**SNPs:** the cis-window chr1:212,525,450–213,531,430 contains 417 SNPs. **414** vary in the training group and were
tested, and clumping kept **47 index SNPs**. The signal peaks at the gene and fades within about 150 kb. Median training
expression rises from about 1 to 5 to 8 with zero, one, and two copies of the lead SNP's T allele.

**Threshold:** p ≤ 0.00001 had the highest validation predictive R².

| Inclusive p cutoff | SNPs | Validation Pearson r | Validation Pearson r² | Validation predictive R² |
|---|---:|---:|---:|---:|
| 1 | 47 | 0.676 | 0.457 | 0.454 |
| 0.05 | 14 | 0.711 | 0.505 | 0.502 |
| 0.01 | 10 | 0.724 | 0.524 | 0.523 |
| 0.001 | 7 | 0.766 | 0.586 | 0.582 |
| **0.00001 (selected)** | **4** | **0.783** | **0.613** | **0.605** |
| Training-mean baseline | 0 | Undefined | Undefined | −0.018 |

Accuracy improved as the cutoff tightened, which suggests that weakly associated SNPs mostly added noise. The model
under-predicts the people with the highest validation expression.

**Final score: four SNPs,** including the lead SNP, which lies inside the gene body.

| SNP | Effect allele | Training beta | Training p |
|---|---|---:|---:|
| rs2279692 | T | 3.830 | 2.92 × 10⁻⁸⁶ |
| rs2492782 | G | 1.300 | 7.40 × 10⁻⁷ |
| rs11119973 | C | −1.753 | 1.10 × 10⁻⁶ |
| rs1317482 | C | 1.435 | 4.12 × 10⁻⁶ |

The frozen equation, fitted on training people only, is **predicted expression = 1.945368 + 0.673416 × genetic score**.

| Original held-out test group (n = 35) | Result |
|---|---:|
| Final predictive R² | **0.740** |
| Pearson r | 0.874 |
| Squared Pearson r | 0.764 |
| Training-mean baseline predictive R² | −0.142 |

**The genetic model beats the baseline in both groups:** 0.605 versus −0.018 in validation, and 0.740 versus −0.142 in
test. A negative baseline means the training mean predicted slightly worse than each group's own mean. Predictions form
small clusters because people with the same four-SNP genotypes receive the same score.

**Takeaways**

- **Strong cis-eQTL evidence can translate into useful prediction.** Genotypes explain most of the variation in
  FLVCR1-AS1 expression among new people.
- **This is the strongest gene on the chromosome, by design,** so it does not describe a typical gene. Effects of a
  top-ranked gene also tend to be overestimated; the calibration slope below 1 reduces the predictions accordingly.
- **Don't over-read the gap between 0.605 and 0.740.** Each group has only about 35 people, so both estimates are
  imprecise. No model was changed after seeing the test result.
- **Prediction does not identify the causal variant.** The score SNPs may tag other variants, and the lead SNP has a
  perfect-LD partner in training (rs963328, with an identical p-value).
- **Genetic scores need a baseline in new people.** An earlier model for LIN9, a gene with weak evidence, could not beat
  the mean; its best validation result was the zero-SNP baseline. FLVCR1-AS1 was chosen by the training-only rule, not by
  prediction performance. The same test people were used to evaluate LIN9, but none of that enters this model. Predicted
  expression can only inform studies of traits and disease for genes whose expression is well predicted.

**Data separation:** the split reproduces the original 275/34/35 groups with seed 42. Betas, SNP eligibility, LD clumps,
and calibration use training people only, validation selected the threshold, and test expression was used only for the
final metrics. An earlier eligibility screen read all 344 people's expression to look for missing or constant values, but
excluded no gene, so it could not affect the selection.

## Week 5: 23andMe and Height Polygenic Scores

The original analysis used a private 23andMe file, and those personal results are not in this repository. The public
notebook runs the same workflow on **NA20792**, one person chosen at random (seed 42) from the 1000 Genomes EUR panel.
These results are that participant's and do not reproduce the private analysis.

| Quantity | Result |
|---|---:|
| Reference group | 488 people (the panel without NA20792) |
| GWAS rows / passed QC | 529,504 / 499,617 |
| Harmonized with the panel | 138,480 |
| Clump index SNPs | 67,385 |
| Scoring SNPs (p ≤ 0.001) | 8,699 |
| Reference mean / SD | 0.875 / 0.869 |
| NA20792 score | 2.652 |
| NA20792 percentile / z-score | 97.34 / 2.05 |

NA20792's score is high within this reference group, but the GWAS is simulated, so this says nothing about their height.

- **The weights are simulated.** The assigned height file comes from the PRS tutorial, which identifies its GWAS as
  simulated. The scores are an educational comparison, not a validated estimate of anyone's height predisposition.
- **The score cannot be read as height or as disease probability.** It has no observed-height calibration, so it cannot be
  converted to centimeters. A percentile is a rank within the 489-person EUR reference, not within the general population.
- **Everyone is scored on the same SNPs, weights, and effect alleles,** so the chip's incomplete coverage reduces all
  scores to the same comparable subset. The p ≤ 0.001 cutoff was fixed in advance and never tuned.
- **Clumping reduces redundant LD information.** It does not prove that a retained SNP is causal or fully independent.

## Week 6: LocusZoom and Linkage Disequilibrium

### Gradescope Answers

**Q1. How many significant cis-eQTLs on chromosome 22?**
**1,316 SNP–gene associations** have p < 0.05 / 280,657 = **1.7815 × 10⁻⁷**, from 280,657 chromosome 22 SNP–gene tests
across 633 genes. They involve **72 unique genes** and 957 unique SNP IDs.

The LD answers below use **FAM118A** (ENSG00000100376.7), cis-window chr22:45,204,849–46,237,836 on GRCh37 (gene body
45,704,849–45,737,836 ±500 kb), with **lead SNP rs104664** at chr22:45,711,854 (p = 3.26 × 10⁻⁸⁷, −log10(p) = 86.49).
Counts include the 635 tested non-lead SNPs. The lead SNP is excluded, and its self-LD is r² = 1.

| Question | Variants | Bonferroni-significant |
|---|---:|---:|
| **Q2.** LD r² > 0.8 with the lead variant | **10** | 10 |
| **Q3.** LD r² > 0.5 with the lead variant | **10** | 10 |
| **Q4.** LD r² > 0.1 with the lead variant | **48** | 19 |

**LocusZoom plot:** `data/week6_results/FAM118A_locuszoom.png`, also shown in the notebook.

### Takeaways

- **Q2 and Q3 match because of a real gap.** No SNP falls between 0.5 and 0.8; the values jump from 0.3078 to 0.8647.
- **LD and significance are separate.** All 10 high-LD SNPs are significant, while 29 of the 38 SNPs with r² between 0.1
  and 0.5 are not.
- **Nearby SNPs can show similar signals because of LD.** If one variant affects expression, SNPs correlated with it will
  also appear associated, which is why association strength tends to fall along with r² in the plot.
- **High LD does not show that a SNP is causal.** rs104664 has the smallest p-value, but 10 other SNPs have r² > 0.8 with
  it and very strong signals (p between about 10⁻⁸⁶ and 10⁻⁶⁵). This data cannot separate them, and the true causal
  variant may not be in this panel. The lead SNP is the best marker in this analysis, not a proven cause.
- **A second signal is neither shown nor ruled out.** The 8 significant SNPs with r² ≤ 0.1 sit near the gene with p around
  10⁻⁸. Very strong signals can carry over even through weak LD, so they may simply tag the main signal. We did not run a
  conditional analysis.
- **LD was fully computed.** All 636 reference SNPs in the window were tested, none has missing genotypes, and PLINK
  returned a valid r² for all 635 non-lead SNPs. The colors are unphased r² from the same 344 people as the regressions.
