# Plan for the December 4 Final Submission

Nothing in this file has been done yet unless it says so.

## Part 2: Disease Polygenic Scores

**First decision:** where the effect sizes come from. The local GWAS files have only Z-scores
([details](gwas_sumstats_inventory.md)). Confirm the choice with the instructor before building on it.

For each eligible disease, reusing `src/prs.py`:

1. Harmonize GWAS variants with the 1000 Genomes panel.
2. Clump with r² 0.1 and 250 kb, as in Week 5.
3. Apply a fixed p-value threshold, justified in advance. Without independent phenotype data, it cannot be optimized.
4. Save the SNP weights and the number of SNPs kept at each step.
5. Score all eligible reference people. Score the private individual **locally only**.
6. Calculate the percentile and z-score, plot the distribution, and save a disease-level summary.

Private scores, plots, and summaries stay in `data/`, which Git ignores. Public versions use the reference group only.

**Interpretation:** a percentile is a rank within the reference group. It is not a disease probability or a
relative risk, and no absolute risk will be calculated without a calibrated model.

## Relationships Between Disease Scores

Using the same reference people:

- **Correlations:** pairwise Pearson correlations, shown as a heatmap.
- **Scatterplots:** selected disease pairs, each with its correlation.
- **Standardized comparison:** each score standardized with the reference mean and SD, then compared across diseases.
- **Interpretation:** consider shared genetic architecture, LD, ancestry, and overlap between GWAS cohorts. A
  correlation alone is not evidence of a causal relationship or a shared mechanism.
- **Optional sensitivity analysis:** other justified thresholds or SNP coverage levels, labeled as sensitivity
  analyses. They will not be tuned on any individual's phenotype.

## Ancestry Principal Components

The local genotype panel is EUR only, so PCs here describe **structure within this reference panel**. They cannot
assign anyone to a worldwide ancestry group.

- **Method:** LD-prune the panel, then run `plink2 --pca`.
- **Uses:** plot the structure, and check how each disease score correlates with the top PCs.
- **eQTLs:** a PC-adjusted scan would be a separate, labeled analysis. The original no-covariate results stay as they are.
- **Scores:** scores will not be residualized automatically, and PCs do not fix every ancestry-related limitation.
- **Privacy:** no ancestry inference or PC analysis will be run or published for the private individual.

## Paper and Repository

- Write the paper, using the numbers in `quarter1_part1_results.md` and the Part 2 outputs.
- Add citations for the datasets, GWAS studies, methods, and readings.
- Add Part 2 to `scripts/run_full_analysis.py` and `README.md`.
- Check the paper against the instructor's rubric.

## Questions for the Instructor

1. Should Part 2 use original GWAS effect sizes, or are Z-based approximate weights acceptable?
2. May the demo's small subset of real expression values be published in a public repository?
3. How exactly was the supplied expression file processed? Its values include negatives, so they are not raw RPKM.
