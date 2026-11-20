# December 4 submission: status and remaining work

Updated November 17, 2026. Part 1 is preserved; Part 2, exploratory analysis, Week 7,
Conda and Docker demonstration workflows are implemented and verified.

- Part 2 uses original RA, T2D and schizophrenia effects, never the local LDSC Z-scores
  as beta weights. Selection, coverage and limitations are recorded in the new methods/results notes.
- The same 489 EUR people contribute to correlations, SNP overlap and PC1–PC2 sensitivity
  analyses. No private ancestry results are generated or published.
- Week 7 uses checksum-verified ENCFF421XIL, lifted to GRCh37 before intersecting the
  saved Part 1 associations and GENCODE v12 genes. Negative lead-SNP overlap findings
  are retained alongside RPS26.
- The demo and automated tests run in a fresh Conda environment and Linux Docker.
- Paper support is in `docs/quarter1_part2_week7_methods_results.md`,
  `results/figure_inventory.tsv` and `docs/paper_figures.tex`. Overleaf was not modified.

Remaining work for the paper: integrate the reviewed captions and results into the course
LaTeX template, cite the original data/software, and check the instructor's rubric. Keep the
SNPs-per-gene histogram supplementary. Clinical risk estimation is outside this analysis;
score ranks are not disease probabilities. RA/T2D targeted-array coverage, ancestry,
uncalibrated effects and the original no-covariate eQTL model belong in the limitations.

The supplied expression values include negatives. The original preprocessing details still
need clarification from the course documentation or instructor; do not describe them as raw
RPKM. Original study archives and extracted GWAS effects remain local under their source
terms. Review the release checklist and generated file allowlist before any staging, commit
or push; none was performed during this work.
