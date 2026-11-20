# Expanded Quarter 1 analysis: methods and actual results

Generated from saved artifacts by `scripts/write_prediction_notes.py`. This report separates
association discovery, expression prediction and disease PRS. No completed discovery outputs
were rebuilt. Full public-data inputs are required to reproduce the expanded analysis.

## Week 4 audit: preserve the noncoding demonstration

FLVCR1-AS1 (ENSG00000198468.2) is annotated as lincRNA (GENCODE v12); it is excluded from the
protein-coding pilot. The existing seed-42 split has 275 training, 34 validation and 35 test
people. Training-only SNP effects, LD clumps (r²=0.1, 250 kb), and score calibration reproduce
the saved files. The five inclusive cutoffs are 1, 0.05, 0.01, 0.001 and 0.00001.
Validation selected p≤1e-05 with 4 SNPs:
predictive R²=0.604984. The saved test predictive R² is
0.739890, Pearson r=0.874170, and Pearson r²=
0.764173. The training-mean test baseline is
-0.142318. Historical test predictions were read to audit metrics,
not newly scored. The files cannot prove how many historical notebook executions occurred.

## Protein-coding prediction protocol

The eligible population is the 14,661 protein-coding genes actually analyzed in discovery;
only their identities and annotation were reused for selection. Discovery p-values and
full-sample betas were not used to build or rank predictive models. Reuse of genotype
exports does not reuse expression effects. A public genotype cache stores one chromosome
at a time as int8 dosages. All score dosages are aligned to the training counted allele;
missing calls in selected scoring SNPs are errors rather than held-out imputation.

A seed-20261120 feature-stratified candidate pool contains 660 genes across
22 autosomes, spanning gene length, local SNP density and
training expression variability. RPS26 was included as a user-proposed candidate before
looking at validation/test prediction performance. Training-only cis-regressions then
classified candidate evidence using min p: strong <10⁻⁸; moderate [10⁻⁸,10⁻³); weak ≥10⁻³.
The frozen pilot has 150 genes across 22 chromosomes
(52 moderate, 44 strong, 54 weak), with the strongest training candidate retained as the single-gene
case study. Quotas were 50 per evidence class; the pool held only
44 strong candidates, so the pre-coded fallback filled the remainder from
the other tiers. Strata are shuffled within each chromosome and visited one chromosome at a
time, so every tier spreads across the genome. Sampling is reproducible; it is not a random
genome-wide benchmark, and the deliberately balanced evidence classes change its average.

This pilot replaces `pilot150`, which is kept frozen and unedited. Its sampler visited strata
in string-sorted chromosome order ('1', '10', '11', ...), so the partial pass that fills each
50-gene quota drew the moderate and weak tiers almost entirely from chromosomes 1, 10, 11 and
12, confounding chromosome with evidence class. The fix was mechanical and fixed before any
new test expression was read; the corrected pilot's cohort was inspected for coverage and
then evaluated once. `pilot150` results (median predictive R²=-0.056;
36.0% positive) are not used in this report.

All genes reuse the original 275/34/35 participant split. For each gene, additive SNP
effects and PLINK2 LD clumping use training people only. For each of the five thresholds,
the sum of beta × allele dosage is calibrated by a linear regression of training expression
on the training genetic score. Validation selects the nonempty, nonconstant genetic model
with maximum predictive R², preferring the larger p cutoff on a tie. Zero-SNP candidates
are labeled baseline models and excluded from genetic threshold selection. A nonempty
model can still perform worse than the training-mean baseline; that outcome is retained.

The cohort, weights, calibration and threshold choices were frozen and hashed before test
expression evaluation. Test performance is predictive R²=1−SSE/SST, with SST centered at
the test expression mean. Pearson r and r² are separate statistics. Baseline predictions
are the training expression mean. Test outcomes are never used to adjust these choices.
Per-gene predictions/metrics are cached; repeating evaluation reads saved results. The
historical Week 4 test split has previously been used for a different gene; for these new
genes it remains excluded from fitting/selection. The supplied expression was already
PEER-normalized before this split; upstream normalization was not relearned within each
split. Results describe this processed phenotype, not independent raw-expression replication.

## Single protein-coding example

RPS26 (ENSG00000197728.5) was the strongest training-only candidate (min p=
2.462e-114), with 299 eligible cis-SNPs and
42 clump index SNPs. Validation selected p≤0.01 and
10 scoring SNPs. Validation predictive R²=0.606747;
held-out predictive R²=0.528043; Pearson r=0.759589;
Pearson r²=0.576976. The training-mean test baseline is
-0.000027. Its frozen model, weights and five-row validation table are
copied to `results/expression_prediction/pilot150_v2/case_study/`.

## Multi-gene predictive performance

- Cohort: 150 genes; 150 nonempty genetic models; 0 without a genetic model.
- Mean test predictive R²: -0.127580; median: -0.030812; minimum: -13.019680.
- Positive predictive R²: 42.67%.
- Outperforming training-mean baseline: 50.00%.
- Models with predictive R²>0.1 / >0.3 / >0.5: 29 / 6 / 2.

The mean is sensitive to single genes: EPCAM has predictive
R²=-13.02 (validation R²=-3.16; test expression SD far below
training SD), so the median is the more representative summary. It was not removed.

Negative R² values were not removed. A negative predictive R² means larger SSE than a
constant prediction at the test-set mean; the deployable training-mean baseline is also
reported separately. Stronger training associations correlate with higher held-out R²
(Spearman r=0.465). These comparisons
are exploratory across related genes and the chosen pilot, not causal effects. Feature
association tables include BH q-values across the listed exploratory features; they do not
account for all gene–gene dependencies. With only 35 test people, per-gene R² is noisy.

## Gene features and accessible chromatin

Saved analyses examine eligible/clumped/selected SNP counts, gene length, local density,
training expression SD, within-gene Bonferroni-significant training pairs and chromosome.
The training significance count uses 0.05 divided by that gene's eligible SNP count,
not the original genome-wide discovery correction (Spearman r with predictive R²=
0.517). Gene/chromosome groups have small
sample counts and are descriptive. The pilot represents 22
of 22 autosomes with 4–11 genes each, too few for per-chromosome inference.
Strong-evidence genes are unevenly spread because few chromosomes contributed strong
candidates, so chromosome remains partly tied to evidence tier.

ATAC annotation reuses lifted GRCh37 B-cell peaks from Week 7. We record accessible
fractions among training-eligible and selected SNPs, significant accessible training
SNPs, and gene-body/cis-window peak counts. Accessible eligible-SNP fraction versus
predictive R² has Spearman r=0.048;
the count of significant accessible training SNPs has r=
0.403. Counts can reflect
association strength, gene length, LD and local density. Residual sensitivity controls
log gene length, log eligible SNP count, log training SD and training signal strength.
Neither overlap nor these associations establish a regulatory or causal mechanism;
primary B-cell accessibility is not identical to LCL accessibility.

## Genome-wide execution plan: not launched

The measured pilot fit was 15.83 seconds and final evaluation
1.34 seconds on this machine, with the genotype cache already built (the
superseded pilot, which built it, took 27.96 seconds). A linear projection
for 14,661 genes is 28.0 minutes, with a rough planning
range of 14.0–55.9 minutes.
The shared genotype cache is 452.2 MB; projected model and
association intermediates add approximately 919.0 MB.
This is an extrapolation, not a full-run benchmark. Peak resident RAM was not measured.

The workflow caches regressions and per-gene models, processes chromosomes sequentially,
and supports `--max-models` checkpoints. Full run commands, only after reviewing pilot
checks/resource estimates:

```bash
python scripts/run_expression_prediction.py --run genomewide --all-genes --max-models 500 --stage fit
# Repeat that fit command until all models are frozen; resume skips completed models.
python scripts/run_expression_prediction.py --run genomewide --stage evaluate
```

No full 14,661-gene prediction run was started. Original discovery tables and Week 4
outputs are unchanged. Exclusions and invalid/nonempty-model status are saved explicitly.

## Disease PRS extension: separate weights and outcomes

Disease scores use published GWAS log-odds effects, not the expression betas above.
Existing three-disease Pearson/Spearman correlations, SNP overlap and reference PCA are
reused. PC1–PC2 residual correlations remain exploratory, not a guarantee that ancestry
bias is eliminated. Fixed thresholds 1, 0.05, 0.01, 0.001 and 0.00001 were compared in
the same 489 public EUR individuals using the same pre-specified LD settings. Full
clumping followed by thresholding reproduces the original primary p≤0.001 SNP sets
and primary scores within PLINK text precision. No threshold was chosen using Tiffany's
outcomes or her percentile. These new sensitivities are reference-only; personal primary
results were not overwritten.

Threshold rank correlations and mean absolute percentile changes are saved by disease.
Ranks can change under different fixed weight sets; there is no independently validated
best threshold in this analysis.

## Disease SNP/eQTL overlap: descriptive, not colocalization

Only SNPs tested in discovery enter the overlap percentage denominator. The background
is the disease GWAS's successfully harmonized SNPs that were also tested for eQTLs,
rather than arbitrary genome positions. Selected PRS SNPs are LD-clumped while the
background is not; LD, targeted-array ascertainment and tested-gene coverage prevent a
simple enrichment or causal interpretation. SNP-level counts and SNP–gene pair counts
are kept separate. Variants outside tested regions are not evidence of no regulation.

- Rheumatoid arthritis: 62/221 tested scoring SNPs have ≥1 BH-significant eQTL (28.05%); separately 109 BH-significant SNP–gene pairs. Background overlap 15.53%.
- Type 2 diabetes: 47/240 tested scoring SNPs have ≥1 BH-significant eQTL (19.58%); separately 66 BH-significant SNP–gene pairs. Background overlap 15.02%.
- Schizophrenia: 323/3032 tested scoring SNPs have ≥1 BH-significant eQTL (10.65%); separately 450 BH-significant SNP–gene pairs. Background overlap 8.16%.

## Reproduce the pilot without altering discovery

```bash
python scripts/verify_week4_prediction.py
python scripts/run_expression_prediction.py --stage fit
python scripts/run_expression_prediction.py --stage evaluate
python scripts/verify_prediction_experiment.py
python scripts/analyze_expression_prediction.py
python scripts/expand_prs_exploration.py
python scripts/write_prediction_notes.py
python scripts/update_research_notebook.py
python scripts/make_paper_support.py
python scripts/export_paper_figures.py --include-private
```

Uses the repository Conda/Docker configuration. The prediction expansion added no packages;
bcftools, bedtools and ucsc-liftover are in `environment.yml` for the Week 7 ATAC annotation
that Analysis 5 reuses. `pilot150_v2` was fitted and evaluated in a fresh environment built from
`environment.yml` (Python 3.12.16). The superseded `pilot150` did not record its environment;
re-running its fit stage in that environment reproduced its pool, cohort, weights and
thresholds exactly. The full inputs are mounted for full-data runs; the small public demo remains separate. No private genotype file is modified or published.
