"""Generate methods/results notes for the frozen pilot and expanded PRS exploration."""
import json
import sys
from pathlib import Path
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from src.prediction import RESULTS, PRIMARY_RUN
from src.utils import REPO_ROOT

def main():
    out=RESULTS/PRIMARY_RUN
    summary=json.loads((out/'performance_summary.json').read_text())
    audit=json.loads((RESULTS/'week4_audit/verification.json').read_text())
    resource=json.loads((out/'genomewide_resource_estimate.json').read_text())
    results=pd.read_csv(out/'gene_prediction_results.tsv',sep='\t')
    case=results.loc[results.Case_Study].iloc[0]
    features=pd.read_csv(out/'feature_performance_associations.tsv',sep='\t').set_index('Feature')
    candidates=pd.read_csv(out/'candidate_genes.tsv',sep='\t')
    candidate_features=pd.read_csv(out/'candidate_training_features.tsv',sep='\t')
    superseded=json.loads((RESULTS/'pilot150/performance_summary.json').read_text())
    overlap=pd.read_csv(REPO_ROOT/'results/part2/exploration/disease_eqtl_overlap.tsv',sep='\t')
    sensitivity=pd.read_csv(REPO_ROOT/'results/part2/exploration/threshold_rank_stability.tsv',sep='\t')
    text=f'''# Expanded Quarter 1 analysis: methods and actual results

Generated from saved artifacts by `scripts/write_prediction_notes.py`. This report separates
association discovery, expression prediction and disease PRS. No completed discovery outputs
were rebuilt. Full public-data inputs are required to reproduce the expanded analysis.

## Week 4 audit: preserve the noncoding demonstration

FLVCR1-AS1 ({audit['Gene_ID']}) is annotated as {audit['Biotype']}; it is excluded from the
protein-coding pilot. The existing seed-42 split has 275 training, 34 validation and 35 test
people. Training-only SNP effects, LD clumps (r²=0.1, 250 kb), and score calibration reproduce
the saved files. The five inclusive cutoffs are 1, 0.05, 0.01, 0.001 and 0.00001.
Validation selected p≤{audit['Selected_Threshold']:g} with {audit['SNPs']} SNPs:
predictive R²={audit['Validation_Predictive_R2']:.6f}. The saved test predictive R² is
{audit['Test_Predictive_R2']:.6f}, Pearson r={audit['Test_Pearson_r']:.6f}, and Pearson r²=
{audit['Test_Pearson_r_squared']:.6f}. The training-mean test baseline is
{audit['Test_Baseline_R2']:.6f}. Historical test predictions were read to audit metrics,
not newly scored. The files cannot prove how many historical notebook executions occurred.

## Protein-coding prediction protocol

The eligible population is the 14,661 protein-coding genes actually analyzed in discovery;
only their identities and annotation were reused for selection. Discovery p-values and
full-sample betas were not used to build or rank predictive models. Reuse of genotype
exports does not reuse expression effects. A public genotype cache stores one chromosome
at a time as int8 dosages. All score dosages are aligned to the training counted allele;
missing calls in selected scoring SNPs are errors rather than held-out imputation.

A seed-20261120 feature-stratified candidate pool contains {len(candidates)} genes across
{candidates.Gene_Chr.nunique()} autosomes, spanning gene length, local SNP density and
training expression variability. RPS26 was included as a user-proposed candidate before
looking at validation/test prediction performance. Training-only cis-regressions then
classified candidate evidence using min p: strong <10⁻⁸; moderate [10⁻⁸,10⁻³); weak ≥10⁻³.
The frozen pilot has {len(results)} genes across {results.Gene_Chr.nunique()} chromosomes
({', '.join(f'{n} {k.lower()}' for k,n in results.Evidence.value_counts().sort_index().items())}), with the strongest training candidate retained as the single-gene
case study. Quotas were 50 per evidence class; the pool held only
{int(candidate_features.Evidence.eq('Strong').sum())} strong candidates, so the pre-coded fallback filled the remainder from
the other tiers. Strata are shuffled within each chromosome and visited one chromosome at a
time, so every tier spreads across the genome. Sampling is reproducible; it is not a random
genome-wide benchmark, and the deliberately balanced evidence classes change its average.

This pilot replaces `pilot150`, which is kept frozen and unedited. Its sampler visited strata
in string-sorted chromosome order ('1', '10', '11', ...), so the partial pass that fills each
50-gene quota drew the moderate and weak tiers almost entirely from chromosomes 1, 10, 11 and
12, confounding chromosome with evidence class. The fix was mechanical and fixed before any
new test expression was read; the corrected pilot's cohort was inspected for coverage and
then evaluated once. `pilot150` results (median predictive R²={superseded['Median_Test_Predictive_R2']:.3f};
{superseded['Percent_Positive_R2']:.1f}% positive) are not used in this report.

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

{case.Gene} ({case.Gene_ID}) was the strongest training-only candidate (min p=
{case.Training_Min_P:.3e}), with {int(case.Eligible_SNPs)} eligible cis-SNPs and
{int(case.Clumped_SNPs)} clump index SNPs. Validation selected p≤{case.Threshold:g} and
{int(case.SNPs)} scoring SNPs. Validation predictive R²={case.Validation_Predictive_R2:.6f};
held-out predictive R²={case.Test_Predictive_R2:.6f}; Pearson r={case.Test_Pearson_r:.6f};
Pearson r²={case.Test_Pearson_r_squared:.6f}. The training-mean test baseline is
{case.Baseline_Test_R2:.6f}. Its frozen model, weights and five-row validation table are
copied to `results/expression_prediction/pilot150_v2/case_study/`.

## Multi-gene predictive performance

- Cohort: {summary['Cohort_Genes']} genes; {summary['Nonempty_Models_Evaluated']} nonempty genetic models; {summary['No_Genetic_Model']} without a genetic model.
- Mean test predictive R²: {summary['Mean_Test_Predictive_R2']:.6f}; median: {summary['Median_Test_Predictive_R2']:.6f}; minimum: {summary['Minimum_Test_Predictive_R2']:.6f}.
- Positive predictive R²: {summary['Percent_Positive_R2']:.2f}%.
- Outperforming training-mean baseline: {summary['Percent_Beating_Training_Mean']:.2f}%.
- Models with predictive R²>0.1 / >0.3 / >0.5: {summary['R2_Above_0_1']} / {summary['R2_Above_0_3']} / {summary['R2_Above_0_5']}.

The mean is sensitive to single genes: {(worst:=results.nsmallest(1,'Test_Predictive_R2').iloc[0]).Gene} has predictive
R²={worst.Test_Predictive_R2:.2f} (validation R²={worst.Validation_Predictive_R2:.2f}; test expression SD far below
training SD), so the median is the more representative summary. It was not removed.

Negative R² values were not removed. A negative predictive R² means larger SSE than a
constant prediction at the test-set mean; the deployable training-mean baseline is also
reported separately. Stronger training associations correlate with higher held-out R²
(Spearman r={features.loc['Training_Signal_MinusLog10P','Spearman_R']:.3f}). These comparisons
are exploratory across related genes and the chosen pilot, not causal effects. Feature
association tables include BH q-values across the listed exploratory features; they do not
account for all gene–gene dependencies. With only 35 test people, per-gene R² is noisy.

## Gene features and accessible chromatin

Saved analyses examine eligible/clumped/selected SNP counts, gene length, local density,
training expression SD, within-gene Bonferroni-significant training pairs and chromosome.
The training significance count uses 0.05 divided by that gene's eligible SNP count,
not the original genome-wide discovery correction (Spearman r with predictive R²=
{features.loc['Training_Gene_Bonferroni_Hits','Spearman_R']:.3f}). Gene/chromosome groups have small
sample counts and are descriptive. The pilot represents {(chromosomes:=pd.read_csv(out/'performance_by_chromosome.tsv',sep='\t')).Genes.gt(0).sum()}
of 22 autosomes with {chromosomes.Genes.min()}–{chromosomes.Genes.max()} genes each, too few for per-chromosome inference.
Strong-evidence genes are unevenly spread because few chromosomes contributed strong
candidates, so chromosome remains partly tied to evidence tier.

ATAC annotation reuses lifted GRCh37 B-cell peaks from Week 7. We record accessible
fractions among training-eligible and selected SNPs, significant accessible training
SNPs, and gene-body/cis-window peak counts. Accessible eligible-SNP fraction versus
predictive R² has Spearman r={features.loc['Eligible_ATAC_Fraction','Spearman_R']:.3f};
the count of significant accessible training SNPs has r=
{features.loc['Training_Significant_ATAC_SNPs','Spearman_R']:.3f}. Counts can reflect
association strength, gene length, LD and local density. Residual sensitivity controls
log gene length, log eligible SNP count, log training SD and training signal strength.
Neither overlap nor these associations establish a regulatory or causal mechanism;
primary B-cell accessibility is not identical to LCL accessibility.

## Genome-wide execution plan: not launched

The measured pilot fit was {resource['Pilot_Fit_Seconds']:.2f} seconds and final evaluation
{resource['Pilot_Test_Seconds']:.2f} seconds on this machine, with the genotype cache already built (the
superseded pilot, which built it, took {json.loads((RESULTS/'pilot150/fit_runtime.json').read_text())['Fit_Seconds']:.2f} seconds). A linear projection
for 14,661 genes is {resource['Projected_Full_Run_Minutes']:.1f} minutes, with a rough planning
range of {resource['Planning_Range_Minutes'][0]:.1f}–{resource['Planning_Range_Minutes'][1]:.1f} minutes.
The shared genotype cache is {resource['Genotype_Cache_MB']:.1f} MB; projected model and
association intermediates add approximately {resource['Projected_Models_And_Associations_MB']:.1f} MB.
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

'''
    for row in overlap.itertuples():
        text+=f'- {row.Disease}: {row.Score_BH_eQTL_SNPs}/{row.Score_SNPs_Tested_For_eQTL} tested scoring SNPs have ≥1 BH-significant eQTL ({row.Score_BH_Percent_Of_Tested:.2f}%); separately {row.BH_Pairs} BH-significant SNP–gene pairs. Background overlap {row.Background_BH_Percent:.2f}%.\n'
    text+='\n## Reproduce the pilot without altering discovery\n\n```bash\npython scripts/verify_week4_prediction.py\npython scripts/run_expression_prediction.py --stage fit\npython scripts/run_expression_prediction.py --stage evaluate\npython scripts/verify_prediction_experiment.py\npython scripts/analyze_expression_prediction.py\npython scripts/expand_prs_exploration.py\npython scripts/write_prediction_notes.py\npython scripts/update_research_notebook.py\npython scripts/make_paper_support.py\npython scripts/export_paper_figures.py --include-private\n```\n\nUses the repository Conda/Docker configuration. The prediction expansion added no packages;\nbcftools, bedtools and ucsc-liftover are in `environment.yml` for the Week 7 ATAC annotation\nthat Analysis 5 reuses. `pilot150_v2` was fitted and evaluated in a fresh environment built from\n`environment.yml` (Python 3.12.16). The superseded `pilot150` did not record its environment;\nre-running its fit stage in that environment reproduced its pool, cohort, weights and\nthresholds exactly. The full inputs are mounted for full-data runs; the small public demo remains separate. No private genotype file is modified or published.\n'
    (REPO_ROOT/'docs/expanded_prediction_methods_results.md').write_text(text)
    print('Wrote expanded methods/results notes from actual frozen outputs.')

if __name__=='__main__':main()
