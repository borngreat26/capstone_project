"""Add public prediction/PRS extensions without executing completed notebook sections."""
import copy
import sys
from pathlib import Path
import nbformat
from nbclient import NotebookClient
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from src.utils import REPO_ROOT

TAG='expanded-research'
SETUP="""from pathlib import Path
import json
import pandas as pd
from IPython.display import display, Image, Markdown
research_root = Path.cwd() if (Path.cwd() / 'src').is_dir() else Path.cwd().parent
prediction_root = research_root / 'results/expression_prediction'
pilot = prediction_root / 'pilot150_v2'
"""

EXPRESSION_CELLS=[
('markdown',r'''## Part 1: Single-Gene Expression Prediction

Association discovery and prediction answer different questions. Discovery tested nearby
SNP–gene associations in the full matched sample. Prediction now estimates **new training-only
effects**, selects a model with validation expression, and assesses it on held-out people.
The completed discovery effects are not used as weights for these test-evaluated predictors.

The shared seed-42 split is 275 training, 34 validation and 35 test people. Training people
alone determine SNP betas, LD clumps (r²=0.1, 250 kb) and expression-score calibration.
Validation selects among inclusive p-value thresholds **1, 0.05, 0.01, 0.001 and 0.00001**.
The selected model is frozen before test outcomes are evaluated.

**Predictive R² = 1 − SSE/SST** measures prediction error, using the observed evaluation-group
mean for SST. It can be negative. **Pearson r²** measures squared correlation and can remain
high despite poor calibration. The training-mean baseline is evaluated separately.

The Week 4 FLVCR1-AS1 demonstration is preserved. GENCODE v12 calls it lincRNA, so it is
excluded from the protein-coding pilot. Its saved results are audited read-only below;
the historical test set is not rescored for this audit.'''),
('code',SETUP+'''week4_audit = json.loads((prediction_root / 'week4_audit/verification.json').read_text())
display(pd.Series({key: week4_audit[key] for key in ['Gene', 'Biotype', 'Split_N',
    'Selected_Threshold', 'SNPs', 'Validation_Predictive_R2', 'Test_Predictive_R2',
    'Test_Pearson_r', 'Test_Pearson_r_squared', 'Test_Baseline_R2',
    'Training_Only_Effects_Verified', 'Training_Only_LD_Clumps_Reproduced',
    'Training_Calibration_Verified', 'Original_Files_Unchanged']}, name='Verified result').to_frame())'''),
('markdown', '''### Protein-coding example selected by training evidence

The example below is the strongest training-only association in the predefined candidate
pool, retained before validation/test performance is seen. RPS26 was a user-proposed
candidate; it still ranked strongest in the training-only scan. The observed-versus-predicted
panel uses only the 35 test people. The second panel shows validation threshold selection,
not test-based tuning. A zero-SNP baseline is never labeled a genetic predictor.'''),
('code', '''prediction_results = pd.read_csv(pilot / 'gene_prediction_results.tsv', sep='\t')
coding_example = prediction_results.loc[prediction_results['Case_Study']].iloc[0]
display(coding_example[['Gene', 'Gene_ID', 'Eligible_SNPs', 'Clumped_SNPs', 'SNPs',
    'Threshold', 'Validation_Predictive_R2', 'Test_Predictive_R2', 'Test_Pearson_r',
    'Test_Pearson_r_squared', 'Baseline_Test_R2']].to_frame('Result'))
display(pd.read_csv(pilot / 'case_study/validation.tsv', sep='\t')
        [['Threshold', 'SNPs', 'Predictive_R2', 'Pearson_r_squared', 'Status']])
display(Image(filename=str(pilot / 'figures/single_gene_prediction.png'), width=1000))'''),
('markdown', '''## Part 1: Genome-Wide Expression Prediction Exploration

This is a **150-gene pilot**, not a completed genome-wide prediction run. A reproducible
660-gene candidate pool spans all 22 autosomes, gene length, local SNP density and training
expression variability. Training min p defines strong (<10⁻⁸), moderate ([10⁻⁸,10⁻³)) and
weak (≥10⁻³) evidence. The pilot has 44 strong (every strong candidate in the pool),
52 moderate and 54 weak genes across all 22 autosomes. Genes were not selected using
validation or test prediction outcomes. This corrected pilot (`pilot150_v2`) replaces
`pilot150`, whose sampler visited chromosomes in string order and drew its moderate and weak
genes almost entirely from chromosomes 1, 10, 11 and 12; that run is kept frozen but not
reported. One gene (EPCAM, R² = −13.0) dominates the mean, so the median is the better
summary. This stratified
pilot is deliberately enriched by evidence class; its average is not a genome-wide estimate.

All genes share the same split. Full-sample eQTL betas are not reused. Genotype exports and
training-only computations are cached; per-gene models are frozen before final evaluation.
Negative test R² values and failures to beat the training-mean baseline are retained.'''),
('code', '''pilot_summary = json.loads((pilot / 'performance_summary.json').read_text())
display(pd.Series(pilot_summary, name='Held-out result').to_frame())
display(pd.read_csv(pilot / 'performance_by_training_evidence.tsv', sep='\t').round(4))
display(Image(filename=str(pilot / 'figures/predictive_r2_distribution.png'), width=1000))'''),
('markdown', '''### What is associated with prediction accuracy?

Training association strength, SNP counts, gene length, density and training expression
variability are explored below. Every association-strength value is from **training only**.
Selected SNP count is influenced by the validation-selected threshold; its relationship
with test accuracy is not a causal effect of adding SNPs. Feature associations and BH
q-values are exploratory across related genes, not independent confirmatory findings.'''),
('code', '''feature_associations = pd.read_csv(pilot / 'feature_performance_associations.tsv', sep='\t')
display(feature_associations.round(5))
display(Image(filename=str(pilot / 'figures/association_strength_vs_prediction.png'), width=850))
display(Image(filename=str(pilot / 'figures/snp_counts_vs_prediction.png'), width=1100))
display(Image(filename=str(pilot / 'figures/gene_characteristics_vs_prediction.png'), width=1100))'''),
('markdown', '''### Prediction and B-cell accessible chromatin

Week 7's GRCh37 accessibility annotations are reused for training-eligible and selected
SNPs. The count of significant training variants uses a **within-gene Bonferroni** family,
distinct from the original genome-wide discovery correction. Counts can be confounded
by gene length, SNP density, LD and training signal strength. An additional residual
sensitivity controls measured locus/expression characteristics. Neither overlap nor a
performance association demonstrates a regulatory mechanism; primary B cells differ from LCLs.'''),
('code', '''display(Image(filename=str(pilot / 'figures/atac_vs_prediction.png'), width=1000))
display(pd.Series(json.loads((pilot / 'atac_adjusted_sensitivity.json').read_text()),
                  name='Exploratory adjusted comparison').to_frame())'''),
('markdown', '''### Verification and the optional full-genome run

The model/threshold manifest was written before test evaluation. Rerunning evaluation
loads saved results rather than revisiting model selection. The supplied expression was
already normalized before this split; upstream PEER normalization was not refitted within
training. With 35 test people, per-gene R² estimates are noisy and are not clinical evidence.

The resource estimate below is a conservative extrapolation, not a full-run benchmark.
**No full 14,661-gene prediction run has been launched.** The CLI processes chromosomes
and supports `--max-models` batches and per-gene resume checkpoints. Run commands and
limitations are in `docs/expanded_prediction_methods_results.md`.'''),
('code', '''prediction_checks = pd.read_csv(pilot / 'verification_checks.tsv', sep='\t')
print('Models passing every saved verification check:',
      int(prediction_checks.drop(columns='Gene_ID').all(axis=1).sum()), '/', len(prediction_checks))
display(pd.Series(json.loads((pilot / 'genomewide_resource_estimate.json').read_text()),
                  name='Planning estimate').to_frame())'''),
]

PRS_CELLS=[
('markdown', '''## Part 2: Disease PRS Exploratory Analysis — Thresholds and eQTL Overlap

The earlier three-disease correlations, shared-SNP results and LD-pruned PCA are retained.
Disease weights come from published GWASs; they are separate from the newly estimated
expression-prediction betas. This extension uses **public reference people only** and does
not alter Tiffany's primary results.

### Fixed-threshold sensitivity

P+T scores were calculated at the five predefined cutoffs. The heatmaps compare reference
rankings between thresholds using Spearman correlations. No threshold was selected to
make a target person's rank high or low. The primary p≤0.001 SNP sets/scores reproduce
the previous analysis. Rank stability does not identify a clinically optimal threshold.'''),
('code',SETUP+'''prs_exploration = research_root / 'results/part2/exploration'
display(pd.read_csv(prs_exploration / 'threshold_summary.tsv', sep='\t')
        [['Disease', 'Threshold', 'Scoring_SNPs', 'Reference_N', 'Primary']])
display(Image(filename=str(prs_exploration / 'figures/threshold_rank_stability.png'), width=1100))'''),
('markdown', '''### Do disease-scoring SNPs also appear in significant cis-eQTL associations?

The table distinguishes unique SNPs from SNP–gene pairs. Percentages use scoring SNPs
that were actually tested for eQTLs. The comparison background is the disease GWAS's
harmonized variants that were tested for eQTLs. Untested SNPs are not evidence of no effect.

Selected score SNPs were LD-clumped while the background was not; LD, gene coverage,
array ascertainment and local density can influence overlap. These are descriptive
variant intersections, **not formal colocalization or evidence of causal pleiotropy**.'''),
('code', '''score_eqtl_overlap = pd.read_csv(prs_exploration / 'disease_eqtl_overlap.tsv', sep='\t')
display(score_eqtl_overlap[['Disease', 'Score_SNPs', 'Score_SNPs_Tested_For_eQTL',
    'Score_BH_eQTL_SNPs', 'Score_Bonferroni_eQTL_SNPs', 'BH_Pairs', 'Bonferroni_Pairs',
    'Score_BH_Percent_Of_Tested', 'Background_BH_Percent']].round(3))
display(Image(filename=str(prs_exploration / 'figures/disease_eqtl_overlap.png'), width=950))'''),
]


def make_cells(entries,prefix):
    result=[]
    for i,(kind,source) in enumerate(entries):
        if kind=='code':compile(source,f'{prefix}-{i}','exec')
        factory=nbformat.v4.new_markdown_cell if kind=='markdown' else nbformat.v4.new_code_cell
        cell=factory(source,metadata={'tags':[TAG]});cell.id=f'{prefix}-{i}';result.append(cell)
    return result


def update_notebooks():
    example=nbformat.read(REPO_ROOT/'notebooks/quarter1_project_public.ipynb',as_version=4)
    expression=make_cells(EXPRESSION_CELLS,'expression-extension')
    prs=make_cells(PRS_CELLS,'prs-extension')
    report=nbformat.v4.new_notebook(cells=expression+prs,metadata=copy.deepcopy(example.metadata))
    NotebookClient(report,timeout=120,kernel_name='python3',record_timing=False,
                   resources={'metadata':{'path':str(REPO_ROOT)}}).execute()
    expression=report.cells[:len(expression)];prs=report.cells[len(expression):]
    for name in ['quarter1_project.ipynb','quarter1_project_public.ipynb']:
        path=REPO_ROOT/'notebooks'/name;notebook=nbformat.read(path,as_version=4)
        before=copy.deepcopy(notebook.cells[:25])
        notebook.cells=[c for c in notebook.cells if TAG not in c.metadata.get('tags',[])]
        part2=next(i for i,c in enumerate(notebook.cells) if c.cell_type=='markdown' and c.source.startswith('## Part 2'))
        notebook.cells[part2:part2]=copy.deepcopy(expression)
        annotation=next(i for i,c in enumerate(notebook.cells) if c.cell_type=='markdown' and 'Functional Annotation of cis-eQTLs' in c.source)
        notebook.cells[annotation:annotation]=copy.deepcopy(prs)
        for cell in notebook.cells:
            if cell.get('id')=='paper-figure-collection-intro':
                if 'Expanded prediction and PRS figures are included.' not in cell.source:
                    cell.source += '\n\nExpanded prediction and PRS figures are included. The complete local handoff is also copied to `paper/figures/`.'
            if cell.get('id')=='paper-figure-collection-export':
                export=nbformat.v4.new_notebook(cells=[copy.deepcopy(cell)],metadata=copy.deepcopy(notebook.metadata))
                NotebookClient(export,timeout=120,kernel_name='python3',record_timing=False,
                               resources={'metadata':{'path':str(REPO_ROOT)}}).execute()
                cell.update(export.cells[0])
        assert notebook.cells[:25]==before
        nbformat.validate(notebook);nbformat.write(notebook,path)
    print('Added executed expression-prediction and expanded PRS sections; completed sections were not executed.')

if __name__=='__main__':update_notebooks()
