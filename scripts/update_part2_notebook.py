"""Refresh the working Part 2 notebook and a separate reference-only public copy.

Only Part 2 report cells and the figure-export cell are executed. Part 1 and the
functional annotation section remain untouched. No Git changes are staged.
"""
import copy
import sys
from pathlib import Path

import nbformat
from nbclient import NotebookClient

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.utils import REPO_ROOT

SETUP = '''from pathlib import Path
import sys
import pandas as pd
from IPython.display import display, Image, Markdown
report_root = Path.cwd() if (Path.cwd() / 'src').is_dir() else Path.cwd().parent
if str(report_root) not in sys.path:
    sys.path.insert(0, str(report_root))
part2 = report_root / 'results/part2'
'''

LOCAL_CELLS = [
    ('markdown', '''### Tiffany's disease PRSs — the main comparison

This working notebook includes the actual local target comparison. The public copy is
`quarter1_project_public.ipynb`. The raw genotype files are not printed or embedded.

**A crucial comparison rule:** score Tiffany and every reference person using the same
compatible, nonmissing SNPs for each disease. Her chip covers only a subset of the
full reference scoring variants. Therefore, the histograms below rescore the reference
people on her shared SNP set. Her score must not be placed on a histogram using the
larger full-reference SNP set.

The table reports her raw score, reference percentile, reference-standardized score
(z-score), and the retained SNP coverage. A percentile is a relative rank, not disease
probability. No thresholds were selected using her disease outcomes.'''),
    ('code', SETUP + '''from scripts.private_part2_report import build_private_report
target_summary, coverage_sensitivity = build_private_report()
target_report = report_root / 'data/private_part2/report'
display(target_summary[['Disease', 'Shared_Scoring_SNPs', 'Full_Scoring_SNPs',
                        'SNPs_Retained_Percent', 'Reference_Individuals', 'Target_PRS',
                        'Reference_Mean', 'Reference_SD', 'Target_Z_Score',
                        'Target_Percentile']].round(4))'''),
    ('markdown', '''### Figure 1 — Where Tiffany lies in each disease score distribution

Each histogram contains the same EUR reference individuals. The **red line marks
Tiffany**; its legend gives her empirical percentile. The z-score reports how many
reference standard deviations her score is above or below the mean. Each panel states
its shared SNP count. These are genetic-score distributions, not distributions of
measured disease outcomes or calibrated disease probabilities.'''),
    ('code', "display(Image(filename=str(target_report / 'target_prs_distributions.png'), width=1100))"),
    ('markdown', '''### Figure 2 — Tiffany's standardized scores across diseases

This plot compares her **relative position** within each disease's matched reference
distribution. Zero is that disease's reference mean; positive values are above it and
negative values are below it. Percentile labels are calculated separately for each
disease. A larger z-score does not mean a larger absolute chance of that disease.'''),
    ('code', '''display(Image(filename=str(target_report / 'target_standardized_comparison.png'), width=900))
interpretation = []
for row in target_summary.itertuples(index=False):
    direction = 'above' if row.Target_Z_Score > 0 else 'below' if row.Target_Z_Score < 0 else 'at'
    interpretation.append(
        f"- **{row.Disease}:** Tiffany is at reference percentile "
        f"{row.Target_Percentile:.1f}, with z={row.Target_Z_Score:.2f} "
        f"({direction} the reference mean), using {row.Shared_Scoring_SNPs} shared SNPs.")
display(Markdown('\\n'.join(interpretation)))'''),
    ('markdown', '''### Figure 3 — How much does limited SNP coverage change the scores?

This sensitivity analysis uses **reference people only**. Each dot compares one
person's full-reference PRS with that same person's score using the subset available
for Tiffany, standardized separately. Agreement with the dashed diagonal means the
relative position is similar under the two sets. Correlations and percentile shifts
quantify sensitivity to coverage; they do not measure clinical prediction accuracy.

Tiffany has no known full-set score because the missing genotypes are unavailable.
We do not fill them in, compare her subset score against full-set scores, or infer
what her full-set percentile would have been. The retained-SNP fraction alone does
not describe how much predictive information remains.'''),
    ('code', '''display(coverage_sensitivity.round(4))
display(Image(filename=str(target_report / 'target_coverage_sensitivity.png'), width=1100))
coverage_notes = []
for row in coverage_sensitivity.itertuples(index=False):
    coverage_notes.append(
        f"- **{row.Disease}:** {row.SNPs_Retained_Percent:.1f}% of selected SNPs retained; "
        f"full-versus-shared Pearson r={row.Pearson_Full_vs_Shared:.2f}. "
        f"The average absolute reference percentile change is "
        f"{row.Mean_Absolute_Percentile_Shift:.1f} percentage points.")
display(Markdown('\\n'.join(coverage_notes)))'''),
]

EXPLORATORY_CELLS = [
    ('markdown', '''### Exploratory reference analysis — disease PRS correlations

**This was the second plot in the earlier notebook.** It compares the full-reference
disease scores across the same reference people; it is not a plot of Tiffany's scores.
Each square is a Pearson correlation between two disease scores. Positive values mean
people with higher scores for one disease tend to have higher scores for the other;
negative values mean the opposite. Values near zero indicate little linear association.
The diagonal is 1 because each score is compared with itself.

The table also reports Spearman rank correlations and PC1–PC2-adjusted correlations.
SNP overlap, LD, population structure and GWAS effect estimates can influence these
relationships. A correlation does not show that one disease causes another.'''),
    ('code', '''display(pd.read_csv(part2 / 'pairwise_correlations.tsv', sep='\\t').round(4))
display(Image(filename=str(part2 / 'figures/prs_correlations.png'), width=800))
display(pd.read_csv(part2 / 'snp_overlap.tsv', sep='\\t').round(4))'''),
    ('markdown', '''### Exploratory reference analysis — population PCA

**This was the third plot in the earlier notebook.** Each dot is a reference person,
not a disease. PC1 and PC2 summarize genetic differences in LD-pruned autosomal
genotypes; nearby points have similar positions along these two axes. The percentages
on the axes show the fraction of genotype variation explained by each component.

This panel describes structure **within the EUR reference panel**, not worldwide
ancestry. Clusters are unlabeled because validated population labels were unavailable.
Tiffany is not projected onto these PCs. The PRS–PC associations and PC-adjusted score
correlations are exploratory checks, not a correction that eliminates ancestry bias.'''),
    ('code', '''display(Image(filename=str(part2 / 'pca/reference_PCA.png'), width=800))
display(pd.read_csv(part2 / 'prs_pc_associations.tsv', sep='\\t').round(4))
display(Image(filename=str(part2 / 'figures/prs_pc_adjusted_correlations.png'), width=800))'''),
]

PUBLIC_CELLS = [
    ('markdown', '''### Public reference PRSs

This copy displays reference results only. Tiffany's values, personal plots and local
report cells are available in the working `quarter1_project.ipynb`, not this public copy.
The full-reference scores below use each disease's complete selected panel SNP set;
the working target comparisons instead use matched shared SNP subsets.'''),
    ('code', SETUP + '''display(pd.read_csv(part2 / 'selected_diseases.tsv', sep='\\t'))
display(pd.read_csv(part2 / 'prs_summary.tsv', sep='\\t'))
display(Image(filename=str(part2 / 'figures/prs_distributions.png'), width=1000))'''),
]


def cells_from(entries, prefix, private=False):
    cells = []
    for i, (kind, source) in enumerate(entries):
        # The entry strings above retain literal escape sequences for executable notebook code.
        if kind == 'code':
            source = source.replace("'\\\\t'", "'\\t'").replace("'\\\\n'", "'\\n'")
            compile(source, f'{prefix}-{i}', 'exec')
        make = nbformat.v4.new_markdown_cell if kind == 'markdown' else nbformat.v4.new_code_cell
        tags = ['part2-week7-report'] + (['private-results'] if private else [])
        cell = make(source, metadata={'tags': tags})
        cell.id = f'{prefix}-{i}'
        cells.append(cell)
    return cells


def refresh_notebooks():
    local_path = REPO_ROOT / 'notebooks/quarter1_project.ipynb'
    public_path = REPO_ROOT / 'notebooks/quarter1_project_public.ipynb'
    original = nbformat.read(local_path, as_version=4)
    start = next(i for i, c in enumerate(original.cells)
                 if c.cell_type == 'markdown' and c.source.startswith('## Part 2')
                 and 'Disease Polygenic Scores and Reference Structure' in c.source)
    end = next(i for i in range(start + 1, len(original.cells))
               if original.cells[i].cell_type == 'markdown'
               and 'Functional Annotation of cis-eQTLs' in original.cells[i].source)
    intro = copy.deepcopy(original.cells[start])
    preserved_extensions = [copy.deepcopy(c) for c in original.cells[start+1:end]
                            if 'expanded-research' in c.metadata.get('tags', [])]
    intro.source = intro.source.replace('Private comparisons remain in ignored local files.',
                                        'The working notebook displays the local target comparisons below.')
    for private in [False, True]:
        notebook = copy.deepcopy(original)
        entries = LOCAL_CELLS if private else PUBLIC_CELLS
        report = cells_from(entries, 'target-report' if private else 'reference-report', private)
        report += cells_from(EXPLORATORY_CELLS, 'reference-exploration')
        extension = nbformat.v4.new_notebook(cells=report, metadata=copy.deepcopy(original.metadata))
        NotebookClient(extension, timeout=180, kernel_name='python3', record_timing=False,
                       resources={'metadata': {'path': str(REPO_ROOT)}}).execute()
        notebook.cells[start:end] = [copy.deepcopy(intro), *extension.cells, *copy.deepcopy(preserved_extensions)]
        notebook.metadata['local_private_analysis'] = private
        for cell in notebook.cells:
            if cell.get('id') == 'paper-figure-collection-intro':
                folder = 'results/paper_figures/' if private else 'results/paper_figures_public/'
                cell.source = ('## Paper Figures — One Folder for LaTeX\n\n'
                               f'Use **`{folder}`** for this notebook. '
                               + ('This local collection includes Tiffany\'s marked plots and local '
                                  'methods/results notes. It is ignored by Git.' if private else
                                  'This public collection contains reference-only figures.')
                               + '\n\nThe following cell refreshes existing figure copies, captions, '
                                 'source inventory and LaTeX examples without rerunning analyses.')
            if cell.get('id') == 'paper-figure-collection-export':
                cell.source = (SETUP + 'from scripts.export_paper_figures import export_paper_figures\n'
                               f'paper_collection = export_paper_figures(include_private={private})\n'
                               'print(paper_collection)')
                cell.metadata['tags'] = ['part2-week7-report'] + (['private-results'] if private else [])
                execution = nbformat.v4.new_notebook(cells=[copy.deepcopy(cell)],
                                                     metadata=copy.deepcopy(original.metadata))
                NotebookClient(execution, timeout=120, kernel_name='python3', record_timing=False,
                               resources={'metadata': {'path': str(REPO_ROOT)}}).execute()
                cell.update(execution.cells[0])
        if not private:
            assert not any('private-results' in c.metadata.get('tags', []) for c in notebook.cells)
        assert notebook.cells[:start] == original.cells[:start]
        nbformat.validate(notebook)
        nbformat.write(notebook, local_path if private else public_path)
    print('Updated the working target analysis and a separate public reference-only notebook.')


if __name__ == '__main__':
    refresh_notebooks()
