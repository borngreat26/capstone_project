"""Collect Quarter 1 figures for a LaTeX paper, with an explicit local-target mode.

Copies existing figures without regenerating them. The default is reference-only;
--include-private adds the target figures to an ignored local working collection.
Run from the repository root: python scripts/export_paper_figures.py
"""
import argparse
import csv
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.utils import REPO_ROOT

FIGURE_FOLDERS = {
    'results/part1/figures': 'part1',
    'results/part2/figures': 'part2',
    'results/part2/pca': 'part2',
    'results/week7/figures': 'week7',
    'results/expression_prediction/pilot150_v2/figures': 'prediction',
    'results/part2/exploration/figures': 'prs_exploration',
}


def latex_examples(rows, folder):
    lines = [r'% Requires \usepackage{graphicx} in the paper preamble.',
             '% These are separate figure examples; select the ones appropriate for the paper.']
    pdf_stems = {Path(row['Filename']).stem for row in rows
                 if Path(row['Filename']).suffix == '.pdf'}
    for row in rows:
        path = Path(row['Filename'])
        if path.suffix == '.png' and path.stem in pdf_stems:
            continue
        caption = row['Caption'].replace('≤', r'$\leq$').replace('²', r'$^2$')
        caption = caption.replace('±', r'$\pm$').replace('−', '-').replace('%', r'\%')
        label = path.stem.replace('_', '-')
        lines.extend([r'\begin{figure}[htbp]', r'\centering',
                      r'\includegraphics[width=\linewidth]{' + folder + '/' + path.name + '}',
                      r'\caption{' + caption + '}', r'\label{fig:' + label + '}',
                      r'\end{figure}', r'% In text: Figure~\ref{fig:' + label + '}.', ''])
    return '\n'.join(lines) + '\n'


def export_paper_figures(include_private=False):
    inventory = REPO_ROOT / 'results/figure_inventory.tsv'
    with inventory.open() as handle:
        original = list(csv.DictReader(handle, delimiter='\t'))
    if include_private:
        private_inventory = REPO_ROOT / 'data/private_part2/report/figure_inventory.tsv'
        if not private_inventory.exists():
            raise FileNotFoundError('Run scripts/private_part2_report.py before exporting target figures')
        with private_inventory.open() as handle:
            original.extend(csv.DictReader(handle, delimiter='\t'))
    if not original:
        raise ValueError('The figure inventory is empty; run scripts/make_paper_support.py.')

    # Validate every source before copying anything. Only inventoried public outputs qualify.
    planned = []
    for row in original:
        relative = Path(row['Filename'])
        parent = relative.parent.as_posix()
        folders = dict(FIGURE_FOLDERS)
        if include_private:
            folders['data/private_part2/report'] = 'personal'
        if parent not in folders or relative.suffix not in ('.png', '.pdf', '.svg'):
            raise ValueError('Figure source is outside the public figure folders')
        source = REPO_ROOT / relative
        if not source.is_file() or source.is_symlink():
            raise ValueError('Figure source is missing or is a symlink: ' + str(relative))
        name = folders[parent] + '_' + relative.name
        planned.append((source, name, row))
    if len({name for _, name, _ in planned}) != len(planned):
        raise ValueError('Two figures would have the same collection filename')

    output = REPO_ROOT / ('results/paper_figures' if include_private else 'results/paper_figures_public')
    output.mkdir(parents=True, exist_ok=True)
    rows = []
    for source, name, row in planned:
        shutil.copy2(source, output / name)
        rows.append({**row, 'Filename': name, 'Original_Filename': row['Filename'],
                     'Personal_Result': row['Filename'].startswith('data/private_part2/')})
    with (output / 'figure_inventory.tsv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)

    # The shareable copy assumes this directory is uploaded as "paper_figures" in Overleaf.
    (output / 'figures.tex').write_text(latex_examples(rows, output.name))
    # Keep the existing repository-level LaTeX support pointed at the same collection.
    if not include_private:
        (REPO_ROOT / 'docs/paper_figures.tex').write_text(
            latex_examples(rows, 'results/paper_figures_public'))

    notes = [
        (REPO_ROOT / 'docs/quarter1_part1_methods_results_notes.md', 'part1_methods_results.md'),
        (REPO_ROOT / 'docs/quarter1_part2_week7_methods_results.md', 'part2_week7_methods_results.md'),
        (REPO_ROOT / 'docs/expanded_prediction_methods_results.md', 'expanded_prediction_methods_results.md'),
    ]
    for source, name in notes:
        if source.exists():
            shutil.copy2(source, output / name)
    if include_private:
        shutil.copy2(REPO_ROOT / 'data/private_part2/report/methods_results.md',
                     output / 'personal_methods_results.md')

    count = len({Path(row['Filename']).stem for row in rows})
    scope = ('LOCAL WORKING COLLECTION: includes personal target results. Do not stage or publish '
             'this folder as a public repository artifact.' if include_private else
             'Public reference-only collection: personal target figures and results are excluded.')
    (output / 'README.md').write_text(f'''# Quarter 1 paper figures

{scope}

Give this entire folder to the agent writing the LaTeX paper. It contains {count}
distinct figures ({len(rows)} PNG/PDF/SVG files) from Part 1, Part 2 and Week 7.
The analysis notebooks can refresh it without rerunning any regressions or scores.

- `figure_inventory.tsv`: filenames, original locations, data sources, suggested captions,
  paper sections and interpretation limitations. `part1_`, `part2_` and `week7_` identify
  the analysis; the same stem in PNG and PDF is one figure in two formats.
- `figures.tex`: ready-to-adapt figure environments and cross-reference examples.
  Upload this directory as `{output.name}/` alongside the paper's main `.tex` file.
  Select the relevant examples rather than including every figure automatically.
- `part1_methods_results.md` and `part2_week7_methods_results.md`: public writing notes
  copied from `docs/`; relative links in these notes refer to the original repository.

Use PDFs where provided. New analysis PNGs are 320 dpi; preserved Part 1 PNGs retain
their original resolution. No image was resized or regenerated for this collection.
The SNPs-tested-per-gene histogram and exploratory sensitivity figures are supplementary
candidates. Captions and notes explain which results are associations or relative scores.

The disease PRSs use published disease GWAS effects, not our expression eQTL betas.
Raw PRSs, z-scores and percentiles are not calibrated disease probabilities. Primary
B-cell ATAC overlap does not prove a causal eQTL or a regulatory link to a gene.
Raw genotypes and original GWAS effects are excluded in both modes. Personal target figures
and a local results note are included only in the explicitly requested local collection.

Refresh after changing the source figures:

```bash
python scripts/make_paper_support.py
{'python scripts/private_part2_report.py' if include_private else ''}
{'python scripts/export_paper_figures.py --include-private' if include_private else ''}
```

To refresh copies using the existing inventory only:

```bash
python scripts/export_paper_figures.py{' --include-private' if include_private else ''}
```
''')
    result = {'Distinct_Figures': count, 'Figure_Files': len(rows),
              'Folder': str(output.relative_to(REPO_ROOT)) + '/'}
    if include_private:
        # The requested local handoff keeps earlier output locations intact.
        handoff = REPO_ROOT / 'paper/figures'
        handoff.mkdir(parents=True, exist_ok=True)
        for path in output.iterdir():
            if path.is_file():shutil.copy2(path, handoff / path.name)
        (handoff / 'figures.tex').write_text(latex_examples(rows, 'figures'))
        result['Paper_Handoff'] = 'paper/figures/'
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--include-private', action='store_true',
                        help='Include target plots in the ignored local working collection')
    print(export_paper_figures(include_private=parser.parse_args().include_private))
