"""Create notebooks/weekly_tasks_public.ipynb, a publishable copy of the private weekly notebook.

Weeks 1-4 and 6 use public data and keep their outputs. The private Week 5 section is left out entirely and
replaced by public_week5_cells.py, which runs the same workflow on a public participant. Those cells are only
re-executed when their code changes. The original notebook is not changed.
"""
import copy
import os
import re
import sys
from pathlib import Path

import nbformat
from nbclient import NotebookClient

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from public_week5_cells import CELLS
from src.utils import REPO_ROOT

SOURCE = REPO_ROOT / 'quater1_weekly_tasks.ipynb'
OUTPUT = REPO_ROOT / 'notebooks/weekly_tasks_public.ipynb'
CODE_REPLACEMENTS = [
    (r"'\./plink2_binary'", 'PLINK2'),
]
SETUP = """# Run from the repository root, with PLINK2 located by src/utils.py.
import os, sys
from pathlib import Path
if Path.cwd().name == 'notebooks':
    os.chdir('..')
sys.path.insert(0, '.')
from src.utils import plink2_path
PLINK2 = plink2_path()"""
INTRO = """# DSC 180A — Weekly Tasks (Public Copy)

The weekly course tasks that led to the Quarter 1 Project. Outputs are from the original run. Findings are in
[weekly_tasks_results.md](../weekly_tasks_results.md).

**The original Week 5 used private 23andMe genotype data, which is not in this repository.** Week 5 here runs the same
workflow on a public 1000 Genomes participant instead, so its results are that participant's. To reproduce results
without the full downloads, use the demo (see `README.md`)."""


def week_5_range(cells):
    start = next(i for i, cell in enumerate(cells) if cell.source.startswith('# Week 5'))
    end = next(i for i, cell in enumerate(cells) if cell.source.startswith('# Week 6'))
    return start, end


def build_week_5_cells():
    cells = []
    for number, (kind, source) in enumerate(CELLS):
        cell = nbformat.v4.new_markdown_cell(source) if kind == 'markdown' else nbformat.v4.new_code_cell(source)
        cell['id'] = f'public-week5-{number}'
        cells.append(cell)
    return cells


def add_outputs(cells):
    """Reuse the saved outputs if the code has not changed; otherwise run the cells from the repository root."""
    if OUTPUT.exists():
        saved = nbformat.read(OUTPUT, 4).cells
        start, end = week_5_range(saved)
        previous = saved[start:end]
        if [c.source for c in previous] == [c.source for c in cells] and any(c.get('outputs') for c in previous):
            return previous
    print('Running the public Week 5 cells...')
    os.environ.pop('MPLBACKEND', None)
    notebook = nbformat.v4.new_notebook(cells=copy.deepcopy(cells), metadata=nbformat.read(SOURCE, 4).metadata)
    NotebookClient(notebook, timeout=None, kernel_name='python3', record_timing=False,
                   resources={'metadata': {'path': str(REPO_ROOT)}}).execute()
    return notebook.cells


def main():
    notebook = nbformat.read(SOURCE, 4)
    start, end = week_5_range(notebook.cells)
    week_5 = add_outputs(build_week_5_cells())
    notebook.cells[start:end] = week_5
    for cell in notebook.cells:
        if cell.cell_type == 'code':
            for old, new in CODE_REPLACEMENTS:
                cell.source = re.sub(old, new, cell.source, flags=re.S)
        else:
            cell.source = cell.source.replace('(weekly_tasks_results.md)', '(../weekly_tasks_results.md)')
    intro, setup = nbformat.v4.new_markdown_cell(INTRO), nbformat.v4.new_code_cell(SETUP)
    intro['id'], setup['id'] = 'public-intro', 'public-setup'
    notebook.cells[:0] = [intro, setup]
    nbformat.validate(notebook)
    OUTPUT.parent.mkdir(exist_ok=True)
    nbformat.write(notebook, OUTPUT)
    print('Wrote', OUTPUT.relative_to(REPO_ROOT))


if __name__ == '__main__':
    main()
