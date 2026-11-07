"""Create notebooks/weekly_tasks_public.ipynb, a publishable copy of the private weekly notebook.

Weeks 1-4 and 6 use public data and keep their outputs. Week 5 outputs are cleared, because Week 5 uses
private 23andMe data. Personal results are not in the notebook's text. The original notebook is not changed.
"""
import re
import sys
from pathlib import Path

import nbformat

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from src.utils import REPO_ROOT

SOURCE = REPO_ROOT / 'quater1_weekly_tasks.ipynb'
OUTPUT = REPO_ROOT / 'notebooks/weekly_tasks_public.ipynb'
# Patterns rather than literal values, so that this published script does not contain personal numbers either.
CODE_REPLACEMENTS = [
    (r"'\./plink2_binary'", 'PLINK2'),
    (r"shutil\.which\('bcftools'\) or '[^']*bcftools'", "shutil.which('bcftools')"),
    (r'display\(tiffany_raw\.head\(3\)\)\n', ''),
    (r"expected5 = \{.*?\}\n", ''),
    (r"    'Results match previously reported values':.*?\n", ''),
    (r"'All representations have [\d,]+ sites': (.*?) == \d+,", r"'All representations have the same sites': \1,"),
    (r"'Everyone scored on the same [\d,]+ unique SNPs': len\(final5\) == \d+ == ", "'Everyone scored on the same unique SNPs': len(final5) == "),
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

**Week 5 used private 23andMe genotype data, which is not in this repository.** Its code and methods are
shown, but its outputs and personal results are removed, and it cannot run without that file. This notebook is a
record of the course tasks. To reproduce results, use `notebooks/quarter1_project.ipynb` or the demo (see `README.md`)."""


def main():
    notebook = nbformat.read(SOURCE, 4)
    in_week_5 = False
    for cell in notebook.cells:
        if cell.cell_type == 'markdown' and cell.source.startswith('# Week'):
            in_week_5 = cell.source.startswith('# Week 5')
        if cell.cell_type == 'code':
            for old, new in CODE_REPLACEMENTS:
                cell.source = re.sub(old, new, cell.source, flags=re.S)
            if in_week_5:
                cell.outputs, cell.execution_count = [], None
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
