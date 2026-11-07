"""Run the full Part 1 genome-wide cis-eQTL analysis and copy its summary results to results/part1/.

Usage, from the repository root:  python scripts/run_full_analysis.py
It needs the full public inputs described in README.md and takes a few minutes. Part 2 is not implemented yet.
"""
import argparse
import os
import shutil
import sys
from pathlib import Path

import nbformat
from nbclient import NotebookClient

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from src.utils import REPO_ROOT, plink2_path

NOTEBOOK = REPO_ROOT / 'notebooks/quarter1_project.ipynb'
WORKING_RESULTS = REPO_ROOT / 'data/quarter1_results/part1'
PUBLISHED_RESULTS = REPO_ROOT / 'results/part1'
REQUIRED_INPUTS = ['data/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt',
                   *[f'data/LDREF/1000G.EUR.{c}.{ext}' for c in range(1, 23) for ext in ('bed', 'bim', 'fam')]]
# The notebook saves its figures directly to results/part1/figures/.
PUBLISHED_FILES = ['table_B_gene_summary.tsv', 'table_C_chromosome_summary.tsv', 'correction_summary.tsv',
                   'excluded_genes.tsv', 'minor_allele_count_summary.tsv', 'top_gene_robustness.tsv',
                   'verification_checks.tsv']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--copy-only', action='store_true', help='Only copy existing results to results/part1/')
    args = parser.parse_args()
    if not args.copy_only:
        missing = [name for name in REQUIRED_INPUTS if not (REPO_ROOT / name).exists()]
        if missing:
            sys.exit(f'Missing {len(missing)} input file(s), starting with {missing[0]}. See "Data" in README.md.')
        print('Using PLINK2 at:', plink2_path())
        os.environ.pop('MPLBACKEND', None)  # Let the notebook save its figures inline.
        notebook = nbformat.read(NOTEBOOK, 4)
        NotebookClient(notebook, timeout=None, kernel_name='python3', record_timing=False, resources={'metadata': {'path': str(NOTEBOOK.parent)}}).execute()
        nbformat.write(notebook, NOTEBOOK)
    PUBLISHED_RESULTS.mkdir(parents=True, exist_ok=True)
    for name in PUBLISHED_FILES:
        shutil.copy2(WORKING_RESULTS / name, PUBLISHED_RESULTS / name)
    print(f'Part 1 finished. Summary results: {PUBLISHED_RESULTS.relative_to(REPO_ROOT)}; '
          f'full association table: {(WORKING_RESULTS / "table_A_genome_wide_cis_eqtl.tsv.gz").relative_to(REPO_ROOT)}')


if __name__ == '__main__':
    main()
