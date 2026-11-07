"""Checks that src/ reproduces the verified earlier weeks. They need the full local data and skip without it.

These tests print no individual-level values. Run from the repository root:
    python -m unittest tests.test_local_full_data -v
"""
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from src.prs import clump, gwas_qc, harmonize_to_reference, read_bim, score_bfile
from src.utils import REPO_ROOT

DATA = REPO_ROOT / 'data'
WEEK4, WEEK5 = DATA / 'week4_results', DATA / 'week5_results'
AUTOSOMES = range(1, 23)


@unittest.skipUnless((DATA / 'Height.gwas.txt.gz').exists() and (WEEK5 / 'height_reference_harmonized.tsv.gz').exists()
                     and (DATA / 'LDREF/1000G.EUR.1.bed').exists(), 'Needs the full local Week 5 inputs and results')
class Week5ReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        gwas = pd.read_csv(DATA / 'Height.gwas.txt.gz', sep='\t')
        gwas['BETA'] = np.log(gwas['OR'])
        panel = pd.concat([read_bim(DATA / f'LDREF/1000G.EUR.{c}') for c in AUTOSOMES], ignore_index=True)
        cls.harmonized, cls.counts = harmonize_to_reference(gwas_qc(gwas), panel)

    def test_harmonization_matches_week_5(self):
        saved = pd.read_csv(WEEK5 / 'height_reference_harmonized.tsv.gz', sep='\t')
        self.assertEqual(self.harmonized['SNP'].tolist(), saved['SNP'].tolist())
        self.assertEqual(self.harmonized['A1'].tolist(), saved['A1'].tolist())
        np.testing.assert_allclose(self.harmonized['BETA'], saved['BETA'], rtol=1e-9)

    def test_clumping_and_threshold_match_week_5(self):
        indexes = []
        with tempfile.TemporaryDirectory() as folder:
            for chrom in AUTOSOMES:
                table = self.harmonized.loc[self.harmonized['CHR'].eq(chrom)]
                indexes += clump(table, DATA / f'LDREF/1000G.EUR.{chrom}', Path(folder) / f'chr{chrom}')
        clumped = self.harmonized.loc[self.harmonized['SNP'].isin(indexes)]
        selected = clumped.loc[clumped['P'].le(0.001)]
        saved = pd.read_csv(WEEK5 / 'height_clumped_p001.tsv', sep='\t')
        self.assertEqual(len(clumped), 67388)
        self.assertEqual(set(selected['SNP']), set(saved['SNP']))
        self.assertEqual(len(selected), 8700)

    def test_reference_scores_match_week_5(self):
        weights = pd.read_csv(WEEK5 / 'height_common_variants.tsv', sep='\t')
        saved = pd.read_csv(WEEK5 / 'reference_height_prs.tsv', sep='\t', dtype={'IID': str}).set_index('IID').iloc[:, 0]
        total = pd.Series(0.0, index=saved.index)
        with tempfile.TemporaryDirectory() as folder:
            for chrom in AUTOSOMES:
                chosen = weights.loc[weights['CHR'].eq(chrom)]
                total = total.add(score_bfile(DATA / f'LDREF/1000G.EUR.{chrom}', chosen, Path(folder) / f'chr{chrom}'))
        # PLINK prints each chromosome's score to six significant digits, so sums can differ by rounding.
        np.testing.assert_allclose(total.loc[saved.index], saved, atol=1e-6)


@unittest.skipUnless((WEEK4 / 'train.keep').exists(), 'Needs the local Week 4 sample split')
class Week4SplitTests(unittest.TestCase):
    def test_prediction_groups_do_not_overlap(self):
        groups = [set(pd.read_csv(WEEK4 / f'{name}.keep', sep='\t', header=None, dtype=str)[1]) for name in ('train', 'validation', 'test')]
        self.assertFalse(groups[0] & groups[1] or groups[0] & groups[2] or groups[1] & groups[2])
        self.assertEqual([len(g) for g in groups], [275, 34, 35])


if __name__ == '__main__':
    unittest.main()
