"""Scientific and smoke tests that need only the public demo data.

Run from the repository root:  python -m unittest discover -s tests -v
"""
import hashlib
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from src.eqtl import add_multiple_testing_columns, fit_cis, summarize_genes
from src.locuszoom import merge_exons
from src.prs import gwas_qc, harmonize_to_reference, read_bim, relative_position, select_weights
from src.utils import REPO_ROOT

DEMO = REPO_ROOT / 'data/demo'
EXPECTED_OUTPUTS = ['demo_eqtl_associations.tsv', 'demo_eqtl_gene_summary.tsv', 'demo_eqtl_top_gene.png',
                    'demo_gwas_harmonized.tsv', 'demo_prs_distribution.png', 'demo_prs_summary.tsv',
                    'demo_reference_prs.tsv', 'demo_snp_counts.tsv', 'demo_weights.tsv']


def read_dosages(bfile, snps):
    """Decode a .bed file by hand: copies of the .bim column 5 allele for each person and SNP."""
    bim = read_bim(bfile)
    people = len(pd.read_csv(f'{bfile}.fam', sep=r'\s+', header=None))
    data = Path(f'{bfile}.bed').read_bytes()
    assert data[:3] == bytes([0x6c, 0x1b, 0x01]), 'Expected a SNP-major PLINK .bed file'
    width = (people + 3) // 4
    copies = {0b00: 2.0, 0b01: np.nan, 0b10: 1.0, 0b11: 0.0}
    dosages = {}
    for snp in snps:
        row = bim.index[bim['SNP'].eq(snp)][0]
        block = data[3 + row * width:3 + (row + 1) * width]
        dosages[snp] = [copies[(block[i // 4] >> (2 * (i % 4))) & 0b11] for i in range(people)]
    return dosages


class HarmonizationTests(unittest.TestCase):
    panel = pd.DataFrame({'CHR': [1, 1, 1, 1], 'SNP': ['same', 'swapped', 'strand', 'moved'], 'CM': 0.0,
                          'BP': [100, 200, 300, 400], 'R1': ['A', 'A', 'A', 'A'], 'R2': ['G', 'G', 'G', 'G']})

    def test_alleles_and_positions(self):
        gwas = pd.DataFrame({'CHR': [1, 1, 1, 1, 1], 'BP': [100, 200, 300, 999, 500],
                             'SNP': ['same', 'swapped', 'strand', 'moved', 'absent'],
                             'A1': ['A', 'G', 'T', 'A', 'A'], 'A2': ['G', 'A', 'C', 'G', 'G'],
                             'BETA': [0.1, 0.2, 0.3, 0.4, 0.5], 'P': 0.01})
        harmonized, counts = harmonize_to_reference(gwas, self.panel)
        result = harmonized.set_index('SNP')
        self.assertEqual(sorted(result.index), ['same', 'strand', 'swapped'])
        self.assertEqual((result.loc['swapped', 'A1'], result.loc['swapped', 'BETA']), ('G', 0.2))
        self.assertEqual((result.loc['strand', 'A1'], result.loc['strand', 'A2'], result.loc['strand', 'BETA']), ('A', 'G', 0.3))
        self.assertEqual((counts['Not in reference panel'], counts['Different position']), (1, 1))

    def test_incompatible_alleles_are_removed(self):
        gwas = pd.DataFrame({'CHR': [1], 'BP': [100], 'SNP': ['same'], 'A1': ['A'], 'A2': ['C'], 'BETA': [0.1], 'P': [0.01]})
        harmonized, counts = harmonize_to_reference(gwas, self.panel)
        self.assertTrue(harmonized.empty)
        self.assertEqual(counts['Incompatible alleles'], 1)

    def test_qc_removes_ambiguous_and_duplicate_snps(self):
        gwas = pd.DataFrame({'SNP': ['a', 'a', 'b', 'c', 'd'], 'A1': ['A', 'A', 'A', 'C', 'A'], 'A2': ['G', 'G', 'T', 'G', 'A'],
                             'P': [0.1, 0.2, 0.1, 0.1, 0.1]})
        passed = gwas_qc(gwas)
        self.assertEqual(passed['SNP'].tolist(), ['a'])
        self.assertEqual(passed['P'].tolist(), [0.1])


class ScoringFunctionTests(unittest.TestCase):
    def test_threshold_is_inclusive_and_needs_an_index_snp(self):
        table = pd.DataFrame({'CHR': 1, 'BP': [1, 2, 3, 4], 'SNP': list('abcd'), 'P': [0.0005, 0.001, 0.002, 0.0001]})
        self.assertEqual(select_weights(table, ['a', 'b', 'c'], 0.001)['SNP'].tolist(), ['a', 'b'])

    def test_percentile_and_z_score(self):
        position = relative_position(pd.Series([1.0, 2.0, 3.0, 4.0]), 3.0)
        self.assertEqual(position['Target_Percentile'], 62.5)
        self.assertAlmostEqual(position['Target_Z_Score'], 0.5 / pd.Series([1.0, 2.0, 3.0, 4.0]).std(ddof=1))


class EqtlFunctionTests(unittest.TestCase):
    def test_results_follow_people_not_file_order(self):
        rng = np.random.default_rng(1)
        ids = [f'p{i}' for i in range(40)]
        genotypes = rng.integers(0, 3, 40).astype(float)
        y = pd.Series(genotypes * 2 + rng.normal(size=40), index=ids)
        columns = {i: f'{i}_{i}' for i in ids}
        variants = pd.DataFrame({'CHR': [1], 'SNP': ['s'], 'POS': [1], 'COUNTED': ['A'], 'ALT': ['G'],
                                 **{columns[i]: [g] for i, g in zip(ids, genotypes)}})
        first = fit_cis(variants, ids, y, columns)
        second = fit_cis(variants, ids[::-1], y.sample(frac=1, random_state=2), columns)
        self.assertAlmostEqual(first['BETA'].iloc[0], second['BETA'].iloc[0])
        self.assertAlmostEqual(first['BETA'].iloc[0], np.polyfit(genotypes, y.to_numpy(), 1)[0])

    def test_corrections_use_the_whole_family(self):
        table = pd.DataFrame({'Gene_ID': ['g1', 'g1', 'g2', 'g2'], 'SNP': list('abcd'), 'P': [0.001, 0.02, 0.03, 0.9]})
        corrected = add_multiple_testing_columns(table)
        np.testing.assert_allclose(corrected['P_Bonferroni'], [0.004, 0.08, 0.12, 1.0])
        np.testing.assert_allclose(corrected['P_BH'], [0.004, 0.04, 0.04, 0.9])
        self.assertEqual(summarize_genes(corrected).loc['g1', 'BH_Hits'], 2)


class DemoRunTests(unittest.TestCase):
    """Runs the real demo twice, then checks what it produced."""

    @classmethod
    def setUpClass(cls):
        cls.folder = tempfile.TemporaryDirectory()
        cls.first, cls.second = Path(cls.folder.name) / 'first', Path(cls.folder.name) / 'second'
        for output in (cls.first, cls.second):
            subprocess.run([sys.executable, str(REPO_ROOT / 'scripts/run_demo.py'), '--output', str(output)],
                           check=True, capture_output=True, text=True, cwd=REPO_ROOT)
        cls.weights = pd.read_csv(cls.first / 'demo_weights.tsv', sep='\t')
        cls.gwas = pd.read_csv(DEMO / 'demo_gwas.tsv', sep='\t')
        cls.reference_scores = pd.read_csv(cls.first / 'demo_reference_prs.tsv', sep='\t', dtype={'IID': str}).set_index('IID')['PRS']
        cls.summary = pd.read_csv(cls.first / 'demo_prs_summary.tsv', sep='\t').set_index('Quantity')['Value']

    @classmethod
    def tearDownClass(cls):
        cls.folder.cleanup()

    def test_expected_outputs_exist(self):
        for name in EXPECTED_OUTPUTS:
            self.assertGreater((self.first / name).stat().st_size, 0, name)

    def test_results_are_reproducible(self):
        for name in EXPECTED_OUTPUTS:
            if name.endswith('.tsv'):
                digest = [hashlib.sha256((folder / name).read_bytes()).hexdigest() for folder in (self.first, self.second)]
                self.assertEqual(digest[0], digest[1], name)

    def test_demo_individual_is_excluded_from_reference(self):
        individual = pd.read_csv(DEMO / 'demo_individual.fam', sep=r'\s+', header=None, dtype=str)[1].tolist()
        self.assertEqual(len(individual), 1)
        self.assertNotIn(individual[0], self.reference_scores.index)

    def test_fixed_threshold_and_only_valid_snps_in_the_score(self):
        self.assertTrue(self.weights['P'].le(0.001).all())
        cases = self.gwas.drop_duplicates('SNP').set_index('SNP').loc[self.weights['SNP'], 'Demo_Case']
        self.assertTrue(cases.isin(['As simulated', 'Opposite strand']).all())

    def test_plink_scores_match_hand_counted_dosages(self):
        bim = read_bim(DEMO / 'demo_reference').set_index('SNP')
        for bfile, expected in [(DEMO / 'demo_reference', self.reference_scores.to_numpy()),
                                (DEMO / 'demo_individual', np.array([self.summary['Target_PRS']]))]:
            dosages = read_dosages(bfile, self.weights['SNP'])
            manual = np.zeros(len(expected))
            for row in self.weights.itertuples():
                column_5 = np.array(dosages[row.SNP])
                manual += row.BETA * (column_5 if row.A1 == bim.loc[row.SNP, 'R1'] else 2 - column_5)
            np.testing.assert_allclose(manual, expected, atol=1e-5)

    def test_score_points_the_right_way_on_held_out_people(self):
        self.assertGreater(self.summary['Reference_Correlation_With_Simulated_Phenotype'], 0)

    def test_eqtl_sample_alignment_and_grch37_positions(self):
        fam = pd.read_csv(DEMO / 'demo_eqtl.fam', sep=r'\s+', header=None, dtype=str)[1].tolist()
        self.assertEqual(fam, pd.read_csv(DEMO / 'demo_expression.tsv', sep='\t', nrows=0).columns[1:].tolist())
        summary = pd.read_csv(self.first / 'demo_eqtl_gene_summary.tsv', sep='\t')
        self.assertEqual(summary['Gene_Symbol'].iloc[0], 'FAM118A')
        associations = pd.read_csv(self.first / 'demo_eqtl_associations.tsv', sep='\t')
        lead = associations.loc[associations['Gene_ID'].eq(summary['Gene_ID'].iloc[0])].sort_values('P').iloc[0]
        # This lead SNP's GRCh37 position; in GRCh38 it would be different.
        self.assertEqual((lead['SNP'], lead['POS'], lead['N']), ('rs104664', 45711854, 344))


class LocusZoomTests(unittest.TestCase):
    def test_merge_exons_combines_overlapping_and_adjacent_exons_only(self):
        rows = pd.DataFrame({3: [100, 150, 201, 400], 4: [160, 200, 250, 450]})
        self.assertEqual(merge_exons(rows), [[100, 250], [400, 450]])


class PublicationSafetyTests(unittest.TestCase):
    @unittest.skipUnless((REPO_ROOT / '.git').exists(), 'Needs a Git repository')
    def test_no_private_files_would_be_published(self):
        sys.path.insert(0, str(REPO_ROOT / 'scripts'))
        from check_public_files import find_problems
        _, problems = find_problems()
        self.assertEqual(problems, [])


if __name__ == '__main__':
    unittest.main()
