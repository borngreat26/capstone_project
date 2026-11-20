"""Consistency checks on saved public results, skipped in the tiny Docker/demo-only image."""
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]

@unittest.skipUnless((ROOT/'results/part2/prs_summary.tsv').exists(),'Full saved public results not included')
class PublicResultsTests(unittest.TestCase):
    def test_score_tables_match_reference_summary(self):
        summary=pd.read_csv(ROOT/'results/part2/prs_summary.tsv',sep='\t').set_index('Disease')
        selected=pd.read_csv(ROOT/'results/part2/selected_diseases.tsv',sep='\t')
        z=pd.read_csv(ROOT/'results/part2/reference_standardized_prs.tsv',sep='\t',index_col='IID')
        for row in selected.itertuples():
            raw=pd.read_csv(ROOT/'results/part2'/row.Slug/'reference_scores.tsv',sep='\t',index_col='IID').PRS
            self.assertEqual(set(raw.index),set(z.index))
            self.assertEqual(len(raw),int(summary.loc[row.Disease,'Individuals']))
            self.assertAlmostEqual(raw.mean(),summary.loc[row.Disease,'Mean'])
            np.testing.assert_allclose(z[row.Disease],((raw-raw.mean())/raw.std(ddof=1)).loc[z.index],atol=1e-12)
    def test_pc_residuals_are_orthogonal_to_covariates(self):
        residual=pd.read_csv(ROOT/'results/part2/pc_adjusted_residuals.tsv',sep='\t',index_col='IID')
        pcs=pd.read_csv(ROOT/'results/part2/pca/reference_pcs.tsv',sep='\t',index_col='IID').loc[residual.index,['PC1','PC2']]
        design=np.column_stack([np.ones(len(pcs)),pcs.to_numpy()])
        np.testing.assert_allclose(design.T@residual.to_numpy(),0,atol=1e-10)
    def test_week7_unique_snps_and_associations_stay_distinct(self):
        table=pd.read_csv(ROOT/'results/week7/snp_overlap_summary.tsv',sep='\t')
        self.assertEqual(table.SNP_Gene_Pairs.iloc[0],6597508)
        self.assertTrue(table.Unique_SNPs.le(table.SNP_Gene_Pairs).all())
        self.assertTrue(table.Overlapping_SNPs.le(table.Overlapping_Pairs).all())
        np.testing.assert_allclose(table.Percent,100*table.Overlapping_SNPs/table.Unique_SNPs)
if __name__=='__main__':unittest.main()
