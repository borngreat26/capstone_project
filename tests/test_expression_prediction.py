"""Tests for prediction metrics, threshold selection and held-out isolation."""
import tempfile
import unittest
from collections import namedtuple
from pathlib import Path
from unittest.mock import patch
import numpy as np
import pandas as pd
from src.prediction import metrics, best_threshold, training_associations, fit_gene


class PredictionTests(unittest.TestCase):
    def test_squared_correlation_is_not_predictive_r2(self):
        result = metrics([1,2,3], [3,2,1])
        self.assertAlmostEqual(result['Pearson_r_squared'], 1.)
        self.assertAlmostEqual(result['Predictive_R2'], -3.)
        self.assertTrue(np.isnan(metrics([1,2,3], [2,2,2])['Pearson_r']))

    def test_threshold_uses_predictive_r2_and_excludes_baseline(self):
        table = pd.DataFrame({'Threshold':[1,.05,.01], 'SNPs':[4,2,0],
                              'Predictive_R2':[.2,.4,.6], 'Pearson_r_squared':[.9,.5,np.nan],
                              'Status':['Genetic score','Genetic score','Training-mean baseline (0 SNPs)']})
        self.assertEqual(best_threshold(table).Threshold, .05)
        self.assertIsNone(best_threshold(table.loc[table.SNPs.eq(0)]))

    def test_test_data_cannot_change_fitting_or_calibration(self):
        rng = np.random.default_rng(19)
        ids = [f'p{i}' for i in range(40)]
        groups = {'train':ids[:28], 'validation':ids[28:34], 'test':ids[34:]}
        g = rng.integers(0,3,size=(4,40))
        frame = pd.DataFrame({'CHR':1,'SNP':['s1','s2','s3','s4'], 'POS':[10,20,30,40],
                              'COUNTED':'A','ALT':'G', **{p+'_'+p:g[:,i] for i,p in enumerate(ids)}})
        expression = pd.Series(2*g[0] + .3*g[1] + rng.normal(size=40), index=ids)
        poisoned = frame.copy()
        for p in groups['test']: poisoned[p+'_'+p] = 2-poisoned[p+'_'+p]
        a = training_associations(frame, groups['train'], expression)
        poisoned_expression = expression.copy()
        poisoned_expression.loc[groups['test']] = 1e12
        b = training_associations(poisoned, groups['train'], poisoned_expression)
        pd.testing.assert_frame_equal(a,b)
        gene = namedtuple('Gene','Gene_ID Gene_Symbol Gene_Chr')('synthetic','Synthetic',1)
        with tempfile.TemporaryDirectory() as directory, patch('src.prediction.clump', return_value=a.SNP.tolist()):
            first = fit_gene(gene,frame,a,expression.loc[groups['train']].to_numpy(),
                             expression.loc[groups['validation']].to_numpy(),groups,Path(directory)/'a')
            second = fit_gene(gene,poisoned,b,expression.loc[groups['train']].to_numpy(),
                              expression.loc[groups['validation']].to_numpy(),groups,Path(directory)/'b')
            self.assertEqual(first,second)
            third = fit_gene(gene,frame,a,expression.loc[groups['train']].to_numpy(),
                             expression.loc[groups['validation']].to_numpy()*10,groups,Path(directory)/'c')
            v1 = pd.read_csv(Path(directory)/'a/validation.tsv',sep='\t')
            v3 = pd.read_csv(Path(directory)/'c/validation.tsv',sep='\t')
            np.testing.assert_allclose(v1[['Intercept','Slope']],v3[['Intercept','Slope']])


if __name__ == '__main__': unittest.main()
