"""Boundary and nesting checks independent of the full dataset."""
import unittest
import numpy as np
import pandas as pd
from src.annotation import snp_counts,interval_counts

class BEDTests(unittest.TestCase):
    def setUp(self):
        self.peaks=pd.DataFrame({'Chrom':['1','1','1'],'Start':[100,110,200],'End':[200,120,201]})
    def test_one_based_snp_boundaries(self):
        np.testing.assert_array_equal(snp_counts(self.peaks,'chr1',[100,101,110,111,120,121,200,201,202]),[0,1,1,2,2,1,1,1,0])
    def test_half_open_gene_boundaries_and_nested_peaks(self):
        np.testing.assert_array_equal(interval_counts(self.peaks,'1',[0,100,120,200,201],[100,120,200,201,202]),[0,2,1,1,0])
    def test_absent_chromosome(self):
        np.testing.assert_array_equal(snp_counts(self.peaks,'2',[101,200]),[0,0])
if __name__=='__main__':unittest.main()
