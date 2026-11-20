"""Audit the saved Week 4 model read-only; do not rescore the historical test set."""
import json
import re
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import linregress
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.prediction import (THRESHOLDS, RESULTS, WORK, metrics, split_ids, load_expression,
                            training_associations, scores_from_frame, digest)
from src.prs import clump
from src.utils import REPO_ROOT

def main():
    source = REPO_ROOT/'data/week4_results'; output = RESULTS/'week4_audit'
    output.mkdir(parents=True, exist_ok=True)
    keep = ['FLVCR1-AS1_validation_metrics.tsv', 'FLVCR1-AS1_test_metrics.tsv',
            'FLVCR1-AS1_frozen_model.tsv', 'FLVCR1-AS1_final_weights.tsv',
            'FLVCR1-AS1_training_eqtl.tsv', 'FLVCR1-AS1_test_predictions.tsv']
    before = {name: digest(source/name) for name in keep}
    groups = split_ids(); expression = load_expression()
    frozen = pd.read_csv(source/'FLVCR1-AS1_frozen_model.tsv', sep='\t').iloc[0]
    gene = frozen.Gene_ID
    annotation = next(line for line in (source/'chr1.gencode12.GRCh37.gtf').open()
                      if f'gene_id "{gene}"' in line and '\tgene\t' in line)
    biotype = re.search(r'gene_type "([^"]+)"', annotation).group(1)
    assert biotype != 'protein_coding'
    variants = pd.read_csv(source/'FLVCR1-AS1_training_cis.traw', sep=r'\s+')
    y = expression.loc[gene, groups['train']].astype(float)
    recomputed = training_associations(variants, groups['train'], y).sort_values('SNP')
    saved = pd.read_csv(source/'FLVCR1-AS1_training_eqtl.tsv', sep='\t').sort_values('SNP')
    assert recomputed.SNP.tolist() == saved.SNP.tolist()
    np.testing.assert_allclose(recomputed.BETA, saved.BETA, rtol=1e-9)
    np.testing.assert_allclose(recomputed.P, saved.P, rtol=1e-8, atol=1e-15)
    wd = WORK/'week4_audit'; wd.mkdir(parents=True, exist_ok=True)
    indexes = clump(saved, REPO_ROOT/'data/LDREF/1000G.EUR.1', wd/'training',
                    keep=source/'train.keep')
    old = pd.read_csv(source/'FLVCR1-AS1_training.clumps', sep=r'\s+').ID
    assert set(indexes) == set(old)
    validation = pd.read_csv(source/'FLVCR1-AS1_validation_metrics.tsv', sep='\t')
    np.testing.assert_array_equal(validation.Threshold, THRESHOLDS)
    for i, row in validation.iterrows():
        weights = saved.loc[saved.SNP.isin(indexes) & saved.P.le(row.Threshold)]
        calibration = linregress(scores_from_frame(variants, weights, groups['train']), y)
        np.testing.assert_allclose([calibration.intercept, calibration.slope],
                                   [row.Calibration_intercept, row.Calibration_slope], rtol=1e-5)
    best = validation.sort_values(['Predictive_R2', 'Threshold'], ascending=[False, False]).iloc[0]
    assert best.Threshold == frozen.Threshold and best.SNPs == frozen.SNPs
    predictions = pd.read_csv(source/'FLVCR1-AS1_test_predictions.tsv', sep='\t')
    test = pd.read_csv(source/'FLVCR1-AS1_test_metrics.tsv', sep='\t', index_col=0).Value
    for quantity, value in metrics(predictions.Observed, predictions.Predicted).items():
        np.testing.assert_allclose(value, test[quantity], atol=1e-12)
    # Read saved test predictions only; no new test scoring or model selection is performed.
    result = {'Gene': frozen.Symbol, 'Gene_ID': gene, 'Protein_Coding_Cohort_Eligible': False,
              'Biotype': biotype + ' (GENCODE v12)', 'Thresholds': list(THRESHOLDS),
              'Split_N': {k: len(v) for k,v in groups.items()}, 'Training_Only_Effects_Verified': True,
              'Training_Only_LD_Clumps_Reproduced': True, 'Training_Calibration_Verified': True,
              'Selected_Threshold': float(frozen.Threshold), 'SNPs': int(frozen.SNPs),
              'Validation_Predictive_R2': float(best.Predictive_R2),
              'Test_Predictive_R2': float(test.Predictive_R2),
              'Test_Pearson_r': float(test.Pearson_r), 'Test_Pearson_r_squared': float(test.Pearson_r_squared),
              'Test_Baseline_R2': float(test.Baseline_predictive_R2),
              'Historical_Test_Audit': 'Read saved predictions; no new test scoring. Historical execution count cannot be proved from files alone.',
              'Original_Files_Unchanged': before == {name: digest(source/name) for name in keep}}
    (output/'verification.json').write_text(json.dumps(result, indent=2)+'\n')
    validation.to_csv(output/'validation_thresholds.tsv', sep='\t', index=False)
    print(json.dumps(result, indent=2))

if __name__ == '__main__': main()
