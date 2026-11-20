"""Verify frozen pilot outputs without selecting models or reevaluating test predictions."""
import json
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import linregress
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from src.prediction import (RESULTS,WORK,PRIMARY_RUN,digest,metrics,best_threshold,split_ids,load_expression,
                            chromosome_data,gene_data,scores_from_frame)
from src.utils import REPO_ROOT

def main():
    run=sys.argv[sys.argv.index('--run')+1] if '--run' in sys.argv else PRIMARY_RUN
    out=RESULTS/run;work=WORK/run;groups=split_ids()
    result=pd.read_csv(out/'gene_prediction_results.tsv',sep='\t')
    frozen=json.loads((out/'frozen_experiment.json').read_text())
    assert digest(out/'cohort.tsv')==frozen['Cohort_SHA256']
    assert digest(REPO_ROOT/'data/quarter1_results/part1/table_A_genome_wide_cis_eqtl.tsv.gz')==json.loads((out/'plan.json').read_text())['Discovery_Table_A_SHA256']
    annotation=pd.read_csv(REPO_ROOT/'data/quarter1_results/part1/protein_coding_genes.tsv',sep='\t')
    assert set(result.Gene_ID)<=set(annotation.Gene_ID)
    expression=load_expression();rows=[];validation_tables=[]
    for chrom,subset in result.groupby('Gene_Chr'):
        variants,g,columns=chromosome_data(chrom)
        for gene in subset.itertuples():
            folder=out/'models'/gene.Gene_ID
            assert digest(folder/'frozen_model.json')==frozen['Models'][gene.Gene_ID]
            model=json.loads((folder/'frozen_model.json').read_text())
            weights=pd.read_csv(folder/'weights.tsv',sep='\t')
            assert digest(folder/'weights.tsv')==model['Weights_SHA256']
            validation=pd.read_csv(folder/'validation.tsv',sep='\t')
            validation_tables.append(validation.assign(Gene_ID=gene.Gene_ID,Gene=gene.Gene))
            best=best_threshold(validation)
            if best is not None:
                assert best.Threshold==model['Threshold'] and best.SNPs==len(weights)
                assert weights.P.le(best.Threshold).all() and weights.N.eq(275).all()
                clumps=pd.read_csv(folder/'training_clump.clumps',sep=r'\s+')
                assert set(weights.SNP)<=set(clumps.ID)
                log=(folder/'training_clump.log').read_text()
                assert '--keep' in log and 'train.keep' in log and '275 samples remaining' in log
                frame=gene_data(gene,variants,g,columns)
                train_y=expression.loc[gene.Gene_ID,groups['train']].to_numpy(float)
                train_score=scores_from_frame(frame,weights,groups['train'])
                fit=linregress(train_score,train_y)
                np.testing.assert_allclose([fit.intercept,fit.slope],[model['Intercept'],model['Slope']],rtol=1e-10,atol=1e-10)
                # Independent single-SNP refits on 275 people, separate from vectorized fit_cis.
                for weight in weights.head(3).itertuples():
                    dosage=frame.set_index('SNP').loc[weight.SNP,[p+'_'+p for p in groups['train']]].to_numpy(float)
                    slope=linregress(dosage,train_y).slope
                    np.testing.assert_allclose(slope,weight.BETA,rtol=1e-9,atol=1e-10)
            predictions=pd.read_csv(work/'test_predictions'/f'{gene.Gene_ID}.tsv',sep='\t')
            assert predictions.IID.tolist()==groups['test']
            if len(weights):
                np.testing.assert_allclose(metrics(predictions.Observed,predictions.Predicted)['Predictive_R2'],gene.Test_Predictive_R2,atol=1e-12)
            np.testing.assert_allclose(metrics(predictions.Observed,predictions.Training_Mean_Baseline)['Predictive_R2'],gene.Baseline_Test_R2,atol=1e-12)
            rows.append({'Gene_ID':gene.Gene_ID,'Training_Only_Refits_And_Calibration':True,
                         'Training_Only_LD_Keep_Verified':True,'Validation_Selected_Threshold':True,
                         'Frozen_Model_Unchanged':True,'Saved_Test_Metrics_Reproduced':True})
    pd.DataFrame(rows).to_csv(out/'verification_checks.tsv',sep='\t',index=False)
    combined=pd.concat(validation_tables,ignore_index=True)
    combined.to_csv(out/'validation_all_thresholds.tsv',sep='\t',index=False)
    features=pd.read_csv(out/'candidate_training_features.tsv',sep='\t')
    features.loc[features.Exclusion.notna()].to_csv(out/'candidate_exclusions.tsv',sep='\t',index=False)
    result.loc[result.SNPs.eq(0),['Gene_ID','Gene','Status']].to_csv(out/'unmodeled_cohort_genes.tsv',sep='\t',index=False)
    print(f'Verified {len(rows)} frozen models, training-only effects/calibration/LD, validation selection, and saved test metrics. Discovery SHA256 unchanged.')

if __name__=='__main__':main()
