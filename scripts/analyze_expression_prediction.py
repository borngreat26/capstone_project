"""Summarize frozen held-out predictions; exploratory analyses do not retune models."""
import argparse
import json
import shutil
import sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr, false_discovery_control
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from src.prediction import RESULTS, WORK, PRIMARY_RUN, metrics, digest
from src.utils import REPO_ROOT
from scripts.run_week7 import save

R2_FLOOR=-1.  # display only: genes below this are drawn as edge markers, never dropped


def mark_clipped(ax, x, y, names):
    """Keep the R² axis readable when one gene is far below the rest; label what is clipped."""
    x,y,names=np.asarray(x,float),np.asarray(y,float),np.asarray(names)
    low=y<R2_FLOOR
    top=max(np.nanmax(y),0)+.08
    ax.set_ylim(R2_FLOOR-.12,top)
    if low.any():
        ax.scatter(x[low],np.full(low.sum(),R2_FLOOR-.06),marker='v',s=40,color='#bf6540',zorder=5,clip_on=False)
        for xi,yi,n in zip(x[low],y[low],names[low]):
            ax.annotate(f'{n} R²={yi:.1f}',(xi,R2_FLOOR-.06),xytext=(5,3),textcoords='offset points',fontsize=7,color='#bf6540')


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--run',default=PRIMARY_RUN);args=ap.parse_args()
    out=RESULTS/args.run;figs=out/'figures';figs.mkdir(exist_ok=True)
    result=pd.read_csv(out/'gene_prediction_results.tsv',sep='\t')
    tested=result.loc[result.SNPs.gt(0)&result.Test_Predictive_R2.notna()].copy()
    summary={'Cohort_Genes':len(result),'Nonempty_Models_Evaluated':len(tested),
             'No_Genetic_Model':int(result.SNPs.eq(0).sum()),
             'Mean_Test_Predictive_R2':tested.Test_Predictive_R2.mean(),
             'Median_Test_Predictive_R2':tested.Test_Predictive_R2.median(),
             'Minimum_Test_Predictive_R2':tested.Test_Predictive_R2.min(),
             'Percent_Positive_R2':100*tested.Test_Predictive_R2.gt(0).mean(),
             'Percent_Beating_Training_Mean':100*tested.Test_Beats_Baseline.mean(),
             'R2_Above_0_1':int(tested.Test_Predictive_R2.gt(.1).sum()),
             'R2_Above_0_3':int(tested.Test_Predictive_R2.gt(.3).sum()),
             'R2_Above_0_5':int(tested.Test_Predictive_R2.gt(.5).sum())}
    (out/'performance_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    result.groupby('Evidence').agg(Genes=('Gene_ID','size'),Mean_R2=('Test_Predictive_R2','mean'),
        Median_R2=('Test_Predictive_R2','median'),Percent_Positive=('Test_Predictive_R2',lambda x:100*x.gt(0).mean()),
        Percent_Beating_Baseline=('Test_Beats_Baseline',lambda x:100*x.mean())).to_csv(out/'performance_by_training_evidence.tsv',sep='\t')

    # ATAC SNP lookup was made from tested discovery variants, not from held-out performance.
    links=pd.read_csv(REPO_ROOT/'results/week7/snp_peak_links.tsv.gz',sep='\t',usecols=['SNP_ID'])
    accessible=set(links.SNP_ID)
    annotations=pd.read_csv(REPO_ROOT/'results/week7/gene_annotation.tsv',sep='\t').set_index('Gene_ID')
    for i,row in result.iterrows():
        association=pd.read_csv(WORK/args.run/'training_associations'/f'{row.Gene_ID}.tsv.gz',sep='\t')
        selected=pd.read_csv(out/'models'/row.Gene_ID/'weights.tsv',sep='\t')
        sig=association.loc[association.P.lt(.05/max(len(association),1))]
        result.loc[i,'Eligible_ATAC_SNPs']=int(association.SNP.isin(accessible).sum())
        result.loc[i,'Eligible_ATAC_Fraction']=association.SNP.isin(accessible).mean()
        result.loc[i,'Training_Significant_ATAC_SNPs']=int(sig.SNP.isin(accessible).sum())
        result.loc[i,'Selected_ATAC_SNPs']=int(selected.SNP.isin(accessible).sum())
        result.loc[i,'Selected_ATAC_Fraction']=selected.SNP.isin(accessible).mean() if len(selected) else np.nan
        result.loc[i,'Gene_Body_ATAC_Peaks']=annotations.loc[row.Gene_ID,'Body_Peaks']
        result.loc[i,'Cis_Window_ATAC_Peaks']=annotations.loc[row.Gene_ID,'Cis_Window_Peaks']
    result['Training_Signal_MinusLog10P']=-np.log10(result.Training_Min_P.clip(lower=1e-300))
    result.to_csv(out/'gene_prediction_features.tsv',sep='\t',index=False)
    features=['Training_Signal_MinusLog10P','Eligible_SNPs','Clumped_SNPs','SNPs','Gene_Length',
              'SNPs_Per_Mb','Training_Expression_SD','Training_Gene_Bonferroni_Hits',
              'Eligible_ATAC_Fraction','Training_Significant_ATAC_SNPs','Selected_ATAC_Fraction',
              'Gene_Body_ATAC_Peaks','Cis_Window_ATAC_Peaks']
    rows=[]
    for feature in features:
        subset=result[[feature,'Test_Predictive_R2']].dropna()
        if subset[feature].nunique()>1:
            r,p=spearmanr(subset[feature],subset.Test_Predictive_R2)
            rows.append({'Feature':feature,'Genes':len(subset),'Spearman_R':r,'P':p})
    associations=pd.DataFrame(rows);associations['BH_Q_Exploratory']=false_discovery_control(associations.P)
    associations.to_csv(out/'feature_performance_associations.tsv',sep='\t',index=False)
    # A descriptive ATAC sensitivity controlling measured locus/expression characteristics.
    controls=np.column_stack([np.ones(len(result)),np.log1p(result.Gene_Length),
                              np.log1p(result.Eligible_SNPs),np.log1p(result.Training_Expression_SD),
                              result.Training_Signal_MinusLog10P])
    good=np.isfinite(controls).all(axis=1)&result.Test_Predictive_R2.notna()&result.Eligible_ATAC_Fraction.notna()
    x=controls[good];y=result.loc[good,'Test_Predictive_R2'].to_numpy();a=result.loc[good,'Eligible_ATAC_Fraction'].to_numpy()
    yr=y-x@np.linalg.lstsq(x,y,rcond=None)[0];ar=a-x@np.linalg.lstsq(x,a,rcond=None)[0]
    partial=pearsonr(yr,ar)
    (out/'atac_adjusted_sensitivity.json').write_text(json.dumps({'Genes':int(good.sum()),
        'Residual_Pearson_R':float(partial.statistic),'Controls':'log gene length, log eligible SNP count, log training expression SD, training -log10(min p)',
        'Interpretation':'Exploratory residual association; not proof of regulation or causal mediation. No calibrated inferential p-value claimed.'},indent=2)+'\n')
    # List every autosome so chromosomes absent from the frozen pilot are visible, not silently dropped.
    by_chromosome=result.groupby('Gene_Chr').agg(Genes=('Gene_ID','size'),Mean_R2=('Test_Predictive_R2','mean'),
                                  Median_R2=('Test_Predictive_R2','median')).reindex(range(1,23))
    by_chromosome['Genes']=by_chromosome.Genes.fillna(0).astype(int);by_chromosome.index.name='Gene_Chr'
    by_chromosome.to_csv(out/'performance_by_chromosome.tsv',sep='\t')
    case=result.loc[result.Case_Study].iloc[0]
    case_folder=out/'case_study';case_folder.mkdir(exist_ok=True)
    for name in ['frozen_model.json','weights.tsv','validation.tsv']:
        shutil.copy2(out/'models'/case.Gene_ID/name,case_folder/name)
    prediction=pd.read_csv(WORK/args.run/'test_predictions'/f'{case.Gene_ID}.tsv',sep='\t')
    np.testing.assert_allclose(metrics(prediction.Observed,prediction.Predicted)['Predictive_R2'],case.Test_Predictive_R2)

    plt.rcParams.update({'font.size':10,'pdf.fonttype':42})
    fig,axes=plt.subplots(1,2,figsize=(10,4.4))
    limits=[min(prediction.Observed.min(),prediction.Predicted.min()),max(prediction.Observed.max(),prediction.Predicted.max())]
    axes[0].scatter(prediction.Observed,prediction.Predicted,s=26,color='#287d8e')
    axes[0].plot(limits,limits,'--',color='gray');axes[0].set(xlabel='Observed normalized expression',ylabel='Predicted normalized expression',
        title=f'{case.Gene}: n={int(case.Test_N)} held-out people\nPredictive R²={case.Test_Predictive_R2:.3f}; Pearson r²={case.Test_Pearson_r_squared:.3f}')
    validation=pd.read_csv(case_folder/'validation.tsv',sep='\t')
    axes[1].plot(range(5),validation.Predictive_R2,'o-',color='#287d8e')
    axes[1].axhline(case.Baseline_Validation_R2,ls='--',color='gray',label='Training-mean baseline')
    chosen=validation.index[validation.Threshold.eq(case.Threshold)][0]
    axes[1].scatter([chosen],[case.Validation_Predictive_R2],s=85,color='#bf6540',zorder=4,label='Selected by validation')
    axes[1].set(xticks=range(5),xticklabels=[f'{p:g}' for p in validation.Threshold],xlabel='Fixed p-value threshold',
                ylabel='Validation predictive R²',title=f'{case.Gene}: threshold selection\n{int(case.SNPs)} selected SNPs');axes[1].legend(frameon=False,fontsize=8)
    save(fig,figs/'single_gene_prediction')
    fig,axes=plt.subplots(1,2,figsize=(10,4))
    clipped=tested.loc[tested.Test_Predictive_R2.lt(R2_FLOOR)]
    axes[0].hist(tested.Test_Predictive_R2.clip(lower=R2_FLOOR),bins=25,color='#287d8e',edgecolor='white');axes[0].axvline(0,color='gray',ls='--')
    axes[0].set(xlabel=f'Held-out predictive R² (values below {R2_FLOOR:g} pooled in the leftmost bar)',ylabel='Genes',title=f'Nonempty genetic models (n={len(tested)})')
    if len(clipped):
        axes[0].annotate('Pooled: '+', '.join(f'{g} ({v:.1f})' for g,v in zip(clipped.Gene,clipped.Test_Predictive_R2)),
                         (.02,.95),xycoords='axes fraction',va='top',fontsize=8,color='#bf6540')
    for label,g in result.groupby('Evidence'):
        values=np.sort(g.Test_Predictive_R2.dropna());axes[1].plot(values,np.arange(1,len(values)+1)/len(values),label=f'{label} (n={len(values)})')
    axes[1].set(xlim=(R2_FLOOR-.05,tested.Test_Predictive_R2.max()+.05),xlabel=f'Held-out predictive R² (axis starts at {R2_FLOOR:g}; all genes counted)',ylabel='Cumulative fraction',title='Stratified by training-only evidence');axes[1].legend(frameon=False)
    save(fig,figs/'predictive_r2_distribution')
    colors={'Strong':'#bf6540','Moderate':'#287d8e','Weak':'#8e9ca5'}
    fig,ax=plt.subplots(figsize=(7,4.5))
    for label,g in result.groupby('Evidence'):
        ax.scatter(g.Training_Signal_MinusLog10P,g.Test_Predictive_R2,label=label,s=25,color=colors[label],alpha=.75)
    r=associations.set_index('Feature').loc['Training_Signal_MinusLog10P','Spearman_R']
    mark_clipped(ax,result.Training_Signal_MinusLog10P,result.Test_Predictive_R2,result.Gene)
    ax.axhline(0,color='gray',ls='--');ax.set(xlabel='−log10(minimum training cis-eQTL p-value)',ylabel='Held-out predictive R²',title=f'Training association strength versus prediction (Spearman r={r:.2f})');ax.legend(frameon=False)
    save(fig,figs/'association_strength_vs_prediction')
    fig,axes=plt.subplots(1,3,figsize=(12,4))
    for ax,feature,label in zip(axes,['Eligible_SNPs','Clumped_SNPs','SNPs'],['Eligible cis-SNPs','LD clump index SNPs','Selected scoring SNPs']):
        ax.scatter(result[feature],result.Test_Predictive_R2,s=20,color='#287d8e',alpha=.7)
        mark_clipped(ax,result[feature],result.Test_Predictive_R2,result.Gene)
        ax.axhline(0,color='gray',ls='--');ax.set(xlabel=label,ylabel='Held-out predictive R²',title=label)
    save(fig,figs/'snp_counts_vs_prediction')
    fig=plt.figure(figsize=(12,8));grid=fig.add_gridspec(2,3)
    spearman=associations.set_index('Feature').Spearman_R
    panels=[('Gene_Length','Gene length (bp)','log'),('SNPs_Per_Mb','Panel cis-SNPs per Mb','linear'),
            ('Training_Expression_SD','Training expression SD','log'),
            ('Training_Gene_Bonferroni_Hits','Within-gene Bonferroni-significant training SNPs','symlog')]
    for i,(feature,label,scale) in enumerate(panels):
        ax=fig.add_subplot(grid[i//3,i%3])
        ax.scatter(result[feature],result.Test_Predictive_R2,s=20,color='#287d8e',alpha=.7)
        ax.set(xlabel=label,ylabel='Held-out predictive R²',xscale=scale,title=f'Spearman r={spearman[feature]:.2f}')
        mark_clipped(ax,result[feature],result.Test_Predictive_R2,result.Gene)
        ax.axhline(0,color='gray',ls='--')
    ax=fig.add_subplot(grid[1,1:])
    rng=np.random.default_rng(0)
    for chrom in by_chromosome.index:
        values=result.loc[result.Gene_Chr.eq(chrom),'Test_Predictive_R2']
        ax.scatter(chrom+rng.uniform(-.2,.2,len(values)),values,s=14,color='#287d8e',alpha=.7)
        if len(values):ax.hlines(values.median(),chrom-.35,chrom+.35,color='#bf6540',lw=2)
    mark_clipped(ax,result.Gene_Chr,result.Test_Predictive_R2,result.Gene)
    ax.axhline(0,color='gray',ls='--')
    ax.set(xticks=by_chromosome.index,xticklabels=[f'{c}\n{n}' for c,n in by_chromosome.Genes.items()],
           xlabel='Chromosome (top) and pilot genes on it (bottom); bar = median',ylabel='Held-out predictive R²',
           title=f'By chromosome: {int(by_chromosome.Genes.gt(0).sum())} of 22 autosomes represented (descriptive)')
    ax.tick_params(axis='x',labelsize=7)
    fig.tight_layout()
    save(fig,figs/'gene_characteristics_vs_prediction')
    fig,axes=plt.subplots(1,2,figsize=(10,4))
    for label,g in result.groupby('Evidence'):
        axes[0].scatter(100*g.Eligible_ATAC_Fraction,g.Test_Predictive_R2,s=25,color=colors[label],label=label,alpha=.7)
    mark_clipped(axes[0],100*result.Eligible_ATAC_Fraction,result.Test_Predictive_R2,result.Gene)
    axes[0].set(xlabel='Training-eligible cis-SNPs overlapping ATAC peaks (%)',ylabel='Held-out predictive R²',title='B-cell accessibility overlap (GRCh37)');axes[0].legend(frameon=False)
    axes[1].scatter(result.Cis_Window_ATAC_Peaks,result.Test_Predictive_R2,s=23,color='#287d8e',alpha=.7)
    mark_clipped(axes[1],result.Cis_Window_ATAC_Peaks,result.Test_Predictive_R2,result.Gene)
    axes[1].set(xlabel='Peaks overlapping gene body ±500 kb',ylabel='Held-out predictive R²',title='Proximity annotation, not causal evidence')
    save(fig,figs/'atac_vs_prediction')

    timing=json.loads((out/'fit_runtime.json').read_text());test_timing=json.loads((out/'test_evaluation.json').read_text())
    models=out/'models';model_bytes=sum(p.stat().st_size for p in models.rglob('*') if p.is_file())
    association_bytes=sum(p.stat().st_size for p in (WORK/args.run/'training_associations').glob('*'))
    ratio=14661/len(result)
    # Deliberately conservative: scale the entire pilot fit time by model count, including
    # startup/cache costs. This is a projection, not a measured full-run benchmark.
    seconds=(timing['Fit_Seconds']+test_timing['Test_Evaluation_Seconds'])*ratio
    estimate={'Pilot_Fit_Seconds':timing['Fit_Seconds'],'Pilot_Test_Seconds':test_timing['Test_Evaluation_Seconds'],
              'Genes_In_Full_Target':14661,'Projected_Full_Run_Minutes':seconds/60,
              'Planning_Range_Minutes':[seconds/120,seconds/30],
              'Genotype_Cache_MB':timing['Genotype_Cache_Bytes']/1e6,
              'Peak_RAM_Measured':False,
              'Memory_Strategy':'One chromosome at a time; int8 dosage cache with memory mapping. Peak resident RAM was not instrumented.',
              'Projected_Models_And_Associations_MB':(model_bytes*ratio+association_bytes*14661/len(pd.read_csv(out/'candidate_genes.tsv',sep='\t')))/1e6,
              'Method':'Linear extrapolation; 0.5x–2x planning range. LD density and machine performance can change runtime/storage.',
              'Full_Genome_Run_Launched':False,
              'Fit_Command_After_Review':'python scripts/run_expression_prediction.py --run genomewide --all-genes --max-models 500 --stage fit',
              'Evaluate_Command_After_Freeze':'python scripts/run_expression_prediction.py --run genomewide --stage evaluate'}
    (out/'genomewide_resource_estimate.json').write_text(json.dumps(estimate,indent=2)+'\n')
    print(json.dumps(summary,indent=2));print(case[['Gene','Test_Predictive_R2','Test_Pearson_r_squared','Baseline_Test_R2']].to_string())

if __name__=='__main__':main()
