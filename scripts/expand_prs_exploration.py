"""Fixed-threshold sensitivity and exploratory disease-SNP/eQTL overlap.

Uses saved harmonized GWASs and public reference data only. No target outcomes or
target percentiles are used to choose thresholds; primary p<=0.001 is unchanged.
"""
import itertools
import sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from src.prediction import THRESHOLDS
from src.prs import clump
from src.utils import REPO_ROOT, run_plink
from scripts.run_week7 import save


def threshold_sensitivity():
    out=REPO_ROOT/'results/part2/exploration';work=REPO_ROOT/'data/part2_sensitivity'
    out.mkdir(parents=True,exist_ok=True);work.mkdir(exist_ok=True)
    if (out/'threshold_scores.tsv').exists():
        return
    selected=pd.read_csv(REPO_ROOT/'results/part2/selected_diseases.tsv',sep='\t')
    results=[];summaries=[];correlations=[]
    for study in selected.itertuples():
        wd=work/study.Slug;wd.mkdir(exist_ok=True)
        harmonized=pd.read_csv(REPO_ROOT/'data/part2_work'/study.Slug/'harmonized.tsv.gz',sep='\t')
        total=None;counts=np.zeros(len(THRESHOLDS),int);primary_ids=set()
        for chrom,table in harmonized.groupby('CHR'):
            prefix=wd/f'chr{int(chrom)}'
            bfile=REPO_ROOT/f'data/LDREF/1000G.EUR.{int(chrom)}'
            if Path(str(prefix)+'.clumps').exists():
                ids=pd.read_csv(str(prefix)+'.clumps',sep=r'\s+').ID
            else:
                ids=clump(table,bfile,prefix,r2=.1,kb=250)
            independent=table.loc[table.SNP.isin(ids)].copy()
            weights=independent[['SNP','A1']].copy()
            for j,cutoff in enumerate(THRESHOLDS):
                weights[f'W{j}']=independent.BETA.where(independent.P.le(cutoff),0.)
                counts[j]+=int(independent.P.le(cutoff).sum())
            primary_ids.update(independent.loc[independent.P.le(.001),'SNP'])
            weights.to_csv(str(prefix)+'_threshold_weights.tsv',sep='\t',index=False)
            run_plink(['--bfile',bfile,'--score',str(prefix)+'_threshold_weights.tsv',
                       '1','2','header-read','no-mean-imputation','cols=+scoresums',
                       '--score-col-nums','3-7','--out',str(prefix)+'_scores'])
            scored=pd.read_csv(str(prefix)+'_scores.sscore',sep=r'\s+',dtype={'IID':str}).set_index('IID')
            assert scored.ALLELE_CT.eq(2*len(weights)).all()
            frame=scored[[f'W{j}_SUM' for j in range(len(THRESHOLDS))]]
            if total is None:total=frame
            else:
                assert set(total.index)==set(frame.index)
                total=total+frame.reindex(total.index)
        primary=pd.read_csv(REPO_ROOT/'results/part2'/study.Slug/'weights.tsv',sep='\t')
        assert primary_ids==set(primary.SNP), 'Primary clump/threshold SNP set changed'
        saved=pd.read_csv(REPO_ROOT/'results/part2'/study.Slug/'reference_scores.tsv',sep='\t',index_col='IID').PRS
        np.testing.assert_allclose(total['W3_SUM'].loc[saved.index],saved,atol=1e-4,rtol=1e-5)
        for j,cutoff in enumerate(THRESHOLDS):
            score=total[f'W{j}_SUM'];sd=score.std(ddof=1)
            z=(score-score.mean())/sd if sd>0 else score*np.nan
            percentile=100*(score.rank(method='average')-.5)/len(score) if counts[j] else score*np.nan
            results.append(pd.DataFrame({'Disease':study.Disease,'Threshold':cutoff,'IID':score.index,
                                         'PRS':score.to_numpy(),'Z':z.to_numpy(),
                                         'Reference_Percentile':percentile.to_numpy()}))
            summaries.append({'Disease':study.Disease,'Threshold':cutoff,'Scoring_SNPs':int(counts[j]),
                              'Reference_N':len(score),'Mean':score.mean(),'SD':sd,'Primary':cutoff==.001})
        for a,b in itertools.combinations(range(len(THRESHOLDS)),2):
            x,y=total.iloc[:,a],total.iloc[:,b]
            shift=100*(x.rank(method='average')-y.rank(method='average')).abs()/len(x)
            correlations.append({'Disease':study.Disease,'Threshold_A':THRESHOLDS[a],'Threshold_B':THRESHOLDS[b],
                                 'Reference_N':len(x),'SNPs_A':int(counts[a]),'SNPs_B':int(counts[b]),
                                 'Pearson_R':x.corr(y),'Spearman_R':x.corr(y,method='spearman'),
                                 'Mean_Absolute_Percentile_Shift':shift.mean()})
        print('Threshold sensitivity complete:',study.Disease,flush=True)
    pd.concat(results,ignore_index=True).to_csv(out/'threshold_scores.tsv',sep='\t',index=False)
    pd.DataFrame(summaries).to_csv(out/'threshold_summary.tsv',sep='\t',index=False)
    pd.DataFrame(correlations).to_csv(out/'threshold_rank_stability.tsv',sep='\t',index=False)


def eqtl_overlap():
    out=REPO_ROOT/'results/part2/exploration'
    if (out/'disease_eqtl_overlap.tsv').exists():return
    selected=pd.read_csv(REPO_ROOT/'results/part2/selected_diseases.tsv',sep='\t')
    sets={};pair_counts={};backgrounds={}
    for study in selected.itertuples():
        sets[study.Disease]=set(pd.read_csv(REPO_ROOT/'results/part2'/study.Slug/'weights.tsv',sep='\t').SNP)
        backgrounds[study.Disease]=set(pd.read_csv(REPO_ROOT/'data/part2_work'/study.Slug/'harmonized.tsv.gz',sep='\t',usecols=['SNP']).SNP)
        pair_counts[study.Disease]=dict(All_Pairs=0,BH_Pairs=0,Bonferroni_Pairs=0)
    tested=set();bh=set();bonf=set();first=True
    import gzip
    with gzip.open(out/'disease_score_eqtl_links.tsv.gz','wt') as handle:
        for chunk in pd.read_csv(REPO_ROOT/'data/quarter1_results/part1/table_A_genome_wide_cis_eqtl.tsv.gz',
                                 sep='\t',usecols=['SNP_ID','Gene_ID','Gene_Symbol','P_BH','P_Bonferroni'],chunksize=250000):
            tested.update(chunk.SNP_ID);bh.update(chunk.loc[chunk.P_BH.lt(.05),'SNP_ID'])
            bonf.update(chunk.loc[chunk.P_Bonferroni.lt(.05),'SNP_ID'])
            for disease,snps in sets.items():
                chosen=chunk.loc[chunk.SNP_ID.isin(snps)].copy()
                pair_counts[disease]['All_Pairs']+=len(chosen)
                pair_counts[disease]['BH_Pairs']+=int(chosen.P_BH.lt(.05).sum())
                pair_counts[disease]['Bonferroni_Pairs']+=int(chosen.P_Bonferroni.lt(.05).sum())
                if len(chosen):
                    chosen.insert(0,'Disease',disease)
                    chosen.to_csv(handle,sep='\t',index=False,header=first);first=False
    rows=[]
    for disease,snps in sets.items():
        eligible=snps&tested; background=backgrounds[disease]&tested
        rows.append({'Disease':disease,'Score_SNPs':len(snps),'Score_SNPs_Tested_For_eQTL':len(eligible),
                     'Score_BH_eQTL_SNPs':len(eligible&bh),'Score_Bonferroni_eQTL_SNPs':len(eligible&bonf),
                     'Score_BH_Percent_Of_Tested':100*len(eligible&bh)/len(eligible),
                     'Score_Bonferroni_Percent_Of_Tested':100*len(eligible&bonf)/len(eligible),
                     'Harmonized_Tested_Background_SNPs':len(background),
                     'Background_BH_eQTL_SNPs':len(background&bh),
                     'Background_BH_Percent':100*len(background&bh)/len(background),**pair_counts[disease]})
    pd.DataFrame(rows).to_csv(out/'disease_eqtl_overlap.tsv',sep='\t',index=False)


def figures():
    out=REPO_ROOT/'results/part2/exploration';figs=out/'figures';figs.mkdir(exist_ok=True)
    stability=pd.read_csv(out/'threshold_rank_stability.tsv',sep='\t')
    fig,axes=plt.subplots(1,3,figsize=(13,4),squeeze=False)
    for ax,(disease,g) in zip(axes[0],stability.groupby('Disease',sort=False)):
        matrix=np.eye(5)
        for row in g.itertuples():
            a,b=THRESHOLDS.index(row.Threshold_A),THRESHOLDS.index(row.Threshold_B)
            matrix[a,b]=matrix[b,a]=row.Spearman_R
        im=ax.imshow(matrix,vmin=-1,vmax=1,cmap='RdBu_r')
        labels=[f'{p:g}' for p in THRESHOLDS]
        ax.set(xticks=range(5),yticks=range(5),xticklabels=labels,yticklabels=labels,title=disease,
               xlabel='Fixed p-value threshold',ylabel='Fixed p-value threshold')
        for i in range(5):
            for j in range(5):ax.text(j,i,f'{matrix[i,j]:.2f}',ha='center',va='center',fontsize=8,
                                      color='white' if abs(matrix[i,j])>.55 else 'black')
    fig.subplots_adjust(left=.06,right=.85,bottom=.16,top=.87,wspace=.45)
    color_axis=fig.add_axes([.90,.22,.012,.58])
    fig.colorbar(im,cax=color_axis,label='Spearman reference rank correlation')
    fig.savefig(figs/'threshold_rank_stability.png',dpi=320);fig.savefig(figs/'threshold_rank_stability.pdf');plt.close(fig)
    overlap=pd.read_csv(out/'disease_eqtl_overlap.tsv',sep='\t')
    fig,ax=plt.subplots(figsize=(8,4));x=np.arange(len(overlap))
    ax.bar(x-.18,overlap.Score_BH_Percent_Of_Tested,.36,label='Selected PRS SNPs',color='#287d8e')
    ax.bar(x+.18,overlap.Background_BH_Percent,.36,label='Harmonized GWAS/tested background',color='#9aa9b4')
    ax.set(xticks=x,xticklabels=overlap.Disease,ylabel='SNPs with ≥1 BH-significant eQTL (%)',
           title='Exploratory disease-score SNP / cis-eQTL overlap (GRCh37)')
    ax.legend(frameon=False);save(fig,figs/'disease_eqtl_overlap')


def main():
    threshold_sensitivity();eqtl_overlap();figures()
    print('Saved public threshold sensitivity and descriptive variant overlap. Primary PRSs unchanged.')

if __name__=='__main__':main()
