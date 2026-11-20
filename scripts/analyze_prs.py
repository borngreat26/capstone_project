"""Explore public disease scores on the same people; PC adjustment is a sensitivity analysis."""
import itertools,sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr,spearmanr
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from src.utils import REPO_ROOT
from scripts.run_week7 import save


def heatmap(matrix,title,path,vmin=-1,vmax=1):
    fig,ax=plt.subplots(figsize=(6,5));im=ax.imshow(matrix,vmin=vmin,vmax=vmax,cmap='RdBu_r')
    ax.set(xticks=range(len(matrix)),yticks=range(len(matrix)),xticklabels=matrix.columns,yticklabels=matrix.index,title=title)
    plt.setp(ax.get_xticklabels(),rotation=30,ha='right')
    for i in range(len(matrix)):
        for j in range(len(matrix)):ax.text(j,i,f'{matrix.iloc[i,j]:.3f}',ha='center',va='center',color='white' if abs(matrix.iloc[i,j])>.55 else 'black')
    fig.colorbar(im,ax=ax,shrink=.75);save(fig,path)

def main():
    root=REPO_ROOT/'results/part2';figs=root/'figures';figs.mkdir(exist_ok=True)
    z=pd.read_csv(root/'reference_standardized_prs.tsv',sep='\t',index_col='IID')
    if z.isna().any().any() or not z.index.is_unique:raise ValueError('Scores must be complete with unique people')
    selected=pd.read_csv(root/'selected_diseases.tsv',sep='\t');summary=pd.read_csv(root/'prs_summary.tsv',sep='\t').set_index('Disease')
    sets={};raw={}
    for folder in root.iterdir():
        if (folder/'weights.tsv').exists():
            weights=pd.read_csv(folder/'weights.tsv',sep='\t');ref=pd.read_csv(folder/'reference_scores.tsv',sep='\t',index_col='IID').PRS
            name=selected.set_index('Slug').loc[folder.name,'Disease']
            sets[name]=set(weights.SNP);raw[name]=ref
    pcs=pd.read_csv(root/'pca/reference_pcs.tsv',sep='\t',dtype={'IID':str}).set_index('IID')
    if not set(z.index)<=set(pcs.index):raise ValueError('PCA sample mismatch')
    x=np.column_stack([np.ones(len(z)),pcs.loc[z.index,['PC1','PC2']].to_numpy()])
    residual=pd.DataFrame(z.to_numpy()-x@np.linalg.lstsq(x,z.to_numpy(),rcond=None)[0],index=z.index,columns=z.columns)
    residual.to_csv(root/'pc_adjusted_residuals.tsv',sep='\t')
    pear=z.corr(method='pearson');spearman=z.corr(method='spearman');adj=residual.corr()
    for name,m in [('pearson',pear),('spearman',spearman),('pc_adjusted_pearson',adj)]:m.to_csv(root/f'{name}_correlations.tsv',sep='\t')
    pairs=[];overlaps=[];jac=pd.DataFrame(np.eye(len(z.columns)),index=z.columns,columns=z.columns)
    for a,b in itertools.combinations(z.columns,2):
        shared=len(sets[a]&sets[b]);union=len(sets[a]|sets[b]);n1,n2=len(sets[a]),len(sets[b])
        pairs.append(dict(Disease_A=a,Disease_B=b,Individuals=len(z),SNPs_A=n1,SNPs_B=n2,Pearson_R=pear.loc[a,b],Pearson_P=pearsonr(z[a],z[b]).pvalue,Spearman_R=spearman.loc[a,b],Spearman_P=spearmanr(z[a],z[b]).pvalue,PC1_PC2_Adjusted_R=adj.loc[a,b]))
        overlaps.append(dict(Disease_A=a,Disease_B=b,SNPs_A=n1,SNPs_B=n2,Shared=shared,Fraction_A=shared/n1,Fraction_B=shared/n2,Jaccard=shared/union))
        jac.loc[a,b]=jac.loc[b,a]=shared/union
    pd.DataFrame(pairs).to_csv(root/'pairwise_correlations.tsv',sep='\t',index=False);pd.DataFrame(overlaps).to_csv(root/'snp_overlap.tsv',sep='\t',index=False)
    pcrows=[]
    for disease in z:
        for pc in ['PC1','PC2']:
            r,p=pearsonr(z[disease],pcs.loc[z.index,pc]);pcrows.append(dict(Disease=disease,PC=pc,Individuals=len(z),Pearson_R=r,P=p))
    pd.DataFrame(pcrows).to_csv(root/'prs_pc_associations.tsv',sep='\t',index=False)
    plt.rcParams.update({'font.size':10,'pdf.fonttype':42})
    heatmap(pear,'Disease score correlations (Pearson)',figs/'prs_correlations')
    heatmap(adj,'Exploratory correlations adjusted for PC1–PC2',figs/'prs_pc_adjusted_correlations')
    heatmap(jac,'Selected SNP overlap (Jaccard)',figs/'snp_overlap',0,1)
    fig,axes=plt.subplots(1,len(z.columns),figsize=(4*len(z.columns),3.7),squeeze=False)
    for ax,disease in zip(axes[0],z.columns):
        ax.hist(raw[disease],bins=24,color='#287d8e',edgecolor='white');ax.set(title=f'{disease}\n{len(sets[disease])} SNPs; n={len(z)}',xlabel='Raw score (log-odds weighted sum)',ylabel='Reference individuals')
    save(fig,figs/'prs_distributions')
    fig,ax=plt.subplots(figsize=(7,4));ax.boxplot([z[c] for c in z],tick_labels=z.columns,showfliers=True);ax.set(ylabel='PRS z-score (reference mean and sample SD)',title='Standardized public reference scores');save(fig,figs/'standardized_prs')
    if pairs:
        fig,axes=plt.subplots(1,len(pairs),figsize=(5*len(pairs),4),squeeze=False)
        for ax,row in zip(axes[0],pairs):
            a,b=row['Disease_A'],row['Disease_B'];ax.scatter(z[a],z[b],s=12,alpha=.5,color='#287d8e');ax.set(xlabel=a+' (z)',ylabel=b+' (z)',title=f"Pearson r={row['Pearson_R']:.3f}; n={len(z)}")
        save(fig,figs/'pairwise_scatterplots')
    print(pd.DataFrame(pairs).to_string(index=False));print(pd.DataFrame(overlaps).to_string(index=False))
if __name__=='__main__':main()
