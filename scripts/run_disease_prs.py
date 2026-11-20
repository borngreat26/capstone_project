"""Fixed-threshold public disease PRSs from verified GRCh37 original effects.

Reuses Week 5 P+T functions. Private comparisons use a separate shared SNP subset
and write only to ignored data/private_part2/. No Z-to-beta conversions.
"""
import argparse,hashlib,json,sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from src.prs import read_bim,gwas_qc,harmonize_to_reference,clump,select_weights,score_bfile,relative_position
from src.utils import REPO_ROOT,run_plink

STUDIES={
 'schizophrenia':dict(name='Schizophrenia',file='GCST006803.txt.gz',accession='CLOZUK+PGC2 2018',effect='OR',format='walters',n=105318,source='Pardiñas et al. 2018; PMID 29483656',reason='Author-documented A1 effect allele, A2 other allele and ORs; European schizophrenia meta-analysis adds a psychiatric disease. GRCh37 coordinates validated against independent panel rsID anchors; original archive MD5 verified. Noncommercial use; original effects are not redistributed.'),
 'rheumatoid_arthritis':dict(name='Rheumatoid arthritis',file='GCST005569.Build37.f.tsv.gz',accession='GCST005569',effect='OR',n=47580,source='Eyre et al. 2012; PMID 23143596',reason='Documented effect/other alleles, ORs, GRCh37 positions; immune disease. Immunochip ascertainment limits genome-wide coverage.'),
 'type2_diabetes':dict(name='Type 2 diabetes',file='GCST007517.Build37.f.tsv.gz',accession='GCST007517',effect='BETA',n=298957,source='Mahajan et al. 2018; PMID 29632382',reason='European BMI-unadjusted exome-array analysis with log-odds betas and GRCh37 positions; metabolic disease. Targeted array limits genome-wide coverage.'),
}

def load_study(study,panel):
    path=REPO_ROOT/'data/gwas_original'/study['file']
    original_count=None;prefilter_removed=0
    if study.get('format')=='walters':
        expected=Path(str(path)+'.md5')
        md5path=path.parent/'GCST006803.md5'
        if hashlib.md5(path.read_bytes()).hexdigest()!=md5path.read_text().split()[0]:
            raise ValueError('Author schizophrenia archive MD5 mismatch')
        pieces=[];original_count=0
        panel_ids=panel.set_index('SNP')[['CHR','BP']]
        checked=matched=0
        for chunk in pd.read_csv(path,sep='\t',chunksize=200000):
            original_count+=len(chunk)
            chunk=chunk.rename(columns={'OR':'odds_ratio'})
            chunk['SNP']=chunk.SNP.str.extract(r'^(rs\d+)(?::|$)',expand=False)
            # Validate coordinates against independently identified rsIDs in the GRCh37 panel.
            joined=chunk.merge(panel_ids,left_on='SNP',right_index=True,suffixes=('','_REF'))
            checked+=len(joined);matched+=int((joined.CHR.eq(joined.CHR_REF)&joined.BP.eq(joined.BP_REF)).sum())
            chunk=chunk.loc[chunk.SNP.isin(panel_ids.index)]
            pieces.append(chunk)
        if checked<10000 or matched/checked<.99:
            raise ValueError('Schizophrenia coordinates do not validate against GRCh37 rsID anchors')
        raw=pd.concat(pieces,ignore_index=True);prefilter_removed=original_count-len(raw)
        raw=raw.rename(columns={'CHR':'chromosome','BP':'base_pair_location','SNP':'variant_id','A1':'effect_allele','A2':'other_allele','P':'p_value'})
    else:
        raw=pd.read_csv(path,sep='\t',dtype={'chromosome':str},low_memory=False)
    d=raw.rename(columns={'chromosome':'CHR','base_pair_location':'BP','variant_id':'SNP','effect_allele':'A1','other_allele':'A2','p_value':'P','beta':'BETA'}).copy()
    d['CHR']=d.CHR.astype(str).str.strip().str.removeprefix('chr')
    d['CHR']=pd.to_numeric(d.CHR,errors='coerce');d['BP']=pd.to_numeric(d.BP,errors='coerce')
    for c in ['A1','A2']:d[c]=d[c].astype(str).str.upper()
    if study['effect']=='OR':
        d['BETA']=np.log(pd.to_numeric(raw.odds_ratio,errors='coerce').where(raw.odds_ratio.gt(0)))
    # Assign rsIDs only from unique GRCh37 panel positions; never invent an allele or orientation.
    position=panel.drop_duplicates(['CHR','BP'],keep=False)[['CHR','BP','SNP']].rename(columns={'SNP':'Panel_SNP'})
    d=d.merge(position,on=['CHR','BP'],how='left',validate='many_to_one')
    d['SNP']=d.SNP.astype('object')
    missing=d.SNP.isna() | d.SNP.astype(str).eq('.')
    d.loc[missing,'SNP']=d.loc[missing,'Panel_SNP']
    counts={'Original_GWAS':original_count if original_count is not None else len(raw),'Nonpanel_Rows_Prefiltered':prefilter_removed,'Missing_ID_After_Position_Match':int(d.SNP.isna().sum())}
    d=d.loc[d.SNP.notna()].copy()
    valid=d.CHR.between(1,22)&d.BP.ge(1)&np.isfinite(pd.to_numeric(d.BETA,errors='coerce'))
    counts['Invalid_Position_Or_Effect']=int((~valid).sum());d=d.loc[valid]
    dup=d.SNP.duplicated(keep=False);counts['Duplicate_ID_Rows_Removed']=int(dup.sum());d=d.loc[~dup]
    qc=gwas_qc(d);counts['Allele_Pvalue_QC_Removed']=len(d)-len(qc)
    h,hcounts=harmonize_to_reference(qc,panel);counts.update(hcounts)
    if study.get('format')=='walters': counts.update(Build_Anchor_rsIDs=checked,GRCh37_Anchor_Matches=matched)
    return h,counts

def score_by_chrom(weights,work,private_bfile=None):
    total=None
    for chrom,g in weights.groupby('CHR'):
        c=int(chrom);bfile=private_bfile or REPO_ROOT/f'data/LDREF/1000G.EUR.{c}'
        scores=score_bfile(bfile,g,work/f'score_chr{c}')
        if not scores.index.is_unique:raise ValueError('Nonunique person IDs')
        if total is None:total=scores
        else:
            if set(total.index)!=set(scores.index):raise ValueError('Reference individuals differ across chromosomes')
            total=total+scores.reindex(total.index)
    return total

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--skip-private',action='store_true');args=ap.parse_args()
    out=REPO_ROOT/'results/part2';work=REPO_ROOT/'data/part2_work';out.mkdir(parents=True,exist_ok=True);work.mkdir(parents=True,exist_ok=True)
    panel=pd.concat([read_bim(REPO_ROOT/f'data/LDREF/1000G.EUR.{c}') for c in range(1,23)],ignore_index=True)
    if panel.SNP.duplicated().any():raise ValueError('Panel rsIDs not unique')
    summaries=[];selected=[];scores=[];private_results=[];private_panels=[]
    for slug,study in sorted(STUDIES.items(), key=lambda item: item[0]=='schizophrenia'):
        wd=work/slug;wd.mkdir(exist_ok=True);dest=out/slug;dest.mkdir(exist_ok=True)
        h,counts=load_study(study,panel);weights=[];nclump=0
        for chrom,g in h.groupby('CHR'):
            # SNPs with P>0.001 cannot become index SNPs passing threshold; omit before clumping.
            eligible=g.loc[g.P.le(.001)]
            if eligible.empty:continue
            ids=clump(eligible,REPO_ROOT/f'data/LDREF/1000G.EUR.{int(chrom)}',wd/f'chr{int(chrom)}',r2=.1,kb=250)
            nclump+=len(ids);weights.append(select_weights(eligible,ids,.001))
        if not weights:raise ValueError(f'No eligible weights for {slug}')
        weights=pd.concat(weights,ignore_index=True)
        ref=score_by_chrom(weights,wd)
        weights.to_csv(dest/'weights.tsv',sep='\t',index=False)
        h.to_csv(wd/'harmonized.tsv.gz',sep='\t',index=False)
        counts.update(After_Clumping_And_Threshold=len(weights),Successfully_Scored=len(weights))
        pd.Series(counts,name='SNPs').rename_axis('Stage').to_csv(dest/'snp_counts.tsv',sep='\t')
        ref.rename('PRS').rename_axis('IID').to_csv(dest/'reference_scores.tsv',sep='\t')
        z=(ref-ref.mean())/ref.std(ddof=1)
        scores.append(z.rename(study['name']))
        summaries.append(dict(Disease=study['name'],Original_SNPs=counts['Original_GWAS'],Harmonized_SNPs=len(h),Scored_SNPs=len(weights),
            Individuals=len(ref),Mean=ref.mean(),SD=ref.std(ddof=1),Median=ref.median(),Q1=ref.quantile(.25),Q3=ref.quantile(.75),IQR=ref.quantile(.75)-ref.quantile(.25),P_Threshold=.001,Clump_R2=.1,Clump_KB=250,Assembly='GRCh37'))
        selected.append(dict(Slug=slug,Disease=study['name'],Accession=study['accession'],Source=study['source'],File=study['file'],Effect='log(OR)' if study['effect']=='OR' else 'log-odds beta',N=study['n'],Ancestry='European',Assembly='GRCh37: verified panel rsID anchors' if study.get('format')=='walters' else 'GRCh37 (formatted file; not GRCh38 .h file)',Justification=study['reason'],SHA256=hashlib.sha256((REPO_ROOT/'data/gwas_original'/study['file']).read_bytes()).hexdigest()))
        # Existing converted file stays untouched. Remove missing calls in a new local-only copy.
        private=REPO_ROOT/'data/week5_results/Tiffany_converted'
        if not args.skip_private and Path(f'{private}.bed').exists():
            pdest=REPO_ROOT/'data/private_part2'/slug;pdest.mkdir(parents=True,exist_ok=True)
            clean=pdest/'complete_calls'
            run_plink(['--bfile',private,'--geno','0','--make-bed','--out',clean])
            pv=read_bim(clean);common,_=harmonize_to_reference(weights[['CHR','BP','SNP','A1','A2','BETA','P']],pv)
            if len(common):
                public_work=pdest/'reference';public_work.mkdir(exist_ok=True)
                matched=score_by_chrom(common,public_work);target=score_bfile(clean,common,pdest/'target')
                if len(target)!=1:raise ValueError('Expected one private target')
                position=relative_position(matched,float(target.iloc[0]));position['Shared_Scoring_SNPs']=len(common)
                private_results.append({'Disease':study['name'],**position})
                private_panels.append((study['name'],matched,float(target.iloc[0]),len(common)))
                pd.Series(position,name='Value').rename_axis('Quantity').to_csv(pdest/'summary.tsv',sep='\t')
                from src.prs import plot_prs_distribution
                plot_prs_distribution(matched,float(target.iloc[0]),study['name']+' — common-SNP comparison','Private target',pdest/'distribution.png')
        print(study['name'],f'{len(h)} harmonized; {len(weights)} scored SNPs; {len(ref)} public people')
    pd.DataFrame(summaries).to_csv(out/'prs_summary.tsv',sep='\t',index=False)
    pd.DataFrame(selected).to_csv(out/'selected_diseases.tsv',sep='\t',index=False)
    pd.concat(scores,axis=1,join='inner').rename_axis('IID').to_csv(out/'reference_standardized_prs.tsv',sep='\t')
    if private_results:
        private_root=REPO_ROOT/'data/private_part2'
        pd.DataFrame(private_results).to_csv(private_root/'combined_summary.tsv',sep='\t',index=False)
        import matplotlib.pyplot as plt
        fig,axes=plt.subplots(1,len(private_panels),figsize=(4*len(private_panels),4),squeeze=False)
        for ax,(name,reference,target,n) in zip(axes[0],private_panels):
            ax.hist(reference,bins=24,color='#287d8e',edgecolor='white')
            ax.axvline(target,color='#bf6540',label='Private target')
            ax.set(title=f'{name}\n{n} common SNPs',xlabel='Raw common-SNP score',ylabel='Reference individuals');ax.legend(frameon=False)
        fig.tight_layout();fig.savefig(private_root/'private_prs_distributions.png',dpi=320);fig.savefig(private_root/'private_prs_distributions.pdf');plt.close(fig)
        ranked=pd.DataFrame(private_results).sort_values('Target_Z_Score')
        fig,ax=plt.subplots(figsize=(7,4));ax.barh(ranked.Disease,ranked.Target_Z_Score,color='#287d8e');ax.set(xlabel='Private target reference-standardized score',title='Common-SNP scores; relative position only')
        fig.tight_layout();fig.savefig(private_root/'private_standardized_comparison.png',dpi=320);fig.savefig(private_root/'private_standardized_comparison.pdf');plt.close(fig)

if __name__=='__main__':main()
