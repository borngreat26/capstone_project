"""Annotate saved Part 1 associations; never rerun eQTL regressions.

Requires the original ENCODE BED, official JSON metadata, UCSC chain and liftOver.
"""
import argparse,gzip,hashlib,json,subprocess,sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from src.annotation import read_peaks,snp_counts,interval_counts,matching_peaks
from src.utils import REPO_ROOT


def save(fig,path):
    fig.tight_layout()
    fig.savefig(path.with_suffix('.png'),dpi=320)
    fig.savefig(path.with_suffix('.pdf'))
    plt.close(fig)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--liftover',default='liftOver');args=ap.parse_args()
    data=REPO_ROOT/'data/encode';out=REPO_ROOT/'results/week7';figs=out/'figures'
    figs.mkdir(parents=True,exist_ok=True)
    meta=json.loads((data/'ENCFF421XIL.metadata.json').read_text())
    exp=json.loads((data/'ENCSR603LVR.metadata.json').read_text())
    assert meta['accession']=='ENCFF421XIL' and meta['dataset']=='/experiments/ENCSR603LVR/'
    assert meta['file_format_type']=='narrowPeak' and exp['assay_term_name']=='ATAC-seq'
    original=data/'ENCFF421XIL.bed'
    content_md5=hashlib.md5(original.read_bytes()).hexdigest()
    if content_md5 != meta['content_md5sum']: raise ValueError('ENCODE uncompressed content MD5 mismatch')
    peaks=read_peaks(original)
    before=len(peaks)
    if meta['assembly']=='GRCh38':
        # Lift BED4 identifiers only. Scores/annotations join back unchanged; original summit offsets
        # remain annotations in the original assembly and are NOT used as lifted summit coordinates.
        bed=data/'ENCFF421XIL.lift_input.bed'
        lift=pd.DataFrame({'c':'chr'+peaks.Chrom,'s':peaks.Start,'e':peaks.End,'id':np.arange(len(peaks))})
        lift.to_csv(bed,sep='\t',header=False,index=False)
        dest=data/'ENCFF421XIL.GRCh37.bed4';unmapped=data/'ENCFF421XIL.unmapped.bed'
        subprocess.run([args.liftover,'-minMatch=0.95',str(bed),str(data/'hg38ToHg19.over.chain.gz'),str(dest),str(unmapped)],check=True)
        mapped=pd.read_csv(dest,sep='\t',header=None,names=['Chrom','Start','End','id'])
        if mapped.id.duplicated().any(): raise ValueError('Nonunique lifted peak')
        peaks=peaks.iloc[mapped.id].reset_index(drop=True).copy()
        peaks[['Original_Chrom','Original_Start','Original_End']]=peaks[['Chrom','Start','End']].to_numpy()
        peaks[['Chrom','Start','End']]=mapped[['Chrom','Start','End']].to_numpy()
        peaks['Chrom']=peaks.Chrom.str.removeprefix('chr')
        if 'SummitOffset' in peaks: peaks=peaks.rename(columns={'SummitOffset':'Original_SummitOffset'})
    elif meta['assembly'] not in ('GRCh37','hg19'): raise ValueError('Unsupported/unverified assembly')
    peaks.to_csv(out/'peaks_GRCh37.tsv.gz',sep='\t',index=False)
    provenance=dict(Experiment=exp['accession'],File=meta['accession'],Biosample=meta['biosample_ontology']['term_name'],
        Assay=exp['assay_term_name'],Format='BED6+4 narrowPeak',Original_Assembly=meta['assembly'],Analysis_Assembly='GRCh37',
        URL='https://www.encodeproject.org'+meta['href'],Official_Compressed_MD5=meta['md5sum'],Verified_Content_MD5=content_md5,
        Local_BED_SHA256=hashlib.sha256(original.read_bytes()).hexdigest(),Original_Download_Date='Unknown: file existed before this work',
        Metadata_Retrieved='2026-11-17',Input_Peaks=before,Mapped_Peaks=len(peaks),Lost_Peaks=before-len(peaks),
        Method='UCSC liftOver hg38ToHg19; -minMatch=0.95, unique mapping; original annotations retained',
        Chain_SHA256=hashlib.sha256((data/'hg38ToHg19.over.chain.gz').read_bytes()).hexdigest())
    (out/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    # GENCODE v12 gene entries are one-based inclusive: convert [start,end] to [start-1,end).
    genes=[]
    with gzip.open(REPO_ROOT/'data/quarter1_results/part1/gencode.v12.annotation.gtf.gz','rt') as handle:
        for line in handle:
            if line.startswith('#'): continue
            f=line.rstrip().split('\t')
            if f[2]!='gene' or f[0].removeprefix('chr') not in list(map(str,range(1,23))): continue
            attr=dict(x.strip().split(' ',1) for x in f[8].strip(';').split(';') if x.strip())
            attr={k:v.strip('"') for k,v in attr.items()}
            if attr.get('gene_type')!='protein_coding': continue
            genes.append([attr['gene_id'],attr['gene_name'],f[0].removeprefix('chr'),int(f[3]),int(f[4])])
    genes=pd.DataFrame(genes,columns=['Gene_ID','Gene_Symbol','Chrom','Start_1based','End_1based'])
    for chrom,group in genes.groupby('Chrom'):
        genes.loc[group.index,'Body_Peaks']=interval_counts(peaks,chrom,group.Start_1based-1,group.End_1based)
        genes.loc[group.index,'Cis_Window_Peaks']=interval_counts(peaks,chrom,np.maximum(0,group.Start_1based-1-500000),group.End_1based+500000)
    for c in ['Body_Peaks','Cis_Window_Peaks']: genes[c]=genes[c].astype(int)
    genes.to_csv(out/'gene_annotation.tsv',sep='\t',index=False)
    leads=pd.read_csv(REPO_ROOT/'results/part1/table_B_gene_summary.tsv',sep='\t').head(4)
    snp_sets={key:{} for key in ['All tested','BH significant','Bonferroni significant']}
    pair_counts={key:[0,0] for key in snp_sets};loci=[];total=0
    source=REPO_ROOT/'data/quarter1_results/part1/table_A_genome_wide_cis_eqtl.tsv.gz'
    overlap_path=out/'overlapping_associations.tsv.gz'
    with gzip.open(overlap_path,'wt') as overlap_handle:
        first=True
        for chunk in pd.read_csv(source,sep='\t',chunksize=250000,dtype={'SNP_Chr':str}):
            total+=len(chunk);counts=np.zeros(len(chunk),dtype=int)
            for chrom,g in chunk.groupby('SNP_Chr'):
                counts[chunk.index.get_indexer(g.index)]=snp_counts(peaks,chrom,g.SNP_Position)
            chunk['ATAC_Peaks']=counts
            hit=chunk.loc[chunk.ATAC_Peaks.gt(0)]
            hit.to_csv(overlap_handle,sep='\t',index=False,header=first);first=False
            for key,mask in [('All tested',np.ones(len(chunk),dtype=bool)),('BH significant',chunk.P_BH.lt(.05)),('Bonferroni significant',chunk.P_Bonferroni.lt(.05))]:
                sub=chunk.loc[mask];pair_counts[key][0]+=len(sub);pair_counts[key][1]+=int(sub.ATAC_Peaks.gt(0).sum())
                # Identity includes SNP ID and genomic location to expose inconsistencies.
                for r in sub[['SNP_ID','SNP_Chr','SNP_Position','ATAC_Peaks']].drop_duplicates().itertuples(index=False,name=None):
                    snp_sets[key][r[:3]]=r[3]
            loci.append(chunk.loc[chunk.Gene_ID.isin(leads.Gene_ID)])
    assert total==6597508, 'Unexpected Table A size: verify the Part 1 source'
    summary=[]
    for key,values in snp_sets.items():
        hits=sum(v>0 for v in values.values())
        summary.append(dict(Set=key,Unique_SNPs=len(values),Overlapping_SNPs=hits,Percent=100*hits/len(values),
            SNP_Gene_Pairs=pair_counts[key][0],Overlapping_Pairs=pair_counts[key][1]))
    summary=pd.DataFrame(summary);summary.to_csv(out/'snp_overlap_summary.tsv',sep='\t',index=False)
    # Per-SNP/peak links carry actual lifted coordinates and preserved accessibility scores.
    links=[]
    bins={}
    for peak in peaks.itertuples():
        for b in range(int(peak.Start)//10000,(int(peak.End)-1)//10000+1):
            bins.setdefault((peak.Chrom,b),[]).append(peak)
    for (snp,chrom,pos),count in snp_sets['All tested'].items():
        if count:
            for p in bins.get((chrom,(pos-1)//10000),[]):
                if not p.Start < pos or not p.End > pos-1: continue
                links.append([snp,chrom,pos,p.Start,p.End,p.Name,p.Score,p.Signal])
    pd.DataFrame(links,columns=['SNP_ID','Chrom','SNP_Position','Peak_Start_0based','Peak_End_exclusive','Peak_Name','Score','Signal']).to_csv(out/'snp_peak_links.tsv.gz',sep='\t',index=False)
    loci=pd.concat(loci);cases=[]
    for gene,g in loci.groupby('Gene_ID',sort=False):
        lead=g.loc[g.P_Value.idxmin()];ann=genes.set_index('Gene_ID').loc[gene]
        cases.append(dict(Gene_ID=gene,Gene=lead.Gene_Symbol,Lead_SNP=lead.SNP_ID,Lead_P=lead.P_Value,
            Lead_Overlaps=bool(lead.ATAC_Peaks>0),BH_SNPs_Overlapping=int(g.loc[g.P_BH.lt(.05)&g.ATAC_Peaks.gt(0),'SNP_ID'].nunique()),
            Bonferroni_SNPs_Overlapping=int(g.loc[g.P_Bonferroni.lt(.05)&g.ATAC_Peaks.gt(0),'SNP_ID'].nunique()),
            Body_Peaks=ann.Body_Peaks,Cis_Window_Peaks=ann.Cis_Window_Peaks))
    pd.DataFrame(cases).to_csv(out/'case_studies.tsv',sep='\t',index=False)
    plt.rcParams.update({'font.size':10,'pdf.fonttype':42})
    fig,ax=plt.subplots(figsize=(6.6,4));bars=ax.bar(summary.Set,summary.Percent,color=['#6c8596','#287d8e','#bf6540'])
    ax.bar_label(bars,fmt='%.2f%%',padding=3);ax.set(ylabel='Unique SNPs overlapping peaks (%)',title='B-cell ATAC-seq overlap (GRCh37)');ax.set_ylim(0,summary.Percent.max()*1.2)
    save(fig,figs/'A_snp_overlap')
    r=loci.loc[loci.Gene_Symbol.eq('RPS26')].sort_values('SNP_Position');lead=r.loc[r.P_Value.idxmin()]
    fig,axes=plt.subplots(3,1,figsize=(8,6),sharex=True,gridspec_kw={'height_ratios':[4,1,1]})
    for overlap,color,label in [(False,'#acb3b8','No peak overlap'),(True,'#287d8e','ATAC peak overlap')]:
        points=r.loc[r.ATAC_Peaks.gt(0).eq(overlap)]
        axes[0].scatter(points.SNP_Position/1e6,-np.log10(points.P_Value),s=15,color=color,label=label)
    axes[0].scatter([lead.SNP_Position/1e6],[-np.log10(lead.P_Value)],marker='*',s=130,c='#bf6540',label=lead.SNP_ID)
    axes[0].set(ylabel='−log10(p)',title='RPS26 regional cis-eQTL associations');axes[0].legend(frameon=False,fontsize=9)
    for p in matching_peaks(peaks,str(lead.SNP_Chr),r.SNP_Position.min()-1,r.SNP_Position.max()).itertuples():
        axes[1].plot([p.Start/1e6,p.End/1e6],[.5,.5],lw=6,color='#287d8e')
    axes[1].set(ylabel='ATAC',ylim=(0,1),yticks=[])
    axes[2].plot([lead.Gene_Start/1e6,lead.Gene_End/1e6],[.5,.5],lw=7,color='#bf6540');axes[2].set(ylabel='RPS26',ylim=(0,1),yticks=[],xlabel=f'Chromosome {lead.SNP_Chr} position (Mb, GRCh37)')
    save(fig,figs/'B_RPS26_accessibility')
    fig,ax=plt.subplots(figsize=(6,4));values=[int(genes.Body_Peaks.gt(0).sum()),int(genes.Cis_Window_Peaks.gt(0).sum())]
    bars=ax.bar(['Gene body','Gene body ±500 kb'],values,color=['#287d8e','#6c8596']);ax.bar_label(bars,padding=3);ax.set(ylabel='Protein-coding genes with ≥1 peak',title=f'GENCODE v12 autosomal genes (n={len(genes):,})',ylim=(0,len(genes)*1.15))
    save(fig,figs/'C_gene_overlap')
    print(summary.to_string(index=False));print(pd.DataFrame(cases).to_string(index=False));print(json.dumps(provenance,indent=2))
if __name__=='__main__':main()
