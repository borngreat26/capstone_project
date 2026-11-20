"""Autosomal LD-pruned PCA of the public EUR reference panel; no private genotypes."""
import argparse,sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from src.utils import REPO_ROOT,run_plink

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--rebuild',action='store_true');args=ap.parse_args()
    work=REPO_ROOT/'data/part2_work/pca';work.mkdir(parents=True,exist_ok=True)
    out=REPO_ROOT/'results/part2/pca';out.mkdir(parents=True,exist_ok=True)
    ps=[];expected=None
    for chrom in range(1,23):
        source=REPO_ROOT/f'data/LDREF/1000G.EUR.{chrom}'
        fam=pd.read_csv(f'{source}.fam',sep=r'\s+',header=None,dtype=str)[[0,1]]
        if expected is None: expected=fam
        if not fam.equals(expected):raise ValueError('Reference sample identities/order differ across chromosomes')
        prefix=work/f'chr{chrom}'
        if args.rebuild or not Path(f'{prefix}.prune.in').exists():
            run_plink(['--bfile',source,'--autosome','--maf','.05','--geno','.02','--hwe','1e-6','--indep-pairwise','200','50','.2','--out',prefix])
        if args.rebuild or not Path(f'{prefix}.pgen').exists():
            run_plink(['--bfile',source,'--extract',f'{prefix}.prune.in','--make-pgen','--out',prefix])
        ps.append(str(prefix))
    merge=work/'merge.list';merge.write_text('\n'.join(ps)+'\n');combined=work/'autosomes'
    if args.rebuild or not Path(f'{combined}.pgen').exists():
        run_plink(['--pmerge-list',merge,'--make-pgen','--out',combined])
    # Request all sample eigenvalues: their sum, rather than only the top ten,
    # is the denominator for explained variance in PLINK's variance-standardized GRM.
    prefix=work/'pca'
    run_plink(['--pfile',combined,'--pca',str(len(expected)-1),'--out',prefix])
    pcs=pd.read_csv(f'{prefix}.eigenvec',sep=r'\s+',dtype={'IID':str})
    eig=np.loadtxt(f'{prefix}.eigenval');var=pd.DataFrame({'PC':np.arange(1,len(eig)+1),'Eigenvalue':eig,'Explained_Fraction':eig/eig.sum()})
    pcs=pcs[[c for c in pcs.columns if not c.startswith('PC') or int(c[2:])<=10]]
    pcs.to_csv(out/'reference_pcs.tsv',sep='\t',index=False);var.to_csv(out/'explained_variance.tsv',sep='\t',index=False)
    counts={str(c):len(Path(f'{p}.prune.in').read_text().splitlines()) for c,p in zip(range(1,23),ps)}
    pd.Series(counts,name='LD_Pruned_SNPs').rename_axis('Chromosome').to_csv(out/'variant_counts.tsv',sep='\t')
    fig,ax=plt.subplots(figsize=(6,5));ax.scatter(pcs.PC1,pcs.PC2,s=14,alpha=.7,color='#287d8e')
    ax.set(title=f'Within-EUR reference structure (n={len(pcs)})',xlabel=f'PC1 ({100*var.Explained_Fraction.iloc[0]:.2f}%)',ylabel=f'PC2 ({100*var.Explained_Fraction.iloc[1]:.2f}%)')
    fig.tight_layout();fig.savefig(out/'reference_PCA.png',dpi=320);fig.savefig(out/'reference_PCA.pdf');plt.close(fig)
    print(f'{len(pcs)} public reference people, {sum(counts.values())} LD-pruned SNPs');print(var.head().to_string(index=False))
if __name__=='__main__':main()
