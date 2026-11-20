"""Reinspect local TCSC files without deriving PRS weights from Z-scores."""
import sys
from pathlib import Path
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.utils import REPO_ROOT

def main():
    rows=[]
    output=REPO_ROOT/'results/part2/gwas_inventory.tsv'
    previous=pd.read_csv(output,sep='\t') if output.exists() else pd.DataFrame()
    cleanup_path=REPO_ROOT/'results/part2/unused_gwas_cleanup.tsv'
    removed=set(pd.read_csv(cleanup_path,sep='\t').Filename) if cleanup_path.exists() else set()
    for path in sorted((REPO_ROOT/'data/sumstats').glob('*.sumstats.gz')):
        n=0; ns=[]; columns=None
        for chunk in pd.read_csv(path, sep=r'\s+', chunksize=250000):
            columns=list(chunk.columns); n+=len(chunk)
            if 'N' in chunk: ns.extend([chunk.N.min(),chunk.N.max()])
        trait=path.name.removesuffix('.sumstats.gz').removeprefix('PASS_')
        rows.append(dict(Filename=path.name,Phenotype=trait.replace('_',' '),
            Original_Source='Author/year suggested by filename only; not independently identified' if any(c.isdigit() for c in trait) else 'Not documented in local file',
            Variants=n,Columns=','.join(columns),Genome_Assembly='Not documented; rsIDs alone do not verify build',
            Effect_Allele='A1 convention unverified against original study',Other_Allele='A2',
            Effect_Size='Z-score only' if 'Z' in columns and not {'BETA','OR'}&set(columns) else 'Inspect original source',
            P_Available='P' in columns,N_Min=min(ns) if ns else None,N_Max=max(ns) if ns else None,
            Ancestry='Not documented in local file',Direct_PRS_Eligible=False,
            Reason='No documented beta or OR; Z must not be used as beta',Local_Status='Available'))
    from scripts.run_disease_prs import STUDIES
    from scripts.fetch_part2_week7 import FILES
    for study in STUDIES.values():
        path=REPO_ROOT/'data/gwas_original'/study['file']
        if not path.exists(): continue
        d=pd.read_csv(path,sep='\t',nrows=0)
        n=sum(len(chunk) for chunk in pd.read_csv(path,sep='\t',usecols=[0],dtype=str,chunksize=250000))
        rows.append(dict(Filename=path.name,Phenotype=study['name'],Original_Source=study['source'],Variants=n,
            Columns=','.join(d.columns),Genome_Assembly='GRCh37; validated panel rsID anchors' if study.get('format')=='walters' else 'GRCh37; EBI formatted Build37 file',Effect_Allele='A1' if study.get('format')=='walters' else 'effect_allele',Other_Allele='A2' if study.get('format')=='walters' else 'other_allele',
            Effect_Size='OR; log transformed for scoring' if study['effect']=='OR' else 'log-odds beta',P_Available=True,
            N_Min=study['n'],N_Max=study['n'],
            Ancestry='European; GWAS Catalog study metadata',Direct_PRS_Eligible=True,
            Reason=study['reason'],Download_URL=FILES['data/gwas_original/'+study['file']],Download_Date='2026-11-17',Local_Status='Available; active PRS input'))
    live={row['Filename'] for row in rows}
    if not previous.empty:
        for record in previous.loc[~previous.Filename.isin(live)].to_dict('records'):
            record['Local_Status']='Removed unused Z-only input (2026-11-20)' if record['Filename'] in removed else 'Not present locally; historical inventory retained'
            rows.append(record)
    out=REPO_ROOT/'results/part2' ;out.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(rows).to_csv(out/'gwas_inventory.tsv',sep='\t',index=False)
    print(f'Inventory: {len(live)} present inputs; {len(rows)-len(live)} historical entries retained. No Z-only file is eligible for PRS.')
if __name__=='__main__': main()
