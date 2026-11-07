"""Inventory the local GWAS summary statistics, to see which can support a polygenic score.

Usage:  python scripts/inspect_sumstats.py
Reads data/sumstats/ and data/LDREF/, and writes docs/gwas_sumstats_inventory.tsv. No scores are computed.
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from src.prs import AMBIGUOUS_PAIRS, COMPLEMENT
from src.utils import REPO_ROOT

panel = pd.concat([pd.read_csv(f'{REPO_ROOT}/data/LDREF/1000G.EUR.{c}.bim', sep=r'\s+', header=None, usecols=[1, 4, 5],
                               names=['SNP', 'R1', 'R2']) for c in range(1, 23)])
heritability = pd.read_csv(REPO_ROOT / 'data/sumstats/alltraits.txt', sep='\t').set_index('Trait_Identifier')
rows = []
for path in sorted((REPO_ROOT / 'data/sumstats').glob('*.sumstats.gz')):
    trait = path.name.removesuffix('.sumstats.gz')
    table = pd.read_csv(path, sep='\t', dtype={'A1': str, 'A2': str})
    columns = table.columns.tolist()
    joined = table.merge(panel, on='SNP')
    direct = ((joined['A1'] == joined['R1']) & (joined['A2'] == joined['R2'])) | ((joined['A1'] == joined['R2']) & (joined['A2'] == joined['R1']))
    c1, c2 = joined['A1'].map(COMPLEMENT), joined['A2'].map(COMPLEMENT)
    flipped = ((c1 == joined['R1']) & (c2 == joined['R2'])) | ((c1 == joined['R2']) & (c2 == joined['R1']))
    rows.append({
        'Trait_File': trait, 'Columns': ','.join(columns), 'A2_Listed_Before_A1': columns.index('A2') < columns.index('A1'),
        'Rows': len(table), 'Median_N': table['N'].median(), 'Missing_Z_Share': round(table['Z'].isna().mean(), 4),
        'Has_Beta_Or_OR': bool({'BETA', 'OR'} & set(columns)), 'Has_SE': 'SE' in columns, 'Has_P': 'P' in columns,
        'Has_Position': bool({'CHR', 'BP', 'POS'} & set(columns)), 'SNPs_In_Panel': len(joined),
        'Panel_Coverage_Pct': round(100 * len(joined) / len(panel), 1),
        'Allele_Match_Pct': round(100 * (direct | flipped).mean(), 2) if len(joined) else None,
        'Ambiguous_Share_In_Panel': round((joined['A1'] + joined['A2']).isin(AMBIGUOUS_PAIRS).mean(), 4),
        'h2g_From_alltraits': heritability['h2g'].get(trait), 'Listed_In_alltraits': trait in heritability.index,
    })
inventory = pd.DataFrame(rows)
inventory.to_csv(REPO_ROOT / 'docs/gwas_sumstats_inventory.tsv', sep='\t', index=False)
print(len(inventory), 'files | with beta or OR:', int(inventory['Has_Beta_Or_OR'].sum()), '| with p-values:',
      int(inventory['Has_P'].sum()), '| with positions:', int(inventory['Has_Position'].sum()))
