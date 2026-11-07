"""Build the small public demo dataset in data/demo/. Maintainers only: graders do not need to run this.

Needs the full local inputs (the LDREF panel, the expression file, and the Part 1 gene table). It reads
only public data. The GWAS is simulated, and every random choice uses SEED.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import linregress

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from src.prs import AMBIGUOUS_PAIRS, COMPLEMENT
from src.utils import REPO_ROOT, run_plink

SEED = 42
CHROM, REGION_START, REGION_END = 22, 44_500_000, 47_000_000
CAUSAL_SNPS, SIMULATED_HERITABILITY = 15, 0.5

panel = REPO_ROOT / 'data/LDREF' / f'1000G.EUR.{CHROM}'
expression_path = REPO_ROOT / 'data/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt'
gene_table_path = REPO_ROOT / 'data/quarter1_results/part1/protein_coding_genes.tsv'
demo = REPO_ROOT / 'data/demo'
work = demo / 'work_tmp'
work.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(SEED)

# Split the panel: one demo individual, a simulated GWAS group, and a separate reference group.
fam = pd.read_csv(f'{panel}.fam', sep=r'\s+', header=None, dtype=str)
shuffled = rng.permutation(fam[1].to_numpy())
demo_id, discovery_ids, reference_ids = shuffled[0], shuffled[1:245], shuffled[245:]
expression_header = pd.read_csv(expression_path, sep='\t', nrows=0).columns[4:]
matched_ids = [i for i in fam[1] if i in set(expression_header)]


def extract(ids, name, folder=demo):
    keep = work / f'{name}.keep'
    pd.DataFrame({'FID': ids, 'IID': ids}).to_csv(keep, sep='\t', index=False, header=False)
    run_plink(['--bfile', panel, '--keep', keep, '--chr', CHROM, '--from-bp', REGION_START, '--to-bp', REGION_END,
               '--snps-only', 'just-acgt', '--max-alleles', '2', '--make-bed', '--out', folder / name])


extract([demo_id], 'demo_individual')
extract(reference_ids, 'demo_reference')
extract(matched_ids, 'demo_eqtl')
extract(discovery_ids, 'discovery', work)
run_plink(['--bfile', work / 'discovery', '--export', 'Av', '--out', work / 'discovery'])
run_plink(['--bfile', demo / 'demo_reference', '--export', 'Av', '--out', work / 'reference'])
discovery = pd.read_csv(work / 'discovery.traw', sep=r'\s+')
reference = pd.read_csv(work / 'reference.traw', sep=r'\s+')
discovery_g = discovery.iloc[:, 6:].to_numpy(float)
reference_g = reference.iloc[:, 6:].to_numpy(float)

# Simulate a phenotype from real genotypes, then run a real per-SNP GWAS in the discovery group.
frequency = discovery_g.mean(axis=1) / 2
common = np.nonzero((frequency > 0.05) & (frequency < 0.95))[0]
causal = rng.choice(common, CAUSAL_SNPS, replace=False)
effects = rng.normal(size=CAUSAL_SNPS)
genetic_discovery = discovery_g[causal].T @ effects
scale = np.sqrt(SIMULATED_HERITABILITY / genetic_discovery.var())
noise_sd = np.sqrt(1 - SIMULATED_HERITABILITY)
phenotype = scale * genetic_discovery + rng.normal(scale=noise_sd, size=discovery_g.shape[1])
varies = discovery_g.std(axis=1) > 0
fit = linregress(discovery_g[varies], phenotype[None, :], axis=1)
gwas = pd.DataFrame({'CHR': discovery.loc[varies, 'CHR'], 'BP': discovery.loc[varies, 'POS'],
                     'SNP': discovery.loc[varies, 'SNP'], 'A1': discovery.loc[varies, 'COUNTED'],
                     'A2': discovery.loc[varies, 'ALT'], 'N': discovery_g.shape[1], 'BETA': fit.slope,
                     'SE': fit.stderr, 'P': fit.pvalue,
                     'MAF': np.minimum(frequency[varies], 1 - frequency[varies])}).reset_index(drop=True)

# Report half of the SNPs with the other allele as A1; the beta changes sign, so the meaning is the same.
swap = rng.random(len(gwas)) < 0.5
gwas.loc[swap, ['A1', 'A2']] = gwas.loc[swap, ['A2', 'A1']].to_numpy()
gwas.loc[swap, 'BETA'] *= -1

# Add a few rows that harmonization must fix or remove. Only SNPs with weak signals are altered,
# and none of them is ambiguous.
gwas['Demo_Case'] = 'As simulated'
safe = gwas.index[~(gwas['A1'] + gwas['A2']).isin(AMBIGUOUS_PAIRS) & gwas['P'].gt(0.05)].to_numpy()
changed = rng.choice(safe, 11, replace=False)
strand, moved, wrong_alleles = changed[:6], changed[6:8], changed[8:11]
gwas.loc[strand, 'A1'] = gwas.loc[strand, 'A1'].map(COMPLEMENT)
gwas.loc[strand, 'A2'] = gwas.loc[strand, 'A2'].map(COMPLEMENT)
gwas.loc[strand, 'Demo_Case'] = 'Opposite strand'
gwas.loc[moved, 'BP'] += 1000
gwas.loc[moved, 'Demo_Case'] = 'Wrong position'
# Complementing only A2 gives a nonambiguous pair (A/G becomes A/C) that fits neither strand of the panel.
gwas.loc[wrong_alleles, 'A2'] = gwas.loc[wrong_alleles, 'A2'].map(COMPLEMENT)
gwas.loc[wrong_alleles, 'Demo_Case'] = 'Incompatible alleles'
not_in_panel = gwas.loc[rng.choice(safe, 3, replace=False)].assign(
    SNP=[f'demo_not_in_panel_{i}' for i in range(1, 4)], Demo_Case='Not in panel')
duplicates = gwas.loc[rng.choice(safe, 2, replace=False)].assign(Demo_Case='Duplicate row')
gwas = pd.concat([gwas, not_in_panel, duplicates], ignore_index=True).sort_values(['BP', 'SNP'])
gwas.to_csv(demo / 'demo_gwas.tsv', sep='\t', index=False)

# Simulated phenotype in the separate reference group, used only to confirm the score points the right way.
genetic_reference = scale * reference_g[causal].T @ effects
pd.DataFrame({'IID': reference.columns[6:].str.split('_', n=1).str[1],
              'Simulated_Phenotype': genetic_reference + rng.normal(scale=noise_sd, size=reference_g.shape[1])
              }).to_csv(demo / 'demo_reference_simulated_phenotype.tsv', sep='\t', index=False)

# Miniature eQTL data: chromosome 22 protein-coding genes whose whole cis-window is inside the region.
genes = pd.read_csv(gene_table_path, sep='\t')
genes = genes.loc[genes['Gene_Chr'].astype(str).eq(str(CHROM)) & genes['Reason'].ne('No expression measurement')].copy()
genes['Cis_Start'] = (genes['Gene_Start'] - 500000).clip(lower=1)
genes['Cis_End'] = genes['Gene_End'] + 500000
genes = genes.loc[genes['Cis_Start'].ge(REGION_START) & genes['Cis_End'].le(REGION_END)]
genes[['Gene_ID', 'Gene_Symbol', 'Gene_Chr', 'Gene_Start', 'Gene_End', 'Cis_Start', 'Cis_End']].to_csv(
    demo / 'demo_genes.tsv', sep='\t', index=False)
expression = pd.read_csv(expression_path, sep='\t', usecols=['TargetID', *matched_ids]).rename(columns={'TargetID': 'Gene_ID'})
expression = expression.loc[expression['Gene_ID'].isin(genes['Gene_ID']), ['Gene_ID', *matched_ids]]
expression.round(6).to_csv(demo / 'demo_expression.tsv', sep='\t', index=False)

for path in work.iterdir():
    path.unlink()
work.rmdir()
for path in demo.glob('*.log'):
    path.unlink()
print('Demo individual:', demo_id, '| reference:', len(reference_ids), '| simulated GWAS people:', len(discovery_ids))
print('Region SNPs:', len(discovery), '| GWAS rows:', len(gwas), '| eQTL genes:', len(genes))
print(gwas['Demo_Case'].value_counts().to_string())
print('Total size (KB):', sum(p.stat().st_size for p in demo.iterdir()) // 1024)
