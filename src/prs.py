"""Polygenic score workflow (P+T): QC, harmonization, clumping, scoring, and relative position.

Extracted from the verified Week 5 height analysis, for use in the demo and in Part 2.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.utils import run_plink

COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
AMBIGUOUS_PAIRS = ['AT', 'TA', 'CG', 'GC']


def gwas_qc(gwas, min_info=0.8, min_maf=0.01):
    """Keep unique SNP IDs with valid p-values and unambiguous DNA alleles; apply INFO/MAF filters if present."""
    unique = gwas.drop_duplicates('SNP', keep='first')
    passed = (unique['P'].between(0, 1) & unique['A1'].isin(list('ACGT')) & unique['A2'].isin(list('ACGT'))
              & unique['A1'].ne(unique['A2']) & ~(unique['A1'] + unique['A2']).isin(AMBIGUOUS_PAIRS))
    if 'INFO' in unique:
        passed &= unique['INFO'].gt(min_info)
    if 'MAF' in unique:
        passed &= unique['MAF'].gt(min_maf)
    return unique.loc[passed].copy()


def read_bim(bfile):
    """Read a PLINK .bim file; R1 and R2 are the two alleles in the reference panel."""
    return pd.read_csv(f'{bfile}.bim', sep=r'\s+', header=None, names=['CHR', 'SNP', 'CM', 'BP', 'R1', 'R2'])


def harmonize_to_reference(gwas, reference_variants):
    """Match GWAS SNPs to the panel by ID, position, and alleles.

    Needs CHR, BP, SNP, A1, A2, BETA, and P. Nonambiguous opposite-strand SNPs are relabeled with
    their complements. The beta keeps its sign, because A1 still names the same biological allele.
    Returns the harmonized table and counts of what happened to each input SNP.
    """
    joined = gwas.merge(reference_variants, on='SNP', suffixes=('_GWAS', ''))
    direct = ((joined['A1'].eq(joined['R1']) & joined['A2'].eq(joined['R2']))
              | (joined['A1'].eq(joined['R2']) & joined['A2'].eq(joined['R1'])))
    c1, c2 = joined['A1'].map(COMPLEMENT), joined['A2'].map(COMPLEMENT)
    reverse = ((c1.eq(joined['R1']) & c2.eq(joined['R2'])) | (c1.eq(joined['R2']) & c2.eq(joined['R1'])))
    same_position = joined['CHR_GWAS'].eq(joined['CHR']) & joined['BP_GWAS'].eq(joined['BP'])
    harmonized = joined.loc[(direct | reverse) & same_position].copy()
    harmonized['A1'] = np.where(direct.loc[harmonized.index], harmonized['A1'], c1.loc[harmonized.index])
    harmonized['A2'] = np.where(direct.loc[harmonized.index], harmonized['A2'], c2.loc[harmonized.index])
    counts = {'Input SNPs after QC': len(gwas), 'Not in reference panel': len(gwas) - len(joined),
              'Different position': int((~same_position).sum()),
              'Incompatible alleles': int((same_position & ~(direct | reverse)).sum()),
              'Complemented to match the panel strand': int((~direct & reverse & same_position).sum()),
              'Harmonized SNPs': len(harmonized)}
    return harmonized[['CHR', 'BP', 'SNP', 'A1', 'A2', 'BETA', 'P', 'R1', 'R2']], counts


def clump(harmonized, bfile, work_prefix, r2=0.1, kb=250):
    """Return the index SNPs from PLINK2 LD clumping of one chromosome's harmonized GWAS."""
    work_prefix = Path(work_prefix)
    gwas_path = Path(f'{work_prefix}_gwas.tsv')
    harmonized.to_csv(gwas_path, sep='\t', index=False)
    run_plink(['--bfile', bfile, '--clump', gwas_path, '--clump-id-field', 'SNP', '--clump-p-field', 'P',
               '--clump-a1-field', 'A1', '--clump-p1', '1', '--clump-p2', '1', '--clump-r2', r2,
               '--clump-kb', kb, '--out', work_prefix])
    return pd.read_csv(f'{work_prefix}.clumps', sep=r'\s+')['ID'].tolist()


def select_weights(harmonized, index_snps, p_threshold):
    """Keep clump index SNPs at or below a fixed p-value threshold."""
    selected = harmonized.loc[harmonized['SNP'].isin(index_snps) & harmonized['P'].le(p_threshold)]
    return selected.sort_values(['CHR', 'BP']).reset_index(drop=True)


def score_bfile(bfile, weights, work_prefix):
    """Sum beta x effect-allele dosage per person. Missing calls are an error, never mean-imputed."""
    weights_path = Path(f'{work_prefix}_weights.tsv')
    weights[['SNP', 'A1', 'BETA']].to_csv(weights_path, sep='\t', index=False)
    run_plink(['--bfile', bfile, '--score', weights_path, '1', '2', '3', 'header-read', 'no-mean-imputation',
               'cols=+scoresums', '--out', work_prefix])
    scored = pd.read_csv(f'{work_prefix}.sscore', sep=r'\s+', dtype={'IID': str}).set_index('IID')
    assert scored['ALLELE_CT'].eq(2 * len(weights)).all(), 'Some scoring SNPs were missing or had missing calls.'
    return scored['BETA_SUM'].rename('PRS')


def relative_position(reference_scores, target_score):
    """Target's percentile (half-credit for ties) and z-score (sample SD) within the reference group."""
    mean, sd = reference_scores.mean(), reference_scores.std(ddof=1)
    below, ties = reference_scores.lt(target_score).sum(), reference_scores.eq(target_score).sum()
    return {'Reference_Individuals': len(reference_scores), 'Reference_Mean': mean, 'Reference_SD': sd,
            'Target_PRS': target_score, 'Target_Percentile': 100 * (below + 0.5 * ties) / len(reference_scores),
            'Target_Z_Score': (target_score - mean) / sd}


def plot_prs_distribution(reference_scores, target_score, title, target_label, output_path):
    """Histogram of reference scores with a vertical line at the target's score."""
    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.hist(reference_scores, bins=24, color='C0', edgecolor='white', linewidth=0.5)
    ax.axvline(target_score, color='C3', linewidth=2, label=target_label)
    ax.set(title=title, xlabel='Polygenic score (sum of beta x effect-allele dosage)', ylabel='Reference individuals')
    ax.legend(frameon=False)
    ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
