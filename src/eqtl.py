"""cis-eQTL mapping: additive SNP-expression regressions within gene body +/- 500 kb.

Extracted from the Week 4 pipeline. The Part 1 notebook and the public demo both use these functions.
"""
import gzip
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, linregress

from src.utils import run_plink


def fit_cis(variants, ids, y, sample_columns):
    """Week 3's additive regressions, one SNP per row; use only the supplied people.

    People are looked up by ID, so genotypes and expression stay paired whatever order the files use.
    """
    ids = [i for i in ids if np.isfinite(y.loc[i])]
    g = variants[[sample_columns[i] for i in ids]].to_numpy(float)
    n = np.isfinite(g).sum(axis=1)
    ordinary = variants['COUNTED'].isin(list('ACGT')) & variants['ALT'].isin(list('ACGT'))
    ordinary &= variants['COUNTED'].ne(variants['ALT'])
    keep = ordinary.to_numpy() & (n >= 3) & (np.nanstd(g, axis=1) > 0)
    result = variants.loc[keep, ['CHR', 'SNP', 'POS', 'COUNTED']].copy()
    if result.empty:
        return result.assign(BETA=[], SE=[], P=[], N=[])
    fit = linregress(g[keep], y.loc[ids].to_numpy(dtype=float)[None, :], axis=1, nan_policy='omit')
    return result.assign(BETA=fit.slope, SE=fit.stderr, P=fit.pvalue, N=n[keep])


def map_chromosome(chrom, genes, expression, ids, sample_columns, bfile, keep_file, output_dir, rebuild=False):
    """Run fit_cis for each eligible gene on one chromosome and save the results.

    `genes` needs Gene_Chr, Cis_Start, and Cis_End. Results are written to a temporary file first,
    so a half-written file is never mistaken for a finished chromosome.
    """
    prefix = Path(output_dir) / f'chr{chrom}'
    output_path = Path(f'{prefix}_cis_eqtl.tsv.gz')
    if output_path.exists() and not rebuild:
        return output_path
    run_plink(['--bfile', bfile, '--keep', keep_file, '--snps-only', 'just-acgt',
               '--min-alleles', '2', '--max-alleles', '2', '--export', 'Av', '--out', prefix])
    chromosome_genes = genes.loc[genes['Gene_Chr'].eq(str(chrom))]
    partial_path = Path(f'{prefix}_partial.tsv.gz')
    written = False
    with gzip.open(partial_path, 'wt', compresslevel=1) as output:
        for block in pd.read_csv(f'{prefix}.traw', sep=r'\s+', chunksize=2000):
            assert block['CHR'].eq(chrom).all() and set(sample_columns.values()).issubset(block.columns)
            overlapping = chromosome_genes.loc[(chromosome_genes.Cis_Start <= block.POS.max())
                                               & (chromosome_genes.Cis_End >= block.POS.min())]
            for gene, region in overlapping.iterrows():
                nearby = block.loc[block['POS'].between(region.Cis_Start, region.Cis_End)]
                if nearby.empty:
                    continue
                result = fit_cis(nearby, ids, expression.loc[gene, ids], sample_columns)
                if not result.empty:
                    result.insert(0, 'Gene_ID', gene)
                    result.to_csv(output, sep='\t', index=False, header=not written)
                    written = True
    os.replace(partial_path, output_path)
    return output_path


def add_multiple_testing_columns(associations):
    """Add BH and Bonferroni p-values, treating every row as one testing family."""
    associations = associations.copy()
    associations['P_BH'] = false_discovery_control(associations['P'], method='bh')
    associations['P_Bonferroni'] = (associations['P'] * len(associations)).clip(upper=1)
    return associations


def summarize_genes(associations):
    """One row per gene: SNPs tested, hits at each level, and the smallest p-value."""
    return associations.groupby('Gene_ID').agg(SNPs_Tested=('SNP', 'size'),
        Nominal_Hits=('P', lambda p: int((p < 0.05).sum())),
        BH_Hits=('P_BH', lambda p: int((p < 0.05).sum())),
        Bonferroni_Hits=('P_Bonferroni', lambda p: int((p < 0.05).sum())),
        Min_P_Value=('P', 'min'))
