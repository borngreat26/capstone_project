"""Run the public demo: a miniature cis-eQTL scan and a P+T polygenic score on a small chr22 region.

Usage, from the repository root:  python scripts/run_demo.py

This reproduces the computational workflow only. The genotypes and expression values are real public data,
but the GWAS is simulated. The results are not Tiffany's results or a disease-risk estimate.
"""
import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from src.eqtl import add_multiple_testing_columns, map_chromosome, summarize_genes
from src.prs import (clump, gwas_qc, harmonize_to_reference, plot_prs_distribution, read_bim, relative_position,
                     score_bfile, select_weights)
from src.utils import REPO_ROOT, plink2_path

DEMO_DATA = REPO_ROOT / 'data/demo'
P_THRESHOLD = 0.001  # Fixed in advance, the same cutoff as Week 5. It is never tuned on the demo individual.


def run_eqtl_demo(output, work):
    """Same mapping functions as Part 1, on a few genes and the people with both data types."""
    genes = pd.read_csv(DEMO_DATA / 'demo_genes.tsv', sep='\t', dtype={'Gene_Chr': str}).set_index('Gene_ID')
    expression = pd.read_csv(DEMO_DATA / 'demo_expression.tsv', sep='\t').set_index('Gene_ID')
    fam = pd.read_csv(DEMO_DATA / 'demo_eqtl.fam', sep=r'\s+', header=None, dtype=str)
    ids = fam[1].tolist()
    assert set(ids) == set(expression.columns), 'Genotype and expression people must match.'
    fam[[0, 1]].to_csv(work / 'eqtl.keep', sep='\t', index=False, header=False)
    sample_columns = dict(zip(ids, fam[0] + '_' + fam[1]))
    path = map_chromosome(22, genes, expression, ids, sample_columns, DEMO_DATA / 'demo_eqtl', work / 'eqtl.keep',
                          work, rebuild=True)
    associations = add_multiple_testing_columns(pd.read_csv(path, sep='\t'))
    summary = summarize_genes(associations).join(genes['Gene_Symbol']).sort_values('Min_P_Value')
    associations.to_csv(output / 'demo_eqtl_associations.tsv', sep='\t', index=False)
    summary.to_csv(output / 'demo_eqtl_gene_summary.tsv', sep='\t')

    top_gene = summary.index[0]
    locus = associations.loc[associations['Gene_ID'].eq(top_gene)]
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.scatter(locus['POS'] / 1e6, -np.log10(locus['P']), s=14, color='C0')
    ax.axhline(-np.log10(0.05 / len(associations)), linestyle='--', color='gray', label='Bonferroni (demo tests only)')
    ax.set(title=f"Demo cis-eQTL scan: {summary.loc[top_gene, 'Gene_Symbol']} ({len(locus)} SNPs tested)",
           xlabel='Chromosome 22 position (Mb, GRCh37)', ylabel='-log10(p-value)')
    ax.legend(frameon=False)
    ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    fig.savefig(output / 'demo_eqtl_top_gene.png', dpi=180)
    plt.close(fig)
    print(f'cis-eQTL demo: {len(summary)} genes, {len(associations):,} SNP-gene tests, '
          f"strongest gene {summary.loc[top_gene, 'Gene_Symbol']} (p = {summary['Min_P_Value'].iloc[0]:.2e})")


def run_prs_demo(output, work):
    """Same P+T functions planned for Part 2: QC, harmonize, clump, threshold, score, and compare."""
    reference_bfile, individual_bfile = DEMO_DATA / 'demo_reference', DEMO_DATA / 'demo_individual'
    gwas = pd.read_csv(DEMO_DATA / 'demo_gwas.tsv', sep='\t')
    passed_qc = gwas_qc(gwas)
    harmonized, counts = harmonize_to_reference(passed_qc, read_bim(reference_bfile))
    index_snps = clump(harmonized, reference_bfile, work / 'clump')
    weights = select_weights(harmonized, index_snps, P_THRESHOLD)
    assert len(weights), 'No SNPs passed the fixed p-value threshold.'

    reference_scores = score_bfile(reference_bfile, weights, work / 'reference')
    individual_scores = score_bfile(individual_bfile, weights, work / 'individual')
    assert len(individual_scores) == 1 and not set(individual_scores.index) & set(reference_scores.index), \
        'The demo individual must not be in the reference group.'
    position = relative_position(reference_scores, float(individual_scores.iloc[0]))

    phenotype = pd.read_csv(DEMO_DATA / 'demo_reference_simulated_phenotype.tsv', sep='\t', dtype={'IID': str}).set_index('IID')
    correlation = np.corrcoef(reference_scores, phenotype['Simulated_Phenotype'].loc[reference_scores.index])[0, 1]

    stages = pd.Series({'GWAS rows': len(gwas), **counts, 'Clump index SNPs': len(index_snps),
                        f'Index SNPs with p <= {P_THRESHOLD}': len(weights)}, name='SNPs')
    stages.rename_axis('Stage').to_csv(output / 'demo_snp_counts.tsv', sep='\t')
    harmonized.to_csv(output / 'demo_gwas_harmonized.tsv', sep='\t', index=False)
    weights[['SNP', 'A1', 'BETA', 'P']].to_csv(output / 'demo_weights.tsv', sep='\t', index=False)
    reference_scores.rename_axis('IID').to_csv(output / 'demo_reference_prs.tsv', sep='\t')
    summary = pd.Series({**position, 'P_Threshold': P_THRESHOLD, 'Scoring_SNPs': len(weights),
                         'Reference_Correlation_With_Simulated_Phenotype': correlation}, name='Value')
    summary.rename_axis('Quantity').to_csv(output / 'demo_prs_summary.tsv', sep='\t')
    plot_prs_distribution(reference_scores, position['Target_PRS'],
        f"Demo score on a simulated GWAS: {len(weights)} SNPs, chr22 region, {len(reference_scores)} EUR reference people",
        f"demo_individual (percentile {position['Target_Percentile']:.1f}, z = {position['Target_Z_Score']:.2f})",
        output / 'demo_prs_distribution.png')
    print(stages.to_string())
    print(f"PRS demo: {len(weights)} SNPs | demo_individual percentile {position['Target_Percentile']:.1f}, "
          f"z = {position['Target_Z_Score']:.2f} | score vs simulated phenotype r = {correlation:.2f}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=REPO_ROOT / 'results/demo', help='Folder for the demo outputs')
    args = parser.parse_args()
    work = args.output / 'work'
    work.mkdir(parents=True, exist_ok=True)
    print('Using PLINK2 at:', plink2_path())
    run_eqtl_demo(args.output, work)
    run_prs_demo(args.output, work)
    print(f'Finished. Outputs are in {args.output}')
    print('Workflow demonstration only: simulated GWAS, small region, not a disease-risk estimate.')


if __name__ == '__main__':
    main()
