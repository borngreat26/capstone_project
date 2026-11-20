"""Build the local target comparison from previously computed scores.

Writes only under ignored data/private_part2/. No association scan, clumping, or
genotype conversion is rerun. The source private genotype files are read only.
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.prs import read_bim, relative_position
from src.utils import REPO_ROOT


def save_figure(fig, output):
    fig.tight_layout()
    fig.savefig(output.with_suffix('.png'), dpi=320)
    fig.savefig(output.with_suffix('.pdf'))
    plt.close(fig)


def manual_target_score(bfile, weights):
    """Independently decode selected BED dosages without exporting any genotypes."""
    bim = read_bim(bfile).reset_index().set_index('SNP')
    fam = pd.read_csv(str(bfile) + '.fam', sep=r'\s+', header=None)
    if len(fam) != 1:
        raise ValueError('Expected one target individual')
    total = 0.0
    with Path(str(bfile) + '.bed').open('rb') as handle:
        if handle.read(3) != bytes([0x6c, 0x1b, 0x01]):
            raise ValueError('Expected SNP-major PLINK BED')
        for row in weights.itertuples():
            variant = bim.loc[row.SNP]
            handle.seek(3 + int(variant['index']))
            dosage = {0: 2., 1: np.nan, 2: 1., 3: 0.}[handle.read(1)[0] & 3]
            if not np.isfinite(dosage) or row.A1 not in (variant.R1, variant.R2):
                raise ValueError('Missing/incompatible scoring allele')
            total += row.BETA * (dosage if row.A1 == variant.R1 else 2 - dosage)
    return total


def build_private_report():
    root = REPO_ROOT
    private = root / 'data/private_part2'
    output = private / 'report'
    output.mkdir(parents=True, exist_ok=True)
    selected = pd.read_csv(root / 'results/part2/selected_diseases.tsv', sep='\t')
    previous = pd.read_csv(private / 'combined_summary.tsv', sep='\t').set_index('Disease')
    summaries, coverage, panels, verification = [], [], [], []

    for study in selected.itertuples():
        source = private / study.Slug
        weights = pd.read_csv(source / 'target_weights.tsv', sep='\t')
        reference = None
        reference_weights = []
        for path in sorted((source / 'reference').glob('score_chr*.sscore')):
            frame = pd.read_csv(path, sep=r'\s+', dtype={'IID': str}).set_index('IID')
            chromosome_weights = pd.read_csv(path.with_name(path.stem + '_weights.tsv'), sep='\t')
            reference_weights.append(chromosome_weights)
            if not frame.index.is_unique or not frame.ALLELE_CT.eq(2 * len(chromosome_weights)).all():
                raise ValueError('Reference sample or scoring coverage mismatch')
            score = frame.BETA_SUM
            if reference is None:
                reference = score
            else:
                if set(reference.index) != set(score.index):
                    raise ValueError('Reference sample sets differ across chromosomes')
                reference = reference + score.reindex(reference.index)
        if reference is None:
            raise FileNotFoundError('Run scripts/run_disease_prs.py to compute target/reference scores')
        joined_weights = pd.concat(reference_weights, ignore_index=True)
        if weights.SNP.duplicated().any() or joined_weights.SNP.duplicated().any():
            raise ValueError('Scoring SNPs must be unique')
        pd.testing.assert_frame_equal(
            weights.set_index('SNP')[['A1', 'BETA']].sort_index(),
            joined_weights.set_index('SNP')[['A1', 'BETA']].sort_index(),
            check_exact=False, atol=1e-12, rtol=1e-12)
        target = pd.read_csv(source / 'target.sscore', sep=r'\s+')
        if len(target) != 1 or target.ALLELE_CT.iloc[0] != 2 * len(weights):
            raise ValueError('Target scoring coverage mismatch')
        if str(target.IID.iloc[0]) in reference.index:
            raise ValueError('The target must not be included in the reference individuals')
        target_score = float(target.BETA_SUM.iloc[0])
        independent = manual_target_score(source / 'complete_calls', weights)
        np.testing.assert_allclose(independent, target_score, atol=1e-4, rtol=1e-5)
        position = relative_position(reference, target_score)
        if not np.isfinite(position['Target_Z_Score']):
            raise ValueError('Reference score variance must be positive')
        for quantity, value in position.items():
            np.testing.assert_allclose(value, previous.loc[study.Disease, quantity], atol=1e-8, rtol=1e-5)

        full = pd.read_csv(root / 'results/part2' / study.Slug / 'reference_scores.tsv',
                           sep='\t', dtype={'IID': str}).set_index('IID').PRS
        full_weights = pd.read_csv(root / 'results/part2' / study.Slug / 'weights.tsv', sep='\t')
        if set(full.index) != set(reference.index) or not set(weights.SNP) <= set(full_weights.SNP):
            raise ValueError('Full/common score comparison does not use the same reference people')
        full = full.reindex(reference.index)
        common_z = (reference - reference.mean()) / reference.std(ddof=1)
        full_z = (full - full.mean()) / full.std(ddof=1)
        full_percentile = 100 * (full.rank(method='average') - .5) / len(full)
        common_percentile = 100 * (reference.rank(method='average') - .5) / len(reference)
        shift = (full_percentile - common_percentile).abs()
        row = {'Disease': study.Disease, **position,
               'Full_Scoring_SNPs': len(full_weights), 'Shared_Scoring_SNPs': len(weights),
               'SNPs_Retained_Percent': 100 * len(weights) / len(full_weights),
               'Reference_Median': reference.median(), 'Reference_Q1': reference.quantile(.25),
               'Reference_Q3': reference.quantile(.75)}
        summaries.append(row)
        coverage.append({'Disease': study.Disease, 'Reference_Individuals': len(reference),
                         'Full_SNPs': len(full_weights), 'Shared_SNPs': len(weights),
                         'SNPs_Retained_Percent': row['SNPs_Retained_Percent'],
                         'Pearson_Full_vs_Shared': pearsonr(full, reference).statistic,
                         'Spearman_Full_vs_Shared': spearmanr(full, reference).statistic,
                         'Mean_Absolute_Percentile_Shift': shift.mean(),
                         'Maximum_Absolute_Percentile_Shift': shift.max()})
        verification.append({'Disease': study.Disease, 'Reference_Matches_Saved_Summary': True,
                             'Same_SNP_Count_For_Target_And_Reference': True,
                             'Same_SNP_Allele_Beta_Set_For_Target_And_Reference': True,
                             'Target_BED_Dosage_Check': True,
                             'Target_Absolute_Score_Error': abs(independent - target_score)})
        panels.append((study.Disease, reference, target_score, position, full_z, common_z, len(weights)))
        reference.rename('PRS').rename_axis('IID').to_csv(source / 'reference_scores.tsv', sep='\t')

    summary = pd.DataFrame(summaries)
    coverage = pd.DataFrame(coverage)
    summary.to_csv(output / 'target_summary.tsv', sep='\t', index=False)
    coverage.to_csv(output / 'coverage_sensitivity.tsv', sep='\t', index=False)
    pd.DataFrame(verification).to_csv(output / 'verification.tsv', sep='\t', index=False)

    plt.rcParams.update({'font.size': 10, 'pdf.fonttype': 42})
    fig, axes = plt.subplots(1, len(panels), figsize=(4.5 * len(panels), 4.3), squeeze=False)
    for ax, (name, reference, score, position, _, _, count) in zip(axes[0], panels):
        ax.hist(reference, bins=24, color='#287d8e', edgecolor='white', alpha=.85)
        ax.axvline(score, color='#b84334', lw=2.3,
                   label=f"Tiffany: percentile {position['Target_Percentile']:.1f}")
        ax.set(title=f'{name}\n{count} shared SNPs; {len(reference)} EUR reference people',
               xlabel='Raw PRS (log-odds weighted sum)', ylabel='Reference individuals')
        ax.text(.97, .83, f"z = {position['Target_Z_Score']:.2f}", transform=ax.transAxes,
                ha='right', color='#b84334')
        ax.legend(frameon=False, fontsize=9, loc='upper right')
        ax.set_ylim(0, ax.get_ylim()[1] * 1.18)
        ax.spines[['top', 'right']].set_visible(False)
    save_figure(fig, output / 'target_prs_distributions')

    ranked = summary.sort_values('Target_Z_Score')
    fig, ax = plt.subplots(figsize=(8, 4.2))
    ax.axvline(0, color='#666666', lw=1, label='Reference mean')
    for x in [-1, 1]:
        ax.axvline(x, color='#cccccc', ls='--', lw=.8)
    ax.scatter(ranked.Target_Z_Score, np.arange(len(ranked)), color='#b84334', s=65, zorder=3)
    for i, row in enumerate(ranked.itertuples()):
        ax.annotate(f'percentile {row.Target_Percentile:.1f}', (row.Target_Z_Score, i),
                    xytext=(8, 8), textcoords='offset points', fontsize=9)
    margin = max(2., float(ranked.Target_Z_Score.abs().max()) + .9)
    ax.set(yticks=np.arange(len(ranked)), yticklabels=ranked.Disease,
           xlim=(-margin, margin), ylim=(-.6, len(ranked) - .35),
           xlabel="Tiffany's PRS z-score within each disease's shared-SNP reference",
           title='Relative genetic-score position across diseases')
    ax.legend(frameon=False, loc='lower right')
    ax.spines[['top', 'right']].set_visible(False)
    save_figure(fig, output / 'target_standardized_comparison')

    fig, axes = plt.subplots(1, len(panels), figsize=(4.4 * len(panels), 4), squeeze=False)
    for ax, (name, _, _, _, full_z, common_z, count) in zip(axes[0], panels):
        r = coverage.set_index('Disease').loc[name, 'Pearson_Full_vs_Shared']
        ax.scatter(full_z, common_z, s=12, color='#287d8e', alpha=.55)
        ax.plot([-4, 4], [-4, 4], color='#999999', ls='--', lw=1)
        ax.set(title=f'{name}\nShared subset: {count} SNPs; Pearson r={r:.2f}',
               xlabel='Full-set reference PRS (z)', ylabel='Shared-subset reference PRS (z)')
        ax.spines[['top', 'right']].set_visible(False)
    save_figure(fig, output / 'target_coverage_sensitivity')

    captions = [
        ('target_prs_distributions', 'Part 2 results',
         'Disease PRS distributions in the EUR reference panel, with Tiffany marked. Each target/reference comparison uses the identical compatible, nonmissing SNP subset.',
         'Uncalibrated scores and reference percentiles are not disease probabilities. SNP coverage differs by disease.'),
        ('target_standardized_comparison', 'Part 2 results',
         "Tiffany's reference-standardized disease scores and empirical percentiles using matched shared-SNP scores.",
         'A higher z-score across diseases does not mean a higher absolute probability of that disease.'),
        ('target_coverage_sensitivity', 'Part 2 supplementary',
         'Full-set versus target-compatible subset scores in the same EUR reference people, standardized separately. Dashed lines indicate equal standardized scores.',
         "These points are reference people. Tiffany's full-set score is unknown; this analysis describes sensitivity to SNP coverage, not clinical accuracy."),
    ]
    inventory = []
    for stem, section, caption, limitations in captions:
        for extension in ['png', 'pdf']:
            inventory.append({'Filename': str((output / f'{stem}.{extension}').relative_to(root)),
                              'Purpose': stem.replace('_', ' '),
                              'Data_Source': 'Saved target-compatible PRSs and full public reference PRSs',
                              'Paper_Section': section, 'Caption': caption, 'Limitations': limitations})
    pd.DataFrame(inventory).to_csv(output / 'figure_inventory.tsv', sep='\t', index=False)
    notes = ['# Tiffany: local Part 2 results', '',
             'Generated from saved scores; original genotype files are read only.',
             'Percentiles use half credit for ties; z-scores use reference sample SD.', '',
             '## Target comparisons', '']
    for row in summary.itertuples():
        notes.append(f'- {row.Disease}: raw score {row.Target_PRS:.6g}; reference percentile '
                     f'{row.Target_Percentile:.2f}; z={row.Target_Z_Score:.3f}; '
                     f'{row.Shared_Scoring_SNPs}/{row.Full_Scoring_SNPs} selected SNPs retained '
                     f'({row.SNPs_Retained_Percent:.1f}%).')
    notes.extend(['', '## Coverage sensitivity in the reference individuals', ''])
    for row in coverage.itertuples():
        notes.append(f'- {row.Disease}: full/shared Pearson r={row.Pearson_Full_vs_Shared:.3f}, '
                     f'Spearman r={row.Spearman_Full_vs_Shared:.3f}; mean absolute percentile '
                     f'change {row.Mean_Absolute_Percentile_Shift:.2f} points.')
    notes.extend(['', 'These are relative genetic-score positions, not disease probabilities or '
                  'validated clinical predictions. No disease outcomes were used to choose thresholds. '
                  'The full-set target score cannot be reconstructed from missing genotypes. '
                  'Reference sensitivity does not measure clinical prediction accuracy. '
                  'Ancestry, study ascertainment, limited targeted-array GWAS coverage, and LD '
                  'also limit comparisons. No private ancestry inference was performed.', ''])
    (output / 'methods_results.md').write_text('\n'.join(notes))
    return summary, coverage


if __name__ == '__main__':
    summary, _ = build_private_report()
    print(f'Built local target comparisons for {len(summary)} diseases; no genotype records exported.')
