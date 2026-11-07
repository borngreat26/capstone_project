"""The public notebook's Week 5: the height polygenic score workflow, run on a public 1000 Genomes participant.

make_public_notebook.py puts these cells in place of the private Week 5 section. They use only public data.
"""

CELLS = [
('markdown', """# Week 5 — Height Polygenic Score for a Public Participant

The original Week 5 analysis used a private 23andMe file, which is not published. This public version runs the same steps on **NA20792**, one person chosen at random (seed 42) from the public 1000 Genomes EUR panel, the same `demo_individual` as the demo. **These results are this participant's and do not reproduce the private analysis.**

The [height GWAS](https://choishingwan.github.io/PRS-Tutorial/base/) is **simulated**, so this is an educational comparison, not an estimate of anyone's height.

### Task 1: The Participant's Genotypes in PLINK and VCF

Split the panel into the participant and a 488-person reference group, so the participant is not compared with themselves. There is no 23andMe file here, so the conversion step becomes a VCF export from PLINK."""),
('code', """from pathlib import Path
import gzip
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm, gaussian_kde
from src.prs import clump, gwas_qc, harmonize_to_reference, read_bim, relative_position, score_bfile, select_weights
from src.utils import run_plink

w5 = Path('data/week5_public_results')
work5 = w5 / 'work'
work5.mkdir(parents=True, exist_ok=True)
chromosomes5 = range(1, 23)
participant_id = pd.read_csv('data/demo/demo_individual.fam', sep=r'\\s+', header=None, dtype=str).iloc[0, 1]
panel_fam5 = pd.read_csv('data/LDREF/1000G.EUR.1.fam', sep=r'\\s+', header=None, dtype=str)
assert panel_fam5[1].eq(participant_id).sum() == 1
keep5 = work5 / 'participant.keep'
panel_fam5.loc[panel_fam5[1].eq(participant_id), [0, 1]].to_csv(keep5, sep='\\t', header=False, index=False)

for chrom in chromosomes5:
    panel = f'data/LDREF/1000G.EUR.{chrom}'
    if not (work5 / f'reference_chr{chrom}.bed').exists():
        run_plink(['--bfile', panel, '--keep', keep5, '--make-bed', '--out', work5 / f'participant_chr{chrom}'])
        run_plink(['--bfile', panel, '--keep', keep5, '--export', 'vcf', 'id-paste=iid', '--out', work5 / f'participant_chr{chrom}'])
        run_plink(['--bfile', panel, '--remove', keep5, '--make-bed', '--out', work5 / f'reference_chr{chrom}'])
print('Participant:', participant_id, '| Reference people:', len(panel_fam5) - 1)"""),
('markdown', """### Exporting and Checking the VCF

Combine the chromosomes into one VCF, then confirm that every call matches the PLINK `.bed` file."""),
('code', """columns5 = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'GT']
alt_copies5 = {'0/0': 0, '0/1': 1, '1/0': 1, '1/1': 2}
bed_to_alt5 = np.array([2, -1, 1, 0])  # .bed codes: 0 = two copies of .bim column 5 (VCF ALT), 2 = one, 3 = none, 1 = missing
genotype_counts5 = pd.Series(0, index=[0, 1, 2])
with gzip.open(w5 / 'participant_autosomal.vcf.gz', 'wt') as merged:
    for chrom in chromosomes5:
        path = work5 / f'participant_chr{chrom}.vcf'
        if chrom == 1:
            merged.writelines(line for line in path.open() if line.startswith('#'))
        calls = pd.read_csv(path, sep='\\t', comment='#', header=None, names=columns5, dtype={'GT': str})
        calls.to_csv(merged, sep='\\t', header=False, index=False)
        vcf_alt = calls['GT'].map(alt_copies5).to_numpy()
        bed_alt = bed_to_alt5[np.fromfile(work5 / f'participant_chr{chrom}.bed', dtype=np.uint8, offset=3) & 3]
        bim = read_bim(work5 / f'participant_chr{chrom}')
        assert (vcf_alt == bed_alt).all() and (bed_alt >= 0).all()
        assert calls['ID'].eq(bim['SNP']).all() and calls['ALT'].eq(bim['R1']).all() and calls['REF'].eq(bim['R2']).all()
        genotype_counts5 = genotype_counts5.add(pd.Series(vcf_alt).value_counts(), fill_value=0).astype(int)

genotype_labels5 = ['No copies of the ALT allele', 'One copy', 'Two copies']
fig, ax = plt.subplots(figsize=(8, 2.6))
bars = ax.barh(genotype_labels5, genotype_counts5, color='C0')
ax.bar_label(bars, labels=[f'{n:,} ({n / genotype_counts5.sum():.1%})' for n in genotype_counts5], padding=4)
ax.invert_yaxis()
ax.set_xlim(0, genotype_counts5.max() * 1.3)
ax.set_title(f'{participant_id}: genotypes at {genotype_counts5.sum():,} SNPs (VCF and PLINK agree)')
ax.xaxis.set_visible(False)
ax.spines[['top', 'right', 'bottom']].set_visible(False)
fig.tight_layout()
fig.savefig(w5 / 'participant_genotype_summary.png', dpi=180)
plt.show()"""),
('markdown', """### Task 2: Height GWAS Summary Statistics

**A1 is the effect allele**, and the [tutorial](https://choishingwan.github.io/PRS-Tutorial/plink/) uses **BETA = ln(OR)**. Base QC follows the tutorial: **INFO > 0.8, MAF > 0.01**, valid alleles, unique SNP IDs, and no A/T or C/G SNPs."""),
('code', """height_gwas = pd.read_csv('data/Height.gwas.txt.gz', sep='\\t')
assert height_gwas['OR'].gt(0).all()
assert np.allclose(2 * norm.sf(abs(np.log(height_gwas['OR']) / height_gwas['SE'])), height_gwas['P'], rtol=0, atol=1e-5)
height_gwas['BETA'] = np.log(height_gwas['OR'])
height_passed = gwas_qc(height_gwas)
display(height_gwas.head(5))

manhattan5 = height_gwas.drop_duplicates('SNP')[['CHR', 'BP', 'P']].sort_values(['CHR', 'BP'])
offsets5 = manhattan5.groupby('CHR')['BP'].max().cumsum().shift(fill_value=0)
x5 = (manhattan5['BP'] + manhattan5['CHR'].map(offsets5)) / 1e6
fig, ax = plt.subplots(figsize=(13, 4))
ax.scatter(x5, -np.log10(manhattan5['P']), s=3, c=np.where(manhattan5['CHR'] % 2, 'C0', 'gray'))
ax.axhline(-np.log10(0.001), linestyle='--', color='black', linewidth=1, label='Score cutoff, p = 0.001')
middles5 = x5.groupby(manhattan5['CHR']).median()
ax.set_xticks(middles5, middles5.index, fontsize=8)
ax.set(xlabel='Chromosome', ylabel='−log10(p)', title=f'Assigned height GWAS (simulated): {len(manhattan5):,} unique SNPs')
ax.margins(x=0.01)
ax.spines[['top', 'right']].set_visible(False)
ax.legend(frameon=False, loc='upper right')
fig.tight_layout()
fig.savefig(w5 / 'height_gwas_manhattan.png', dpi=180)
plt.show()"""),
('markdown', """### Task 3: Clumping and Thresholding

Match GWAS SNPs to the panel by ID, position, and alleles, then clump within 250 kb at **r² > 0.1** using the reference group's genotypes. Index SNPs with **p ≤ 0.001** are kept. This cutoff is fixed and never tuned."""),
('code', """reference_variants5 = pd.concat([read_bim(work5 / f'reference_chr{c}') for c in chromosomes5], ignore_index=True)
height_harmonized, harmonize_counts5 = harmonize_to_reference(height_passed, reference_variants5)
index_snps5 = []
for chrom in chromosomes5:
    index_snps5 += clump(height_harmonized.loc[height_harmonized['CHR'].eq(chrom)], work5 / f'reference_chr{chrom}', work5 / f'clump_chr{chrom}')
weights5 = select_weights(height_harmonized, index_snps5, 0.001)
assert weights5['SNP'].is_unique and weights5['P'].le(0.001).all()
weights5.to_csv(w5 / 'height_weights_p001.tsv', sep='\\t', index=False)

stages5 = pd.Series({'GWAS rows': len(height_gwas), 'Passed QC': len(height_passed), 'Harmonized with the panel': len(height_harmonized),
                     'Clump index SNPs': len(index_snps5), 'Index SNPs with p ≤ 0.001': len(weights5)})
display(stages5.to_frame('SNPs'))
by_chromosome5 = weights5.groupby('CHR').size().reindex(chromosomes5, fill_value=0)
fig, (stage_ax, chromosome_ax) = plt.subplots(1, 2, figsize=(14, 3.8))
bars = stage_ax.barh(stages5.index, stages5, color='C0')
stage_ax.bar_label(bars, fmt='{:,.0f}', padding=4)
stage_ax.set_xscale('log')
stage_ax.set_xlim(1, stages5.max() * 8)
stage_ax.invert_yaxis()
stage_ax.set_title('SNPs remaining at each stage (log scale)')
stage_ax.xaxis.set_visible(False)
stage_ax.spines['bottom'].set_visible(False)
chromosome_ax.bar(by_chromosome5.index, by_chromosome5, color='C0')
chromosome_ax.set_xticks(by_chromosome5.index, by_chromosome5.index, fontsize=9)
chromosome_ax.set(xlabel='Chromosome', ylabel='Scoring SNPs', title='Scoring SNPs by chromosome')
for ax in (stage_ax, chromosome_ax):
    ax.spines[['top', 'right']].set_visible(False)
fig.tight_layout()
fig.savefig(w5 / 'snp_filtering.png', dpi=180)
plt.show()"""),
('markdown', """### Task 4: Computing the Height PRS

The score is $PRS_i=\\sum_j\\beta_jG_{ij}$, summed across chromosomes. Everyone gets the same SNPs, weights, and effect alleles, and missing calls are an error rather than being filled in.

The z-score uses the reference mean and sample SD. The percentile is a rank in this **488-person EUR reference**. It is not a disease probability, and neither value converts to height."""),
('code', """reference_prs5 = pd.Series(0.0, index=panel_fam5.loc[panel_fam5[1].ne(participant_id), 1], name='PRS')
participant_prs5 = 0.0
for chrom in chromosomes5:
    weights = weights5.loc[weights5['CHR'].eq(chrom)]
    if len(weights):
        reference_prs5 = reference_prs5.add(score_bfile(work5 / f'reference_chr{chrom}', weights, work5 / f'score_reference_chr{chrom}'))
        participant_prs5 += float(score_bfile(work5 / f'participant_chr{chrom}', weights, work5 / f'score_participant_chr{chrom}').iloc[0])
assert reference_prs5.notna().all()
position5 = relative_position(reference_prs5, participant_prs5)
summary5 = pd.Series({'Scoring SNPs': len(weights5), **position5})
summary5.to_csv(w5 / 'participant_prs_summary.tsv', sep='\\t', header=['Value'])
reference_prs5.rename_axis('IID').to_csv(w5 / 'reference_height_prs.tsv', sep='\\t')
display(summary5.to_frame('Value'))"""),
('markdown', """### Independent Verification of the Scores

These checks do not rely on PLINK's scoring: effect alleles are counted directly from the `.bed` bytes for everyone. PLINK's six-digit output can differ by a tiny rounding amount."""),
('code', """def hand_scores(bfile, weights):
    '''Sum beta x effect-allele copies from the raw .bed bytes.'''
    people = len(pd.read_csv(f'{bfile}.fam', sep=r'\\s+', header=None))
    bim = read_bim(bfile)
    row = pd.Series(bim.index, index=bim['SNP'])
    bed = Path(f'{bfile}.bed').read_bytes()
    width = (people + 3) // 4
    total = np.zeros(people)
    for snp, a1, beta in zip(weights['SNP'], weights['A1'], weights['BETA']):
        packed = np.frombuffer(bed, np.uint8, width, 3 + row[snp] * width)
        codes = np.stack([(packed >> shift) & 3 for shift in (0, 2, 4, 6)], axis=1).ravel()[:people]
        assert not (codes == 1).any()
        first_allele = np.select([codes == 0, codes == 2], [2, 1], default=0)  # copies of .bim column 5
        total += beta * (first_allele if a1 == bim.loc[row[snp], 'R1'] else 2 - first_allele)
    return total

hand_reference5 = pd.Series(0.0, index=reference_prs5.index)
hand_participant5 = 0.0
for chrom, weights in weights5.groupby('CHR'):
    hand_reference5 += hand_scores(work5 / f'reference_chr{chrom}', weights)
    hand_participant5 += hand_scores(work5 / f'participant_chr{chrom}', weights)[0]
reference_difference5 = (hand_reference5 - reference_prs5).abs().max()
participant_difference5 = abs(hand_participant5 - participant_prs5)
assert reference_difference5 < 1e-5 and participant_difference5 < 1e-5
assert participant_id not in reference_prs5.index and len(reference_prs5) == 488
print(f'Largest difference from PLINK: reference {reference_difference5:.2e} | participant {participant_difference5:.2e}')"""),
('markdown', """### Height PRS Distribution

Reference scores with a kernel density estimate, and the participant's score on the same scale."""),
('code', """mean5, sd5 = position5['Reference_Mean'], position5['Reference_SD']
fig5, ax5 = plt.subplots(figsize=(10, 5.5))
ax5.hist(reference_prs5, bins=24, density=True, alpha=0.5, color='steelblue', edgecolor='white',
         label=f'1000 Genomes EUR reference (n={len(reference_prs5)})')
grid5 = np.linspace(min(reference_prs5.min(), participant_prs5) - 0.3 * sd5, max(reference_prs5.max(), participant_prs5) + 0.3 * sd5, 400)
ax5.plot(grid5, gaussian_kde(reference_prs5)(grid5), color='steelblue', label='Empirical KDE')
ax5.axvline(participant_prs5, color='crimson', linewidth=2.2, label=f'{participant_id} PRS = {participant_prs5:.5f}')
ax5.axvline(mean5, color='black', linestyle='--', label=f'Reference mean = {mean5:.5f}')
ax5.set(xlabel='Height polygenic score (PRS)', ylabel='Density', title=f'Height PRS Distribution: Public Participant {participant_id}')
ax5.text(0.02, 0.97, f"{participant_id}: percentile {position5['Target_Percentile']:.2f}\\nz-score {position5['Target_Z_Score']:.3f}\\n{len(weights5):,} SNPs; p ≤ 0.001",
         transform=ax5.transAxes, va='top', bbox={'facecolor': 'white', 'alpha': 0.9, 'edgecolor': 'none'})
ax5.legend(loc='upper right', fontsize=9)
fig5.tight_layout(rect=(0, 0.06, 1, 1))
fig5.text(0.5, 0.025, 'Educational comparison using simulated-height GWAS weights; scores are not height in cm.', ha='center', fontsize=9)
fig5.savefig(w5 / 'participant_height_prs_distribution.png', dpi=180)
plt.show()"""),
('markdown', "**Findings:** see Week 5 in [weekly_tasks_results.md](../weekly_tasks_results.md)."),
]
