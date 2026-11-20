"""Append/update independent result-reading sections without executing Part 1 again."""
import sys
from pathlib import Path
import nbformat
from nbclient import NotebookClient
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from src.utils import REPO_ROOT

TAG='part2-week7-report'

PROJECT_CONNECTIONS = r'''### How the analyses connect: two different kinds of beta

**Part 1 asks which nearby SNPs are associated with a gene's expression.** For each
eligible protein-coding gene across chromosomes 1–22, we fit a separate linear regression
for each SNP in the gene body ±500 kb:

$$E_{ig}=\alpha_{jg}+\beta^{\mathrm{eQTL}}_{jg}G_{ij}+\epsilon_{ijg}.$$

Here, $G_{ij}$ is person $i$'s dosage of the counted allele at SNP $j$, and $E_{ig}$ is
the supplied normalized expression for gene $g$. Each eQTL beta measures an expression
change per additional counted allele. The scan produces SNP–gene effects, p-values and
genome-wide multiple-testing corrections. It does not estimate disease effects.

**Week 4 used training-derived eQTL betas for an expression prediction exercise.** After
SNP selection, an expression score $S_{ig}=\sum_j\hat\beta^{\mathrm{eQTL}}_{jg}G_{ij}$
was calibrated on training data and evaluated with separate validation/test individuals.
This explains the remembered connection between our regressions and a weighted score.
That score predicts expression; those weights are not used for the Part 2 disease PRSs.

**Part 2 asks where a person's disease-associated genetic score falls in a reference group.**
Published disease GWASs already estimated allele effects using disease phenotypes. We
reuse documented log-odds betas (or convert documented odds ratios using $\beta=\ln(OR)$),
align the effect alleles to the genotypes, clump correlated SNPs and retain p≤0.001:

$$\mathrm{PRS}_{id}=\sum_{j\in J_d}\hat\beta^{\mathrm{GWAS}}_{jd}G^{\mathrm{effect}}_{ij}.$$

The genotype dosage is 0, 1 or 2 copies of the named effect allele. These weights come from
original RA, T2D and schizophrenia studies, not the original Z-only TCSC files and not our
Part 1 expression regressions. We do not fit a disease model on the reference individuals,
because disease outcomes are not available for this analysis.

Part 2 returns a raw PRS, reference mean/SD, a z-score and a reference percentile for a
comparison individual. We also examine cross-disease score correlations, selected-SNP
overlap, population PCs and changes after PC adjustment. A percentile is a rank among
the comparison people; neither a percentile nor the raw sum is a calibrated probability
of disease. Private comparisons use the same compatible SNP subset for target and reference.

**Week 7 adds biological context to Part 1.** It checks whether eQTL SNPs and gene regions
overlap B-cell accessible chromatin. It does not change PRS weights or establish a causal
chain from a SNP through expression to disease. Demonstrating that chain would require
additional analyses and evidence beyond this project.

### Disease PRS data and scoring choices'''

PAPER_FIGURE_INTRO = '''## Paper Figures — One Folder for LaTeX

All inventoried public Quarter 1 figures are collected in **`results/paper_figures/`**,
with captions, data sources, interpretation limits, writing notes and LaTeX examples.
Give that entire folder to the agent preparing the paper. Prefixes identify Part 1,
Part 2 and Week 7; PNG and PDF copies with the same stem represent one figure.

The cell below copies existing outputs and does not rerun regressions, PRS calculations
or figure generation. Original analysis figures remain intact. Private target figures
are excluded. Run `python scripts/make_paper_support.py` after changing the source figures
to refresh both the inventory and this collection.'''

PAPER_FIGURE_CODE = '''from pathlib import Path
import sys
paper_root = Path.cwd() if (Path.cwd() / 'src').is_dir() else Path.cwd().parent
if str(paper_root) not in sys.path:
    sys.path.insert(0, str(paper_root))
from scripts.export_paper_figures import export_paper_figures
paper_collection = export_paper_figures()
print(f"Paper figures: {paper_collection['Distinct_Figures']} distinct figures, "
      f"{paper_collection['Figure_Files']} image files")
print('Folder:', paper_collection['Folder'])
print('Captions and sources: figure_inventory.tsv; LaTeX examples: figures.tex')'''

def main():
    path=REPO_ROOT/'notebooks/quarter1_project.ipynb';nb=nbformat.read(path,as_version=4)
    nb.cells=[c for c in nb.cells if TAG not in c.metadata.get('tags',[])]
    cells=[nbformat.v4.new_markdown_cell('## Part 2 — Disease Polygenic Scores and Reference Structure\n\n'
        + PROJECT_CONNECTIONS + '''

The original TCSC files contained Z-scores rather than usable PRS betas and were removed
after inventorying at the user's request. The verified original
GWASs and selection reasons are listed below. All analyses use a fixed p≤0.001 cutoff,
EUR LD clumping (r²=0.1, 250 kb), GRCh37 coordinates, and allele-aligned log-odds weights.
Raw scores, reference-standardized scores, and percentiles describe different quantities;
none is a calibrated disease probability. Private comparisons remain in ignored local files.

The public PCA uses QC-filtered, LD-pruned autosomal genotypes. PCs describe within-EUR
structure; clusters are unlabeled because validated population metadata were not supplied.
PC adjustment is an exploratory sensitivity analysis. Correlations and SNP overlap do not
establish causality, pleiotropy, or a shared disease mechanism.

Run `scripts/run_disease_prs.py`, `scripts/run_pca.py` and `scripts/analyze_prs.py` from the
repository root to recreate the results. These cells read saved outputs and do not refit Part 1.'''),
    nbformat.v4.new_code_cell('''from pathlib import Path
import pandas as pd
from IPython.display import display, Image
report_root = Path.cwd() if (Path.cwd() / 'src').exists() else Path.cwd().parent
for name in ['selected_diseases.tsv', 'prs_summary.tsv', 'pairwise_correlations.tsv',
             'snp_overlap.tsv', 'prs_pc_associations.tsv']:
    report_path = report_root / 'results/part2' / name
    if report_path.exists():
        print(name)
        display(pd.read_csv(report_path, sep='\\t'))
    else:
        print('Run the Part 2 scripts to generate', name)
for name in ['figures/prs_distributions.png', 'figures/prs_correlations.png', 'pca/reference_PCA.png']:
    report_path = report_root / 'results/part2' / name
    if report_path.exists(): display(Image(filename=str(report_path), width=760))'''),
    nbformat.v4.new_markdown_cell('''## Week 7 — Functional Annotation of cis-eQTLs Using ENCODE ATAC-seq

Accessible chromatin is DNA available for regulatory protein binding. ATAC-seq measures
transposase-accessible regions. [Roadmap Epigenomics (2015)](https://pubmed.ncbi.nlm.nih.gov/25693563/)
and [ENCODE (2020)](https://www.nature.com/articles/s41586-020-2493-4) provide the motivation
for annotating association signals with cell-context epigenomic measurements.

ENCFF421XIL (ENCSR603LVR) contains primary B-cell ATAC-seq peaks. B cells are relevant
to B-cell-derived LCL expression, but primary B cells and EBV-transformed LCLs are distinct
experimental contexts. The original BED6+4 peaks use GRCh38. We verify the official
uncompressed checksum, then use UCSC liftOver to GRCh37 before comparing with Part 1.

A one-based SNP p overlaps a BED interval [start,end) when start≤p−1<end. We separately
count unique SNPs and SNP–gene associations for all tested, BH-significant, and
Bonferroni-significant associations. GENCODE v12 gene-body and gene-body ±500 kb window
overlaps are separate annotations; the gene denominator includes unexpressed genes.

The case studies are the four strongest Part 1 genes, selected regardless of overlap.
RPS26's lead-SNP result and the negative findings for other lead SNPs are reported below.
Overlap is evidence of genomic co-location, not a causal variant or regulatory connection.
Differences in overlap percentages are descriptive; no formal enrichment claim is made.
Run `python scripts/run_week7.py` to reproduce this section from saved Part 1 associations.'''),
    nbformat.v4.new_code_cell('''import json
week7 = report_root / 'results/week7'
if (week7 / 'provenance.json').exists():
    provenance = json.loads((week7 / 'provenance.json').read_text())
    print(f"Peaks: {provenance['Input_Peaks']:,} original; {provenance['Mapped_Peaks']:,} lifted; {provenance['Lost_Peaks']:,} lost")
    for name in ['snp_overlap_summary.tsv', 'case_studies.tsv']:
        print(name)
        display(pd.read_csv(week7 / name, sep='\\t'))
    annotation = pd.read_csv(week7 / 'gene_annotation.tsv', sep='\\t')
    display(pd.DataFrame({'Annotated_genes': [len(annotation)],
                          'Gene_body_overlap': [int(annotation.Body_Peaks.gt(0).sum())],
                          'Cis_window_overlap': [int(annotation.Cis_Window_Peaks.gt(0).sum())]}))
    for name in ['A_snp_overlap.png', 'B_RPS26_accessibility.png', 'C_gene_overlap.png']:
        display(Image(filename=str(week7 / 'figures' / name), width=760))
else:
    print('Run scripts/fetch_part2_week7.py and scripts/run_week7.py first.')''')]
    cells.extend([nbformat.v4.new_markdown_cell(PAPER_FIGURE_INTRO),
                  nbformat.v4.new_code_cell(PAPER_FIGURE_CODE)])
    for cell in cells:cell.metadata['tags']=[TAG]
    extension=nbformat.v4.new_notebook(cells=cells,metadata=nb.metadata)
    NotebookClient(extension,timeout=120,kernel_name='python3',record_timing=False,resources={'metadata':{'path':str(REPO_ROOT)}}).execute()
    nb.cells.extend(extension.cells);nbformat.write(nb,path)
    from scripts.update_part2_notebook import refresh_notebooks
    refresh_notebooks()
    if (REPO_ROOT/'results/expression_prediction/pilot150_v2/performance_summary.json').exists():
        from scripts.update_research_notebook import update_notebooks
        update_notebooks()
    print('Updated only appended report sections; original Part 1 cells were not executed.')
if __name__=='__main__':main()
