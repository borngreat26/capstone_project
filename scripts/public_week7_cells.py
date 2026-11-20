"""Add the public Week 7 report without executing earlier weekly tasks.

Run: python scripts/public_week7_cells.py
The public-notebook generator uses the same cells and reuses unchanged outputs.
"""
import copy
import re
import sys
from pathlib import Path

import nbformat
from nbclient import NotebookClient

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.utils import REPO_ROOT

OUTPUT = REPO_ROOT / 'notebooks/weekly_tasks_public.ipynb'
TAG = 'public-week7'
CELLS = [
    ('markdown', '''# Week 7 — Functional Annotation of cis-eQTLs Using ENCODE ATAC-seq

Can the cis-eQTL variants and genes from Part 1 be found in experimentally measured
open chromatin? We reuse the completed genome-wide associations and annotate their
positions; no eQTL regressions are rerun here.

**Background.** Accessible chromatin is DNA available to regulatory proteins. ATAC-seq
uses transposase insertion to measure accessibility. The
[Roadmap Epigenomics study (2015)](https://pubmed.ncbi.nlm.nih.gov/25693563/)
and [ENCODE study (2020)](https://www.nature.com/articles/s41586-020-2493-4)
show how epigenomic measurements help identify candidate regulatory regions.

Our expression measurements come from lymphoblastoid cell lines (LCLs), derived from
B cells. Primary B-cell ATAC-seq is biologically relevant, but primary B cells and
EBV-transformed LCLs can have different accessibility patterns. Genomic overlap alone
does not establish a causal variant or a regulatory connection to a gene.'''),
    ('markdown', '''### Step 1: Verify the ENCODE File and Genome Assembly

The requested file is [ENCFF421XIL](https://www.encodeproject.org/files/ENCFF421XIL/)
from [ENCSR603LVR](https://www.encodeproject.org/experiments/ENCSR603LVR/):
primary B-cell ATAC-seq pseudoreplicated peaks in **BED6+4 narrowPeak** format.
The extra columns contain signal, significance and summit annotations.

The original peaks are **GRCh38**, while Part 1 and GENCODE v12 use **GRCh37**.
The existing BED's uncompressed MD5 was checked against official ENCODE metadata,
then UCSC liftOver with the `hg38ToHg19` chain and `-minMatch=0.95` converted the
intervals before overlap analysis. Original signal annotations were preserved;
original summit offsets are not treated as lifted summit coordinates.

These cells read the saved public results in `results/week7/`. To recreate the full
annotation from the saved Part 1 Table A, run these commands from the repository root:

```bash
python scripts/fetch_part2_week7.py
python scripts/run_week7.py
```

Existing input files are retained. Download provenance, assembly conversion and
overlap calculations are implemented in those scripts and `src/annotation.py`.'''),
    ('code', r'''from pathlib import Path
import json
import pandas as pd
from IPython.display import display, Image, Markdown

w7_root = Path.cwd()
if not (w7_root / 'src').is_dir():
    w7_root = w7_root.parent
w7 = w7_root / 'results/week7'
w7_required = ['provenance.json', 'snp_overlap_summary.tsv', 'gene_annotation.tsv',
               'case_studies.tsv', 'figures/A_snp_overlap.png',
               'figures/B_RPS26_accessibility.png', 'figures/C_gene_overlap.png']
if any(not (w7 / name).exists() for name in w7_required):
    raise FileNotFoundError('Recreate the public results with scripts/run_week7.py first.')

w7_provenance = json.loads((w7 / 'provenance.json').read_text())
display(pd.Series({key: w7_provenance[key] for key in
    ['Experiment', 'File', 'Biosample', 'Assay', 'Format', 'Original_Assembly',
     'Analysis_Assembly', 'Verified_Content_MD5', 'Method', 'Input_Peaks',
     'Mapped_Peaks', 'Lost_Peaks']}, name='Verified metadata').to_frame())
assert w7_provenance['Analysis_Assembly'] == 'GRCh37'
assert (w7_provenance['Mapped_Peaks'] + w7_provenance['Lost_Peaks']
        == w7_provenance['Input_Peaks'])
print(f"LiftOver retained {w7_provenance['Mapped_Peaks']:,} of "
      f"{w7_provenance['Input_Peaks']:,} peaks; {w7_provenance['Lost_Peaks']:,} were lost.")'''),
    ('markdown', '''### Step 2: Overlap Unique SNPs and SNP–Gene Associations

BED intervals are **zero-based, end-exclusive**: `[start, end)`. A one-based SNP
position `p` is the interval `[p−1, p)`, so it overlaps when `start ≤ p−1 < end`.
For example, BED `[100, 101)` contains one-based position 101, not position 100 or 102.

We analyze all tested cis-SNPs, SNPs in BH-significant associations, and SNPs in
Bonferroni-significant associations. Both significance definitions use the original
**genome-wide adjusted p < 0.05** across the Part 1 testing family. A SNP is counted
once per set even if it occurs in several associations. The association counts are
reported separately; overlap with multiple peaks does not multiply either count.'''),
    ('code', r'''w7_overlap = pd.read_csv(w7 / 'snp_overlap_summary.tsv', sep='\t')
assert w7_overlap['Overlapping_SNPs'].le(w7_overlap['Unique_SNPs']).all()
assert w7_overlap['Overlapping_Pairs'].le(w7_overlap['SNP_Gene_Pairs']).all()
assert ((w7_overlap['Percent']
         - 100 * w7_overlap['Overlapping_SNPs'] / w7_overlap['Unique_SNPs']).abs()
        < 1e-10).all()
display(w7_overlap.assign(Percent=w7_overlap['Percent'].map('{:.3f}%'.format)))
display(Image(filename=str(w7 / 'figures/A_snp_overlap.png'), width=850))'''),
    ('markdown', '''### Step 3: Annotate Gene Bodies and cis-Windows Separately

GENCODE v12 gene coordinates are one-based and inclusive. A gene from `start` to
`end` becomes BED `[start−1, end)`. We count peaks overlapping this gene body and,
separately, peaks overlapping the gene body extended by 500 kb on each side.

The denominator here is **all annotated autosomal protein-coding genes**, including
genes without testable expression in Part 1. A peak in a large cis-window is a
proximity annotation; it is not evidence that the peak regulates that gene.'''),
    ('code', r'''w7_genes = pd.read_csv(w7 / 'gene_annotation.tsv', sep='\t')
w7_gene_counts = pd.DataFrame({
    'Region': ['Gene body', 'Gene body ±500 kb'],
    'Annotated_genes': [len(w7_genes), len(w7_genes)],
    'Genes_with_peak': [int(w7_genes['Body_Peaks'].gt(0).sum()),
                        int(w7_genes['Cis_Window_Peaks'].gt(0).sum())]
})
w7_gene_counts['Percent'] = (100 * w7_gene_counts['Genes_with_peak']
                             / w7_gene_counts['Annotated_genes'])
display(w7_gene_counts.assign(Percent=w7_gene_counts['Percent'].map('{:.2f}%'.format)))
display(Image(filename=str(w7 / 'figures/C_gene_overlap.png'), width=800))'''),
    ('markdown', '''### Step 4: Examine RPS26 and Other Strong eQTL Genes

The case studies are the four genes with the smallest Part 1 lead-association p-values,
selected **regardless of accessibility overlap**. RPS26 is the primary example.
The table distinguishes lead-SNP overlap, overlapping significant unique SNPs, and
peak counts in gene bodies versus cis-windows.

The regional figure shows existing RPS26 association p-values, its lead SNP,
lifted ATAC-seq peak intervals and the gene body on the same GRCh37 coordinate axis.
It is a Matplotlib regional association plot, not output from LocusZoom software.'''),
    ('code', r'''w7_cases = pd.read_csv(w7 / 'case_studies.tsv', sep='\t').sort_values('Lead_P')
display(w7_cases.assign(Lead_P=w7_cases['Lead_P'].map('{:.3e}'.format)))
display(w7_genes.loc[w7_genes['Gene_ID'].isin(w7_cases['Gene_ID'])])
display(Image(filename=str(w7 / 'figures/B_RPS26_accessibility.png'), width=950))'''),
    ('markdown', '''### Findings and Interpretation

The summary below is generated from the saved tables. Higher overlap percentages
among significant SNPs are descriptive, not a formal enrichment test. LD, allele
frequency, local SNP density and distance from genes can influence these percentages.
Neither overlap nor a nearby accessible peak identifies a causal regulatory mechanism.'''),
    ('code', r'''w7_findings = []
for row in w7_overlap.itertuples(index=False):
    w7_findings.append(
        f"- **{row.Set}:** {row.Overlapping_SNPs:,}/{row.Unique_SNPs:,} unique SNPs "
        f"overlap peaks ({row.Percent:.2f}%); separately, "
        f"{row.Overlapping_Pairs:,}/{row.SNP_Gene_Pairs:,} SNP–gene pairs overlap.")
for row in w7_cases.itertuples(index=False):
    outcome = 'overlaps' if row.Lead_Overlaps else 'does not overlap'
    w7_findings.append(
        f"- **{row.Gene}:** lead SNP {row.Lead_SNP} {outcome} an ATAC-seq peak; "
        f"{row.BH_SNPs_Overlapping} BH-significant and "
        f"{row.Bonferroni_SNPs_Overlapping} Bonferroni-significant unique SNPs overlap.")
display(Markdown('\n'.join(w7_findings)))'''),
]


def week_7_range(cells):
    """Find Week 7, preserving any later weeks if they are present."""
    start = next((i for i, cell in enumerate(cells)
                  if cell.cell_type == 'markdown'
                  and re.match(r'^#{1,2} Week 7(?:\s|$)', cell.source)), len(cells))
    end = next((i for i in range(start + 1, len(cells))
                if cells[i].cell_type == 'markdown'
                and re.match(r'^#{1,2} Week \d+(?:\s|$)', cells[i].source)), len(cells))
    return start, end


def add_week_7(notebook, previous=None):
    cells = []
    for number, (kind, source) in enumerate(CELLS):
        create = (nbformat.v4.new_markdown_cell if kind == 'markdown'
                  else nbformat.v4.new_code_cell)
        cell = create(source, metadata={'tags': [TAG]})
        cell.id = f'public-week7-{number}'
        cells.append(cell)
    if previous is None and OUTPUT.exists():
        previous = nbformat.read(OUTPUT, as_version=4)
    saved = []
    if previous is not None:
        start, end = week_7_range(previous.cells)
        saved = previous.cells[start:end]
    if ([c.source for c in saved] == [c.source for c in cells]
            and any(c.get('outputs') for c in saved)):
        cells = copy.deepcopy(saved)
    else:
        report = nbformat.v4.new_notebook(cells=cells, metadata=copy.deepcopy(notebook.metadata))
        NotebookClient(report, timeout=120, kernel_name='python3', record_timing=False,
                       resources={'metadata': {'path': str(REPO_ROOT)}}).execute()
        cells = report.cells
    start, end = week_7_range(notebook.cells)
    notebook.cells[start:end] = cells


def main():
    notebook = nbformat.read(OUTPUT, as_version=4)
    add_week_7(notebook, previous=copy.deepcopy(notebook))
    nbformat.validate(notebook)
    nbformat.write(notebook, OUTPUT)
    print('Updated public Week 7; earlier weeks were not executed or changed.')


if __name__ == '__main__':
    main()
