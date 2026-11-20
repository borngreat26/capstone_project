# Genetic Variation, Gene Expression, and Polygenic Scores

DSC 180A capstone, Quarter 1 Project. This repository studies how common genetic variants relate to gene
expression and to disease-associated traits.

- **Part 1: Genome-wide cis-eQTL analysis.** Find associations between SNPs and the expression of nearby
  protein-coding genes, using 1000 Genomes genotypes and lymphoblastoid cell line (LCL) expression data.
- **Part 2: Disease polygenic risk scores (PRS).** Use disease GWAS summary statistics to build genetic scores for
  1000 Genomes individuals, and compare one person's score with the reference distribution.

**To try it quickly,** go to [Run the Demo](#run-the-demo). It needs no private data or large downloads.

## Current Progress

| Component | Status |
|---|---|
| Weeks 1–6 course tasks | Done ([public copy](notebooks/weekly_tasks_public.ipynb), [results](weekly_tasks_results.md)) |
| Part 1: genome-wide cis-eQTLs | **Done:** 14,661 protein-coding genes, 6,597,508 SNP–gene tests ([results](quarter1_part1_results.md)) |
| Public demo and tests | Done |
| Part 2: disease PRS | Completed for three verified original GWASs; see the results and limitations below |
| Relationships between scores, ancestry PCs | Completed on 489 public EUR individuals |
| Week 7: functional annotation | Completed after verified GRCh38-to-GRCh37 conversion |
| Research paper | Planned for December 4 ([plan](docs/final_submission_plan.md)) |

The genotype panel and expression file provide 14,661 testable genes, not the roughly 20,000 protein-coding genes
in the genome. The [Part 1 results](quarter1_part1_results.md) explain the difference.

## Repository Contents

```text
├── README.md
├── environment.yml                 Python packages and PLINK2
├── quarter1_part1_results.md       Part 1 findings and interpretation
├── weekly_tasks_results.md         Weeks 1–6 findings and takeaways
├── notebooks/
│   ├── quarter1_project.ipynb      Local working analysis, including personal target results
│   ├── quarter1_project_public.ipynb  Sanitized reference-only report
│   └── weekly_tasks_public.ipynb   Weeks 1–6, with Week 5 run on a public participant
├── src/
│   ├── eqtl.py                     cis-eQTL regressions and multiple-testing correction
│   ├── locuszoom.py                Builds the Week 6 LocusZoom page
│   ├── prs.py                      QC, harmonization, clumping, scoring, percentile and z-score
│   └── utils.py                    Paths and finding PLINK2
├── scripts/
│   ├── run_demo.py                 Runs the public demo
│   ├── run_full_analysis.py        Runs Part 1 on the full data
│   ├── build_demo_data.py          How the demo data were made (maintainers only)
│   ├── inspect_sumstats.py         Inventory of the GWAS files
│   ├── render_locuszoom.py         Renders the Week 6 LocusZoom page to a PNG (needs Playwright)
│   ├── make_public_notebook.py     Creates the public weekly notebook (maintainers only)
│   └── public_week5_cells.py       The public notebook's Week 5, on a public participant
├── data/demo/                      Small public demo dataset, with its own README
├── results/part1/                  Part 1 summary tables and figures
├── results/demo/                   Demo outputs, created when you run it
├── tests/                          Automated checks
└── docs/                           GWAS inventory, final-submission plan, release checklist
```

## Installation

You need [Conda](https://docs.conda.io/en/latest/miniconda.html). From the repository root:

```bash
conda env create -f environment.yml
conda activate capstone
```

This installs Python 3.12, constrained scientific packages, PLINK2, BCFtools, bedtools and UCSC liftOver from Conda-forge/Bioconda. SciPy 1.16 is required for the existing vectorized regression API. The final environment was tested in a fresh temporary prefix.

**BCFtools:** included for genotype conversion utilities. The public demo needs no private conversion. SAMtools is only an optional helper for preparing FASTA references in the legacy private workflow.

**PLINK2:** the code looks for it in this order: the `PLINK2_BIN` environment variable, `plink2` on your `PATH`, then a
local `./plink2_binary`. To use your own [PLINK2 download](https://www.cog-genomics.org/plink/2.0/):

```bash
export PLINK2_BIN=/path/to/plink2
```

**Tested on:** fresh Conda on macOS arm64 and a Linux Docker container.
Windows and Linux x86_64 have not been independently tested in this run. The Docker image
uses Linux packages and never copies the Apple Silicon executable.

## Run the Demo

From the repository root, with the environment active:

```bash
python scripts/run_demo.py
```

It takes a few seconds and writes to `results/demo/`:

| Output | What it shows |
|---|---|
| `demo_eqtl_associations.tsv` | SNP–gene results for 12 genes, with BH and Bonferroni p-values |
| `demo_eqtl_gene_summary.tsv` | One row per gene |
| `demo_eqtl_top_gene.png` | Association plot for the strongest gene |
| `demo_snp_counts.tsv` | SNPs kept at each step |
| `demo_gwas_harmonized.tsv` | GWAS after QC and harmonization |
| `demo_weights.tsv` | Final SNPs, effect alleles, and weights |
| `demo_reference_prs.tsv` | Score for each reference person |
| `demo_prs_summary.tsv` | The demo individual's score, percentile, and z-score |
| `demo_prs_distribution.png` | Reference distribution, with a line at the demo individual's score |

To run the tests:

```bash
python -m unittest discover -s tests -v
```

Tests that need the full data skip themselves when it is missing.

**What the demo is:** a public 1000 Genomes person (`demo_individual`, chosen with a fixed seed and excluded from
the reference group) goes through the same functions as the real analysis, on a small region of chromosome 22. The
genotypes and expression values are real public data, but **the GWAS is simulated**. It shows that the computation
works. It says nothing about anyone's disease risk, and its outputs differ from the private analysis by design.
See [data/demo/README.md](data/demo/README.md).

## Data

Everything below stays out of Git. Store each file at the path shown, relative to the repository root.

| Data | Used by | Location | Source and access |
|---|---|---|---|
| 1000 Genomes EUR genotypes (GRCh37, 489 people) | Part 1, Weeks 1–6 | `data/LDREF/1000G.EUR.{1..22}.{bed,bim,fam}` | Public: [LDREF.tar.bz2](https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2) |
| LCL expression | Part 1, Weeks 2–6 | `data/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt` | Public GEUVADIS study: [E-GEUV-1](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-GEUV-1). Supplied through the course. Decompress it if your copy ends in `.gz` |
| GENCODE v12 annotation (GRCh37) | Part 1 | `data/quarter1_results/part1/gencode.v12.annotation.gtf.gz` | Public: [GENCODE](https://www.gencodegenes.org/human/release_12.html). The notebook downloads it |
| Height GWAS (simulated, from the PRS tutorial) | Week 5 | `data/Height.gwas.txt.gz` | [Course link](https://drive.google.com/file/d/1RWjk49QNZj9zvJHc9X_wyZ51fdy6xQjv/view) |
| GRCh37 reference FASTA | Week 5 | `data/reference/human_g1k_v37.autosomes.fa` | Public: [1000 Genomes](https://1000genomes.s3.amazonaws.com/technical/reference/human_g1k_v37.fasta.gz), chromosomes 1–22 only |
| Active disease GWAS effects | Part 2 | `data/gwas_original/` | Original RA/T2D GWAS Catalog files and author-hosted schizophrenia file; `scripts/fetch_part2_week7.py` records exact URLs |
| Personal 23andMe genotypes | Week 5 | Local only | **Private and not available.** See [Privacy](#privacy-and-limitations) |

To get the genotype panel:

```bash
mkdir -p data
curl -L -o data/LDREF.tar.bz2 https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2
tar -xjf data/LDREF.tar.bz2 -C data
```

**Conventions:** all genotype coordinates and gene annotations use GRCh37. Run commands from the repository root.
Random seeds are fixed: 42 for the Week 4 split and for the demo.

## Run the Full Analysis

**Part 1** needs only the public genotype panel and expression file:

```bash
python scripts/run_full_analysis.py
```

This runs the Part 1 cells of `notebooks/quarter1_project_public.ipynb` (3–5 minutes, about 1.3 GB of working files). Full results go
to `data/quarter1_results/part1/`, including Table A with all 6.6 million associations. Summary tables are
copied to `results/part1/`:

- **Tables:** `table_B_gene_summary.tsv`, `table_C_chromosome_summary.tsv`, `correction_summary.tsv`,
  `excluded_genes.tsv`, `minor_allele_count_summary.tsv`, `top_gene_robustness.tsv`, and `verification_checks.tsv`.
- **Figures:** the notebook saves all 7 plots to `results/part1/figures/`, numbered in notebook order:
  `01_genes_kept_at_each_filtering_step.png`
  `02_multiple_testing_results.png`
  `03_significant_cis_eqtls_by_chromosome.png`
  `04_snps_tested_per_gene.png`
  `05_manhattan_strongest_cis_eqtl_per_gene.png`
  `06_RPS26_association_and_expression_by_genotype.png`
  `07_RPS26_locuszoom.png`

**Weekly tasks:** `notebooks/weekly_tasks_public.ipynb` is a record of the course tasks, with their original outputs. It
is not a one-command reproduction, because Weeks 3, 4, and 6 read small GENCODE extracts that were prepared by hand.
Part 1 above covers the same methods genome-wide and needs no such files.

**Week 5 in public:** the original comparison used a private 23andMe file, so it cannot be reproduced from this
repository. The public notebook runs the same workflow on NA20792, a randomly chosen public participant, using the full
panel and the height GWAS. Those results are that participant's, and the small demo uses the same person.

**Part 2 and Week 7:** use the commands below. They reuse the saved Part 1 results and do not rerun its regressions.

## Privacy and Limitations

- **Personal genotypes are excluded.** Week 5 used a personal 23andMe file, and its owner asked that it not be
  published. The raw file, converted files, individual genotype tables, and personal results stay on the
  analyst's computer. `.gitignore` excludes them.
- **The demo uses a public reference individual** so anyone can run the workflow without private data. That person
  is not a stand-in for the private individual's results, so the outputs differ.
- **A PRS percentile is not an absolute disease risk.** It is a rank within one reference group. Turning a score
  into a probability or relative risk requires a model calibrated on real outcomes.
- **The reference group and SNP coverage affect interpretation.** The panel is 489 European-ancestry people, and
  scores may not transfer to other ancestries. Consumer genotyping chips measure only some GWAS variants, which
  limits a score even when everyone is scored on the same SNPs.
- **cis-eQTL results are associations.** Linkage disequilibrium means a significant SNP may only tag the true
  causal variant, and the model has no covariates for ancestry or technical factors. See the
  [Part 1 results](quarter1_part1_results.md).

## References

**Course**
- [DSC 180A Quarter 1 Project](https://dsc-capstone.org/2026-27/assignments/projects/q1/)
- [Genetic risk prediction capstone site](https://tiffanyamariuta.github.io/capstone-genetic-risk-prediction/)
- [PRS tutorial (Choi et al.)](https://choishingwan.github.io/PRS-Tutorial/)

**Data**
- [1000 Genomes Project](https://www.internationalgenome.org/)
- [GEUVADIS expression, E-GEUV-1](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-GEUV-1)
- [GENCODE release 12](https://www.gencodegenes.org/human/release_12.html)
- [TCSC summary statistics](https://github.com/TiffanyAmariuta/TCSC/tree/main/sumstats)

**Software**
- [PLINK 2.0](https://www.cog-genomics.org/plink/2.0/)
- [BCFtools and SAMtools](https://www.htslib.org/)

## Part 2, PCA and Week 7

The original 90 TCSC downloads are LD Score regression inputs: they have Z-scores, not
PRS effect coefficients. No Z-score was used as a beta. The expanded inventory is
[`results/part2/gwas_inventory.tsv`](results/part2/gwas_inventory.tsv). The selected original
studies and their limitations are in [`selected_diseases.tsv`](results/part2/selected_diseases.tsv).
RA and T2D use targeted arrays with limited, nonrandom coverage. Schizophrenia uses the
author-hosted Pardiñas 2018 CLOZUK+PGC2 meta-analysis with documented ORs/effect alleles,
a verified original MD5, and GRCh37 positions checked against 1,147,331 independent panel
rsID anchors. These scores describe disease-associated genetic variation and are not
validated clinical risk predictors. The author portal does not explicitly label the schizophrenia
build; its GRCh37 assignment is empirically verified, rather than assumed.

From the repository root, with Conda active:

```bash
python scripts/inventory_part2.py
python scripts/fetch_part2_week7.py
python scripts/run_disease_prs.py
python scripts/run_pca.py
python scripts/analyze_prs.py
python scripts/run_week7.py
python scripts/make_paper_support.py
python scripts/update_analysis_notebook.py
```

The 90 unused Z-only TCSC archives were removed at the user's request on November 20, 2026,
freeing about 1.06 GB. Their historical inventory and deletion checksums are retained in
`results/part2/gwas_inventory.tsv` and `unused_gwas_cleanup.tsv`; they are not required or
downloaded for the active analysis. The inventory script preserves those historical rows.
The selected-disease analysis only needs the original files fetched by `fetch_part2_week7.py`
and the EUR panel. The fetcher
preserves existing inputs. For a public-only run use `run_disease_prs.py --skip-private`.
If the existing converted private target file is present, the default run writes its common-SNP
comparisons exclusively to ignored `data/private_part2/`; those files must stay local.
Original genotype files are read only. Original schizophrenia effects are for noncommercial
research and cannot be redistributed without meeting the author conditions; the archive and
all extracted GWAS weights are ignored by Git and excluded from Docker. Public and private-comparison scoring sets are separate,
and every person compared within either set uses that set's identical nonmissing variants.

PCA uses all 22 autosomes, MAF≥0.05, missingness≤0.02, HWE p≥10⁻⁶, and LD pruning
(200-variant windows, 50-variant steps, r²=0.2). PC1–PC2 adjustment is exploratory;
PCs describe structure within the EUR panel, not global ancestry. No population labels
were inferred from sample IDs. Main results are under `results/part2/`: raw scores,
reference z-scores, count audits, summary statistics, Pearson/Spearman correlations,
SNP-set overlap, PCA variance fractions and PC-adjusted correlations.

Week 7 reads the completed Part 1 Table A and GENCODE v12 GTF. It verifies the existing
ENCFF421XIL content checksum against official metadata and converts GRCh38 peaks to GRCh37
using UCSC liftOver (`-minMatch=0.95`) and the UCSC hg38ToHg19 chain. An existing file is
not downloaded again. The original BED and its metadata remain under ignored `data/encode/`.
A SNP at one-based p is the half-open interval [p−1,p). Gene-body and ±500 kb overlaps
are separate annotations. Peak counts, losses, checksums, unique-SNP counts and association
counts are saved under `results/week7/`. Detailed association and peak-link files remain local;
summary tables and figures can be reviewed for publication.

**For the LaTeX paper, give the writing agent the entire
[`results/paper_figures/`](results/paper_figures/) folder.** This ignored local collection includes
the personal target plots as well as Part 1, Part 2 and Week 7 figures, caption/source inventory, methods/results notes and ready-to-adapt
LaTeX examples. The final cell of `notebooks/quarter1_project.ipynb` refreshes the collection
from existing outputs. It does not rerun the analyses. The same refresh is available with
`python scripts/export_paper_figures.py --include-private`. Run
`python scripts/private_part2_report.py` to refresh target comparisons from saved scores.
The default exporter and `scripts/make_paper_support.py` refresh only
`results/paper_figures_public/`, which excludes personal target results.

`quarter1_project.ipynb` is the local working notebook; `quarter1_project_public.ipynb`
is the sanitized reference-only copy. Run `python scripts/update_part2_notebook.py` to
refresh both report sections without rerunning Part 1. The working file was previously
tracked by Git: its ignore rule alone cannot protect it. It is excluded from the release
allowlist, and the privacy checker blocks a Git index containing personal report cells.
Do not use `git add -a` for a public release; removing that local file from tracking can be
reviewed separately. No staging or history changes were made automatically.

New figures are 320 dpi PNG plus vector PDF. Existing Part 1 figures were preserved at their
original resolution. [`results/figure_inventory.tsv`](results/figure_inventory.tsv) lists figures,
sources, captions and limitations. [`docs/paper_figures.tex`](docs/paper_figures.tex) supplies
LaTeX figure environments and `Figure~\ref{...}` examples. No Overleaf document was changed.
[`docs/quarter1_part2_week7_methods_results.md`](docs/quarter1_part2_week7_methods_results.md)
contains results generated from the actual tables.

## Docker and reproducibility verification

The expanded training-only expression-prediction workflow is described in
[`docs/expanded_prediction_methods_results.md`](docs/expanded_prediction_methods_results.md).
It preserves Week 4's noncoding FLVCR1-AS1 example and adds a frozen 150-gene
protein-coding pilot. Its gene-level performance, negative R² values, baselines,
training-feature associations and resource estimates are under
`results/expression_prediction/pilot150_v2/` (the earlier `pilot150` draw is kept frozen but
superseded; see the notes). Per-gene fitting intermediates remain local.
The full 14,661-gene prediction run has **not** been launched.

```bash
python scripts/verify_week4_prediction.py
python scripts/run_expression_prediction.py --stage fit
python scripts/run_expression_prediction.py --stage evaluate
python scripts/verify_prediction_experiment.py
python scripts/analyze_expression_prediction.py
python scripts/expand_prs_exploration.py
python scripts/write_prediction_notes.py
python scripts/make_paper_support.py
python scripts/update_research_notebook.py
```

Training effects, LD clumping and calibration use only the 275 training people;
34 validation people choose the cutoff; 35 test people evaluate frozen models.
Prediction R² is `1-SSE/SST`, distinct from Pearson r². Repeat evaluation reuses saved
metrics. `--run` names isolate immutable experiments, and `--max-models` supports batch
checkpoints. A genome-wide run requires the separate explicit `--all-genes` option;
review the pilot resource estimate before using it. This is a stratified pilot, not a
genome-wide performance benchmark.

The new disease PRS extension compares reference rankings across five fixed thresholds
and describes SNP-level/pair-level eQTL overlap against harmonized/tested backgrounds.
Primary PRSs and personal results are unchanged. It is not formal colocalization.

The working notebook contains separate single-gene prediction, multi-gene prediction
and disease PRS exploration sections. Its final export also copies the complete local
paper handoff to **`paper/figures/`**, preserving existing figure locations. That directory
is ignored because it can include personal target plots; the reference-only public
collection remains `results/paper_figures_public/`.

```bash
docker build -t capstone-q1 .
docker run --rm capstone-q1
docker run --rm capstone-q1 python -m unittest discover -s tests -v
```

The image uses a Linux Conda environment. Explicit COPY allowlists and `.dockerignore`
exclude full genotypes, GWAS archives, private files, credentials, notebooks and local binaries.
The default demo writes inside the disposable container. Mount only the intended output/data
directories when retaining results. To run a full public-data script, mount local public inputs
at `/workspace/data` and an output directory at `/workspace/results`, then pass the relevant
`python scripts/...` command. Do not mount private data into a public demonstration.

The same demo runs eQTL regressions, P+T scoring, a PRS distribution and synthetic ATAC
interval overlaps. It contains a public 1000 Genomes individual excluded from its reference
group; simulated GWAS weights and synthetic peaks are labeled. It does not reproduce
personal numerical results.

For the full public analysis, acquire the panel, GEUVADIS expression and GENCODE v12 inputs
listed in the Data table above, fetch the selected original GWASs and ENCODE resources, then
run Part 1 only if its saved Table A is absent. The remaining scripts reuse that table. The
notebook report updater executes only the appended result-reading cells, preserving Part 1.
Full data and public demo outputs are distinct; private comparisons cannot be reproduced
without the owner's local files. Review [`docs/release_checklist.md`](docs/release_checklist.md)
before staging any files.
