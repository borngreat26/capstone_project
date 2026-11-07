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
| Part 2: multiple-disease PRS | **Not started.** The shared code and height workflow are ready, but disease weights are not ([why](docs/gwas_sumstats_inventory.md)) |
| Relationships between scores, ancestry PCs | Planned |
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
│   ├── quarter1_project.ipynb      Part 1: genome-wide cis-eQTL analysis
│   └── weekly_tasks_public.ipynb   Weeks 1–6, with private Week 5 results removed
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
│   └── check_public_files.py       Scans files for private data before publishing
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

This installs Python 3.14.10, the exact packages used here, and PLINK2 from Bioconda.

**PLINK2:** the code looks for it in this order: the `PLINK2_BIN` environment variable, `plink2` on your `PATH`, then a
local `./plink2_binary`. To use your own [PLINK2 download](https://www.cog-genomics.org/plink/2.0/):

```bash
export PLINK2_BIN=/path/to/plink2
```

**Tested on:** macOS on Apple Silicon, with PLINK2 v2.0.0-a.7.13 (Bioconda) and v2.0.0-b.1-dev, which gave identical
results. Linux should work through the same Conda environment, but has not been tested. Windows has not been tested. BCFtools and SAMtools (1.25 here) are only needed for the private Week 5
conversion, so they are not installed by default.

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
| GWAS summary statistics | Part 2 | `data/sumstats/*.sumstats.gz` | Public: [TCSC repository](https://github.com/TiffanyAmariuta/TCSC/tree/main/sumstats) |
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

This runs `notebooks/quarter1_project.ipynb` (3–5 minutes, about 1.3 GB of working files). Full results go
to `data/quarter1_results/part1/`, including Table A with all 6.6 million associations. Summary tables and
figures are copied to `results/part1/`:

- **Tables:** `table_B_gene_summary.tsv`, `table_C_chromosome_summary.tsv`, `correction_summary.tsv`,
  `excluded_genes.tsv`, `minor_allele_count_summary.tsv`, `top_gene_robustness.tsv`, and `verification_checks.tsv`.
- **Figures:** `gene_filtering.png`, `multiple_testing_summary.png`, `significant_by_chromosome.png`,
  `snps_tested_per_gene.png`, `genome_wide_lead_cis_eqtl.png`, `top_gene_locus.png`, and `RPS26_locuszoom.png`.

**Weekly tasks:** `notebooks/weekly_tasks_public.ipynb` is a record of the course tasks, with their original outputs. It
is not a one-command reproduction, because Weeks 3, 4, and 6 read small GENCODE extracts that were prepared by hand.
Part 1 above covers the same methods genome-wide and needs no such files.

**What cannot run without private data:** Week 5's comparison needs the private 23andMe file, so it cannot be
reproduced from this repository. The demo shows the same workflow on a public individual.

**Part 2** has not been implemented yet.

## Privacy and Limitations

- **Personal genotypes are excluded.** Week 5 used a personal 23andMe file, and its owner asked that it not be
  published. The raw file, converted files, individual genotype tables, and personal results stay on the
  analyst's computer. `.gitignore` excludes them, and `scripts/check_public_files.py` scans for them.
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
