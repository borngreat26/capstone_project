# GWAS Summary Statistics: What the Local Files Support

`scripts/inspect_sumstats.py` checked all 90 `.sumstats.gz` files in `data/sumstats/`, from the
[TCSC repository](https://github.com/TiffanyAmariuta/TCSC/tree/main/sumstats). The results for each file are in
[gwas_sumstats_inventory.tsv](gwas_sumstats_inventory.tsv). These files are not copied into this repository.

## Main Finding: These Files Cannot Be Used as PRS Weights Directly

**All 90 files have Z-scores and sample sizes only. None has a beta, odds ratio, standard error, p-value,
allele frequency, chromosome, or position.** Z-scores are not genetic effect sizes, because a Z-score also
depends on sample size and allele frequency. No PRS has been calculated from these files.

| Property | Result across the 90 files |
|---|---|
| Columns | `SNP, A1, A2, N, CHISQ, Z` (68 files); `SNP, A1, A2, Z, N` (9); `SNP, A1, A2, N, Z` (7); `SNP, A2, A1, Z, N` (6) |
| SNP IDs | rsIDs |
| Genome build | Not stated, and no positions to check. Matching to the GRCh37 panel has to use rsIDs |
| Rows | 352,680 to 1,201,412 |
| Share of panel SNPs covered | 29.6% to 99.7% (median about 99%) |
| Alleles matching the panel, among SNPs in both | 100% |
| Strand-ambiguous SNPs (A/T or C/G) | None among SNPs in the panel |
| Missing Z-scores | None |

**Handle with care:**

- **Allele order:** six files list `A2` before `A1`, so columns must be read by name, never by position.
- **Effect direction:** `A1` is presumed to be the allele the Z-score refers to, which is the LD Score
  regression convention. The files do not state this, so check it against the original study.
- **Details not in the files:** ancestry, case and control counts, and usage limits have to come from each
  original publication. This inventory did not verify them.

## Options for Part 2

| Option | How it works | Trade-offs |
|---|---|---|
| **A. Original GWAS downloads (recommended)** | Get betas or odds ratios, standard errors, p-values, and positions from each published study | The most defensible choice, and it matches the Week 5 workflow. Takes more work, and access varies by study |
| **B. Approximate from Z** | $\hat\beta \approx Z / \sqrt{2p(1-p)(N+Z^2)}$, with $p$ from a reference panel and p-values from Z | Works with local files. Effects are on a standardized scale, depend on the reference frequencies, and need effective N for disease traits. Needs instructor approval |

P-values can be calculated exactly from Z, so clumping is possible with these files. The missing piece is the weights.

## Suggested Disease Traits

These files have good panel coverage and are named for diseases. Counts come from the inventory. The final choice
depends on finding effect estimates and on the instructor's advice.

| Disease | File | Median N | Panel coverage |
|---|---|---:|---:|
| Schizophrenia | `PASS_Schizophrenia_Pardinas2018` | 94,437 | 82.7% |
| Type 2 diabetes | `PASS_Type_2_Diabetes` | 54,337 | 81.4% |
| Rheumatoid arthritis | `PASS_Rheumatoid_Arthritis` | 38,242 | 94.5% |
| Inflammatory bowel disease | `PASS_IBD_deLange2017` | 59,957 | 96.7% |
| Alzheimer's disease | `PASS_Alzheimers_deRojas2021` | 409,435 | 99.7% |
| Atrial fibrillation | `PASS_AtrialFibrillation_Nielsen2018` | 1,030,836 | 99.7% |
| Major depression | `PASS_MDD_Howard2019` | 380,296 | 97.4% |
| Bipolar disorder | `PASS_BIP_Stahl2019` | 46,582 | 99.5% |

Schizophrenia and type 2 diabetes cover less of the panel, so more of their selected SNPs may be lost.
Median N is not an effective sample size for a case-control study.
