# Quarter 2 Feasibility: Can the Existing GWAS Data Support the Proposed Study?

**Project:** Integrating cis-eQTL Associations and Chromatin Accessibility to Prioritize Disease-Associated Regulatory Genes
**Audit date:** 2026-11-20. **Scope:** read-only audit. No Quarter 1 file was modified, nothing was downloaded or rerun genome-wide, and no individual-level genotypes are reported.

**How the numbers were produced.** A read-only profiling script read every GWAS file in full, the GRCh37 1000 Genomes EUR panel (`data/LDREF/1000G.EUR.*.bim`), the Part 1 cis-eQTL table, and the lifted-over B-cell ATAC peaks. Outputs:

- [`results/q2_planning/gwas_feasibility_table.tsv`](../results/q2_planning/gwas_feasibility_table.tsv): one row per GWAS
- [`results/q2_planning/gwas_lead_loci_500kb.tsv`](../results/q2_planning/gwas_lead_loci_500kb.tsv): lead SNP of each genome-wide significant locus
- [`results/q2_planning/gwas_locus_density.png`](../results/q2_planning/gwas_locus_density.png): SNP density at those loci

Confidence labels used throughout: **Confirmed** means read from the file, or recomputed from it. **Documented** means stated in a local metadata or source file. **Not available** means absent locally.

---

## 1. Executive Summary

- **Four GWAS files exist locally.** Three are real disease GWASs: rheumatoid arthritis (RA), type 2 diabetes (T2D) and schizophrenia (SCZ). The fourth is the simulated Week 5 height file from the PRS tutorial.
- The 90 TCSC/LDSC `.sumstats.gz` files had **Z-scores only**. They were deleted on 2026-11-20 and are not candidates.
- **All three disease files have real effect estimates and standard errors.** Effect allele, other allele, p-values and GRCh37 positions are present in each.
- **The Quarter 1 PRS pipeline did not use Z-scores as betas.** It took ORs to log(OR) and used log-odds betas as given. Effect alleles align correctly with PLINK dosages. I found no correctness errors, only a few minor QC gaps (Section 3).
- **Recommended primary disease: rheumatoid arthritis** (`GCST005569`).
  - It is the only immune disease available, and it is B-cell relevant.
  - It has beta and SE in GRCh37.
  - Its Immunochip design gives dense SNP coverage at RA loci: 8 of 12 genome-wide significant loci have ≥50 SNPs shared with the eQTL table within ±250 kb.
  - Its weaknesses: the genome-wide coverage is narrow (only 4.6% of eQTL-tested SNPs are covered), there are only 12 significant loci, and allele frequencies are absent.
- **Optional secondary: schizophrenia, as a tissue-mismatched contrast.** It is dense and complete (96.2% of eQTL-tested SNPs covered; 136 loci), but it is biologically unrelated to B cells. T2D is not suitable for colocalization, because its exome-array SNPs are too sparse (median of 6 shared SNPs per locus).
- **Formal colocalization (`coloc.abf`) is possible in principle, but not yet runnable.**
  - The data have the required inputs: beta, varbeta (SE²), the trait type, and matchable SNP IDs and alleles. The eQTL side also needs either `sdY` or MAF plus N, and both can be derived from local files.
  - R, `coloc` and `susieR` are **not installed**. Installing them is the first required step.
- **Recommendation: Option A, scoped narrowly.**
  - Question A (eQTL/ATAC enrichment) is the main quantitative result.
  - Question B runs `coloc.abf` at the RA Immunochip loci, plus ATAC annotation and gene ranking.
  - Fall back to Option B if the coloc preprocessing checks fail.

---

## 2. Existing GWAS Dataset Inventory

### 2.1 Where the files are

| Path | Status |
|---|---|
| `data/gwas_original/GCST005569.Build37.f.tsv.gz` (+ `.meta.yaml`, `.readme.txt`) | Active; RA |
| `data/gwas_original/GCST007517.Build37.f.tsv.gz` (+ `.meta.yaml`) | Active; T2D |
| `data/gwas_original/GCST006803.txt.gz` (+ `.md5`, `.source.html`) | Active; SCZ |
| `data/Height.gwas.txt.gz` | Week 5 teaching file, simulated (`README.md` line 128) |
| `data/sumstats/alltraits.txt`, `AllPairs_Cov.txt` | LDSC h² and rg metadata only; no SNP data |
| 90 × `PASS_*/UKB_*.sumstats.gz` | **Deleted 2026-11-20** (`data/sumstats/README.md`, `results/part2/unused_gwas_cleanup.tsv`). Columns were `SNP, A1, A2, (N, CHISQ,) Z`: no beta, SE, P or position |

Other files whose names contain "gwas" are not disease GWASs:

- `data/part2_work/*/chr*_gwas.tsv`: PRS clumping intermediates
- `data/prediction_work/week4_audit/training_gwas.tsv`: expression-prediction training
- `data/demo/demo_gwas.tsv`: demo
- `data/week6_results/FAM118A_summary_statistics.tsv`: eQTL

### 2.2 Verified properties

| Property | Rheumatoid arthritis | Type 2 diabetes | Schizophrenia | Height |
|---|---|---|---|---|
| Study | Eyre et al. 2012, Immunochip; PMID 23143596 | Mahajan et al. 2018, exome array, BMI-unadjusted; PMID 29632382 | Pardiñas et al. 2018, CLOZUK+PGC2; PMID 29483656 | Simulated (PRS tutorial) |
| Source URL (from `results/part2/gwas_inventory.tsv`) | EBI FTP, GCST005569 `harmonised/…-Build37.f.tsv.gz` | EBI FTP, GCST007517 `harmonised/…-build37.f.tsv.gz` | `walters.psycm.cf.ac.uk/clozuk_pgc2.meta.sumstats.txt.gz` (MD5 verified) | Course Google Drive link |
| Format | gzip TSV, GWAS Catalog "formatted" file (not the GRCh38 `.h` file) | gzip TSV, GWAS Catalog formatted file | gzip TSV, author file | gzip TSV |
| Variants (rows) | 129,273 | 133,871 | 8,171,061 | 529,504 |
| Chromosomes | 1–22, X, Y | 1–22, 23 | 1–22, 23 | 1–22 |
| Assembly | **GRCh37, confirmed**: 50,503/50,503 panel rsIDs at identical GRCh37 positions | **GRCh37, confirmed**: 18,645 chr:pos matches to the panel, 100% allele-compatible. No rsIDs | **GRCh37, confirmed**: 1,147,331/1,147,331 panel rsIDs at identical positions | GRCh37: 100% rsID position match |
| Effect allele | `effect_allele` | `effect_allele` | `A1` (documented as effect allele; lower-case in file) | `A1` |
| Other allele | `other_allele` | `other_allele` | `A2` | `A2` |
| Beta | `beta`: confirmed equal to ln(`odds_ratio`) (max diff 1e-16) | `beta` (log-odds) | Not present; derive as ln(OR) | Not present; ln(OR) used |
| Odds ratio | `odds_ratio`, with 95% CI | Column present but **100% NA** | `OR` | `OR` (quantitative trait coded as OR in the tutorial) |
| SE | `standard_error`: 4 NA, **151 rows = 0** (all with OR = 1, P ≈ 1) | `standard_error` | `SE`: documented as SE of log OR | `SE` |
| Z-score | Not present (derivable as beta/SE) | Not present | Not present | Not present |
| P-value | `p_value`: Wald P reproduced from beta/SE (median \|Δlog10P\| < 1e-15). MHC minimum 9.9e-324 (underflow) | `p_value`: Wald agreement median \|Δlog10P\| = 0.11, from rounding of beta/SE to 4 decimals | `P`: median \|Δlog10P\| = 0.0003 | `P` |
| Allele frequency | `effect_allele_frequency`: **100% NA** | `effect_allele_frequency`: 649 rows at 0 or 1 | `Freq.A1` (1000G EUR, *not* study frequency; 963 rows at 0 or 1) | `MAF` |
| Sample size | No column. Total 47,580 documented (`.meta.yaml`); **cases/controls not available locally** | Per-SNP `Neff` 381–133,611; total 298,957 documented; **cases/controls not available locally** | No column. **Documented: 105,318 = 40,675 cases + 64,643 controls** (`.source.html`) | `N` = 388,028 |
| Trait type | Case-control | Case-control | Case-control | Quantitative (simulated) |
| Ancestry | European (documented) | European (documented) | European (documented) | n/a |
| Duplicate IDs | 962 rows (481 rsIDs), all X/Y pseudoautosomal pairs | No IDs (`variant_id` 100% NA); 4 duplicated chr:pos:allele rows | 1,096 rows; 1,412 duplicated chr:pos rows; 112 duplicated chr:pos:alleles | 4 |
| Non-SNV rows | 345 | 1 | 492,376 indels | 0 |
| Strand-ambiguous SNVs (A/T, C/G) | 15,613 | 15,626 | 1,178,211 | 29,887 |
| p < 5×10⁻⁸ rows; loci (500 kb merge, autosomes) | 2,353; **12** (1 in MHC) | 174; **57** | 18,237; **136** | 9,065; 1,011 |

**Metadata caveat:** both `.meta.yaml` files say `genome_assembly: GRCh38`. That describes the GWAS Catalog *harmonised* `.h.tsv.gz` file, not the local formatted `Build37` file. The GRCh37 label above comes from direct position checks, not from that YAML field.

**Licence:** the SCZ source page restricts use to non-commercial research and prohibits redistribution. `results/part2/*/weights.tsv` and `data/*` are git-ignored, which is consistent with that restriction.

---

## 3. Quarter 1 PRS Pipeline Audit

Code: [`scripts/run_disease_prs.py`](../scripts/run_disease_prs.py) and [`src/prs.py`](../src/prs.py).

| # | Question | Finding |
|---|---|---|
| 1 | GWAS files used | The three `data/gwas_original/` files (`STUDIES`, `run_disease_prs.py:16-20`). The Z-only TCSC files were never scored |
| 2 | Diseases | Rheumatoid arthritis, type 2 diabetes, schizophrenia |
| 3 | Effect-size column | RA: `odds_ratio`. SCZ: `OR`. T2D: `beta` |
| 4 | Scale | RA and SCZ: OR converted to log(OR). T2D: log-odds beta, used as is |
| 5 | OR transformed? | Yes: `np.log(odds_ratio)` (`run_disease_prs.py:52-53`). For RA this equals the file's own `beta` |
| 6 | Z-scores used? | **No.** No code path converts Z to a weight (`grep` of `scripts/` and `src/`; `inventory_part2.py` only labels Z-only files) |
| 7 | Allele alignment | GWAS SNPs are matched to the panel by rsID **and** position. Alleles are accepted if they match directly or after complementing. On a strand flip, the effect allele label is complemented but its beta keeps its sign. Weights are written as `SNP, A1 (GWAS effect allele), BETA`, and PLINK `--score` counts A1 explicitly. **Verified in this audit:** all 229 RA weights and all 250 T2D weights have A1 equal to the original effect allele (or its complement) and the original effect value |
| 8 | Strand-ambiguous SNPs | Dropped (`gwas_qc`, `src/prs.py:17-27`), along with indels and non-ACGT alleles |
| 9 | Genome build | GRCh37 throughout. T2D rows without rsIDs received an rsID only from a *unique* GRCh37 panel position. SCZ coordinates were validated against 1,147,331 panel rsID anchors (100%) |
| 10 | Clumping and threshold | PLINK2 `--clump`, r² 0.1, 250 kb, 1000G EUR LD (489 people). Only SNPs with P ≤ 0.001 were clumped and kept. Sensitivity thresholds (1, 0.05, 0.01, 1e-5) are in `results/part2/exploration/threshold_summary.tsv` |
| 11 | SNPs kept | RA 129,273 → 50,503 harmonized → **229** scored. T2D 133,871 → 18,645 → **250**. SCZ 8,171,061 → 1,147,331 → **3,556** (`results/part2/*/snp_counts.tsv`) |
| 12 | Standardization | Each raw score was z-scored on the 489 1000G EUR reference individuals (mean, sample SD; `run_disease_prs.py:100-101`). The private target was placed relative to a separately rescored common-SNP reference, so the comparison is valid |

**Verdict: no Z-score-as-beta error.** Minor observations; none invalidates the Q1 results, and none was changed:

- `gwas_qc` applies MAF and INFO filters only when columns literally named `MAF`/`INFO` exist. So SCZ's `Freq.A1` and the INFO > 0.9 companion list were not used. RA and T2D had no INFO or usable MAF, so no MAF filter was applied to them either.
- RA's 151 rows with SE = 0 are harmless for PRS, which does not use SE, and all have P ≈ 1. They **must** be dropped for colocalization.
- T2D lost 115,226 of 133,871 rows because they lack rsIDs and their positions are not in the HapMap3-like panel. That is expected for an exome array.
- Clumping LD came from all 489 reference people, which include the 344 eQTL individuals. That is appropriate for PRS.

---

## 4. Genome Assembly and Variant Harmonization

| Resource | Assembly | Evidence |
|---|---|---|
| cis-eQTL table and genotypes (`data/LDREF/1000G.EUR.*`) | GRCh37 | Part 1 verification check "SNP positions match the reference panel"; GENCODE v12 |
| RA, T2D, SCZ GWAS (local files) | GRCh37 | Section 2.2 position checks |
| ENCODE B-cell ATAC peaks `ENCFF421XIL` | GRCh38 → **lifted to GRCh37** | `data/encode/ENCFF421XIL.GRCh37.bed4`. 100,364 of 100,655 peaks mapped (`results/week7/provenance.json`) |

**No further liftover is needed.** Everything is already in GRCh37.

Harmonization rules to use in Quarter 2 (the same rules as Q1, with additions for coloc):

1. Match by chromosome + position + allele pair. rsID is secondary; T2D has none.
2. Drop A/T and C/G SNPs. RA has no allele frequencies, so they cannot be resolved.
3. Drop indels, X/Y/23, duplicated positions, and SE ≤ 0 or missing.
4. Orient every GWAS beta to the eQTL **counted allele**, flipping the sign when the GWAS effect allele is the eQTL's other allele (after complementing where needed).

---

## 5. Disease-by-Disease Feasibility Comparison

| Disease | SNP coverage (eQTL-tested SNPs matched) | Beta available | SE available | Alleles available | N available | PRS usable | Coloc usable | Recommendation |
|---|---|---|---|---|---|---|---|---|
| Rheumatoid arthritis | 46,440 (4.6%); dense at Immunochip loci: 8/12 loci ≥50 SNPs in ±250 kb | Yes (`beta` = ln OR) | Yes (155 NA/0 to drop) | Yes; no frequency | Total only (47,580); cases/controls not local | Yes (used) | **After preprocessing**, `coloc.abf` at dense loci | **Primary** |
| Schizophrenia | 963,883 (96.2%); 117/136 loci ≥50 SNPs | Derive ln(OR) | Yes | Yes + 1000G frequency | Yes: 40,675 / 64,643 | Yes (used) | Technically yes, after preprocessing | Secondary contrast only |
| Type 2 diabetes | 17,363 (1.7%); 1/57 loci ≥50 SNPs (median 6) | Yes | Yes | Yes + frequency | Per-SNP Neff; cases/controls not local | Yes (used) | **No**: too sparse for regional coloc | Exclude from coloc |
| Height (simulated) | 121,926 (12.2%) | ln(OR) | Yes | Yes | Yes | Teaching only | Not meaningful | Exclude |

**Descriptive ATAC overlap** (not LD-aware; correlated SNPs are counted repeatedly, so this is *not* an enrichment test):

| Disease | Significant GWAS rows in B-cell peaks | All rows in peaks |
|---|---|---|
| RA | 5.4% | 2.9% |
| SCZ | 3.6% | 1.2% |
| T2D | 2.9% | 3.3% |

**Biological relevance:**

- **RA** is an autoimmune disease with an established B-cell component (autoantibodies; anti-CD20 therapy is used clinically). Its significant loci include *PTPN22* (rs2476601, P = 9×10⁻⁶²), *STAT4*, *TNFAIP3*, *IL6R*, *TYK2*, *PADI4*, *AFF3*, *RBPJ*, *ANKRD55* (nearest-gene labels, by position) and MHC.
- **SCZ** is a brain disorder. LCL and B-cell regulatory data are a weak tissue match, so any shared signal would need cautious interpretation.
- **T2D** is driven by islets, liver, adipose and muscle. Its relevance to B cells is low.

---

## 6. GWAS–cis-eQTL Overlap Feasibility

**eQTL table** ([`data/quarter1_results/part1/table_A_genome_wide_cis_eqtl.tsv.gz`](../data/quarter1_results/part1/table_A_genome_wide_cis_eqtl.tsv.gz)):

- 6,597,508 SNP–gene pairs; 1,001,779 unique SNPs; 14,661 genes
- N = 344 in every row; no missing values; no duplicate pairs
- 101,916 pairs significant at BH < 0.05; 20,172 at Bonferroni < 0.05

**Overlap at significant GWAS loci:**

| | RA | SCZ | T2D |
|---|---|---|---|
| Significant GWAS SNPs that are eQTL-tested | 319 | 2,305 | 103 |
| … that are also BH-significant eQTL SNPs | 136 | 900 | 14 |
| Loci where a significant GWAS SNP is itself a BH eQTL SNP | 2 of 12 (*PTPN22* region, MHC) | 29 of 136 | 7 of 57 |

These counts show that the two kinds of signal *co-occur*. **They do not show a shared causal variant.** In LD-dense regions such as the MHC, co-occurrence is expected by chance.

**Panel limitation:** the eQTL SNPs are the ~1.0M HapMap3-like SNPs in `data/LDREF` and are not imputed. A true causal variant may be absent, in which case colocalization relies on tag SNPs.

---

## 7. Colocalization Data Requirements

Source: the [coloc data article](https://chr1swallace.github.io/coloc/articles/a02_data.html) and `check_dataset`.

- **Each dataset:** one trait in one region, with `type` always given. The preferred inputs are `beta` and `varbeta` (SE², not SE).
- **`type = "cc"`, beta/varbeta form:** the official minimal example has no `s` (case fraction); `s` is required only for the p-value/MAF form.
- **`type = "quant"`:** needs `sdY`, or else `N` + `MAF` so that coloc can estimate it.
- **SNP matching:** SNPs are matched by ID, and only shared SNPs contribute. The region must be **densely covered and not p-value thresholded**, and the alleles must be aligned.
- **Model assumption:** `coloc.abf` assumes at most one causal variant per trait in the region. `coloc.susie` relaxes this, but needs an LD matrix aligned to the same alleles.

| Requirement | RA GWAS | eQTL (Part 1) | Status |
|---|---|---|---|
| SNP ID / position | rsID + GRCh37 position | `SNP_ID`, `SNP_Chr`, `SNP_Position` | ✅ |
| Effect / other allele | `effect_allele` / `other_allele` | `Counted_Allele` is in table A. **The other allele is not in table A** but is the `ALT` column of `chromosomes/chr*.traw` (= `.bim` column 5). Verified on chr22: the counted allele is always `.bim` column 6 and is usually the major allele (median counted frequency 0.78) | ✅ after a join |
| beta, varbeta | `beta`, `standard_error`² (drop 155 rows) | `Beta`, `SE`² | ✅ after preprocessing |
| P | `p_value` | `P_Value` | ✅ |
| Type | cc | quant | ✅ |
| N | Total 47,580 (documented) | 344 | ✅ |
| `s` (cases/total) | **Not available locally** | n/a | Not needed for the beta/varbeta form; obtain it for completeness and for reporting |
| MAF | **Not available** (column all NA) | Not in table A; **recoverable** from the `.traw` dosages or `plink2 --freq --keep matched.keep` | eQTL ✅ after a cheap computation; GWAS MAF not needed for coloc.abf with beta/varbeta |
| sdY (eQTL) | n/a | Expression is PEER-normalized and **not standardized**. Compute each gene's SD in the 344 people from `GD462…resk10.txt`, or let coloc estimate it from MAF + N | ✅ after preprocessing |
| LD (coloc.susie only) | 1000G EUR reference; out-of-sample for the RA cohorts | In-sample: the 344 people are in `data/LDREF` | ⚠️ Possible; LD mismatch with the multi-cohort RA meta-analysis is a risk |
| Software | — | — | ❌ **R, `coloc` and `susieR` are not installed** (`Rscript` not found) |

**Classification:**

| Dataset | Classification |
|---|---|
| RA | **Supported after documented preprocessing** (`coloc.abf`), at the ~8 dense Immunochip loci. Not supported genome-wide |
| SCZ | Supported after preprocessing. Biologically weak for this tissue |
| T2D | **Not supported** without different summary statistics (a genome-wide imputed T2D GWAS) |
| Height | Not applicable (simulated) |

**What coloc does not establish:** a high PP.H4 indicates that one shared variant is *probable* under the model's priors. It is not proof of causality, and sharing significant SNPs alone is not colocalization.

---

## 8. Recommended Primary Disease

**Primary: rheumatoid arthritis**, file `data/gwas_original/GCST005569.Build37.f.tsv.gz`.

1. **Biological relevance.** RA is an autoimmune disease with B-cell involvement, and its loci include genes active in lymphocytes (*PTPN22*, *STAT4*, *TNFAIP3*, *IL6R*, *TYK2*). LCL eQTLs and primary B-cell ATAC-seq are a plausible regulatory context.
2. **Dataset.** Eyre et al. 2012 Immunochip meta-analysis, GWAS Catalog formatted GRCh37 file.
3. **Completeness.** Beta, SE, both alleles, P and GRCh37 position are all present. Allele frequency and per-SNP N are missing. Coverage is targeted, not genome-wide.
4. **Match to the eQTL data.** 46,440 SNPs are shared, 100% position- and allele-compatible. 8 of 12 significant loci have ≥50 shared SNPs within ±250 kb.
5. **ENCODE context.** It is appropriate. The peaks come from a primary B cell (CL:0000236), adult male donor, pseudoreplicated peaks lifted to GRCh37. LCLs are EBV-transformed B cells, so the cell types are related but not identical.
6. **Formal coloc.** Possible with `coloc.abf` at the dense loci once R/coloc is installed. Exclude the MHC locus (chr6:28–34 Mb), where LD and HLA effects make single-variant colocalization unreliable.
7. **Additional data needed.**
   - Required: R + `coloc` (and optionally `susieR`).
   - Recommended: RA case/control counts from Eyre et al. 2012 or the GWAS Catalog record (to report `s`), and RA allele frequencies (to rescue A/T and C/G SNPs). Both are optional for `coloc.abf`.
8. **Limitations.**
   - Only 12 significant loci (11 outside the MHC), so the gene ranking will be small.
   - Immunochip density varies by locus.
   - The eQTL sample is small (N = 344) and not imputed.
   - LCLs are not primary B cells.
   - ENCODE accessibility comes from one donor.
   - Ancestry is European only.

**Optional secondary: schizophrenia** (`GCST006803.txt.gz`), as a *tissue-mismatch contrast*. SCZ data are dense and complete, and a lower B-cell accessibility signal than RA would be informative. Do not frame SCZ hits in LCLs as mechanistic.

**Not recommended: T2D.** The exome array is too sparse for regional analysis, and the trait is not B-cell related.

---

## 9. Proposed Quarter 2 Analysis

**Recommendation: Option A, scoped narrowly, with Option B as a pre-declared fallback.**

**Question A (main quantitative result).** Test whether significant cis-eQTL SNPs fall in B-cell ATAC peaks more often than matched nonsignificant SNPs.

- Matching: MAF, distance to TSS, and LD score or number of LD partners. Week 7 found 4.8% of Bonferroni eQTL SNPs in peaks vs 1.6% of all tested SNPs (`results/week7/snp_overlap_summary.tsv`), but that comparison is not yet matched.
- Use resampling or logistic regression.

**Question B (RA).**

1. Define regions around RA lead SNPs at p < 5×10⁻⁸, excluding the MHC.
2. For each gene whose cis window overlaps a region, run `coloc.abf` on the shared, allele-aligned SNPs:
   - RA: `type = "cc"`, beta, varbeta
   - eQTL: `type = "quant"`, beta, varbeta, N = 344, plus `sdY` or MAF
3. Report PP.H3 and PP.H4, with prior sensitivity (`coloc::sensitivity`).
4. Optional: `coloc.susie` using in-sample LD for the eQTL side.
5. Annotate credible-set SNPs with B-cell ATAC peaks.
6. Rank genes by PP.H4, eQTL strength and accessibility.

**Fall back to Option B if any of these happen:**

- R/coloc cannot be installed.
- Allele alignment checks fail.
- Fewer than ~5 loci have adequate shared density.

In that case, report descriptive locus overlap and LD-aware interpretation, and state explicitly that **overlap does not establish a shared causal variant**.

**Why not Option C?** The data for formal coloc at RA loci exist, so downgrading GWAS integration is unnecessary. However, because RA offers only about 11 usable loci, Question A should still carry the main statistical weight.

---

## 10. Missing Data and Remaining Questions

| Item | Status | How to resolve |
|---|---|---|
| R, `coloc`, `susieR` | Not installed | Install (needs your approval; this is a download) |
| RA cases/controls (`s`) | Not in local files; only total N = 47,580 | Eyre et al. 2012 / GWAS Catalog GCST005569 record |
| RA allele frequencies | Column all NA | Use 1000G EUR frequencies from `data/LDREF` as a proxy, or drop A/T and C/G SNPs |
| eQTL other allele and MAF | Not in table A | Join from `.traw`/`.bim`; `plink2 --freq` on 344 people (cheap) |
| eQTL `sdY` | Expression not standardized | Per-gene SD from the expression file, or estimate from MAF + N |
| Exact RA sample definition (47,580 vs the paper's discovery set) | Not verified | Check the paper |
| T2D cases/controls | Not local (per-SNP Neff only) | Only needed if T2D is used, which is not recommended |
| Causal variants missing from the HapMap3-like panel | Inherent | Acknowledge as a limitation; imputation is out of scope |

---

## 11. Recommended Next Steps

1. Approve installing R + `coloc` (+ `susieR`) into the `capstone` environment.
2. Look up the RA case/control counts and the exact sample composition in Eyre et al. 2012.
3. Write `scripts/prepare_coloc_inputs.py` (new outputs only):
   - Take the RA regions from `results/q2_planning/gwas_lead_loci_500kb.tsv`, excluding the MHC.
   - Join eQTL rows with the `.traw` ALT allele and computed MAF.
   - Align alleles to the counted allele; drop ambiguous SNPs and SE ≤ 0.
   - Save one input table per locus and gene.
4. Run `coloc.abf` and `check_alignment` on 1–2 pilot loci (e.g. *PTPN22*, *TNFAIP3*) before scaling up.
5. In parallel, finish Question A's matched-background enrichment, which is independent of the GWAS work.
6. Update the proposal:
   - Primary disease: RA (Immunochip GRCh37 file).
   - Secondary contrast: SCZ (optional).
   - State the Option A scope and the Option B fallback criteria.
