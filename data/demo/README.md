# Demo Data

A small dataset that lets anyone run this project's workflow without private data or the full downloads.
**It only demonstrates the computation. Its scores are not disease-risk estimates, and they are not the
results of the private Week 5 analysis.** `scripts/build_demo_data.py` created every file here (random seed 42).

| File | Contents | Source |
|---|---|---|
| `demo_individual.{bed,bim,fam}` | One person, chosen at random from the panel | Real public genotypes |
| `demo_reference.{bed,bim,fam}` | 244 reference people, not including the demo individual | Real public genotypes |
| `demo_gwas.tsv` | Summary statistics for 1,548 SNPs, plus test rows | **Simulated** (see below) |
| `demo_reference_simulated_phenotype.tsv` | Simulated trait for the reference people | **Simulated** |
| `demo_eqtl.{bed,bim,fam}` | The 344 people who also have expression data | Real public genotypes |
| `demo_expression.tsv` | Expression for 12 genes, rounded to 6 decimals | Real public data |
| `demo_genes.tsv` | The 12 protein-coding genes and their cis-windows | GENCODE v12 |

**Region:** all SNPs are on chromosome 22, 44.5–47.0 Mb (GRCh37), with their real alleles and positions.

**Sources:** genotypes are a subset of the [1000 Genomes EUR LD reference](https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2).
Expression comes from the course's GEUVADIS file (ArrayExpress [E-GEUV-1](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-GEUV-1)).
Genes come from [GENCODE v12](https://www.gencodegenes.org/human/release_12.html).

## The Simulated GWAS

This GWAS is simulated, not published, so its weights say nothing about any real trait.

1. The panel's other 244 people form a discovery group, separate from the reference group and the demo
   individual. Their genotypes are not in this folder.
2. A trait is simulated from their real genotypes, with 15 causal SNPs and 50% heritability.
3. Each SNP gets a real linear regression on that trait, which gives `BETA`, `SE`, and `P`.

- **`A1`** is the effect allele.
- **`BETA`** is the change in the simulated trait per copy of `A1`. It is not an odds ratio.
- **`N`** is 244. **`MAF`** is measured in the discovery group.

## Test Rows

`Demo_Case` labels rows changed on purpose so that harmonization has something to catch. Real GWAS files have no such column.

| `Demo_Case` | Rows | What should happen |
|---|---:|---|
| `Opposite strand` | 6 | Relabeled with complementary alleles, and kept |
| `Wrong position` | 2 | Removed |
| `Incompatible alleles` | 3 | Removed |
| `Not in panel` | 3 | Removed |
| `Duplicate row` | 2 | Only the first copy is kept |
