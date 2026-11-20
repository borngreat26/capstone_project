# Quarter 1 paper figures

Public reference-only collection: personal target figures and results are excluded.

Give this entire folder to the agent writing the LaTeX paper. It contains 25
distinct figures (43 PNG/PDF/SVG files) from Part 1, Part 2 and Week 7.
The analysis notebooks can refresh it without rerunning any regressions or scores.

- `figure_inventory.tsv`: filenames, original locations, data sources, suggested captions,
  paper sections and interpretation limitations. `part1_`, `part2_` and `week7_` identify
  the analysis; the same stem in PNG and PDF is one figure in two formats.
- `figures.tex`: ready-to-adapt figure environments and cross-reference examples.
  Upload this directory as `paper_figures_public/` alongside the paper's main `.tex` file.
  Select the relevant examples rather than including every figure automatically.
- `part1_methods_results.md` and `part2_week7_methods_results.md`: public writing notes
  copied from `docs/`; relative links in these notes refer to the original repository.

Use PDFs where provided. New analysis PNGs are 320 dpi; preserved Part 1 PNGs retain
their original resolution. No image was resized or regenerated for this collection.
The SNPs-tested-per-gene histogram and exploratory sensitivity figures are supplementary
candidates. Captions and notes explain which results are associations or relative scores.

The disease PRSs use published disease GWAS effects, not our expression eQTL betas.
Raw PRSs, z-scores and percentiles are not calibrated disease probabilities. Primary
B-cell ATAC overlap does not prove a causal eQTL or a regulatory link to a gene.
Raw genotypes and original GWAS effects are excluded in both modes. Personal target figures
and a local results note are included only in the explicitly requested local collection.

Refresh after changing the source figures:

```bash
python scripts/make_paper_support.py


```

To refresh copies using the existing inventory only:

```bash
python scripts/export_paper_figures.py
```
