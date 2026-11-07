# Privacy and Reproducibility Checklist

Run this before every push:

```bash
python scripts/check_public_files.py
python -m unittest discover -s tests
git status --short
```

## Privacy

- [x] Git history reviewed: the only earlier commit contains `.gitignore` and `README.md`.
- [x] `.gitignore` excludes personal genotype files everywhere in the repository, not only in `data/`.
- [x] The local weekly notebook is ignored. `notebooks/weekly_tasks_public.ipynb` is the version to publish.
- [x] Public Week 5: the private section is left out and replaced by the same workflow on a public participant.
- [x] The scanner checks every file Git would publish, including personal values read from local results at run time.
- [x] Demo genotypes contain only people from the public 1000 Genomes panel.
- [x] Owner reviewed the items below that needed a decision.

## Must Stay Local

- The raw 23andMe file and anything in `data/week5_results/`, including converted VCF, BCF, and PLINK files, call tables,
  scores, summaries, logs, and figures.
- `quater1_weekly_tasks.ipynb`, whose Week 5 outputs contain personal results, and its backups in `data/notebook_backups/`.

## Excluded Because of Size or Redistribution

- `data/LDREF/`, the full expression file, `data/reference/`, `data/sumstats/`, `data/Height.gwas.txt.gz`, and other `data/` folders.
- `plink2_binary` and `plink2/`: install PLINK2 as the README describes.

## Owner Decisions

Decided by the repository owner on 2026-11-06:

- **Demo expression values:** publishing the 12 genes of real GEUVADIS data in `data/demo/` is accepted.
- **The name "Tiffany" in the public notebook:** stays in the Week 5 text only. It contains no genotypes or results.
- **No `LICENSE` file:** none will be added. Without one, others have no formal permission to reuse the code.

## Reproducibility

- [x] `environment.yml` lists the versions from the working environment, and a fresh environment was built from it.
- [x] PLINK2 is configurable, with no personal absolute paths.
- [x] From a copy of only the publishable files, in a fresh environment: the demo ran, the tests passed, and Part 1
  reproduced its tables byte for byte.
- [x] Moving code into `src/` left all Part 1 result tables byte-identical.
- [x] Random seeds are fixed, and the demo's outputs are identical across runs.
- [ ] Linux and Windows are untested.
- [ ] The expression file has no confirmed direct download link; it was supplied through the course.
