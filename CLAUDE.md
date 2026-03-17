# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Purpose

Bioinformatics analysis exploring the intersection of:
- **Hot vs cold tumors** (immune-inflamed vs immune-excluded) using TCGA pan-cancer RNA-seq
- **p16/CDKN2A** as a marker of cellular senescence and SASP (senescence-associated secretory phenotype)
- Whether senescence status in the tumor microenvironment correlates with immune phenotype and patient survival

## Setup

**Python** (data acquisition + ssGSEA scoring):
```bash
pip install -r requirements.txt
jupyter lab
```

**R** (statistical testing + survival analysis):
```bash
Rscript r_packages.R   # one-time package install
```

## Notebook Execution Order

Run notebooks sequentially — each produces parquet files consumed by the next:

| Notebook | Language | Output | Runtime |
|---|---|---|---|
| `01_data_download.ipynb` | Python | Raw data in `data/` | ~30 min (network) |
| `02_immune_classification.ipynb` | Python | `data/sample_immune_labels.parquet` | ~5 min |
| `03_p16_senescence.ipynb` | Python | `data/sample_senescence.parquet` | ~20–30 min (ssGSEA) |
| `04_overlap_analysis.Rmd` | R | Figures in `figures/` | ~5 min |
| `05_survival.Rmd` | R | Survival figures + `data/cox_results.csv` | ~5 min |
| `06_figures.ipynb` | Python | Composite figures | ~10 min (UMAP) |

Knit the `.Rmd` notebooks with `rmarkdown::render("notebooks/04_overlap_analysis.Rmd")` or via RStudio.

**Note:** The ssGSEA step in notebook 03 caches its output to `data/ssgsea_senescence_scores.parquet` — delete this file to recompute.

## Architecture

```
src/signatures.py        # All curated gene lists (IFN-γ, SASP, senescence effectors, MSigDB set names)
src/utils.py             # Shared functions: log2tpm, signature_score, classify_hot_cold
notebooks/               # Analysis notebooks (01–06; 01–03, 06 Python; 04–05 R Markdown)
r_packages.R             # One-time R package install script
data/                    # Gitignored raw + intermediate data (parquet files bridge Python → R)
figures/                 # Gitignored output plots
```

**Key design decisions:**
- Hot/cold classification uses **IFN-γ 6-gene signature** (median split); CD8+ T cell fraction from Thorsson et al. 2018 CIBERSORT is used as secondary validation
- CDKN2A genomic status (deleted/mutated) is pulled from cBioPortal TCGA PanCancer Atlas 2018 study (`tcga_pan_can_atlas_2018`)
- ssGSEA senescence scores use MSigDB sets: `FRIDMAN_SENESCENCE_UP`, `REACTOME_CELLULAR_SENESCENCE`, plus custom SASP and effector lists defined in `src/signatures.py`
- **Python handles data acquisition and ssGSEA scoring (notebooks 01–03); R handles all statistical testing and survival analysis (notebooks 04–05).** Data passes between layers via parquet files read with `arrow`.
- Survival analysis classifies p16 high/low **within each tumor type** (not globally) to account for baseline CDKN2A expression differences across cancer types
- Cox model uses `strata(tumor_type)` to allow per-type baseline hazards without one-hot encoding all tumor types

## Data Sources

All data is publicly available and downloaded at runtime:
- **UCSC Xena** — TCGA pan-cancer RNA-seq (Xena hub S3 bucket)
- **cBioPortal** — REST API at `cbioportal.org/api`

## Pre-specified Primary Hypotheses

The following are the **confirmatory** hypotheses for this analysis. All other tests are exploratory.

| # | Hypothesis | Test | Notebook |
|---|---|---|---|
| H1 | Pan-cancer Spearman ρ(CDKN2A_expr, IFN-γ score) > 0 | Spearman, BH-FDR | 04 §2 |
| H2 | SASP score is higher in hot vs cold tumors (pan-cancer) | Mann-Whitney U | 04 §5 |
| H3 | CDKN2A deletion enriched in cold tumors (pan-cancer) | Chi-squared | 04 §4 |
| H4 | p16-high hot tumors have better OS than p16-low cold tumors | Log-rank (4-group KM) | 05 §4 |
| H5 | CDKN2A_expr is an independent predictor of OS (age-adjusted, strata tumor type) | Cox PH (HR, 95% CI) | 05 §6 |

**Multiple testing:** BH FDR is applied within each analysis. Across H1–H5, a Bonferroni threshold of α = 0.01 is used for the five primary tests. All other p-values reported are exploratory and should be interpreted accordingly.

**Expected effect sizes (for sanity-checking results):**
- H1: ρ ~ 0.10–0.25 (weak-moderate positive; p16 promotes SASP → immune infiltration)
- H2: SASP NES higher in hot, p < 0.01
- H3: CDKN2A deletion rate higher in cold tumors, p < 0.01
- H4: Hot/p16-high best OS, cold/p16-low worst, global log-rank p < 0.01
- H5: CDKN2A_expr HR < 1.0 (higher expression → lower hazard, p < 0.05)

## Tumor Purity Notes

TCGA bulk RNA-seq mixes tumor cells with stroma and infiltrating immune cells. Purity estimates (TCGA ABSOLUTE, from Aran et al. 2015) are downloaded in notebook 01 and used as a covariate in Cox models (notebook 05). Spearman correlations involving immune scores should be interpreted with awareness that high-stroma samples inflate immune signals.

## New Data Sources (added)

- **TCGA ABSOLUTE purity** — UCSC Xena (Aran et al. 2015 consensus estimates)
- **Thorsson immune subtypes (C1–C6)** — UCSC Xena (`Subtype_Immune_Model_Based.txt.gz`)
- **TMB (nonsynonymous mutations/Mb)** — cBioPortal per-study clinical attributes
- **DSS / PFI survival endpoints** — already in `tcga_clinical_survival.tsv` (unused until notebook 05 update)
