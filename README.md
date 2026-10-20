# Hot/Cold Tumor × p16/CDKN2A — TCGA Pan-Cancer Analysis

Bioinformatics analysis exploring whether cellular senescence status (p16/CDKN2A) correlates with immune phenotype and patient survival across TCGA pan-cancer RNA-seq data.

**Central question:** Do p16-high tumors have more immune-inflamed microenvironments, and does that translate to better survival?

**GitHub Pages:** [kuod.github.io/p16-hot-cold](https://kuod.github.io/p16-hot-cold/)

## Background

- **Hot vs cold tumors** — immune-inflamed vs immune-excluded phenotypes, classified by a 6-gene IFN-γ signature (adapted from Ayers et al. 2017: IFNG, STAT1, CCL5, CXCL9, CXCL10, IDO1; CCL5 replaces HLA-DRA); samples are split into pan-cancer IFN-γ **tertiles** (cold / intermediate / hot)
- **p16/CDKN2A** — CDK4/6 inhibitor and marker of cellular senescence; senescent cells secrete pro-inflammatory SASP cytokines (IL-6, CXCL8, CCL2) that may recruit T cells
- **Hypothesis:** p16 → SASP → T cell recruitment → hot tumor phenotype → better survival

## Quick Start

```bash
pip install -r requirements.txt
streamlit run app.py
```

Run the unit tests with `pytest tests`.

The dashboard runs with partial data — tabs show info banners for files not yet generated.

To run the full pipeline, see [start-here.md](start-here.md).

## Pipeline

| Step | Notebook | Language | Runtime |
|---|---|---|---|
| 1 | `01_data_download.ipynb` | Python | ~30 min |
| 2 | `02_immune_classification.ipynb` | Python | ~5 min |
| 3 | `03_p16_senescence.ipynb` | Python | ~20–30 min (ssGSEA) |
| 4 | `04_overlap_analysis.Rmd` | R | ~5 min |
| 5 | `05_survival.Rmd` | R | ~5 min |
| 6 | `06_figures.ipynb` | Python | ~10 min |

R notebooks knit via:
```bash
Rscript -e 'rmarkdown::render("notebooks/04_overlap_analysis.Rmd")'
Rscript -e 'rmarkdown::render("notebooks/05_survival.Rmd")'
```

## Pre-specified Hypotheses

| # | Hypothesis | Test |
|---|---|---|
| H1 | Pan-cancer ρ(CDKN2A_expr, IFN-γ score) > 0 | Spearman (single pan-cancer `cor.test`; BH-FDR only on exploratory tables) |
| H2 | SASP score higher in hot vs cold tumors | Three-group Kruskal-Wallis (tertiles) on SASP score, plus group medians |
| H3 | CDKN2A deletion enriched in cold tumors | Chi-squared (pooled) + Mantel-Haenszel (stratified by tumor type) |
| H4 | p16-high/hot tumors have better OS than p16-low/cold | Log-rank (6-group KM: tertile × p16; unstratified and stratified by tumor type) |
| H5 | CDKN2A_expr is an independent predictor of OS | Cox PH (HR, 95% CI) |

Bonferroni α = 0.01 across H1–H5. All other tests are exploratory (BH FDR α = 0.05).

## Data Sources

All data is public and downloaded at runtime. Tumor purity (ABSOLUTE) comes from the NCI GDC and is included as a covariate in the Cox models (notebook 05; the SASP and IFN-γ purity sensitivity analyses are in notebook 04):

- **UCSC Xena** — TCGA pan-cancer RNA-seq (RSEM log2(norm_count+1), EB++ AdjustPANCAN)
- **cBioPortal** — CDKN2A/TP53/RB1/ATM/PTEN CNA + mutation status (TCGA PanCancer Atlas 2018)
- **Thorsson et al. 2018** — CIBERSORT immune fractions and immune subtype (C1–C6)
- **Taylor et al. 2018** — TCGA PanCanAtlas ABSOLUTE tumor purity/ploidy calls (NCI GDC); about 95.8% of the 9,702 retained samples have a purity value

## Architecture

```
src/signatures.py    # Curated gene lists (IFN-γ, SASP, senescence effectors, MSigDB set names)
src/utils.py         # Shared utilities: log2tpm, signature_score, classify_hot_cold
notebooks/           # Analysis pipeline (01–03, 06 Python; 04–05 R Markdown)
app.py               # Streamlit dashboard
r_packages.R         # One-time R package install
data/                # Gitignored intermediate data (parquet files bridge Python → R)
figures/             # Gitignored output plots
```

## Dashboard

```bash
streamlit run app.py
```

| Tab | Contents |
|---|---|
| Overview | Sample counts, hot tumor fraction by cancer type, age distributions |
| Immune Classification | IFN-γ score histogram, CD8+ fraction, IFN-γ vs CD8 scatter |
| Senescence Scores | Per-score violin plots, alteration rates, Thorsson immune subtype cross-tabs |
| Correlations | Interactive scatter (any score vs IFN-γ / CD8 fraction), per-type ρ |
| Survival | Cox forest plot, per-tumor-type log-rank table, interactive Kaplan-Meier explorer |
| Pathway | Interactive p16/CDK4-6/Rb/E2F pathway diagram |
| Figures | Static PNG gallery |
| Hypotheses | Pre-specified hypothesis tests with live results |
| Summary | Narrative synthesis of all findings |

## Disk Space

Run `clean.sh` to remove regenerable files (~370 MB):

```bash
./clean.sh           # dry run
./clean.sh --force   # delete
```

## License

See [LICENSE](LICENSE).
