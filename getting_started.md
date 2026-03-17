# Getting Started: Visualizing Results

There are three ways to explore results depending on how far through the pipeline you are.

---

## Option A — Interactive dashboard (recommended)

The Streamlit app at `app.py` is the primary interface. It works with partial data and shows info banners for anything not yet generated.

```bash
streamlit run app.py
```

Opens at `http://localhost:8501`. Use the **sidebar** to filter tumor types and navigate between sections:

| Tab | What you see |
|---|---|
| **Overview** | Sample counts, hot tumor fraction by cancer type, age distributions |
| **Immune Classification** | IFN-γ score histogram, CD8+ fraction boxplots, IFN-γ vs CD8 scatter |
| **Senescence Scores** | Per-score violin plots (hot vs cold), alteration rates by gene, CDKN2A by tumor type, Thorsson immune subtype cross-tabs |
| **Correlations** | Interactive scatter (any senescence score vs IFN-γ / CD8 fraction), Spearman ρ, per-type log-rank p-values |
| **Survival** | Cox forest plot (HR ± 95% CI), per-tumor-type log-rank table, interactive Kaplan-Meier explorer (endpoint × stratification × follow-up slider) |
| **Pathway** | Interactive p16/CDK4-6/Rb/E2F pathway diagram with SASP outputs and cancer escape routes |
| **Figures** | All static PNGs from the pipeline displayed as a gallery |

The KM explorer in the Survival tab lets you switch between OS / DSS / PFI endpoints and stratify by hot/cold, p16 status, or the 4-group combination — no rerunning notebooks needed.

---

## Option B — Static figures

All plots are written to `figures/` during notebook execution. Key outputs:

| File | Source | Content |
|---|---|---|
| `ifng_score_distribution.png` | notebook 02 | IFN-γ score histogram with hot/cold split |
| `hot_fraction_by_tumor_type.png` | notebook 02 | Bar chart of hot fraction per cancer type |
| `cdkn2a_by_tumor_type.png` | notebook 03 | CDKN2A expression boxplots across 33 types |
| `cdkn2a_vs_ifng_scatter.png` | notebook 04 | Pan-cancer CDKN2A vs IFN-γ scatter with Spearman ρ |
| `cdkn2a_ifng_corr_by_type.png` | notebook 04 | Per-tumor-type Spearman ρ heatmap |
| `spearman_heatmap_pancancer.png` | notebook 04 | Multi-score × IFN-γ correlation heatmap |
| `sasp_score_hot_cold.png` | notebook 04 | SASP NES violin: hot vs cold (H2) |
| `km_hot_cold.png` | notebook 05 | KM curves: hot vs cold OS |
| `km_4group.png` | notebook 05 | KM curves: hot/cold × p16-high/low 4-group (H4) |
| `cox_forest_plot.png` | notebook 05 | Cox HR forest plot with 95% CI (H5) |
| `fig1_senescence_heatmap.png` | notebook 06 | Composite: senescence score heatmap |
| `fig2_umap.png` | notebook 06 | UMAP of samples colored by immune/senescence features |
| `fig3_correlation_matrix.png` | notebook 06 | Correlation matrix across all scores |
| `fig4_survival_composite.png` | notebook 06 | Four-panel survival composite figure |

Open any PNG directly or view them all at once in the **Figures** tab of the dashboard.

---

## Option C — R Markdown HTML reports

Notebooks 04 and 05 knit to self-contained HTML files that include all statistical output (p-values, BH FDR results, `cox.zph()` PH checks):

```
notebooks/04_overlap_analysis.html   # H1–H3 statistical tests
notebooks/05_survival.html           # H4–H5, KM plots, Cox tables
```

Open these in any browser — no server needed. To regenerate them:

```r
rmarkdown::render("notebooks/04_overlap_analysis.Rmd")
rmarkdown::render("notebooks/05_survival.Rmd")
```

Or from the terminal:

```bash
Rscript -e 'rmarkdown::render("notebooks/04_overlap_analysis.Rmd")'
Rscript -e 'rmarkdown::render("notebooks/05_survival.Rmd")'
```

---

## Prerequisites

If you haven't run the pipeline yet, you need the intermediate parquet files before the dashboard or figures will populate. See `start-here.md` for the full step-by-step setup. Minimum steps to get the dashboard working:

1. Run notebooks 01–03 (Python) to produce `data/sample_immune_labels.parquet` and `data/sample_senescence.parquet`
2. Knit notebooks 04–05 (R) to produce `data/cox_results.csv` and `data/logrank_per_type.csv`
3. Optionally run notebook 06 (Python) to generate the UMAP and composite figures

The dashboard runs fine with partial data — any tab missing its source files will display an explanatory banner rather than crashing.
