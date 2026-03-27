# Start Here

Complete these steps in order. Each step depends on the previous one.

---

## Step 1 — Install dependencies

```bash
# Python
pip install -r requirements.txt

# R (run once)
Rscript r_packages.R
```

---

## Step 2 — Download raw data (`01_data_download.ipynb`)

Open Jupyter and run notebook 01:

```bash
jupyter lab
```

Everything downloads automatically — UCSC Xena RNA-seq, clinical survival table,
Thorsson et al. 2018 CIBERSORT immune fractions (NCI GDC), and cBioPortal CNA/mutation
data. No manual steps required. Expect ~30 minutes on a typical connection.

---

## Step 3 — Classify hot/cold tumors (`02_immune_classification.ipynb`)

Run notebook 02. Produces:

- `data/sample_immune_labels.parquet`
- `figures/hot_fraction_by_tumor_type.png`
- `figures/ifng_score_distribution.png`

~5 minutes.

---

## Step 4 — Score p16/senescence (`03_p16_senescence.ipynb`)

Run notebook 03. The ssGSEA step is the slow part (~20–30 min) and caches to
`data/ssgsea_senescence_scores.parquet` — subsequent runs are fast. Produces:

- `data/sample_senescence.parquet`
- `figures/cdkn2a_by_tumor_type.png`

---

## Step 5 — Statistical analysis in R

Knit both R Markdown notebooks. In R or RStudio:

```r
rmarkdown::render("notebooks/04_overlap_analysis.Rmd")
rmarkdown::render("notebooks/05_survival.Rmd")
```

Or from the terminal:

```bash
Rscript -e 'rmarkdown::render("notebooks/04_overlap_analysis.Rmd")'
Rscript -e 'rmarkdown::render("notebooks/05_survival.Rmd")'
```

Notebook 04 produces figures in `figures/`. Notebook 05 produces figures plus
`data/cox_results.csv` and `data/logrank_per_type.csv`.

---

## Step 6 — Summary figures (`06_figures.ipynb`)

Run notebook 06 in Jupyter. Produces composite figures including UMAP (~10 min).

---

## Step 7 — Launch the dashboard

```bash
streamlit run app.py
```

Opens at `http://localhost:8501`. All 6 tabs are populated once Steps 2–6 are
complete. The app runs fine with partial data — tabs show info banners for any
files not yet generated.

---

## Key hypotheses being tested

| Hypothesis | Test | Notebook |
|---|---|---|
| SASP-high → hot (SASP drives TIL infiltration) | Mann-Whitney U, Spearman ρ | 04 |
| CDKN2A deletion → cold (loss of SASP signaling) | Chi-squared contingency | 04 |
| Relationship is cancer-type specific | Per-type Spearman, BH FDR | 04 |
| Hot tumors → better survival | KM + log-rank | 05 |
| p16-high/hot combination → best prognosis | 4-group KM + log-rank | 05 |
| IFN-γ and CDKN2A independently predict survival | Cox PH with strata | 05 |

---

## Freeing disk space

Run `clean.sh` to remove all regenerable files (~370 MB). A dry-run is shown by default:

```bash
./clean.sh           # preview what would be deleted
./clean.sh --force   # actually delete
```

Everything removed is recreated by re-running the pipeline (Steps 2–6 above).
The largest item is the raw TCGA expression matrix (~316 MB, ~30 min to redownload).

---

## Verification checklist

Before trusting results, confirm:

- [ ] `nrow` of both parquet files matches (inner join should retain ~10–11k samples)
- [ ] Pan-cancer Spearman ρ (CDKN2A vs IFN-γ) from R notebook 04 matches the Python
      scatter in the dashboard within ±0.01
- [ ] Hot tumors show better prognosis in known immunotherapy-responsive types
      (SKCM, BLCA) — sanity-checks biology
- [ ] `cox.zph()` output in notebook 05 shows non-significant p-values (PH assumption holds)
- [ ] No tumor type contributes >20% of total samples (check Overview tab in dashboard)
