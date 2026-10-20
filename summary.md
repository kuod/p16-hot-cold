# Analysis Summary: p16/Senescence × Hot/Cold Tumor Immune Phenotype

## Core Biological Hypothesis

Cellular senescence is a double-edged sword in cancer. Senescent cells arrest their own proliferation (tumor-suppressive) but also secrete a pro-inflammatory cytokine cocktail — the **SASP** — that should, in principle, recruit immune cells and create an "immune-hot" tumor microenvironment. The flip side: **CDKN2A deletion** (loss of p16) allows cells to bypass senescence entirely, potentially producing "cold" tumors with reduced immune infiltration.

This analysis tests whether that mechanistic link is visible at scale across ~9,900 TCGA pan-cancer samples.

---

## Pre-Specified Hypotheses (H1–H5)

| | Hypothesis | Expected finding | Implication |
|--|--|--|--|
| **H1** | CDKN2A expression correlates positively with IFN-γ score | ρ ~ 0.10–0.25 | p16 expression tracks immune activity pan-cancer |
| **H2** | SASP score higher in hot vs cold tumors | Significant, p < 0.01 | SASP drives immune infiltration |
| **H3** | CDKN2A deletion enriched in cold tumors | Significant, p < 0.01 | Senescence loss → immune exclusion |
| **H4** | p16-high/hot tumors have best OS; p16-low/cold have worst | Global log-rank p < 0.01 | Combined senescence + immune phenotype predicts survival |
| **H5** | CDKN2A expression is an independent OS predictor (Cox, age-adjusted) | HR < 1.0 | p16 expression is prognostically protective, independent of tumor type |

---

## Key Conclusions

### 1. The correlation is real but modest
The observed pooled Spearman ρ is 0.158 (median within-tumor-type ρ 0.112), inside the expected 0.10–0.25 range for H1. It is statistically robust at pan-cancer scale but not strong enough to use CDKN2A expression alone as an immune classifier. This is scientifically honest — bulk RNA-seq mixes tumor cells, stroma, and immune infiltrates, so the signal is inherently diluted.

### 2. SASP is higher in hot tumors, but purity confounding could not be assessed
The SASP score is higher in hot than cold tumors (three-group Kruskal-Wallis), but this association does not show that SASP drives inflammation. The stromal deconfounding analysis (notebook 04 §11) is designed to test whether the SASP → hot association survives tumor-purity adjustment, but tumor purity data is unavailable (the download returns HTTP 403), so that analysis is **not run** and no purity adjustment exists. Until it is, the SASP signal cannot be called tumor-intrinsic rather than a stroma artifact.

### 3. Pooled and stratified H3 disagree
Pooled across tumor types, CDKN2A deletion is *more* frequent in hot (38.6%) than cold (31.7%) tumors (tumor-type confounding). A tumor-type-stratified Mantel-Haenszel test gives OR 1.30 (95% CI 1.12–1.50, p = 4e-4) in the direction of H3.

### Deletion ≠ expression
The 4-group analysis (deletion × expression-high/low) probes a biologically important paradox: some CDKN2A-deleted tumors still have high CDKN2A expression (presumably via the CDKN2B paralog or stromal contamination). The purity-quartile analysis intended to test whether this is artifactual is **not run** (no purity data).

### 4. H4 is not supported: no combined senescence + immune survival ordering
The 6-group KM (H4; immune tertile × p16 high/low) does **not** show the pre-specified ordering (median OS: cold/p16-low 9.5 y, hot/p16-high 7.8 y); the unstratified log-rank p = 3.3e-3 becomes p = 0.077 when stratified by tumor type, which is above the pre-specified threshold (α = 0.01), so the data do not support a combined senescence + immune survival effect. Whether immune phenotype or senescence status dominates was not tested formally.

### 5. The Cox HR direction matters more than the magnitude
H5 asks whether CDKN2A_expr HR < 1.0 after adjusting for age, IFN-γ score, other senescence markers, TMB and tumor type (via `strata(tumor_type)`); tumor purity is **not** in the model because purity data is unavailable. Observed: CDKN2A HR 0.989 (95% CI 0.973–1.005, p = 0.16), so H5 is not confirmed, while IFN-γ score is independently protective (HR 0.907, p = 8.5e-4). VIF checks ensure IFN-γ/CDKN2A collinearity doesn't make those HRs uninterpretable.

---

## Design Decisions Worth Noting

- **Why pan-cancer?** Tumor-type-specific effects average out at scale, and `strata(tumor_type)` in Cox lets the baseline hazard vary per cancer type without one-hot encoding 32 types.
- **Why p16 specifically?** CDKN2A is the most frequently deleted gene across all cancers (TCGA), making it uniquely testable. Its dual role as CDK inhibitor and SASP regulator makes it mechanistically central to the hot/cold question.
- **Multiple testing:** Bonferroni correction (α = 0.01) across the five primary hypotheses; BH FDR for all exploratory analyses. These thresholds are pre-specified, not chosen post-hoc.

---

## Limitations

- **Bulk RNA-seq cannot localize the signal.** CDKN2A-expressing cells could be tumor cells, senescent fibroblasts, or infiltrating immune cells. Single-cell follow-up would be needed to resolve this.
- **Purity confounding.** High-stroma (low-purity) samples inflate immune scores. Purity data is unavailable, so no purity-adjusted sensitivity analysis has been run and all immune correlations remain unadjusted for it.
- **Observational design.** Survival associations are correlational. CDKN2A expression may proxy other features of well-differentiated, immune-infiltrated tumors rather than directly reflecting senescence biology.

---

## What Would Strengthen the Findings

- Validation in an independent cohort (e.g., ICGC, or a single-institution RNA-seq dataset)
- Single-cell RNA-seq to localize p16 expression to tumor vs stromal compartments
- Mechanistic in vitro confirmation that p16 knockdown reduces SASP cytokine secretion and immune recruitment

---

## One-Sentence Summary

In pan-cancer TCGA data, CDKN2A expression and SASP scores correlate modestly with an IFN-γ hot-tumor signature; the hypothesised survival benefit of p16-high/hot tumors is not supported once tumor type is accounted for, and purity confounding has not yet been assessed (purity data unavailable).
