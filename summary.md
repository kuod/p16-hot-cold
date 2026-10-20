# Analysis Summary: p16/Senescence × Hot/Cold Tumor Immune Phenotype

## Core Biological Hypothesis

Cellular senescence is a double-edged sword in cancer. Senescent cells arrest their own proliferation (tumor-suppressive) but also secrete a pro-inflammatory cytokine cocktail — the **SASP** — that should, in principle, recruit immune cells and create an "immune-hot" tumor microenvironment. The flip side: **CDKN2A deletion** (loss of p16) allows cells to bypass senescence entirely, potentially producing "cold" tumors with reduced immune infiltration.

This analysis tests whether that mechanistic link is visible at scale across 9,702 TCGA pan-cancer samples.

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
The observed pooled Spearman ρ is 0.150 in the full cohort (p = 1.3e-49; 0.165 in the confirmatory set of 25 tumor types with n ≥ 100), with a median within-tumor-type ρ of 0.112 (highest in kidney chromophobe 0.44, thyroid 0.40 and ovarian 0.36; lowest in DLBCL −0.19, stomach −0.12 and cholangiocarcinoma −0.06), inside the expected 0.10–0.25 range for H1. It is statistically robust at pan-cancer scale but not strong enough to use CDKN2A expression alone as an immune classifier. This is scientifically honest — bulk RNA-seq mixes tumor cells, stroma, and immune infiltrates, so the signal is inherently diluted.

### 2. SASP is higher in hot tumors, and the difference survives purity adjustment (attenuated)
The SASP score is higher in hot than cold tumors (three-group Kruskal-Wallis p = 8.5e-259; median −0.071 hot vs −0.097 intermediate vs −0.166 cold; hot > cold in 24 of the 26 tumor types with at least 5 tumors per group, a descriptive count rather than a formal test), but this association does not show that SASP drives inflammation. SASP correlates negatively with ABSOLUTE purity (ρ = −0.35), so stroma inflates it. The stromal deconfounding analysis (notebook 04 §11), which now runs because purity is available (ABSOLUTE, TCGA PanCanAtlas via the GDC, Taylor et al. 2018; 9,295 of 9,702 samples, 95.8%, median 0.64), shows that after residualising SASP on purity the SASP–IFN-γ correlation falls from ρ = 0.37 to 0.24 and the hot > cold difference remains (Kruskal-Wallis p = 3.2e-90; adjusted medians 0.006 hot, 0.000 intermediate, −0.043 cold). So the signal is not purely a stroma artifact, but it is attenuated, and the adjustment is a pooled linear residualisation that does not account for tumor type.

### 3. Pooled and stratified H3 disagree
Pooled across tumor types, CDKN2A deletion is *more* frequent in hot (39.2%) than cold (31.3%) tumors (tumor-type confounding). A tumor-type-stratified Mantel-Haenszel test gives OR 1.25 (95% CI 1.08–1.44, p = 2.6e-3, cold > hot in 21 of 29 strata) in the direction of H3, so H3 is still supported but the effect is weaker than before; the pooled direction (hot > cold) remains confounded by tumor type.

### Deletion ≠ expression
The 4-group analysis (deletion × expression-high/low) probes a biologically important paradox: some CDKN2A-deleted tumors still have high CDKN2A expression (presumably via the CDKN2B paralog or stromal contamination). The purity-quartile analysis (notebook 04 §9b, `figures/cdkn2a_expr_vs_del_by_purity.png`) now runs with the ABSOLUTE purity values: deleted tumors have similar or slightly lower CDKN2A expression than intact ones in every purity quartile, with wide overlap, and high-expression deleted samples persist in the highest-purity quartile, which argues against a purely low-purity artifact. This is a visual comparison only (no formal test).

### 4. H4 is not supported: no combined senescence + immune survival ordering
The 6-group KM (H4; immune tertile × p16 high/low) does **not** show the pre-specified ordering (median OS: cold/p16-low 10.8 y, hot/p16-high 7.2 y; cold/p16-high 7.3 y, hot/p16-low 8.1 y, intermediate/p16-high 6.5 y, intermediate/p16-low 7.4 y); the unstratified log-rank p = 1.5e-4 becomes p = 0.058 when stratified by tumor type, which is above the pre-specified threshold (α = 0.01), so the data do not support a combined senescence + immune survival effect. The three-group (cold/intermediate/hot) log-rank is p = 4.3e-3. Whether immune phenotype or senescence status dominates was not tested formally.

### 5. The Cox HR direction matters more than the magnitude
H5 asks whether CDKN2A_expr HR < 1.0 after adjusting for age, IFN-γ score, other senescence markers, TMB and tumor type (via `strata(tumor_type)`); ABSOLUTE tumor purity is now also in the model (n = 8,878 complete cases for purity and TMB, 2,479 events). Observed: CDKN2A HR 0.993 (95% CI 0.977–1.009, p = 0.40), so H5 is not confirmed. The IFN-γ score is nominally associated with longer survival (HR 0.930, 0.873–0.991, p = 0.024), which does not meet the α = 0.01 used for the primary hypotheses; it was 0.907 in the earlier model without purity (a different sample, so not a like-for-like comparison). Purity itself had HR 1.29 (1.007–1.651, p = 0.044). The exploratory CDKN2A × IFN-γ interaction (also adjusted for age, purity and TMB, stratified by tumor type, same n) has HR 0.961 (0.945–0.978, p = 5.8e-6). The Fine-Gray model for cancer-specific death (not stratified by tumor type, n = 8,590) gives CDKN2A sub-distribution HR 1.050 (1.028–1.072, p = 6.3e-6), above 1 and so of opposite sign to the non-significant all-cause Cox estimate; the models differ in endpoint, stratification and sample, so this is hypothesis-generating. Per-tumor-type three-group log-rank tests pass FDR < 0.05 in LGG, KIRP and PAAD (adjusted p 2.3e-14, 3.5e-5, 0.017; LUAD 0.059 and UVM 0.060 do not; LGG has only 4 hot tumors). VIF checks ensure IFN-γ/CDKN2A collinearity doesn't make those HRs uninterpretable.

---

## Design Decisions Worth Noting

- **Why pan-cancer?** Tumor-type-specific effects average out at scale, and `strata(tumor_type)` in Cox lets the baseline hazard vary per cancer type without one-hot encoding 32 types.
- **Why p16 specifically?** CDKN2A is the most frequently deleted gene across all cancers (TCGA), making it uniquely testable. Its dual role as CDK inhibitor and SASP regulator makes it mechanistically central to the hot/cold question.
- **Multiple testing:** Bonferroni correction (α = 0.01) across the five primary hypotheses; BH FDR for all exploratory analyses. These thresholds are pre-specified, not chosen post-hoc.

---

## Limitations

- **Bulk RNA-seq cannot localize the signal.** CDKN2A-expressing cells could be tumor cells, senescent fibroblasts, or infiltrating immune cells. Single-cell follow-up would be needed to resolve this.
- **Purity confounding.** High-stroma (low-purity) samples inflate immune scores. ABSOLUTE purity (95.8% coverage) is now in the Cox and Fine-Gray models and in a purity-adjusted SASP / IFN-γ sensitivity analysis, but the H1 correlation, the six-group survival analysis and the deletion tests are not purity-adjusted, and purity estimates are themselves imperfect. Stromal confounding is reduced, not eliminated.
- **Observational design.** Survival associations are correlational. CDKN2A expression may proxy other features of well-differentiated, immune-infiltrated tumors rather than directly reflecting senescence biology.

---

## What Would Strengthen the Findings

- Validation in an independent cohort (e.g., ICGC, or a single-institution RNA-seq dataset)
- Single-cell RNA-seq to localize p16 expression to tumor vs stromal compartments
- Mechanistic in vitro confirmation that p16 knockdown reduces SASP cytokine secretion and immune recruitment

---

## One-Sentence Summary

In pan-cancer TCGA data, CDKN2A expression and SASP scores correlate modestly with an IFN-γ hot-tumor signature; the hypothesised survival benefit of p16-high/hot tumors is not supported once tumor type is accounted for, and the SASP association survives, in attenuated form, adjustment for ABSOLUTE tumor purity.
