# scRNA-seq Follow-up: Localizing the CDKN2A/SASP Signal

The central unanswered question from the bulk analysis is **which cell type** in the TME
expresses p16 and drives SASP. Bulk RNA-seq cannot resolve this — these datasets can.

---

## Priority 1: Localize CDKN2A/SASP to a Cell Type

### Pan-cancer stromal / CAF atlases
| Study | Tumor types | Cells | Key use |
|---|---|---|---|
| Qian et al. 2020 *(Cell Research)* | 10 types, paired tumor/normal | ~200k | Which compartment (epithelial/stromal/immune) is CDKN2A highest? |
| Elyada et al. 2019 *(Cancer Discovery)* | Pancreatic | — | iCAFs vs myCAFs; iCAFs are major SASP producers and likely p16-high |
| Dominguez et al. 2020 *(Cell)* | Multi-tissue fibroblasts | — | Contains a senescent-like fibroblast cluster (CDKN2A/CDKN1A/IGFBP7) — direct validation target |

### Senescence-specific
- **Saul et al. 2022** *(Nature Aging)* — the SenMayo gene set paper; deposited scRNA-seq with annotated senescent clusters. Score your `SENESCENCE_EFFECTORS` and `SASP_GENES` lists against their clusters as a direct cross-validation.
- Search CellxGene for **"senescence"** — datasets from aged lung and liver contexts contain senescence-enriched populations with established marker annotations.

---

## Priority 2: Immune Exclusion Mechanism (cold tumor / CDKN2A-deleted arm)

| Study | Tumor type | Key use |
|---|---|---|
| Jerby-Arnon et al. 2018 *(Cell)* | Melanoma | Tumor cell exclusion program — check anti-correlation with CDKN2A in malignant cells |
| Pelka et al. 2021 *(Cell)* | Colorectal (~370k cells) | Detailed myeloid/T cell states + spatial component; CRC has strong hot/cold biology |
| Gavish et al. 2023 *(Nature)* | Pan-cancer, 24 types | Pan-cancer malignant cell programs — look for a senescent malignant meta-program |

---

## Priority 3: Validate H2/H3 at Single-Cell Resolution

| Study | Cells | Key use |
|---|---|---|
| Salcher et al. 2022 — TCIA *(Cancer Cell)* | ~500k immune, 13 types | Do SASP scores in myeloid/stromal clusters track CD8 infiltration per tumor type? Mirrors bulk H2 |
| Zheng et al. 2021 *(Science)* | ~300k T cells, pan-cancer | Does your IFN-γ signature map onto their "hot" T cell states vs exclusion states? |

---

## CellxGene Search Strategy

At **cellxgene.cziscience.com**:

1. Filter **disease = "cancer"** + **tissue = [top tumor types from per-type log-rank results in notebook 05]**
2. Sort by cell count — prioritize >50k cells for stable rare-cluster detection
3. For the senescence question, filter **cell type = "fibroblast"** — CAFs are the most plausible bulk SASP source and are well-annotated in most CellxGene cancer collections
4. HTAN collections (pancreatic, colorectal, lung) are the most mature and have stromal annotations detailed enough to score SASP/senescence signatures

---

## The Specific Experiment to Run

Given a dataset with annotated cell types:

**1. Score cells with your existing signatures**
```python
import scanpy as sc

for name, genes in [("SASP", SASP_GENES), ("senescence", SENESCENCE_EFFECTORS)]:
    sc.tl.score_genes(adata, genes, score_name=name)
```
Uses the same mean z-score logic as the bulk pipeline — directly comparable.

**2. Identify which cell type carries the signal**
```python
sc.pl.violin(adata, ["SASP", "CDKN2A"], groupby="cell_type")
```
Expect CAFs and potentially senescent epithelial cells to top the list.

**3. Test whether stromal senescence drives the bulk correlation**

Per sample, correlate:
- **fraction of SASP-high CAFs** (scRNA-seq) vs **bulk IFN-γ score** (from notebook 02)

If this per-sample correlation is strong, stromal senescence is the likely driver of the H1 bulk signal — not tumor-cell-intrinsic p16.

**4. Validate H3 at cell-type resolution**

Compare CDKN2A expression in malignant cells between samples classified as hot vs cold by your IFN-γ threshold. If the difference disappears in malignant cells but persists in CAFs, the deletion signal in bulk is a stromal phenomenon.

---

## Interpretation Framework

| Finding | Conclusion |
|---|---|
| CDKN2A/SASP highest in CAFs | Bulk signal is stromal; p16 in fibroblasts drives immune recruitment |
| CDKN2A/SASP highest in malignant cells | Tumor-cell-intrinsic senescence drives hot phenotype |
| CDKN2A/SASP highest in myeloid cells | Macrophage senescence (inflammaging) drives IFN-γ signal |
| Signal diffuse across cell types | Bulk correlation reflects tissue-level co-occurrence, not a single driver |

The stromal hypothesis (CAFs) is the most biologically plausible given SASP biology and would be the strongest story for a follow-up paper.
