"""
Curated gene signatures used throughout the analysis.
All gene symbols are HGNC-approved.
"""

# Standard hot-tumor / IFN-gamma proxy (Ayers et al. 2017, 6-gene version)
IFNG_SIGNATURE = ["IFNG", "STAT1", "CCL5", "CXCL9", "CXCL10", "IDO1"]

# Senescence-associated secretory phenotype (SASP) core genes
SASP_GENES = [
    "IL6", "CXCL8",  # IL8 official symbol is CXCL8
    "MMP3", "MMP1", "MMP9", "MMP13",
    "VEGFA", "CXCL1", "CCL2", "CCL7",
    "IL1A", "IL1B", "CXCL2", "CXCL3",
    "IGFBP3", "IGFBP7",
]

# Core senescence effectors / cell-cycle arrest
SENESCENCE_EFFECTORS = [
    "CDKN2A",   # p16^INK4a
    "CDKN1A",   # p21^CIP1/WAF1
    "TP53",
    "RB1",
    "CDKN2B",   # p15^INK4b
    "GLB1",     # beta-galactosidase (SA-β-gal proxy)
    "LMNB1",    # loss = senescence marker
    "H2AFX",    # gamma-H2AX surrogate (DNA damage); HGNC: H2AFX, alias H2AX
]

# Immune checkpoint / T cell exhaustion
CHECKPOINT_GENES = [
    "PDCD1",    # PD-1
    "CD274",    # PD-L1
    "CTLA4",
    "LAG3",
    "HAVCR2",   # TIM-3
    "TIGIT",
]

# CD8+ cytotoxic T cell markers
CD8_TCELL_GENES = ["CD8A", "CD8B", "GZMB", "PRF1", "IFNG"]

# Individual expression markers tracked alongside CDKN2A_expr
SENESCENCE_EXPRESSION_MARKERS = ["CDKN1A", "SERPINE1", "GDF15", "LMNB1"]

# Immunosuppressive arm of SASP — pro-tumorigenic, T-cell suppressive
SASP_SUPPRESSIVE = [
    "IL10",      # anti-inflammatory cytokine, suppresses CD8+ T cells
    "TGFB1",     # master immune suppressor; also induces p15 (CDKN2B)
    "TGFB2",
    "CXCL12",    # SDF-1; drives T-cell exclusion in cold tumors
    "VEGFA",     # pro-angiogenic + immune exclusion; also in SASP_GENES (pro-inflammatory role)
                 # dual role: SASP_GENES scores its senescence-promoting arm;
                 # SASP_SUPPRESSIVE scores its immunosuppressive/exclusion arm
    "PTGS2",     # COX-2 / prostaglandin E2 — immunosuppressive
    "IDO1",      # tryptophan catabolism → T cell anergy
]

# Minimal 3-gene senescence marker score (harder to explain by stromal contamination)
# p16 expression + p21 expression + LMNB1 loss (downregulated in senescence)
SENESCENCE_MINIMAL_SCORE = ["CDKN2A", "CDKN1A", "LMNB1"]

# Senescence pathway genes to pull from cBioPortal CNA + mutation profiles
# Format: (gene_name, entrez_id, fetch_cna, fetch_mut)
# Covers: cell-cycle brakes (INK4/CIP), DNA-damage sensors, and p53 pathway
SENESCENCE_GENOMIC_LOCI = [
    ("tp53",   7157, True, True),   # p53 — central senescence effector
    ("rb1",    5925, True, True),   # Rb — downstream of p16
    ("atm",    472,  True, True),   # DNA damage sensor → p53 activation
    ("pten",   5728, True, True),   # PTEN loss → PI3K/AKT → senescence bypass
    ("cdkn1a", 1026, True, False),  # p21^CIP1 — p53 target, CDK inhibitor
    ("cdkn2b", 1030, True, False),  # p15^INK4b — paralog of p16
    ("cdkn2c", 1031, True, False),  # p18^INK4c — INK4 family
]

# Note: SASP_SUPPRESSIVE above uses custom genes (not MSigDB) and will be scored
# via ssGSEA by passing the gene list directly (not via MSIGDB_SETS lookup).

# MSigDB gene set names to query via gseapy
MSIGDB_SETS = {
    "senescence_up": "FRIDMAN_SENESCENCE_UP",
    "reactome_senescence": "REACTOME_CELLULAR_SENESCENCE",
    "hallmark_ifng": "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "hallmark_inflammatory": "HALLMARK_INFLAMMATORY_RESPONSE",
    "hallmark_tnfa": "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "senmayo":         "SAUL_SEN_MAYO",            # Saul et al. 2022 Nature Aging
    "hallmark_p53":    "HALLMARK_P53_PATHWAY",
    "gobp_senescence": "GOBP_CELLULAR_SENESCENCE",  # C5 GO biological process
    "senescence_dn":   "FRIDMAN_SENESCENCE_DN",     # downregulated in senescence
}
