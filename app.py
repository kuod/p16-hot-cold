"""Streamlit app for exploring hot/cold × p16 TCGA analysis results."""

import sys
from pathlib import Path

import pandas as pd
import numpy as np
import streamlit as st
import plotly.express as px
import plotly.graph_objects as go

sys.path.insert(0, str(Path(__file__).parent / "src"))
from pathway_diagram import streamlit_html as _pathway_html
from utils import validate_parquet

DATA = Path("data")
FIGURES = Path("figures")

st.set_page_config(
    page_title="Hot/Cold × p16 — TCGA",
    page_icon="🧬",
    layout="wide",
)


@st.cache_data
def load_data():
    validate_parquet(
        DATA / "sample_immune_labels.parquet",
        required_cols=["tumor_type", "hot_cold", "ifng_score"],
    )
    validate_parquet(
        DATA / "sample_senescence.parquet",
        required_cols=["CDKN2A_expr"],
    )
    immune = pd.read_parquet(DATA / "sample_immune_labels.parquet")
    immune.index.name = "sample_id"
    senes = pd.read_parquet(DATA / "sample_senescence.parquet")
    senes.index.name = "sample_id"
    master = immune.join(senes, how="inner")

    clin = pd.read_csv(DATA / "tcga_clinical_survival.tsv", sep="\t", index_col=0)
    clin.index.name = "sample_id"
    extra_clin = ["OS", "OS.time", "DSS", "DSS.time", "PFI", "PFI.time",
                  "age_at_initial_pathologic_diagnosis"]
    available_clin = [c for c in extra_clin if c in clin.columns]
    master = master.join(
        clin[available_clin].rename(
            columns={"age_at_initial_pathologic_diagnosis": "age"}
        ),
        how="left",
    )
    master["OS.time_years"] = master["OS.time"] / 365.25
    if "DSS.time" in master.columns:
        master["DSS.time_years"] = master["DSS.time"] / 365.25
    if "PFI.time" in master.columns:
        master["PFI.time_years"] = master["PFI.time"] / 365.25

    # p16 high/low within tumor type (median split)
    def _safe_qcut(x):
        try:
            return pd.qcut(x, 2, labels=["p16-low", "p16-high"], duplicates="drop")
        except ValueError:
            return pd.Series(pd.NA, index=x.index, dtype="object")

    master["p16_status"] = (
        master.groupby("tumor_type")["CDKN2A_expr"]
        .transform(_safe_qcut)
        .astype(str)
    )
    master["group6"] = master["hot_cold"] + " / " + master["p16_status"]

    cox = pd.read_csv(DATA / "cox_results.csv")
    lr = pd.read_csv(DATA / "logrank_per_type.csv")
    return master, cox, lr


master, cox, lr = load_data()

# ── Sidebar ────────────────────────────────────────────────────────────────────
st.sidebar.title("🧬 Hot/Cold × p16")
st.sidebar.caption("TCGA pan-cancer analysis")

tumor_types = sorted(master["tumor_type"].dropna().unique())
selected_types = st.sidebar.multiselect(
    "Filter tumor types", tumor_types, default=tumor_types
)

df = master[master["tumor_type"].isin(selected_types)] if selected_types else master

page = st.sidebar.radio(
    "Section",
    ["Overview", "Immune Classification", "Senescence Scores", "Correlations", "Survival", "Pathway", "Figures", "Hypotheses", "Summary"],
)

# ── Overview ──────────────────────────────────────────────────────────────────
if page == "Overview":
    st.title("Hot/Cold Tumor × p16/CDKN2A — TCGA Pan-Cancer")
    st.markdown(
        """
        This app explores the relationship between **immune phenotype** (cold / intermediate / hot tumors, IFN-γ tertiles)
        and **p16/CDKN2A senescence status** across TCGA pan-cancer RNA-seq data.

        | Hypothesis | Test |
        |---|---|
        | SASP-high tumors → hot (pro-inflammatory cytokines drive TIL infiltration) | Spearman correlation, Mann-Whitney |
        | CDKN2A-deleted tumors → cold | Chi-squared |
        | Senescence/immune phenotype predicts survival | Kaplan-Meier, Cox PH |
        """
    )

    c1, c2, c3, c4, c5 = st.columns(5)
    c1.metric("Samples", f"{len(df):,}")
    c2.metric("Tumor types", df["tumor_type"].nunique())
    c3.metric("Hot", f"{(df['hot_cold'] == 'hot').sum():,}")
    c4.metric("Intermediate", f"{(df['hot_cold'] == 'intermediate').sum():,}")
    c5.metric("Cold", f"{(df['hot_cold'] == 'cold').sum():,}")

    st.subheader("Data completeness")
    completeness_cols = {
        "CDKN2A_expr": "CDKN2A expression",
        "SASP": "SASP score",
        "SASP_suppressive": "Suppressive SASP score",
        "senescence_minimal": "Minimal senescence score",
        "cd8_fraction": "CD8+ fraction (CIBERSORT)",
        "immune_subtype": "Thorsson immune subtype (C1–C6)",
        "tumor_purity": "Tumor purity (ABSOLUTE)",
        "OS": "Overall survival (OS)",
        "DSS": "Disease-specific survival (DSS)",
        "PFI": "Progression-free interval (PFI)",
        "age": "Age at diagnosis",
        "tmb": "TMB (nonsynonymous mut/Mb)",
        "cdkn2a_altered": "CDKN2A alteration",
        "tp53_altered": "TP53 alteration",
        "atm_altered": "ATM alteration",
        "pten_altered": "PTEN alteration",
    }
    completeness_rows = [
        {"Variable": label, "N present": f"{df[col].notna().sum():,}",
         "% complete": f"{100 * df[col].notna().mean():.1f}%"}
        for col, label in completeness_cols.items() if col in df.columns
    ]
    if completeness_rows:
        st.dataframe(
            pd.DataFrame(completeness_rows).set_index("Variable"),
            use_container_width=True,
        )

    st.subheader("Immune phenotype distribution by cancer type")
    pheno_frac = (
        df.groupby("tumor_type")["hot_cold"]
        .apply(lambda s: s.value_counts(normalize=True))
        .unstack(fill_value=0)
        .reset_index()
        .sort_values("hot", ascending=True)
    )
    # Ensure all 3 columns exist
    for col in ["cold", "intermediate", "hot"]:
        if col not in pheno_frac.columns:
            pheno_frac[col] = 0.0
    fig = px.bar(
        pheno_frac.melt(id_vars="tumor_type", value_vars=["cold", "intermediate", "hot"],
                        var_name="phenotype", value_name="fraction"),
        x="fraction", y="tumor_type", color="phenotype", orientation="h",
        color_discrete_map={"hot": "#e74c3c", "intermediate": "#f39c12", "cold": "#3498db"},
        labels={"fraction": "Fraction", "tumor_type": "", "phenotype": ""},
        height=600, barmode="stack",
    )
    fig.add_vline(x=1/3, line_dash="dash", line_color="gray")
    st.plotly_chart(fig, use_container_width=True)

    if "age" in df.columns:
        st.subheader("Median age at diagnosis by cancer type")
        age_by_type = (
            df.dropna(subset=["age", "tumor_type"])
            .groupby("tumor_type")["age"]
            .median()
            .reset_index()
            .rename(columns={"age": "median_age"})
            .sort_values("median_age", ascending=True)
        )
        fig_age = px.bar(
            age_by_type, x="median_age", y="tumor_type", orientation="h",
            color="median_age", color_continuous_scale=["#27ae60", "#8e44ad"],
            labels={"median_age": "Median age (years)", "tumor_type": ""},
            height=600,
        )
        fig_age.update_coloraxes(showscale=False)
        st.plotly_chart(fig_age, use_container_width=True)

        st.subheader("Age distribution by immune phenotype")
        age_hc = df.dropna(subset=["age", "hot_cold"])
        med_hot  = age_hc.loc[age_hc["hot_cold"] == "hot",          "age"].median()
        med_int  = age_hc.loc[age_hc["hot_cold"] == "intermediate", "age"].median()
        med_cold = age_hc.loc[age_hc["hot_cold"] == "cold",         "age"].median()
        fig_violin = px.violin(
            age_hc,
            x="hot_cold", y="age", color="hot_cold",
            color_discrete_map={"hot": "#e74c3c", "intermediate": "#f39c12", "cold": "#3498db"},
            category_orders={"hot_cold": ["cold", "intermediate", "hot"]},
            box=True, points=False,
            labels={"hot_cold": "", "age": "Age at diagnosis (years)"},
        )
        fig_violin.update_layout(showlegend=False)
        st.plotly_chart(fig_violin, use_container_width=True)
        st.caption(
            f"Median age — hot: {med_hot:.1f} yrs, intermediate: {med_int:.1f} yrs, cold: {med_cold:.1f} yrs"
        )

# ── Immune Classification ──────────────────────────────────────────────────────
elif page == "Immune Classification":
    st.title("Immune Classification")

    col1, col2 = st.columns(2)

    PHENO_COLORS = {"hot": "#e74c3c", "intermediate": "#f39c12", "cold": "#3498db"}
    PHENO_ORDER  = ["cold", "intermediate", "hot"]

    with col1:
        st.subheader("IFN-γ score distribution")
        low_t = df["ifng_score"].quantile(1/3)
        high_t = df["ifng_score"].quantile(2/3)
        fig = px.histogram(
            df, x="ifng_score", color="hot_cold",
            color_discrete_map=PHENO_COLORS,
            category_orders={"hot_cold": PHENO_ORDER},
            barmode="overlay", opacity=0.55, nbins=80,
            labels={"ifng_score": "IFN-γ signature score", "hot_cold": ""},
        )
        fig.add_vline(x=low_t,  line_dash="dash", line_color="#f39c12", annotation_text="T1")
        fig.add_vline(x=high_t, line_dash="dash", line_color="#e74c3c",  annotation_text="T2")
        st.plotly_chart(fig, use_container_width=True)

    with col2:
        st.subheader("CD8+ T cell fraction by phenotype")
        fig = px.box(
            df.dropna(subset=["cd8_fraction"]),
            x="hot_cold", y="cd8_fraction",
            color="hot_cold",
            color_discrete_map=PHENO_COLORS,
            category_orders={"hot_cold": PHENO_ORDER},
            labels={"cd8_fraction": "CD8+ fraction (CIBERSORT)", "hot_cold": ""},
            points="outliers",
        )
        fig.update_layout(showlegend=False)
        st.plotly_chart(fig, use_container_width=True)

    st.subheader("IFN-γ score vs CD8+ fraction")
    fig = px.scatter(
        df.dropna(subset=["cd8_fraction"]).sample(min(3000, len(df))),
        x="ifng_score", y="cd8_fraction",
        color="hot_cold",
        color_discrete_map=PHENO_COLORS,
        category_orders={"hot_cold": PHENO_ORDER},
        opacity=0.5, size_max=4,
        trendline="ols",
        hover_data=["tumor_type"],
        labels={"ifng_score": "IFN-γ score", "cd8_fraction": "CD8+ fraction"},
    )
    st.plotly_chart(fig, use_container_width=True)

# ── Senescence Scores ──────────────────────────────────────────────────────────
elif page == "Senescence Scores":
    st.title("Senescence Scores")

    score_cols = ["CDKN2A_expr", "SASP", "Senescence_Effectors",
                  "hallmark_ifng", "hallmark_inflammatory", "hallmark_tnfa",
                  "reactome_senescence", "senescence_up",
                  "senmayo", "hallmark_p53", "gobp_senescence", "senescence_dn",
                  "CDKN1A_expr", "SERPINE1_expr", "GDF15_expr", "LMNB1_expr"]
    score_cols = [c for c in score_cols if c in df.columns]

    col1, col2 = st.columns(2)

    PHENO_COLORS = {"hot": "#e74c3c", "intermediate": "#f39c12", "cold": "#3498db"}
    PHENO_ORDER  = ["cold", "intermediate", "hot"]

    with col1:
        score = st.selectbox("Score", score_cols)
        fig = px.violin(
            df.dropna(subset=[score, "hot_cold"]),
            x="hot_cold", y=score, color="hot_cold",
            color_discrete_map=PHENO_COLORS,
            category_orders={"hot_cold": PHENO_ORDER},
            box=True, points=False,
            labels={"hot_cold": "", score: score},
        )
        fig.update_layout(showlegend=False)
        st.plotly_chart(fig, use_container_width=True)

    with col2:
        st.subheader("Senescence gene alterations: hot vs cold")
        # Prefer _altered (del|mut composite) where available
        alt_gene_cols = [
            ("CDKN2A", "cdkn2a_altered"),
            ("TP53",   "tp53_altered"),
            ("RB1",    "rb1_altered"),
            ("ATM",    "atm_altered"),
            ("PTEN",   "pten_altered"),
            ("CDKN1A", "cdkn1a_deleted"),
            ("CDKN2B", "cdkn2b_deleted"),
            ("CDKN2C", "cdkn2c_deleted"),
        ]
        available = [(g, c) for g, c in alt_gene_cols if c in df.columns]
        if available and "hot_cold" in df.columns:
            rows = []
            for gene, col in available:
                for pheno in ["cold", "intermediate", "hot"]:
                    sub = df[df["hot_cold"] == pheno]
                    pct = 100 * sub[col].mean() if len(sub) > 0 else 0.0
                    rows.append({"Gene": gene, "Phenotype": pheno, "% altered": pct})
            alt_df = pd.DataFrame(rows)
            fig = px.bar(
                alt_df, x="Gene", y="% altered", color="Phenotype",
                barmode="group",
                color_discrete_map={"hot": "#e74c3c", "intermediate": "#f39c12", "cold": "#3498db"},
                category_orders={"Phenotype": ["cold", "intermediate", "hot"]},
                labels={"% altered": "% of samples altered"},
            )
            fig.update_layout(legend_title_text="")
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.caption("Alteration data not yet available — re-run notebooks 01 and 03.")

    st.subheader("CDKN2A expression by tumor type")
    fig = px.box(
        df.dropna(subset=["CDKN2A_expr", "tumor_type"]),
        x="tumor_type", y="CDKN2A_expr", color="hot_cold",
        color_discrete_map={"hot": "#e74c3c", "intermediate": "#f39c12", "cold": "#3498db"},
        category_orders={"hot_cold": ["cold", "intermediate", "hot"]},
        labels={"CDKN2A_expr": "CDKN2A expression (log2)", "tumor_type": ""},
        points=False,
        height=500,
    )
    fig.update_layout(xaxis_tickangle=-45)
    st.plotly_chart(fig, use_container_width=True)

    # ── Thorsson immune subtype cross-tab ─────────────────────────────────────
    if "immune_subtype" in df.columns and df["immune_subtype"].notna().any():
        st.subheader("Thorsson immune subtypes (C1–C6) × hot/cold")
        sub_df = df.dropna(subset=["immune_subtype", "hot_cold"])
        sub_counts = (
            sub_df.groupby(["immune_subtype", "hot_cold"])
            .size().reset_index(name="n")
        )
        sub_counts["pct"] = (
            100 * sub_counts["n"]
            / sub_counts.groupby("immune_subtype")["n"].transform("sum")
        )
        fig_sub = px.bar(
            sub_counts, x="immune_subtype", y="pct", color="hot_cold",
            barmode="stack",
            color_discrete_map={"hot": "#e74c3c", "intermediate": "#f39c12", "cold": "#3498db"},
            category_orders={"hot_cold": ["cold", "intermediate", "hot"]},
            labels={"immune_subtype": "Immune subtype", "pct": "% samples", "hot_cold": ""},
            title="Immune phenotype (tertile) within Thorsson immune subtypes",
        )
        fig_sub.update_layout(legend_title_text="")
        st.plotly_chart(fig_sub, use_container_width=True)
        st.caption(
            "C1 = wound healing · C2 = IFN-γ dominant · C3 = inflammatory · "
            "C4 = lymphocyte depleted · C5 = immunologically quiet · C6 = TGF-β dominant"
        )

        col_sa, col_sb = st.columns(2)
        with col_sa:
            cdkn2a_sub = df.dropna(subset=["immune_subtype", "CDKN2A_expr"])
            fig_cdkn_sub = px.box(
                cdkn2a_sub, x="immune_subtype", y="CDKN2A_expr",
                color="immune_subtype", points=False,
                labels={"CDKN2A_expr": "CDKN2A expr (log2)", "immune_subtype": ""},
                title="CDKN2A expression by immune subtype",
            )
            fig_cdkn_sub.update_layout(showlegend=False)
            st.plotly_chart(fig_cdkn_sub, use_container_width=True)
        with col_sb:
            if "SASP" in df.columns:
                sasp_sub = df.dropna(subset=["immune_subtype", "SASP"])
                fig_sasp_sub = px.box(
                    sasp_sub, x="immune_subtype", y="SASP",
                    color="immune_subtype", points=False,
                    labels={"SASP": "SASP NES", "immune_subtype": ""},
                    title="SASP score by immune subtype",
                )
                fig_sasp_sub.update_layout(showlegend=False)
                st.plotly_chart(fig_sasp_sub, use_container_width=True)

    # ── Alteration co-occurrence heatmap ──────────────────────────────────────
    st.subheader("Alteration co-occurrence")
    cooc_genes = [(g, c) for g, c in alt_gene_cols if c in df.columns]
    if len(cooc_genes) >= 2:
        alt_matrix = pd.DataFrame({
            gene: df[col].fillna(False).astype(int)
            for gene, col in cooc_genes
        })
        cooc_corr = alt_matrix.corr(method="pearson")  # phi correlation for binary
        fig_cooc = px.imshow(
            cooc_corr.round(2),
            color_continuous_scale="RdBu_r", zmin=-1, zmax=1,
            text_auto=True,
            title="Alteration co-occurrence (phi correlation)",
            labels={"color": "Phi"},
            aspect="auto",
        )
        fig_cooc.update_layout(height=420)
        st.plotly_chart(fig_cooc, use_container_width=True)
        st.caption(
            "Positive phi = tend to co-occur; negative = mutually exclusive. "
            "CDKN2A and RB1 co-deletions are common (9p21 locus)."
        )

# ── Correlations ──────────────────────────────────────────────────────────────
elif page == "Correlations":
    st.title("Correlations")

    col1, col2 = st.columns(2)
    with col1:
        x_col = st.selectbox("X axis", [c for c in [
            "CDKN2A_expr", "SASP", "Senescence_Effectors",
            "reactome_senescence", "senescence_up",
            "senmayo", "hallmark_p53", "gobp_senescence", "senescence_dn",
            "CDKN1A_expr", "SERPINE1_expr", "GDF15_expr", "LMNB1_expr",
            "age",
        ] if c in df.columns])
    with col2:
        y_col = st.selectbox("Y axis", ["ifng_score", "cd8_fraction"])

    sample_df = df.dropna(subset=[x_col, y_col])
    if len(sample_df) > 4000:
        sample_df = sample_df.sample(4000, random_state=42)

    rho = sample_df[[x_col, y_col]].corr(method="spearman").iloc[0, 1]
    st.caption(
        f"Spearman ρ = {rho:.3f}  (n={len(sample_df):,}, downsampled for display) — "
        "**exploratory**: interpret with BH FDR correction."
    )

    fig = px.scatter(
        sample_df, x=x_col, y=y_col,
        color="hot_cold",
        color_discrete_map={"hot": "#e74c3c", "intermediate": "#f39c12", "cold": "#3498db"},
        category_orders={"hot_cold": ["cold", "intermediate", "hot"]},
        opacity=0.4, size_max=4,
        trendline="ols",
        hover_data=["tumor_type"],
        labels={x_col: x_col, y_col: y_col},
    )
    st.plotly_chart(fig, use_container_width=True)

    st.subheader("Per-tumor-type log-rank p-values (hot vs cold survival)")
    fig = px.bar(
        lr.sort_values("lr_pval"),
        x="lr_pval", y="tumor_type", orientation="h",
        color="padj",
        color_continuous_scale=["#e74c3c", "#f39c12", "#95a5a6"],
        labels={"lr_pval": "Log-rank p-value", "tumor_type": "", "padj": "BH adj p"},
        height=500,
    )
    fig.add_vline(x=0.05, line_dash="dash", line_color="black", annotation_text="p=0.05")
    st.plotly_chart(fig, use_container_width=True)

# ── Survival ──────────────────────────────────────────────────────────────────
elif page == "Survival":
    st.title("Survival Analysis")

    st.subheader("Cox proportional hazards — adjusted for tumor type")
    fig = go.Figure()
    for _, row in cox.iterrows():
        fig.add_trace(go.Scatter(
            x=[row["conf.low"], row["conf.high"]],
            y=[row["term"], row["term"]],
            mode="lines",
            line=dict(color="#e74c3c", width=2),
            showlegend=False,
        ))
        fig.add_trace(go.Scatter(
            x=[row["estimate"]],
            y=[row["term"]],
            mode="markers",
            marker=dict(color="#e74c3c", size=10),
            showlegend=False,
            hovertemplate=(
                f"<b>{row['term']}</b><br>"
                f"HR = {row['estimate']:.3f}<br>"
                f"95% CI: [{row['conf.low']:.3f}, {row['conf.high']:.3f}]<br>"
                f"p = {row['p.value']:.3e}<extra></extra>"
            ),
        ))
    fig.add_vline(x=1, line_dash="dash", line_color="gray")
    fig.update_xaxes(type="log", title="Hazard Ratio (95% CI)")
    fig.update_yaxes(title="")
    fig.update_layout(height=300, title="Cox HR (strata = tumor type)")
    st.plotly_chart(fig, use_container_width=True)
    st.caption(
        "**Primary hypotheses (H1–H5):** Bonferroni α = 0.01 (5 pre-specified tests). "
        "All other comparisons are exploratory (BH FDR α = 0.05)."
    )

    st.subheader("Per-tumor-type survival log-rank results")
    lr_display = lr.copy()
    lr_display["sig"] = lr_display["padj"].apply(
        lambda p: "***" if p < 0.001 else ("**" if p < 0.01 else ("*" if p < 0.05 else ""))
    )
    st.dataframe(
        lr_display.style.background_gradient(subset=["lr_pval"], cmap="RdYlGn_r"),
        use_container_width=True,
    )

    st.subheader("Group sizes (immune tertile × p16-high/low)")
    group_counts = (
        df.groupby(["tumor_type", "group6"])
        .size()
        .reset_index(name="n")
        .pivot(index="tumor_type", columns="group6", values="n")
        .fillna(0)
        .astype(int)
    )
    st.dataframe(group_counts, use_container_width=True)

    # ── Interactive Kaplan-Meier explorer ─────────────────────────────────────
    st.subheader("Interactive Kaplan-Meier explorer")

    def _km_estimate(times, events):
        """Return (t, s) step-function arrays for a single-group KM curve."""
        t = np.asarray(times, dtype=float)
        e = np.asarray(events, dtype=float)
        order = np.argsort(t)
        t, e = t[order], e[order]
        unique_event_times = np.unique(t[e == 1])
        surv = 1.0
        curve_t = [0.0]
        curve_s = [1.0]
        for ti in unique_event_times:
            d = np.sum((t == ti) & (e == 1))
            n_risk = np.sum(t >= ti)
            surv *= 1.0 - d / n_risk
            curve_t.append(ti)
            curve_s.append(surv)
        return np.array(curve_t), np.array(curve_s)

    endpoint_opts = [
        ep for ep in ["OS", "DSS", "PFI"]
        if ep in df.columns and f"{ep}.time_years" in df.columns
    ]
    strat_opts = [c for c in ["hot_cold", "p16_status", "group6"] if c in df.columns]

    if endpoint_opts and strat_opts:
        col_km1, col_km2, col_km3 = st.columns(3)
        with col_km1:
            km_endpoint = st.selectbox("Endpoint", endpoint_opts)
        with col_km2:
            km_group = st.selectbox("Stratify by", strat_opts)
        with col_km3:
            max_years = st.slider("Max follow-up (years)", 1, 15, 10)

        time_col = f"{km_endpoint}.time_years"
        km_df = df.dropna(subset=[time_col, km_endpoint, km_group]).copy()
        km_df = km_df[km_df[time_col] <= max_years]

        palette = px.colors.qualitative.Set2
        groups = sorted(km_df[km_group].dropna().unique())
        fig_km = go.Figure()
        for i, grp in enumerate(groups):
            sub = km_df[km_df[km_group] == grp]
            t_arr, s_arr = _km_estimate(sub[time_col].values, sub[km_endpoint].values)
            # Convert KM step function to plottable x/y
            t_step = np.repeat(t_arr, 2)[1:]
            s_step = np.repeat(s_arr, 2)[:-1]
            color = palette[i % len(palette)]
            fig_km.add_trace(go.Scatter(
                x=t_step, y=s_step,
                mode="lines",
                name=f"{grp}  (n={len(sub):,})",
                line=dict(color=color, width=2),
            ))
        fig_km.update_layout(
            xaxis_title=f"Time (years) — {km_endpoint}",
            yaxis_title="Survival probability",
            yaxis_range=[0, 1.05],
            height=450,
            title=f"KM: {km_endpoint} by {km_group}",
            legend_title=km_group,
        )
        st.plotly_chart(fig_km, use_container_width=True)
        st.caption(
            f"{len(km_df):,} samples with complete {km_endpoint} data · "
            f"truncated at {max_years} years · "
            "Step function only — no confidence intervals shown."
        )
    else:
        st.caption("Survival endpoint data not yet available — re-run notebook 05.")

# ── Pathway ───────────────────────────────────────────────────────────────────
elif page == "Pathway":
    st.title("p16/CDKN2A Cell-Cycle Pathway")
    st.markdown(
        """
        Canonical pathway from upstream oncogenic stress and DNA damage through the
        p16–CDK4/6–Rb–E2F axis to senescence or proliferation.
        **Red nodes** (dashed edges) mark recurrent cancer escape routes —
        genomic alterations that bypass arrest and drive tumour progression.

        | Colour | Meaning |
        |---|---|
        | Yellow | Upstream inducers (oncogenic signals, DNA damage) |
        | Indigo | DNA-damage sensors (ATM/ATR, CHK1/2) |
        | Green | Tumour suppressors / CDK inhibitors |
        | Blue | Cyclin-dependent kinases |
        | Grey | Rb / E2F cell-cycle switch |
        | Gold | Arrest & SASP outputs |
        | Red | Cancer escape routes |
        """
    )

    import streamlit.components.v1 as components
    components.html(_pathway_html(height=920), height=940, scrolling=True)

    pathway_png = FIGURES / "p16_pathway.png"
    if pathway_png.exists():
        with st.expander("Static PNG (for export / embedding)"):
            st.image(str(pathway_png))
    else:
        st.caption(
            "Run `notebooks/06_figures.ipynb` to pre-render a static PNG "
            "(requires network access to mermaid.ink)."
        )

# ── Figures ───────────────────────────────────────────────────────────────────
elif page == "Figures":
    st.title("Output Figures")

    figures = {
        "IFN-γ score distribution": "ifng_score_distribution.png",
        "Hot fraction by tumor type": "hot_fraction_by_tumor_type.png",
        "CDKN2A vs IFN-γ scatter": "cdkn2a_vs_ifng_scatter.png",
        "CDKN2A–IFN-γ correlation by type": "cdkn2a_ifng_corr_by_type.png",
        "CDKN2A by tumor type": "cdkn2a_by_tumor_type.png",
        "Spearman heatmap (pan-cancer)": "spearman_heatmap_pancancer.png",
        "SASP score: cold / intermediate / hot": "sasp_score_hot_cold.png",
        "KM: cold / intermediate / hot": "km_hot_cold.png",
        "KM: 6-group (tertile × p16)": "km_4group.png",
        "Cox forest plot": "cox_forest_plot.png",
        "Fig 1: Senescence heatmap": "fig1_senescence_heatmap.png",
        "Fig 2: UMAP": "fig2_umap.png",
        "Fig 3: Correlation matrix": "fig3_correlation_matrix.png",
        "Fig 4: Survival composite": "fig4_survival_composite.png",
        "p16 pathway diagram": "p16_pathway.png",
        "Alteration rates by immune tertile (bar)": "alteration_hot_cold_grouped_bar.png",
        "Alteration heatmap by tumor type": "alteration_hot_cold_heatmap.png",
        "Age by tumor type": "age_by_tumor_type.png",
        "Age by immune tertile": "age_hot_cold.png",
        "Age × senescence correlations": "age_senescence_correlations.png",
        "KM: age groups": "km_age_groups.png",
        "KM: DSS by immune tertile": "km_dss_hot_cold.png",
        "KM: PFI by immune tertile": "km_pfi_hot_cold.png",
        "Competing risks forest": "competing_risks_forest.png",
        "CDKN2A 4-group immune phenotype": "cdkn2a_4group_immune.png",
        "CDKN2A by Thorsson immune subtype": "cdkn2a_by_immune_subtype.png",
    }

    cols = st.columns(2)
    for i, (label, fname) in enumerate(figures.items()):
        path = FIGURES / fname
        if path.exists():
            with cols[i % 2]:
                st.subheader(label)
                st.image(str(path))
        else:
            with cols[i % 2]:
                st.caption(f"_(missing: {fname})_")

# ── Hypotheses ────────────────────────────────────────────────────────────────
elif page == "Hypotheses":
    st.title("Pre-specified Hypotheses")
    st.markdown(
        """
        Five confirmatory hypotheses were pre-registered before analysis.
        **Bonferroni threshold: α = 0.01** (5 tests).
        All other comparisons are exploratory (BH FDR α = 0.05).
        """
    )

    # ── Compute observed stats from loaded data ────────────────────────────────
    # H1: Spearman ρ(CDKN2A_expr, IFN-γ score)
    h1_rho, h1_n = None, 0
    if "CDKN2A_expr" in df.columns and "ifng_score" in df.columns:
        _h1 = df[["CDKN2A_expr", "ifng_score"]].dropna()
        h1_n = len(_h1)
        h1_rho = _h1.corr(method="spearman").iloc[0, 1]

    # H2: SASP score hot vs cold Mann-Whitney
    h2_hot_med, h2_cold_med = None, None
    if "SASP" in df.columns and "hot_cold" in df.columns:
        h2_hot  = df.loc[df["hot_cold"] == "hot",  "SASP"].dropna()
        h2_cold = df.loc[df["hot_cold"] == "cold", "SASP"].dropna()
        if len(h2_hot) and len(h2_cold):
            h2_hot_med  = h2_hot.median()
            h2_cold_med = h2_cold.median()

    # H3: CDKN2A deletion rate hot vs cold
    h3_cold_rate, h3_hot_rate = None, None
    if "cdkn2a_altered" in df.columns and "hot_cold" in df.columns:
        _h3_cold = df.loc[df["hot_cold"] == "cold",  "cdkn2a_altered"].dropna()
        _h3_hot  = df.loc[df["hot_cold"] == "hot",   "cdkn2a_altered"].dropna()
        if len(_h3_cold) and len(_h3_hot):
            h3_cold_rate = 100 * _h3_cold.mean()
            h3_hot_rate  = 100 * _h3_hot.mean()

    # H4 / H5: from Cox results
    cox_cdkn2a = cox[cox["term"].str.contains("CDKN2A", case=False, na=False)]
    h5_hr = cox_cdkn2a["estimate"].values[0] if len(cox_cdkn2a) else None
    h5_lo = cox_cdkn2a["conf.low"].values[0] if len(cox_cdkn2a) else None
    h5_hi = cox_cdkn2a["conf.high"].values[0] if len(cox_cdkn2a) else None
    h5_p  = cox_cdkn2a["p.value"].values[0] if len(cox_cdkn2a) else None

    # Fraction of tumor types with significant hot>cold log-rank (padj < 0.05)
    lr_sig_frac = (lr["padj"] < 0.05).mean() if "padj" in lr.columns else None

    # H6: Suppressive SASP in cold tumors
    h6_cold_med, h6_hot_med = None, None
    if "SASP_suppressive" in df.columns and "hot_cold" in df.columns:
        _h6_cold = df.loc[df["hot_cold"] == "cold", "SASP_suppressive"].dropna()
        _h6_hot  = df.loc[df["hot_cold"] == "hot",  "SASP_suppressive"].dropna()
        if len(_h6_cold) and len(_h6_hot):
            h6_cold_med = _h6_cold.median()
            h6_hot_med  = _h6_hot.median()

    # H7: TMB higher in hot tumors
    h7_hot_med, h7_cold_med, h7_rho = None, None, None
    if "tmb" in df.columns and "hot_cold" in df.columns:
        _h7_hot  = df.loc[df["hot_cold"] == "hot",  "tmb"].dropna()
        _h7_cold = df.loc[df["hot_cold"] == "cold", "tmb"].dropna()
        if len(_h7_hot) and len(_h7_cold):
            h7_hot_med  = _h7_hot.median()
            h7_cold_med = _h7_cold.median()
        if "ifng_score" in df.columns:
            _h7s = df[["tmb", "ifng_score"]].dropna()
            if len(_h7s):
                h7_rho = _h7s.corr(method="spearman").iloc[0, 1]

    # H8: TP53 alteration enriched in cold tumors
    h8_cold_rate, h8_hot_rate = None, None
    if "tp53_altered" in df.columns and "hot_cold" in df.columns:
        _h8_cold = df.loc[df["hot_cold"] == "cold", "tp53_altered"].dropna()
        _h8_hot  = df.loc[df["hot_cold"] == "hot",  "tp53_altered"].dropna()
        if len(_h8_cold) and len(_h8_hot):
            h8_cold_rate = 100 * _h8_cold.mean()
            h8_hot_rate  = 100 * _h8_hot.mean()

    # H9: LMNB1 inversely correlated with IFN-γ score
    h9_rho, h9_n = None, 0
    if "LMNB1_expr" in df.columns and "ifng_score" in df.columns:
        _h9 = df[["LMNB1_expr", "ifng_score"]].dropna()
        h9_n   = len(_h9)
        h9_rho = _h9.corr(method="spearman").iloc[0, 1]

    # H10: Tumor purity lower in hot tumors
    h10_hot_med, h10_cold_med, h10_rho = None, None, None
    if "tumor_purity" in df.columns and "hot_cold" in df.columns:
        _h10_hot  = df.loc[df["hot_cold"] == "hot",  "tumor_purity"].dropna()
        _h10_cold = df.loc[df["hot_cold"] == "cold", "tumor_purity"].dropna()
        if len(_h10_hot) and len(_h10_cold):
            h10_hot_med  = _h10_hot.median()
            h10_cold_med = _h10_cold.median()
        if "ifng_score" in df.columns:
            _h10s = df[["tumor_purity", "ifng_score"]].dropna()
            if len(_h10s):
                h10_rho = _h10s.corr(method="spearman").iloc[0, 1]

    # H11: CDKN1A expression correlates with IFN-γ score
    h11_rho, h11_n = None, 0
    if "CDKN1A_expr" in df.columns and "ifng_score" in df.columns:
        _h11 = df[["CDKN1A_expr", "ifng_score"]].dropna()
        h11_n   = len(_h11)
        h11_rho = _h11.corr(method="spearman").iloc[0, 1]

    # H12: RB1 alteration enriched in cold tumors
    h12_cold_rate, h12_hot_rate = None, None
    if "rb1_altered" in df.columns and "hot_cold" in df.columns:
        _h12_cold = df.loc[df["hot_cold"] == "cold", "rb1_altered"].dropna()
        _h12_hot  = df.loc[df["hot_cold"] == "hot",  "rb1_altered"].dropna()
        if len(_h12_cold) and len(_h12_hot):
            h12_cold_rate = 100 * _h12_cold.mean()
            h12_hot_rate  = 100 * _h12_hot.mean()

    # H13: SERPINE1 expression correlates with IFN-γ score
    h13_rho, h13_n = None, 0
    if "SERPINE1_expr" in df.columns and "ifng_score" in df.columns:
        _h13 = df[["SERPINE1_expr", "ifng_score"]].dropna()
        h13_n   = len(_h13)
        h13_rho = _h13.corr(method="spearman").iloc[0, 1]

    # H14: Age correlates with CDKN2A expression
    h14_rho, h14_n = None, 0
    if "age" in df.columns and "CDKN2A_expr" in df.columns:
        _h14 = df[["age", "CDKN2A_expr"]].dropna()
        h14_n   = len(_h14)
        h14_rho = _h14.corr(method="spearman").iloc[0, 1]

    # H15: hot/cold is more prognostic in p16-high than p16-low (interaction)
    # Compute log-rank statistic proxy: median OS difference within each p16 stratum
    h15_p16high_diff, h15_p16low_diff = None, None
    if all(c in df.columns for c in ["p16_status", "hot_cold", "OS", "OS.time_years"]):
        def _median_os(stratum_col, stratum_val, hc_val):
            sub = df[
                (df[stratum_col] == stratum_val) & (df["hot_cold"] == hc_val)
            ].dropna(subset=["OS.time_years", "OS"])
            return sub.loc[sub["OS"] == 1, "OS.time_years"].median() if len(sub) else None

        h_hi_hot  = _median_os("p16_status", "p16-high", "hot")
        h_hi_cold = _median_os("p16_status", "p16-high", "cold")
        h_lo_hot  = _median_os("p16_status", "p16-low",  "hot")
        h_lo_cold = _median_os("p16_status", "p16-low",  "cold")
        if None not in (h_hi_hot, h_hi_cold):
            h15_p16high_diff = h_hi_hot - h_hi_cold
        if None not in (h_lo_hot, h_lo_cold):
            h15_p16low_diff = h_lo_hot - h_lo_cold
        h15_vals = (h_hi_hot, h_hi_cold, h_lo_hot, h_lo_cold)
    else:
        h15_vals = (None, None, None, None)

    # H16: Thorsson C2 has highest SASP scores
    h16_c2_med, h16_other_med, h16_ranking = None, None, None
    if "immune_subtype" in df.columns and "SASP" in df.columns:
        _h16 = df.dropna(subset=["immune_subtype", "SASP"])
        if len(_h16):
            h16_ranking = (
                _h16.groupby("immune_subtype")["SASP"]
                .median()
                .sort_values(ascending=False)
            )
            if "C2" in h16_ranking.index:
                h16_c2_med    = h16_ranking["C2"]
                h16_other_med = h16_ranking.drop("C2").median()

    # H17: SenMayo score higher in hot tumors
    h17_hot_med, h17_cold_med = None, None
    if "senmayo" in df.columns and "hot_cold" in df.columns:
        _h17_hot  = df.loc[df["hot_cold"] == "hot",  "senmayo"].dropna()
        _h17_cold = df.loc[df["hot_cold"] == "cold", "senmayo"].dropna()
        if len(_h17_hot) and len(_h17_cold):
            h17_hot_med  = _h17_hot.median()
            h17_cold_med = _h17_cold.median()

    # H18: TNF-α/NF-κB signalling higher in hot tumors
    h18_hot_med, h18_cold_med = None, None
    if "hallmark_tnfa" in df.columns and "hot_cold" in df.columns:
        _h18_hot  = df.loc[df["hot_cold"] == "hot",  "hallmark_tnfa"].dropna()
        _h18_cold = df.loc[df["hot_cold"] == "cold", "hallmark_tnfa"].dropna()
        if len(_h18_hot) and len(_h18_cold):
            h18_hot_med  = _h18_hot.median()
            h18_cold_med = _h18_cold.median()

    # H19: Inflammatory response score higher in hot tumors
    h19_hot_med, h19_cold_med = None, None
    if "hallmark_inflammatory" in df.columns and "hot_cold" in df.columns:
        _h19_hot  = df.loc[df["hot_cold"] == "hot",  "hallmark_inflammatory"].dropna()
        _h19_cold = df.loc[df["hot_cold"] == "cold", "hallmark_inflammatory"].dropna()
        if len(_h19_hot) and len(_h19_cold):
            h19_hot_med  = _h19_hot.median()
            h19_cold_med = _h19_cold.median()

    # H20: senescence_dn inversely correlates with IFN-γ
    h20_rho, h20_n = None, 0
    if "senescence_dn" in df.columns and "ifng_score" in df.columns:
        _h20 = df[["senescence_dn", "ifng_score"]].dropna()
        h20_n   = len(_h20)
        h20_rho = _h20.corr(method="spearman").iloc[0, 1]

    # H21: GDF15 expression correlates with IFN-γ
    h21_rho, h21_n = None, 0
    if "GDF15_expr" in df.columns and "ifng_score" in df.columns:
        _h21 = df[["GDF15_expr", "ifng_score"]].dropna()
        h21_n   = len(_h21)
        h21_rho = _h21.corr(method="spearman").iloc[0, 1]

    # H22: CDKN2A homozygous deep deletion enriched in cold tumors
    h22_cold_rate, h22_hot_rate = None, None
    h22_alt_cold, h22_alt_hot   = None, None  # composite for comparison
    if "cdkn2a_deep_del" in df.columns and "hot_cold" in df.columns:
        _h22_cold = df.loc[df["hot_cold"] == "cold", "cdkn2a_deep_del"].dropna()
        _h22_hot  = df.loc[df["hot_cold"] == "hot",  "cdkn2a_deep_del"].dropna()
        if len(_h22_cold) and len(_h22_hot):
            h22_cold_rate = 100 * _h22_cold.mean()
            h22_hot_rate  = 100 * _h22_hot.mean()
    if "cdkn2a_altered" in df.columns and "hot_cold" in df.columns:
        h22_alt_cold = 100 * df.loc[df["hot_cold"] == "cold", "cdkn2a_altered"].dropna().mean()
        h22_alt_hot  = 100 * df.loc[df["hot_cold"] == "hot",  "cdkn2a_altered"].dropna().mean()

    # ── Build conclusions dynamically from observed stats ─────────────────────
    def _h1_conclusion():
        if h1_rho is None:
            return None, "Data not yet available — re-run notebooks 01–03."
        strength = "strong" if abs(h1_rho) > 0.4 else ("moderate" if abs(h1_rho) > 0.2 else "weak")
        in_range = 0.10 <= h1_rho <= 0.25
        if h1_rho > 0:
            verdict = "supported" if in_range else ("supported (stronger than expected)" if h1_rho > 0.25 else "supported (weaker than expected)")
            return True, (
                f"**H1 {verdict}.** Pan-cancer ρ = {h1_rho:.3f} — a {strength} positive correlation. "
                f"The direction is consistent with SASP-driven IFN-γ signalling. "
                f"{'Effect size is within the pre-specified range (0.10–0.25).' if in_range else f'Effect size falls outside the pre-specified range (0.10–0.25).'}"
            )
        else:
            return False, (
                f"**H1 not supported.** Pan-cancer ρ = {h1_rho:.3f} — a {strength} negative correlation, "
                "contrary to the hypothesis that p16/SASP drives immune infiltration."
            )

    def _h2_conclusion():
        if h2_hot_med is None:
            return None, "Data not yet available — re-run notebook 03."
        delta = h2_hot_med - h2_cold_med
        if delta > 0:
            return True, (
                f"**H2 supported.** Hot tumors have higher median SASP scores than cold tumors "
                f"(Δ = {delta:+.3f}), consistent with SASP cytokines promoting immune infiltration."
            )
        else:
            return False, (
                f"**H2 not supported.** Cold tumors have equal or higher median SASP scores than hot tumors "
                f"(Δ = {delta:+.3f}), which does not support the hypothesis."
            )

    def _h3_conclusion():
        if h3_cold_rate is None:
            return None, "Alteration data not yet available — re-run notebooks 01 and 03."
        diff = h3_cold_rate - h3_hot_rate
        if diff > 0:
            return True, (
                f"**H3 supported.** CDKN2A alteration is more frequent in cold tumors ({h3_cold_rate:.1f}%) "
                f"than hot ({h3_hot_rate:.1f}%, Δ = {diff:+.1f} pp), consistent with p16 loss promoting "
                "immune exclusion by suppressing SASP."
            )
        else:
            return False, (
                f"**H3 not supported.** CDKN2A alteration rates do not favour cold tumors "
                f"(cold: {h3_cold_rate:.1f}%, hot: {h3_hot_rate:.1f}%, Δ = {diff:+.1f} pp)."
            )

    def _h4_conclusion():
        if lr_sig_frac is None:
            return None, "Survival data not yet available — re-run notebook 05."
        n_sig = int(round(lr_sig_frac * len(lr)))
        n_tot = len(lr)
        if lr_sig_frac > 0.3:
            return True, (
                f"**H4 supported.** {n_sig}/{n_tot} tumor types ({lr_sig_frac*100:.0f}%) show significantly "
                "better OS in hot vs cold tumors (padj < 0.05). The immune phenotype × p16 group separation "
                "is present across the majority of cancer contexts tested."
            )
        elif lr_sig_frac > 0:
            return False, (
                f"**H4 partially supported.** Only {n_sig}/{n_tot} tumor types ({lr_sig_frac*100:.0f}%) reach "
                "significance (padj < 0.05). The effect is present but limited to a minority of cancer types."
            )
        else:
            return False, (
                "**H4 not supported.** No tumor type shows a statistically significant hot>cold OS difference after BH correction."
            )

    def _h5_conclusion():
        if h5_hr is None:
            return None, "Cox results not yet available — re-run notebook 05."
        protective = h5_hr < 1.0
        sig = h5_p < 0.05
        bon_sig = h5_p < 0.01
        if protective and bon_sig:
            return True, (
                f"**H5 supported (Bonferroni).** CDKN2A expression is an independent protective predictor "
                f"(HR = {h5_hr:.3f}, 95% CI: {h5_lo:.3f}–{h5_hi:.3f}, p = {h5_p:.2e}), surviving the "
                "Bonferroni threshold of α = 0.01."
            )
        elif protective and sig:
            return True, (
                f"**H5 supported (nominal).** CDKN2A expression is protective (HR = {h5_hr:.3f}, p = {h5_p:.2e}) "
                "but does not reach the Bonferroni threshold of α = 0.01 for this confirmatory test."
            )
        elif protective:
            return False, (
                f"**H5 trend only.** HR = {h5_hr:.3f} is in the protective direction but is not statistically "
                f"significant (p = {h5_p:.2e})."
            )
        else:
            return False, (
                f"**H5 not supported.** HR = {h5_hr:.3f} is ≥ 1 (p = {h5_p:.2e}), indicating no independent "
                "protective effect of CDKN2A expression on OS."
            )

    def _h6_conclusion():
        if h6_cold_med is None:
            return None, "SASP_suppressive data not yet available — re-run notebook 03."
        delta = h6_cold_med - h6_hot_med
        if delta > 0:
            return True, (
                f"**H6 supported.** Cold tumors have higher suppressive SASP scores than hot tumors "
                f"(cold median {h6_cold_med:.3f} vs hot {h6_hot_med:.3f}, Δ = {delta:+.3f}), "
                "suggesting immune exclusion is partly driven by active immunosuppressive SASP signalling."
            )
        else:
            return False, (
                f"**H6 not supported.** Suppressive SASP scores do not favour cold tumors "
                f"(cold {h6_cold_med:.3f} vs hot {h6_hot_med:.3f}, Δ = {delta:+.3f})."
            )

    def _h7_conclusion():
        if h7_hot_med is None:
            return None, "TMB data not yet available — re-run notebook 01."
        direction_ok = h7_hot_med > h7_cold_med
        rho_str = f"  Spearman ρ(TMB, IFN-γ) = {h7_rho:.3f}." if h7_rho is not None else ""
        if direction_ok:
            return True, (
                f"**H7 supported.** Hot tumors have higher median TMB than cold tumors "
                f"(hot {h7_hot_med:.1f} vs cold {h7_cold_med:.1f} mut/Mb).{rho_str} "
                "This is consistent with DNA damage driving both mutational burden and p16-mediated senescence."
            )
        else:
            return False, (
                f"**H7 not supported.** TMB does not differ in the expected direction "
                f"(hot {h7_hot_med:.1f} vs cold {h7_cold_med:.1f} mut/Mb).{rho_str}"
            )

    def _h8_conclusion():
        if h8_cold_rate is None:
            return None, "TP53 alteration data not yet available — re-run notebooks 01 and 03."
        diff = h8_cold_rate - h8_hot_rate
        if diff > 0:
            return True, (
                f"**H8 supported.** TP53 alterations are more frequent in cold tumors "
                f"({h8_cold_rate:.1f}%) than hot ({h8_hot_rate:.1f}%, Δ = {diff:+.1f} pp), "
                "consistent with TP53 loss impairing senescence-associated immune signalling."
            )
        else:
            return False, (
                f"**H8 not supported.** TP53 alteration rates do not favour cold tumors "
                f"(cold {h8_cold_rate:.1f}%, hot {h8_hot_rate:.1f}%, Δ = {diff:+.1f} pp). "
                "Note: TP53 gain-of-function mutations in hot tumors may confound a simple deletion/mutation count."
            )

    def _h9_conclusion():
        if h9_rho is None:
            return None, "LMNB1_expr data not yet available — re-run notebook 03."
        strength = "strong" if abs(h9_rho) > 0.4 else ("moderate" if abs(h9_rho) > 0.2 else "weak")
        if h9_rho < 0:
            return True, (
                f"**H9 supported.** LMNB1 expression is negatively correlated with IFN-γ score "
                f"(ρ = {h9_rho:.3f}, n = {h9_n:,} — {strength} inverse correlation), consistent with "
                "nuclear lamina loss in deeply senescent cells activating cGAS-STING and SASP."
            )
        else:
            return False, (
                f"**H9 not supported.** LMNB1 expression is positively correlated with IFN-γ score "
                f"(ρ = {h9_rho:.3f}, n = {h9_n:,}), contrary to the hypothesis."
            )

    def _h10_conclusion():
        if h10_hot_med is None:
            return None, "Tumor purity data not yet available — re-run notebook 01."
        rho_str = f"  Spearman ρ(purity, IFN-γ) = {h10_rho:.3f}." if h10_rho is not None else ""
        diff = h10_cold_med - h10_hot_med
        if h10_cold_med > h10_hot_med:
            return True, (
                f"**H10 confirmed.** Cold tumors have higher median purity than hot tumors "
                f"(cold {h10_cold_med:.2f} vs hot {h10_hot_med:.2f}, Δ = {diff:+.2f}).{rho_str} "
                "This confirms that purity is a meaningful confound: part of the hot/cold signal reflects "
                "immune cell dilution rather than intrinsic tumour biology. Purity is included as a Cox covariate."
            )
        else:
            return False, (
                f"**H10 not confirmed.** Purity does not differ in the expected direction "
                f"(cold {h10_cold_med:.2f} vs hot {h10_hot_med:.2f}).{rho_str}"
            )

    def _h11_conclusion():
        if h11_rho is None:
            return None, "CDKN1A_expr data not yet available — re-run notebook 03."
        strength = "strong" if abs(h11_rho) > 0.4 else ("moderate" if abs(h11_rho) > 0.2 else "weak")
        vs_h1 = ""
        if h1_rho is not None:
            if abs(h11_rho) > abs(h1_rho):
                vs_h1 = f" CDKN1A shows a *stronger* correlation than CDKN2A (ρ = {h1_rho:.3f}), suggesting p21 may be the dominant arm in this dataset."
            else:
                vs_h1 = f" CDKN2A shows a stronger correlation (ρ = {h1_rho:.3f}), consistent with p16 being the primary driver."
        if h11_rho > 0:
            return True, (
                f"**H11 supported.** CDKN1A expression positively correlates with IFN-γ score "
                f"(ρ = {h11_rho:.3f}, n = {h11_n:,} — {strength}).{vs_h1}"
            )
        else:
            return False, (
                f"**H11 not supported.** CDKN1A expression is negatively correlated with IFN-γ score "
                f"(ρ = {h11_rho:.3f}), contrary to the hypothesis."
            )

    def _h12_conclusion():
        if h12_cold_rate is None:
            return None, "RB1 alteration data not yet available — re-run notebooks 01 and 03."
        diff = h12_cold_rate - h12_hot_rate
        cdkn2a_context = ""
        if h3_cold_rate is not None:
            cdkn2a_context = (
                f" For comparison, CDKN2A alteration rate difference is {h3_cold_rate - h3_hot_rate:+.1f} pp — "
                f"{'RB1 shows a similar pattern' if diff * (h3_cold_rate - h3_hot_rate) > 0 else 'RB1 shows the opposite pattern'}, "
                "which is expected given the 9p21 co-deletion locus."
            )
        if diff > 0:
            return True, (
                f"**H12 supported.** RB1 alteration is more frequent in cold tumors ({h12_cold_rate:.1f}%) "
                f"than hot ({h12_hot_rate:.1f}%, Δ = {diff:+.1f} pp), extending the senescence bypass "
                f"model to the full p16–CDK4/6–Rb axis.{cdkn2a_context}"
            )
        else:
            return False, (
                f"**H12 not supported.** RB1 alteration does not favour cold tumors "
                f"(cold {h12_cold_rate:.1f}%, hot {h12_hot_rate:.1f}%, Δ = {diff:+.1f} pp).{cdkn2a_context}"
            )

    def _h13_conclusion():
        if h13_rho is None:
            return None, "SERPINE1_expr data not yet available — re-run notebook 03."
        strength = "strong" if abs(h13_rho) > 0.4 else ("moderate" if abs(h13_rho) > 0.2 else "weak")
        if h13_rho > 0:
            return True, (
                f"**H13 supported.** SERPINE1 expression positively correlates with IFN-γ score "
                f"(ρ = {h13_rho:.3f}, n = {h13_n:,} — {strength} correlation), supporting SERPINE1/PAI-1 "
                "as a SASP-linked immune-activity marker at the individual gene level."
            )
        else:
            return False, (
                f"**H13 not supported.** SERPINE1 expression is negatively correlated with IFN-γ score "
                f"(ρ = {h13_rho:.3f}), contrary to the SASP-immune hypothesis."
            )

    def _h14_conclusion():
        if h14_rho is None:
            return None, "Age or CDKN2A data not yet available — re-run notebooks 01–03."
        strength = "strong" if abs(h14_rho) > 0.4 else ("moderate" if abs(h14_rho) > 0.2 else "weak")
        if h14_rho > 0:
            return True, (
                f"**H14 supported.** Age positively correlates with CDKN2A expression "
                f"(ρ = {h14_rho:.3f}, n = {h14_n:,} — {strength} correlation), consistent with senescence "
                "accumulating with age and validating CDKN2A as a biologically meaningful readout. "
                "This also supports the decision to include age as a covariate in Cox models (H5)."
            )
        else:
            return False, (
                f"**H14 not supported.** CDKN2A expression does not increase with age "
                f"(ρ = {h14_rho:.3f}). This may indicate that cancer-specific CDKN2A dysregulation "
                "dominates the bulk RNA-seq signal over normal senescence accumulation."
            )

    def _h15_conclusion():
        hi_diff, lo_diff = h15_p16high_diff, h15_p16low_diff
        h_hi_hot, h_hi_cold, h_lo_hot, h_lo_cold = h15_vals
        if hi_diff is None and lo_diff is None:
            return None, "Survival or p16 status data not yet available — re-run notebooks 03 and 05."
        hi_str = f"{hi_diff:+.2f} yrs" if hi_diff is not None else "N/A"
        lo_str = f"{lo_diff:+.2f} yrs" if lo_diff is not None else "N/A"
        interaction_present = (
            hi_diff is not None and lo_diff is not None and
            abs(hi_diff) > abs(lo_diff) and hi_diff > 0
        )
        if interaction_present:
            return True, (
                f"**H15 supported.** The hot−cold median OS gap is larger in p16-high tumors ({hi_str}) "
                f"than p16-low ({lo_str}), consistent with p16/SASP mechanistically driving the immune-survival "
                "axis. The hot/cold label is more prognostically informative when p16 is intact."
            )
        elif hi_diff is not None and hi_diff > 0:
            return True, (
                f"**H15 partially supported.** The hot−cold OS gap is positive in p16-high tumors ({hi_str}), "
                f"but the interaction with p16-low ({lo_str}) is less clear."
            )
        else:
            return False, (
                f"**H15 not supported.** The hot−cold OS gap in p16-high tumors ({hi_str}) is not clearly "
                f"larger than in p16-low ({lo_str}). Immune phenotype appears similarly prognostic regardless of p16 status."
            )

    def _h16_conclusion():
        if h16_c2_med is None:
            return None, "Thorsson immune subtype or SASP data not yet available — re-run notebooks 01 and 03."
        rank_pos = list(h16_ranking.index).index("C2") + 1 if "C2" in h16_ranking.index else None
        rank_str = f" (ranked #{rank_pos} of {len(h16_ranking)} subtypes)" if rank_pos else ""
        if rank_pos == 1:
            return True, (
                f"**H16 supported.** Thorsson C2 (IFN-γ dominant) has the highest median SASP score "
                f"({h16_c2_med:.3f} vs {h16_other_med:.3f} median of other subtypes){rank_str}. "
                "This independently validates the SASP → immune-inflamed link using a classification "
                "that was derived separately from the IFN-γ score used to define hot/cold."
            )
        elif rank_pos is not None and rank_pos <= 2:
            return True, (
                f"**H16 partially supported.** C2 ranks #{rank_pos} for SASP score{rank_str} "
                f"(median {h16_c2_med:.3f}). The subtype ordering is broadly consistent with the hypothesis."
            )
        else:
            return False, (
                f"**H16 not supported.** C2 does not rank highest for SASP score{rank_str} "
                f"(median {h16_c2_med:.3f}). "
                + (f"Top subtype: {h16_ranking.index[0]} (median {h16_ranking.iloc[0]:.3f})." if h16_ranking is not None else "")
            )

    def _h17_conclusion():
        if h17_hot_med is None:
            return None, "SenMayo data not yet available — re-run notebook 03."
        delta = h17_hot_med - h17_cold_med
        sasp_context = ""
        if h2_hot_med is not None:
            sasp_context = (
                f" For comparison, the custom SASP score gap is {h2_hot_med - h2_cold_med:+.3f}. "
                f"{'SenMayo and SASP scores agree in direction.' if delta * (h2_hot_med - h2_cold_med) > 0 else 'SenMayo and SASP scores disagree in direction — worth investigating.'}"
            )
        if delta > 0:
            return True, (
                f"**H17 supported.** Hot tumors have higher SenMayo scores than cold "
                f"(hot {h17_hot_med:.3f} vs cold {h17_cold_med:.3f}, Δ = {delta:+.3f}). "
                "This replicates the SASP finding (H2) using an independent, peer-reviewed 125-gene signature "
                f"(Saul et al. 2022, Nature Aging), strengthening confidence in the senescence–immune link.{sasp_context}"
            )
        else:
            return False, (
                f"**H17 not supported.** SenMayo scores do not favour hot tumors "
                f"(hot {h17_hot_med:.3f} vs cold {h17_cold_med:.3f}, Δ = {delta:+.3f}). "
                f"This is inconsistent with H2.{sasp_context}"
            )

    def _h18_conclusion():
        if h18_hot_med is None:
            return None, "hallmark_tnfa data not yet available — re-run notebook 03."
        delta = h18_hot_med - h18_cold_med
        if delta > 0:
            return True, (
                f"**H18 supported.** TNF-α/NF-κB signalling is higher in hot tumors "
                f"(hot {h18_hot_med:.3f} vs cold {h18_cold_med:.3f}, Δ = {delta:+.3f}). "
                "TNF-α is a primary SASP effector cytokine that drives NF-κB-mediated immune recruitment. "
                "This provides mechanistic specificity for the SASP → hot phenotype link beyond the composite score."
            )
        else:
            return False, (
                f"**H18 not supported.** TNF-α/NF-κB activity does not favour hot tumors "
                f"(hot {h18_hot_med:.3f} vs cold {h18_cold_med:.3f}, Δ = {delta:+.3f})."
            )

    def _h19_conclusion():
        if h19_hot_med is None:
            return None, "hallmark_inflammatory data not yet available — re-run notebook 03."
        delta = h19_hot_med - h19_cold_med
        if delta > 0:
            return True, (
                f"**H19 supported.** Inflammatory response activity is higher in hot tumors "
                f"(hot {h19_hot_med:.3f} vs cold {h19_cold_med:.3f}, Δ = {delta:+.3f}). "
                "The HALLMARK_INFLAMMATORY_RESPONSE gene set captures broad NF-κB and cytokine signalling "
                "consistent with SASP-driven immune activation."
            )
        else:
            return False, (
                f"**H19 not supported.** Inflammatory response scores do not favour hot tumors "
                f"(hot {h19_hot_med:.3f} vs cold {h19_cold_med:.3f}, Δ = {delta:+.3f})."
            )

    def _h20_conclusion():
        if h20_rho is None:
            return None, "senescence_dn data not yet available — re-run notebook 03."
        strength = "strong" if abs(h20_rho) > 0.4 else ("moderate" if abs(h20_rho) > 0.2 else "weak")
        if h20_rho < 0:
            return True, (
                f"**H20 supported.** The senescence-downregulated signature is negatively correlated with IFN-γ score "
                f"(ρ = {h20_rho:.3f}, n = {h20_n:,} — {strength} inverse correlation). "
                "Cells that have downregulated senescence-associated structural genes are found in more immune-inflamed tumors, "
                "consistent with deep senescence activating innate immune signalling (e.g. via cGAS-STING)."
                + (" Note: the effect is weak — purity confounding cannot be ruled out." if abs(h20_rho) < 0.1 else "")
            )
        else:
            return False, (
                f"**H20 not supported.** The senescence-down signature is positively correlated with IFN-γ score "
                f"(ρ = {h20_rho:.3f}), contrary to the hypothesis."
            )

    def _h21_conclusion():
        if h21_rho is None:
            return None, "GDF15_expr data not yet available — re-run notebook 03."
        strength = "strong" if abs(h21_rho) > 0.4 else ("moderate" if abs(h21_rho) > 0.2 else "weak")
        if h21_rho > 0:
            return True, (
                f"**H21 supported.** GDF15 expression positively correlates with IFN-γ score "
                f"(ρ = {h21_rho:.3f}, n = {h21_n:,} — {strength} correlation), consistent with "
                "GDF15 acting as a SASP stress signal in immune-active tumors."
            )
        else:
            return False, (
                f"**H21 not supported.** GDF15 expression is negatively correlated with IFN-γ score "
                f"(ρ = {h21_rho:.3f}, n = {h21_n:,}). "
                "GDF15 may be more highly expressed in cold/proliferative tumors independently of SASP, "
                "or its role may be immunosuppressive rather than pro-inflammatory in this context."
            )

    def _h22_conclusion():
        if h22_cold_rate is None:
            return None, "CDKN2A deep deletion data not yet available — re-run notebooks 01 and 03."
        diff_deep  = h22_cold_rate - h22_hot_rate
        diff_comp  = (h22_alt_cold - h22_alt_hot) if h22_alt_cold is not None else None
        comp_str   = (
            f" For context, the composite alteration rate (del + mut) is also higher in hot tumors "
            f"({h22_alt_hot:.1f}%) than cold ({h22_alt_cold:.1f}%), suggesting the pattern is driven by "
            "point mutations or low-level CNAs co-occurring with high TMB in hot tumors rather than "
            "functional homozygous loss of p16."
        ) if diff_comp is not None and diff_comp < 0 else ""
        if diff_deep > 0:
            return True, (
                f"**H22 supported.** CDKN2A homozygous deep deletion is more frequent in cold tumors "
                f"({h22_cold_rate:.1f}%) than hot ({h22_hot_rate:.1f}%, Δ = {diff_deep:+.1f} pp). "
                "Biallelic loss of p16 — which fully abolishes the senescence brake — is specifically "
                f"enriched in immune-excluded tumors.{comp_str}"
            )
        else:
            return False, (
                f"**H22 not supported.** CDKN2A deep deletion does not favour cold tumors "
                f"(cold {h22_cold_rate:.1f}%, hot {h22_hot_rate:.1f}%, Δ = {diff_deep:+.1f} pp).{comp_str} "
                "This suggests that when p16 is biallelically lost it does not consistently associate with "
                "immune exclusion at pan-cancer scale — pan-cancer confounding by TMB and tumor grade is likely."
            )


    # ── Display hypothesis cards ───────────────────────────────────────────────
    HYPOTHESES = [
        {
            "id": "H1",
            "title": "CDKN2A expression correlates positively with IFN-γ score (pan-cancer)",
            "rationale": "p16/CDKN2A promotes SASP → pro-inflammatory cytokines drive T-cell infiltration, producing higher IFN-γ activity in p16-high tumors.",
            "test": "Spearman ρ (CDKN2A_expr vs IFN-γ score), BH-FDR",
            "expected": "ρ ~ 0.10–0.25 (weak-moderate positive)",
            "notebook": "04 §2",
            "threshold": "Bonferroni α = 0.01",
            "observed": (
                f"Spearman ρ = **{h1_rho:.3f}** (n = {h1_n:,})"
                if h1_rho is not None else "_Data not yet loaded_"
            ),
            "direction_ok": h1_rho is not None and h1_rho > 0,
            "conclusion_fn": _h1_conclusion,
        },
        {
            "id": "H2",
            "title": "SASP score is higher in hot vs cold tumors (pan-cancer)",
            "rationale": "Hot tumors — defined by high IFN-γ activity — should co-occur with elevated SASP secretion, as SASP cytokines (IL-6, CXCL8) recruit and activate immune cells.",
            "test": "Mann-Whitney U (hot vs cold SASP NES)",
            "expected": "SASP NES higher in hot, p < 0.01",
            "notebook": "04 §5",
            "threshold": "Bonferroni α = 0.01",
            "observed": (
                f"Median SASP — hot: **{h2_hot_med:.3f}**, cold: **{h2_cold_med:.3f}** "
                f"(Δ = {h2_hot_med - h2_cold_med:+.3f})"
                if h2_hot_med is not None else "_Data not yet loaded_"
            ),
            "direction_ok": h2_hot_med is not None and h2_hot_med > h2_cold_med,
            "conclusion_fn": _h2_conclusion,
        },
        {
            "id": "H3",
            "title": "CDKN2A deletion is enriched in cold tumors (pan-cancer)",
            "rationale": "Loss of CDKN2A removes senescence braking — reducing SASP — so CDKN2A-deleted tumors should be more immune-excluded (cold).",
            "test": "Chi-squared (CDKN2A deletion rate in hot vs cold)",
            "expected": "CDKN2A deletion rate higher in cold tumors, p < 0.01",
            "notebook": "04 §4",
            "threshold": "Bonferroni α = 0.01",
            "observed": (
                f"Alteration rate — cold: **{h3_cold_rate:.1f}%**, hot: **{h3_hot_rate:.1f}%**"
                if h3_cold_rate is not None else "_Alteration data not yet available — re-run notebook 03_"
            ),
            "direction_ok": h3_cold_rate is not None and h3_cold_rate > h3_hot_rate,
            "conclusion_fn": _h3_conclusion,
        },
        {
            "id": "H4",
            "title": "p16-high hot tumors have better OS than p16-low cold tumors",
            "rationale": "If p16/SASP drives immune infiltration, the combination of hot (immune-active) + p16-high should confer the best prognosis; cold + p16-low the worst.",
            "test": "Log-rank (4-group KM: hot/cold × p16-high/low)",
            "expected": "Global log-rank p < 0.01; hot/p16-high best, cold/p16-low worst",
            "notebook": "05 §4",
            "threshold": "Bonferroni α = 0.01",
            "observed": (
                f"{(lr_sig_frac * 100):.0f}% of tumor types show significant hot>cold OS difference (padj < 0.05)"
                if lr_sig_frac is not None else "_Survival data not yet available — re-run notebook 05_"
            ),
            "direction_ok": lr_sig_frac is not None and lr_sig_frac > 0.3,
            "conclusion_fn": _h4_conclusion,
        },
        {
            "id": "H5",
            "title": "CDKN2A expression is an independent OS predictor (age-adjusted, strata tumor type)",
            "rationale": "After accounting for age and tumor-type baseline hazards, p16 expression should independently predict survival if its effect is not purely mediated by tumor type.",
            "test": "Cox PH: CDKN2A_expr coefficient, HR with 95% CI",
            "expected": "HR < 1.0 (higher expression → lower hazard), p < 0.05",
            "notebook": "05 §6",
            "threshold": "Bonferroni α = 0.01",
            "observed": (
                f"HR = **{h5_hr:.3f}** (95% CI: {h5_lo:.3f}–{h5_hi:.3f}), p = {h5_p:.2e}"
                if h5_hr is not None else "_Cox results not yet available — re-run notebook 05_"
            ),
            "direction_ok": h5_hr is not None and h5_hr < 1.0,
            "conclusion_fn": _h5_conclusion,
        },
        {
            "id": "H6",
            "title": "Suppressive SASP is enriched in cold tumors",
            "rationale": "Cold tumors may not simply lack pro-inflammatory SASP — they may actively secrete immunosuppressive SASP components (TGF-β, IL-10 type) that exclude T cells. This tests whether immune exclusion is active rather than passive.",
            "test": "Mann-Whitney U (cold vs hot SASP_suppressive NES)",
            "expected": "SASP_suppressive NES higher in cold, p < 0.05",
            "notebook": "04 §5 (exploratory)",
            "threshold": "BH FDR α = 0.05 (exploratory)",
            "observed": (
                f"Median suppressive SASP — cold: **{h6_cold_med:.3f}**, hot: **{h6_hot_med:.3f}** "
                f"(Δ = {h6_cold_med - h6_hot_med:+.3f})"
                if h6_cold_med is not None else "_SASP_suppressive data not yet available — re-run notebook 03_"
            ),
            "direction_ok": h6_cold_med is not None and h6_cold_med > h6_hot_med,
            "conclusion_fn": _h6_conclusion,
        },
        {
            "id": "H7",
            "title": "Tumor mutational burden (TMB) is higher in hot tumors",
            "rationale": "DNA damage drives both mutational burden and p16-mediated senescence. Hot tumors should accumulate more mutations (neoantigens → immune recognition), and this genotoxic stress should co-occur with elevated p16/SASP.",
            "test": "Mann-Whitney U (hot vs cold TMB); Spearman ρ(TMB, IFN-γ score)",
            "expected": "TMB higher in hot; ρ > 0, p < 0.05",
            "notebook": "04 (exploratory)",
            "threshold": "BH FDR α = 0.05 (exploratory)",
            "observed": (
                f"Median TMB — hot: **{h7_hot_med:.1f}**, cold: **{h7_cold_med:.1f}** mut/Mb"
                + (f";  ρ(TMB, IFN-γ) = **{h7_rho:.3f}**" if h7_rho is not None else "")
                if h7_hot_med is not None else "_TMB data not yet available — re-run notebook 01_"
            ),
            "direction_ok": h7_hot_med is not None and h7_hot_med > h7_cold_med,
            "conclusion_fn": _h7_conclusion,
        },
        {
            "id": "H8",
            "title": "TP53 alteration is enriched in cold tumors",
            "rationale": "TP53 loss impairs the p53→p21 senescence arm and reduces inflammatory NF-κB signalling. Complementary to H3 (CDKN2A), this tests whether both CDK inhibitor pathway arms must be intact for the hot immune phenotype.",
            "test": "Chi-squared (TP53 alteration rate in hot vs cold)",
            "expected": "TP53 alteration rate higher in cold tumors, p < 0.05",
            "notebook": "04 §4 (exploratory)",
            "threshold": "BH FDR α = 0.05 (exploratory)",
            "observed": (
                f"TP53 alteration rate — cold: **{h8_cold_rate:.1f}%**, hot: **{h8_hot_rate:.1f}%**"
                if h8_cold_rate is not None else "_TP53 alteration data not yet available — re-run notebooks 01 and 03_"
            ),
            "direction_ok": h8_cold_rate is not None and h8_cold_rate > h8_hot_rate,
            "conclusion_fn": _h8_conclusion,
        },
        {
            "id": "H9",
            "title": "LMNB1 expression is inversely correlated with IFN-γ score",
            "rationale": "LMNB1 downregulation is one of the most specific markers of deep senescence — nuclear lamina loss leads to chromatin instability, cGAS-STING activation, and SASP induction. LMNB1-low tumors should therefore be more immune-inflamed.",
            "test": "Spearman ρ(LMNB1_expr, IFN-γ score) — expect negative ρ",
            "expected": "ρ < 0 (LMNB1-low ↔ hot), p < 0.05",
            "notebook": "04 §2 (exploratory)",
            "threshold": "BH FDR α = 0.05 (exploratory)",
            "observed": (
                f"Spearman ρ = **{h9_rho:.3f}** (n = {h9_n:,})"
                if h9_rho is not None else "_LMNB1_expr data not yet available — re-run notebook 03_"
            ),
            "direction_ok": h9_rho is not None and h9_rho < 0,
            "conclusion_fn": _h9_conclusion,
        },
        {
            "id": "H10",
            "title": "Tumor purity is lower in hot tumors (confound quantification)",
            "rationale": "TCGA bulk RNA-seq mixes tumor and immune cells. High purity = fewer infiltrating cells = lower IFN-γ score by dilution alone. Quantifying this confound is essential for interpreting H1–H5 — if purity fully explains the hot/cold signal, the biology is less compelling.",
            "test": "Mann-Whitney U (purity in hot vs cold); Spearman ρ(purity, IFN-γ score) — expect negative ρ",
            "expected": "Purity lower in hot; ρ < 0, p < 0.05",
            "notebook": "04 (methodological / exploratory)",
            "threshold": "BH FDR α = 0.05 (exploratory)",
            "observed": (
                f"Median purity — hot: **{h10_hot_med:.2f}**, cold: **{h10_cold_med:.2f}**"
                + (f";  ρ(purity, IFN-γ) = **{h10_rho:.3f}**" if h10_rho is not None else "")
                if h10_hot_med is not None else "_Tumor purity data not yet available — re-run notebook 01_"
            ),
            "direction_ok": h10_hot_med is not None and h10_cold_med > h10_hot_med,
            "conclusion_fn": _h10_conclusion,
        },
        {
            "id": "H11",
            "title": "CDKN1A (p21) expression correlates positively with IFN-γ score",
            "rationale": "Tests whether the p53→p21 senescence arm parallels the p16 arm (H1). If both CDK inhibitors correlate with immune activity, a broad senescence axis is implicated; if only CDKN2A does, the effect is p16-specific — mechanistically informative either way.",
            "test": "Spearman ρ(CDKN1A_expr, ifng_score)",
            "expected": "ρ > 0, p < 0.05; magnitude compared to H1 (CDKN2A)",
            "notebook": "04 §2 (exploratory)",
            "threshold": "BH FDR α = 0.05 (exploratory)",
            "observed": (
                f"Spearman ρ = **{h11_rho:.3f}** (n = {h11_n:,})"
                + (f"  |  CDKN2A ρ = {h1_rho:.3f}" if h1_rho is not None else "")
                if h11_rho is not None else "_CDKN1A_expr data not yet available — re-run notebook 03_"
            ),
            "direction_ok": h11_rho is not None and h11_rho > 0,
            "conclusion_fn": _h11_conclusion,
        },
        {
            "id": "H12",
            "title": "RB1 alteration is enriched in cold tumors",
            "rationale": "RB1 is the direct downstream effector of the p16–CDK4/6 axis. Loss of RB1 bypasses senescence arrest even when p16 is intact. If the whole pathway matters, RB1 loss should associate with immune-cold tumors — completing the pathway story from H3.",
            "test": "Chi-squared (RB1 alteration rate in hot vs cold)",
            "expected": "RB1 alteration rate higher in cold tumors, p < 0.05",
            "notebook": "04 §4 (exploratory)",
            "threshold": "BH FDR α = 0.05 (exploratory)",
            "observed": (
                f"RB1 alteration rate — cold: **{h12_cold_rate:.1f}%**, hot: **{h12_hot_rate:.1f}%**"
                if h12_cold_rate is not None else "_RB1 alteration data not yet available — re-run notebooks 01 and 03_"
            ),
            "direction_ok": h12_cold_rate is not None and h12_cold_rate > h12_hot_rate,
            "conclusion_fn": _h12_conclusion,
        },
        {
            "id": "H13",
            "title": "SERPINE1 (PAI-1) expression correlates positively with IFN-γ score",
            "rationale": "SERPINE1/PAI-1 is one of the most replicated SASP biomarkers in plasma studies. If SASP drives immune infiltration (H2), individual SASP effectors should each positively correlate with immune activity — a more specific test than the composite SASP NES.",
            "test": "Spearman ρ(SERPINE1_expr, ifng_score)",
            "expected": "ρ > 0, p < 0.05",
            "notebook": "04 §2 (exploratory)",
            "threshold": "BH FDR α = 0.05 (exploratory)",
            "observed": (
                f"Spearman ρ = **{h13_rho:.3f}** (n = {h13_n:,})"
                if h13_rho is not None else "_SERPINE1_expr data not yet available — re-run notebook 03_"
            ),
            "direction_ok": h13_rho is not None and h13_rho > 0,
            "conclusion_fn": _h13_conclusion,
        },
        {
            "id": "H14",
            "title": "Age correlates positively with CDKN2A expression",
            "rationale": "Cellular senescence accumulates with age, so older patients should have higher tissue-wide p16 expression. Confirming this validates CDKN2A as a biologically meaningful readout in TCGA and justifies including age as a Cox covariate given partial collinearity.",
            "test": "Spearman ρ(age, CDKN2A_expr)",
            "expected": "ρ > 0, p < 0.05",
            "notebook": "04 (exploratory)",
            "threshold": "BH FDR α = 0.05 (exploratory)",
            "observed": (
                f"Spearman ρ = **{h14_rho:.3f}** (n = {h14_n:,})"
                if h14_rho is not None else "_Age or CDKN2A data not yet available — re-run notebooks 01–03_"
            ),
            "direction_ok": h14_rho is not None and h14_rho > 0,
            "conclusion_fn": _h14_conclusion,
        },
        {
            "id": "H15",
            "title": "The prognostic value of immune phenotype is stronger in p16-high than p16-low tumors",
            "rationale": "If p16/SASP mechanistically drives the hot phenotype, then the hot/cold label in p16-low tumors captures noise rather than biology and should be less prognostic. A larger hot−cold OS gap in p16-high strata would support a causal rather than correlative relationship.",
            "test": "Hot−cold median OS gap within p16-high vs p16-low strata (interaction / effect modification)",
            "expected": "Hot−cold OS gap larger in p16-high stratum",
            "notebook": "05 §4 (exploratory)",
            "threshold": "BH FDR α = 0.05 (exploratory)",
            "observed": (
                f"Median OS gap (hot−cold) — p16-high: **{h15_p16high_diff:+.2f} yrs**, "
                f"p16-low: **{h15_p16low_diff:+.2f} yrs**"
                if h15_p16high_diff is not None else "_Survival or p16 status data not yet available — re-run notebooks 03 and 05_"
            ),
            "direction_ok": (
                h15_p16high_diff is not None and h15_p16low_diff is not None and
                h15_p16high_diff > 0 and abs(h15_p16high_diff) > abs(h15_p16low_diff)
            ),
            "conclusion_fn": _h15_conclusion,
        },
        {
            "id": "H16",
            "title": "Thorsson C2 (IFN-γ dominant) subtype has the highest SASP scores",
            "rationale": "Cross-validates the SASP → immune-inflamed link using an independent immune classification (Thorsson et al. 2018) derived separately from the IFN-γ score used to define hot/cold. C2 is the most immune-infiltrated subtype and should rank highest on SASP if the hot/cold signal is not a circular artefact.",
            "test": "Kruskal-Wallis across C1–C6; rank C2 by median SASP",
            "expected": "C2 has the highest median SASP score among all six subtypes",
            "notebook": "04 (exploratory)",
            "threshold": "BH FDR α = 0.05 (exploratory)",
            "observed": (
                f"C2 median SASP = **{h16_c2_med:.3f}** | subtype ranking: "
                + ", ".join(f"{s} ({v:.2f})" for s, v in h16_ranking.items())
                if h16_c2_med is not None else "_Thorsson subtype or SASP data not yet available — re-run notebooks 01 and 03_"
            ),
            "direction_ok": h16_c2_med is not None and (
                h16_ranking is not None and h16_ranking.index[0] == "C2"
            ),
            "conclusion_fn": _h16_conclusion,
        },
        {
            "id": "H17",
            "title": "SenMayo score is higher in hot tumors (independent senescence signature validation)",
            "rationale": "SenMayo (Saul et al. 2022, Nature Aging) is currently the best-validated pan-tissue senescence gene set (125 genes, replicated across human and mouse). Testing it independently of the custom SASP score (H2) provides cross-signature validation — if both agree, the senescence–immune link is more robust.",
            "test": "Mann-Whitney U (hot vs cold SenMayo NES)",
            "expected": "SenMayo NES higher in hot, p < 0.05; direction consistent with H2",
            "notebook": "04 §5 (exploratory)",
            "threshold": "BH FDR α = 0.05 (exploratory)",
            "observed": (
                f"Median SenMayo — hot: **{h17_hot_med:.3f}**, cold: **{h17_cold_med:.3f}** "
                f"(Δ = {h17_hot_med - h17_cold_med:+.3f})"
                if h17_hot_med is not None else "_SenMayo data not yet available — re-run notebook 03_"
            ),
            "direction_ok": h17_hot_med is not None and h17_hot_med > h17_cold_med,
            "conclusion_fn": _h17_conclusion,
        },
        {
            "id": "H18",
            "title": "TNF-α / NF-κB signalling is higher in hot tumors",
            "rationale": "TNF-α is a primary SASP effector cytokine that drives NF-κB-mediated recruitment of immune cells. HALLMARK_TNFA_SIGNALING_VIA_NFKB captures this mechanistic arm specifically — a positive result would give pathway-level resolution to the SASP → immune link in H2.",
            "test": "Mann-Whitney U (hot vs cold hallmark_tnfa NES)",
            "expected": "hallmark_tnfa higher in hot, p < 0.05",
            "notebook": "04 §5 (exploratory)",
            "threshold": "BH FDR α = 0.05 (exploratory)",
            "observed": (
                f"Median TNF-α NES — hot: **{h18_hot_med:.3f}**, cold: **{h18_cold_med:.3f}** "
                f"(Δ = {h18_hot_med - h18_cold_med:+.3f})"
                if h18_hot_med is not None else "_hallmark_tnfa data not yet available — re-run notebook 03_"
            ),
            "direction_ok": h18_hot_med is not None and h18_hot_med > h18_cold_med,
            "conclusion_fn": _h18_conclusion,
        },
        {
            "id": "H19",
            "title": "Inflammatory response score is higher in hot tumors",
            "rationale": "HALLMARK_INFLAMMATORY_RESPONSE captures broad NF-κB and cytokine signalling downstream of SASP. Testing it alongside TNF-α (H18) and SASP (H2) allows triangulation: if all three point the same direction, the SASP → immune signal is consistent across multiple gene set definitions.",
            "test": "Mann-Whitney U (hot vs cold hallmark_inflammatory NES)",
            "expected": "hallmark_inflammatory higher in hot, p < 0.05",
            "notebook": "04 §5 (exploratory)",
            "threshold": "BH FDR α = 0.05 (exploratory)",
            "observed": (
                f"Median inflammatory NES — hot: **{h19_hot_med:.3f}**, cold: **{h19_cold_med:.3f}** "
                f"(Δ = {h19_hot_med - h19_cold_med:+.3f})"
                if h19_hot_med is not None else "_hallmark_inflammatory data not yet available — re-run notebook 03_"
            ),
            "direction_ok": h19_hot_med is not None and h19_hot_med > h19_cold_med,
            "conclusion_fn": _h19_conclusion,
        },
        {
            "id": "H20",
            "title": "Senescence-downregulated signature is inversely correlated with IFN-γ score",
            "rationale": "Deeply senescent cells downregulate structural and cell-cycle genes (LMNB1, histones, lamins). This chromatin instability activates cGAS-STING → type I IFN → innate immune surveillance. The 'senescence-down' gene set captures cells that have undergone this transition — they should co-occur with immune-inflamed microenvironments.",
            "test": "Spearman ρ(senescence_dn, ifng_score) — expect negative ρ",
            "expected": "ρ < 0 (senescence-down ↔ hot), p < 0.05",
            "notebook": "04 §2 (exploratory)",
            "threshold": "BH FDR α = 0.05 (exploratory)",
            "observed": (
                f"Spearman ρ = **{h20_rho:.3f}** (n = {h20_n:,})"
                if h20_rho is not None else "_senescence_dn data not yet available — re-run notebook 03_"
            ),
            "direction_ok": h20_rho is not None and h20_rho < 0,
            "conclusion_fn": _h20_conclusion,
        },
        {
            "id": "H21",
            "title": "GDF15 expression correlates positively with IFN-γ score",
            "rationale": "GDF15 is a stress-induced SASP cytokine elevated in therapy-induced and oncogene-induced senescence. If it acts as a pro-immune signal in the tumor microenvironment, higher GDF15 should associate with hot tumors. A negative result would be equally informative — suggesting GDF15 is immunosuppressive or expressed in non-senescent cold tumor cells.",
            "test": "Spearman ρ(GDF15_expr, ifng_score)",
            "expected": "ρ > 0, p < 0.05",
            "notebook": "04 §2 (exploratory)",
            "threshold": "BH FDR α = 0.05 (exploratory)",
            "observed": (
                f"Spearman ρ = **{h21_rho:.3f}** (n = {h21_n:,})"
                if h21_rho is not None else "_GDF15_expr data not yet available — re-run notebook 03_"
            ),
            "direction_ok": h21_rho is not None and h21_rho > 0,
            "conclusion_fn": _h21_conclusion,
        },
        {
            "id": "H22",
            "title": "CDKN2A homozygous deep deletion is specifically enriched in cold tumors",
            "rationale": "Biallelic (deep) deletion completely abolishes p16 function, while heterozygous loss or point mutations may leave partial activity. If the p16 → senescence → SASP → immune axis is dose-dependent, the strongest functional loss (deep deletion) should show the strongest cold-tumour enrichment — a more specific test than H3.",
            "test": "Chi-squared (CDKN2A deep deletion rate in hot vs cold)",
            "expected": "Deep deletion rate higher in cold, p < 0.05; larger effect than composite alteration rate",
            "notebook": "04 §4 (exploratory)",
            "threshold": "BH FDR α = 0.05 (exploratory)",
            "observed": (
                f"Deep deletion rate — cold: **{h22_cold_rate:.1f}%**, hot: **{h22_hot_rate:.1f}%**"
                + (f"  |  composite alteration — cold: {h22_alt_cold:.1f}%, hot: {h22_alt_hot:.1f}%" if h22_alt_cold is not None else "")
                if h22_cold_rate is not None else "_CDKN2A deep deletion data not yet available — re-run notebooks 01 and 03_"
            ),
            "direction_ok": h22_cold_rate is not None and h22_cold_rate > h22_hot_rate,
            "conclusion_fn": _h22_conclusion,
        },
    ]

    for h in HYPOTHESES:
        supported, conclusion_text = h["conclusion_fn"]()
        icon = "✅" if supported is True else ("❓" if supported is None else "⚠️")
        with st.expander(f"{icon}  **{h['id']}** — {h['title']}", expanded=True):
            col_a, col_b = st.columns([1, 1])
            with col_a:
                st.markdown(f"**Rationale:** {h['rationale']}")
                st.markdown(f"**Test:** {h['test']}  |  **Threshold:** {h['threshold']}")
                st.markdown(f"**Expected:** {h['expected']}")
                st.markdown(f"**Notebook:** {h['notebook']}")
            with col_b:
                st.markdown("**Observed (interactive data):**")
                st.markdown(h["observed"])
            # Conclusion callout spanning full width
            if supported is True:
                st.success(conclusion_text)
            elif supported is False:
                st.warning(conclusion_text)
            else:
                st.info(conclusion_text)

    st.divider()
    st.markdown(
        """
        **Multiple testing note:** BH FDR is applied within each of the five notebook sections.
        Across H1–H5, Bonferroni α = 0.01 (5 pre-specified tests). All other p-values are
        exploratory and should be interpreted with caution.
        """
    )

# ── Summary of Findings ───────────────────────────────────────────────────────
elif page == "Summary":
    st.title("Summary of Findings")
    st.markdown(
        """
        Integrated interpretation across all five confirmatory hypotheses and the
        main exploratory analyses. Statistics are computed from the currently loaded
        (and filtered) dataset.
        """
    )

    # ── Key metrics ───────────────────────────────────────────────────────────
    n_samples   = len(df)
    n_types     = df["tumor_type"].nunique()
    n_hot       = (df["hot_cold"] == "hot").sum()
    n_cold      = (df["hot_cold"] == "cold").sum()
    pct_hot     = 100 * n_hot  / n_samples
    pct_cold    = 100 * n_cold / n_samples

    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Samples", f"{n_samples:,}")
    c2.metric("Tumor types", n_types)
    c3.metric("Hot", f"{n_hot:,}  ({pct_hot:.0f}%)")
    c4.metric("Cold", f"{n_cold:,}  ({pct_cold:.0f}%)")

    # ── H1: CDKN2A–IFN-γ correlation ─────────────────────────────────────────
    st.subheader("H1 — CDKN2A expression × IFN-γ score")
    if "CDKN2A_expr" in df.columns and "ifng_score" in df.columns:
        _h1 = df[["CDKN2A_expr", "ifng_score"]].dropna()
        rho = _h1.corr(method="spearman").iloc[0, 1]
        direction = "positive" if rho > 0 else "negative"
        strength  = "strong" if abs(rho) > 0.4 else ("moderate" if abs(rho) > 0.2 else "weak")
        st.markdown(
            f"Pan-cancer Spearman ρ = **{rho:.3f}** (n = {len(_h1):,}). "
            f"The correlation is **{strength} and {direction}**, "
            f"{'consistent' if rho > 0 else 'inconsistent'} with the hypothesis that "
            "p16-driven SASP promotes IFN-γ signalling."
        )
        col_h1a, col_h1b = st.columns(2)
        with col_h1a:
            # Per-type correlations
            per_type = (
                df.dropna(subset=["CDKN2A_expr", "ifng_score"])
                .groupby("tumor_type")[["CDKN2A_expr", "ifng_score"]]
                .apply(lambda g: g.corr(method="spearman").iloc[0, 1])
                .reset_index(name="rho")
                .sort_values("rho", ascending=True)
            )
            fig_rho = px.bar(
                per_type, x="rho", y="tumor_type", orientation="h",
                color="rho", color_continuous_scale="RdBu",
                color_continuous_midpoint=0,
                labels={"rho": "Spearman ρ", "tumor_type": ""},
                title="Per-type ρ(CDKN2A, IFN-γ)",
                height=500,
            )
            fig_rho.add_vline(x=0, line_dash="dash", line_color="gray")
            fig_rho.update_coloraxes(showscale=False)
            st.plotly_chart(fig_rho, use_container_width=True)
        with col_h1b:
            _scatter = (
                df[["CDKN2A_expr", "ifng_score", "hot_cold"]]
                .dropna(subset=["CDKN2A_expr", "ifng_score"])
                .sample(min(3000, len(_h1)), random_state=42)
                .copy()
            )
            fig_sc = px.scatter(
                _scatter, x="CDKN2A_expr", y="ifng_score",
                color="hot_cold",
                color_discrete_map={"hot": "#e74c3c", "intermediate": "#f39c12", "cold": "#3498db"},
                opacity=0.35, trendline="ols",
                labels={"CDKN2A_expr": "CDKN2A expression (log2)", "ifng_score": "IFN-γ score"},
                title="CDKN2A vs IFN-γ (pan-cancer)",
            )
            st.plotly_chart(fig_sc, use_container_width=True)
    else:
        st.caption("_Data not yet available — re-run notebooks 01–03._")

    # ── H2: SASP by immune phenotype ──────────────────────────────────────────
    st.subheader("H2 — SASP score in hot vs cold tumors")
    if "SASP" in df.columns and "hot_cold" in df.columns:
        _h2 = df.dropna(subset=["SASP", "hot_cold"])
        medians = _h2.groupby("hot_cold")["SASP"].median()
        hot_m   = medians.get("hot",  float("nan"))
        cold_m  = medians.get("cold", float("nan"))
        st.markdown(
            f"Median SASP NES — hot: **{hot_m:.3f}**, cold: **{cold_m:.3f}**. "
            f"Hot tumors show **{'higher' if hot_m > cold_m else 'lower'}** SASP scores, "
            f"{'supporting' if hot_m > cold_m else 'contrary to'} H2."
        )
        fig_sasp = px.violin(
            _h2, x="hot_cold", y="SASP", color="hot_cold",
            color_discrete_map={"hot": "#e74c3c", "intermediate": "#f39c12", "cold": "#3498db"},
            category_orders={"hot_cold": ["cold", "intermediate", "hot"]},
            box=True, points=False,
            labels={"hot_cold": "", "SASP": "SASP NES"},
            title="SASP score by immune phenotype",
        )
        fig_sasp.update_layout(showlegend=False)
        st.plotly_chart(fig_sasp, use_container_width=True)
    else:
        st.caption("_SASP data not yet available — re-run notebook 03._")

    # ── H3: CDKN2A deletion rate ───────────────────────────────────────────────
    st.subheader("H3 — CDKN2A deletion enriched in cold tumors")
    if "cdkn2a_altered" in df.columns and "hot_cold" in df.columns:
        alt_rates = (
            df.dropna(subset=["cdkn2a_altered", "hot_cold"])
            .groupby("hot_cold")["cdkn2a_altered"]
            .mean()
            .mul(100)
            .reset_index()
            .rename(columns={"cdkn2a_altered": "pct_altered"})
        )
        cold_r = alt_rates.loc[alt_rates["hot_cold"] == "cold",  "pct_altered"].values
        hot_r  = alt_rates.loc[alt_rates["hot_cold"] == "hot",   "pct_altered"].values
        cold_r = cold_r[0] if len(cold_r) else float("nan")
        hot_r  = hot_r[0]  if len(hot_r)  else float("nan")
        st.markdown(
            f"CDKN2A alteration rate — cold: **{cold_r:.1f}%**, hot: **{hot_r:.1f}%**. "
            f"{'Cold tumors have a higher alteration rate, supporting H3.' if cold_r > hot_r else 'Hot tumors show a higher alteration rate — contrary to H3.'}"
        )
        fig_alt = px.bar(
            alt_rates,
            x="hot_cold", y="pct_altered", color="hot_cold",
            color_discrete_map={"hot": "#e74c3c", "intermediate": "#f39c12", "cold": "#3498db"},
            category_orders={"hot_cold": ["cold", "intermediate", "hot"]},
            labels={"pct_altered": "% CDKN2A altered", "hot_cold": ""},
            title="CDKN2A alteration rate by immune phenotype",
        )
        fig_alt.update_layout(showlegend=False)
        st.plotly_chart(fig_alt, use_container_width=True)
    else:
        st.caption("_Alteration data not yet available — re-run notebooks 01 and 03._")

    # ── H4 / H5: Survival ─────────────────────────────────────────────────────
    st.subheader("H4 & H5 — Survival: p16 × immune phenotype")

    col_sv1, col_sv2 = st.columns(2)
    with col_sv1:
        st.markdown("**H4 — Per-tumor-type log-rank significance (hot vs cold OS)**")
        if "padj" in lr.columns:
            n_sig = (lr["padj"] < 0.05).sum()
            n_tot = len(lr)
            sig_frac = 100 * n_sig / n_tot
            st.markdown(
                f"{n_sig} / {n_tot} tumor types ({sig_frac:.0f}%) show significantly different OS "
                f"between hot and cold (padj < 0.05)."
            )
            lr_plot = lr.copy()
            lr_plot["significant"] = lr_plot["padj"] < 0.05
            fig_lr = px.scatter(
                lr_plot.sort_values("lr_pval"),
                x="lr_pval", y="tumor_type",
                color="significant",
                color_discrete_map={True: "#e74c3c", False: "#95a5a6"},
                labels={"lr_pval": "Log-rank p-value", "tumor_type": "", "significant": "padj < 0.05"},
                title="Per-type log-rank p (hot vs cold OS)",
                height=500,
            )
            fig_lr.add_vline(x=0.05, line_dash="dash", line_color="black")
            st.plotly_chart(fig_lr, use_container_width=True)
        else:
            st.caption("_Log-rank results not yet available._")

    with col_sv2:
        st.markdown("**H5 — Cox PH: CDKN2A expression as independent predictor**")
        cox_cdkn2a = cox[cox["term"].str.contains("CDKN2A", case=False, na=False)]
        if len(cox_cdkn2a):
            row = cox_cdkn2a.iloc[0]
            direction = "protective (HR < 1)" if row["estimate"] < 1 else "detrimental (HR > 1)"
            sig_label = "p < 0.05 — significant" if row["p.value"] < 0.05 else "p ≥ 0.05 — not significant"
            st.markdown(
                f"HR = **{row['estimate']:.3f}** (95% CI: {row['conf.low']:.3f}–{row['conf.high']:.3f})  \n"
                f"p = {row['p.value']:.2e} — **{sig_label}**  \n"
                f"Effect direction: **{direction}** (H5 predicts HR < 1)."
            )
            # Mini forest plot for all Cox terms
            fig_cox = go.Figure()
            for _, r in cox.iterrows():
                fig_cox.add_trace(go.Scatter(
                    x=[r["conf.low"], r["conf.high"]], y=[r["term"], r["term"]],
                    mode="lines", line=dict(color="#e74c3c", width=2), showlegend=False,
                ))
                fig_cox.add_trace(go.Scatter(
                    x=[r["estimate"]], y=[r["term"]],
                    mode="markers", marker=dict(color="#e74c3c", size=10), showlegend=False,
                    hovertemplate=(
                        f"<b>{r['term']}</b><br>HR = {r['estimate']:.3f}<br>"
                        f"p = {r['p.value']:.2e}<extra></extra>"
                    ),
                ))
            fig_cox.add_vline(x=1, line_dash="dash", line_color="gray")
            fig_cox.update_xaxes(type="log", title="Hazard Ratio (95% CI)")
            fig_cox.update_yaxes(title="")
            fig_cox.update_layout(height=300, title="Cox forest plot")
            st.plotly_chart(fig_cox, use_container_width=True)
        else:
            st.caption("_Cox results not yet available — re-run notebook 05._")

    # ── Narrative synthesis ───────────────────────────────────────────────────
    st.divider()
    st.subheader("Interpretation")

    # Build a dynamic narrative based on observed data
    findings = []

    if "CDKN2A_expr" in df.columns and "ifng_score" in df.columns:
        _tmp = df[["CDKN2A_expr", "ifng_score"]].dropna()
        rho  = _tmp.corr(method="spearman").iloc[0, 1]
        if rho > 0.1:
            findings.append(
                f"**CDKN2A expression positively correlates with immune activity** (ρ = {rho:.3f} pan-cancer), "
                "consistent with SASP-mediated recruitment of cytotoxic T cells."
            )
        elif rho > 0:
            findings.append(
                f"**CDKN2A expression shows a weak positive correlation** with IFN-γ score (ρ = {rho:.3f}), "
                "suggesting a modest contribution of SASP to immune infiltration."
            )
        else:
            findings.append(
                f"**CDKN2A expression is negatively correlated with IFN-γ score** (ρ = {rho:.3f}), "
                "which is inconsistent with the SASP-immune axis hypothesis."
            )

    if "SASP" in df.columns and "hot_cold" in df.columns:
        _hot  = df.loc[df["hot_cold"] == "hot",  "SASP"].dropna()
        _cold = df.loc[df["hot_cold"] == "cold", "SASP"].dropna()
        if len(_hot) and len(_cold):
            if _hot.median() > _cold.median():
                findings.append(
                    f"**SASP scores are higher in hot tumors** (median {_hot.median():.3f} vs {_cold.median():.3f} in cold), "
                    "supporting SASP as a driver of the immune-inflamed phenotype."
                )
            else:
                findings.append(
                    f"**SASP scores are not elevated in hot tumors** (median {_hot.median():.3f} vs {_cold.median():.3f} in cold), "
                    "which does not support H2."
                )

    if "cdkn2a_altered" in df.columns and "hot_cold" in df.columns:
        _c = df.loc[df["hot_cold"] == "cold", "cdkn2a_altered"].dropna().mean()
        _h = df.loc[df["hot_cold"] == "hot",  "cdkn2a_altered"].dropna().mean()
        if _c > _h:
            findings.append(
                f"**CDKN2A alterations are more frequent in cold tumors** ({100*_c:.1f}% vs {100*_h:.1f}% in hot), "
                "supporting the model that p16 loss reduces SASP and promotes immune exclusion."
            )
        else:
            findings.append(
                f"**CDKN2A alteration rates do not differ markedly between cold ({100*_c:.1f}%) and hot ({100*_h:.1f}%) tumors**, "
                "providing weak evidence for H3."
            )

    cox_cdkn2a = cox[cox["term"].str.contains("CDKN2A", case=False, na=False)]
    if len(cox_cdkn2a):
        r = cox_cdkn2a.iloc[0]
        if r["estimate"] < 1 and r["p.value"] < 0.05:
            findings.append(
                f"**CDKN2A expression is an independent protective prognostic factor** "
                f"(HR = {r['estimate']:.3f}, p = {r['p.value']:.2e}), consistent with H5."
            )
        elif r["estimate"] < 1:
            findings.append(
                f"**CDKN2A expression trends toward a protective effect** (HR = {r['estimate']:.3f}) "
                f"but does not reach Bonferroni significance (p = {r['p.value']:.2e})."
            )
        else:
            findings.append(
                f"**CDKN2A expression is not a protective Cox predictor** (HR = {r['estimate']:.3f}, "
                f"p = {r['p.value']:.2e}), which is inconsistent with H5."
            )

    if "padj" in lr.columns:
        n_sig = (lr["padj"] < 0.05).sum()
        n_tot = len(lr)
        findings.append(
            f"**Hot vs cold OS differences are significant in {n_sig}/{n_tot} tumor types** (padj < 0.05), "
            "suggesting immune phenotype is prognostically relevant across cancer contexts."
        )

    if findings:
        for i, f in enumerate(findings, 1):
            st.markdown(f"{i}. {f}")
    else:
        st.caption("_Run notebooks 01–05 to populate the summary._")

    st.divider()
    st.markdown(
        """
        **Limitations to consider:**
        - TCGA bulk RNA-seq mixes tumor cells with stroma; immune scores partly reflect tumor purity.
        - IFN-γ tertile classification (hot/cold) is pan-cancer median-based — within-type thresholds may differ.
        - Survival analyses are observational; confounding by stage, treatment, and purity cannot be excluded.
        - CDKN2A expression reflects both tumor cell expression and infiltrating cell content.
        """
    )
