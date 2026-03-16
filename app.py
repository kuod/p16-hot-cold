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
    master["group4"] = master["hot_cold"] + " / " + master["p16_status"]

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
    ["Overview", "Immune Classification", "Senescence Scores", "Correlations", "Survival", "Pathway", "Figures"],
)

# ── Overview ──────────────────────────────────────────────────────────────────
if page == "Overview":
    st.title("Hot/Cold Tumor × p16/CDKN2A — TCGA Pan-Cancer")
    st.markdown(
        """
        This app explores the relationship between **immune phenotype** (hot vs cold tumors)
        and **p16/CDKN2A senescence status** across TCGA pan-cancer RNA-seq data.

        | Hypothesis | Test |
        |---|---|
        | SASP-high tumors → hot (pro-inflammatory cytokines drive TIL infiltration) | Spearman correlation, Mann-Whitney |
        | CDKN2A-deleted tumors → cold | Chi-squared |
        | Senescence/immune phenotype predicts survival | Kaplan-Meier, Cox PH |
        """
    )

    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Samples", f"{len(df):,}")
    c2.metric("Tumor types", df["tumor_type"].nunique())
    c3.metric("Hot", f"{(df['hot_cold'] == 'hot').sum():,}")
    c4.metric("Cold", f"{(df['hot_cold'] == 'cold').sum():,}")

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

    st.subheader("Hot tumor fraction by cancer type")
    hot_frac = (
        df.groupby("tumor_type")["hot_cold"]
        .apply(lambda s: (s == "hot").mean())
        .reset_index()
        .rename(columns={"hot_cold": "hot_fraction"})
        .sort_values("hot_fraction", ascending=True)
    )
    fig = px.bar(
        hot_frac, x="hot_fraction", y="tumor_type", orientation="h",
        color="hot_fraction", color_continuous_scale=["#3498db", "#e74c3c"],
        labels={"hot_fraction": "Fraction hot", "tumor_type": ""},
        height=600,
    )
    fig.add_vline(x=0.5, line_dash="dash", line_color="gray")
    fig.update_coloraxes(showscale=False)
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
        med_hot  = age_hc.loc[age_hc["hot_cold"] == "hot",  "age"].median()
        med_cold = age_hc.loc[age_hc["hot_cold"] == "cold", "age"].median()
        fig_violin = px.violin(
            age_hc, x="hot_cold", y="age", color="hot_cold",
            color_discrete_map={"hot": "#e74c3c", "cold": "#3498db"},
            box=True, points=False,
            labels={"hot_cold": "", "age": "Age at diagnosis (years)"},
        )
        fig_violin.update_layout(showlegend=False)
        st.plotly_chart(fig_violin, use_container_width=True)
        st.caption(
            f"Median age — hot: {med_hot:.1f} years, cold: {med_cold:.1f} years"
        )

# ── Immune Classification ──────────────────────────────────────────────────────
elif page == "Immune Classification":
    st.title("Immune Classification")

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("IFN-γ score distribution")
        fig = px.histogram(
            df, x="ifng_score", color="hot_cold",
            color_discrete_map={"hot": "#e74c3c", "cold": "#3498db"},
            barmode="overlay", opacity=0.65, nbins=80,
            labels={"ifng_score": "IFN-γ signature score", "hot_cold": ""},
        )
        fig.add_vline(x=df["ifng_score"].median(), line_dash="dash", line_color="black",
                      annotation_text="median")
        st.plotly_chart(fig, use_container_width=True)

    with col2:
        st.subheader("CD8+ T cell fraction by phenotype")
        fig = px.box(
            df.dropna(subset=["cd8_fraction"]),
            x="hot_cold", y="cd8_fraction",
            color="hot_cold",
            color_discrete_map={"hot": "#e74c3c", "cold": "#3498db"},
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
        color_discrete_map={"hot": "#e74c3c", "cold": "#3498db"},
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

    with col1:
        score = st.selectbox("Score", score_cols)
        fig = px.violin(
            df.dropna(subset=[score, "hot_cold"]),
            x="hot_cold", y=score, color="hot_cold",
            color_discrete_map={"hot": "#e74c3c", "cold": "#3498db"},
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
                for pheno in ["hot", "cold"]:
                    sub = df[df["hot_cold"] == pheno]
                    pct = 100 * sub[col].mean() if len(sub) > 0 else 0.0
                    rows.append({"Gene": gene, "Phenotype": pheno, "% altered": pct})
            alt_df = pd.DataFrame(rows)
            fig = px.bar(
                alt_df, x="Gene", y="% altered", color="Phenotype",
                barmode="group",
                color_discrete_map={"hot": "#e74c3c", "cold": "#3498db"},
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
        color_discrete_map={"hot": "#e74c3c", "cold": "#3498db"},
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
            color_discrete_map={"hot": "#e74c3c", "cold": "#3498db"},
            labels={"immune_subtype": "Immune subtype", "pct": "% samples", "hot_cold": ""},
            title="Hot/cold fraction within Thorsson immune subtypes",
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
        color_discrete_map={"hot": "#e74c3c", "cold": "#3498db"},
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

    st.subheader("Group sizes (hot/cold × p16-high/low)")
    group_counts = (
        df.groupby(["tumor_type", "group4"])
        .size()
        .reset_index(name="n")
        .pivot(index="tumor_type", columns="group4", values="n")
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
    strat_opts = [c for c in ["hot_cold", "p16_status", "group4"] if c in df.columns]

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
        "SASP score: hot vs cold": "sasp_score_hot_cold.png",
        "KM: hot vs cold": "km_hot_cold.png",
        "KM: 4-group": "km_4group.png",
        "Cox forest plot": "cox_forest_plot.png",
        "Fig 1: Senescence heatmap": "fig1_senescence_heatmap.png",
        "Fig 2: UMAP": "fig2_umap.png",
        "Fig 3: Correlation matrix": "fig3_correlation_matrix.png",
        "Fig 4: Survival composite": "fig4_survival_composite.png",
        "p16 pathway diagram": "p16_pathway.png",
        "Alteration rates: hot vs cold (bar)": "alteration_hot_cold_grouped_bar.png",
        "Alteration heatmap by tumor type": "alteration_hot_cold_heatmap.png",
        "Age by tumor type": "age_by_tumor_type.png",
        "Age: hot vs cold": "age_hot_cold.png",
        "Age × senescence correlations": "age_senescence_correlations.png",
        "KM: age groups": "km_age_groups.png",
        "KM: DSS hot vs cold": "km_dss_hot_cold.png",
        "KM: PFI hot vs cold": "km_pfi_hot_cold.png",
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
