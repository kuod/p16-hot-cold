"""Statistics behind the dashboard hypothesis cards (numpy/scipy/pandas only).

These reproduce the tests used in notebooks 04/05 so card verdicts agree with them:
  H3  tumor-type-stratified Cochran-Mantel-Haenszel (R mantelhaen.test, exact=FALSE)
  H4  6-group log-rank, unstratified and stratified by tumor type (R survdiff)
"""
import numpy as np
import pandas as pd
from scipy import stats


def cmh_odds_ratio(tables, correct=True):
    """Mantel-Haenszel common OR, Robins-Breslow-Greenland 95% CI, CMH p-value.

    tables: iterable of 2x2 arrays [[a, b], [c, d]] (one per stratum).
    Strata with n < 2 are skipped. Matches R mantelhaen.test(exact=FALSE).
    """
    R = S = PR = PS_QR = QS = 0.0
    A = EA = V = 0.0
    k = 0
    for t in tables:
        a, b, c, d = (float(x) for x in np.asarray(t).ravel())
        n = a + b + c + d
        if n < 2:
            continue
        k += 1
        r, s = a * d / n, b * c / n
        p, q = (a + d) / n, (b + c) / n
        R += r; S += s
        PR += p * r; PS_QR += p * s + q * r; QS += q * s
        A += a
        r1, c1 = a + b, a + c
        EA += r1 * c1 / n
        V += r1 * (c + d) * c1 * (b + d) / (n * n * (n - 1))
    if k == 0 or R == 0 or S == 0 or V == 0:
        return {"or": np.nan, "ci_low": np.nan, "ci_high": np.nan, "p": np.nan, "n_strata": k}
    or_mh = R / S
    se = np.sqrt(PR / (2 * R ** 2) + PS_QR / (2 * R * S) + QS / (2 * S ** 2))
    lo, hi = np.exp(np.log(or_mh) + np.array([-1, 1]) * 1.959963984540054 * se)
    dev = abs(A - EA) - (0.5 if correct else 0.0)
    chi2 = max(dev, 0.0) ** 2 / V
    return {"or": or_mh, "ci_low": lo, "ci_high": hi,
            "p": float(stats.chi2.sf(chi2, 1)), "n_strata": k}


def stratified_deletion_test(df, flag_col="cdkn2a_deleted", group_col="hot_cold",
                             strata_col="tumor_type"):
    """H3: cold-vs-hot odds of deletion, stratified by tumor type (notebook 04 §4).

    OR > 1 means deletion is more common in cold. Also returns the number of
    strata where the cold deletion rate exceeds the hot one, and pooled rates.
    """
    d = df[[flag_col, group_col, strata_col]].dropna()
    d = d[d[group_col].isin(["cold", "hot"])].copy()
    d[flag_col] = d[flag_col].astype(bool)
    tables, cold_gt_hot = [], 0
    for _, g in d.groupby(strata_col):
        # notebook keeps only strata with both groups and both deletion states
        if g[group_col].nunique() < 2 or g[flag_col].nunique() < 2:
            continue
        cold, hot = g[g[group_col] == "cold"][flag_col], g[g[group_col] == "hot"][flag_col]
        tables.append([[cold.sum(), hot.sum()], [(~cold).sum(), (~hot).sum()]])
        cold_gt_hot += int(cold.mean() > hot.mean())
    res = cmh_odds_ratio(tables)
    res["cold_gt_hot_strata"] = cold_gt_hot
    pc, ph = d.loc[d[group_col] == "cold", flag_col], d.loc[d[group_col] == "hot", flag_col]
    res["pooled_cold_pct"] = 100 * pc.mean() if len(pc) else np.nan
    res["pooled_hot_pct"] = 100 * ph.mean() if len(ph) else np.nan
    return res


def logrank_kgroup(time, event, group, strata=None):
    """k-group log-rank chi-square test, optionally stratified (R survdiff).

    Returns (chi2, df, p). Stratified: O-E and variance summed over strata.
    """
    t = np.asarray(time, float); e = np.asarray(event, int)
    g = pd.Categorical(np.asarray(group)); codes = g.codes
    k_all = len(g.categories)
    s = np.zeros(len(t), int) if strata is None else pd.Categorical(np.asarray(strata)).codes
    OE = np.zeros(k_all); Vt = np.zeros((k_all, k_all))
    for sv in np.unique(s):
        m = s == sv
        ts, es, gs = t[m], e[m], codes[m]
        for tt in np.unique(ts[es == 1]):
            at_risk = ts >= tt
            n_k = np.bincount(gs[at_risk], minlength=k_all).astype(float)
            n = n_k.sum()
            dead = (ts == tt) & (es == 1)
            d_k = np.bincount(gs[dead], minlength=k_all).astype(float)
            d = d_k.sum()
            OE += d_k - d * n_k / n
            if n > 1:
                f = d * (n - d) / (n * n * (n - 1))
                Vt += f * (n * np.diag(n_k) - np.outer(n_k, n_k))
    # drop groups with no subjects / no variance (R: groups with n == 0)
    present = np.bincount(codes, minlength=k_all) > 0
    idx = np.where(present)[0]
    df_ = len(idx) - 1
    if df_ < 1:
        return np.nan, 0, np.nan
    oe, V = OE[idx][:-1], Vt[np.ix_(idx, idx)][:-1, :-1]
    chi2 = float(oe @ np.linalg.pinv(V) @ oe)
    return chi2, df_, float(stats.chi2.sf(chi2, df_))


def km_median(time, event):
    """Kaplan-Meier median survival time (NaN if never reached)."""
    t = np.asarray(time, float); e = np.asarray(event, int)
    surv = 1.0
    for tt in np.unique(t[e == 1]):
        n = (t >= tt).sum(); d = ((t == tt) & (e == 1)).sum()
        surv *= 1 - d / n
        if surv <= 0.5:
            return float(tt)
    return np.nan


def add_group6(df, time_col="OS.time_years", event_col="OS", expr_col="CDKN2A_expr",
               strata_col="tumor_type"):
    """Notebook 05 §2: restrict to samples with OS, split CDKN2A within tumor type
    (dplyr::ntile(2): rank ties 'first', extra sample goes to the low half), drop
    types with a single level, build 'hot_cold / p16-status' labels in `group6`."""
    d = df.dropna(subset=[time_col, event_col, expr_col, strata_col, "hot_cold"]).copy()
    r = d.groupby(strata_col)[expr_col].rank(method="first")
    n = d.groupby(strata_col)[expr_col].transform("size")
    d["p16_status"] = np.where(np.floor(2 * (r - 1) / n) + 1 == 2, "p16-high", "p16-low")
    d = d[d.groupby(strata_col)["p16_status"].transform("nunique") == 2]
    d["group6"] = d["hot_cold"].astype(str) + " / " + d["p16_status"]
    return d


def six_group_survival(df, group_col="group6", time_col="OS.time_years", event_col="OS",
                       strata_col="tumor_type", min_n=30, min_events=5):
    """H4 (notebook 05 §4): 6-group log-rank, unstratified and stratified + KM medians.

    `df` should already carry `group6` (see add_group6).
    """
    d = df.dropna(subset=[group_col, time_col, event_col]).copy()
    d = d[~d[group_col].astype(str).str.contains("nan")]
    d = d[d[time_col] > 0] if False else d
    g = d.groupby(group_col)
    keep = g[event_col].transform("size").ge(min_n) & g[event_col].transform("sum").ge(min_events)
    d = d[keep]
    out = {"n": len(d), "groups": sorted(d[group_col].unique())}
    _, dfree, out["p_unstrat"] = logrank_kgroup(d[time_col], d[event_col], d[group_col])
    _, _, out["p_strat"] = logrank_kgroup(d[time_col], d[event_col], d[group_col], d[strata_col])
    out["df"] = dfree
    out["medians"] = {k: km_median(x[time_col], x[event_col]) for k, x in d.groupby(group_col)}
    return out


def h4_ordering_observed(medians, best="hot / p16-high", worst="cold / p16-low"):
    """Claimed ordering: best group has the longest median OS, worst the shortest."""
    v = {k: m for k, m in medians.items() if not np.isnan(m)}
    if best not in v or worst not in v:
        return False
    return v[best] == max(v.values()) and v[worst] == min(v.values())


def within_type_spearman(df, x, y, strata_col="tumor_type", min_n=10):
    """Median of per-tumor-type Spearman rho (H1 summary)."""
    d = df[[x, y, strata_col]].dropna()
    rhos = [stats.spearmanr(g[x], g[y])[0] for _, g in d.groupby(strata_col) if len(g) >= min_n]
    return float(np.nanmedian(rhos)) if rhos else np.nan
