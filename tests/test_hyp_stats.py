import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import hyp_stats as hs

DATA = Path("/Users/david/git/cold-tumor-p16/data")


def test_cmh_single_stratum_matches_pearson_relation():
    t = [[10, 5], [4, 12]]
    r = hs.cmh_odds_ratio([t])
    assert r["or"] == pytest.approx(6.0)
    # CMH (Yates) = Yates-corrected Pearson chi2 * (n-1)/n
    n = 31
    chi2 = stats.chi2_contingency(np.array(t), correction=True)[0]
    assert r["p"] == pytest.approx(stats.chi2.sf(chi2 * (n - 1) / n, 1), rel=1e-6)
    assert r["ci_low"] < 6.0 < r["ci_high"]


def test_cmh_simpson_paradox_reverses_pooled_direction():
    # within each type cold has the higher deletion rate, but the type with
    # high deletion is mostly hot -> pooled direction is reversed.
    rows = []
    # (type, n_cold, deleted_cold, n_hot, deleted_hot)
    for typ, nc, dc, nh, dh in [("A", 200, 110, 1000, 500), ("B", 1000, 60, 200, 10)]:
        rows += [(typ, "cold", i < dc) for i in range(nc)]
        rows += [(typ, "hot", i < dh) for i in range(nh)]
    d = pd.DataFrame(rows, columns=["tumor_type", "hot_cold", "cdkn2a_deleted"])
    r = hs.stratified_deletion_test(d)
    assert r["pooled_cold_pct"] < r["pooled_hot_pct"]
    assert r["or"] > 1
    assert r["cold_gt_hot_strata"] == 2


def test_cmh_degenerate_returns_nan():
    assert np.isnan(hs.cmh_odds_ratio([])["or"])


def _surv(n=120, seed=1):
    rng = np.random.default_rng(seed)
    g = rng.integers(0, 2, n)
    t = rng.exponential(1 + g, n)
    e = (rng.random(n) < 0.8).astype(int)
    return t, e, g


def test_logrank_two_group_matches_scipy():
    t, e, g = _surv()
    chi2, df, p = hs.logrank_kgroup(t, e, g)
    assert df == 1
    if hasattr(stats, "logrank"):
        def cd(k):
            return stats.CensoredData(uncensored=t[(g == k) & (e == 1)], right=t[(g == k) & (e == 0)])
        res = stats.logrank(cd(0), cd(1))
        assert p == pytest.approx(res.pvalue, rel=1e-6)


def test_logrank_constant_strata_equals_unstratified():
    t, e, g = _surv()
    a = hs.logrank_kgroup(t, e, g)
    b = hs.logrank_kgroup(t, e, g, strata=np.zeros(len(t)))
    assert a[2] == pytest.approx(b[2])


def test_stratified_logrank_removes_baseline_confounding():
    # no group effect within strata, but the strata differ hugely in survival
    # and group composition -> unstratified significant, stratified not.
    rng = np.random.default_rng(3)
    rows = []
    for s, scale, p1 in [("fast", 0.5, 0.8), ("slow", 5.0, 0.2)]:
        for _ in range(300):
            rows.append((s, int(rng.random() < p1), rng.exponential(scale), 1))
    d = pd.DataFrame(rows, columns=["s", "g", "t", "e"])
    assert hs.logrank_kgroup(d.t, d.e, d.g)[2] < 1e-6
    assert hs.logrank_kgroup(d.t, d.e, d.g, d.s)[2] > 0.01


def test_km_median_and_not_reached():
    assert hs.km_median([1, 2, 3, 4], [1, 1, 1, 1]) == 2.0
    assert np.isnan(hs.km_median([1, 2, 3, 4], [1, 0, 0, 0]))


def test_ordering_check():
    m = {"hot / p16-high": 9.0, "cold / p16-low": 3.0, "hot / p16-low": 5.0}
    assert hs.h4_ordering_observed(m)
    assert not hs.h4_ordering_observed({"hot / p16-high": 7.8, "cold / p16-low": 9.5, "x": 7.0})


def test_add_group6_splits_within_type_after_dropping_missing_os():
    d = pd.DataFrame({
        "tumor_type": ["A"] * 5 + ["B"],
        "hot_cold": ["hot"] * 6,
        "CDKN2A_expr": [1, 2, 3, 4, 5, 9.0],
        "OS": [1, 1, 1, 1, np.nan, 1],
        "OS.time_years": [1.0] * 6,
    })
    g = hs.add_group6(d)
    assert "B" not in set(g.tumor_type)  # single-sample type has one level -> dropped
    assert (g["p16_status"].value_counts() == 2).all()
    assert set(g["group6"]) == {"hot / p16-low", "hot / p16-high"}


@pytest.mark.skipif(not (DATA / "sample_immune_labels.parquet").exists(),
                    reason="gitignored analysis data not available")
def test_documented_numbers_on_real_data():
    imm = pd.read_parquet(DATA / "sample_immune_labels.parquet")
    sen = pd.read_parquet(DATA / "sample_senescence.parquet")
    m = imm.join(sen, how="inner")
    r = hs.stratified_deletion_test(m)
    assert r["or"] == pytest.approx(1.297, abs=1e-3)
    assert (r["ci_low"], r["ci_high"]) == pytest.approx((1.124, 1.496), abs=1e-3)
    assert r["p"] == pytest.approx(3.6e-4, rel=0.05)
    assert (r["cold_gt_hot_strata"], r["n_strata"]) == (21, 29)
    clin = pd.read_csv(DATA / "tcga_clinical_survival.tsv", sep="\t", index_col=0)
    m = m.join(clin[["OS", "OS.time"]], how="left")
    m["OS.time_years"] = m["OS.time"] / 365.25
    h4 = hs.six_group_survival(hs.add_group6(m))
    assert h4["p_unstrat"] == pytest.approx(3.3e-3, rel=0.05)
    assert h4["p_strat"] == pytest.approx(0.077, abs=0.002)
    assert h4["medians"]["cold / p16-low"] == pytest.approx(9.50, abs=0.01)
    assert h4["medians"]["hot / p16-high"] == pytest.approx(7.83, abs=0.01)
    assert not hs.h4_ordering_observed(h4["medians"])
