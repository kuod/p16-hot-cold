import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from utils import (classify_hot_cold, normalize_barcode, safe_qcut,
                   signature_score, spearman_table, validate_parquet, zscore_df)


def test_normalize_barcode():
    assert normalize_barcode("TCGA-AB-1234-01A-11R-A29S-07") == "TCGA-AB-1234-01A"
    assert normalize_barcode("TCGA-AB") == "TCGA-AB"


def test_zscore_constant_gene_is_zero():
    df = pd.DataFrame([[1.0, 2.0, 3.0], [5.0, 5.0, 5.0]], index=["a", "const"])
    z = zscore_df(df)
    assert (z.loc["const"] == 0).all()
    assert z.loc["a"].mean() == pytest.approx(0)


def test_signature_score_and_no_genes():
    expr = pd.DataFrame([[1.0, 2.0, 3.0], [3.0, 2.0, 1.0]], index=["A", "B"], columns=list("xyz"))
    assert signature_score(expr, ["A"]).tolist() == pytest.approx([-1, 0, 1])
    with pytest.raises(ValueError):
        signature_score(expr, ["NOPE"])


@pytest.mark.parametrize("n", [2, 3])
def test_classify_preserves_nan(n):
    s = pd.Series([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, np.nan])
    out = classify_hot_cold(s, n)
    assert pd.isna(out.iloc[-1])
    assert out.iloc[0] == "cold" and out.iloc[5] == "hot"


def test_classify_invalid_groups():
    with pytest.raises(ValueError):
        classify_hot_cold(pd.Series([1.0, 2.0]), 4)


def test_safe_qcut_degenerate():
    out = safe_qcut(pd.Series([1.0, 1.0, 1.0]), 2, labels=["lo", "hi"])
    assert out.isna().all()


def test_spearman_table_too_few_pairs():
    x = pd.DataFrame({"a": [1.0, 2.0, np.nan]})
    y = pd.DataFrame({"b": [2.0, np.nan, np.nan]})
    assert np.isnan(spearman_table(x, y).loc["a", "b"][0])


def test_validate_parquet_errors(tmp_path):
    with pytest.raises(FileNotFoundError):
        validate_parquet(tmp_path / "nope.parquet", ["a"])
    p = tmp_path / "d.parquet"
    pd.DataFrame({"a": range(5)}).to_parquet(p)
    with pytest.raises(ValueError):
        validate_parquet(p, ["b"], min_rows=1)
    with pytest.raises(ValueError):
        validate_parquet(p, ["a"], min_rows=100)
