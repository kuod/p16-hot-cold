"""
Shared utility functions for the hot/cold × p16 analysis.
"""

import numpy as np
import pandas as pd
from scipy import stats


def normalize_barcode(barcode: str) -> str:
    """
    Normalize a TCGA barcode to the first 4 hyphen-separated components
    (e.g. TCGA-AB-1234-01A-... → TCGA-AB-1234-01).
    Handles both full 7-component barcodes and already-shortened ones.
    """
    parts = str(barcode).split("-")
    return "-".join(parts[:4]) if len(parts) >= 4 else barcode


def zscore_df(df: pd.DataFrame) -> pd.DataFrame:
    """Z-score each gene (row) across samples."""
    return df.sub(df.mean(axis=1), axis=0).div(df.std(axis=1), axis=0)


def log2tpm(tpm: pd.DataFrame, pseudocount: float = 1.0) -> pd.DataFrame:
    """log2(TPM + pseudocount) transform."""
    return np.log2(tpm + pseudocount)


def signature_score(expr: pd.DataFrame, genes: list[str]) -> pd.Series:
    """
    Simple mean z-score signature scoring.

    Parameters
    ----------
    expr : genes × samples TPM or log-TPM DataFrame
    genes : gene list defining the signature

    Returns
    -------
    Per-sample signature score (Series indexed by sample barcode).
    """
    present = [g for g in genes if g in expr.index]
    missing = set(genes) - set(present)
    if missing:
        print(f"  Warning: {len(missing)} genes not found: {missing}")
    z = zscore_df(expr)
    return z.loc[present].mean(axis=0)


def classify_hot_cold(scores: pd.Series, n_groups: int = 2) -> pd.Series:
    """
    Median-split (n_groups=2) or tertile (n_groups=3) classification
    of a continuous immune/IFN-γ score into hot/cold labels.
    """
    if n_groups == 2:
        threshold = scores.median()
        return scores.map(lambda x: "hot" if x >= threshold else "cold")
    elif n_groups == 3:
        low, high = scores.quantile([1 / 3, 2 / 3])
        def _label(x):
            if x <= low:
                return "cold"
            elif x >= high:
                return "hot"
            return "intermediate"
        return scores.map(_label)
    raise ValueError("n_groups must be 2 or 3")


def validate_parquet(path, required_cols: list[str], min_rows: int = 100) -> pd.DataFrame:
    """
    Load a parquet file and validate it has the required columns and minimum row count.
    Raises ValueError with a descriptive message on failure.
    Returns the loaded DataFrame on success.
    """
    from pathlib import Path
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(
            f"Expected parquet not found: {path}\n"
            f"Re-run the upstream notebook to generate it."
        )
    df = pd.read_parquet(path)
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(
            f"{path.name}: missing expected columns: {missing}\n"
            f"Available columns: {list(df.columns)}"
        )
    if len(df) < min_rows:
        raise ValueError(
            f"{path.name}: only {len(df)} rows (expected ≥ {min_rows}). "
            f"Upstream notebook may have failed."
        )
    print(f"  ✓ {path.name}: {len(df):,} rows, {len(df.columns)} cols")
    return df


def safe_qcut(series: "pd.Series", q: int, labels=None) -> "pd.Series":
    """
    Like pd.qcut(series, q, labels=labels) but returns NaN (not raises) when
    the series has too few unique values to form q bins (e.g. all-deleted
    CDKN2A in a tiny cohort).
    """
    import pandas as pd
    try:
        return pd.qcut(series, q, labels=labels, duplicates="drop")
    except ValueError:
        return pd.Series([pd.NA] * len(series), index=series.index, dtype="object")


def spearman_table(x_df: pd.DataFrame, y_df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute Spearman correlation between every column in x_df and every
    column in y_df, returning a DataFrame of (rho, pvalue) tuples shaped
    (x_cols × y_cols).

    Both DataFrames must share the same index (sample barcodes).
    """
    results = []
    for xc in x_df.columns:
        row = {}
        for yc in y_df.columns:
            mask = x_df[xc].notna() & y_df[yc].notna()
            rho, pval = stats.spearmanr(x_df.loc[mask, xc], y_df.loc[mask, yc])
            row[yc] = (round(rho, 4), round(pval, 4))
        results.append(pd.Series(row, name=xc))
    return pd.DataFrame(results)
