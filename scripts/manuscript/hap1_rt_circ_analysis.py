#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from _mudata_utils import (
    CIRC_MODALITY_CANDIDATES,
    RNA_MODALITY_CANDIDATES,
    RT_MODALITY_CANDIDATES,
    common_obs_names,
    correlation_rows,
    fail,
    find_modality,
    get_or_compute_circ_metric,
    get_or_compute_rna_metric,
    get_or_compute_rt_metric,
    load_mudata,
    write_table,
)


def build_cell_table(path: Path) -> pd.DataFrame:
    mdata = load_mudata(path)
    rna = mdata.mod[find_modality(mdata, RNA_MODALITY_CANDIDATES)]
    circ = mdata.mod[find_modality(mdata, CIRC_MODALITY_CANDIDATES)]
    rt = mdata.mod[find_modality(mdata, RT_MODALITY_CANDIDATES)]
    shared = common_obs_names(rna, circ, rt)
    if not shared:
        fail("No shared cells across RNA, circRNA, and RT modalities.")

    frame = pd.DataFrame(index=pd.Index(shared, name="cell_id"))
    frame["circRNA_count"] = get_or_compute_circ_metric(mdata, circ, "circRNA_count").reindex(shared)
    frame["circRNA_total_support"] = get_or_compute_circ_metric(mdata, circ, "circRNA_total_support").reindex(shared)
    frame["detected_genes"] = get_or_compute_rna_metric(mdata, rna, "detected_genes").reindex(shared)
    frame["frac_rt_pos"] = get_or_compute_rt_metric(mdata, rt, "frac_rt_pos").reindex(shared)
    frame["frac_rt_neg"] = get_or_compute_rt_metric(mdata, rt, "frac_rt_neg").reindex(shared)
    return frame.reset_index()


def ols_table(frame: pd.DataFrame) -> pd.DataFrame:
    cols = ["circRNA_count", "detected_genes", "frac_rt_pos"]
    subset = frame[cols].dropna()
    n = int(len(subset))
    if n < 4:
        fail("OLS requires at least 4 complete cells for circRNA_count ~ detected_genes + frac_rt_pos.")
    y = subset["circRNA_count"].to_numpy(dtype=float)
    x = np.column_stack(
        [
            np.ones(n, dtype=float),
            subset["detected_genes"].to_numpy(dtype=float),
            subset["frac_rt_pos"].to_numpy(dtype=float),
        ]
    )
    terms = ["Intercept", "detected_genes", "frac_rt_pos"]
    if np.linalg.matrix_rank(x) < x.shape[1]:
        fail("OLS design matrix is rank deficient; check detected_genes and frac_rt_pos variation.")
    beta, _, _, _ = np.linalg.lstsq(x, y, rcond=None)
    fitted = x @ beta
    resid = y - fitted
    df_resid = n - x.shape[1]
    rss = float(np.sum(resid**2))
    tss = float(np.sum((y - y.mean()) ** 2))
    sigma2 = rss / df_resid
    cov = sigma2 * np.linalg.inv(x.T @ x)
    se = np.sqrt(np.diag(cov))
    t_values = beta / se
    p_values = 2 * stats.t.sf(np.abs(t_values), df=df_resid)
    r_squared = 1.0 - rss / tss if tss > 0 else np.nan
    return pd.DataFrame(
        {
            "term": terms,
            "estimate": beta,
            "std_error": se,
            "t": t_values,
            "p": p_values,
            "n": n,
            "df_resid": df_resid,
            "r_squared": r_squared,
        }
    )


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Analyze HAP1 scRR replication-timing state versus circRNA burden."
    )
    parser.add_argument("--input", required=True, type=Path, help="HAP1 RNA+circ+RT .h5mu file.")
    parser.add_argument("--correlations-out", required=True, type=Path, help="Output correlation TSV/CSV.")
    parser.add_argument("--ols-out", required=True, type=Path, help="Output OLS coefficient TSV/CSV.")
    parser.add_argument(
        "--scatter-out",
        type=Path,
        default=None,
        help="Optional per-cell metrics TSV/CSV for scatter plotting.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    frame = build_cell_table(args.input)
    pairs = [
        ("frac_rt_pos", "circRNA_count"),
        ("frac_rt_pos", "circRNA_total_support"),
        ("frac_rt_neg", "circRNA_count"),
        ("frac_rt_neg", "circRNA_total_support"),
        ("detected_genes", "circRNA_count"),
    ]
    write_table(pd.DataFrame(correlation_rows(frame, pairs)), args.correlations_out)
    write_table(ols_table(frame), args.ols_out)
    if args.scatter_out is not None:
        write_table(frame, args.scatter_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
