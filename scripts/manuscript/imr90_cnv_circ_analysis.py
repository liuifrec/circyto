#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from _mudata_utils import (
    CIRC_MODALITY_CANDIDATES,
    CNV_MODALITY_CANDIDATES,
    common_obs_names,
    correlation_rows,
    dense_matrix,
    fail,
    feature_coordinates,
    find_modality,
    get_or_compute_circ_metric,
    get_or_compute_cnv_metric,
    load_mudata,
    split_host_genes,
    write_table,
)


def build_cell_table(mdata, circ, cnv) -> pd.DataFrame:
    shared = common_obs_names(circ, cnv)
    if not shared:
        fail("No shared cells across circRNA and CNV modalities.")
    frame = pd.DataFrame(index=pd.Index(shared, name="cell_id"))
    frame["circRNA_count"] = get_or_compute_circ_metric(mdata, circ, "circRNA_count").reindex(shared)
    frame["circRNA_total_support"] = get_or_compute_circ_metric(mdata, circ, "circRNA_total_support").reindex(shared)
    frame["frac_non_diploid"] = get_or_compute_cnv_metric(mdata, cnv, "frac_non_diploid").reindex(shared)
    frame["frac_loss"] = get_or_compute_cnv_metric(mdata, cnv, "frac_loss").reindex(shared)
    frame["frac_gain"] = get_or_compute_cnv_metric(mdata, cnv, "frac_gain").reindex(shared)
    return frame.reset_index()


def host_gene_enrichment(
    circ,
    cell_table: pd.DataFrame,
    *,
    high_quantile: float,
    low_quantile: float,
    metric: str,
) -> pd.DataFrame:
    if "host_gene" not in circ.var.columns:
        fail("circ.var is missing `host_gene`; run host-gene annotation or repair first.")
    if metric not in cell_table.columns:
        fail(f"Missing CNV metric `{metric}` for CNV-high/CNV-low comparison.")
    high_cutoff = float(cell_table[metric].quantile(high_quantile))
    low_cutoff = float(cell_table[metric].quantile(low_quantile))
    high_cells = cell_table.loc[cell_table[metric] >= high_cutoff, "cell_id"].astype(str).tolist()
    low_cells = cell_table.loc[cell_table[metric] <= low_cutoff, "cell_id"].astype(str).tolist()
    if not high_cells or not low_cells:
        fail("CNV-high/CNV-low groups are empty; adjust quantile thresholds.")

    circ_cells = list(map(str, circ.obs_names))
    high_idx = [circ_cells.index(cell) for cell in high_cells if cell in circ_cells]
    low_idx = [circ_cells.index(cell) for cell in low_cells if cell in circ_cells]
    matrix = dense_matrix(circ)
    host_to_cols: dict[str, list[int]] = {}
    for col_idx, value in enumerate(circ.var["host_gene"]):
        for gene in split_host_genes(value):
            host_to_cols.setdefault(gene, []).append(col_idx)

    rows: list[dict[str, object]] = []
    for host_gene, cols in sorted(host_to_cols.items()):
        high_values = matrix[np.ix_(high_idx, cols)].sum(axis=1)
        low_values = matrix[np.ix_(low_idx, cols)].sum(axis=1)
        high_mean = float(np.mean(high_values)) if len(high_values) else np.nan
        low_mean = float(np.mean(low_values)) if len(low_values) else np.nan
        rows.append(
            {
                "host_gene": host_gene,
                "n_circRNAs": int(len(cols)),
                "cnv_metric": metric,
                "high_group_cutoff": high_cutoff,
                "low_group_cutoff": low_cutoff,
                "n_high_cells": int(len(high_idx)),
                "n_low_cells": int(len(low_idx)),
                "high_mean_support": high_mean,
                "low_mean_support": low_mean,
                "difference": high_mean - low_mean,
                "log2_fc_pseudocount1": float(np.log2((high_mean + 1.0) / (low_mean + 1.0))),
                "high_detected_cells": int(np.count_nonzero(high_values > 0)),
                "low_detected_cells": int(np.count_nonzero(low_values > 0)),
            }
        )
    return pd.DataFrame(rows).sort_values(["difference", "host_gene"], ascending=[False, True])


def local_cnv_at_circ_loci(circ, cnv, *, host_genes: set[str] | None = None) -> pd.DataFrame:
    from scipy import stats

    circ_coords = feature_coordinates(circ, kind="circ")
    cnv_coords = feature_coordinates(cnv, kind="cnv")
    if circ_coords is None or cnv_coords is None:
        return pd.DataFrame(
            [
                {
                    "status": "missing_coordinates",
                    "message": "circRNA or CNV genomic coordinates were not available.",
                }
            ]
        )
    circ_var = circ.var.copy().reset_index(drop=True)
    if "host_gene" not in circ_var.columns:
        circ_var["host_gene"] = ""
    circ_coords = circ_coords.reset_index(drop=True)
    cnv_coords = cnv_coords.reset_index(drop=True)
    cnv_matrix = dense_matrix(cnv)
    circ_matrix = dense_matrix(circ)
    shared = common_obs_names(circ, cnv)
    circ_cell_idx = [list(map(str, circ.obs_names)).index(cell) for cell in shared]
    cnv_cell_idx = [list(map(str, cnv.obs_names)).index(cell) for cell in shared]

    rows: list[dict[str, object]] = []
    for _, coord in circ_coords.iterrows():
        circ_idx = int(coord["feature_position"])
        genes = split_host_genes(circ_var.loc[circ_idx, "host_gene"])
        if host_genes and not (genes & host_genes):
            continue
        same_chrom = cnv_coords[cnv_coords["chrom"].astype(str) == str(coord["chrom"])]
        hits = same_chrom[(same_chrom["start"] <= coord["midpoint"]) & (same_chrom["end"] > coord["midpoint"])]
        if hits.empty:
            continue
        cnv_idx = int(hits.iloc[0]["feature_position"])
        circ_support = circ_matrix[circ_cell_idx, circ_idx].astype(float)
        local_cnv = cnv_matrix[cnv_cell_idx, cnv_idx].astype(float)
        if len(shared) >= 3 and len(np.unique(local_cnv)) >= 2 and len(np.unique(circ_support)) >= 2:
            pearson = stats.pearsonr(local_cnv, circ_support)
            pearson_r = float(pearson.statistic)
            pearson_p = float(pearson.pvalue)
        else:
            pearson_r = np.nan
            pearson_p = np.nan
        state_summary: dict[str, float] = {}
        for state in sorted(set(local_cnv.tolist())):
            state_summary[str(state)] = float(np.mean(circ_support[local_cnv == state]))
        rows.append(
            {
                "status": "ok",
                "circ_id": str(circ.var_names[circ_idx]),
                "host_gene": ";".join(sorted(genes)),
                "chrom": coord["chrom"],
                "start": int(coord["start"]),
                "end": int(coord["end"]),
                "cnv_bin_id": str(cnv.var_names[cnv_idx]),
                "n_cells": int(len(shared)),
                "mean_local_copy_number": float(np.mean(local_cnv)),
                "mean_circ_support": float(np.mean(circ_support)),
                "pearson_r": pearson_r,
                "pearson_p": pearson_p,
                "mean_support_by_copy_number": json.dumps(state_summary, sort_keys=True),
            }
        )
    if not rows:
        return pd.DataFrame(
            [
                {
                    "status": "no_locus_matches",
                    "message": "No circRNA loci overlapped available CNV bins.",
                }
            ]
        )
    return pd.DataFrame(rows)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Analyze IMR90 scRR CNV burden and circRNA programs.")
    parser.add_argument("--input", required=True, type=Path, help="IMR90 RNA+circ+CNV .h5mu file.")
    parser.add_argument("--correlations-out", required=True, type=Path, help="Output global correlation TSV/CSV.")
    parser.add_argument("--enrichment-out", required=True, type=Path, help="Output CNV-high vs CNV-low host-gene table.")
    parser.add_argument("--cell-metrics-out", type=Path, default=None, help="Optional per-cell metrics TSV/CSV.")
    parser.add_argument("--local-cnv-out", type=Path, default=None, help="Optional local CNV-at-circRNA-locus TSV/CSV.")
    parser.add_argument("--high-quantile", type=float, default=0.75, help="CNV-high quantile threshold.")
    parser.add_argument("--low-quantile", type=float, default=0.25, help="CNV-low quantile threshold.")
    parser.add_argument("--cnv-metric", default="frac_non_diploid", help="CNV burden metric for high/low grouping.")
    parser.add_argument(
        "--local-host-genes",
        default="COL1A1,FN1,VIM,COL6A3,P4HB",
        help="Comma-separated host genes for optional local CNV analysis; use empty string for all.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    mdata = load_mudata(args.input)
    circ = mdata.mod[find_modality(mdata, CIRC_MODALITY_CANDIDATES)]
    cnv = mdata.mod[find_modality(mdata, CNV_MODALITY_CANDIDATES)]
    cell_table = build_cell_table(mdata, circ, cnv)
    pairs = [
        ("frac_non_diploid", "circRNA_count"),
        ("frac_non_diploid", "circRNA_total_support"),
        ("frac_loss", "circRNA_count"),
        ("frac_gain", "circRNA_count"),
    ]
    write_table(pd.DataFrame(correlation_rows(cell_table, pairs)), args.correlations_out)
    enrichment = host_gene_enrichment(
        circ,
        cell_table,
        high_quantile=args.high_quantile,
        low_quantile=args.low_quantile,
        metric=args.cnv_metric,
    )
    write_table(enrichment, args.enrichment_out)
    if args.cell_metrics_out is not None:
        write_table(cell_table, args.cell_metrics_out)
    if args.local_cnv_out is not None:
        genes = {gene.strip() for gene in args.local_host_genes.split(",") if gene.strip()}
        write_table(local_cnv_at_circ_loci(circ, cnv, host_genes=genes or None), args.local_cnv_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
