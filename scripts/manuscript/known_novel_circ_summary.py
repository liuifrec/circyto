#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from _mudata_utils import (
    CIRC_MODALITY_CANDIDATES,
    dataset_name_from_path,
    dense_matrix,
    find_modality,
    first_existing_column,
    host_gene_recovery,
    load_mudata,
    nonempty_string_mask,
    write_table,
)


KNOWN_STATUS_COLUMNS = [
    "known_status",
    "annotation_status",
    "circ_annotation_status",
    "known_novel",
    "novelty",
    "status",
]
KNOWN_ID_COLUMNS = [
    "known_circ_id",
    "known_circRNA_id",
    "circatlas_id",
    "circAtlas_id",
    "annotation_id",
    "matched_id",
]


def infer_known_status(circ) -> pd.Series:
    var = circ.var.copy()
    bool_col = first_existing_column(var, ["is_known", "known"])
    if bool_col is not None:
        values = var[bool_col]
        if values.dtype == bool:
            return values.map({True: "known", False: "novel"})
        lowered = values.fillna("").astype(str).str.lower()
        return lowered.isin({"true", "1", "yes", "known"}).map({True: "known", False: "novel"})

    status_col = first_existing_column(var, KNOWN_STATUS_COLUMNS)
    if status_col is not None:
        lowered = var[status_col].fillna("").astype(str).str.lower()
        status = pd.Series("unknown", index=var.index, dtype=object)
        status[lowered.str.contains("novel", na=False)] = "novel"
        status[lowered.str.contains("known|match|exact|database", na=False)] = "known"
        return status

    id_col = first_existing_column(var, KNOWN_ID_COLUMNS)
    if id_col is not None:
        return nonempty_string_mask(var[id_col]).map({True: "known", False: "novel"})

    return pd.Series("unknown", index=var.index, dtype=object)


def summarize_one(path: Path, dataset_name: str | None = None) -> pd.DataFrame:
    mdata = load_mudata(path)
    circ = mdata.mod[find_modality(mdata, CIRC_MODALITY_CANDIDATES)]
    status = infer_known_status(circ).reset_index(drop=True)
    support = dense_matrix(circ).sum(axis=0).astype(float)
    has_host = (
        nonempty_string_mask(circ.var["host_gene"]).to_numpy()
        if "host_gene" in circ.var.columns
        else np.zeros(circ.n_vars, dtype=bool)
    )
    rows = []
    for label, indices in status.groupby(status).groups.items():
        idx = np.asarray(list(indices), dtype=int)
        values = support[idx]
        rows.append(
            {
                "dataset": dataset_name or dataset_name_from_path(path),
                "known_status": label,
                "n_circRNAs": int(len(idx)),
                "support_min": float(np.min(values)) if len(values) else np.nan,
                "support_q25": float(np.quantile(values, 0.25)) if len(values) else np.nan,
                "support_median": float(np.median(values)) if len(values) else np.nan,
                "support_mean": float(np.mean(values)) if len(values) else np.nan,
                "support_q75": float(np.quantile(values, 0.75)) if len(values) else np.nan,
                "support_max": float(np.max(values)) if len(values) else np.nan,
                "host_gene_annotated": int(has_host[idx].sum()) if len(idx) else 0,
                "host_gene_recovery": float(has_host[idx].mean()) if len(idx) else 0.0,
            }
        )
    annotated, total, recovery = host_gene_recovery(circ)
    rows.append(
        {
            "dataset": dataset_name or dataset_name_from_path(path),
            "known_status": "all",
            "n_circRNAs": int(total),
            "support_min": float(np.min(support)) if total else np.nan,
            "support_q25": float(np.quantile(support, 0.25)) if total else np.nan,
            "support_median": float(np.median(support)) if total else np.nan,
            "support_mean": float(np.mean(support)) if total else np.nan,
            "support_q75": float(np.quantile(support, 0.75)) if total else np.nan,
            "support_max": float(np.max(support)) if total else np.nan,
            "host_gene_annotated": annotated,
            "host_gene_recovery": recovery,
        }
    )
    return pd.DataFrame(rows)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Summarize known versus novel circRNAs in one or more .h5mu files.")
    parser.add_argument("h5mu", nargs="+", type=Path, help="Input MuData .h5mu file(s).")
    parser.add_argument(
        "--dataset-name",
        action="append",
        default=None,
        help="Optional dataset label. Repeat once per input, in the same order.",
    )
    parser.add_argument("--out", required=True, type=Path, help="Output TSV/CSV path.")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    if args.dataset_name and len(args.dataset_name) != len(args.h5mu):
        raise SystemExit("[ERROR] --dataset-name must be repeated once per input .h5mu file.")
    frames = [
        summarize_one(path, args.dataset_name[i] if args.dataset_name else None)
        for i, path in enumerate(args.h5mu)
    ]
    write_table(pd.concat(frames, ignore_index=True), args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
