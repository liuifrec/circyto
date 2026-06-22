#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from _mudata_utils import (
    CIRC_MODALITY_CANDIDATES,
    common_obs_names,
    dataset_name_from_path,
    find_modality,
    host_gene_recovery,
    host_gene_source_counts,
    json_dumps,
    load_mudata,
    write_table,
)


def summarize_one(path: Path, dataset_name: str | None = None) -> dict[str, object]:
    mdata = load_mudata(path)
    modalities = list(mdata.mod.keys())
    circ_name = find_modality(mdata, CIRC_MODALITY_CANDIDATES, required=False)
    modality_cells = {name: int(mdata.mod[name].n_obs) for name in modalities}
    modality_features = {name: int(mdata.mod[name].n_vars) for name in modalities}
    pairwise_overlap: dict[str, int] = {}
    for i, left in enumerate(modalities):
        for right in modalities[i + 1 :]:
            pairwise_overlap[f"{left}&{right}"] = len(common_obs_names(mdata.mod[left], mdata.mod[right]))
    all_overlap = len(common_obs_names(*(mdata.mod[name] for name in modalities))) if modalities else 0
    all_union = len(set().union(*(set(map(str, mdata.mod[name].obs_names)) for name in modalities))) if modalities else 0

    row: dict[str, object] = {
        "dataset": dataset_name or dataset_name_from_path(path),
        "path": str(path),
        "modalities": ";".join(map(str, modalities)),
        "modality_cells": json_dumps(modality_cells),
        "modality_features": json_dumps(modality_features),
        "all_modality_cell_overlap": int(all_overlap),
        "all_modality_cell_union": int(all_union),
        "pairwise_cell_overlap": json_dumps(pairwise_overlap),
        "circRNA_count": 0,
        "host_gene_annotated": 0,
        "host_gene_total": 0,
        "host_gene_recovery": 0.0,
        "host_gene_source_counts": "{}",
    }
    if circ_name is not None:
        circ = mdata.mod[circ_name]
        annotated, total, recovery = host_gene_recovery(circ)
        row.update(
            {
                "circRNA_count": int(circ.n_vars),
                "host_gene_annotated": annotated,
                "host_gene_total": total,
                "host_gene_recovery": recovery,
                "host_gene_source_counts": json_dumps(host_gene_source_counts(circ)),
            }
        )
    return row


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Summarize modality inventory and host-gene annotation recovery for one or more .h5mu files."
    )
    parser.add_argument("h5mu", nargs="+", type=Path, help="Input MuData .h5mu file(s).")
    parser.add_argument(
        "--dataset-name",
        action="append",
        default=None,
        help="Optional dataset label. Repeat once per input, in the same order.",
    )
    parser.add_argument("--out", type=Path, default=None, help="Output TSV/CSV path. Defaults to stdout.")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    if args.dataset_name and len(args.dataset_name) != len(args.h5mu):
        parser.error("--dataset-name must be repeated once per input .h5mu file.")
    rows = [
        summarize_one(path, args.dataset_name[i] if args.dataset_name else None)
        for i, path in enumerate(args.h5mu)
    ]
    write_table(pd.DataFrame(rows), args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
