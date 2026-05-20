#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import pandas as pd


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _read_cell_index(path: Path) -> list[str]:
    if not path.exists():
        return []
    return [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def build_summary(workdir: Path, *, top_n: int = 10) -> dict[str, Any]:
    rna_dir = workdir / "rna"
    gene_counts_path = rna_dir / "gene_counts.tsv"
    feature_table_path = rna_dir / "gene_feature_table.tsv"
    rna_summary_path = rna_dir / "rna_import_summary.json"
    cell_index_path = workdir / "matrix" / "cell_index.txt"

    gene_counts_exists = gene_counts_path.exists()
    feature_table_exists = feature_table_path.exists()
    rna_summary_exists = rna_summary_path.exists()

    df = pd.read_csv(gene_counts_path, sep="\t", keep_default_na=False) if gene_counts_exists else None
    rna_summary = _load_json(rna_summary_path) if rna_summary_exists else {}

    cell_columns: list[str] = []
    top_genes: list[dict[str, Any]] = []
    cell_total_counts: list[dict[str, Any]] = []
    cell_id_match: dict[str, Any] | None = None

    if df is not None:
        cell_columns = [column for column in df.columns if column not in {"gene_id", "gene_name"}]
        numeric = df[cell_columns].apply(pd.to_numeric, errors="coerce") if cell_columns else pd.DataFrame()
        if cell_columns:
            gene_totals = numeric.sum(axis=1)
            ranked_gene_idx = gene_totals.sort_values(ascending=False).index.tolist()[:top_n]
            top_genes = [
                {
                    "gene_id": str(df.iloc[idx]["gene_id"]),
                    "gene_name": str(df.iloc[idx]["gene_name"]),
                    "total_count": int(gene_totals.iloc[idx]),
                }
                for idx in ranked_gene_idx
            ]

            totals_by_cell = {str(column): int(numeric[column].sum()) for column in cell_columns}
            ranked_cells = sorted(totals_by_cell.items(), key=lambda item: (item[1], item[0]))
            cell_total_counts = [{"cell_id": cell_id, "total_count": total} for cell_id, total in ranked_cells]

            matrix_cell_ids = _read_cell_index(cell_index_path)
            if matrix_cell_ids:
                matches = cell_columns == matrix_cell_ids
                cell_id_match = {
                    "matrix_cell_index_present": True,
                    "matches": matches,
                    "rna_only": sorted(set(cell_columns) - set(matrix_cell_ids)),
                    "matrix_only": sorted(set(matrix_cell_ids) - set(cell_columns)),
                }
            else:
                cell_id_match = {
                    "matrix_cell_index_present": False,
                    "matches": None,
                    "rna_only": list(cell_columns),
                    "matrix_only": [],
                }

    return {
        "workdir": str(workdir.resolve()),
        "gene_counts_exists": gene_counts_exists,
        "gene_feature_table_exists": feature_table_exists,
        "rna_import_summary_exists": rna_summary_exists,
        "n_genes": int(rna_summary.get("n_genes", df.shape[0] if df is not None else 0) or 0),
        "n_cells": int(rna_summary.get("n_cells", len(cell_columns)) or 0),
        "total_counts_sum": int(rna_summary.get("total_counts_sum", 0) or 0),
        "assigned_templates": int(rna_summary.get("assigned_templates", 0) or 0),
        "ambiguous_templates_excluded": int(rna_summary.get("ambiguous_templates_excluded", 0) or 0),
        "unassigned_templates": int(rna_summary.get("unassigned_templates", 0) or 0),
        "top_expressed_genes": top_genes,
        "lowest_total_rna_cells": cell_total_counts[:top_n],
        "highest_total_rna_cells": list(reversed(cell_total_counts[-top_n:])),
        "cell_ids_match_matrix": cell_id_match,
    }


def _print_text(summary: dict[str, Any]) -> None:
    print(f"Workdir: {summary['workdir']}")
    print(
        f"RNA outputs: gene_counts={summary['gene_counts_exists']} "
        f"gene_feature_table={summary['gene_feature_table_exists']} "
        f"rna_import_summary={summary['rna_import_summary_exists']}"
    )
    print(
        f"Counts: n_genes={summary['n_genes']} n_cells={summary['n_cells']} "
        f"total_counts_sum={summary['total_counts_sum']}"
    )
    print(
        f"Templates: assigned={summary['assigned_templates']} "
        f"ambiguous={summary['ambiguous_templates_excluded']} "
        f"unassigned={summary['unassigned_templates']}"
    )
    print("")
    print("Top expressed genes:")
    if summary["top_expressed_genes"]:
        for item in summary["top_expressed_genes"]:
            print(f"- {item['gene_id']} {item['gene_name']}: total_count={item['total_count']}")
    else:
        print("- none")
    print("")
    print("Lowest total RNA count cells:")
    if summary["lowest_total_rna_cells"]:
        for item in summary["lowest_total_rna_cells"]:
            print(f"- {item['cell_id']}: total_count={item['total_count']}")
    else:
        print("- none")
    print("")
    print("Highest total RNA count cells:")
    if summary["highest_total_rna_cells"]:
        for item in summary["highest_total_rna_cells"]:
            print(f"- {item['cell_id']}: total_count={item['total_count']}")
    else:
        print("- none")
    print("")
    print(f"Cell IDs match matrix/cell_index.txt: {summary['cell_ids_match_matrix']}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Read-only inspection of add-rna-profile outputs.")
    parser.add_argument("--workdir", required=True, type=Path, help="Completed workflow directory containing rna/ outputs.")
    parser.add_argument("--json", action="store_true", help="Emit JSON instead of human-readable text.")
    args = parser.parse_args(argv)

    try:
        summary = build_summary(args.workdir)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    if args.json:
        print(json.dumps(summary, indent=2, sort_keys=True))
    else:
        _print_text(summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
