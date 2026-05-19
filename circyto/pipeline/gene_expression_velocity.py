from __future__ import annotations

import os
import json
from pathlib import Path
from typing import Any

import pandas as pd


VALID_GENE_EXPRESSION_METHODS = {"none", "featurecounts", "velocyto"}
VALID_VELOCITY_LAYERS = {"none", "velocyto"}
REGENERABLE_ALIGNMENT_SUFFIXES = (
    ".sam",
    ".bam",
    ".bai",
    ".junction",
)
REGENERABLE_ALIGNMENT_NAMES = (
    "Aligned.out.sam",
    "Chimeric.out.junction",
    "Unmapped.out.mate1",
    "Unmapped.out.mate2",
    "bwa_rescue.sam",
)

VALID_VELOCITY_LAYOUT_FILES = {
    "barcodes.tsv",
    "features.tsv",
    "spliced.mtx",
    "unspliced.mtx",
}


def validate_full_length_future_options(
    *,
    export_mudata: bool,
    gene_counts: Path | None,
    gene_expression_method: str,
    velocity_layers: str,
    cleanup_intermediates: bool,
    dry_run: bool,
) -> None:
    if gene_expression_method not in VALID_GENE_EXPRESSION_METHODS:
        raise ValueError(
            f"Unsupported --gene-expression-method: {gene_expression_method}. "
            f"Choose from: {', '.join(sorted(VALID_GENE_EXPRESSION_METHODS))}"
        )
    if velocity_layers not in VALID_VELOCITY_LAYERS:
        raise ValueError(
            f"Unsupported --velocity-layers: {velocity_layers}. "
            f"Choose from: {', '.join(sorted(VALID_VELOCITY_LAYERS))}"
        )
    if velocity_layers == "velocyto" and gene_expression_method != "velocyto":
        raise ValueError("--velocity-layers velocyto requires --gene-expression-method velocyto")

    if gene_counts is not None and gene_expression_method == "velocyto":
        raise NotImplementedError(
            "--gene-counts TSV import currently supports only normalized gene-count table staging, "
            "not velocyto-based RNA layer import."
        )

    if export_mudata or gene_expression_method != "none" or velocity_layers != "none":
        raise NotImplementedError(
            "Gene-expression and velocity-compatible outputs are not implemented yet for "
            "`circyto workflow full-length-circrna`. Current behavior remains circ-only h5ad. "
            "See docs/gene_expression_velocity_integration.md for the staged design."
        )

    if cleanup_intermediates and not dry_run:
        raise NotImplementedError(
            "--cleanup-intermediates execution is not implemented yet. "
            "Use it only with --dry-run to inspect the planned cleanup policy."
        )


def validate_feature_id_uniqueness(feature_ids: list[str], *, label: str) -> None:
    normalized = [str(feature_id).strip() for feature_id in feature_ids]
    duplicates = sorted({feature_id for feature_id in normalized if feature_id and normalized.count(feature_id) > 1})
    if duplicates:
        raise ValueError(
            f"Duplicate {label} IDs detected: {', '.join(duplicates[:5])}. "
            f"{label} identifiers must be unique."
        )


def validate_cell_id_consistency(
    circ_cell_ids: list[str],
    rna_cell_ids: list[str],
    *,
    circ_label: str = "circ",
    rna_label: str = "rna",
) -> None:
    circ = [str(cell_id).strip() for cell_id in circ_cell_ids if str(cell_id).strip()]
    rna = [str(cell_id).strip() for cell_id in rna_cell_ids if str(cell_id).strip()]
    if not circ:
        raise ValueError(f"{circ_label} cell ID set is empty")
    if not rna:
        raise ValueError(f"{rna_label} cell ID set is empty")
    if circ != rna:
        only_circ = sorted(set(circ) - set(rna))
        only_rna = sorted(set(rna) - set(circ))
        raise ValueError(
            f"Cell ID mismatch between {circ_label} and {rna_label}. "
            f"Only in {circ_label}: {only_circ[:5]}. Only in {rna_label}: {only_rna[:5]}."
        )


def validate_gene_expression_table_schema(path: Path) -> dict[str, Any]:
    df = pd.read_csv(path, sep="\t", keep_default_na=False)
    required = ["gene_id", "gene_name"]
    missing = [column for column in required if column not in df.columns]
    if missing:
        raise ValueError(
            f"{path} is missing required gene-expression columns: {', '.join(missing)}. "
            "Expected at minimum: gene_id, gene_name, and one or more cell columns."
        )
    if df.empty:
        raise ValueError(f"{path} has no gene-expression rows")
    cell_columns = [column for column in df.columns if column not in {"gene_id", "gene_name"}]
    if not cell_columns:
        raise ValueError(f"{path} has no cell columns after gene_id and gene_name")
    validate_feature_id_uniqueness(df["gene_id"].astype(str).tolist(), label="gene")
    validate_feature_id_uniqueness([str(column) for column in cell_columns], label="cell")
    numeric = df[cell_columns].apply(pd.to_numeric, errors="coerce")
    if numeric.isna().any().any():
        raise ValueError(
            f"{path} contains non-numeric gene-count values in one or more cell columns."
        )
    return {
        "path": str(path.resolve()),
        "n_genes": int(df.shape[0]),
        "n_cells": int(len(cell_columns)),
        "cell_ids": [str(column) for column in cell_columns],
        "feature_id_column": "gene_id",
        "feature_name_column": "gene_name",
    }


def import_gene_counts_table(
    *,
    path: Path,
    expected_cell_ids: list[str],
    outdir: Path,
) -> dict[str, Any]:
    summary = validate_gene_expression_table_schema(path)
    df = pd.read_csv(path, sep="\t", keep_default_na=False)
    cell_columns = [column for column in df.columns if column not in {"gene_id", "gene_name"}]
    validate_cell_id_consistency(expected_cell_ids, [str(column) for column in cell_columns], circ_label="circ", rna_label="rna")

    numeric = df[cell_columns].apply(pd.to_numeric, errors="raise")
    normalized_counts = pd.concat([df[["gene_id", "gene_name"]].copy(), numeric], axis=1)
    feature_table = df[["gene_id", "gene_name"]].copy()

    outdir.mkdir(parents=True, exist_ok=True)
    gene_counts_out = outdir / "gene_counts.tsv"
    feature_out = outdir / "gene_feature_table.tsv"
    summary_out = outdir / "rna_import_summary.json"
    normalized_counts.to_csv(gene_counts_out, sep="\t", index=False)
    feature_table.to_csv(feature_out, sep="\t", index=False)

    payload = {
        **summary,
        "cell_ids": [str(column) for column in cell_columns],
        "output_gene_counts": str(gene_counts_out.resolve()),
        "output_gene_feature_table": str(feature_out.resolve()),
        "feature_table_columns": ["gene_id", "gene_name"],
        "count_table_columns": ["gene_id", "gene_name", *[str(column) for column in cell_columns]],
        "total_counts_sum": float(numeric.to_numpy().sum()),
    }
    summary_out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    payload["output_rna_import_summary"] = str(summary_out.resolve())
    return payload


def validate_velocity_layers_schema(path: Path) -> dict[str, Any]:
    if not path.exists() or not path.is_dir():
        raise FileNotFoundError(f"Velocity layer directory not found: {path}")
    names = {child.name for child in path.iterdir()}
    missing = sorted(VALID_VELOCITY_LAYOUT_FILES - names)
    if missing:
        raise ValueError(
            f"{path} is missing required velocity-layer files: {', '.join(missing)}. "
            "Expected at minimum: barcodes.tsv, features.tsv, spliced.mtx, unspliced.mtx."
        )
    barcodes = [
        line.strip()
        for line in (path / "barcodes.tsv").read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if not barcodes:
        raise ValueError(f"{path / 'barcodes.tsv'} has no cell IDs")
    feature_rows = [
        line.rstrip("\n").split("\t")
        for line in (path / "features.tsv").read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if not feature_rows:
        raise ValueError(f"{path / 'features.tsv'} has no feature rows")
    gene_ids = [row[0].strip() for row in feature_rows if row and row[0].strip()]
    if len(gene_ids) != len(feature_rows):
        raise ValueError(f"{path / 'features.tsv'} contains blank gene_id entries")
    validate_feature_id_uniqueness(gene_ids, label="gene")
    validate_feature_id_uniqueness(barcodes, label="cell")
    return {
        "path": str(path.resolve()),
        "n_cells": int(len(barcodes)),
        "n_genes": int(len(gene_ids)),
        "has_ambiguous": (path / "ambiguous.mtx").exists(),
        "cell_ids": barcodes,
        "feature_ids": gene_ids,
    }


def build_cleanup_plan(*, outdir: Path) -> dict[str, Any]:
    candidates: list[dict[str, Any]] = []
    preserved = [
        "workflow_summary.json",
        "detector summaries",
        "matrix/",
        "anndata/",
        "mudata/",
        "qc/",
        "manifests/",
        "logs/",
        "provenance",
    ]
    for subdir_name in ("align", "ciri3"):
        subdir = outdir / subdir_name
        if not subdir.exists():
            continue
        for root, _, filenames in os.walk(subdir):
            for filename in filenames:
                path = Path(root) / filename
                if (
                    path.suffix.lower() in REGENERABLE_ALIGNMENT_SUFFIXES
                    or filename in REGENERABLE_ALIGNMENT_NAMES
                ):
                    try:
                        size = int(path.stat().st_size)
                    except OSError:
                        size = 0
                    candidates.append({"path": str(path.resolve()), "bytes": size})
    candidates.sort(key=lambda item: (-int(item["bytes"]), str(item["path"])))
    return {
        "enabled": True,
        "mode": "planned-only",
        "would_delete_only_regenerable_alignment_intermediates": True,
        "preserve": preserved,
        "candidate_count": len(candidates),
        "candidate_bytes": int(sum(int(item["bytes"]) for item in candidates)),
        "delete_candidates": candidates[:20],
    }
