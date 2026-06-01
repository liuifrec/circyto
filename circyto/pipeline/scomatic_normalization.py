from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Sequence

import pandas as pd

from circyto.pipeline.workflow_reporting import (
    apply_standard_provenance,
    utc_now_iso,
    write_json,
)


SCOMATIC_CANDIDATE_COLUMNS = [
    "variant_id",
    "cell_id",
    "chrom",
    "pos",
    "ref",
    "alt",
    "gene",
    "filter_status",
    "candidate_variant_class",
    "read_support",
    "vaf",
    "caller",
]

TERMINOLOGY_NOTE = (
    "SComatic outputs are treated as RNA-derived candidate variant signals, "
    "not validated somatic mutations."
)

_REQUIRED_FIELD_ALIASES = {
    "cell_id": ("cell_id", "cell", "cellid", "barcode", "sample", "sample_id", "cb", "index"),
    "chrom": ("chrom", "chr", "chromosome", "#chrom"),
    "pos": ("pos", "position", "start"),
    "ref": ("ref", "reference", "ref_allele", "reference_allele"),
    "alt": ("alt", "alternative", "alt_allele", "alternate_allele"),
    "gene": ("gene", "gene_name", "host_gene"),
    "filter_status": ("filter_status", "filter", "scomatic_filter", "status"),
    "read_support": (
        "read_support",
        "support",
        "alt_reads",
        "alt_count",
        "read_count_alt",
        "dp_alt",
        "n_alt",
        "alt_depth",
    ),
    "vaf": ("vaf", "af", "allele_fraction", "allele_frequency", "variant_allele_fraction"),
}

_OPTIONAL_FIELD_ALIASES = {
    "variant_id": ("variant_id", "id", "variant", "mutation_id"),
    "candidate_variant_class": ("candidate_variant_class", "variant_class", "class"),
    "caller": ("caller",),
}

_DEFAULT_VALUES = {
    "candidate_variant_class": "RNA-derived candidate variant signal",
    "caller": "SComatic",
}


def _normalize_column_name(value: object) -> str:
    return str(value).strip().lower().replace(" ", "_").replace("-", "_")


def _column_lookup(df: pd.DataFrame) -> dict[str, str]:
    lookup: dict[str, str] = {}
    for column in df.columns:
        lookup.setdefault(_normalize_column_name(column), str(column))
    return lookup


def _find_column(df: pd.DataFrame, aliases: Sequence[str]) -> str | None:
    lookup = _column_lookup(df)
    for alias in aliases:
        match = lookup.get(_normalize_column_name(alias))
        if match is not None:
            return match
    return None


def _read_delimited_table(path: Path) -> pd.DataFrame:
    sep = "," if path.suffix.lower() == ".csv" else "\t"
    try:
        return pd.read_csv(path, sep=sep, keep_default_na=False)
    except Exception as exc:
        raise ValueError(f"Could not read SComatic output table {path}: {exc}") from exc


def _resolve_required_columns(path: Path, df: pd.DataFrame) -> dict[str, str]:
    resolved: dict[str, str] = {}
    missing: list[str] = []
    for field, aliases in _REQUIRED_FIELD_ALIASES.items():
        column = _find_column(df, aliases)
        if column is None:
            missing.append(field)
        else:
            resolved[field] = column
    if missing:
        available = ", ".join(str(column) for column in df.columns) or "<none>"
        raise ValueError(
            f"Could not normalize SComatic output {path}: missing required SComatic "
            f"columns for fields: {', '.join(missing)}. Available columns: {available}"
        )
    return resolved


def _numeric_series(path: Path, df: pd.DataFrame, column: str, field: str) -> pd.Series:
    try:
        return pd.to_numeric(df[column], errors="raise")
    except Exception as exc:
        raise ValueError(
            f"Could not normalize SComatic output {path}: column {column!r} for "
            f"field {field!r} must be numeric."
        ) from exc


def _integer_series(path: Path, df: pd.DataFrame, column: str, field: str) -> pd.Series:
    values = _numeric_series(path, df, column, field)
    non_integer = values.dropna().map(float).map(lambda value: not value.is_integer())
    if bool(non_integer.any()):
        raise ValueError(
            f"Could not normalize SComatic output {path}: column {column!r} for "
            f"field {field!r} must contain integer values."
        )
    return values.astype(int)


def _string_series(df: pd.DataFrame, column: str) -> pd.Series:
    return df[column].astype(str).str.strip()


def _generated_variant_ids(normalized: pd.DataFrame) -> pd.Series:
    return (
        normalized["chrom"].astype(str)
        + ":"
        + normalized["pos"].astype(str)
        + ":"
        + normalized["ref"].astype(str)
        + ">"
        + normalized["alt"].astype(str)
        + ":"
        + normalized["cell_id"].astype(str)
    )


def normalize_scomatic_output_table(path: Path) -> tuple[pd.DataFrame, dict[str, Any]]:
    raw = _read_delimited_table(path)
    if raw.empty:
        raise ValueError(f"Could not normalize SComatic output {path}: table is empty.")

    required_columns = _resolve_required_columns(path, raw)
    normalized = pd.DataFrame(
        {
            "cell_id": _string_series(raw, required_columns["cell_id"]),
            "chrom": _string_series(raw, required_columns["chrom"]),
            "pos": _integer_series(path, raw, required_columns["pos"], "pos"),
            "ref": _string_series(raw, required_columns["ref"]),
            "alt": _string_series(raw, required_columns["alt"]),
            "gene": _string_series(raw, required_columns["gene"]),
            "filter_status": _string_series(raw, required_columns["filter_status"]),
            "read_support": _integer_series(path, raw, required_columns["read_support"], "read_support"),
            "vaf": _numeric_series(path, raw, required_columns["vaf"], "vaf").astype(float),
        }
    )

    for field in ("cell_id", "chrom", "ref", "alt"):
        if bool(normalized[field].astype(str).str.strip().eq("").any()):
            raise ValueError(
                f"Could not normalize SComatic output {path}: required field "
                f"{field!r} contains empty values."
            )

    optional_columns: dict[str, str | None] = {}
    for field, aliases in _OPTIONAL_FIELD_ALIASES.items():
        optional_columns[field] = _find_column(raw, aliases)

    generated_variant_ids = _generated_variant_ids(normalized)
    variant_id_col = optional_columns["variant_id"]
    if variant_id_col is None:
        normalized["variant_id"] = generated_variant_ids
    else:
        variant_ids = _string_series(raw, variant_id_col)
        normalized["variant_id"] = variant_ids.where(variant_ids.ne(""), generated_variant_ids)

    for field in ("candidate_variant_class", "caller"):
        column = optional_columns[field]
        default = _DEFAULT_VALUES[field]
        if column is None:
            normalized[field] = default
        else:
            values = _string_series(raw, column)
            normalized[field] = values.where(values.ne(""), default)

    summary = {
        "path": str(path.resolve()),
        "rows": int(raw.shape[0]),
        "columns": [str(column) for column in raw.columns],
        "resolved_required_columns": required_columns,
        "resolved_optional_columns": optional_columns,
    }
    return normalized[SCOMATIC_CANDIDATE_COLUMNS], summary


def _load_cell_annotation(path: Path, candidate_cells: set[str]) -> dict[str, Any]:
    df = _read_delimited_table(path)
    cell_col = _find_column(df, ("Index", "cell_id", "cell", "CB"))
    cell_type_col = _find_column(df, ("Cell_type", "cell_type"))
    if cell_col is None or cell_type_col is None:
        available = ", ".join(str(column) for column in df.columns) or "<none>"
        raise ValueError(
            f"Cell annotation table {path} must contain Index and Cell_type columns. "
            f"Available columns: {available}"
        )
    cell_ids = df[cell_col].astype(str).str.strip().tolist()
    duplicate_cells = sorted({cell for cell in cell_ids if cell_ids.count(cell) > 1})
    if duplicate_cells:
        preview = ", ".join(duplicate_cells[:10])
        raise ValueError(f"Cell annotation table {path} contains duplicate cell IDs: {preview}")
    annotation_cells = set(cell_ids)
    missing_cells = sorted(candidate_cells - annotation_cells)
    if missing_cells:
        preview = ", ".join(missing_cells[:10])
        raise ValueError(
            f"Cell annotation table {path} is missing SComatic candidate cell IDs: {preview}"
        )
    return {
        "path": str(path.resolve()),
        "rows": int(df.shape[0]),
        "cell_column": cell_col,
        "cell_type_column": cell_type_col,
        "n_candidate_cells_present": int(len(candidate_cells)),
    }


def _load_provenance_metadata(path: Path | None) -> dict[str, Any] | None:
    if path is None:
        return None
    text = path.read_text(encoding="utf-8")
    payload: dict[str, Any] = {
        "path": str(path.resolve()),
    }
    try:
        payload["format"] = "json"
        payload["content"] = json.loads(text) if text.strip() else None
    except json.JSONDecodeError:
        payload["format"] = "text"
        payload["content"] = text
    return payload


def normalize_scomatic_results(
    *,
    scomatic_output_paths: Sequence[Path],
    outdir: Path,
    cell_annotation_path: Path | None = None,
    provenance_metadata_path: Path | None = None,
) -> dict[str, Any]:
    if not scomatic_output_paths:
        raise ValueError("At least one --scomatic-output file is required.")

    normalized_tables: list[pd.DataFrame] = []
    source_summaries: list[dict[str, Any]] = []
    for path in scomatic_output_paths:
        if not path.exists():
            raise FileNotFoundError(f"SComatic output file not found: {path}")
        normalized, source_summary = normalize_scomatic_output_table(path)
        normalized_tables.append(normalized)
        source_summaries.append(source_summary)

    candidate_df = pd.concat(normalized_tables, ignore_index=True)
    candidate_cells = set(candidate_df["cell_id"].astype(str))
    cell_annotation = (
        _load_cell_annotation(cell_annotation_path, candidate_cells)
        if cell_annotation_path is not None
        else None
    )
    provenance_metadata = _load_provenance_metadata(provenance_metadata_path)

    outdir.mkdir(parents=True, exist_ok=True)
    candidate_out = outdir / "scomatic_candidate_summary.tsv"
    candidate_df.to_csv(candidate_out, sep="\t", index=False)

    summary_path = outdir / "normalize_scomatic_results_summary.json"
    summary = {
        "mode": "read-only-normalization",
        "description": "Normalized external SComatic output files; SComatic was not run.",
        "n_input_files": int(len(scomatic_output_paths)),
        "n_input_rows": int(sum(source["rows"] for source in source_summaries)),
        "n_candidates": int(candidate_df.shape[0]),
        "n_cells": int(candidate_df["cell_id"].nunique()),
        "cell_ids": sorted(candidate_cells),
        "source_files": source_summaries,
        "cell_annotation": cell_annotation,
        "provenance_metadata": provenance_metadata,
        "schema": SCOMATIC_CANDIDATE_COLUMNS,
        "terminology_note": TERMINOLOGY_NOTE,
        "output_scomatic_candidate_summary": str(candidate_out.resolve()),
        "output_normalize_scomatic_results_summary": str(summary_path.resolve()),
    }
    summary = apply_standard_provenance(
        summary,
        command_name="circyto normalize-scomatic-results",
        workflow_type="scomatic-normalization",
        protocol="external-scomatic",
        read_layout="unknown",
        genome_fasta=None,
        gtf=None,
        detector_backend="SComatic",
        started_at=utc_now_iso(),
        completed_at=utc_now_iso(),
        elapsed_seconds=0.0,
    )
    write_json(summary_path, summary)
    return summary
