from __future__ import annotations

from dataclasses import dataclass
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
    "RNA-derived candidate variant signals without orthogonal DNA validation."
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

_SCOMATIC_CALLING_REQUIRED_ALIASES = {
    "chrom": ("#CHROM", "CHROM", "chrom"),
    "pos": ("Start", "POS", "Position"),
    "ref": ("REF",),
    "alt": ("ALT",),
    "filter_status": ("FILTER",),
    "cell_types": ("Cell_types", "Cell_type", "cell_types", "cell_type"),
    "read_support": ("Bc", "Read_count_ALT", "read_support"),
    "vaf": ("VAF",),
}

_SCOMATIC_CALLING_OPTIONAL_ALIASES = {
    "gene": ("Gene", "gene", "gene_name", "host_gene"),
    "cell_type_filter": ("Cell_type_Filter", "cell_type_filter"),
}

_MISSING_VALUES = {"", ".", "NA", "NaN", "nan", "None", "none"}


@dataclass(frozen=True)
class _TablePreamble:
    sep: str
    skiprows: int
    metadata_lines: list[str]
    header_columns: list[str]
    header_line_number: int


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


def _inspect_table_preamble(path: Path) -> _TablePreamble:
    sep = "," if path.suffix.lower() == ".csv" else "\t"
    metadata_lines: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle):
            stripped = line.rstrip("\n\r")
            if not stripped:
                continue
            if sep == "\t" and stripped.startswith("##"):
                metadata_lines.append(stripped)
                continue
            return _TablePreamble(
                sep=sep,
                skiprows=line_number,
                metadata_lines=metadata_lines,
                header_columns=stripped.split(sep),
                header_line_number=line_number + 1,
            )
    return _TablePreamble(
        sep=sep,
        skiprows=0,
        metadata_lines=metadata_lines,
        header_columns=[],
        header_line_number=0,
    )


def _read_delimited_table(path: Path, preamble: _TablePreamble | None = None) -> pd.DataFrame:
    preamble = preamble or _inspect_table_preamble(path)
    try:
        return pd.read_csv(path, sep=preamble.sep, keep_default_na=False, skiprows=preamble.skiprows)
    except Exception as exc:
        raise ValueError(f"Could not read SComatic output table {path}: {exc}") from exc


def _header_lookup(columns: Sequence[str]) -> dict[str, str]:
    lookup: dict[str, str] = {}
    for column in columns:
        lookup.setdefault(_normalize_column_name(column), str(column))
    return lookup


def _find_header_column(columns: Sequence[str], aliases: Sequence[str]) -> str | None:
    lookup = _header_lookup(columns)
    for alias in aliases:
        match = lookup.get(_normalize_column_name(alias))
        if match is not None:
            return match
    return None


def _is_scomatic_calling_header(columns: Sequence[str]) -> bool:
    if not columns:
        return False
    return all(
        _find_header_column(columns, aliases) is not None
        for aliases in _SCOMATIC_CALLING_REQUIRED_ALIASES.values()
    )


def _infer_scomatic_calling_stage(path: Path) -> str:
    name = path.name.lower()
    if "step1" in name:
        return "step1"
    if "step2" in name:
        return "step2"
    return "calling"


def _is_missing(value: object) -> bool:
    return str(value).strip() in _MISSING_VALUES


def _split_scomatic_groups(value: object, *, delimiter: str = ",") -> list[str]:
    if _is_missing(value):
        return []
    return [part.strip() for part in str(value).split(delimiter)]


def _value_for_group(path: Path, values: list[str], index: int, field: str) -> str:
    if index < len(values):
        return values[index]
    if len(values) == 1:
        return values[0]
    raise ValueError(
        f"Could not normalize SComatic calling table {path}: field {field!r} "
        f"has {len(values)} comma-separated values but expected value {index + 1}."
    )


def _numeric_text(path: Path, value: str, field: str) -> float:
    try:
        return float(value)
    except Exception as exc:
        raise ValueError(
            f"Could not normalize SComatic calling table {path}: value {value!r} "
            f"for field {field!r} must be numeric."
        ) from exc


def _integer_text(path: Path, value: str, field: str) -> int:
    numeric = _numeric_text(path, value, field)
    if not numeric.is_integer():
        raise ValueError(
            f"Could not normalize SComatic calling table {path}: value {value!r} "
            f"for field {field!r} must contain an integer value."
        )
    return int(numeric)


def _read_scomatic_calling_header(path: Path, preamble: _TablePreamble) -> tuple[dict[str, str], dict[str, str | None]]:
    required: dict[str, str] = {}
    missing: list[str] = []
    for field, aliases in _SCOMATIC_CALLING_REQUIRED_ALIASES.items():
        column = _find_header_column(preamble.header_columns, aliases)
        if column is None:
            missing.append(field)
        else:
            required[field] = column
    if missing:
        available = ", ".join(preamble.header_columns) or "<none>"
        raise ValueError(
            f"Could not normalize SComatic calling table {path}: missing required "
            f"SComatic step columns for fields: {', '.join(missing)}. "
            f"Available columns: {available}"
        )
    optional = {
        field: _find_header_column(preamble.header_columns, aliases)
        for field, aliases in _SCOMATIC_CALLING_OPTIONAL_ALIASES.items()
    }
    return required, optional


def _normalize_scomatic_calling_row(
    *,
    path: Path,
    row: dict[str, str],
    required_columns: dict[str, str],
    optional_columns: dict[str, str | None],
    caller: str,
) -> list[dict[str, Any]]:
    alt_groups = _split_scomatic_groups(row[required_columns["alt"]])
    cell_type_groups = _split_scomatic_groups(row[required_columns["cell_types"]])
    if not alt_groups or not cell_type_groups:
        return []

    read_support_groups = _split_scomatic_groups(row[required_columns["read_support"]])
    vaf_groups = _split_scomatic_groups(row[required_columns["vaf"]])
    if not read_support_groups or not vaf_groups:
        return []

    cell_type_filter_col = optional_columns.get("cell_type_filter")
    cell_type_filter_groups = (
        _split_scomatic_groups(row[cell_type_filter_col])
        if cell_type_filter_col is not None
        else []
    )
    gene_col = optional_columns.get("gene")
    gene = str(row[gene_col]).strip() if gene_col is not None else "."
    if _is_missing(gene):
        gene = "."

    chrom = str(row[required_columns["chrom"]]).strip()
    pos = _integer_text(path, str(row[required_columns["pos"]]).strip(), "pos")
    ref = str(row[required_columns["ref"]]).strip()
    row_filter = str(row[required_columns["filter_status"]]).strip()
    if _is_missing(row_filter):
        row_filter = "."

    records: list[dict[str, Any]] = []
    for cell_type_index, cell_type in enumerate(cell_type_groups):
        if _is_missing(cell_type):
            continue
        alt_group = _value_for_group(path, alt_groups, cell_type_index, "ALT")
        if _is_missing(alt_group):
            continue
        read_support_group = _value_for_group(path, read_support_groups, cell_type_index, "Bc")
        vaf_group = _value_for_group(path, vaf_groups, cell_type_index, "VAF")
        cell_type_filter = (
            _value_for_group(path, cell_type_filter_groups, cell_type_index, "Cell_type_Filter")
            if cell_type_filter_groups
            else "."
        )

        alts = _split_scomatic_groups(alt_group, delimiter="|")
        read_support_values = _split_scomatic_groups(read_support_group, delimiter="|")
        vaf_values = _split_scomatic_groups(vaf_group, delimiter="|")
        for alt_index, alt in enumerate(alts):
            if _is_missing(alt):
                continue
            read_support = _value_for_group(path, read_support_values, alt_index, "Bc")
            vaf = _value_for_group(path, vaf_values, alt_index, "VAF")
            filter_status = row_filter
            if filter_status == "." and not _is_missing(cell_type_filter):
                filter_status = cell_type_filter
            records.append(
                {
                    "variant_id": f"{chrom}:{pos}:{ref}>{alt}:{cell_type}",
                    "cell_id": cell_type,
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "gene": gene,
                    "filter_status": filter_status,
                    "candidate_variant_class": _DEFAULT_VALUES["candidate_variant_class"],
                    "read_support": _integer_text(path, read_support, "read_support"),
                    "vaf": _numeric_text(path, vaf, "vaf"),
                    "caller": caller,
                }
            )
    return records


def normalize_scomatic_calling_table(
    path: Path,
    *,
    preamble: _TablePreamble | None = None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    preamble = preamble or _inspect_table_preamble(path)
    required_columns, optional_columns = _read_scomatic_calling_header(path, preamble)
    header = preamble.header_columns
    stage = _infer_scomatic_calling_stage(path)
    caller = f"SComatic-{stage}"

    records: list[dict[str, Any]] = []
    data_rows = 0
    candidate_rows = 0
    with path.open("r", encoding="utf-8") as handle:
        for _ in range(preamble.header_line_number):
            next(handle, None)
        for line_number, line in enumerate(handle, start=preamble.header_line_number + 1):
            stripped = line.rstrip("\n\r")
            if not stripped:
                continue
            elements = stripped.split(preamble.sep)
            if len(elements) != len(header):
                raise ValueError(
                    f"Could not normalize SComatic calling table {path}: line {line_number} "
                    f"has {len(elements)} columns but header has {len(header)} columns."
                )
            data_rows += 1
            row = dict(zip(header, elements))
            row_records = _normalize_scomatic_calling_row(
                path=path,
                row=row,
                required_columns=required_columns,
                optional_columns=optional_columns,
                caller=caller,
            )
            if row_records:
                candidate_rows += 1
                records.extend(row_records)

    normalized = pd.DataFrame.from_records(records, columns=SCOMATIC_CANDIDATE_COLUMNS)
    summary = {
        "path": str(path.resolve()),
        "rows": int(data_rows),
        "candidate_rows": int(candidate_rows),
        "normalized_candidate_rows": int(normalized.shape[0]),
        "columns": [str(column) for column in header],
        "parser": "scomatic-calling-step",
        "scomatic_stage": stage,
        "metadata_line_count": int(len(preamble.metadata_lines)),
        "metadata_lines": preamble.metadata_lines,
        "header_line_number": int(preamble.header_line_number),
        "resolved_required_columns": required_columns,
        "resolved_optional_columns": optional_columns,
        "cell_id_semantics": (
            "Real SComatic step1/step2 outputs are cell-type-level calls; "
            "cell_id is populated from the SComatic Cell_types column."
        ),
    }
    return normalized, summary


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
    preamble = _inspect_table_preamble(path)
    if _is_scomatic_calling_header(preamble.header_columns):
        return normalize_scomatic_calling_table(path, preamble=preamble)

    raw = _read_delimited_table(path, preamble=preamble)
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
        "parser": "generic-scomatic-candidate-table",
        "metadata_line_count": int(len(preamble.metadata_lines)),
        "metadata_lines": preamble.metadata_lines,
        "header_line_number": int(preamble.header_line_number),
        "resolved_required_columns": required_columns,
        "resolved_optional_columns": optional_columns,
    }
    return normalized[SCOMATIC_CANDIDATE_COLUMNS], summary


def _load_cell_annotation(
    path: Path,
    candidate_cells: set[str],
    *,
    candidate_cell_semantics: str = "cell_id",
) -> dict[str, Any]:
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

    if candidate_cell_semantics == "cell_type":
        validation_col = cell_type_col
        validation_values = set(df[cell_type_col].astype(str).str.strip())
        missing_label = "SComatic candidate cell types"
    else:
        validation_col = cell_col
        validation_values = set(cell_ids)
        missing_label = "SComatic candidate cell IDs"

    missing_cells = sorted(candidate_cells - validation_values)
    if missing_cells:
        preview = ", ".join(missing_cells[:10])
        raise ValueError(
            f"Cell annotation table {path} is missing {missing_label}: {preview}"
        )
    return {
        "path": str(path.resolve()),
        "rows": int(df.shape[0]),
        "cell_column": cell_col,
        "cell_type_column": cell_type_col,
        "validated_against": validation_col,
        "candidate_cell_semantics": candidate_cell_semantics,
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
    parser_types = {str(source.get("parser", "")) for source in source_summaries}
    candidate_cell_semantics = (
        "cell_type"
        if parser_types and parser_types <= {"scomatic-calling-step"}
        else "cell_id"
    )
    cell_annotation = (
        _load_cell_annotation(
            cell_annotation_path,
            candidate_cells,
            candidate_cell_semantics=candidate_cell_semantics,
        )
        if cell_annotation_path is not None
        else None
    )
    provenance_metadata = _load_provenance_metadata(provenance_metadata_path)

    outdir.mkdir(parents=True, exist_ok=True)
    candidate_out = outdir / "scomatic_candidate_summary.tsv"
    candidate_df.to_csv(candidate_out, sep="\t", index=False)

    summary_path = outdir / "normalize_scomatic_results_summary.json"
    candidate_provenance_path = Path(str(candidate_out) + ".provenance.json")
    summary = {
        "mode": "read-only-normalization",
        "description": "Normalized external SComatic output files; SComatic was not run.",
        "n_input_files": int(len(scomatic_output_paths)),
        "n_input_rows": int(sum(source["rows"] for source in source_summaries)),
        "n_candidates": int(candidate_df.shape[0]),
        "n_cells": int(candidate_df["cell_id"].nunique()),
        "cell_ids": sorted(candidate_cells),
        "candidate_cell_semantics": candidate_cell_semantics,
        "source_files": source_summaries,
        "cell_annotation": cell_annotation,
        "provenance_metadata": provenance_metadata,
        "schema": SCOMATIC_CANDIDATE_COLUMNS,
        "terminology_note": TERMINOLOGY_NOTE,
        "output_scomatic_candidate_summary": str(candidate_out.resolve()),
        "output_scomatic_candidate_summary_provenance": str(candidate_provenance_path.resolve()),
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
    write_json(candidate_provenance_path, summary)
    return summary
