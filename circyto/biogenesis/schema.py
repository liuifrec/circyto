from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np
import pandas as pd


BIOGENESIS_SCHEMA_VERSION = "1.0"

PROVENANCE_FIELDS = (
    "dataset_id",
    "detector",
    "reference_genome",
    "reference_annotation",
    "workflow_uuid",
)


@dataclass(frozen=True)
class TableSchema:
    name: str
    version: str
    required_fields: tuple[str, ...]
    optional_fields: tuple[str, ...]


@dataclass(frozen=True)
class ValidatedBiogenesisBundle:
    candidates: pd.DataFrame
    cell_contexts: pd.DataFrame
    observations: pd.DataFrame


class BiogenesisSchemaError(ValueError):
    """Raised when a biogenesis table violates its versioned data contract."""


CANDIDATE_SCHEMA = TableSchema(
    name="circRNA_candidates",
    version=BIOGENESIS_SCHEMA_VERSION,
    required_fields=(
        "schema_version",
        "circ_id",
        "dataset_id",
        "chrom",
        "start",
        "end",
        "strand",
        "coordinate_system",
        "detector",
        "reference_genome",
        "reference_annotation",
        "workflow_uuid",
        "label_status",
    ),
    optional_fields=(
        "host_gene",
        "host_gene_id",
        "donor_window",
        "acceptor_window",
        "exon_count",
        "circ_length",
        "upstream_exon_length",
        "downstream_exon_length",
        "upstream_intron_length",
        "downstream_intron_length",
        "upstream_repeat_family",
        "downstream_repeat_family",
        "repeat_pairing_score",
        "donor_splice_site",
        "acceptor_splice_site",
        "donor_splice_score",
        "acceptor_splice_score",
        "known_status",
        "known_database",
        "known_database_id",
    ),
)

CELL_CONTEXT_SCHEMA = TableSchema(
    name="cell_contexts",
    version=BIOGENESIS_SCHEMA_VERSION,
    required_fields=(
        "schema_version",
        "cell_id",
        "dataset_id",
        "donor_id",
        "protocol",
        "detector",
        "reference_genome",
        "reference_annotation",
        "workflow_uuid",
        "total_rna_counts",
        "detected_genes",
        "circRNA_count",
        "circRNA_total_support",
    ),
    optional_fields=(
        "condition",
        "cell_type",
        "batch_id",
        "library_id",
        "read_layout",
        "strandedness",
        "assigned_reads",
        "mitochondrial_fraction",
        "ribosomal_fraction",
        "alignment_status",
        "detector_status",
    ),
)

OBSERVATION_SCHEMA = TableSchema(
    name="cell_circ_observations",
    version=BIOGENESIS_SCHEMA_VERSION,
    required_fields=(
        "schema_version",
        "cell_id",
        "circ_id",
        "dataset_id",
        "donor_id",
        "protocol",
        "detector",
        "reference_genome",
        "reference_annotation",
        "workflow_uuid",
        "count",
        "detected",
        "bsj_support",
        "candidate_label_status",
        "detectability_offset",
    ),
    optional_fields=(
        "host_gene_expression",
        "label_confidence",
        "count_uncertainty",
        "detection_probability",
        "observation_weight",
    ),
)

_LABEL_STATUSES = {"positive", "unlabelled"}
_KNOWN_STATUSES = {"known", "novel", "unknown"}
_COORDINATE_SYSTEMS = {"0-based-half-open", "1-based-closed"}
_CELL_NUMERIC_PREFIXES = (
    "host_expr__",
    "rbp_program__",
    "splicing_program__",
    "latent__",
    "multimodal__",
)


def _frame_copy(frame: pd.DataFrame, schema: TableSchema) -> pd.DataFrame:
    if not isinstance(frame, pd.DataFrame):
        raise BiogenesisSchemaError(
            f"{schema.name}: expected a pandas DataFrame, got {type(frame).__name__}"
        )
    if frame.empty:
        raise BiogenesisSchemaError(f"{schema.name}: table must contain at least one record")
    missing = [field for field in schema.required_fields if field not in frame.columns]
    if missing:
        raise BiogenesisSchemaError(
            f"{schema.name}: missing required columns: {', '.join(missing)}"
        )
    return frame.copy()


def _example_rows(mask: pd.Series) -> list[int]:
    return [int(position) for position in np.flatnonzero(mask.to_numpy())[:5]]


def _validate_version(frame: pd.DataFrame, schema: TableSchema) -> None:
    values = frame["schema_version"].astype("string").fillna("").str.strip()
    bad = values != schema.version
    if bad.any():
        found = sorted(set(values.loc[bad].tolist()))
        raise BiogenesisSchemaError(
            f"{schema.name}.schema_version: expected '{schema.version}', "
            f"found {found}; row positions {_example_rows(bad)}"
        )
    frame["schema_version"] = schema.version


def _normalize_required_strings(
    frame: pd.DataFrame,
    columns: Iterable[str],
    *,
    table_name: str,
) -> None:
    for column in columns:
        values = frame[column].astype("string").fillna("").str.strip()
        bad = values == ""
        if bad.any():
            raise BiogenesisSchemaError(
                f"{table_name}.{column}: values must be non-empty; "
                f"row positions {_example_rows(bad)}"
            )
        frame[column] = values.astype(str)


def _normalize_optional_strings(frame: pd.DataFrame, columns: Iterable[str]) -> None:
    for column in columns:
        if column in frame.columns:
            frame[column] = frame[column].astype("string").fillna("").str.strip().astype(str)


def _normalize_numeric(
    frame: pd.DataFrame,
    column: str,
    *,
    table_name: str,
    required: bool,
    integer: bool = False,
    minimum: float | None = None,
    maximum: float | None = None,
) -> None:
    if column not in frame.columns:
        return
    raw = frame[column]
    numeric = pd.to_numeric(raw, errors="coerce")
    invalid = numeric.isna() & raw.astype("string").fillna("").str.strip().ne("")
    if required:
        invalid |= numeric.isna()
    invalid |= numeric.notna() & ~np.isfinite(numeric.astype(float))
    if integer:
        invalid |= numeric.notna() & (numeric.astype(float) % 1 != 0)
    if minimum is not None:
        invalid |= numeric.notna() & (numeric < minimum)
    if maximum is not None:
        invalid |= numeric.notna() & (numeric > maximum)
    if invalid.any():
        constraints = ["numeric"]
        if integer:
            constraints.append("integral")
        if minimum is not None:
            constraints.append(f">= {minimum:g}")
        if maximum is not None:
            constraints.append(f"<= {maximum:g}")
        raise BiogenesisSchemaError(
            f"{table_name}.{column}: expected {' and '.join(constraints)} values; "
            f"row positions {_example_rows(invalid)}"
        )
    if integer:
        frame[column] = numeric.astype("int64" if required else "Int64")
    else:
        frame[column] = numeric.astype(float)


def _validate_enum(
    frame: pd.DataFrame,
    column: str,
    allowed: set[str],
    *,
    table_name: str,
    allow_empty: bool = False,
) -> None:
    if column not in frame.columns:
        return
    values = frame[column].astype("string").fillna("").str.strip().str.lower()
    bad = ~values.isin(allowed | ({""} if allow_empty else set()))
    if bad.any():
        found = sorted(set(values.loc[bad].tolist()))
        raise BiogenesisSchemaError(
            f"{table_name}.{column}: expected one of {sorted(allowed)}, "
            f"found {found}; row positions {_example_rows(bad)}"
        )
    frame[column] = values.astype(str)


def _validate_unique(
    frame: pd.DataFrame,
    columns: list[str],
    *,
    table_name: str,
) -> None:
    duplicate = frame.duplicated(columns, keep=False)
    if duplicate.any():
        examples = frame.loc[duplicate, columns].head(5).to_dict("records")
        raise BiogenesisSchemaError(
            f"{table_name}: duplicate key {tuple(columns)}; examples {examples}"
        )


def _normalize_boolean(frame: pd.DataFrame, column: str, *, table_name: str) -> None:
    values = frame[column].astype("string").fillna("").str.strip().str.lower()
    mapping = {
        "true": True,
        "1": True,
        "yes": True,
        "false": False,
        "0": False,
        "no": False,
    }
    bad = ~values.isin(mapping)
    if bad.any():
        found = sorted(set(values.loc[bad].tolist()))
        raise BiogenesisSchemaError(
            f"{table_name}.{column}: expected boolean true/false values, "
            f"found {found}; row positions {_example_rows(bad)}"
        )
    frame[column] = values.map(mapping).astype(bool)


def validate_candidate_records(frame: pd.DataFrame) -> pd.DataFrame:
    """Validate and normalize circRNA candidate records for schema v1.0."""
    schema = CANDIDATE_SCHEMA
    out = _frame_copy(frame, schema)
    _validate_version(out, schema)
    _normalize_required_strings(
        out,
        (
            "circ_id",
            "dataset_id",
            "chrom",
            "strand",
            "coordinate_system",
            "detector",
            "reference_genome",
            "reference_annotation",
            "workflow_uuid",
            "label_status",
        ),
        table_name=schema.name,
    )
    _normalize_optional_strings(
        out,
        (
            "host_gene",
            "host_gene_id",
            "donor_window",
            "acceptor_window",
            "upstream_repeat_family",
            "downstream_repeat_family",
            "donor_splice_site",
            "acceptor_splice_site",
            "known_status",
            "known_database",
            "known_database_id",
        ),
    )
    _normalize_numeric(
        out, "start", table_name=schema.name, required=True, integer=True, minimum=0
    )
    _normalize_numeric(
        out, "end", table_name=schema.name, required=True, integer=True, minimum=1
    )
    for column in (
        "exon_count",
        "circ_length",
        "upstream_exon_length",
        "downstream_exon_length",
        "upstream_intron_length",
        "downstream_intron_length",
    ):
        _normalize_numeric(
            out, column, table_name=schema.name, required=False, integer=True, minimum=0
        )
    for column in (
        "repeat_pairing_score",
        "donor_splice_score",
        "acceptor_splice_score",
    ):
        _normalize_numeric(
            out, column, table_name=schema.name, required=False, minimum=0
        )
    bad_interval = out["end"] <= out["start"]
    if bad_interval.any():
        raise BiogenesisSchemaError(
            f"{schema.name}: end must be greater than start; "
            f"row positions {_example_rows(bad_interval)}"
        )
    _validate_enum(out, "strand", {"+", "-"}, table_name=schema.name)
    _validate_enum(
        out, "coordinate_system", _COORDINATE_SYSTEMS, table_name=schema.name
    )
    invalid_one_based_start = (
        (out["coordinate_system"] == "1-based-closed") & (out["start"] < 1)
    )
    if invalid_one_based_start.any():
        raise BiogenesisSchemaError(
            f"{schema.name}: 1-based-closed coordinates require start >= 1; "
            f"row positions {_example_rows(invalid_one_based_start)}"
        )
    _validate_enum(out, "label_status", _LABEL_STATUSES, table_name=schema.name)
    _validate_enum(
        out,
        "known_status",
        _KNOWN_STATUSES,
        table_name=schema.name,
        allow_empty=True,
    )
    _validate_unique(out, ["circ_id"], table_name=schema.name)
    return out


def validate_cell_context_records(frame: pd.DataFrame) -> pd.DataFrame:
    """Validate and normalize per-cell modelling context records for schema v1.0."""
    schema = CELL_CONTEXT_SCHEMA
    out = _frame_copy(frame, schema)
    _validate_version(out, schema)
    _normalize_required_strings(
        out,
        (
            "cell_id",
            "dataset_id",
            "donor_id",
            "protocol",
            "detector",
            "reference_genome",
            "reference_annotation",
            "workflow_uuid",
        ),
        table_name=schema.name,
    )
    _normalize_optional_strings(
        out,
        (
            "condition",
            "cell_type",
            "batch_id",
            "library_id",
            "read_layout",
            "strandedness",
            "alignment_status",
            "detector_status",
        ),
    )
    for column in (
        "total_rna_counts",
        "detected_genes",
        "circRNA_count",
        "circRNA_total_support",
    ):
        _normalize_numeric(
            out, column, table_name=schema.name, required=True, integer=True, minimum=0
        )
    _normalize_numeric(
        out,
        "assigned_reads",
        table_name=schema.name,
        required=False,
        integer=True,
        minimum=0,
    )
    for column in ("mitochondrial_fraction", "ribosomal_fraction"):
        _normalize_numeric(
            out,
            column,
            table_name=schema.name,
            required=False,
            minimum=0,
            maximum=1,
        )
    for column in out.columns:
        if column.startswith(_CELL_NUMERIC_PREFIXES):
            _normalize_numeric(
                out, column, table_name=schema.name, required=False
            )
    _validate_unique(out, ["cell_id"], table_name=schema.name)
    return out


def validate_observation_records(frame: pd.DataFrame) -> pd.DataFrame:
    """Validate and normalize cell-by-circRNA observations for schema v1.0."""
    schema = OBSERVATION_SCHEMA
    out = _frame_copy(frame, schema)
    _validate_version(out, schema)
    _normalize_required_strings(
        out,
        (
            "cell_id",
            "circ_id",
            "dataset_id",
            "donor_id",
            "protocol",
            "detector",
            "reference_genome",
            "reference_annotation",
            "workflow_uuid",
            "candidate_label_status",
        ),
        table_name=schema.name,
    )
    for column in ("count", "bsj_support"):
        _normalize_numeric(
            out, column, table_name=schema.name, required=True, minimum=0
        )
    _normalize_numeric(
        out, "detectability_offset", table_name=schema.name, required=True
    )
    _normalize_numeric(
        out,
        "host_gene_expression",
        table_name=schema.name,
        required=False,
        minimum=0,
    )
    _normalize_numeric(
        out,
        "label_confidence",
        table_name=schema.name,
        required=False,
        minimum=0,
        maximum=1,
    )
    _normalize_numeric(
        out,
        "count_uncertainty",
        table_name=schema.name,
        required=False,
        minimum=0,
    )
    _normalize_numeric(
        out,
        "detection_probability",
        table_name=schema.name,
        required=False,
        minimum=0,
        maximum=1,
    )
    _normalize_numeric(
        out,
        "observation_weight",
        table_name=schema.name,
        required=False,
        minimum=0,
    )
    _normalize_boolean(out, "detected", table_name=schema.name)
    _validate_enum(
        out,
        "candidate_label_status",
        _LABEL_STATUSES,
        table_name=schema.name,
    )
    has_signal = (out["count"] > 0) | (out["bsj_support"] > 0)
    inconsistent = out["detected"] != has_signal
    if inconsistent.any():
        raise BiogenesisSchemaError(
            f"{schema.name}.detected: must equal whether count or bsj_support is positive; "
            f"row positions {_example_rows(inconsistent)}"
        )
    _validate_unique(out, ["cell_id", "circ_id"], table_name=schema.name)
    return out


def _validate_linked_values(
    observations: pd.DataFrame,
    reference: pd.DataFrame,
    *,
    key: str,
    columns: Iterable[str],
    reference_name: str,
) -> None:
    indexed = reference.set_index(key)
    for column in columns:
        expected = observations[key].map(indexed[column])
        mismatch = observations[column].astype(str) != expected.astype(str)
        if mismatch.any():
            examples = observations.loc[mismatch, [key, column]].head(5).to_dict("records")
            raise BiogenesisSchemaError(
                f"biogenesis_bundle: observations.{column} does not match "
                f"{reference_name}.{column} for linked {key}; examples {examples}"
            )


def validate_biogenesis_bundle(
    candidates: pd.DataFrame,
    cell_contexts: pd.DataFrame,
    observations: pd.DataFrame,
) -> ValidatedBiogenesisBundle:
    """Validate all three tables and their foreign-key/provenance relationships."""
    validated_candidates = validate_candidate_records(candidates)
    validated_contexts = validate_cell_context_records(cell_contexts)
    validated_observations = validate_observation_records(observations)

    unknown_circs = sorted(
        set(validated_observations["circ_id"]) - set(validated_candidates["circ_id"])
    )
    if unknown_circs:
        raise BiogenesisSchemaError(
            f"biogenesis_bundle: observations reference unknown circ_id values: "
            f"{unknown_circs[:5]}"
        )
    unknown_cells = sorted(
        set(validated_observations["cell_id"]) - set(validated_contexts["cell_id"])
    )
    if unknown_cells:
        raise BiogenesisSchemaError(
            f"biogenesis_bundle: observations reference unknown cell_id values: "
            f"{unknown_cells[:5]}"
        )

    _validate_linked_values(
        validated_observations,
        validated_candidates,
        key="circ_id",
        columns=PROVENANCE_FIELDS,
        reference_name=CANDIDATE_SCHEMA.name,
    )
    _validate_linked_values(
        validated_observations,
        validated_contexts,
        key="cell_id",
        columns=(
            "dataset_id",
            "donor_id",
            "protocol",
            "detector",
            "reference_genome",
            "reference_annotation",
            "workflow_uuid",
        ),
        reference_name=CELL_CONTEXT_SCHEMA.name,
    )
    candidate_labels = validated_candidates.set_index("circ_id")["label_status"]
    expected_labels = validated_observations["circ_id"].map(candidate_labels)
    mismatch = (
        validated_observations["candidate_label_status"].astype(str)
        != expected_labels.astype(str)
    )
    if mismatch.any():
        examples = validated_observations.loc[
            mismatch, ["circ_id", "candidate_label_status"]
        ].head(5).to_dict("records")
        raise BiogenesisSchemaError(
            "biogenesis_bundle: observation candidate_label_status does not match "
            f"candidate label_status; examples {examples}"
        )

    return ValidatedBiogenesisBundle(
        candidates=validated_candidates,
        cell_contexts=validated_contexts,
        observations=validated_observations,
    )
