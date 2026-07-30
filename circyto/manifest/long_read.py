from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from circyto.paths import resolve_manifest_path


LONG_READ_SCHEMA_VERSION = "circyto.long_read_single_cell.v1"
BIOLOGICAL_INTERPRETATION_BOUNDARY = (
    "experimental_long_read_interoperability_not_circrna_validation"
)

LONG_READ_REQUIRED_COLUMNS = [
    "schema_version",
    "cell_id",
    "long_read_fastq",
    "protocol",
    "sequencing_platform",
    "archive_library_selection",
    "library_preparation_summary",
    "molecule_type",
    "barcode_status",
    "source_accession",
    "dataset_id",
    "biological_interpretation_boundary",
]

VALID_MOLECULE_TYPES = {"cdna", "direct_rna"}
_SAFE_CELL_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


@dataclass(frozen=True)
class LongReadManifestRow:
    cell_id: str
    long_read_fastq: str
    protocol: str
    sequencing_platform: str
    archive_library_selection: str
    library_preparation_summary: str
    molecule_type: str
    barcode_status: str
    source_accession: str
    dataset_id: str
    biological_interpretation_boundary: str = BIOLOGICAL_INTERPRETATION_BOUNDARY
    schema_version: str = LONG_READ_SCHEMA_VERSION

    def to_dict(self) -> dict[str, str]:
        return {
            "schema_version": self.schema_version,
            "cell_id": self.cell_id,
            "long_read_fastq": self.long_read_fastq,
            "protocol": self.protocol,
            "sequencing_platform": self.sequencing_platform,
            "archive_library_selection": self.archive_library_selection,
            "library_preparation_summary": self.library_preparation_summary,
            "molecule_type": self.molecule_type,
            "barcode_status": self.barcode_status,
            "source_accession": self.source_accession,
            "dataset_id": self.dataset_id,
            "biological_interpretation_boundary": self.biological_interpretation_boundary,
        }


def _validate_row(
    row: LongReadManifestRow,
    *,
    path: Path,
    line_number: int,
    validate_files: bool,
) -> LongReadManifestRow:
    prefix = f"{path}:{line_number}"
    if row.schema_version != LONG_READ_SCHEMA_VERSION:
        raise ValueError(
            f"{prefix}: unsupported schema_version={row.schema_version!r}; "
            f"expected {LONG_READ_SCHEMA_VERSION!r}"
        )
    if not row.cell_id:
        raise ValueError(f"{prefix}: empty cell_id")
    if not _SAFE_CELL_ID.fullmatch(row.cell_id):
        raise ValueError(
            f"{prefix}: unsafe cell_id={row.cell_id!r}; use letters, digits, '.', '_' or '-' "
            "and do not include path separators"
        )
    for field_name in (
        "long_read_fastq",
        "protocol",
        "sequencing_platform",
        "archive_library_selection",
        "library_preparation_summary",
        "barcode_status",
        "source_accession",
        "dataset_id",
    ):
        if not getattr(row, field_name).strip():
            raise ValueError(f"{prefix}: empty required field {field_name}")
    if row.molecule_type not in VALID_MOLECULE_TYPES:
        raise ValueError(
            f"{prefix}: invalid molecule_type={row.molecule_type!r}; "
            f"expected one of {sorted(VALID_MOLECULE_TYPES)}"
        )
    if row.biological_interpretation_boundary != BIOLOGICAL_INTERPRETATION_BOUNDARY:
        raise ValueError(
            f"{prefix}: biological_interpretation_boundary must be "
            f"{BIOLOGICAL_INTERPRETATION_BOUNDARY!r}"
        )
    fastq = resolve_manifest_path(path, row.long_read_fastq)
    suffixes = [suffix.lower() for suffix in fastq.suffixes]
    if not (
        suffixes[-1:] in ([".fastq"], [".fq"])
        or suffixes[-2:] in ([".fastq", ".gz"], [".fq", ".gz"])
    ):
        raise ValueError(f"{prefix}: long_read_fastq must end in .fastq, .fq, .fastq.gz or .fq.gz")
    if validate_files and not fastq.is_file():
        raise FileNotFoundError(f"{prefix}: long_read_fastq not found: {fastq}")
    return LongReadManifestRow(
        **{
            **row.to_dict(),
            "long_read_fastq": str(fastq),
            "molecule_type": row.molecule_type.lower(),
        }
    )


def write_long_read_manifest_tsv(rows: Iterable[LongReadManifestRow], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows_list = list(rows)
    seen: set[str] = set()
    validated: list[LongReadManifestRow] = []
    for line_number, row in enumerate(rows_list, start=2):
        checked = _validate_row(
            row,
            path=path,
            line_number=line_number,
            validate_files=False,
        )
        if checked.cell_id in seen:
            raise ValueError(f"{path}:{line_number}: duplicate cell_id={checked.cell_id!r}")
        seen.add(checked.cell_id)
        validated.append(checked)
    if not validated:
        raise ValueError("Long-read manifest must contain at least one row")

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=LONG_READ_REQUIRED_COLUMNS, delimiter="\t")
        writer.writeheader()
        for row in sorted(validated, key=lambda item: item.cell_id):
            payload = row.to_dict()
            fastq = Path(payload["long_read_fastq"])
            if fastq.is_absolute():
                try:
                    payload["long_read_fastq"] = str(fastq.relative_to(path.parent.resolve()))
                except ValueError:
                    pass
            writer.writerow(payload)


def read_long_read_manifest_tsv(
    path: Path,
    *,
    validate_files: bool = True,
) -> list[LongReadManifestRow]:
    if not path.is_file():
        raise FileNotFoundError(f"Long-read manifest not found: {path}")

    rows: list[LongReadManifestRow] = []
    seen: set[str] = set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = list(reader.fieldnames or [])
        missing = [column for column in LONG_READ_REQUIRED_COLUMNS if column not in header]
        if missing:
            raise ValueError(
                f"Long-read manifest missing required columns at {path}: {', '.join(missing)}"
            )
        for line_number, raw in enumerate(reader, start=2):
            row = LongReadManifestRow(
                **{
                    column: (raw.get(column) or "").strip()
                    for column in LONG_READ_REQUIRED_COLUMNS
                }
            )
            row = _validate_row(
                row,
                path=path,
                line_number=line_number,
                validate_files=validate_files,
            )
            if row.cell_id in seen:
                raise ValueError(f"{path}:{line_number}: duplicate cell_id={row.cell_id!r}")
            seen.add(row.cell_id)
            rows.append(row)
    if not rows:
        raise ValueError(f"Long-read manifest contains 0 data rows: {path}")
    return rows


def validate_long_read_manifest_tsv(
    path: Path,
    *,
    strict: bool = False,
) -> tuple[bool, list[str], dict[str, int]]:
    summary = {"cells": 0, "cdna_rows": 0, "direct_rna_rows": 0, "missing_files": 0}
    try:
        rows = read_long_read_manifest_tsv(path, validate_files=strict)
    except Exception as exc:
        message = str(exc)
        if "not found" in message:
            summary["missing_files"] += 1
        return False, [message], summary
    summary["cells"] = len(rows)
    summary["cdna_rows"] = sum(row.molecule_type == "cdna" for row in rows)
    summary["direct_rna_rows"] = sum(row.molecule_type == "direct_rna" for row in rows)
    return True, [], summary
