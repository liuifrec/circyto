from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping

from circyto.paths import resolve_manifest_path


CIRI_LONG_SCHEMA_VERSION = "circyto.ciri_long_input.v1"
CIRI_LONG_PLATFORM = "OXFORD_NANOPORE"
CIRI_LONG_MOLECULE_TYPE = "cDNA"
CIRI_LONG_LIBRARY_PREPARATION = "rolling_circle_reverse_transcription"
CIRI_LONG_INTERPRETATION_BOUNDARY = (
    "full_length_circrna_detection_from_rcrt_library"
)

CIRI_LONG_REQUIRED_COLUMNS = [
    "schema_version",
    "sample_id",
    "reads_path",
    "sequencing_platform",
    "molecule_type",
    "library_preparation",
    "circRNA_enrichment",
    "source_accession",
    "dataset_id",
    "reference_id",
    "reference_build",
    "biological_interpretation_boundary",
]

_SAFE_SAMPLE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
_POLYA_SELECTION = re.compile(
    r"(?:poly\s*\(\s*a\s*\)|poly[-_ ]?a(?:denylated)?|oligo\s*\(\s*d?t\s*\)|oligo[-_ ]?d?t)",
    flags=re.IGNORECASE,
)


@dataclass(frozen=True)
class CiriLongManifestRow:
    sample_id: str
    reads_path: str
    source_accession: str
    dataset_id: str
    reference_id: str
    reference_build: str
    sequencing_platform: str = CIRI_LONG_PLATFORM
    molecule_type: str = CIRI_LONG_MOLECULE_TYPE
    library_preparation: str = CIRI_LONG_LIBRARY_PREPARATION
    circRNA_enrichment: bool = True
    biological_interpretation_boundary: str = CIRI_LONG_INTERPRETATION_BOUNDARY
    schema_version: str = CIRI_LONG_SCHEMA_VERSION
    extra: Mapping[str, str] | None = None

    def to_dict(self) -> dict[str, str]:
        payload = {
            "schema_version": self.schema_version,
            "sample_id": self.sample_id,
            "reads_path": self.reads_path,
            "sequencing_platform": self.sequencing_platform,
            "molecule_type": self.molecule_type,
            "library_preparation": self.library_preparation,
            "circRNA_enrichment": "true" if self.circRNA_enrichment else "false",
            "source_accession": self.source_accession,
            "dataset_id": self.dataset_id,
            "reference_id": self.reference_id,
            "reference_build": self.reference_build,
            "biological_interpretation_boundary": self.biological_interpretation_boundary,
        }
        if self.extra:
            payload.update({str(key): str(value) for key, value in self.extra.items()})
        return payload


def _parse_boolean(value: object, *, field_name: str, prefix: str) -> bool:
    if isinstance(value, bool):
        return value
    normalized = str(value).strip().lower()
    if normalized in {"true", "1", "yes"}:
        return True
    if normalized in {"false", "0", "no"}:
        return False
    raise ValueError(
        f"{prefix}: {field_name} must be true or false; got {value!r}"
    )


def _has_supported_read_suffix(path: Path) -> bool:
    suffixes = [suffix.lower() for suffix in path.suffixes]
    return (
        suffixes[-1:] in ([".fa"], [".fasta"], [".fq"], [".fastq"])
        or suffixes[-2:]
        in (
            [".fa", ".gz"],
            [".fasta", ".gz"],
            [".fq", ".gz"],
            [".fastq", ".gz"],
        )
    )


def _poly_a_selection_evidence(row: CiriLongManifestRow) -> list[str]:
    evidence: list[str] = []
    for key, value in dict(row.extra or {}).items():
        normalized_key = key.strip().lower()
        if normalized_key not in {
            "library_selection",
            "rna_selection",
            "selection_method",
            "priming_method",
            "library_preparation_summary",
        }:
            continue
        if _POLYA_SELECTION.search(str(value)):
            evidence.append(f"{key}={value!r}")
    return evidence


def _validate_row(
    row: CiriLongManifestRow,
    *,
    manifest_path: Path,
    line_number: int,
    validate_files: bool,
) -> CiriLongManifestRow:
    prefix = f"{manifest_path}:{line_number}"
    reserved_extra = sorted(
        set(dict(row.extra or {})).intersection(CIRI_LONG_REQUIRED_COLUMNS)
    )
    if reserved_extra:
        raise ValueError(
            f"{prefix}: extra fields must not override required columns: "
            + ", ".join(reserved_extra)
        )
    if row.schema_version != CIRI_LONG_SCHEMA_VERSION:
        raise ValueError(
            f"{prefix}: unsupported schema_version={row.schema_version!r}; "
            f"expected {CIRI_LONG_SCHEMA_VERSION!r}"
        )
    if not row.sample_id:
        raise ValueError(f"{prefix}: empty sample_id")
    if not _SAFE_SAMPLE_ID.fullmatch(row.sample_id):
        raise ValueError(
            f"{prefix}: unsafe sample_id={row.sample_id!r}; use letters, digits, "
            "'.', '_' or '-' and do not include path separators"
        )
    if row.sequencing_platform != CIRI_LONG_PLATFORM:
        raise ValueError(
            f"{prefix}: sequencing_platform must be {CIRI_LONG_PLATFORM!r}; "
            f"got {row.sequencing_platform!r}"
        )
    if row.molecule_type != CIRI_LONG_MOLECULE_TYPE:
        if row.molecule_type.strip().lower().replace("-", "_") in {
            "direct_rna",
            "directrna",
        }:
            reason = "direct RNA is not a CIRI-long RCRT cDNA input"
        else:
            reason = "unknown or unsupported molecule chemistry"
        raise ValueError(
            f"{prefix}: molecule_type must be {CIRI_LONG_MOLECULE_TYPE!r}; "
            f"got {row.molecule_type!r} ({reason})"
        )
    if row.library_preparation != CIRI_LONG_LIBRARY_PREPARATION:
        raise ValueError(
            f"{prefix}: library_preparation must be "
            f"{CIRI_LONG_LIBRARY_PREPARATION!r}; got {row.library_preparation!r}. "
            "Ordinary ONT cDNA and unknown chemistry are not CIRI-long-compatible."
        )
    if not bool(row.circRNA_enrichment):
        raise ValueError(
            f"{prefix}: circRNA_enrichment must be true for a CIRI-long input"
        )
    poly_a_evidence = _poly_a_selection_evidence(row)
    if poly_a_evidence:
        raise ValueError(
            f"{prefix}: ordinary poly(A)/oligo(dT)-selected cDNA is rejected: "
            + ", ".join(poly_a_evidence)
        )
    if not row.source_accession.strip() and not row.dataset_id.strip():
        raise ValueError(
            f"{prefix}: at least one of source_accession or dataset_id must be populated"
        )
    if not row.reference_id.strip():
        raise ValueError(f"{prefix}: empty reference_id; reference identity is required")
    if not row.reference_build.strip():
        raise ValueError(
            f"{prefix}: empty reference_build; reference identity is required"
        )
    if (
        row.biological_interpretation_boundary
        != CIRI_LONG_INTERPRETATION_BOUNDARY
    ):
        raise ValueError(
            f"{prefix}: biological_interpretation_boundary must be "
            f"{CIRI_LONG_INTERPRETATION_BOUNDARY!r}"
        )
    if not row.reads_path.strip():
        raise ValueError(f"{prefix}: empty reads_path")
    reads_path = resolve_manifest_path(manifest_path, row.reads_path)
    if reads_path.suffix.lower() in {".bam", ".sam", ".cram"}:
        raise ValueError(
            f"{prefix}: reads_path must be raw FASTA/FASTQ, not a generic alignment: "
            f"{reads_path}"
        )
    if not _has_supported_read_suffix(reads_path):
        raise ValueError(
            f"{prefix}: reads_path must end in .fa, .fasta, .fq, .fastq, "
            "or a gzip-compressed form"
        )
    if validate_files and not reads_path.is_file():
        raise FileNotFoundError(f"{prefix}: reads_path not found: {reads_path}")
    return CiriLongManifestRow(
        sample_id=row.sample_id,
        reads_path=str(reads_path),
        source_accession=row.source_accession,
        dataset_id=row.dataset_id,
        reference_id=row.reference_id,
        reference_build=row.reference_build,
        sequencing_platform=row.sequencing_platform,
        molecule_type=row.molecule_type,
        library_preparation=row.library_preparation,
        circRNA_enrichment=True,
        biological_interpretation_boundary=row.biological_interpretation_boundary,
        schema_version=row.schema_version,
        extra=dict(row.extra or {}) or None,
    )


def read_ciri_long_manifest_tsv(
    path: Path,
    *,
    validate_files: bool = True,
) -> list[CiriLongManifestRow]:
    if not path.is_file():
        raise FileNotFoundError(f"CIRI-long manifest not found: {path}")
    rows: list[CiriLongManifestRow] = []
    seen: set[str] = set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = list(reader.fieldnames or [])
        missing = [
            column for column in CIRI_LONG_REQUIRED_COLUMNS if column not in header
        ]
        if missing:
            raise ValueError(
                f"CIRI-long manifest missing required columns at {path}: "
                + ", ".join(missing)
            )
        for line_number, raw in enumerate(reader, start=2):
            prefix = f"{path}:{line_number}"
            enrichment = _parse_boolean(
                raw.get("circRNA_enrichment", ""),
                field_name="circRNA_enrichment",
                prefix=prefix,
            )
            extra = {
                key: (value or "").strip()
                for key, value in raw.items()
                if key not in CIRI_LONG_REQUIRED_COLUMNS and (value or "").strip()
            }
            row = CiriLongManifestRow(
                schema_version=(raw.get("schema_version") or "").strip(),
                sample_id=(raw.get("sample_id") or "").strip(),
                reads_path=(raw.get("reads_path") or "").strip(),
                sequencing_platform=(raw.get("sequencing_platform") or "").strip(),
                molecule_type=(raw.get("molecule_type") or "").strip(),
                library_preparation=(raw.get("library_preparation") or "").strip(),
                circRNA_enrichment=enrichment,
                source_accession=(raw.get("source_accession") or "").strip(),
                dataset_id=(raw.get("dataset_id") or "").strip(),
                reference_id=(raw.get("reference_id") or "").strip(),
                reference_build=(raw.get("reference_build") or "").strip(),
                biological_interpretation_boundary=(
                    raw.get("biological_interpretation_boundary") or ""
                ).strip(),
                extra=extra or None,
            )
            checked = _validate_row(
                row,
                manifest_path=path,
                line_number=line_number,
                validate_files=validate_files,
            )
            if checked.sample_id in seen:
                raise ValueError(
                    f"{path}:{line_number}: duplicate sample_id={checked.sample_id!r}"
                )
            seen.add(checked.sample_id)
            rows.append(checked)
    if not rows:
        raise ValueError(f"CIRI-long manifest contains 0 data rows: {path}")
    return rows


def write_ciri_long_manifest_tsv(
    rows: Iterable[CiriLongManifestRow],
    path: Path,
) -> None:
    rows_list = list(rows)
    if not rows_list:
        raise ValueError("CIRI-long manifest must contain at least one row")
    validated: list[CiriLongManifestRow] = []
    seen: set[str] = set()
    extra_columns: set[str] = set()
    for line_number, row in enumerate(rows_list, start=2):
        checked = _validate_row(
            row,
            manifest_path=path,
            line_number=line_number,
            validate_files=False,
        )
        if checked.sample_id in seen:
            raise ValueError(
                f"{path}:{line_number}: duplicate sample_id={checked.sample_id!r}"
            )
        seen.add(checked.sample_id)
        extra_columns.update(dict(checked.extra or {}))
        validated.append(checked)

    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = CIRI_LONG_REQUIRED_COLUMNS + sorted(extra_columns)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in sorted(validated, key=lambda item: item.sample_id):
            payload = row.to_dict()
            reads_path = Path(payload["reads_path"])
            if reads_path.is_absolute():
                try:
                    payload["reads_path"] = str(
                        reads_path.relative_to(path.parent.resolve())
                    )
                except ValueError:
                    pass
            writer.writerow(payload)


def validate_ciri_long_manifest_tsv(
    path: Path,
    *,
    strict: bool = False,
) -> tuple[bool, list[str], dict[str, object]]:
    summary: dict[str, object] = {
        "schema_version": CIRI_LONG_SCHEMA_VERSION,
        "samples": 0,
        "reference_identities": [],
        "chemistry_gate": "not_evaluated",
        "missing_files": 0,
    }
    try:
        rows = read_ciri_long_manifest_tsv(path, validate_files=strict)
    except Exception as exc:
        message = str(exc)
        if "not found" in message:
            summary["missing_files"] = 1
        summary["chemistry_gate"] = "rejected"
        return False, [message], summary
    summary["samples"] = len(rows)
    summary["reference_identities"] = sorted(
        {f"{row.reference_id}:{row.reference_build}" for row in rows}
    )
    summary["chemistry_gate"] = "accepted_rcrt_circrna_enriched"
    return True, [], summary
