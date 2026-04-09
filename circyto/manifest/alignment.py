from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from circyto.paths import resolve_manifest_path


ALIGNMENT_REQUIRED_COLUMNS = [
    "cell_id",
    "bam",
    "sam",
    "group_id",
    "read_layout",
    "aligner",
    "reference",
    "cache_key",
    "source_manifest",
]


@dataclass(frozen=True)
class AlignmentManifestRow:
    cell_id: str
    bam: str = ""
    sam: str = ""
    group_id: str = ""
    read_layout: str = ""
    aligner: str = ""
    reference: str = ""
    cache_key: str = ""
    source_manifest: str = ""
    extra: Optional[Dict[str, str]] = None

    def to_dict(self) -> Dict[str, str]:
        data = {
            "cell_id": self.cell_id,
            "bam": self.bam,
            "sam": self.sam,
            "group_id": self.group_id,
            "read_layout": self.read_layout,
            "aligner": self.aligner,
            "reference": self.reference,
            "cache_key": self.cache_key,
            "source_manifest": self.source_manifest,
        }
        if self.extra:
            data.update({k: str(v) for k, v in self.extra.items()})
        return data

    @property
    def alignment_path(self) -> str:
        return self.bam or self.sam


def write_alignment_manifest_tsv(rows: Iterable[AlignmentManifestRow], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows_list = list(rows)

    extras: List[str] = []
    seen = set(ALIGNMENT_REQUIRED_COLUMNS)
    for row in rows_list:
        if row.extra:
            for key in row.extra:
                if key not in seen:
                    extras.append(key)
                    seen.add(key)
    fieldnames = ALIGNMENT_REQUIRED_COLUMNS + sorted(extras)

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in sorted(rows_list, key=lambda item: item.cell_id):
            writer.writerow(row.to_dict())


def read_alignment_manifest_tsv(path: Path, *, validate_files: bool = True) -> List[AlignmentManifestRow]:
    if not path.exists():
        raise FileNotFoundError(f"Alignment manifest not found: {path}")

    rows: List[AlignmentManifestRow] = []
    seen_ids: set[str] = set()

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = reader.fieldnames or []
        if "cell_id" not in header:
            raise KeyError(f"Alignment manifest missing required column 'cell_id': {path}")

        for i, raw in enumerate(reader, start=2):
            cell_id = (raw.get("cell_id") or "").strip()
            if not cell_id:
                raise ValueError(f"Empty cell_id at {path}:{i}")
            if cell_id in seen_ids:
                raise ValueError(f"Duplicate cell_id '{cell_id}' at {path}:{i}")
            seen_ids.add(cell_id)

            bam = (raw.get("bam") or "").strip()
            sam = (raw.get("sam") or "").strip()
            if not bam and not sam:
                raise ValueError(f"Alignment manifest row missing bam/sam for cell_id={cell_id} at {path}:{i}")
            if bam and sam:
                raise ValueError(f"Alignment manifest row has both bam and sam for cell_id={cell_id} at {path}:{i}")

            alignment_path = resolve_manifest_path(path, bam or sam)
            if validate_files and not alignment_path.exists():
                raise FileNotFoundError(
                    f"Alignment file not found for cell_id={cell_id}: {alignment_path}"
                )

            extras = {
                key: value
                for key, value in raw.items()
                if key not in ALIGNMENT_REQUIRED_COLUMNS and value not in (None, "")
            }

            rows.append(
                AlignmentManifestRow(
                    cell_id=cell_id,
                    bam=str(alignment_path) if bam else "",
                    sam=str(alignment_path) if sam else "",
                    group_id=(raw.get("group_id") or "").strip(),
                    read_layout=(raw.get("read_layout") or "").strip(),
                    aligner=(raw.get("aligner") or "").strip(),
                    reference=(raw.get("reference") or "").strip(),
                    cache_key=(raw.get("cache_key") or "").strip(),
                    source_manifest=(raw.get("source_manifest") or "").strip(),
                    extra=extras or None,
                )
            )

    if not rows:
        raise ValueError(f"Alignment manifest contains 0 data rows: {path}")
    return rows


def validate_alignment_manifest_tsv(path: Path, strict: bool = False) -> Tuple[bool, List[str], Dict[str, int]]:
    errors: List[str] = []
    summary = {"cells": 0, "bam_rows": 0, "sam_rows": 0, "missing_files": 0, "invalid_layouts": 0}

    try:
        rows = read_alignment_manifest_tsv(path, validate_files=strict)
    except Exception as exc:
        return False, [str(exc)], summary

    for row in rows:
        summary["cells"] += 1
        if row.bam:
            summary["bam_rows"] += 1
        if row.sam:
            summary["sam_rows"] += 1
        if row.read_layout and row.read_layout not in {"single-end", "paired-end"}:
            summary["invalid_layouts"] += 1
            errors.append(f"Invalid read_layout for cell_id={row.cell_id}: {row.read_layout}")
        alignment_path = Path(row.alignment_path)
        if strict and not alignment_path.exists():
            summary["missing_files"] += 1
            errors.append(f"Alignment file not found: {alignment_path}")
        if strict and row.reference:
            ref_path = resolve_manifest_path(path, row.reference)
            if not ref_path.exists():
                summary["missing_files"] += 1
                errors.append(f"Reference file not found for cell_id={row.cell_id}: {ref_path}")

    return len(errors) == 0, errors, summary
