from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple
import csv


REQUIRED_COLUMNS = [
    "cell_id",
    "platform",
    "read1",
    "read2",
    "bam",
    "library_id",
    "n_input_reads",
]


@dataclass(frozen=True)
class ManifestRow:
    cell_id: str
    platform: str
    read1: str = ""
    read2: str = ""
    bam: str = ""
    library_id: str = ""
    n_input_reads: int = 0
    extra: Optional[Dict[str, str]] = None

    def to_dict(self) -> Dict[str, str]:
        d = {
            "cell_id": self.cell_id,
            "platform": self.platform,
            "read1": self.read1,
            "read2": self.read2,
            "bam": self.bam,
            "library_id": self.library_id,
            "n_input_reads": str(int(self.n_input_reads)),
        }
        if self.extra:
            d.update({k: str(v) for k, v in self.extra.items()})
        return d


def write_manifest_tsv(rows: Iterable[ManifestRow], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows_list = list(rows)

    # Stable column set: required first, then extras sorted.
    extras: List[str] = []
    seen = set(REQUIRED_COLUMNS)
    for r in rows_list:
        if r.extra:
            for k in r.extra.keys():
                if k not in seen:
                    extras.append(k)
                    seen.add(k)
    extras = sorted(extras)

    fieldnames = REQUIRED_COLUMNS + extras
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for r in sorted(rows_list, key=lambda x: x.cell_id):
            w.writerow(r.to_dict())


def validate_manifest_tsv(path: Path, strict: bool = False) -> Tuple[bool, List[str], Dict[str, int]]:
    errors: List[str] = []
    summary = {"cells": 0, "fastq_rows": 0, "bam_rows": 0, "missing_files": 0}

    if not path.exists():
        return False, [f"Manifest not found: {path}"], summary

    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        header = r.fieldnames or []
        for col in REQUIRED_COLUMNS:
            if col not in header:
                errors.append(f"Missing required column: {col}")
        if errors:
            return False, errors, summary

        seen_ids = set()
        for i, row in enumerate(r, start=2):
            summary["cells"] += 1

            cell_id = (row.get("cell_id") or "").strip()
            if not cell_id:
                errors.append(f"Line {i}: empty cell_id")
                continue
            if cell_id in seen_ids:
                errors.append(f"Line {i}: duplicate cell_id: {cell_id}")
            seen_ids.add(cell_id)

            platform = (row.get("platform") or "").strip()
            if not platform:
                errors.append(f"Line {i}: empty platform")

            library_id = (row.get("library_id") or "").strip()
            if not library_id:
                errors.append(f"Line {i}: empty library_id")

            n_input_reads = (row.get("n_input_reads") or "").strip()
            try:
                n = int(n_input_reads)
                if n < 0:
                    errors.append(f"Line {i}: n_input_reads < 0")
            except Exception:
                errors.append(f"Line {i}: n_input_reads is not an int: {n_input_reads}")

            read1 = (row.get("read1") or "").strip()
            read2 = (row.get("read2") or "").strip()
            bam = (row.get("bam") or "").strip()

            fastq_mode = bool(read1)
            bam_mode = bool(bam)

            if fastq_mode and bam_mode:
                errors.append(f"Line {i}: both FASTQ and BAM provided; choose one mode")
                continue
            if (not fastq_mode) and (not bam_mode):
                errors.append(f"Line {i}: neither FASTQ nor BAM provided")
                continue

            if fastq_mode:
                summary["fastq_rows"] += 1
                p1 = (path.parent / read1) if not Path(read1).is_absolute() else Path(read1)
                if strict and not p1.exists():
                    summary["missing_files"] += 1
                    errors.append(f"Line {i}: read1 not found: {p1}")
                if read2:
                    p2 = (path.parent / read2) if not Path(read2).is_absolute() else Path(read2)
                    if strict and not p2.exists():
                        summary["missing_files"] += 1
                        errors.append(f"Line {i}: read2 not found: {p2}")

            if bam_mode:
                summary["bam_rows"] += 1
                pb = (path.parent / bam) if not Path(bam).is_absolute() else Path(bam)
                if strict and not pb.exists():
                    summary["missing_files"] += 1
                    errors.append(f"Line {i}: bam not found: {pb}")

    return (len(errors) == 0), errors, summary
