from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Optional

from circyto.detectors.ciri3 import Ciri3Detector
from circyto.manifest.alignment import (
    AlignmentManifestRow,
    read_alignment_manifest_tsv,
    write_alignment_manifest_tsv,
)
from circyto.pipeline.align_manifest import (
    plan_alignment_cache,
    prepare_alignment_cache,
    read_source_manifest,
    run_detector_alignment_manifest,
)
from circyto.pipeline.collect import collect_matrix
from circyto.protocols import get_protocol_preset, normalize_protocol_name


ProgressFn = Callable[[str], None]
VALID_STRANDEDNESS = {"", "forward", "reverse", "unstranded"}


@dataclass(frozen=True)
class RunCiri3Params:
    manifest: Path
    outdir: Path
    genome_fasta: Path
    gtf: Path
    star_index: Path
    protocol: Optional[str] = None
    strandedness: Optional[str] = None
    threads: int = 8
    parallel: int = 1
    chunk_size: int = 1
    dry_run: bool = False
    fail_fast: bool = False
    command_template: Optional[str] = None


def _emit(progress: ProgressFn, message: str) -> None:
    progress(f"[run-ciri3] {message}")


def _read_manifest_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = list(reader.fieldnames or [])
        rows = [{k: "" if v is None else str(v) for k, v in row.items()} for row in reader]
    return header, rows


def _write_manifest_rows(path: Path, header: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in header})


def _normalize_strandedness(value: str | None) -> str:
    text = str(value or "").strip().lower()
    if text not in VALID_STRANDEDNESS:
        supported = ", ".join(sorted(v for v in VALID_STRANDEDNESS if v))
        raise ValueError(f"Unsupported strandedness '{value}'. Supported: {supported}")
    return text


def _resolve_row_protocol(row, default_protocol: str | None) -> str:
    row_protocol = normalize_protocol_name(getattr(row, "protocol", "") or getattr(row, "platform", ""))
    cli_protocol = normalize_protocol_name(default_protocol)
    if row_protocol and cli_protocol and row_protocol != cli_protocol:
        raise ValueError(
            f"Protocol mismatch for sample_id/cell_id={row.cell_id}: manifest={row_protocol} cli={cli_protocol}"
        )
    resolved = row_protocol or cli_protocol
    if not resolved:
        raise ValueError(
            f"Missing protocol for sample_id/cell_id={row.cell_id}. Provide a manifest protocol column or --protocol."
        )
    return get_protocol_preset(resolved).name


def _resolve_row_strandedness(row, default_strandedness: str | None) -> str:
    row_value = _normalize_strandedness(getattr(row, "strandedness", ""))
    cli_value = _normalize_strandedness(default_strandedness)
    if row_value and cli_value and row_value != cli_value:
        raise ValueError(
            f"Strandedness mismatch for sample_id/cell_id={row.cell_id}: manifest={row_value} cli={cli_value}"
        )
    return row_value or cli_value


def _build_subset_rows(
    rows: list[dict[str, str]],
    *,
    keep_cell_ids: set[str],
    protocol_by_cell: dict[str, str],
    strandedness_by_cell: dict[str, str],
) -> list[dict[str, str]]:
    subset: list[dict[str, str]] = []
    for row in rows:
        cell_id = (row.get("cell_id") or row.get("sample_id") or "").strip()
        if cell_id not in keep_cell_ids:
            continue
        payload = dict(row)
        payload["cell_id"] = cell_id
        payload["protocol"] = protocol_by_cell[cell_id]
        payload["strandedness"] = strandedness_by_cell[cell_id]
        subset.append(payload)
    return subset


def _subset_header(header: list[str]) -> list[str]:
    ordered = list(header)
    for column in ("cell_id", "protocol", "strandedness", "read_layout"):
        if column not in ordered:
            ordered.append(column)
    return ordered


def _star_extra_flags(rows) -> str:
    any_gz = any(
        path is not None and path.suffix == ".gz"
        for row in rows
        for path in (row.read1, row.effective_read2)
    )
    parts = ["--genomeDir", "{star_index}"]
    if any_gz:
        parts.extend(["--readFilesCommand", "zcat"])
    return " ".join(parts)


def _dry_run_alignment_rows(
    rows,
    *,
    outdir: Path,
    aligner: str,
    ref_fa: Path,
) -> list[AlignmentManifestRow]:
    planned: list[AlignmentManifestRow] = []
    for row in rows:
        staged_dir = outdir / "planned"
        sam_path = staged_dir / f"{row.cell_id}.sam"
        kwargs = {
            "cell_id": row.cell_id,
            "bam": "",
            "sam": str(sam_path.resolve()),
            "group_id": row.group_id,
            "read_layout": row.read_layout,
            "aligner": aligner,
            "reference": str(ref_fa.resolve()),
            "cache_key": f"dryrun-{row.cell_id}-{aligner}",
            "source_manifest": "",
            "mapper_mode": "1" if aligner == "star" else "0",
            "artifact_bucket": "star" if aligner == "star" else "bwa_mem",
            "sortedness": "unsorted",
        }
        if aligner == "star":
            kwargs["chimeric_junction"] = str((staged_dir / f"{row.cell_id}.Chimeric.out.junction").resolve())
            kwargs["unmapped_mate1"] = str((staged_dir / f"{row.cell_id}.Unmapped.out.mate1").resolve())
            kwargs["unmapped_mate2"] = str((staged_dir / f"{row.cell_id}.Unmapped.out.mate2").resolve())
            kwargs["bwa_sam"] = str((staged_dir / f"{row.cell_id}.bwa_rescue.sam").resolve())
        planned.append(AlignmentManifestRow(**kwargs))
    return planned


def _merge_alignment_manifests(paths: list[Path], output_path: Path) -> Path:
    rows: list[AlignmentManifestRow] = []
    for path in paths:
        rows.extend(read_alignment_manifest_tsv(path, validate_files=False))
    write_alignment_manifest_tsv(rows, output_path)
    return output_path


def _write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _print_dry_run_commands(progress: ProgressFn, plan_path: Path, detector_plan_path: Path) -> None:
    alignment_plan = json.loads(plan_path.read_text(encoding="utf-8"))
    detector_plan = json.loads(detector_plan_path.read_text(encoding="utf-8"))
    for item in alignment_plan.get("command_preview", []):
        _emit(progress, f"STAR/BWA plan {item['cell_id']}: {item['command']}")
    for item in detector_plan.get("command_preview", []):
        _emit(progress, f"CIRI3 plan {item['cell_id']}: {item['command']}")


def _write_alignment_summary(path: Path, *, source_rows, manifest_paths: list[Path], dry_run: bool) -> None:
    payload = {
        "workflow": "run-ciri3",
        "dry_run": dry_run,
        "planned_cells": len(source_rows),
        "read_layout_counts": {
            "paired-end": sum(1 for row in source_rows if row.read_layout == "paired-end"),
            "single-end": sum(1 for row in source_rows if row.read_layout == "single-end"),
        },
        "sub_manifests": [str(item.resolve()) for item in manifest_paths],
    }
    _write_json(path, payload)


def run_ciri3_workflow(params: RunCiri3Params, *, progress: ProgressFn = print) -> None:
    root_dirs = {
        "root": params.outdir,
        "align": params.outdir / "align",
        "ciri3": params.outdir / "ciri3",
        "matrix": params.outdir / "matrix",
        "qc": params.outdir / "qc",
        "logs": params.outdir / "logs",
        "manifests": params.outdir / "manifests",
    }
    for path in root_dirs.values():
        path.mkdir(parents=True, exist_ok=True)

    source_rows = read_source_manifest(params.manifest, validate_files=True)
    header, raw_rows = _read_manifest_rows(params.manifest)
    if len(source_rows) != len(raw_rows):
        raise ValueError("Manifest parsing mismatch between structured and raw readers")

    protocol_by_cell: dict[str, str] = {}
    strandedness_by_cell: dict[str, str] = {}
    paired_cells: set[str] = set()
    single_cells: set[str] = set()
    for row in source_rows:
        protocol_by_cell[row.cell_id] = _resolve_row_protocol(row, params.protocol)
        strandedness_by_cell[row.cell_id] = _resolve_row_strandedness(row, params.strandedness)
        if row.read_layout == "paired-end":
            paired_cells.add(row.cell_id)
        else:
            single_cells.add(row.cell_id)

    subset_header = _subset_header(header)
    subset_manifests: list[tuple[str, Path, list]] = []
    if paired_cells:
        paired_rows = [row for row in source_rows if row.cell_id in paired_cells]
        paired_manifest = root_dirs["manifests"] / "paired_end_manifest.tsv"
        _write_manifest_rows(
            paired_manifest,
            subset_header,
            _build_subset_rows(
                raw_rows,
                keep_cell_ids=paired_cells,
                protocol_by_cell=protocol_by_cell,
                strandedness_by_cell=strandedness_by_cell,
            ),
        )
        subset_manifests.append(("star", paired_manifest, paired_rows))
    if single_cells:
        single_rows = [row for row in source_rows if row.cell_id in single_cells]
        single_manifest = root_dirs["manifests"] / "single_end_manifest.tsv"
        _write_manifest_rows(
            single_manifest,
            subset_header,
            _build_subset_rows(
                raw_rows,
                keep_cell_ids=single_cells,
                protocol_by_cell=protocol_by_cell,
                strandedness_by_cell=strandedness_by_cell,
            ),
        )
        subset_manifests.append(("bwa-mem", single_manifest, single_rows))

    det = Ciri3Detector(command_template=params.command_template)
    manifest_paths: list[Path] = []
    plan_paths: list[Path] = []
    detector_manifest = root_dirs["align"] / "alignment_manifest.tsv"

    if params.dry_run:
        synthetic_rows: list[AlignmentManifestRow] = []
        for aligner, subset_manifest, rows in subset_manifests:
            subdir = root_dirs["align"] / aligner.replace("-", "_")
            extra_flags = _star_extra_flags(rows).format(star_index=str(params.star_index.resolve())) if aligner == "star" else ""
            plan = plan_alignment_cache(
                manifest=subset_manifest,
                outdir=subdir,
                aligner=aligner,
                ref_fa=params.genome_fasta,
                detector_hint="ciri3",
                threads=params.threads,
                parallel=params.parallel,
                chunk_size=params.chunk_size,
                extra_flags=extra_flags,
            )
            if plan.get("errors"):
                raise RuntimeError("Alignment preflight failed: " + "; ".join(plan["errors"]))
            plan_path = subdir / "alignment_prepare_plan.json"
            _write_json(plan_path, plan)
            plan_paths.append(plan_path)
            synthetic_rows.extend(
                _dry_run_alignment_rows(rows, outdir=subdir, aligner=aligner, ref_fa=params.genome_fasta)
            )
        write_alignment_manifest_tsv(synthetic_rows, detector_manifest)
        _write_alignment_summary(
            root_dirs["align"] / "alignment_prepare_summary.json",
            source_rows=source_rows,
            manifest_paths=[detector_manifest],
            dry_run=True,
        )
        run_detector_alignment_manifest(
            detector=det,
            manifest=detector_manifest,
            outdir=root_dirs["ciri3"],
            ref_fa=params.genome_fasta,
            gtf=params.gtf,
            threads=params.threads,
            parallel=params.parallel,
            chunk_size=params.chunk_size,
            dry_run=True,
            fail_fast=params.fail_fast,
        )
        for plan_path in plan_paths:
            _print_dry_run_commands(progress, plan_path, root_dirs["ciri3"] / "detector_alignment_plan.json")
        _write_json(
            root_dirs["logs"] / "run_ciri3_summary.json",
            {
                "workflow": "run-ciri3",
                "protocols": sorted(set(protocol_by_cell.values())),
                "strandedness": sorted({value for value in strandedness_by_cell.values() if value}),
                "manifest": str(params.manifest.resolve()),
                "dry_run": True,
            },
        )
        return

        for aligner, subset_manifest, rows in subset_manifests:
            subdir = root_dirs["align"] / aligner.replace("-", "_")
            extra_flags = _star_extra_flags(rows).format(star_index=str(params.star_index.resolve())) if aligner == "star" else ""
            manifest_path = prepare_alignment_cache(
            manifest=subset_manifest,
            outdir=subdir,
            aligner=aligner,
            ref_fa=params.genome_fasta,
            detector_hint="ciri3",
            threads=params.threads,
            parallel=params.parallel,
            chunk_size=params.chunk_size,
            extra_flags=extra_flags,
            output_format="sam",
            fail_fast=params.fail_fast,
        )
        manifest_paths.append(manifest_path)

    _merge_alignment_manifests(manifest_paths, detector_manifest)
    _write_alignment_summary(
        root_dirs["align"] / "alignment_prepare_summary.json",
        source_rows=source_rows,
        manifest_paths=manifest_paths + [detector_manifest],
        dry_run=False,
    )
    run_detector_alignment_manifest(
        detector=det,
        manifest=detector_manifest,
        outdir=root_dirs["ciri3"],
        ref_fa=params.genome_fasta,
        gtf=params.gtf,
        threads=params.threads,
        parallel=params.parallel,
        chunk_size=params.chunk_size,
        fail_fast=params.fail_fast,
    )
    collect_matrix(
        str(root_dirs["ciri3"]),
        str(root_dirs["matrix"] / "circ_counts.mtx"),
        str(root_dirs["matrix"] / "circ_index.txt"),
        str(root_dirs["matrix"] / "cell_index.txt"),
        min_count_per_cell=1,
    )
    _write_json(
        root_dirs["logs"] / "run_ciri3_summary.json",
        {
            "workflow": "run-ciri3",
            "protocols": sorted(set(protocol_by_cell.values())),
            "strandedness": sorted({value for value in strandedness_by_cell.values() if value}),
            "manifest": str(params.manifest.resolve()),
            "alignment_manifest": str(detector_manifest.resolve()),
            "dry_run": False,
        },
    )
