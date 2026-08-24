from __future__ import annotations

import csv
import json
import re
import shlex
import shutil
import subprocess
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import pandas as pd

from circyto import __version__
from circyto.multimodal.sync import MUDATA_SYNC_POLICY, mudata_from_modalities, read_h5mu, write_h5mu
from circyto.pipeline.annotate_host_gene import normalize_host_gene_annotations
from circyto.pipeline.scomatic_normalization import (
    SCOMATIC_CANDIDATE_COLUMNS,
    TERMINOLOGY_NOTE,
    normalize_scomatic_results,
)
from circyto.pipeline.workflow_reporting import sanitize_frame_for_anndata, write_json

try:
    import anndata as ad  # type: ignore

    HAS_ANNDATA = True
except Exception:
    HAS_ANNDATA = False

try:
    import mudata as mu  # type: ignore

    HAS_MUDATA = True
except Exception:
    HAS_MUDATA = False


SUPPORTED_PROTOCOLS = {
    "smartseq2": "smartseq2",
    "smart-seq2": "smartseq2",
    "smart_seq2": "smartseq2",
    "smartseq3": "smartseq3",
    "smart-seq3": "smartseq3",
    "smart_seq3": "smartseq3",
    "ramda": "ramda",
    "ramda-seq": "ramda",
    "ramda_seq": "ramda",
    "shinramda": "shinramda",
    "shin-ramda": "shinramda",
    "shin_ramda": "shinramda",
    "scrr": "scrr-rna",
    "scrr-rna": "scrr-rna",
    "scrr_rna": "scrr-rna",
    "scrr-rna-arm": "scrr-rna",
}

CELL_ID_ALIASES = ("cell_id", "cell", "sample", "sample_id", "barcode", "Index")
BAM_ALIASES = ("bam", "bam_path", "alignment_bam", "bam_file")
SAM_ALIASES = ("sam", "sam_path", "alignment_sam", "sam_file")
CELL_TYPE_ALIASES = ("cell_type", "Cell_type", "celltype", "cell_label", "cluster", "annotation")


def normalize_full_length_protocol(protocol: str) -> str:
    key = str(protocol).strip().lower()
    if key not in SUPPORTED_PROTOCOLS:
        allowed = ", ".join(sorted(set(SUPPORTED_PROTOCOLS.values())))
        raise ValueError(f"Unsupported SComatic full-length protocol {protocol!r}. Accepted normalized values: {allowed}")
    return SUPPORTED_PROTOCOLS[key]


def _read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", keep_default_na=False)


def _find_column(df: pd.DataFrame, aliases: Sequence[str]) -> str | None:
    lookup = {str(column).strip().lower(): str(column) for column in df.columns}
    for alias in aliases:
        match = lookup.get(str(alias).strip().lower())
        if match is not None:
            return match
    return None


def _require_column(df: pd.DataFrame, aliases: Sequence[str], label: str, source: Path) -> str:
    column = _find_column(df, aliases)
    if column is None:
        raise ValueError(f"{source} is missing a {label} column. Tried: {', '.join(aliases)}")
    return column


def _resolve_path(value: str, *, base: Path) -> Path:
    path = Path(str(value).strip()).expanduser()
    if not path.is_absolute():
        path = base / path
    return path.resolve(strict=False)


def _safe_id(value: str, fallback: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("._-")
    return safe or fallback


def _safe_cell_type(value: str, fallback: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9_-]+", "_", str(value)).strip("_-")
    safe = re.sub(r"_+", "_", safe)
    return safe or fallback


def _load_cell_metadata(path: Path | None) -> tuple[dict[str, dict[str, str]], str | None]:
    if path is None:
        return {}, None
    df = _read_tsv(path)
    cell_col = _require_column(df, CELL_ID_ALIASES, "cell identifier", path)
    cell_type_col = _find_column(df, CELL_TYPE_ALIASES)
    metadata: dict[str, dict[str, str]] = {}
    for _, row in df.iterrows():
        cell_id = str(row[cell_col]).strip()
        if not cell_id:
            continue
        metadata[cell_id] = {str(column): str(row[column]).strip() for column in df.columns}
    return metadata, cell_type_col


def _cell_type_for_row(
    *,
    row: dict[str, str],
    cell_id: str,
    cell_type_column: str | None,
    default_cell_type: str | None,
    metadata_by_cell: dict[str, dict[str, str]],
    metadata_cell_type_column: str | None,
) -> str:
    raw = ""
    if cell_type_column:
        raw = str(row.get(cell_type_column, "")).strip()
        if not raw and cell_id in metadata_by_cell:
            raw = str(metadata_by_cell[cell_id].get(cell_type_column, "")).strip()
    if not raw and metadata_cell_type_column and cell_id in metadata_by_cell:
        raw = str(metadata_by_cell[cell_id].get(metadata_cell_type_column, "")).strip()
    if not raw and default_cell_type:
        raw = str(default_cell_type).strip()
    if not raw:
        raw = cell_id
    return _safe_cell_type(raw, fallback=_safe_cell_type(cell_id, fallback="unknown"))


def _write_annotation_table(rows: list[dict[str, str]], path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["Index", "Cell_type"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({"Index": row["cell_id"], "Cell_type": row["cell_type"]})


def _rewrite_sam_tags(sam_text: str, *, cell_id: str) -> tuple[str, dict[str, int]]:
    total_reads = 0
    missing_cb = 0
    replaced_cb = 0
    missing_nM = 0
    missing_NH = 0
    expected_cb = f"CB:Z:{cell_id}"
    out_lines: list[str] = []
    for raw_line in sam_text.splitlines():
        if raw_line.startswith("@") or not raw_line.strip():
            out_lines.append(raw_line)
            continue
        fields = raw_line.split("\t")
        total_reads += 1
        has_cb = False
        has_nM = False
        has_NH = False
        for idx in range(11, len(fields)):
            if fields[idx].startswith("CB:Z:"):
                has_cb = True
                if fields[idx] != expected_cb:
                    fields[idx] = expected_cb
                    replaced_cb += 1
            elif fields[idx].startswith("nM:"):
                has_nM = True
            elif fields[idx].startswith("NH:"):
                has_NH = True
        if not has_cb:
            fields.append(expected_cb)
            missing_cb += 1
        if not has_nM:
            missing_nM += 1
        if not has_NH:
            missing_NH += 1
        out_lines.append("\t".join(fields))
    return "\n".join(out_lines) + ("\n" if sam_text.endswith("\n") else ""), {
        "total_reads": total_reads,
        "missing_cb_before": missing_cb,
        "replaced_cb": replaced_cb,
        "missing_nM": missing_nM,
        "missing_NH": missing_NH,
    }


def _run_checked(cmd: list[str], *, input_data: bytes | str | None = None, text: bool = False) -> subprocess.CompletedProcess[Any]:
    result = subprocess.run(
        cmd,
        input=input_data,
        capture_output=True,
        text=text,
        check=False,
    )
    if result.returncode != 0:
        stderr = result.stderr if isinstance(result.stderr, str) else result.stderr.decode(errors="replace")
        raise RuntimeError(f"Command failed ({result.returncode}): {shlex.join(cmd)}\n{stderr}")
    return result


def _prepare_one_alignment(
    *,
    samtools: str,
    threads: int,
    alignment_path: Path,
    cell_id: str,
    sorted_bam: Path,
    tag_report: Path,
) -> dict[str, int]:
    view = _run_checked([samtools, "view", "-h", str(alignment_path)], text=True)
    rewritten, tag_counts = _rewrite_sam_tags(str(view.stdout), cell_id=cell_id)
    bam_bytes = _run_checked([samtools, "view", "-b", "-"], input_data=rewritten.encode("utf-8")).stdout
    _run_checked(
        [samtools, "sort", "-@", str(threads), "-o", str(sorted_bam), "-"],
        input_data=bam_bytes,
    )
    _run_checked([samtools, "index", str(sorted_bam)])
    pd.DataFrame([tag_counts]).to_csv(tag_report, sep="\t", index=False)
    return tag_counts


def _parse_one_cell_manifest(
    *,
    alignment_manifest_path: Path,
    outdir: Path,
    cell_metadata_path: Path | None,
    cell_type_column: str | None,
    default_cell_type: str | None,
) -> list[dict[str, str]]:
    manifest_df = _read_tsv(alignment_manifest_path)
    cell_col = _require_column(manifest_df, CELL_ID_ALIASES, "cell identifier", alignment_manifest_path)
    bam_col = _find_column(manifest_df, BAM_ALIASES)
    sam_col = _find_column(manifest_df, SAM_ALIASES)
    if bam_col is None and sam_col is None:
        raise ValueError(f"{alignment_manifest_path} must include a BAM or SAM path column.")
    if cell_type_column and cell_type_column not in manifest_df.columns:
        metadata_df = _read_tsv(cell_metadata_path) if cell_metadata_path is not None else pd.DataFrame()
        if cell_type_column not in metadata_df.columns:
            raise ValueError(f"--cell-type-column {cell_type_column!r} was not found in the manifest or metadata.")

    metadata_by_cell, metadata_cell_type_column = _load_cell_metadata(cell_metadata_path)
    seen_cells: set[str] = set()
    seen_file_ids: dict[str, int] = {}
    rows: list[dict[str, str]] = []
    per_cell_dir = outdir / "per_cell"
    tag_dir = outdir / "tag_reports"
    for idx, raw_row in manifest_df.iterrows():
        row = {str(column): str(raw_row[column]).strip() for column in manifest_df.columns}
        cell_id = str(row[cell_col]).strip()
        if not cell_id:
            raise ValueError(f"Empty cell_id at {alignment_manifest_path}:{idx + 2}")
        if re.search(r"\s", cell_id):
            raise ValueError(f"cell_id contains whitespace, which is unsafe for SComatic CB tags: {cell_id!r}")
        if cell_id in seen_cells:
            raise ValueError(f"Duplicate cell_id in alignment manifest: {cell_id}")
        seen_cells.add(cell_id)

        bam = row.get(bam_col, "") if bam_col is not None else ""
        sam = row.get(sam_col, "") if sam_col is not None else ""
        if bool(bam) == bool(sam):
            raise ValueError(f"Row for cell_id={cell_id} must provide exactly one BAM or SAM path.")
        alignment_path = _resolve_path(bam or sam, base=alignment_manifest_path.parent)
        if not alignment_path.exists():
            raise FileNotFoundError(f"Alignment file not found for cell_id={cell_id}: {alignment_path}")

        cell_type = _cell_type_for_row(
            row=row,
            cell_id=cell_id,
            cell_type_column=cell_type_column,
            default_cell_type=default_cell_type,
            metadata_by_cell=metadata_by_cell,
            metadata_cell_type_column=metadata_cell_type_column,
        )
        base_file_id = _safe_id(cell_id, fallback=f"cell_{idx + 1}")
        count = seen_file_ids.get(base_file_id, 0) + 1
        seen_file_ids[base_file_id] = count
        file_id = base_file_id if count == 1 else f"{base_file_id}_{count}"
        rows.append(
            {
                "file_id": file_id,
                "cell_id": cell_id,
                "alignment_path": str(alignment_path),
                "alignment_kind": "bam" if bam else "sam",
                "cell_type": cell_type,
                "sorted_bam": str(per_cell_dir / f"{file_id}.scomatic.sorted.bam"),
                "tag_report": str(tag_dir / f"{file_id}.tags.tsv"),
            }
        )
    if not rows:
        raise ValueError(f"Alignment manifest contains 0 rows: {alignment_manifest_path}")
    return rows


def _shell_args(value: str | None) -> list[str]:
    return shlex.split(value or "")


def _expand_args(args: Sequence[str], mapping: dict[str, str]) -> list[str]:
    return [str(arg).format(**mapping) for arg in args]


def _shlex_join(cmd: Sequence[str]) -> str:
    return shlex.join([str(item) for item in cmd])


def prepare_scomatic_input(
    *,
    outdir: Path,
    protocol: str,
    alignment_manifest_path: Path | None = None,
    merged_bam_path: Path | None = None,
    cell_metadata_path: Path | None = None,
    reference_fasta_path: Path | None = None,
    sample_id: str = "scomatic",
    threads: int = 1,
    cell_type_column: str | None = None,
    default_cell_type: str | None = None,
    samtools: str = "samtools",
) -> dict[str, Any]:
    normalized_protocol = normalize_full_length_protocol(protocol)
    if threads < 1:
        raise ValueError("--threads must be a positive integer.")
    if (alignment_manifest_path is None) == (merged_bam_path is None):
        raise ValueError("Provide exactly one of --alignment-manifest or --merged-bam.")

    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "per_cell").mkdir(exist_ok=True)
    (outdir / "tag_reports").mkdir(exist_ok=True)
    (outdir / "merged").mkdir(exist_ok=True)

    annotation_tsv = outdir / "cell_annotations.tsv"
    rows_tsv = outdir / "adapter_rows.tsv"
    prepared_bams_tsv = outdir / "scomatic_prepared_bams.tsv"
    merged_bam_out: Path
    tag_summary = {
        "total_reads": 0,
        "missing_cb_before_injection": 0,
        "replaced_cb_not_matching_cell_id": 0,
        "missing_nM": 0,
        "missing_NH": 0,
    }

    if alignment_manifest_path is not None:
        samtools_path = shutil.which(samtools)
        if samtools_path is None:
            raise RuntimeError(f"samtools not found on PATH: {samtools}")
        mode = "one-cell-per-bam"
        rows = _parse_one_cell_manifest(
            alignment_manifest_path=alignment_manifest_path,
            outdir=outdir,
            cell_metadata_path=cell_metadata_path,
            cell_type_column=cell_type_column,
            default_cell_type=default_cell_type,
        )
        pd.DataFrame(rows).to_csv(rows_tsv, sep="\t", index=False)
        _write_annotation_table(rows, annotation_tsv)
        sorted_bams: list[Path] = []
        for row in rows:
            sorted_bam = Path(row["sorted_bam"])
            tag_report = Path(row["tag_report"])
            sorted_bam.parent.mkdir(parents=True, exist_ok=True)
            tag_report.parent.mkdir(parents=True, exist_ok=True)
            counts = _prepare_one_alignment(
                samtools=samtools_path,
                threads=threads,
                alignment_path=Path(row["alignment_path"]),
                cell_id=row["cell_id"],
                sorted_bam=sorted_bam,
                tag_report=tag_report,
            )
            tag_summary["total_reads"] += int(counts["total_reads"])
            tag_summary["missing_cb_before_injection"] += int(counts["missing_cb_before"])
            tag_summary["replaced_cb_not_matching_cell_id"] += int(counts["replaced_cb"])
            tag_summary["missing_nM"] += int(counts["missing_nM"])
            tag_summary["missing_NH"] += int(counts["missing_NH"])
            sorted_bams.append(sorted_bam)

        merged_bam_out = outdir / "merged" / f"{_safe_id(sample_id, fallback='scomatic')}.scomatic.bam"
        if len(sorted_bams) == 1:
            shutil.copyfile(sorted_bams[0], merged_bam_out)
        else:
            _run_checked([samtools_path, "merge", "-@", str(threads), "-f", str(merged_bam_out)] + [str(path) for path in sorted_bams])
        _run_checked([samtools_path, "index", str(merged_bam_out)])
    else:
        mode = "merged-bam"
        assert merged_bam_path is not None
        if not merged_bam_path.exists():
            raise FileNotFoundError(f"Merged BAM not found: {merged_bam_path}")
        metadata_by_cell, metadata_cell_type_column = _load_cell_metadata(cell_metadata_path)
        if not metadata_by_cell:
            raise ValueError("--cell-metadata is required for --merged-bam mode and must contain at least one cell.")
        rows = []
        for idx, cell_id in enumerate(sorted(metadata_by_cell)):
            cell_type = _cell_type_for_row(
                row=metadata_by_cell[cell_id],
                cell_id=cell_id,
                cell_type_column=cell_type_column,
                default_cell_type=default_cell_type,
                metadata_by_cell=metadata_by_cell,
                metadata_cell_type_column=metadata_cell_type_column,
            )
            rows.append(
                {
                    "file_id": _safe_id(cell_id, fallback=f"cell_{idx + 1}"),
                    "cell_id": cell_id,
                    "alignment_path": str(merged_bam_path.resolve()),
                    "alignment_kind": "merged_bam",
                    "cell_type": cell_type,
                    "sorted_bam": "",
                    "tag_report": "",
                }
            )
        pd.DataFrame(rows).to_csv(rows_tsv, sep="\t", index=False)
        _write_annotation_table(rows, annotation_tsv)
        merged_bam_out = merged_bam_path.resolve()

    prepared_bams = pd.DataFrame(
        [
            {
                "sample_id": sample_id,
                "mode": mode,
                "bam": str(merged_bam_out),
                "cell_annotations": str(annotation_tsv.resolve()),
                "protocol": normalized_protocol,
                "reference_fasta": str(reference_fasta_path.resolve()) if reference_fasta_path is not None else "",
            }
        ]
    )
    prepared_bams.to_csv(prepared_bams_tsv, sep="\t", index=False)

    limitations = [
        "SComatic outputs are RNA-derived candidate variant signals unless orthogonal DNA validation exists.",
        "Homogeneous datasets or one-cell-type annotations may produce no Step2 candidate calls.",
        "Without biologically meaningful cell types, Step1/Step2 operate on technical groups rather than validated cell-type contrasts.",
        "Without a matched Panel of Normals, downstream calls should remain exploratory.",
    ]
    summary = {
        "command_name": "circyto prepare-scomatic-input",
        "circyto_version": __version__,
        "description": "Prepared full-length RNA alignments and cell annotations for external SComatic interoperability.",
        "protocol": normalized_protocol,
        "supported_protocols": sorted(set(SUPPORTED_PROTOCOLS.values())),
        "mode": mode,
        "sample_id": sample_id,
        "runs_scomatic": False,
        "scomatic_result_terminology": "RNA-derived candidate variant signals",
        "terminology_note": TERMINOLOGY_NOTE,
        "alignment_manifest": str(alignment_manifest_path.resolve()) if alignment_manifest_path is not None else None,
        "source_merged_bam": str(merged_bam_path.resolve()) if merged_bam_path is not None else None,
        "cell_metadata": str(cell_metadata_path.resolve()) if cell_metadata_path is not None else None,
        "reference_fasta": str(reference_fasta_path.resolve()) if reference_fasta_path is not None else None,
        "outdir": str(outdir.resolve()),
        "n_cells": int(len(rows)),
        "merged_bam": str(merged_bam_out.resolve()),
        "merged_bam_index": str(Path(str(merged_bam_out) + ".bai").resolve()) if mode == "one-cell-per-bam" else None,
        "cell_annotation_tsv": str(annotation_tsv.resolve()),
        "adapter_rows_tsv": str(rows_tsv.resolve()),
        "prepared_bams_tsv": str(prepared_bams_tsv.resolve()),
        "tag_summary": tag_summary,
        "scomatic_tag_contract": {
            "CB": "required by SComatic; injected or normalized to cell_id for one-cell-per-BAM input",
            "nM": "preserved and counted; required when running SComatic with --max_nM",
            "NH": "preserved and counted; required when running SComatic with --max_NH",
        },
        "limitations": limitations,
    }
    summary_path = outdir / "prepare_scomatic_input_summary.json"
    summary["output_summary_json"] = str(summary_path.resolve())
    write_json(summary_path, summary)
    return summary


def run_scomatic(
    *,
    prepared_dir: Path,
    outdir: Path,
    scomatic_dir: Path,
    reference_fasta_path: Path,
    threads: int = 1,
    python_executable: str = "python",
    basecellcounter_script: Path | None = None,
    step1_script: Path | None = None,
    step2_script: Path | None = None,
    basecellcounter_args: str | None = None,
    step1_args: str | None = None,
    step2_args: str | None = None,
    run_step1: bool = False,
    run_step2: bool = False,
    execute: bool = False,
) -> dict[str, Any]:
    if threads < 1:
        raise ValueError("--threads must be a positive integer.")
    prepared_summary_path = prepared_dir / "prepare_scomatic_input_summary.json"
    if not prepared_summary_path.exists():
        raise FileNotFoundError(f"Missing prepared input summary: {prepared_summary_path}")
    prepared_summary = json.loads(prepared_summary_path.read_text(encoding="utf-8"))
    merged_bam = str(prepared_summary.get("merged_bam", ""))
    cell_annotations = str(prepared_summary.get("cell_annotation_tsv", ""))
    if not merged_bam:
        raise ValueError(f"{prepared_summary_path} does not record merged_bam.")

    outdir.mkdir(parents=True, exist_ok=True)
    basecellcounter_out = outdir / "basecellcounter"
    step1_out = outdir / "step1"
    step2_out = outdir / "step2"
    basecellcounter_out.mkdir(exist_ok=True)
    if run_step1:
        step1_out.mkdir(exist_ok=True)
    if run_step2:
        step2_out.mkdir(exist_ok=True)

    base_script = basecellcounter_script or scomatic_dir / "scripts" / "BaseCellCounter" / "BaseCellCounter.py"
    s1_script = step1_script or scomatic_dir / "scripts" / "BaseCellCalling" / "BaseCellCalling.step1.py"
    s2_script = step2_script or scomatic_dir / "scripts" / "BaseCellCalling" / "BaseCellCalling.step2.py"
    mapping = {
        "merged_bam": merged_bam,
        "cell_annotations": cell_annotations,
        "reference_fasta": str(reference_fasta_path.resolve()),
        "outdir": str(outdir.resolve()),
        "basecellcounter_out": str(basecellcounter_out.resolve()),
        "step1_out": str(step1_out.resolve()),
        "step2_out": str(step2_out.resolve()),
        "threads": str(threads),
    }

    base_args = _expand_args(
        _shell_args(basecellcounter_args)
        or ["--bam", "{merged_bam}", "--ref", "{reference_fasta}", "--out_folder", "{basecellcounter_out}", "--nprocs", "{threads}"],
        mapping,
    )
    commands: list[dict[str, Any]] = [
        {
            "stage": "BaseCellCounter",
            "script": str(base_script),
            "script_exists": base_script.exists(),
            "argv": [python_executable, str(base_script)] + base_args,
        }
    ]
    if run_step1:
        s1_args = _expand_args(_shell_args(step1_args) or ["--infile", "{basecellcounter_out}", "--out_folder", "{step1_out}"], mapping)
        commands.append(
            {
                "stage": "Step1",
                "script": str(s1_script),
                "script_exists": s1_script.exists(),
                "argv": [python_executable, str(s1_script)] + s1_args,
            }
        )
    if run_step2:
        s2_args = _expand_args(_shell_args(step2_args) or ["--infile", "{step1_out}", "--out_folder", "{step2_out}"], mapping)
        commands.append(
            {
                "stage": "Step2",
                "script": str(s2_script),
                "script_exists": s2_script.exists(),
                "argv": [python_executable, str(s2_script)] + s2_args,
            }
        )

    plan_path = outdir / "scomatic_run_plan.sh"
    plan_lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "",
        "# Generated by circyto run-scomatic.",
        "# Outputs remain RNA-derived candidate variant signals unless orthogonal DNA validation exists.",
    ]
    plan_lines.extend(_shlex_join(command["argv"]) for command in commands)
    plan_path.write_text("\n".join(plan_lines) + "\n", encoding="utf-8")
    plan_path.chmod(0o755)

    execution_results: list[dict[str, Any]] = []
    if execute:
        for command in commands:
            if not Path(command["script"]).exists():
                raise FileNotFoundError(f"SComatic script not found for {command['stage']}: {command['script']}")
            result = subprocess.run(command["argv"], capture_output=True, text=True, check=False)
            log_base = outdir / f"{str(command['stage']).lower()}_execution"
            (Path(str(log_base) + ".stdout.log")).write_text(result.stdout, encoding="utf-8")
            (Path(str(log_base) + ".stderr.log")).write_text(result.stderr, encoding="utf-8")
            execution_results.append(
                {
                    "stage": command["stage"],
                    "returncode": int(result.returncode),
                    "stdout_log": str(Path(str(log_base) + ".stdout.log").resolve()),
                    "stderr_log": str(Path(str(log_base) + ".stderr.log").resolve()),
                }
            )
            if result.returncode != 0:
                raise RuntimeError(f"SComatic {command['stage']} failed with exit code {result.returncode}.")

    summary = {
        "command_name": "circyto run-scomatic",
        "circyto_version": __version__,
        "description": "Generated or executed an external SComatic command plan for full-length RNA candidate-signal interoperability.",
        "execute": bool(execute),
        "prepared_dir": str(prepared_dir.resolve()),
        "scomatic_dir": str(scomatic_dir.resolve()),
        "reference_fasta": str(reference_fasta_path.resolve()),
        "merged_bam": merged_bam,
        "cell_annotations": cell_annotations,
        "commands": commands,
        "command_plan": str(plan_path.resolve()),
        "execution_results": execution_results,
        "terminology_note": TERMINOLOGY_NOTE,
        "limitations": [
            "External SComatic parameters, reference, and Panel of Normals determine candidate quality.",
            "Without cell-type contrasts or PoN filtering, Step2 may produce no candidates or exploratory-only candidates.",
            "circyto records and imports SComatic outputs as RNA-derived candidate variant signals.",
        ],
    }
    summary_path = outdir / "run_scomatic_summary.json"
    summary["output_summary_json"] = str(summary_path.resolve())
    write_json(summary_path, summary)
    return summary


def import_scomatic(
    *,
    scomatic_output_paths: Sequence[Path],
    outdir: Path,
    cell_annotation_path: Path | None = None,
    provenance_metadata_path: Path | None = None,
) -> dict[str, Any]:
    summary = normalize_scomatic_results(
        scomatic_output_paths=scomatic_output_paths,
        outdir=outdir,
        cell_annotation_path=cell_annotation_path,
        provenance_metadata_path=provenance_metadata_path,
    )
    summary = dict(summary)
    summary["command_name"] = "circyto import-scomatic"
    summary["description"] = "Imported SComatic Step1/Step2 or candidate tables into circyto's candidate-signal schema."
    summary["output_import_scomatic_summary"] = str((outdir / "import_scomatic_summary.json").resolve())
    write_json(outdir / "import_scomatic_summary.json", summary)
    return summary


def _candidate_feature_id(row: pd.Series) -> str:
    return f"{row['chrom']}:{row['pos']}:{row['ref']}>{row['alt']}"


def _join_unique(values: pd.Series) -> str:
    unique = sorted({str(value) for value in values if str(value) != ""})
    return ";".join(unique)


def _candidate_snv_anndata(candidate_df: pd.DataFrame, obs_names: list[str]):
    if not HAS_ANNDATA:
        raise RuntimeError("anndata is required for merge-scomatic")
    missing = [column for column in SCOMATIC_CANDIDATE_COLUMNS if column not in candidate_df.columns]
    if missing:
        raise ValueError(f"SComatic candidate summary is missing required columns: {', '.join(missing)}")

    df = candidate_df.copy()
    df["cell_id"] = df["cell_id"].astype(str)
    df["feature_id"] = df.apply(_candidate_feature_id, axis=1)
    features = sorted(df["feature_id"].astype(str).unique().tolist())
    obs_index = {cell_id: idx for idx, cell_id in enumerate(obs_names)}
    var_index = {feature_id: idx for idx, feature_id in enumerate(features)}
    read_support = np.zeros((len(obs_names), len(features)), dtype=np.float32)
    vaf = np.zeros((len(obs_names), len(features)), dtype=np.float32)
    presence = np.zeros((len(obs_names), len(features)), dtype=np.int8)
    for _, row in df.iterrows():
        cell_id = str(row["cell_id"])
        if cell_id not in obs_index:
            continue
        obs_idx = obs_index[cell_id]
        var_idx = var_index[str(row["feature_id"])]
        read_support[obs_idx, var_idx] += float(row["read_support"])
        vaf[obs_idx, var_idx] = max(vaf[obs_idx, var_idx], float(row["vaf"]))
        presence[obs_idx, var_idx] = 1

    obs = pd.DataFrame({"cell_id": obs_names}, index=obs_names)
    var = (
        df.groupby("feature_id", as_index=False)
        .agg(
            chrom=("chrom", "first"),
            pos=("pos", "first"),
            ref=("ref", "first"),
            alt=("alt", "first"),
            gene=("gene", _join_unique),
            filter_status_values=("filter_status", _join_unique),
            candidate_variant_class=("candidate_variant_class", _join_unique),
            caller=("caller", _join_unique),
        )
        .set_index("feature_id", drop=False)
        .reindex(features)
    )
    var.index.name = None
    adata = ad.AnnData(
        X=read_support,
        obs=sanitize_frame_for_anndata(obs),
        var=sanitize_frame_for_anndata(var),
    )
    adata.layers["vaf"] = vaf
    adata.layers["presence"] = presence
    adata.uns["circyto"] = {
        "modality_name": "candidate_snv",
        "command_name": "circyto merge-scomatic",
        "circyto_version": __version__,
        "X_semantics": "cell x candidate variant read_support from RNA-derived SComatic output",
        "layers": {
            "vaf": "cell x candidate variant maximum VAF",
            "presence": "binary candidate signal presence",
        },
        "terminology_note": TERMINOLOGY_NOTE,
    }
    return adata


def merge_scomatic(
    *,
    input_h5mu: Path,
    scomatic_candidates_path: Path,
    output_h5mu: Path,
    summary_json: Path | None = None,
    allow_partial: bool = False,
    overwrite: bool = False,
) -> dict[str, Any]:
    if not HAS_ANNDATA:
        raise RuntimeError("anndata is required for merge-scomatic")
    if not HAS_MUDATA:
        raise RuntimeError("mudata is not installed; install circyto[mudata] or pip install mudata")
    if output_h5mu.exists() and not overwrite:
        raise ValueError(f"Output already exists: {output_h5mu}. Use --overwrite to replace it.")
    if summary_json is None:
        summary_json = output_h5mu.with_suffix(".summary.json")
    if summary_json.exists() and not overwrite:
        raise ValueError(f"Summary already exists: {summary_json}. Use --overwrite to replace it.")

    mdata = read_h5mu(str(input_h5mu))
    candidate_df = pd.read_csv(scomatic_candidates_path, sep="\t", keep_default_na=False)
    candidate_cells = set(candidate_df["cell_id"].astype(str)) if "cell_id" in candidate_df.columns else set()
    existing_obs = [str(value) for value in mdata.obs_names.astype(str).tolist()]
    missing_cells = sorted(candidate_cells - set(existing_obs))
    if missing_cells and not allow_partial:
        preview = ", ".join(missing_cells[:10])
        raise ValueError(
            "SComatic candidate cell IDs are absent from the input MuData obs: "
            f"{preview}. Use --allow-partial for cell-type-level or partial candidate tables."
        )
    obs_names = existing_obs if not allow_partial else sorted(set(existing_obs) | candidate_cells)
    candidate_adata = _candidate_snv_anndata(candidate_df, obs_names)
    modalities = {modality: adata_obj.copy() for modality, adata_obj in mdata.mod.items()}
    if "circ" in modalities:
        modalities["circ"].var = normalize_host_gene_annotations(modalities["circ"].var)
    modalities["candidate_snv"] = candidate_adata
    merged = mudata_from_modalities(modalities)
    merged.uns.update(mdata.uns)
    merged.uns["circyto_scomatic"] = {
        "command_name": "circyto merge-scomatic",
        "circyto_version": __version__,
        "mudata_sync_policy": MUDATA_SYNC_POLICY,
        "source_h5mu": str(input_h5mu.resolve()),
        "scomatic_candidates": str(scomatic_candidates_path.resolve()),
        "allow_partial": bool(allow_partial),
        "terminology_note": TERMINOLOGY_NOTE,
    }
    output_h5mu.parent.mkdir(parents=True, exist_ok=True)
    write_h5mu(merged, str(output_h5mu))
    summary = {
        "command_name": "circyto merge-scomatic",
        "description": "Merged RNA-derived SComatic candidate signals into MuData as candidate_snv.",
        "source_h5mu": str(input_h5mu.resolve()),
        "scomatic_candidates": str(scomatic_candidates_path.resolve()),
        "output_h5mu": str(output_h5mu.resolve()),
        "output_summary_json": str(summary_json.resolve()),
        "allow_partial": bool(allow_partial),
        "n_input_obs": int(len(existing_obs)),
        "n_candidate_cells": int(len(candidate_cells)),
        "n_missing_candidate_cells": int(len(missing_cells)),
        "missing_candidate_cells": missing_cells,
        "n_obs": int(merged.n_obs),
        "modalities": list(merged.mod.keys()),
        "modality_shapes": {key: [int(value.n_obs), int(value.n_vars)] for key, value in merged.mod.items()},
        "terminology_note": TERMINOLOGY_NOTE,
    }
    write_json(summary_json, summary)
    return summary
