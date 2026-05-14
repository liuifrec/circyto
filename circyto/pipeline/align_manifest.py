from __future__ import annotations

import csv
import hashlib
import json
import os
import shlex
import shutil
import subprocess
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional

from circyto.detectors import DetectorBase, DetectorRunInputs, DetectorResult
from circyto.detectors.base import get_detector_capabilities
from circyto.detectors.ciri3 import Ciri3Detector
from circyto.manifest.alignment import (
    AlignmentManifestRow,
    read_alignment_manifest_tsv,
    write_alignment_manifest_tsv,
)
from circyto.paths import resolve_manifest_path

VALID_READ_LAYOUTS = {"single-end", "paired-end"}
VALID_SUBSETS = {"failed", "missing", "stale", "incomplete", "all-failed-chunks"}


def _alignment_output_has_calls(path: Path, detector_name: str) -> bool:
    if not path.exists() or path.stat().st_size == 0:
        return False
    text = path.read_text(encoding="utf-8", errors="replace").splitlines()
    if detector_name in {"ciri-full", "ciri2", "ciri3"}:
        return len([line for line in text[1:] if line.strip()]) > 0
    return len([line for line in text if line.strip() and not line.startswith("#")]) > 0


def _count_alignment_output_rows(path: Path, detector_name: str) -> int | None:
    if not path.exists() or path.stat().st_size == 0:
        return 0
    text = path.read_text(encoding="utf-8", errors="replace").splitlines()
    if detector_name in {"ciri-full", "ciri2", "ciri3"}:
        return len([line for line in text[1:] if line.strip()])
    return len([line for line in text if line.strip() and not line.startswith("#")])


def _alignment_outcome_category(status: str, *, raw_rows: int | None, normalized_rows: int | None) -> str:
    if status == "failed":
        return "failed"
    if status == "skipped_existing":
        return "skipped-existing"
    if normalized_rows and normalized_rows > 0:
        return "success-non-empty"
    if raw_rows is not None and raw_rows > 0 and (normalized_rows or 0) == 0:
        return "success-normalized-empty"
    return "success-empty"


@dataclass(frozen=True)
class SourceManifestRow:
    cell_id: str
    read1: Optional[Path] = None
    read2: Optional[Path] = None
    bam: Optional[Path] = None
    platform: str = ""
    protocol: str = ""
    strandedness: str = ""
    library_id: str = ""
    group_id: str = ""
    declared_read_layout: str = ""
    extra: Optional[Dict[str, str]] = None

    @property
    def read_layout(self) -> str:
        if self.declared_read_layout in VALID_READ_LAYOUTS:
            return self.declared_read_layout
        return "paired-end" if self.read2 is not None else "single-end"

    @property
    def input_mode(self) -> str:
        return "alignment" if self.bam is not None and self.read1 is None else "fastq"

    @property
    def effective_read2(self) -> Optional[Path]:
        if self.read_layout != "paired-end":
            return None
        return self.read2


def _summary_path(outdir: Path, name: str) -> Path:
    return outdir / name


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _merge_cell_records(
    existing_summary: dict[str, Any] | None,
    new_records: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    merged: dict[str, dict[str, Any]] = {}
    for record in (existing_summary or {}).get("cells", []):
        cell_id = str(record.get("cell_id", "")).strip()
        if cell_id:
            merged[cell_id] = dict(record)
    for record in new_records:
        cell_id = str(record.get("cell_id", "")).strip()
        if cell_id:
            merged[cell_id] = dict(record)
    return [merged[cell_id] for cell_id in sorted(merged)]


def _read_manifest_rows_raw(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = list(reader.fieldnames or [])
        rows = [{k: "" if v is None else str(v) for k, v in row.items()} for row in reader]
    return header, rows


def _write_manifest_rows_raw(path: Path, *, header: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in header})


def _pick_col(row: dict[str, str], keys: tuple[str, ...]) -> Optional[str]:
    for key in keys:
        value = row.get(key)
        if key in row and value is not None and str(value).strip():
            return key
    return None


def _manifest_value(raw: dict[str, str], *keys: str) -> str:
    for key in keys:
        value = raw.get(key)
        if value is None:
            continue
        text = str(value).strip()
        if text:
            return text
    return ""


def _validated_read_layout(raw: dict[str, str], *, path: Path, line_number: int, cell_id: str) -> str:
    read_layout = (raw.get("read_layout") or "").strip()
    if read_layout and read_layout not in VALID_READ_LAYOUTS:
        raise ValueError(f"Invalid read_layout '{read_layout}' for cell_id={cell_id} at {path}:{line_number}")
    return read_layout


def read_source_manifest(path: Path, *, validate_files: bool = True) -> List[SourceManifestRow]:
    if not path.exists():
        raise FileNotFoundError(f"Manifest not found: {path}")

    rows: List[SourceManifestRow] = []
    seen_ids: set[str] = set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        if not fieldnames or ("cell_id" not in fieldnames and "sample_id" not in fieldnames):
            raise KeyError(f"Manifest missing required column 'cell_id' or 'sample_id': {path}")

        for i, raw in enumerate(reader, start=2):
            cell_id = (raw.get("cell_id") or raw.get("sample_id") or "").strip()
            if not cell_id:
                raise ValueError(f"Empty cell_id at {path}:{i}")
            if cell_id in seen_ids:
                raise ValueError(f"Duplicate cell_id '{cell_id}' at {path}:{i}")
            seen_ids.add(cell_id)

            read1_raw = _manifest_value(raw, "read1", "r1", "fastq_1")
            read2_raw = _manifest_value(raw, "read2", "r2", "fastq_2")
            read_layout = _validated_read_layout(raw, path=path, line_number=i, cell_id=cell_id)
            read1 = resolve_manifest_path(path, read1_raw) if read1_raw else None
            read2 = resolve_manifest_path(path, read2_raw) if read2_raw else None
            bam_raw = (raw.get("bam") or "").strip()
            bam = resolve_manifest_path(path, bam_raw) if bam_raw else None

            if read1 is None and bam is None:
                raise ValueError(f"Manifest row for cell_id={cell_id} must provide FASTQ or BAM input")
            if validate_files:
                if read1 is not None and not read1.exists():
                    raise FileNotFoundError(f"Manifest read1 not found for cell_id={cell_id}: {read1}")
                if read2 is not None and not read2.exists():
                    raise FileNotFoundError(f"Manifest read2 not found for cell_id={cell_id}: {read2}")
                if bam is not None and not bam.exists():
                    raise FileNotFoundError(f"Manifest bam not found for cell_id={cell_id}: {bam}")

            extras = {
                key: value
                for key, value in raw.items()
                if key
                not in {
                    "cell_id",
                    "sample_id",
                    "read1",
                    "read2",
                    "r1",
                    "r2",
                    "fastq_1",
                    "fastq_2",
                    "bam",
                    "platform",
                    "protocol",
                    "strandedness",
                    "library_id",
                    "group_id",
                    "read_layout",
                }
                and value not in (None, "")
            }
            rows.append(
                SourceManifestRow(
                    cell_id=cell_id,
                    read1=read1,
                    read2=read2,
                    bam=bam,
                    platform=(raw.get("platform") or "").strip(),
                    protocol=(raw.get("protocol") or raw.get("platform") or "").strip(),
                    strandedness=(raw.get("strandedness") or "").strip(),
                    library_id=(raw.get("library_id") or "").strip(),
                    group_id=(raw.get("group_id") or raw.get("library_id") or "").strip(),
                    declared_read_layout=read_layout,
                    extra=extras or None,
                )
            )

    if not rows:
        raise ValueError(f"Manifest contains 0 data rows: {path}")
    return rows


def _alignment_provenance_path(alignment_path: Path) -> Path:
    return Path(str(alignment_path) + ".provenance.json")


def _detector_output_provenance_path(output_path: Path) -> Path:
    return Path(str(output_path) + ".provenance.json")


def _chunk_dir(outdir: Path) -> Path:
    return outdir / "chunks"


def _chunk_summary_path(outdir: Path, chunk_index: int) -> Path:
    return _chunk_dir(outdir) / f"chunk_{chunk_index:05d}.json"


def _next_chunk_index(outdir: Path) -> int:
    max_index = 0
    for path in _chunk_dir(outdir).glob("chunk_*.json"):
        try:
            max_index = max(max_index, int(path.stem.split("_", 1)[1]))
        except (IndexError, ValueError):
            continue
    return max_index + 1


def _load_json(path: Path) -> dict[str, Any] | None:
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None


def _alignment_cache_key(
    row: SourceManifestRow,
    *,
    aligner: str,
    ref_fa: Path | None,
    detector_hint: str | None,
    threads: int,
    extra_flags: str,
) -> str:
    payload = {
        "cell_id": row.cell_id,
        "read1": str(row.read1.resolve()) if row.read1 else None,
        "read2": str(row.effective_read2.resolve()) if row.effective_read2 else None,
        "bam": str(row.bam.resolve()) if row.bam else None,
        "aligner": aligner,
        "ref_fa": str(ref_fa.resolve()) if ref_fa else None,
        "detector_hint": detector_hint,
        "read_layout": row.read_layout,
        "threads": threads,
        "extra_flags": extra_flags,
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode("utf-8")).hexdigest()[:16]


def _copy_or_link(src: Path, dest: Path, *, link_mode: str) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() or dest.is_symlink():
        dest.unlink()
    if link_mode == "copy":
        shutil.copy2(src, dest)
        return
    try:
        # Use an absolute source so manifest rows derived from dest.resolve()
        # point at the cached artifact rather than a doubly-prefixed relative path.
        os.symlink(src.resolve(), dest)
    except OSError:
        shutil.copy2(src, dest)


def _tool_exists(name: str) -> bool:
    return shutil.which(name) is not None


def _tail_text(path: Path, max_lines: int = 40) -> str:
    try:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
        return ""
    return "\n".join(lines[-max_lines:])


def _shell_join(parts: list[str]) -> str:
    return shlex.join(parts)


def _append_log_line(path: Path, line: str, *, mode: str = "a") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open(mode, encoding="utf-8") as handle:
        handle.write(line.rstrip("\n") + "\n")


def _emit_prepare_progress(*, cell_id: str, stage: str, log_path: Path, detail: str = "") -> None:
    line = f"[circyto] cell={cell_id} stage={stage}"
    if detail:
        line += f" {detail}"
    print(line, flush=True)
    _append_log_line(log_path, line)


def _build_bwa_mem_command(
    *,
    row: SourceManifestRow,
    ref_fa: Path,
    threads: int,
    extra_flags: str,
) -> list[str]:
    cmd = ["bwa", "mem", "-t", str(threads)]
    if extra_flags.strip():
        cmd.extend(shlex.split(extra_flags))
    cmd.extend([str(ref_fa), str(row.read1)])
    if row.effective_read2 is not None:
        cmd.append(str(row.effective_read2))
    return cmd


def _alignment_failure_details(
    *,
    row: SourceManifestRow,
    command: list[str],
    log_path: Path,
    exit_code: int | None,
    stderr_tail_lines: int = 40,
) -> dict[str, Any]:
    return {
        "read_layout": row.read_layout,
        "read1": str(row.read1.resolve()) if row.read1 else None,
        "read2": str(row.effective_read2.resolve()) if row.effective_read2 else None,
        "command": _shell_join(command),
        "exit_code": exit_code,
        "log_path": str(log_path.resolve()),
        "stderr_tail": _tail_text(log_path, max_lines=stderr_tail_lines),
    }


def _format_alignment_failure(prefix: str, details: dict[str, Any], *, extra: str = "") -> str:
    segments = [
        prefix,
        f"read_layout={details.get('read_layout')}",
        f"command={details.get('command')}",
        f"read1={details.get('read1')}",
        f"read2={details.get('read2') or '<none>'}",
        f"exit={details.get('exit_code')}",
        f"log={details.get('log_path')}",
    ]
    if extra:
        segments.append(extra)
    message = "; ".join(segments)
    stderr_tail = str(details.get("stderr_tail") or "").strip()
    if stderr_tail:
        message += f"\n--- log tail ---\n{stderr_tail}"
    return message


def _row_field(row: AlignmentManifestRow, name: str, default: str = "") -> str:
    value = getattr(row, name, "") or ""
    if value:
        return str(value)
    extra = row.extra or {}
    return str(extra.get(name, default) or default)


def _alignment_output_path(cell_cache_dir: Path, output_format: str) -> Path:
    suffix = ".sam" if output_format == "sam" else ".bam"
    return cell_cache_dir / f"alignment{suffix}"


def _effective_output_format(*, detector_hint: str | None, aligner: str, output_format: str) -> str:
    if detector_hint == "ciri3" and aligner in {"bwa-mem", "star"}:
        return "sam"
    return output_format


def _staged_alignment_bucket(*, aligner: str, output_format: str) -> str:
    if aligner == "bwa-mem":
        return "bwa_mem" if output_format == "sam" else "bwa_mem_bam"
    if aligner == "star":
        return "star"
    if aligner == "reuse-existing":
        return "reused_alignment"
    return aligner.replace("-", "_")


def _star_genome_dir(extra_flags: str) -> str | None:
    tokens = shlex.split(extra_flags.strip()) if extra_flags.strip() else []
    for i, token in enumerate(tokens):
        if token == "--genomeDir" and i + 1 < len(tokens):
            return tokens[i + 1]
        if token.startswith("--genomeDir="):
            return token.split("=", 1)[1]
    return os.environ.get("CIRCYTO_STAR_GENOME_DIR")


def _flag_present(tokens: list[str], flag: str) -> bool:
    if flag in tokens:
        return True
    prefix = f"{flag}="
    return any(token.startswith(prefix) for token in tokens)


def _star_read_files_command(extra_flags: str) -> str | None:
    tokens = shlex.split(extra_flags.strip()) if extra_flags.strip() else []
    for i, token in enumerate(tokens):
        if token == "--readFilesCommand" and i + 1 < len(tokens):
            return tokens[i + 1]
        if token.startswith("--readFilesCommand="):
            return token.split("=", 1)[1]
    return None


def _star_tmp_dir_base() -> str | None:
    return os.environ.get("CIRCYTO_STAR_TMPDIR") or os.environ.get("TMPDIR")


def _effective_alignment_extra_flags(*, detector_hint: str | None, aligner: str, extra_flags: str) -> str:
    if detector_hint == "ciri3" and aligner == "star":
        tokens = shlex.split(extra_flags.strip()) if extra_flags.strip() else []
        required_pairs = [
            ("--outSAMtype", ["SAM"]),
            ("--outReadsUnmapped", ["Fastx"]),
            ("--outSJfilterOverhangMin", ["15", "12", "12", "12"]),
            ("--alignSJoverhangMin", ["15"]),
            ("--alignSJDBoverhangMin", ["15"]),
            ("--outFilterMultimapNmax", ["20"]),
            ("--outFilterScoreMin", ["1"]),
            ("--outFilterMatchNmin", ["1"]),
            ("--outFilterMismatchNmax", ["2"]),
            ("--chimSegmentMin", ["15"]),
            ("--chimScoreMin", ["15"]),
            ("--chimJunctionOverhangMin", ["15"]),
        ]
        merged: list[str] = []
        for flag, values in required_pairs:
            if not _flag_present(tokens, flag):
                merged.extend([flag, *values])
        merged.extend(tokens)
        return shlex.join(merged)
    if extra_flags.strip():
        return extra_flags
    if detector_hint == "ciri3" and aligner == "bwa-mem":
        return "-k 15 -T 15"
    return extra_flags


def _validate_prepare_config(
    *,
    rows: List[SourceManifestRow],
    outdir: Path,
    aligner: str,
    detector_hint: str | None,
    ref_fa: Path | None,
    chunk_size: int,
    parallel: int,
    command_template: str | None,
    extra_flags: str,
    output_format: str,
) -> list[str]:
    errors: list[str] = []
    effective_output_format = _effective_output_format(
        detector_hint=detector_hint,
        aligner=aligner,
        output_format=output_format,
    )
    if chunk_size < 1:
        errors.append("chunk_size must be >= 1")
    if parallel < 1:
        errors.append("parallel must be >= 1")
    if output_format not in {"bam", "sam"}:
        errors.append("output_format must be 'bam' or 'sam'")
    if outdir.exists() and (outdir / "alignment_manifest.tsv").exists():
        summary = _load_json(outdir / "alignment_prepare_summary.json")
        if summary and summary.get("aligner") not in {None, aligner}:
            errors.append(
                f"Existing output directory {outdir} was initialized with aligner={summary.get('aligner')}, not {aligner}"
            )
    if aligner == "bwa-mem":
        if ref_fa is None:
            errors.append("bwa-mem alignment requires --ref-fa")
        if not _tool_exists("bwa"):
            errors.append("bwa not found on PATH")
        if effective_output_format == "bam" and not _tool_exists("samtools"):
            errors.append("samtools not found on PATH; BAM output requires samtools")
        for row in rows:
            if row.read1 is None:
                errors.append(f"bwa-mem requires FASTQ input; cell_id={row.cell_id} only provides BAM")
    elif aligner == "star":
        if not _tool_exists("STAR"):
            errors.append("STAR not found on PATH")
        if _star_genome_dir(extra_flags) is None:
            errors.append(
                "STAR alignment requires a genome index. Pass it via "
                '--extra-flags "--genomeDir /path/to/star_index" '
                "or set CIRCYTO_STAR_GENOME_DIR."
            )
        for row in rows:
            if row.read1 is None:
                errors.append(f"star requires FASTQ input; cell_id={row.cell_id} only provides BAM")
            if detector_hint == "ciri3" and row.read_layout != "paired-end":
                errors.append(
                    "STAR + CIRI3 is currently supported only for paired-end data; "
                    f"cell_id={row.cell_id} has read_layout={row.read_layout}. "
                    "Single-end STAR + CIRI3 is not implemented in the current hybrid rescue path. "
                    "Use BWA + CIRI3 for the sentinel/production baseline, or use paired-end STAR + CIRI3."
                )
            if any(path and path.suffix == ".gz" for path in (row.read1, row.effective_read2)):
                if _star_read_files_command(extra_flags) is None:
                    errors.append(
                        "STAR input FASTQs appear to be gzipped; configure decompression with "
                        '--extra-flags "--readFilesCommand zcat" '
                        f"(cell_id={row.cell_id})."
                    )
    elif aligner == "reuse-existing":
        for row in rows:
            if row.bam is None:
                errors.append(f"reuse-existing requires BAM input; cell_id={row.cell_id} has no bam column")
    elif command_template is None:
        errors.append(
            "No built-in alignment recipe for this request. Provide --command-template or use --aligner reuse-existing/bwa-mem/star."
        )
    if detector_hint == "ciri3" and aligner == "reuse-existing":
        errors.append("CIRI3 requires aligner-generated unsorted SAM input; reuse-existing BAM cache is not supported.")
    return errors


def _validate_resume_state(outdir: Path) -> list[str]:
    errors: list[str] = []
    manifest_path = outdir / "alignment_manifest.tsv"
    summary_path = outdir / "alignment_prepare_summary.json"
    if manifest_path.exists() and not summary_path.exists():
        errors.append(f"Bad resume state: {manifest_path} exists but {summary_path} is missing")
    return errors


def _run_shell_template(
    *,
    template: str,
    context: dict[str, Any],
    log_path: Path,
) -> None:
    try:
        cmd = template.format(**context)
    except KeyError as exc:
        raise RuntimeError(f"Alignment template references unsupported placeholder: {exc}") from exc
    with log_path.open("w", encoding="utf-8") as log_handle:
        result = subprocess.run(
            cmd,
            shell=True,
            executable="/bin/bash",
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            check=False,
        )
    if result.returncode != 0:
        tail = _tail_text(log_path)
        suffix = f"\n--- log tail ---\n{tail}" if tail else ""
        raise RuntimeError(f"Alignment command failed; see {log_path}.{suffix}")


def _alignment_command_preview(
    *,
    row: SourceManifestRow,
    aligner: str,
    ref_fa: Path | None,
    threads: int,
    extra_flags: str,
    output_format: str,
    command_template: str | None,
    outdir: Path,
) -> str:
    effective_output_format = _effective_output_format(
        detector_hint=None,
        aligner=aligner,
        output_format=output_format,
    )
    cache_key = _alignment_cache_key(
        row,
        aligner=aligner,
        ref_fa=ref_fa,
        detector_hint=None,
        threads=threads,
        extra_flags=extra_flags,
    )
    cell_cache_dir = outdir / "cache" / cache_key
    out_path = _alignment_output_path(cell_cache_dir, effective_output_format)
    if aligner == "reuse-existing":
        return f"reuse-existing {row.bam} -> {out_path}"
    if aligner == "bwa-mem":
        if ref_fa is None:
            return "<bwa-mem requires ref_fa>"
        cmd = _build_bwa_mem_command(
            row=row,
            ref_fa=ref_fa,
            threads=threads,
            extra_flags=extra_flags,
        )
        if effective_output_format == "bam":
            return _shell_join(cmd) + f" | samtools sort -@ {max(1, threads)} -o {shlex.quote(str(out_path))}"
        return _shell_join(cmd) + f" > {shlex.quote(str(out_path))}"
    if aligner == "star":
        prefix = f"{cell_cache_dir / 'star_run'}/"
        cmd = ["STAR", "--runThreadN", str(threads)]
        if extra_flags.strip():
            cmd.extend(shlex.split(extra_flags))
        if _star_genome_dir(extra_flags) is None:
            cmd.extend(["--genomeDir", "<genomeDir-required>"])
        if row.read1 is not None:
            cmd.extend(["--readFilesIn", str(row.read1)])
            if row.effective_read2 is not None:
                cmd.append(str(row.effective_read2))
        cmd.extend(["--outFileNamePrefix", prefix])
        return " ".join(shlex.quote(part) for part in cmd)
    if command_template:
        context = {
            "cell_id": row.cell_id,
            "read1": row.read1 or "",
            "read2": row.effective_read2 or "",
            "bam": row.bam or "",
            "ref_fa": ref_fa or "",
            "out_path": out_path,
            "threads": threads,
            "extra_flags": extra_flags,
            "read_layout": row.read_layout,
            "log_path": cell_cache_dir / f"{row.cell_id}.align.log",
        }
        return command_template.format(**context)
    return "<no preview available>"


def _run_bwa_mem_alignment(
    *,
    row: SourceManifestRow,
    ref_fa: Path,
    out_path: Path,
    threads: int,
    extra_flags: str,
    output_format: str,
    log_path: Path,
) -> None:
    bwa_cmd = _build_bwa_mem_command(
        row=row,
        ref_fa=ref_fa,
        threads=threads,
        extra_flags=extra_flags,
    )

    tmp_out = out_path.with_suffix(out_path.suffix + ".partial")
    if tmp_out.exists():
        tmp_out.unlink()

    with log_path.open("w", encoding="utf-8") as log_handle:
        log_handle.write(f"$ {_shell_join(bwa_cmd)}\n")
        log_handle.write(f"# read_layout={row.read_layout}\n")
        log_handle.write(f"# read1={row.read1}\n")
        log_handle.write(f"# read2={row.effective_read2 or '<none>'}\n")
        log_handle.flush()
        bwa_proc = subprocess.Popen(
            bwa_cmd,
            stdout=subprocess.PIPE,
            stderr=log_handle,
            text=False,
        )
        if output_format == "bam":
            samtools_cmd = ["samtools", "sort", "-@", str(max(1, threads)), "-o", str(tmp_out)]
            sort_proc = subprocess.Popen(
                samtools_cmd,
                stdin=bwa_proc.stdout,
                stdout=log_handle,
                stderr=log_handle,
                text=False,
            )
            assert bwa_proc.stdout is not None
            bwa_proc.stdout.close()
            bwa_rc = bwa_proc.wait()
            sort_rc = sort_proc.wait()
            if bwa_rc != 0 or sort_rc != 0:
                details = _alignment_failure_details(
                    row=row,
                    command=bwa_cmd,
                    log_path=log_path,
                    exit_code=bwa_rc if bwa_rc != 0 else sort_rc,
                )
                raise RuntimeError(
                    _format_alignment_failure(
                        f"bwa-mem pipeline failed for {row.cell_id}",
                        details,
                        extra=f"samtools_sort={sort_rc}",
                    )
                )
        else:
            assert bwa_proc.stdout is not None
            with tmp_out.open("wb") as out_handle:
                shutil.copyfileobj(bwa_proc.stdout, out_handle)
            bwa_proc.stdout.close()
            bwa_rc = bwa_proc.wait()
            if bwa_rc != 0:
                details = _alignment_failure_details(
                    row=row,
                    command=bwa_cmd,
                    log_path=log_path,
                    exit_code=bwa_rc,
                )
                raise RuntimeError(_format_alignment_failure(f"bwa mem failed for {row.cell_id}", details))

    tmp_out.replace(out_path)


def _run_star_alignment(
    *,
    row: SourceManifestRow,
    out_path: Path,
    chimeric_path: Path,
    unmapped_mate1_path: Path | None,
    unmapped_mate2_path: Path | None,
    bwa_sam_path: Path | None,
    ref_fa: Path | None,
    threads: int,
    extra_flags: str,
    log_path: Path,
) -> None:
    genome_dir = _star_genome_dir(extra_flags)
    if genome_dir is None:
        raise RuntimeError("STAR alignment requires --genomeDir in --extra-flags or CIRCYTO_STAR_GENOME_DIR.")

    star_dir = out_path.parent / "star_run"
    star_dir.mkdir(parents=True, exist_ok=True)
    extra_tokens = shlex.split(extra_flags.strip()) if extra_flags.strip() else []
    skip_next = False
    filtered_tokens: list[str] = []
    for token in extra_tokens:
        if skip_next:
            skip_next = False
            continue
        if token == "--genomeDir":
            skip_next = True
            continue
        if token.startswith("--genomeDir="):
            continue
        filtered_tokens.append(token)
    if not _flag_present(filtered_tokens, "--outSAMtype"):
        filtered_tokens.extend(["--outSAMtype", "SAM"])
    if not _flag_present(filtered_tokens, "--chimOutType"):
        filtered_tokens.extend(["--chimOutType", "Junctions"])

    with tempfile.TemporaryDirectory(prefix="circyto_star_", dir=_star_tmp_dir_base()) as tmp_root_str:
        tmp_root = Path(tmp_root_str)
        run_dir = tmp_root / "run"
        run_dir.mkdir(parents=True, exist_ok=True)
        prefix = f"{run_dir}/"
        out_tmp_dir = tmp_root / "_STARtmp"
        cmd = [
            "STAR",
            "--runThreadN",
            str(max(1, threads)),
            "--genomeDir",
            genome_dir,
            "--readFilesIn",
            str(row.read1),
        ]
        if row.effective_read2 is not None:
            cmd.append(str(row.effective_read2))
        cmd.extend(filtered_tokens)
        cmd.extend(["--outFileNamePrefix", prefix, "--outTmpDir", str(out_tmp_dir)])

        with log_path.open("a", encoding="utf-8") as log_handle:
            log_handle.write(
                f"# STAR temp workspace: {tmp_root}"
                " (override with CIRCYTO_STAR_TMPDIR or TMPDIR)\n"
            )
            log_handle.write(
                "# Running STAR in a local Linux temp directory before copying outputs "
                "back into the alignment cache.\n"
            )
            log_handle.flush()
        _emit_prepare_progress(
            cell_id=row.cell_id,
            stage="star-start",
            log_path=log_path,
            detail=f"genome_dir={genome_dir} threads={max(1, threads)}",
        )
        with log_path.open("a", encoding="utf-8") as log_handle:
            result = subprocess.run(
                cmd,
                stdout=log_handle,
                stderr=subprocess.STDOUT,
                check=False,
                text=True,
            )
        _append_log_line(log_path, f"[circyto] cell={row.cell_id} stage=star-end returncode={result.returncode}")

        for name in ("Log.out", "Log.progress.out", "Log.final.out", "SJ.out.tab"):
            candidate = run_dir / name
            if candidate.exists():
                shutil.copy2(candidate, star_dir / name)

        if result.returncode != 0:
            tail = _tail_text(log_path)
            suffix = f"\n--- log tail ---\n{tail}" if tail else ""
            raise RuntimeError(
                f"STAR alignment failed for {row.cell_id}; see {log_path}. "
                "If you are on WSL or a network-mounted workspace, prefer a local Linux temp directory "
                "for STAR via CIRCYTO_STAR_TMPDIR or by leaving TMPDIR unset so /tmp is used."
                f"{suffix}"
            )

        produced_sam = run_dir / "Aligned.out.sam"
        produced_chimeric = run_dir / "Chimeric.out.junction"
        produced_unmapped_mate1 = run_dir / "Unmapped.out.mate1"
        produced_unmapped_mate2 = run_dir / "Unmapped.out.mate2"
        if not produced_sam.exists():
            raise RuntimeError(f"STAR completed but did not produce {produced_sam}")
        if not produced_chimeric.exists():
            raise RuntimeError(f"STAR completed but did not produce {produced_chimeric}")
        _emit_prepare_progress(cell_id=row.cell_id, stage="star-copy-start", log_path=log_path)
        shutil.copy2(produced_sam, out_path)
        shutil.copy2(produced_chimeric, chimeric_path)
        if unmapped_mate1_path is not None:
            if not produced_unmapped_mate1.exists():
                raise RuntimeError(f"STAR completed but did not produce {produced_unmapped_mate1}")
            shutil.copy2(produced_unmapped_mate1, unmapped_mate1_path)
        if unmapped_mate2_path is not None:
            if not produced_unmapped_mate2.exists():
                raise RuntimeError(f"STAR completed but did not produce {produced_unmapped_mate2}")
            shutil.copy2(produced_unmapped_mate2, unmapped_mate2_path)
        _emit_prepare_progress(
            cell_id=row.cell_id,
            stage="star-copy-end",
            log_path=log_path,
            detail=(
                f"sam_bytes={out_path.stat().st_size} "
                f"chimeric_bytes={chimeric_path.stat().st_size} "
                f"mate1_bytes={unmapped_mate1_path.stat().st_size if unmapped_mate1_path else 0} "
                f"mate2_bytes={unmapped_mate2_path.stat().st_size if unmapped_mate2_path else 0}"
            ),
        )

    if bwa_sam_path is not None:
        if ref_fa is None:
            raise RuntimeError("STAR+CIRI3 hybrid rescue requires ref_fa.")
        if unmapped_mate1_path is None or unmapped_mate2_path is None:
            raise RuntimeError("STAR+CIRI3 hybrid rescue requires both unmapped mates.")
        bwa_cmd = [
            "bwa",
            "mem",
            "-T",
            "19",
            str(ref_fa),
            str(unmapped_mate1_path),
            str(unmapped_mate2_path),
        ]
        bwa_started = time.perf_counter()
        _emit_prepare_progress(
            cell_id=row.cell_id,
            stage="bwa-rescue-start",
            log_path=log_path,
            detail=(
                f"mate1_bytes={unmapped_mate1_path.stat().st_size} "
                f"mate2_bytes={unmapped_mate2_path.stat().st_size}"
            ),
        )
        with log_path.open("a", encoding="utf-8") as log_handle:
            log_handle.write(f"\n$ {shlex.join(bwa_cmd)} > {bwa_sam_path}\n")
            log_handle.flush()
            with bwa_sam_path.open("w", encoding="utf-8") as out_handle:
                result = subprocess.run(
                    bwa_cmd,
                    stdout=out_handle,
                    stderr=log_handle,
                    check=False,
                    text=True,
                )
        if result.returncode != 0:
            tail = _tail_text(log_path)
            suffix = f"\n--- log tail ---\n{tail}" if tail else ""
            raise RuntimeError(f"STAR rescue bwa mem failed for {row.cell_id}; see {log_path}.{suffix}")
        if not bwa_sam_path.exists() or bwa_sam_path.stat().st_size == 0:
            raise RuntimeError(f"STAR rescue bwa mem produced empty output for {row.cell_id}: {bwa_sam_path}")
        _emit_prepare_progress(
            cell_id=row.cell_id,
            stage="bwa-rescue-end",
            log_path=log_path,
            detail=f"seconds={round(time.perf_counter() - bwa_started, 3)} output_bytes={bwa_sam_path.stat().st_size}",
        )


def _flatten_detector_result(result: Any) -> List[DetectorResult]:
    if isinstance(result, DetectorResult):
        return [result]
    if isinstance(result, (list, tuple)):
        flat: List[DetectorResult] = []
        for item in result:
            flat.extend(_flatten_detector_result(item))
        return flat
    raise TypeError(f"Detector returned unsupported result type: {type(result)!r}")


def _run_detector(detector: DetectorBase, inputs: DetectorRunInputs) -> DetectorResult | list[DetectorResult]:
    if inputs.input_mode == "alignment" and hasattr(detector, "run_from_alignment"):
        return getattr(detector, "run_from_alignment")(inputs)
    if inputs.input_mode == "fastq" and hasattr(detector, "run_from_fastq"):
        return getattr(detector, "run_from_fastq")(inputs)
    return detector.run(inputs)


def _build_detector_alignment_stamp(
    detector: DetectorBase,
    *,
    row: AlignmentManifestRow,
    ref_fa: Path | None,
    gtf: Path | None,
    threads: int,
) -> dict[str, Any]:
    return {
        "detector": detector.name,
        "cell_id": row.cell_id,
        "alignment_path": str(Path(row.alignment_path).resolve()),
        "aligner": row.aligner,
        "sortedness": _row_field(row, "sortedness"),
        "mapper_mode": _row_field(row, "mapper_mode"),
        "chimeric_junction": row.chimeric_junction,
        "unmapped_mate1": row.unmapped_mate1,
        "unmapped_mate2": row.unmapped_mate2,
        "bwa_sam": row.bwa_sam,
        "cache_key": row.cache_key,
        "source_manifest": row.source_manifest,
        "read_layout": row.read_layout,
        "ref_fa": str(ref_fa.resolve()) if ref_fa else None,
        "gtf": str(gtf.resolve()) if gtf else None,
        "threads": threads,
        "input_mode": "alignment",
        "execution_mode": "alignment-first",
    }


def _provenance_matches(existing: dict[str, Any] | None, expected: dict[str, Any]) -> bool:
    if existing is None:
        return False
    if existing == expected:
        return True
    return all(existing.get(key) == value for key, value in expected.items()) or all(
        expected.get(key) == value for key, value in existing.items() if key in expected
    )


def _load_provenance(path: Path) -> dict[str, Any] | None:
    return _load_json(path)


def _validate_alignment_rows(
    rows: List[AlignmentManifestRow],
    *,
    detector: DetectorBase | None = None,
    ref_fa: Path | None = None,
    gtf: Path | None = None,
) -> list[str]:
    errors: list[str] = []
    for row in rows:
        alignment_path = Path(row.alignment_path)
        if not alignment_path.exists():
            errors.append(f"Missing alignment file for cell_id={row.cell_id}: {alignment_path}")
        if not row.read_layout:
            errors.append(f"Missing read_layout for cell_id={row.cell_id}")
        elif row.read_layout not in VALID_READ_LAYOUTS:
            errors.append(f"Invalid read_layout for cell_id={row.cell_id}: {row.read_layout}")
        if not row.cache_key:
            errors.append(f"Missing cache_key for cell_id={row.cell_id}")
        if not row.aligner:
            errors.append(f"Missing aligner for cell_id={row.cell_id}")
        if detector is not None:
            caps = get_detector_capabilities(detector)
            if not caps.accepts_alignment:
                errors.append(f"Detector '{detector.name}' does not accept alignment inputs")
        if detector.name == "ciri3":
            ciri = detector if isinstance(detector, Ciri3Detector) else Ciri3Detector()
            ok, ciri_errors, _ = ciri.validate_runtime()
            if not ok:
                errors.extend([f"CIRI3 preflight: {msg}" for msg in ciri_errors])
            if row.bam:
                errors.append(f"CIRI3 requires unsorted SAM input; cell_id={row.cell_id} currently points to BAM")
            sortedness = _row_field(row, "sortedness")
            if sortedness and sortedness != "unsorted":
                errors.append(f"CIRI3 requires unsorted alignment input; cell_id={row.cell_id} has sortedness={sortedness}")
            mapper_mode = _row_field(row, "mapper_mode")
            if mapper_mode == "1":
                if row.aligner != "star":
                    errors.append(f"CIRI3 STAR mode requires aligner=star; cell_id={row.cell_id} has aligner={row.aligner or 'unknown'}")
                if not row.sam:
                    errors.append(f"CIRI3 STAR mode requires STAR-generated aligned SAM; cell_id={row.cell_id}")
                if not row.chimeric_junction:
                    errors.append(f"CIRI3 STAR mode requires chimeric_junction; cell_id={row.cell_id}")
                if not row.bwa_sam:
                    errors.append(f"CIRI3 STAR mode requires bwa_sam rescue alignment; cell_id={row.cell_id}")
                if row.sam and not Path(row.sam).exists():
                    errors.append(f"CIRI3 STAR mode aligned SAM missing for cell_id={row.cell_id}: {row.sam}")
                if row.chimeric_junction and not Path(row.chimeric_junction).exists():
                    errors.append(f"CIRI3 STAR mode chimeric_junction missing for cell_id={row.cell_id}: {row.chimeric_junction}")
                if row.bwa_sam and not Path(row.bwa_sam).exists():
                    errors.append(f"CIRI3 STAR mode bwa_sam missing for cell_id={row.cell_id}: {row.bwa_sam}")
        if detector.name in {"ciri2", "ciri-full"} and gtf is None:
            errors.append(f"Detector '{detector.name}' requires --gtf")
    if ref_fa is not None and not ref_fa.exists():
        errors.append(f"Reference FASTA not found: {ref_fa}")
    if gtf is not None and not gtf.exists():
        errors.append(f"GTF not found: {gtf}")
    return errors


def plan_alignment_cache(
    *,
    manifest: Path,
    outdir: Path,
    aligner: str = "reuse-existing",
    ref_fa: Path | None = None,
    detector_hint: str | None = None,
    threads: int = 8,
    parallel: int = 4,
    sentinel_cells: int = 0,
    chunk_size: int = 25,
    command_template: str | None = None,
    extra_flags: str = "",
    output_format: str = "bam",
    preview_rows: int = 3,
) -> dict[str, Any]:
    rows = read_source_manifest(manifest, validate_files=True)
    effective_extra_flags = _effective_alignment_extra_flags(
        detector_hint=detector_hint,
        aligner=aligner,
        extra_flags=extra_flags,
    )
    effective_output_format = _effective_output_format(
        detector_hint=detector_hint,
        aligner=aligner,
        output_format=output_format,
    )
    errors = _validate_prepare_config(
        rows=rows,
        outdir=outdir,
        aligner=aligner,
        detector_hint=detector_hint,
        ref_fa=ref_fa,
        chunk_size=chunk_size,
        parallel=parallel,
        command_template=command_template,
        extra_flags=effective_extra_flags,
        output_format=output_format,
    )
    errors.extend(_validate_resume_state(outdir))
    chunk_count = (len(rows) + max(1, chunk_size) - 1) // max(1, chunk_size)
    n_fastq = sum(1 for row in rows if row.input_mode == "fastq")
    n_bam = sum(1 for row in rows if row.input_mode == "alignment")
    payload = {
        "manifest": str(manifest.resolve()),
        "outdir": str(outdir.resolve()),
        "aligner": aligner,
        "detector_hint": detector_hint,
        "threads": threads,
        "parallel": parallel,
        "sentinel_cells": sentinel_cells,
        "chunk_size": chunk_size,
        "chunk_count": chunk_count,
        "output_format": effective_output_format,
        "extra_flags": effective_extra_flags,
        "n_rows": len(rows),
        "n_fastq_rows": n_fastq,
        "n_bam_rows": n_bam,
        "command_preview": [
            {
                "cell_id": row.cell_id,
                "command": _alignment_command_preview(
                    row=row,
                    aligner=aligner,
                    ref_fa=ref_fa,
                    threads=threads,
                    extra_flags=effective_extra_flags,
                    output_format=effective_output_format,
                    command_template=command_template,
                    outdir=outdir,
                ),
            }
            for row in rows[: max(0, preview_rows)]
        ],
        "errors": errors,
    }
    return payload


def summarize_alignment_chunks(outdir: Path) -> dict[str, Any]:
    chunk_dir = _chunk_dir(outdir)
    summaries = sorted(chunk_dir.glob("chunk_*.json"))
    chunks: list[dict[str, Any]] = []
    status_counts: dict[str, int] = {}
    for path in summaries:
        payload = _load_json(path)
        if payload is None:
            continue
        chunks.append(payload)
        status = str(payload.get("status", "unknown"))
        status_counts[status] = status_counts.get(status, 0) + 1
    failed_chunks = [chunk["chunk_index"] for chunk in chunks if chunk.get("status") in {"failed", "partial_failure"}]
    failed_cells = sorted(
        {
            cell["cell_id"]
            for chunk in chunks
            for cell in chunk.get("cells", [])
            if cell.get("status") == "failed"
        }
    )
    summary = _load_json(outdir / "alignment_prepare_summary.json") or _load_json(outdir / "detector_run_summary.json") or {}
    return {
        "outdir": str(outdir.resolve()),
        "n_chunks": len(chunks),
        "status_counts": status_counts,
        "failed_chunks": failed_chunks,
        "failed_cells": failed_cells,
        "summary": summary,
        "chunks": chunks,
    }


def export_manifest_subset(
    *,
    manifest: Path,
    run_dir: Path,
    out_path: Path,
    subset: str,
    chunk_index: int | None = None,
) -> Path:
    if subset not in VALID_SUBSETS and chunk_index is None:
        raise ValueError(f"Unsupported subset '{subset}'. Supported: {', '.join(sorted(VALID_SUBSETS))}")
    header, rows = _read_manifest_rows_raw(manifest)
    source_cell_ids = [str(row.get("cell_id", "")).strip() for row in rows]
    summary = _load_json(run_dir / "alignment_prepare_summary.json") or _load_json(run_dir / "detector_run_summary.json") or {}
    chunk_summary = summarize_alignment_chunks(run_dir)
    status_by_cell = {str(cell.get("cell_id")): str(cell.get("status")) for cell in summary.get("cells", [])}
    run_state = summarize_run_state(
        manifest=manifest,
        run_dir=run_dir,
        mode="prepare" if (run_dir / "alignment_prepare_summary.json").exists() else "detector",
    )

    if chunk_index is not None:
        selected = {
            str(cell.get("cell_id"))
            for chunk in chunk_summary["chunks"]
            if int(chunk.get("chunk_index", -1)) == chunk_index
            for cell in chunk.get("cells", [])
        }
    elif subset == "failed":
        selected = {cell_id for cell_id, status in status_by_cell.items() if status == "failed"}
    elif subset == "missing":
        selected = {cell_id for cell_id in source_cell_ids if cell_id and cell_id not in status_by_cell}
    elif subset == "stale":
        selected = set(run_state["stale_cells"])
    elif subset == "incomplete":
        selected = set(run_state["incomplete_cells"])
    else:
        failed_chunks = set(chunk_summary["failed_chunks"])
        selected = {
            str(cell.get("cell_id"))
            for chunk in chunk_summary["chunks"]
            if int(chunk.get("chunk_index", -1)) in failed_chunks
            for cell in chunk.get("cells", [])
        }

    filtered = [row for row in rows if str(row.get("cell_id", "")).strip() in selected]
    _write_manifest_rows_raw(out_path, header=header, rows=filtered)
    return out_path


def summarize_run_state(
    *,
    manifest: Path,
    run_dir: Path,
    mode: str,
) -> dict[str, Any]:
    if mode not in {"prepare", "detector"}:
        raise ValueError("mode must be 'prepare' or 'detector'")
    header, rows = _read_manifest_rows_raw(manifest)
    cell_ids = [str(row.get("cell_id", "")).strip() for row in rows if str(row.get("cell_id", "")).strip()]
    summary = _load_json(run_dir / "alignment_prepare_summary.json") or _load_json(run_dir / "detector_run_summary.json") or {}
    status_by_cell = {str(cell.get("cell_id")): str(cell.get("status")) for cell in summary.get("cells", [])}
    stale_cells: list[str] = []

    if mode == "prepare":
        alignment_manifest = run_dir / "alignment_manifest.tsv"
        if alignment_manifest.exists():
            for row in read_alignment_manifest_tsv(alignment_manifest, validate_files=False):
                alignment_path = Path(row.alignment_path)
                if not alignment_path.exists() or _load_provenance(_alignment_provenance_path(alignment_path)) is None:
                    stale_cells.append(row.cell_id)
    else:
        for cell_id in cell_ids:
            output_path = run_dir / f"{cell_id}.tsv"
            if output_path.exists() and _load_provenance(_detector_output_provenance_path(output_path)) is None:
                stale_cells.append(cell_id)
            if status_by_cell.get(cell_id) in {"success", "skipped_existing"} and not output_path.exists():
                stale_cells.append(cell_id)

    completed = sorted(cell_id for cell_id, status in status_by_cell.items() if status not in {"failed"})
    failed = sorted(cell_id for cell_id, status in status_by_cell.items() if status == "failed")
    missing = sorted(cell_id for cell_id in cell_ids if cell_id not in status_by_cell)
    incomplete = sorted(set(failed) | set(missing) | set(stale_cells))

    return {
        "manifest": str(manifest.resolve()),
        "run_dir": str(run_dir.resolve()),
        "mode": mode,
        "planned_cells": len(cell_ids),
        "completed_cells": len(completed),
        "failed_cells": len(failed),
        "missing_cells": len(missing),
        "stale_cells": sorted(set(stale_cells)),
        "incomplete_cells": incomplete,
        "failed_cell_ids": failed,
        "missing_cell_ids": missing,
        "summary": summary,
    }


def prepare_alignment_cache(
    *,
    manifest: Path,
    outdir: Path,
    aligner: str = "reuse-existing",
    ref_fa: Path | None = None,
    detector_hint: str | None = None,
    threads: int = 8,
    parallel: int = 4,
    sentinel_cells: int = 0,
    chunk_size: int = 25,
    command_template: str | None = None,
    extra_flags: str = "",
    link_mode: str = "symlink",
    index_bam: bool = False,
    output_format: str = "bam",
    dry_run: bool = False,
    fail_fast: bool = False,
) -> Path:
    plan = plan_alignment_cache(
        manifest=manifest,
        outdir=outdir,
        aligner=aligner,
        ref_fa=ref_fa,
        detector_hint=detector_hint,
        threads=threads,
        parallel=parallel,
        sentinel_cells=sentinel_cells,
        chunk_size=chunk_size,
        command_template=command_template,
        extra_flags=extra_flags,
        output_format=output_format,
    )
    if dry_run:
        _write_json(_summary_path(outdir, "alignment_prepare_plan.json"), plan)
        return outdir / "alignment_manifest.tsv"
    if plan["errors"]:
        raise RuntimeError("Alignment preflight failed: " + "; ".join(plan["errors"]))

    rows = read_source_manifest(manifest, validate_files=True)
    effective_extra_flags = _effective_alignment_extra_flags(
        detector_hint=detector_hint,
        aligner=aligner,
        extra_flags=extra_flags,
    )
    effective_output_format = _effective_output_format(
        detector_hint=detector_hint,
        aligner=aligner,
        output_format=output_format,
    )
    outdir.mkdir(parents=True, exist_ok=True)
    cache_dir = outdir / "cache"
    staged_dir = outdir / _staged_alignment_bucket(aligner=aligner, output_format=effective_output_format)
    summary_path = _summary_path(outdir, "alignment_prepare_summary.json")
    manifest_path = outdir / "alignment_manifest.tsv"

    ordered_rows = rows[:]
    if sentinel_cells > 0:
        ordered_rows = rows[:sentinel_cells] + rows[sentinel_cells:]

    effective_parallel = max(1, min(parallel, len(ordered_rows)))
    existing_completed = {
        row.cell_id: row for row in read_alignment_manifest_tsv(manifest_path, validate_files=False)
    } if manifest_path.exists() else {}
    completed_rows: dict[str, AlignmentManifestRow] = dict(existing_completed)
    per_cell_records: List[dict[str, Any]] = []
    failures: List[dict[str, Any]] = []
    started_at = time.perf_counter()

    def _prepare_one(row: SourceManifestRow) -> tuple[AlignmentManifestRow, dict[str, Any]]:
        cache_key = _alignment_cache_key(
            row,
            aligner=aligner,
            ref_fa=ref_fa,
            detector_hint=detector_hint,
            threads=threads,
            extra_flags=effective_extra_flags,
        )
        cell_cache_dir = cache_dir / cache_key
        out_path = _alignment_output_path(cell_cache_dir, effective_output_format)
        staged_path = staged_dir / f"{row.cell_id}{out_path.suffix}"
        chimeric_path = cell_cache_dir / "Chimeric.out.junction"
        staged_chimeric_path = staged_dir / f"{row.cell_id}.Chimeric.out.junction"
        unmapped_mate1_path = cell_cache_dir / "Unmapped.out.mate1" if aligner == "star" and detector_hint == "ciri3" else None
        unmapped_mate2_path = cell_cache_dir / "Unmapped.out.mate2" if aligner == "star" and detector_hint == "ciri3" else None
        bwa_sam_path = cell_cache_dir / "bwa_rescue.sam" if aligner == "star" and detector_hint == "ciri3" else None
        staged_unmapped_mate1_path = staged_dir / f"{row.cell_id}.Unmapped.out.mate1" if unmapped_mate1_path else None
        staged_unmapped_mate2_path = staged_dir / f"{row.cell_id}.Unmapped.out.mate2" if unmapped_mate2_path else None
        staged_bwa_sam_path = staged_dir / f"{row.cell_id}.bwa_rescue.sam" if bwa_sam_path else None
        log_path = cell_cache_dir / f"{row.cell_id}.align.log"
        provenance = {
            "cell_id": row.cell_id,
            "aligner": aligner,
            "detector_hint": detector_hint,
            "cache_key": cache_key,
            "manifest": str(manifest.resolve()),
            "read_layout": row.read_layout,
            "read1": str(row.read1.resolve()) if row.read1 else None,
            "read2": str(row.effective_read2.resolve()) if row.effective_read2 else None,
            "input_bam": str(row.bam.resolve()) if row.bam else None,
            "ref_fa": str(ref_fa.resolve()) if ref_fa else None,
            "threads": threads,
            "extra_flags": effective_extra_flags,
            "output_format": effective_output_format,
        }

        existing_prov = _load_provenance(_alignment_provenance_path(out_path))
        if out_path.exists() and _provenance_matches(existing_prov, provenance):
            _copy_or_link(out_path, staged_path, link_mode=link_mode)
            sortedness = str(existing_prov.get("sortedness", "unknown"))
            mapper_mode = str(existing_prov.get("mapper_mode", ""))
            artifact_bucket = _staged_alignment_bucket(aligner=aligner, output_format=effective_output_format)
            row_kwargs = {
                "chimeric_junction": "",
                "unmapped_mate1": "",
                "unmapped_mate2": "",
                "bwa_sam": "",
            }
            extra = {"reused_alignment": "true"}
            if aligner == "star" and chimeric_path.exists():
                _copy_or_link(chimeric_path, staged_chimeric_path, link_mode=link_mode)
                row_kwargs["chimeric_junction"] = str(staged_chimeric_path.resolve())
            if staged_unmapped_mate1_path is not None and unmapped_mate1_path is not None and unmapped_mate1_path.exists():
                _copy_or_link(unmapped_mate1_path, staged_unmapped_mate1_path, link_mode=link_mode)
                row_kwargs["unmapped_mate1"] = str(staged_unmapped_mate1_path.resolve())
            if staged_unmapped_mate2_path is not None and unmapped_mate2_path is not None and unmapped_mate2_path.exists():
                _copy_or_link(unmapped_mate2_path, staged_unmapped_mate2_path, link_mode=link_mode)
                row_kwargs["unmapped_mate2"] = str(staged_unmapped_mate2_path.resolve())
            if staged_bwa_sam_path is not None and bwa_sam_path is not None and bwa_sam_path.exists():
                _copy_or_link(bwa_sam_path, staged_bwa_sam_path, link_mode=link_mode)
                row_kwargs["bwa_sam"] = str(staged_bwa_sam_path.resolve())
            return (
                AlignmentManifestRow(
                    cell_id=row.cell_id,
                    bam=str(staged_path.resolve()) if staged_path.suffix == ".bam" else "",
                    sam=str(staged_path.resolve()) if staged_path.suffix == ".sam" else "",
                    group_id=row.group_id,
                    read_layout=row.read_layout,
                    aligner=aligner,
                    reference=str(ref_fa.resolve()) if ref_fa else "",
                    cache_key=cache_key,
                    source_manifest=str(manifest.resolve()),
                    mapper_mode=mapper_mode,
                    artifact_bucket=artifact_bucket,
                    sortedness=sortedness,
                    **row_kwargs,
                    extra=extra,
                ),
                {
                    "cell_id": row.cell_id,
                    "status": "reused_cached",
                    "reused_alignment": True,
                    "cache_key": cache_key,
                    "alignment_path": str(staged_path.resolve()),
                    "artifact_bucket": artifact_bucket,
                    "sortedness": sortedness,
                    "mapper_mode": mapper_mode,
                    "read_layout": row.read_layout,
                },
            )

        cell_cache_dir.mkdir(parents=True, exist_ok=True)
        cell_started = time.perf_counter()
        reused_alignment = False
        if log_path.exists():
            log_path.unlink()
        _emit_prepare_progress(
            cell_id=row.cell_id,
            stage="prepare-start",
            log_path=log_path,
            detail=f"aligner={aligner} read_layout={row.read_layout}",
        )
        if aligner == "reuse-existing":
            assert row.bam is not None
            source = row.bam
            _copy_or_link(source, out_path, link_mode="copy")
            reused_alignment = True
        elif aligner == "bwa-mem":
            assert ref_fa is not None
            _run_bwa_mem_alignment(
                row=row,
                ref_fa=ref_fa,
                out_path=out_path,
                threads=threads,
                extra_flags=effective_extra_flags,
                output_format=effective_output_format,
                log_path=log_path,
            )
        elif aligner == "star":
            _run_star_alignment(
                row=row,
                out_path=out_path,
                chimeric_path=chimeric_path,
                unmapped_mate1_path=unmapped_mate1_path,
                unmapped_mate2_path=unmapped_mate2_path,
                bwa_sam_path=bwa_sam_path,
                ref_fa=ref_fa,
                threads=threads,
                extra_flags=effective_extra_flags,
                log_path=log_path,
            )
        else:
            context = {
                "cell_id": row.cell_id,
                "read1": row.read1 or "",
                "read2": row.effective_read2 or "",
                "bam": row.bam or "",
                "ref_fa": ref_fa or "",
                "out_path": out_path,
                "threads": threads,
                "extra_flags": effective_extra_flags,
                "read_layout": row.read_layout,
                "log_path": log_path,
            }
            assert command_template is not None
            _run_shell_template(template=command_template, context=context, log_path=log_path)
            if not out_path.exists():
                raise RuntimeError(f"Alignment command completed but did not produce expected output: {out_path}")

        if index_bam and out_path.suffix == ".bam":
            subprocess.run(["samtools", "index", str(out_path)], check=False)

        provenance.update(
            {
                "reused_alignment": reused_alignment,
                "elapsed_seconds": round(time.perf_counter() - cell_started, 3),
                "alignment_path": str(out_path.resolve()),
                "sortedness": "unsorted" if effective_output_format == "sam" else "sorted",
                "mapper_mode": "1" if aligner == "star" else "0",
            }
        )
        if aligner == "star":
            provenance["chimeric_junction"] = str(chimeric_path.resolve())
        if unmapped_mate1_path is not None:
            provenance["unmapped_mate1"] = str(unmapped_mate1_path.resolve())
        if unmapped_mate2_path is not None:
            provenance["unmapped_mate2"] = str(unmapped_mate2_path.resolve())
        if bwa_sam_path is not None:
            provenance["bwa_sam"] = str(bwa_sam_path.resolve())
        _write_json(_alignment_provenance_path(out_path), provenance)
        _emit_prepare_progress(cell_id=row.cell_id, stage="stage-artifacts-start", log_path=log_path)
        _copy_or_link(out_path, staged_path, link_mode=link_mode)
        sortedness = "unsorted" if effective_output_format == "sam" else "sorted"
        mapper_mode = "1" if aligner == "star" else "0"
        artifact_bucket = _staged_alignment_bucket(aligner=aligner, output_format=effective_output_format)
        row_kwargs = {
            "chimeric_junction": "",
            "unmapped_mate1": "",
            "unmapped_mate2": "",
            "bwa_sam": "",
        }
        extra = {"reused_alignment": str(reused_alignment).lower()}
        if aligner == "star":
            _copy_or_link(chimeric_path, staged_chimeric_path, link_mode=link_mode)
            row_kwargs["chimeric_junction"] = str(staged_chimeric_path.resolve())
        if staged_unmapped_mate1_path is not None and unmapped_mate1_path is not None:
            _copy_or_link(unmapped_mate1_path, staged_unmapped_mate1_path, link_mode=link_mode)
            row_kwargs["unmapped_mate1"] = str(staged_unmapped_mate1_path.resolve())
        if staged_unmapped_mate2_path is not None and unmapped_mate2_path is not None:
            _copy_or_link(unmapped_mate2_path, staged_unmapped_mate2_path, link_mode=link_mode)
            row_kwargs["unmapped_mate2"] = str(staged_unmapped_mate2_path.resolve())
        if staged_bwa_sam_path is not None and bwa_sam_path is not None:
            _copy_or_link(bwa_sam_path, staged_bwa_sam_path, link_mode=link_mode)
            row_kwargs["bwa_sam"] = str(staged_bwa_sam_path.resolve())
        _emit_prepare_progress(
            cell_id=row.cell_id,
            stage="stage-artifacts-end",
            log_path=log_path,
            detail=f"seconds={round(time.perf_counter() - cell_started, 3)}",
        )
        return (
            AlignmentManifestRow(
                cell_id=row.cell_id,
                bam=str(staged_path.resolve()) if staged_path.suffix == ".bam" else "",
                sam=str(staged_path.resolve()) if staged_path.suffix == ".sam" else "",
                group_id=row.group_id,
                read_layout=row.read_layout,
                aligner=aligner,
                reference=str(ref_fa.resolve()) if ref_fa else "",
                cache_key=cache_key,
                source_manifest=str(manifest.resolve()),
                mapper_mode=mapper_mode,
                artifact_bucket=artifact_bucket,
                sortedness=sortedness,
                **row_kwargs,
                extra=extra,
            ),
            {
                "cell_id": row.cell_id,
                "status": "reused_input" if reused_alignment else "aligned",
                "reused_alignment": reused_alignment,
                "cache_key": cache_key,
                "alignment_path": str(staged_path.resolve()),
                "artifact_bucket": artifact_bucket,
                "sortedness": sortedness,
                "mapper_mode": mapper_mode,
                "read_layout": row.read_layout,
                "seconds": round(time.perf_counter() - cell_started, 3),
            },
        )

    chunk_counter = 0
    for chunk_start in range(0, len(ordered_rows), max(1, chunk_size)):
        chunk_counter += 1
        chunk = ordered_rows[chunk_start : chunk_start + max(1, chunk_size)]
        chunk_started = time.perf_counter()
        chunk_records: list[dict[str, Any]] = []
        chunk_failures: list[dict[str, Any]] = []
        with ThreadPoolExecutor(max_workers=min(effective_parallel, len(chunk))) as executor:
            future_map = {executor.submit(_prepare_one, row): row for row in chunk}
            for future in as_completed(future_map):
                row = future_map[future]
                try:
                    manifest_row, record = future.result()
                    completed_rows[manifest_row.cell_id] = manifest_row
                    chunk_records.append(record)
                    per_cell_records.append(record)
                except Exception as exc:
                    cache_key = _alignment_cache_key(
                        row,
                        aligner=aligner,
                        ref_fa=ref_fa,
                        detector_hint=detector_hint,
                        threads=threads,
                        extra_flags=effective_extra_flags,
                    )
                    failure_log_path = cache_dir / cache_key / f"{row.cell_id}.align.log"
                    failure = {
                        "cell_id": row.cell_id,
                        "error": str(exc),
                        "read_layout": row.read_layout,
                        "read1": str(row.read1.resolve()) if row.read1 else None,
                        "read2": str(row.effective_read2.resolve()) if row.effective_read2 else None,
                        "log_path": str(failure_log_path.resolve()),
                        "stderr_tail": _tail_text(failure_log_path),
                        "command": _alignment_command_preview(
                            row=row,
                            aligner=aligner,
                            ref_fa=ref_fa,
                            threads=threads,
                            extra_flags=effective_extra_flags,
                            output_format=effective_output_format,
                            command_template=command_template,
                            outdir=outdir,
                        ),
                    }
                    chunk_failures.append(failure)
                    failures.append(failure)
                    record = {
                        "cell_id": row.cell_id,
                        "status": "failed",
                        "reused_alignment": False,
                        "error": str(exc),
                        "read_layout": row.read_layout,
                        "read1": failure["read1"],
                        "read2": failure["read2"],
                        "log_path": failure["log_path"],
                        "stderr_tail": failure["stderr_tail"],
                        "command": failure["command"],
                    }
                    chunk_records.append(record)
                    per_cell_records.append(record)

        if completed_rows:
            print(f"[circyto] stage=manifest-write-start rows={len(completed_rows)} path={manifest_path}", flush=True)
            write_alignment_manifest_tsv(completed_rows.values(), manifest_path)
            print(f"[circyto] stage=manifest-write-end rows={len(completed_rows)} path={manifest_path}", flush=True)

        chunk_status = "success" if not chunk_failures else ("partial_failure" if len(chunk_failures) < len(chunk) else "failed")
        _write_json(
            _chunk_summary_path(outdir, chunk_counter),
            {
                "chunk_index": chunk_counter,
                "chunk_start": chunk_start,
                "chunk_size": len(chunk),
                "status": chunk_status,
                "elapsed_seconds": round(time.perf_counter() - chunk_started, 3),
                "cells": sorted(chunk_records, key=lambda item: item["cell_id"]),
                "failures": chunk_failures,
            },
        )
        print(
            f"[circyto] align chunk={chunk_counter} cells={len(chunk)} status={chunk_status} "
            f"failures={len(chunk_failures)} manifest={manifest_path}"
        )
        if fail_fast and chunk_failures:
            break

    per_cell_records.sort(key=lambda item: item["cell_id"])
    status_counts: dict[str, int] = {}
    for record in per_cell_records:
        status = str(record["status"])
        status_counts[status] = status_counts.get(status, 0) + 1

    payload = {
        "manifest": str(manifest.resolve()),
        "outdir": str(outdir.resolve()),
        "aligner": aligner,
        "detector_hint": detector_hint,
        "threads": threads,
        "parallel_requested": parallel,
        "parallel_effective": effective_parallel,
        "chunk_size": chunk_size,
        "sentinel_cells": sentinel_cells,
        "n_manifest_rows": len(rows),
        "status_counts": status_counts,
        "elapsed_seconds": round(time.perf_counter() - started_at, 3),
        "cells": per_cell_records,
        "n_chunks": chunk_counter,
        "output_format": effective_output_format,
        "planned_cells": len(rows),
        "completed_cells": status_counts.get("aligned", 0) + status_counts.get("reused_input", 0) + status_counts.get("reused_cached", 0),
        "failed_cells": status_counts.get("failed", 0),
        "failed_cell_ids": sorted(item["cell_id"] for item in failures),
        "chunk_status_counts": summarize_alignment_chunks(outdir)["status_counts"],
        "failed_chunk_indices": summarize_alignment_chunks(outdir)["failed_chunks"],
        "command_template": command_template,
        "continue_on_error": not fail_fast,
        "read_layout_counts": {
            layout: sum(1 for row in rows if row.read_layout == layout) for layout in sorted(VALID_READ_LAYOUTS)
        },
    }
    if failures:
        payload["failures"] = failures
    _write_json(summary_path, payload)

    if failures:
        failed_cells = ", ".join(item["cell_id"] for item in failures[:5])
        suffix = "" if len(failures) <= 5 else f" (+{len(failures) - 5} more)"
        raise RuntimeError(
            f"Alignment preparation completed with failures for {len(failures)}/{len(rows)} cells: "
            f"{failed_cells}{suffix}. Summary: {summary_path}"
        )

    return manifest_path


def run_detector_alignment_manifest(
    *,
    detector: DetectorBase,
    manifest: Path,
    outdir: Path,
    ref_fa: Path | None = None,
    gtf: Path | None = None,
    threads: int = 8,
    parallel: int = 4,
    chunk_size: int = 50,
    sentinel_cells: int = 0,
    dry_run: bool = False,
    fail_fast: bool = False,
) -> List[DetectorResult]:
    rows = read_alignment_manifest_tsv(manifest, validate_files=not dry_run)
    errors = [] if dry_run else _validate_alignment_rows(rows, detector=detector, ref_fa=ref_fa, gtf=gtf)
    if chunk_size < 1:
        errors.append("chunk_size must be >= 1")
    if errors:
        raise RuntimeError("Alignment-detector preflight failed: " + "; ".join(errors))
    outdir.mkdir(parents=True, exist_ok=True)
    if dry_run:
        preview: list[dict[str, Any]] = []
        if isinstance(detector, Ciri3Detector):
            for row in rows[:3]:
                inputs = DetectorRunInputs(
                    cell_id=row.cell_id,
                    bam=Path(row.alignment_path) if row.bam else None,
                    sam=Path(row.alignment_path) if row.sam else None,
                    outdir=outdir,
                    ref_fa=ref_fa,
                    gtf=gtf,
                    threads=threads,
                    input_mode="alignment",
                    read_layout=row.read_layout or None,
                    alignment_group=row.group_id or None,
                    extra={"alignment_manifest_row": row.to_dict()},
                )
                preview.append(
                    {
                        "cell_id": row.cell_id,
                        "command": detector.preview_command(
                            inputs,
                            raw_output=(outdir / f"{row.cell_id}.ciri3_run" / "ciri3_raw.tsv"),
                            run_dir=(outdir / f"{row.cell_id}.ciri3_run"),
                            log_path=(outdir / f"{row.cell_id}.ciri3.log"),
                        ),
                    }
                )
        _write_json(
            outdir / "detector_alignment_plan.json",
            {
                "detector": detector.name,
                "manifest": str(manifest.resolve()),
                "outdir": str(outdir.resolve()),
                "threads": threads,
                "parallel": parallel,
                "chunk_size": chunk_size,
                "sentinel_cells": sentinel_cells,
                "n_manifest_rows": len(rows),
                "command_preview": preview,
            },
        )
        return []

    ordered_rows = rows[:]
    if sentinel_cells > 0:
        ordered_rows = rows[:sentinel_cells] + rows[sentinel_cells:]

    started_at = time.perf_counter()
    caps = get_detector_capabilities(detector)
    effective_parallel = max(1, min(parallel, caps.max_parallel, len(rows)))
    summary_path = _summary_path(outdir, "detector_run_summary.json")
    existing_summary = _load_json(summary_path) or {}
    chunk_index_base = _next_chunk_index(outdir)
    results: List[DetectorResult] = []
    per_cell_records: List[dict[str, Any]] = []
    failures: List[dict[str, Any]] = []

    def _run_one(row: AlignmentManifestRow) -> tuple[List[DetectorResult], dict[str, Any]]:
        cell_started = time.perf_counter()
        alignment_path = Path(row.alignment_path)
        row_extra = row.extra or {}
        row_sortedness = _row_field(row, "sortedness")
        row_mapper_mode = _row_field(row, "mapper_mode")
        expected_output = outdir / f"{row.cell_id}.tsv"
        expected_stamp = _build_detector_alignment_stamp(
            detector,
            row=row,
            ref_fa=ref_fa,
            gtf=gtf,
            threads=threads,
        )
        if expected_output.exists() and _provenance_matches(
            _load_provenance(_detector_output_provenance_path(expected_output)), expected_stamp
        ):
            return (
                [DetectorResult(detector=detector.name, cell_id=row.cell_id, outdir=outdir, tsv_path=expected_output, meta={"skipped_existing": True})],
                {
                    "cell_id": row.cell_id,
                    "status": "skipped_existing",
                    "outcome_category": "skipped-existing",
                    "seconds": 0.0,
                    "tsv_path": str(expected_output),
                    "read_layout": row.read_layout,
                    "execution_mode": "alignment-first",
                    "input_mode": "alignment",
                    "reused_alignment": row_extra.get("reused_alignment") == "true",
                    "detector_backend": detector.name,
                    "alignment_group": row.group_id,
                    "reference_used": str(ref_fa) if ref_fa else None,
                    "input_file_type": "bam" if row.bam else "sam",
                    "input_sortedness": row_sortedness,
                    "mapper_mode": row_mapper_mode,
                    "raw_output_path": None,
                    "raw_row_count": None,
                    "normalized_row_count": _count_alignment_output_rows(expected_output, detector.name),
                },
            )

        inputs = DetectorRunInputs(
            cell_id=row.cell_id,
            bam=alignment_path if row.bam else None,
            sam=alignment_path if row.sam else None,
            outdir=outdir,
            ref_fa=ref_fa,
            gtf=gtf,
            threads=threads,
            input_mode="alignment",
            read_layout=row.read_layout or None,
            alignment_group=row.group_id or None,
            provenance={
                "alignment_manifest": str(manifest.resolve()),
                "cache_key": row.cache_key,
                "aligner": row.aligner,
                "source_manifest": row.source_manifest,
                "reused_alignment": row_extra.get("reused_alignment") == "true",
            },
            extra={"alignment_manifest_row": row.to_dict()},
        )
        raw_result = _run_detector(detector, inputs)
        flat_results = _flatten_detector_result(raw_result)
        primary_result = flat_results[0]
        primary_meta = primary_result.meta or {}
        primary_path = primary_result.tsv_path
        _write_json(_detector_output_provenance_path(primary_path), expected_stamp)
        normalized_rows = _count_alignment_output_rows(primary_path, detector.name)
        raw_output_path = primary_meta.get("raw_output_path")
        raw_rows = _count_alignment_output_rows(Path(raw_output_path), detector.name) if raw_output_path else None
        status = "success" if (normalized_rows or 0) > 0 else "empty"
        return flat_results, {
            "cell_id": row.cell_id,
            "status": status,
            "outcome_category": _alignment_outcome_category(status, raw_rows=raw_rows, normalized_rows=normalized_rows),
            "seconds": round(time.perf_counter() - cell_started, 3),
            "tsv_path": str(primary_path),
            "read_layout": row.read_layout,
            "execution_mode": "alignment-first",
            "input_mode": "alignment",
            "reused_alignment": row_extra.get("reused_alignment") == "true",
            "detector_backend": detector.name,
            "alignment_group": row.group_id,
            "reference_used": str(ref_fa) if ref_fa else None,
            "input_file_type": primary_meta.get("input_file_type") or ("bam" if row.bam else "sam"),
            "input_sortedness": primary_meta.get("input_sortedness") or row_sortedness,
            "mapper_mode": primary_meta.get("mapper_mode") or row_mapper_mode,
            "raw_output_path": raw_output_path,
            "raw_row_count": raw_rows,
            "normalized_row_count": normalized_rows,
        }

    chunk_counter = 0
    for chunk_start in range(0, len(ordered_rows), max(1, chunk_size)):
        chunk_counter += 1
        chunk_index = chunk_index_base + chunk_counter - 1
        chunk = ordered_rows[chunk_start : chunk_start + max(1, chunk_size)]
        chunk_started = time.perf_counter()
        chunk_records: list[dict[str, Any]] = []
        chunk_failures: list[dict[str, Any]] = []
        with ThreadPoolExecutor(max_workers=min(effective_parallel, len(chunk))) as executor:
            future_map = {executor.submit(_run_one, row): row for row in chunk}
            for future in as_completed(future_map):
                row = future_map[future]
                try:
                    cell_results, record = future.result()
                    results.extend(cell_results)
                    chunk_records.append(record)
                    per_cell_records.append(record)
                except Exception as exc:
                    failure = {"cell_id": row.cell_id, "error": str(exc), "read_layout": row.read_layout}
                    chunk_failures.append(failure)
                    failures.append(failure)
                    record = {
                        "cell_id": row.cell_id,
                        "status": "failed",
                        "outcome_category": "failed",
                        "seconds": None,
                        "error": str(exc),
                        "read_layout": row.read_layout,
                        "execution_mode": "alignment-first",
                        "input_mode": "alignment",
                        "reused_alignment": (row.extra or {}).get("reused_alignment") == "true",
                        "detector_backend": detector.name,
                        "alignment_group": row.group_id,
                        "reference_used": str(ref_fa) if ref_fa else None,
                        "input_file_type": "bam" if row.bam else "sam",
                        "input_sortedness": _row_field(row, "sortedness"),
                        "mapper_mode": _row_field(row, "mapper_mode"),
                        "raw_output_path": None,
                        "raw_row_count": None,
                        "normalized_row_count": None,
                    }
                    chunk_records.append(record)
                    per_cell_records.append(record)

        chunk_status = "success" if not chunk_failures else ("partial_failure" if len(chunk_failures) < len(chunk) else "failed")
        _write_json(
            _chunk_summary_path(outdir, chunk_index),
            {
                "chunk_index": chunk_index,
                "chunk_start": chunk_start,
                "chunk_size": len(chunk),
                "status": chunk_status,
                "elapsed_seconds": round(time.perf_counter() - chunk_started, 3),
                "cells": sorted(chunk_records, key=lambda item: item["cell_id"]),
                "failures": chunk_failures,
            },
        )
        print(
            f"[circyto] detector={detector.name} chunk={chunk_index} cells={len(chunk)} "
            f"status={chunk_status} failures={len(chunk_failures)} outdir={outdir}"
        )
        if fail_fast and chunk_failures:
            break

    per_cell_records = _merge_cell_records(existing_summary, per_cell_records)
    status_counts: dict[str, int] = {}
    for record in per_cell_records:
        status = str(record["status"])
        status_counts[status] = status_counts.get(status, 0) + 1

    chunk_snapshot = summarize_alignment_chunks(outdir)
    preserve_existing_scope = bool(existing_summary) and len(rows) < int(existing_summary.get("planned_cells", 0) or 0)
    planned_cells = max(len(per_cell_records), len(rows), int(existing_summary.get("planned_cells", 0) or 0))
    payload = {
        "detector": detector.name,
        "manifest": str(manifest.resolve()) if not preserve_existing_scope else str(existing_summary.get("manifest", manifest.resolve())),
        "outdir": str(outdir.resolve()),
        "threads": threads,
        "parallel_requested": parallel,
        "parallel_effective": effective_parallel,
        "n_manifest_rows": planned_cells if preserve_existing_scope else len(rows),
        "status_counts": status_counts,
        "elapsed_seconds": round(time.perf_counter() - started_at, 3),
        "input_mode": "alignment",
        "execution_mode": "alignment-first",
        "chunk_size": chunk_size,
        "sentinel_cells": sentinel_cells,
        "n_chunks": chunk_snapshot["n_chunks"],
        "cells": per_cell_records,
        "planned_cells": planned_cells,
        "completed_cells": status_counts.get("success", 0) + status_counts.get("skipped_existing", 0),
        "failed_cells": status_counts.get("failed", 0),
        "failed_cell_ids": sorted(record["cell_id"] for record in per_cell_records if str(record.get("status")) == "failed"),
        "failed_chunk_indices": chunk_snapshot["failed_chunks"],
        "chunk_status_counts": chunk_snapshot["status_counts"],
        "read_layout_counts": (
            existing_summary.get("read_layout_counts")
            if preserve_existing_scope
            else {layout: sum(1 for row in rows if row.read_layout == layout) for layout in sorted(VALID_READ_LAYOUTS)}
        ),
        "aligners": existing_summary.get("aligners") if preserve_existing_scope else sorted({row.aligner for row in rows}),
        "references": existing_summary.get("references") if preserve_existing_scope else sorted({row.reference for row in rows if row.reference}),
        "detector_backend": detector.name,
        "continue_on_error": not fail_fast,
    }
    if isinstance(detector, Ciri3Detector):
        payload["command_template"] = detector.resolve_command_template()
    if failures:
        payload["failures"] = failures
    _write_json(summary_path, payload)

    if failures:
        failed_cells = ", ".join(item["cell_id"] for item in failures[:5])
        suffix = "" if len(failures) <= 5 else f" (+{len(failures) - 5} more)"
        raise RuntimeError(
            f"Detector '{detector.name}' completed with failures for {len(failures)}/{len(rows)} alignment rows: "
            f"{failed_cells}{suffix}. Summary: {summary_path}"
        )

    return results
