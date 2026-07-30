from __future__ import annotations

import csv
import hashlib
import json
import re
import shlex
import shutil
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from circyto import __version__
from circyto.manifest.ciri_long import (
    CIRI_LONG_INTERPRETATION_BOUNDARY,
    CIRI_LONG_LIBRARY_PREPARATION,
    CIRI_LONG_SCHEMA_VERSION,
    CiriLongManifestRow,
    read_ciri_long_manifest_tsv,
)
from circyto.pipeline.workflow_reporting import utc_now_iso, write_json


BWA_INDEX_SUFFIXES = (".amb", ".ann", ".bwt", ".pac", ".sa")
OFFICIAL_COORDINATE_SYSTEM = "CIRI-long/CIRI-series 1-based closed"
NORMALIZED_COORDINATE_SYSTEM = "0-based half-open"
COORDINATE_CONVERSION_RULE = (
    "normalized_start=original_start-1; normalized_end=original_end"
)
CHEMISTRY_GATE_DECISION = "accepted_rcrt_circrna_enriched_cdna"

BSJ_SCHEMA_VERSION = "circyto.ciri_long_bsj.v1"
ISOFORM_SCHEMA_VERSION = "circyto.ciri_long_isoform.v1"
EXPRESSION_SCHEMA_VERSION = "circyto.ciri_long_expression.v1"
ISOFORM_USAGE_SCHEMA_VERSION = "circyto.ciri_long_isoform_usage.v1"
READ_ASSIGNMENT_SCHEMA_VERSION = "circyto.ciri_long_read_assignment.v1"

BSJ_COLUMNS = [
    "schema_version",
    "circRNA_id",
    "chromosome",
    "start",
    "end",
    "strand",
    "support",
    "circRNA_type",
    "host_gene_id",
    "host_gene_name",
    "host_gene_type",
    "source_cohort",
    "original_start",
    "original_end",
    "original_coordinate_system",
    "normalized_coordinate_system",
    "coordinate_conversion_rule",
    "official_source",
    "official_feature_type",
    "official_score",
    "splice_site",
    "equivalent_sequence",
    "major_isoform_length",
    "original_attributes_json",
    "original_info_line",
]

ISOFORM_COLUMNS = [
    "schema_version",
    "parent_circRNA_id",
    "isoform_id",
    "isoform_ordinal",
    "exon_block_structure",
    "transcript_length",
    "major_isoform_status",
    "host_gene_id",
    "host_gene_name",
    "host_gene_type",
    "source_cohort",
    "original_exon_block_structure",
    "original_coordinate_system",
    "normalized_coordinate_system",
    "coordinate_conversion_rule",
    "original_circRNA_attributes_json",
]

EXPRESSION_COLUMNS = [
    "schema_version",
    "circRNA_id",
    "sample_id",
    "support",
    "source_cohort",
    "original_value",
]

ISOFORM_USAGE_COLUMNS = [
    "schema_version",
    "parent_circRNA_id",
    "isoform_id",
    "sample_id",
    "isoform_usage",
    "source_cohort",
    "original_value",
]

_SAFE_PREFIX = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
_ATTRIBUTE = re.compile(r"([A-Za-z0-9_.:-]+)\s+\"([^\"]*)\"")


@dataclass(frozen=True)
class ToolVersion:
    executable: str
    version: str
    raw_output: str

    def to_dict(self) -> dict[str, str]:
        return {
            "executable": self.executable,
            "version": self.version,
            "raw_output": self.raw_output,
        }


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _file_identity(path: Path) -> dict[str, Any]:
    resolved = path.resolve()
    return {
        "path": str(resolved),
        "sha256": sha256_file(resolved),
        "bytes": resolved.stat().st_size,
    }


def _resolve_executable(value: str, *, label: str) -> str:
    candidate = Path(value).expanduser()
    if candidate.parent != Path(".") or candidate.is_absolute():
        if not candidate.is_file():
            raise FileNotFoundError(f"{label} executable not found: {candidate}")
        return str(candidate.resolve())
    resolved = shutil.which(value)
    if resolved is None:
        raise FileNotFoundError(f"{label} executable not found on PATH: {value}")
    return str(Path(resolved).resolve())


def _first_nonempty_line(text: str) -> str:
    for line in text.splitlines():
        stripped = line.strip()
        if stripped:
            return stripped
    return "unknown"


def detect_ciri_long_version(executable: str) -> ToolVersion:
    resolved = _resolve_executable(executable, label="CIRI-long")
    result = subprocess.run(
        [resolved, "--version"],
        capture_output=True,
        text=True,
        check=False,
        shell=False,
    )
    output = "\n".join(
        part.strip() for part in (result.stdout, result.stderr) if part and part.strip()
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"`{resolved} --version` failed with exit code {result.returncode}: "
            f"{output or 'no output'}"
        )
    line = _first_nonempty_line(output)
    match = re.search(r"(?i)(?:CIRI-long\s*)?(?:version\s*)?v?(\d+(?:\.\d+)+(?:[-+._A-Za-z0-9]*)?)", line)
    version = match.group(1) if match else line
    return ToolVersion(resolved, version, output)


def detect_bwa_version(executable: str) -> ToolVersion:
    resolved = _resolve_executable(executable, label="BWA")
    result = subprocess.run(
        [resolved],
        capture_output=True,
        text=True,
        check=False,
        shell=False,
    )
    output = "\n".join(
        part.strip() for part in (result.stdout, result.stderr) if part and part.strip()
    )
    match = re.search(r"(?im)^\s*Version:\s*([^\s]+)", output)
    if match:
        version = match.group(1)
    else:
        version = _first_nonempty_line(output)
    if version == "unknown":
        raise RuntimeError(f"Unable to determine BWA version from {resolved}")
    return ToolVersion(resolved, version, output)


def bwa_index_paths(reference_fasta: Path) -> list[Path]:
    reference = reference_fasta.resolve()
    return [Path(f"{reference}{suffix}") for suffix in BWA_INDEX_SUFFIXES]


def check_ciri_long_readiness(
    *,
    ciri_long: str = "CIRI-long",
    bwa: str = "bwa",
    reference_fasta: Path | None = None,
    gtf: Path | None = None,
    circ_annotation: Path | None = None,
) -> dict[str, Any]:
    errors: list[str] = []
    warnings: list[str] = []
    ciri_version: ToolVersion | None = None
    bwa_version: ToolVersion | None = None
    try:
        ciri_version = detect_ciri_long_version(ciri_long)
    except (FileNotFoundError, RuntimeError, OSError) as exc:
        errors.append(str(exc))
    try:
        bwa_version = detect_bwa_version(bwa)
    except (FileNotFoundError, RuntimeError, OSError) as exc:
        errors.append(str(exc))

    reference: dict[str, Any] | None = None
    index_assets: list[dict[str, Any]] = []
    if reference_fasta is not None:
        if not reference_fasta.is_file():
            errors.append(f"Reference FASTA not found: {reference_fasta}")
        else:
            reference = _file_identity(reference_fasta)
            for path in bwa_index_paths(reference_fasta):
                if path.is_file():
                    index_assets.append(_file_identity(path))
                else:
                    errors.append(f"Missing BWA reference index asset: {path}")

    annotation: dict[str, Any] | None = None
    if gtf is not None:
        if not gtf.is_file():
            errors.append(f"GTF annotation not found: {gtf}")
        else:
            annotation = _file_identity(gtf)
    else:
        warnings.append(
            "No GTF supplied; host-gene and exon/isoform annotation may be limited."
        )

    circ_asset: dict[str, Any] | None = None
    if circ_annotation is not None:
        if not circ_annotation.is_file():
            errors.append(
                f"Optional circRNA annotation asset not found: {circ_annotation}"
            )
        else:
            circ_asset = _file_identity(circ_annotation)

    return {
        "ok": not errors,
        "errors": errors,
        "warnings": warnings,
        "tools": {
            "ciri_long": ciri_version.to_dict() if ciri_version else None,
            "bwa": bwa_version.to_dict() if bwa_version else None,
        },
        "reference": reference,
        "bwa_index_assets": index_assets,
        "gtf": annotation,
        "circ_annotation": circ_asset,
    }


def _validate_prefix(value: str, *, label: str) -> str:
    if not _SAFE_PREFIX.fullmatch(value):
        raise ValueError(
            f"Unsafe {label}={value!r}; use letters, digits, '.', '_' or '-'"
        )
    return value


def build_ciri_long_call_argv(
    *,
    executable: str,
    reads_path: Path,
    output_dir: Path,
    reference_fasta: Path,
    prefix: str,
    threads: int,
    gtf: Path | None = None,
    circ_annotation: Path | None = None,
) -> list[str]:
    if threads < 1:
        raise ValueError("threads must be >= 1")
    _validate_prefix(prefix, label="sample prefix")
    argv = [
        executable,
        "call",
        "-i",
        str(reads_path),
        "-o",
        str(output_dir),
        "-r",
        str(reference_fasta),
        "-p",
        prefix,
    ]
    if gtf is not None:
        argv.extend(["-a", str(gtf)])
    if circ_annotation is not None:
        argv.extend(["-c", str(circ_annotation)])
    argv.extend(["-t", str(threads)])
    return argv


def build_ciri_long_collapse_argv(
    *,
    executable: str,
    calls_list: Path,
    output_dir: Path,
    reference_fasta: Path,
    prefix: str,
    threads: int,
    gtf: Path | None = None,
    circ_annotation: Path | None = None,
) -> list[str]:
    if threads < 1:
        raise ValueError("threads must be >= 1")
    _validate_prefix(prefix, label="cohort prefix")
    argv = [
        executable,
        "collapse",
        "-i",
        str(calls_list),
        "-o",
        str(output_dir),
        "-p",
        prefix,
        "-r",
        str(reference_fasta),
    ]
    if gtf is not None:
        argv.extend(["-a", str(gtf)])
    if circ_annotation is not None:
        argv.extend(["-c", str(circ_annotation)])
    argv.extend(["-t", str(threads)])
    return argv


def _require_one_reference(rows: Sequence[CiriLongManifestRow]) -> tuple[str, str]:
    identities = {(row.reference_id, row.reference_build) for row in rows}
    if len(identities) != 1:
        rendered = ", ".join(f"{ref}:{build}" for ref, build in sorted(identities))
        raise ValueError(
            "A CIRI-long call/collapse cohort must use one explicit reference "
            f"identity; found {rendered}"
        )
    return next(iter(identities))


def _official_call_outputs(sample_dir: Path, prefix: str) -> dict[str, str]:
    return {
        "candidate_fasta": str((sample_dir / f"{prefix}.cand_circ.fa").resolve()),
        "call_json": str((sample_dir / f"{prefix}.json").resolve()),
        "call_log": str((sample_dir / f"{prefix}.log").resolve()),
        "low_confidence_fasta": str(
            (sample_dir / f"{prefix}.low_confidence.fa").resolve()
        ),
        "ccs_fasta": str((sample_dir / "tmp" / f"{prefix}.ccs.fa").resolve()),
        "raw_fasta": str((sample_dir / "tmp" / f"{prefix}.raw.fa").resolve()),
        "splice_site_index": str((sample_dir / "tmp" / "ss.idx").resolve()),
    }


def _official_collapse_outputs(collapse_dir: Path, prefix: str) -> dict[str, str]:
    return {
        "info": str((collapse_dir / f"{prefix}.info").resolve()),
        "expression": str((collapse_dir / f"{prefix}.expression").resolve()),
        "isoforms": str((collapse_dir / f"{prefix}.isoforms").resolve()),
        "reads": str((collapse_dir / f"{prefix}.reads").resolve()),
        "collapse_log": str((collapse_dir / f"{prefix}.log").resolve()),
        "splice_site_index": str((collapse_dir / "tmp" / "ss.idx").resolve()),
        "corrected_pickle": str(
            (collapse_dir / "tmp" / f"{prefix}.corrected.pkl").resolve()
        ),
    }


def _output_artifacts(outputs: Mapping[str, str]) -> dict[str, dict[str, Any]]:
    artifacts: dict[str, dict[str, Any]] = {}
    for name, value in outputs.items():
        path = Path(value)
        if path.is_file():
            artifacts[name] = {**_file_identity(path), "exists": True}
        else:
            artifacts[name] = {"path": str(path.resolve()), "exists": False}
    return artifacts


def _base_provenance(
    *,
    stage: str,
    readiness: Mapping[str, Any],
    reference_id: str,
    reference_build: str,
    argv: Sequence[str],
    started_at: str,
    status: str,
    warnings: Iterable[str],
) -> dict[str, Any]:
    return {
        "schema_version": "circyto.ciri_long_provenance.v1",
        "circyto_version": __version__,
        "workflow_type": "ciri_long_rcrt_adapter",
        "stage": stage,
        "status": status,
        "started_at": started_at,
        "completed_at": utc_now_iso(),
        "shell": False,
        "argv": list(argv),
        "argv_shell_rendering_for_display_only": shlex.join(list(argv)),
        "ciri_long": readiness.get("tools", {}).get("ciri_long"),
        "bwa": readiness.get("tools", {}).get("bwa"),
        "reference": {
            **dict(readiness.get("reference") or {}),
            "reference_id": reference_id,
            "reference_build": reference_build,
        },
        "bwa_index_assets": list(readiness.get("bwa_index_assets") or []),
        "gtf": readiness.get("gtf"),
        "circ_annotation": readiness.get("circ_annotation"),
        "chemistry_gate": {
            "decision": CHEMISTRY_GATE_DECISION,
            "required_library_preparation": CIRI_LONG_LIBRARY_PREPARATION,
            "circRNA_enrichment_required": True,
            "biological_interpretation_boundary": CIRI_LONG_INTERPRETATION_BOUNDARY,
        },
        "warnings": list(warnings),
        "rejected_assumptions": [
            "ordinary poly(A)-selected cDNA is not accepted",
            "direct RNA is not accepted",
            "generic minimap2 BAM/SAM input is not accepted",
            "bulk samples are not represented as single cells",
        ],
    }


def _run_and_log(
    argv: Sequence[str],
    *,
    stdout_path: Path,
    stderr_path: Path,
) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(
        list(argv),
        capture_output=True,
        text=True,
        check=False,
        shell=False,
    )
    stdout_path.write_text(result.stdout or "", encoding="utf-8")
    stderr_path.write_text(result.stderr or "", encoding="utf-8")
    return result


def run_ciri_long_call_stage(
    *,
    manifest_path: Path,
    reference_fasta: Path,
    outdir: Path,
    gtf: Path | None = None,
    circ_annotation: Path | None = None,
    threads: int = 1,
    ciri_long: str = "CIRI-long",
    bwa: str = "bwa",
    execute: bool = False,
) -> dict[str, Any]:
    rows = read_ciri_long_manifest_tsv(manifest_path, validate_files=True)
    reference_id, reference_build = _require_one_reference(rows)
    readiness = check_ciri_long_readiness(
        ciri_long=ciri_long,
        bwa=bwa,
        reference_fasta=reference_fasta,
        gtf=gtf,
        circ_annotation=circ_annotation,
    )
    if not readiness["ok"]:
        raise RuntimeError(
            "CIRI-long readiness failed:\n- " + "\n- ".join(readiness["errors"])
        )
    ciri_executable = str(readiness["tools"]["ciri_long"]["executable"])
    root = outdir / "ciri_long"
    call_root = root / "call"
    call_root.mkdir(parents=True, exist_ok=True)
    call_manifest_path = call_root / "call_manifest.tsv"
    call_summary_path = call_root / "call_summary.json"
    resume_planned = False
    if call_manifest_path.exists():
        if execute and call_summary_path.is_file():
            previous = json.loads(call_summary_path.read_text(encoding="utf-8"))
            resume_planned = (
                previous.get("status") == "planned"
                and previous.get("executed") is False
            )
        if not resume_planned:
            raise FileExistsError(
                f"CIRI-long call output already exists: {call_manifest_path}"
            )

    call_rows: list[dict[str, str]] = []
    stage_started = utc_now_iso()
    for row in rows:
        sample_dir = call_root / row.sample_id
        if sample_dir.exists() and any(sample_dir.iterdir()):
            existing = {item.name for item in sample_dir.iterdir()}
            if not (resume_planned and existing <= {"provenance.json"}):
                raise FileExistsError(
                    f"CIRI-long sample output already exists for {row.sample_id}: "
                    f"{sample_dir}"
                )
        sample_dir.mkdir(parents=True, exist_ok=True)
        reads_path = Path(row.reads_path)
        argv = build_ciri_long_call_argv(
            executable=ciri_executable,
            reads_path=reads_path,
            output_dir=sample_dir,
            reference_fasta=reference_fasta.resolve(),
            prefix=row.sample_id,
            threads=threads,
            gtf=gtf.resolve() if gtf else None,
            circ_annotation=circ_annotation.resolve() if circ_annotation else None,
        )
        official_outputs = _official_call_outputs(sample_dir, row.sample_id)
        wrapper_stdout = sample_dir / "wrapper.stdout.log"
        wrapper_stderr = sample_dir / "wrapper.stderr.log"
        started_at = utc_now_iso()
        start_time = time.monotonic()
        status = "planned"
        returncode: int | None = None
        if execute:
            result = _run_and_log(
                argv,
                stdout_path=wrapper_stdout,
                stderr_path=wrapper_stderr,
            )
            returncode = int(result.returncode)
            status = "completed" if result.returncode == 0 else "failed"
            if returncode == 0:
                required = [
                    Path(official_outputs["candidate_fasta"]),
                    Path(official_outputs["low_confidence_fasta"]),
                ]
                missing = [str(path) for path in required if not path.is_file()]
                if missing:
                    status = "incomplete_outputs"
        provenance = _base_provenance(
            stage="call",
            readiness=readiness,
            reference_id=reference_id,
            reference_build=reference_build,
            argv=argv,
            started_at=started_at,
            status=status,
            warnings=readiness["warnings"],
        )
        provenance.update(
            {
                "elapsed_seconds": round(time.monotonic() - start_time, 3),
                "sample_id": row.sample_id,
                "source_accession": row.source_accession,
                "dataset_id": row.dataset_id,
                "input_reads": _file_identity(reads_path),
                "manifest_schema_version": row.schema_version,
                "manifest_path": str(manifest_path.resolve()),
                "official_output_files": _output_artifacts(official_outputs),
                "wrapper_stdout_log": str(wrapper_stdout.resolve()),
                "wrapper_stderr_log": str(wrapper_stderr.resolve()),
                "returncode": returncode,
                "executed": execute,
            }
        )
        provenance_path = sample_dir / "provenance.json"
        write_json(provenance_path, provenance)
        if execute and returncode != 0:
            raise RuntimeError(
                f"CIRI-long call failed for {row.sample_id} with exit code "
                f"{returncode}; see {wrapper_stderr}"
            )
        if execute and status == "incomplete_outputs":
            raise RuntimeError(
                f"CIRI-long call completed for {row.sample_id} but required "
                "official outputs are missing: " + ", ".join(missing)
            )
        call_rows.append(
            {
                "sample_id": row.sample_id,
                "candidate_fasta": official_outputs["candidate_fasta"],
                "low_confidence_fasta": official_outputs["low_confidence_fasta"],
                "call_json": official_outputs["call_json"],
                "call_log": official_outputs["call_log"],
                "provenance": str(provenance_path.resolve()),
                "source_accession": row.source_accession,
                "dataset_id": row.dataset_id,
                "reference_id": row.reference_id,
                "reference_build": row.reference_build,
                "status": status,
            }
        )

    with call_manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(call_rows[0]),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(call_rows)
    summary = {
        "schema_version": "circyto.ciri_long_call_summary.v1",
        "stage": "call",
        "status": "completed" if execute else "planned",
        "executed": execute,
        "started_at": stage_started,
        "completed_at": utc_now_iso(),
        "manifest": str(manifest_path.resolve()),
        "call_manifest": str(call_manifest_path.resolve()),
        "sample_count": len(call_rows),
        "samples": call_rows,
        "readiness": readiness,
    }
    summary_path = call_summary_path
    write_json(summary_path, summary)
    summary["summary_path"] = str(summary_path.resolve())
    return summary


def _read_call_manifest(
    path: Path,
    *,
    validate_files: bool,
) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(f"CIRI-long call manifest not found: {path}")
    required = {
        "sample_id",
        "candidate_fasta",
        "reference_id",
        "reference_build",
    }
    rows: list[dict[str, str]] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = sorted(required - set(reader.fieldnames or []))
        if missing:
            raise ValueError(
                f"CIRI-long call manifest missing columns at {path}: "
                + ", ".join(missing)
            )
        seen: set[str] = set()
        for line_number, raw in enumerate(reader, start=2):
            row = {str(key): (value or "").strip() for key, value in raw.items()}
            sample_id = row["sample_id"]
            _validate_prefix(sample_id, label="sample_id")
            if sample_id in seen:
                raise ValueError(
                    f"{path}:{line_number}: duplicate sample_id={sample_id!r}"
                )
            seen.add(sample_id)
            candidate = Path(row["candidate_fasta"]).expanduser()
            if not candidate.is_absolute():
                candidate = (path.parent / candidate).resolve()
            if re.search(r"\s", str(candidate)):
                raise ValueError(
                    "CIRI-long collapse sample-list paths cannot contain whitespace: "
                    f"{candidate}"
                )
            if validate_files and not candidate.is_file():
                raise FileNotFoundError(
                    f"{path}:{line_number}: candidate_fasta not found: {candidate}"
                )
            row["candidate_fasta"] = str(candidate)
            rows.append(row)
    if not rows:
        raise ValueError(f"CIRI-long call manifest contains 0 rows: {path}")
    return rows


def run_ciri_long_collapse_stage(
    *,
    call_manifest_path: Path,
    reference_fasta: Path,
    outdir: Path,
    prefix: str = "cohort",
    gtf: Path | None = None,
    circ_annotation: Path | None = None,
    threads: int = 1,
    ciri_long: str = "CIRI-long",
    bwa: str = "bwa",
    execute: bool = False,
) -> dict[str, Any]:
    _validate_prefix(prefix, label="cohort prefix")
    call_rows = _read_call_manifest(call_manifest_path, validate_files=execute)
    identities = {
        (row["reference_id"], row["reference_build"]) for row in call_rows
    }
    if len(identities) != 1:
        raise ValueError(
            "CIRI-long call manifest must contain one reference identity; found "
            + ", ".join(f"{ref}:{build}" for ref, build in sorted(identities))
        )
    reference_id, reference_build = next(iter(identities))
    readiness = check_ciri_long_readiness(
        ciri_long=ciri_long,
        bwa=bwa,
        reference_fasta=reference_fasta,
        gtf=gtf,
        circ_annotation=circ_annotation,
    )
    if not readiness["ok"]:
        raise RuntimeError(
            "CIRI-long readiness failed:\n- " + "\n- ".join(readiness["errors"])
        )
    root = outdir / "ciri_long"
    collapse_dir = root / "collapse"
    collapse_dir.mkdir(parents=True, exist_ok=True)
    provenance_path = collapse_dir / "provenance.json"
    summary_path = collapse_dir / "collapse_summary.json"
    resume_planned = False
    if provenance_path.exists():
        if execute and summary_path.is_file():
            previous = json.loads(summary_path.read_text(encoding="utf-8"))
            resume_planned = (
                previous.get("status") == "planned"
                and previous.get("executed") is False
            )
        existing = {item.name for item in collapse_dir.iterdir()}
        allowed = {"calls.list", "provenance.json", "collapse_summary.json"}
        if not (resume_planned and existing <= allowed):
            raise FileExistsError(
                f"CIRI-long collapse output already exists: {provenance_path}"
            )
    calls_list = collapse_dir / "calls.list"
    calls_list.write_text(
        "".join(
            f"{row['sample_id']} {row['candidate_fasta']}\n" for row in call_rows
        ),
        encoding="utf-8",
    )
    ciri_executable = str(readiness["tools"]["ciri_long"]["executable"])
    argv = build_ciri_long_collapse_argv(
        executable=ciri_executable,
        calls_list=calls_list.resolve(),
        output_dir=collapse_dir.resolve(),
        reference_fasta=reference_fasta.resolve(),
        prefix=prefix,
        threads=threads,
        gtf=gtf.resolve() if gtf else None,
        circ_annotation=circ_annotation.resolve() if circ_annotation else None,
    )
    outputs = _official_collapse_outputs(collapse_dir, prefix)
    stdout_log = collapse_dir / "wrapper.stdout.log"
    stderr_log = collapse_dir / "wrapper.stderr.log"
    started_at = utc_now_iso()
    start_time = time.monotonic()
    status = "planned"
    returncode: int | None = None
    if execute:
        result = _run_and_log(
            argv,
            stdout_path=stdout_log,
            stderr_path=stderr_log,
        )
        returncode = int(result.returncode)
        status = "completed" if result.returncode == 0 else "failed"
        if returncode == 0:
            required = [
                Path(outputs["info"]),
                Path(outputs["expression"]),
                Path(outputs["isoforms"]),
                Path(outputs["reads"]),
            ]
            missing = [str(path) for path in required if not path.is_file()]
            if missing:
                status = "incomplete_outputs"
    provenance = _base_provenance(
        stage="collapse",
        readiness=readiness,
        reference_id=reference_id,
        reference_build=reference_build,
        argv=argv,
        started_at=started_at,
        status=status,
        warnings=readiness["warnings"],
    )
    provenance.update(
        {
            "elapsed_seconds": round(time.monotonic() - start_time, 3),
            "cohort_prefix": prefix,
            "call_manifest": _file_identity(call_manifest_path),
            "calls_list": _file_identity(calls_list),
            "source_samples": [
                {
                    "sample_id": row["sample_id"],
                    "candidate_fasta": row["candidate_fasta"],
                    "candidate_fasta_sha256": (
                        sha256_file(Path(row["candidate_fasta"]))
                        if Path(row["candidate_fasta"]).is_file()
                        else None
                    ),
                }
                for row in call_rows
            ],
            "official_output_files": _output_artifacts(outputs),
            "wrapper_stdout_log": str(stdout_log.resolve()),
            "wrapper_stderr_log": str(stderr_log.resolve()),
            "returncode": returncode,
            "executed": execute,
        }
    )
    write_json(provenance_path, provenance)
    if execute and returncode != 0:
        raise RuntimeError(
            f"CIRI-long collapse failed with exit code {returncode}; see {stderr_log}"
        )
    if execute and status == "incomplete_outputs":
        raise RuntimeError(
            "CIRI-long collapse completed but required official outputs are "
            "missing: " + ", ".join(missing)
        )
    summary = {
        "schema_version": "circyto.ciri_long_collapse_summary.v1",
        "stage": "collapse",
        "status": status,
        "executed": execute,
        "cohort_prefix": prefix,
        "sample_count": len(call_rows),
        "calls_list": str(calls_list.resolve()),
        "provenance": str(provenance_path.resolve()),
        "official_output_files": outputs,
    }
    write_json(summary_path, summary)
    summary["summary_path"] = str(summary_path.resolve())
    return summary


def parse_ciri_long_attributes(text: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for match in _ATTRIBUTE.finditer(text):
        attributes[match.group(1)] = match.group(2)
    return attributes


def _parse_positive_coordinate(value: str, *, path: Path, line_number: int) -> int:
    try:
        coordinate = int(value)
    except ValueError as exc:
        raise ValueError(
            f"{path}:{line_number}: coordinate must be an integer; got {value!r}"
        ) from exc
    if coordinate < 1:
        raise ValueError(
            f"{path}:{line_number}: CIRI-long coordinates must be >= 1; got {coordinate}"
        )
    return coordinate


def _normalize_isoform_structure(
    structure: str,
    *,
    path: Path,
    line_number: int,
) -> tuple[list[list[int]], int]:
    blocks: list[list[int]] = []
    transcript_length = 0
    for raw_block in structure.split(","):
        match = re.fullmatch(r"(\d+)-(\d+)", raw_block.strip())
        if not match:
            raise ValueError(
                f"{path}:{line_number}: malformed CIRI-long isoform block "
                f"{raw_block!r} in {structure!r}"
            )
        start = int(match.group(1))
        end = int(match.group(2))
        if start < 1 or end < start:
            raise ValueError(
                f"{path}:{line_number}: invalid 1-based closed exon block "
                f"{raw_block!r}"
            )
        blocks.append([start - 1, end])
        transcript_length += end - start + 1
    return blocks, transcript_length


def parse_ciri_long_info(
    path: Path,
    *,
    cohort: str,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    bsj_rows: list[dict[str, Any]] = []
    isoform_rows: list[dict[str, Any]] = []
    seen_circs: set[str] = set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.rstrip("\n")
            if not stripped:
                continue
            fields = stripped.split("\t")
            if len(fields) != 9:
                raise ValueError(
                    f"{path}:{line_number}: expected 9 GTF-like columns; "
                    f"found {len(fields)}"
                )
            chrom, source, feature_type, start_raw, end_raw, score, strand, _, attr_text = fields
            original_start = _parse_positive_coordinate(
                start_raw, path=path, line_number=line_number
            )
            original_end = _parse_positive_coordinate(
                end_raw, path=path, line_number=line_number
            )
            if original_end < original_start:
                raise ValueError(
                    f"{path}:{line_number}: end precedes start: "
                    f"{original_start}>{original_end}"
                )
            attributes = parse_ciri_long_attributes(attr_text)
            circ_id = attributes.get("circ_id", "").strip()
            if not circ_id:
                raise ValueError(
                    f"{path}:{line_number}: CIRI-long info row lacks circ_id attribute"
                )
            if circ_id in seen_circs:
                raise ValueError(
                    f"{path}:{line_number}: duplicate circ_id={circ_id!r}"
                )
            seen_circs.add(circ_id)
            bsj_rows.append(
                {
                    "schema_version": BSJ_SCHEMA_VERSION,
                    "circRNA_id": circ_id,
                    "chromosome": chrom,
                    "start": original_start - 1,
                    "end": original_end,
                    "strand": strand,
                    "support": score,
                    "circRNA_type": attributes.get("circ_type", ""),
                    "host_gene_id": attributes.get("gene_id", ""),
                    "host_gene_name": attributes.get("gene_name", ""),
                    "host_gene_type": attributes.get("gene_type", ""),
                    "source_cohort": cohort,
                    "original_start": original_start,
                    "original_end": original_end,
                    "original_coordinate_system": OFFICIAL_COORDINATE_SYSTEM,
                    "normalized_coordinate_system": NORMALIZED_COORDINATE_SYSTEM,
                    "coordinate_conversion_rule": COORDINATE_CONVERSION_RULE,
                    "official_source": source,
                    "official_feature_type": feature_type,
                    "official_score": score,
                    "splice_site": attributes.get("splice_site", ""),
                    "equivalent_sequence": attributes.get("equivalent_seq", ""),
                    "major_isoform_length": attributes.get("circ_len", ""),
                    "original_attributes_json": json.dumps(
                        attributes, sort_keys=True, separators=(",", ":")
                    ),
                    "original_info_line": stripped,
                }
            )
            structures = [
                item.strip()
                for item in attributes.get("isoform", "").split("|")
                if item.strip()
            ]
            for ordinal, structure in enumerate(structures, start=1):
                blocks, transcript_length = _normalize_isoform_structure(
                    structure,
                    path=path,
                    line_number=line_number,
                )
                isoform_rows.append(
                    {
                        "schema_version": ISOFORM_SCHEMA_VERSION,
                        "parent_circRNA_id": circ_id,
                        "isoform_id": f"{circ_id}|{structure}",
                        "isoform_ordinal": ordinal,
                        "exon_block_structure": json.dumps(
                            blocks, separators=(",", ":")
                        ),
                        "transcript_length": transcript_length,
                        "major_isoform_status": (
                            "reported_major" if ordinal == 1 else "other_reported_isoform"
                        ),
                        "host_gene_id": attributes.get("gene_id", ""),
                        "host_gene_name": attributes.get("gene_name", ""),
                        "host_gene_type": attributes.get("gene_type", ""),
                        "source_cohort": cohort,
                        "original_exon_block_structure": structure,
                        "original_coordinate_system": OFFICIAL_COORDINATE_SYSTEM,
                        "normalized_coordinate_system": NORMALIZED_COORDINATE_SYSTEM,
                        "coordinate_conversion_rule": COORDINATE_CONVERSION_RULE,
                        "original_circRNA_attributes_json": json.dumps(
                            attributes, sort_keys=True, separators=(",", ":")
                        ),
                    }
                )
    return bsj_rows, isoform_rows


def _read_matrix_long(
    path: Path,
    *,
    id_column: str,
    value_name: str,
    schema_version: str,
    cohort: str,
    include_parent: bool,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = list(reader.fieldnames or [])
        if id_column not in header:
            raise ValueError(
                f"{path}: expected identifier column {id_column!r}; found {header}"
            )
        samples = [column for column in header if column != id_column]
        if not samples:
            return []
        for line_number, raw in enumerate(reader, start=2):
            feature_id = (raw.get(id_column) or "").strip()
            if not feature_id:
                raise ValueError(f"{path}:{line_number}: empty {id_column}")
            for sample_id in samples:
                value = (raw.get(sample_id) or "").strip()
                if value == "":
                    raise ValueError(
                        f"{path}:{line_number}: missing value for sample {sample_id!r}"
                    )
                try:
                    float(value)
                except ValueError as exc:
                    raise ValueError(
                        f"{path}:{line_number}: non-numeric value {value!r} "
                        f"for sample {sample_id!r}"
                    ) from exc
                item = {
                    "schema_version": schema_version,
                    (
                        "isoform_id"
                        if include_parent
                        else "circRNA_id"
                    ): feature_id,
                    "sample_id": sample_id,
                    value_name: value,
                    "source_cohort": cohort,
                    "original_value": value,
                }
                if include_parent:
                    item["parent_circRNA_id"] = feature_id.split("|", 1)[0]
                rows.append(item)
    return rows


def parse_ciri_long_expression(
    path: Path,
    *,
    cohort: str,
) -> list[dict[str, str]]:
    return _read_matrix_long(
        path,
        id_column="circ_ID",
        value_name="support",
        schema_version=EXPRESSION_SCHEMA_VERSION,
        cohort=cohort,
        include_parent=False,
    )


def parse_ciri_long_isoform_usage(
    path: Path,
    *,
    cohort: str,
) -> list[dict[str, str]]:
    return _read_matrix_long(
        path,
        id_column="isoform_ID",
        value_name="isoform_usage",
        schema_version=ISOFORM_USAGE_SCHEMA_VERSION,
        cohort=cohort,
        include_parent=True,
    )


def parse_ciri_long_read_assignments(
    path: Path,
    *,
    cohort: str,
) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        official_columns = list(reader.fieldnames or [])
        required = {"read_id", "circ_id", "sample"}
        missing = sorted(required - set(official_columns))
        if missing:
            raise ValueError(
                f"{path}: CIRI-long reads table missing columns: "
                + ", ".join(missing)
            )
        rows: list[dict[str, str]] = []
        for line_number, raw in enumerate(reader, start=2):
            item = {column: (raw.get(column) or "") for column in official_columns}
            if not item["read_id"].strip():
                raise ValueError(f"{path}:{line_number}: empty read_id")
            original = dict(item)
            item.update(
                {
                    "schema_version": READ_ASSIGNMENT_SCHEMA_VERSION,
                    "source_cohort": cohort,
                    "original_row_json": json.dumps(
                        original, sort_keys=True, separators=(",", ":")
                    ),
                }
            )
            rows.append(item)
    return official_columns, rows


def _write_tsv(
    path: Path,
    rows: Iterable[Mapping[str, Any]],
    *,
    fieldnames: Sequence[str],
) -> int:
    row_list = list(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fieldnames),
            delimiter="\t",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(row_list)
    return len(row_list)


def normalize_ciri_long_outputs(
    *,
    collapse_dir: Path,
    outdir: Path,
    prefix: str = "cohort",
) -> dict[str, Any]:
    _validate_prefix(prefix, label="cohort prefix")
    source_paths = {
        key: Path(value)
        for key, value in _official_collapse_outputs(collapse_dir, prefix).items()
        if key in {"info", "expression", "isoforms", "reads"}
    }
    missing = [str(path) for path in source_paths.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError(
            "Missing required CIRI-long collapse outputs: " + ", ".join(missing)
        )
    bsj_rows, isoform_rows = parse_ciri_long_info(
        source_paths["info"], cohort=prefix
    )
    expression_rows = parse_ciri_long_expression(
        source_paths["expression"], cohort=prefix
    )
    usage_rows = parse_ciri_long_isoform_usage(
        source_paths["isoforms"], cohort=prefix
    )
    read_columns, read_rows = parse_ciri_long_read_assignments(
        source_paths["reads"], cohort=prefix
    )

    normalized_dir = outdir / "ciri_long" / "normalized"
    normalized_dir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "circRNA_bsj": normalized_dir / "circRNA_bsj.tsv",
        "circRNA_isoforms": normalized_dir / "circRNA_isoforms.tsv",
        "circRNA_expression": normalized_dir / "circRNA_expression.tsv",
        "circRNA_isoform_usage": normalized_dir / "circRNA_isoform_usage.tsv",
        "circRNA_read_assignments": normalized_dir
        / "circRNA_read_assignments.tsv",
    }
    for path in outputs.values():
        if path.exists():
            raise FileExistsError(f"Normalized CIRI-long output already exists: {path}")
    counts = {
        "bsj_records": _write_tsv(
            outputs["circRNA_bsj"], bsj_rows, fieldnames=BSJ_COLUMNS
        ),
        "isoform_records": _write_tsv(
            outputs["circRNA_isoforms"],
            isoform_rows,
            fieldnames=ISOFORM_COLUMNS,
        ),
        "expression_records": _write_tsv(
            outputs["circRNA_expression"],
            expression_rows,
            fieldnames=EXPRESSION_COLUMNS,
        ),
        "isoform_usage_records": _write_tsv(
            outputs["circRNA_isoform_usage"],
            usage_rows,
            fieldnames=ISOFORM_USAGE_COLUMNS,
        ),
        "read_assignment_records": _write_tsv(
            outputs["circRNA_read_assignments"],
            read_rows,
            fieldnames=read_columns
            + ["schema_version", "source_cohort", "original_row_json"],
        ),
    }
    summary = {
        "schema_version": "circyto.ciri_long_import_summary.v1",
        "circyto_version": __version__,
        "created_at": utc_now_iso(),
        "cohort_prefix": prefix,
        "source_outputs": {
            key: _file_identity(path) for key, path in source_paths.items()
        },
        "upstream_collapse_provenance": (
            _file_identity(collapse_dir / "provenance.json")
            if (collapse_dir / "provenance.json").is_file()
            else None
        ),
        "normalized_outputs": {
            key: _file_identity(path) for key, path in outputs.items()
        },
        "record_counts": counts,
        "coordinate_policy": {
            "original_coordinate_system": OFFICIAL_COORDINATE_SYSTEM,
            "normalized_coordinate_system": NORMALIZED_COORDINATE_SYSTEM,
            "conversion_rule": COORDINATE_CONVERSION_RULE,
            "official_identifiers_preserved": True,
            "official_representations_preserved": True,
        },
        "entity_separation": {
            "bsj_level": "circRNA_bsj.tsv",
            "isoform_level": "circRNA_isoforms.tsv",
            "expression_level": "circRNA_expression.tsv",
            "isoform_usage_level": "circRNA_isoform_usage.tsv",
            "read_evidence_level": "circRNA_read_assignments.tsv",
        },
        "warnings": [
            "CIRI-long detector output normalization is not independent biological validation.",
            "Bulk samples are not represented as single cells.",
            "Isoform usage is preserved as reported and is not an absolute read count.",
        ],
    }
    summary_path = normalized_dir / "import_summary.json"
    write_json(summary_path, summary)
    summary["summary_path"] = str(summary_path.resolve())
    return summary
