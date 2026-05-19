from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from circyto.pipeline.workflow_reporting import load_json


REQUIRED_PROVENANCE_FIELDS = [
    "command_name",
    "circyto_version",
    "workflow_type",
    "protocol",
    "read_layout",
    "genome_fasta",
    "gtf",
    "detector_backend",
    "started_at",
    "completed_at",
    "elapsed_seconds",
    "hostname",
    "python_version",
    "cleanup_summary",
    "workflow_uuid",
    "completed_stages",
    "skipped_stages",
    "failed_stages",
    "partial_outputs_detected",
]


def _check_json(path: Path, errors: list[str]) -> dict[str, Any] | None:
    if not path.exists():
        errors.append(f"Missing JSON file: {path}")
        return None
    try:
        return load_json(path)
    except json.JSONDecodeError:
        errors.append(f"Corrupted JSON file: {path}")
        return None


def _required_paths_for_workflow(summary: dict[str, Any], workdir: Path) -> dict[str, Path]:
    paths = dict(summary.get("paths", {}))
    required: dict[str, Path] = {"workflow_summary": workdir / "workflow_summary.json"}
    workflow_type = str(summary.get("workflow_type", summary.get("workflow", "")))
    if workflow_type == "smartseq3-ciri3":
        for key in ("selected_manifest", "alignment_manifest", "detector_summary", "matrix", "cell_index", "circ_index"):
            value = paths.get(key)
            if value:
                required[key] = Path(str(value))
    else:
        for key in ("manifest", "alignment_manifest", "detector_summary", "matrix", "cell_index", "circ_index"):
            value = paths.get(key)
            if value:
                required[key] = Path(str(value))
    h5ad = paths.get("h5ad")
    if h5ad:
        required["h5ad"] = Path(str(h5ad))
    return required


def check_workflow_integrity(workdir: Path) -> dict[str, Any]:
    summary_path = workdir / "workflow_summary.json"
    errors: list[str] = []
    warnings: list[str] = []
    details: dict[str, Any] = {
        "workdir": str(workdir.resolve()),
        "required_outputs": {},
        "provenance_files_checked": [],
    }

    summary = _check_json(summary_path, errors)
    if summary is None:
        return {"ok": False, "workdir": str(workdir.resolve()), "errors": errors, "warnings": warnings, "details": details}

    missing_fields = [field for field in REQUIRED_PROVENANCE_FIELDS if field not in summary]
    if missing_fields:
        errors.append(f"workflow_summary.json missing required provenance fields: {', '.join(missing_fields)}")

    required_paths = _required_paths_for_workflow(summary, workdir)
    for label, path in required_paths.items():
        details["required_outputs"][label] = str(path)
        if not path.exists():
            errors.append(f"Missing required output: {label} -> {path}")

    alignment_manifest = None
    manifest_path = required_paths.get("alignment_manifest")
    if manifest_path and manifest_path.exists():
        alignment_manifest = manifest_path
    else:
        for candidate in (
            workdir / "align" / "alignment_manifest.tsv",
            workdir / "align" / "star" / "alignment_manifest.tsv",
            workdir / "align" / "bwa_mem" / "alignment_manifest.tsv",
        ):
            if candidate.exists():
                alignment_manifest = candidate
                break
    if alignment_manifest is None:
        errors.append("No alignment manifest found under completed workflow directory.")

    matrix_path = required_paths.get("matrix")
    if matrix_path and matrix_path.exists():
        try:
            with matrix_path.open("r", encoding="utf-8") as handle:
                next((line for line in handle if line.strip() and not line.startswith("%")), None)
        except OSError:
            errors.append(f"Could not read matrix file: {matrix_path}")

    rna_dir = workdir / "rna"
    if rna_dir.exists():
        for name in ("gene_counts.tsv", "gene_feature_table.tsv", "rna_import_summary.json"):
            path = rna_dir / name
            if not path.exists():
                errors.append(f"RNA directory present but missing {name}")
        rna_summary = _check_json(rna_dir / "rna_import_summary.json", errors)
        if rna_summary is not None:
            if "method" not in rna_summary:
                errors.append("rna_import_summary.json missing method")

    for prov_path in sorted(workdir.rglob("*.provenance.json")):
        details["provenance_files_checked"].append(str(prov_path.relative_to(workdir)))
        _check_json(prov_path, errors)

    partial_outputs = list(summary.get("partial_outputs_detected", []))
    if partial_outputs:
        warnings.append(f"partial_outputs_detected: {', '.join(partial_outputs)}")

    return {
        "ok": len(errors) == 0,
        "workdir": str(workdir.resolve()),
        "workflow_type": summary.get("workflow_type", summary.get("workflow")),
        "workflow_uuid": summary.get("workflow_uuid"),
        "errors": errors,
        "warnings": warnings,
        "details": details,
    }
