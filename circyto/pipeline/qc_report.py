"""Offline, read-only interpretation of existing workflow QC artifacts.

Only report.html and metrics.json are written. No scientific objects, detector
outputs, or source summaries are loaded through pipeline execution code.
"""
from __future__ import annotations

import csv
import hashlib
import html
import io
import json
import os
import re
import tempfile
from collections import Counter
from decimal import Decimal, InvalidOperation
from pathlib import Path
from statistics import median
from typing import Any, Callable

from circyto import __version__


SCHEMA_VERSION = "circyto.qc_report.v1"
REPORT_FEATURE = "qc-report-usability (unreleased; outside the v0.10.0 manuscript baseline)"
SUCCESS = {"success", "empty", "skipped_existing", "aligned", "reused_input", "reused_cached", "completed"}
CELL_COMPLETION = {
    "alignment": {"aligned", "reused_input", "reused_cached"},
    "detector": {"success", "empty", "skipped_existing"},
}
FAILED = {"failed", "error", "partial_failure", "aborted"}
OUTPUTS = {
    "matrix/circ_counts.mtx": "circRNA candidates × cells (Matrix Market)",
    "matrix/circ_index.txt": "Matrix row identities, in order",
    "matrix/cell_index.txt": "Matrix column identities, in order",
    "matrix/circ_feature_table.tsv": "Candidate coordinates and annotations",
    "anndata/circ_counts.h5ad": "AnnData: cells × circRNA candidates (optional)",
    "mudata/circyto_multimodal.h5mu": "MuData: aligned modalities (optional)",
    "qc/cell_qc.tsv": "Per-cell QC, including selected zero-count cells",
    "qc/circ_qc.tsv": "Per-candidate QC and host-gene annotations",
}
SOURCES = ("workflow_summary.json", "qc/cell_qc.tsv", "qc/circ_qc.tsv",
           "ciri3/detector_run_summary.json", "align/alignment_prepare_summary.json",
           "rna/rna_import_summary.json")


def _text(value: Any) -> str:
    """Limit exported free text to scalars; redact paths and network locations."""
    if not isinstance(value, (str, int, float, bool)):
        return "Not available"
    value = str(value)
    value = re.sub(r"(?:https?|file|s3)://[^\s<>\"']+", "[location omitted]", value)
    value = re.sub(r"(?<![\w<])(?:[A-Za-z]:[\\/]|/|~/)[^\s<>\"']+", "[path omitted]", value)
    return "".join(c for c in value if c >= " " or c in "\n\t")


def _escape(value: Any) -> str:
    return html.escape(str(value), quote=True)


def _integer(value: Any) -> int | None:
    if value is None or isinstance(value, bool):
        return None
    try:
        number = Decimal(str(value))
        if number.is_finite() and number >= 0 and number == number.to_integral_value():
            return int(number)
    except (InvalidOperation, ValueError):
        pass
    return None


class _Inputs:
    def __init__(self, root: Path):
        self.root = root
        self.sources: list[dict[str, Any]] = []
        self.warnings: list[str] = []

    def read(self, name: str, *, optional: bool = False) -> str | None:
        source: dict[str, Any] = {"path": name, "available": False}
        self.sources.append(source)
        path = self.root / name
        try:
            if not path.resolve().is_relative_to(self.root):
                raise ValueError("Symlink points outside the workflow directory")
            data = path.read_bytes()
            text = data.decode("utf-8-sig")
        except (OSError, UnicodeError, ValueError) as exc:
            reason = "File is missing" if isinstance(exc, FileNotFoundError) else "File is unreadable or outside the workflow directory"
            source["reason"] = reason
            if not isinstance(exc, FileNotFoundError):
                source["valid"] = False
            if not optional or not isinstance(exc, FileNotFoundError):
                self.warnings.append(f"{name}: {reason}.")
            return None
        source.update(available=True, sha256=hashlib.sha256(data).hexdigest(), bytes=len(data))
        return text

    def invalid(self, name: str, reason: str) -> None:
        self.warnings.append(f"{name}: {reason}.")
        for source in self.sources:
            if source["path"] == name:
                source.update(valid=False, reason=reason)

    def json(self, name: str, *, optional: bool = False) -> dict[str, Any]:
        text = self.read(name, optional=optional)
        if text is None:
            return {}
        try:
            value = json.loads(text)
            if not isinstance(value, dict):
                raise ValueError("Expected an object")
            return value
        except (ValueError, RecursionError):
            self.invalid(name, "Invalid JSON object; values are unavailable")
            return {}

    def table(self, name: str, id_column: str) -> tuple[list[dict[str, str]] | None, set[str]]:
        text = self.read(name)
        if text is None:
            return None, set()
        try:
            reader = csv.DictReader(io.StringIO(text), delimiter="\t", strict=True)
            fields = reader.fieldnames or []
            rows = list(reader)
            ids = [row.get(id_column, "") for row in rows]
            if (id_column not in fields or len(set(fields)) != len(fields)
                    or any(None in row or None in row.values() for row in rows)
                    or any(not value.strip() for value in ids) or len(set(ids)) != len(ids)):
                raise ValueError("Invalid table")
            return rows, set(fields)
        except (csv.Error, ValueError):
            self.invalid(name, f"Malformed TSV or missing/duplicate {id_column}; table metrics are unavailable")
            return None, set()


def _metric(value: Any, unit: str, source: str, definition: str, reason: str | None = None) -> dict[str, Any]:
    return {"value": value, "unit": unit, "source": source, "definition": definition,
            "reason": reason if value is None else None}


def _fraction(numerator: int | None, denominator: int | None, source: str, definition: str) -> dict[str, Any]:
    result = _metric(numerator / denominator if numerator is not None and denominator else None,
                     "fraction", source, definition, "A complete numerator and a nonzero denominator are required.")
    result.update(numerator=numerator, denominator=denominator)
    return result


def _counts(rows: list[dict[str, str]] | None, fields: set[str], column: str) -> list[int | None] | None:
    if rows is None or column not in fields:
        return None
    return [_integer(row[column]) for row in rows]


def _complete(values: list[int | None] | None) -> bool:
    return values is not None and all(value is not None for value in values)


def _status_counts(value: Any) -> dict[str, int]:
    if not isinstance(value, dict):
        return {}
    return {_text(key): count for key, raw in value.items() if (count := _integer(raw)) is not None}


def collect_report_metrics(workdir: Path) -> dict[str, Any]:
    """Read fixed local artifact names, with no dependency on original paths."""
    root = Path(workdir).resolve()
    if not root.is_dir():
        raise ValueError("--workdir must be an existing workflow directory")
    inputs = _Inputs(root)
    summary = inputs.json("workflow_summary.json")
    cells, cell_fields = inputs.table("qc/cell_qc.tsv", "cell_id")
    candidates, circ_fields = inputs.table("qc/circ_qc.tsv", "circ_id")
    detector = inputs.json("ciri3/detector_run_summary.json", optional=True)
    alignment = inputs.json("align/alignment_prepare_summary.json", optional=True)
    rna = inputs.json("rna/rna_import_summary.json", optional=True)
    n_cells = len(cells) if cells is not None else None
    n_candidates = len(candidates) if candidates is not None else None
    counts = _counts(cells, cell_fields, "circRNA_count")
    support = _counts(cells, cell_fields, "total_circRNA_support")
    detected = _counts(candidates, circ_fields, "n_cells_detected")
    circ_support = _counts(candidates, circ_fields, "total_support")
    cell_source = "qc/cell_qc.tsv"
    circ_source = "qc/circ_qc.tsv"
    missing = "Source table/column is missing, malformed, or contains missing/invalid nonnegative integer values."
    metrics = {
        "cells": _metric(n_cells, "cells", cell_source, "Number of unique cell_id rows in cell_qc.tsv.", missing),
        "candidates": _metric(n_candidates, "candidates", circ_source, "Number of unique circ_id rows in circ_qc.tsv.", missing),
        "total_support": _metric(sum(support) if _complete(support) else None, "support", cell_source,
                                 "Sum of total_circRNA_support over every QC cell; detector support units, not molecules or UMIs.", missing),
        "candidate_total_support": _metric(sum(circ_support) if _complete(circ_support) else None, "support", circ_source,
                                           "Independent sum of total_support over every candidate QC row.", missing),
        "median_candidates_per_cell": _metric(median(counts) if _complete(counts) and counts else None, "candidates/cell", cell_source,
                                               "Median circRNA_count over every QC cell, including zero-count cells.", missing),
        "zero_candidate_cells": _metric(counts.count(0) if _complete(counts) else None, "cells", cell_source,
                                        "Number of QC cells with circRNA_count = 0; failed cells are not biological negatives.", missing),
    }
    has_host = candidates is not None and "host_gene" in circ_fields
    blank_hosts = {"", ".", "-", "na", "n/a", "nan", "none", "null", "unknown", "unassigned"}
    annotated = sum(row["host_gene"].strip().lower() not in blank_hosts for row in candidates) if has_host else None
    metrics["host_gene_annotated_candidates"] = _metric(annotated, "candidates", circ_source,
        "Candidate rows with a recorded, non-placeholder host_gene; annotation is not experimental validation.", "host_gene column or candidate table is unavailable.")
    metrics["host_gene_coverage"] = _fraction(annotated, n_candidates, circ_source,
        "Candidate rows with a recorded host_gene / all candidate QC rows. Does not establish annotation completeness or correctness.")

    checks: list[dict[str, Any]] = []

    def check(name: str, left: Any, right: Any) -> None:
        ok = left == right if left is not None and right is not None else None
        checks.append({"check": name, "left": left, "right": right, "ok": ok})
        if ok is False:
            inputs.warnings.append(f"QC consistency mismatch: {name} ({left} versus {right}).")

    check("cell support sum = candidate support sum", metrics["total_support"]["value"], metrics["candidate_total_support"]["value"])
    check("cell candidate counts sum = candidate detected-cell counts sum",
          sum(counts) if _complete(counts) else None, sum(detected) if _complete(detected) else None)
    matrix = summary.get("matrix") if isinstance(summary.get("matrix"), dict) else {}
    check("QC cells = summary matrix cells", n_cells, _integer(matrix.get("n_cells")))
    check("QC candidates = summary matrix candidates", n_candidates, _integer(matrix.get("n_circRNAs")))
    expected_cells = _integer(summary.get("selected_cell_count", summary.get("planned_cells")))
    check("QC cells = selected/planned cells", n_cells, expected_cells)
    prevalence_valid = n_cells is not None and bool(n_cells) and _complete(detected)
    if detected is not None and n_cells is not None and any(value is not None and value > n_cells for value in detected):
        inputs.warnings.append("Candidate n_cells_detected exceeds the QC cell denominator; prevalence is unavailable.")
        prevalence_valid = False
        checks.append({"check": "candidate prevalence bounds", "ok": False})
    if counts is not None and n_candidates is not None and any(value is not None and value > n_candidates for value in counts):
        inputs.warnings.append("A cell candidate count exceeds the number of candidate QC rows.")
        checks.append({"check": "cell candidate count bounds", "ok": False})
    if any(item["ok"] is False for item in checks):
        prevalence_valid = False

    per_cell = []
    for index, row in enumerate(cells or []):
        per_cell.append({"cell_id": row["cell_id"],
                         "candidate_count": counts[index] if counts is not None else None,
                         "total_support": support[index] if support is not None else None,
                         "detector_status": _text(row.get("detector_status", "unknown")),
                         "alignment_status": _text(row.get("alignment_status", "unknown"))})
    per_candidate = []
    for index, row in enumerate(candidates or []):
        value = detected[index] if detected is not None else None
        per_candidate.append({"circ_id": row["circ_id"], "n_cells_detected": value,
                              "total_support": circ_support[index] if circ_support is not None else None,
                              "host_gene": row["host_gene"] if has_host else None,
                              "prevalence": None})

    evidence: list[dict[str, Any]] = []
    for kind, data in (("alignment", alignment), ("detector", detector)):
        data_source = "align/alignment_prepare_summary.json" if kind == "alignment" else "ciri3/detector_run_summary.json"
        for source, raw in (("workflow_summary.json", summary.get(f"{kind}_status_counts")), (data_source, data.get("status_counts"))):
            if raw is not None and (not isinstance(raw, dict) or any(_integer(value) is None for value in raw.values())):
                inputs.invalid(source, f"Invalid {kind} status counts; completion cannot be established")
        records = data.get("cells")
        if isinstance(records, list) and records:
            statuses = Counter(_text(row.get("status", "unknown")) if isinstance(row, dict) else "unknown" for row in records)
            evidence.append({"kind": kind, "source": data_source + " (cell records)", "status_counts": dict(statuses)})
        for source, raw in (("workflow_summary.json", summary.get(f"{kind}_status_counts")),
                            (data_source, data.get("status_counts")),
                            (cell_source, dict(Counter(row[f"{kind}_status"] for row in per_cell)) if f"{kind}_status" in cell_fields else None)):
            status_counts = _status_counts(raw)
            if status_counts:
                evidence.append({"kind": kind, "source": source, "status_counts": status_counts})
        # Some older summaries provide failure counters without status_counts.
        nested = summary.get(kind) if isinstance(summary.get(kind), dict) else {}
        for source, record in (("workflow_summary.json", nested), (kind + " summary", data)):
            failed = _integer(record.get("failed_cells", record.get("failed")))
            if failed:
                evidence.append({"kind": kind, "source": source, "status_counts": {"failed": failed}})

    reasons: list[str] = []
    stages = summary.get("stage_graph", [])
    stages = [row for row in stages if isinstance(row, dict)] if isinstance(stages, list) else []
    stage_rows = [{"stage": _text(row.get("stage", "unknown")), "status": _text(row.get("status", "unknown"))} for row in stages]
    failed = bool(summary.get("failed_stages")) or str(summary.get("status", "")).lower() in FAILED
    failed = failed or any(row["status"].lower() in FAILED for row in stage_rows)
    for record in evidence:
        for status, count in record["status_counts"].items():
            if count and status.lower() in FAILED:
                failed = True
            elif count and status.lower() not in CELL_COMPLETION[record["kind"]]:
                reasons.append(f"{record['kind']} has cells with status {status}.")
    if failed:
        reasons.append("Recorded stage, alignment, or detector failures; QC may describe partial outputs.")
    if not summary:
        reasons.append("workflow_summary.json is missing or invalid; completion cannot be established.")
    elif not summary.get("completed_at"):
        reasons.append("No completion timestamp is recorded; completion cannot be established.")
    if summary.get("status") and str(summary["status"]).lower() not in SUCCESS:
        reasons.append("The summary does not record a completed workflow status.")
    if any(row["status"].lower() not in {"completed", "skipped", "disabled", "delegated"} for row in stage_rows):
        reasons.append("One or more workflow stages are not complete.")
    for kind in ("alignment", "detector"):
        records = [record for record in evidence if record["kind"] == kind]
        if not records:
            reasons.append(f"{kind.capitalize()} completion evidence is unavailable.")
        else:
            totals = {sum(record["status_counts"].values()) for record in records}
            if len(totals) > 1 or (n_cells is not None and any(total != n_cells for total in totals)):
                reasons.append(f"{kind.capitalize()} status totals do not match the QC cell population.")
            distributions = {tuple(sorted((key, value) for key, value in row["status_counts"].items() if value)) for row in records}
            if len(distributions) > 1:
                reasons.append(f"{kind.capitalize()} status distributions disagree between sources.")
    if cells is None or candidates is None:
        reasons.append("Required QC tables are unavailable.")
    if any(source.get("valid") is False for source in inputs.sources if source["path"] != "rna/rna_import_summary.json"):
        reasons.append("One or more source artifacts could not be interpreted reliably.")
    if any(item["ok"] is False for item in checks):
        reasons.append("QC sources are inconsistent; inspect the recorded values before interpretation.")

    # The aggregate candidate TSV cannot be restricted to a subset of cells.
    # Only use its denominator when the entire QC cohort was evaluated. Keep
    # raw QC counts (including failed/missing rows) intact for troubleshooting.
    evaluated_cells = sum(all(row[f"{kind}_status"].lower() in statuses
                              for kind, statuses in CELL_COMPLETION.items()) for row in per_cell)
    cohort_valid = bool(n_cells) and evaluated_cells == n_cells and not reasons and summary.get("dry_run") is not True
    cohort_reason = (
        "Detection rates require consistent completion evidence and completed alignment and detector statuses for every QC cell. "
        "Failed, missing, or unprocessed cells are not negative detections. Aggregate candidate QC cannot be restricted to an evaluated subset."
    )
    metrics["zero_candidate_fraction"] = _fraction(
        metrics["zero_candidate_cells"]["value"], n_cells if cohort_valid else None, cell_source,
        "Zero-count cells / all QC cells, available only when the entire QC cohort was evaluated; not biological absence.")
    if not cohort_valid:
        metrics["zero_candidate_fraction"]["reason"] = cohort_reason
        inputs.warnings.append(cohort_reason)
    prevalence_valid = prevalence_valid and cohort_valid and _complete(counts)
    prevalence_reason = None if prevalence_valid else (
        cohort_reason if not cohort_valid else "Complete, consistent candidate and cell counts are required for prevalence.")
    if prevalence_valid:
        for row in per_candidate:
            row["prevalence"] = row["n_cells_detected"] / n_cells
    outputs = [{"path": name, "description": description,
                "available": (root / name).is_file() and (root / name).resolve().is_relative_to(root)}
               for name, description in OUTPUTS.items()]
    if any(not row["available"] for row in outputs[:3]):
        reasons.append("One or more matrix/index files are missing from this copy of the workflow.")
    if summary.get("partial_outputs_detected"):
        inputs.warnings.append("The workflow summary recorded missing outputs. Inspect the original summary; some intermediates may have been cleaned.")
    state = "failed" if failed else "dry_run" if summary.get("dry_run") is True else "incomplete" if reasons else "completed"
    if summary.get("dry_run") is True:
        reasons.insert(0, "This is a planning-only dry run; detection has not been established.")
    if not summary and not failed:
        state = "unknown"
    recorded_warnings = summary.get("warnings", [])
    if isinstance(recorded_warnings, list):
        inputs.warnings.extend(_text(warning) for warning in recorded_warnings)
    if n_cells == 0:
        inputs.warnings.append("The cell QC table contains no cells; fractions requiring a cell denominator are unavailable.")
    if not prevalence_valid:
        inputs.warnings.append(f"Candidate prevalence is unavailable: {prevalence_reason}")
    for name, values in (("circRNA_count", counts), ("total_circRNA_support", support), ("n_cells_detected", detected), ("total_support", circ_support)):
        if not _complete(values):
            inputs.warnings.append(f"{name}: {missing}")

    provenance = {key: _text(summary[key]) for key in (
        "workflow_type", "workflow", "protocol", "read_layout", "detector_backend", "circyto_version",
        "python_version", "started_at", "completed_at", "workflow_uuid") if summary.get(key) is not None}
    # Reference basenames aid orientation without disclosing machine paths.
    for key in ("genome_fasta", "gtf"):
        value = summary.get(key)
        if isinstance(value, str) and value:
            provenance[key + "_basename"] = _text(value.replace("\\", "/").rsplit("/", 1)[-1])
    return {
        "schema_version": SCHEMA_VERSION, "reporter_version": __version__, "reporter_feature": REPORT_FEATURE,
        "workflow": {"status": state, "reasons": list(dict.fromkeys(reasons)), "stages": stage_rows},
        "metrics": metrics, "per_cell": per_cell, "per_candidate": per_candidate,
        "candidate_prevalence_denominator": n_cells if prevalence_valid else None,
        "candidate_prevalence_reason": prevalence_reason,
        "cell_evaluation": {"completed_status_cells": evaluated_cells if cells is not None else None,
                            "other_status_cells": n_cells - evaluated_cells if n_cells is not None else None,
                            "cohort_eligible_for_rates": cohort_valid},
        "detector_evidence": evidence,
        "consistency_checks": checks, "warnings": list(dict.fromkeys(inputs.warnings)),
        "provenance": provenance, "sources": inputs.sources, "outputs": outputs,
        "optional_summaries": {"rna": {key: _text(rna[key]) for key in ("method", "n_cells", "n_genes") if key in rna}},
    }


def _display(value: Any) -> str:
    if value is None:
        return "Not available"
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value)


def _table(headers: list[str], rows: list[list[Any]], *, limit: int = 100) -> str:
    if not rows:
        return '<p class="muted">Not available: no rows to display.</p>'
    header = "".join(f'<th scope="col">{_escape(value)}</th>' for value in headers)
    body = "".join("<tr>" + "".join(f"<td>{_escape(_display(value))}</td>" for value in row) + "</tr>" for row in rows[:limit])
    note = f'<p class="muted">Showing {min(limit, len(rows))} of {len(rows)} rows. All values are in metrics.json and the source TSV.</p>'
    return f'<div class="table-scroll"><table><thead><tr>{header}</tr></thead><tbody>{body}</tbody></table></div>{note}'


def _histogram(values: list[int | None], title: str, axis: str) -> str:
    if not values or any(value is None for value in values):
        return '<p class="muted">Not available: complete numeric QC values are required.</p>'
    frequencies = Counter(values)
    if len(frequencies) <= 12:
        bins = [(str(value), frequencies[value]) for value in sorted(frequencies)]
    else:
        width = max(1, (max(values) + 11) // 12)
        binned = Counter(value // width for value in values)
        bins = [(f"{index * width}–{(index + 1) * width - 1}", binned[index]) for index in sorted(binned)]
    peak = max(count for _, count in bins)
    spacing = 480 / len(bins)
    parts = [f'<svg viewBox="0 0 560 250" role="img" aria-label="{_escape(title)}"><title>{_escape(title)}</title>',
             '<line x1="48" y1="194" x2="535" y2="194" stroke="#ccd6df"/>']
    for index, (label, count) in enumerate(bins):
        x = 52 + index * spacing
        height = 135 * count / peak
        parts.append(f'<rect x="{x:.2f}" y="{194-height:.2f}" width="{spacing*.68:.2f}" height="{height:.2f}" rx="3" fill="#157d86"><title>{_escape(label)}: {count}</title></rect>')
        parts.append(f'<text x="{x+spacing*.34:.2f}" y="{185-height:.2f}" text-anchor="middle">{count}</text><text x="{x+spacing*.34:.2f}" y="214" text-anchor="middle">{_escape(label)}</text>')
    parts.append(f'<text x="290" y="243" text-anchor="middle">{_escape(axis)}</text></svg>')
    return "".join(parts)


def render_report(payload: dict[str, Any]) -> str:
    """Render escaped HTML with embedded CSS/SVG and no scripts or web assets."""
    from circyto.pipeline.qc_report_style import STYLE

    metrics = payload["metrics"]
    workflow = payload["workflow"]
    provenance = payload["provenance"]
    status = workflow["status"]
    cards = []
    for key, label in (("cells", "QC cells"), ("candidates", "circRNA candidates"), ("total_support", "Total support"),
                       ("median_candidates_per_cell", "Median candidates / cell"), ("zero_candidate_fraction", "Zero-count cells"),
                       ("host_gene_coverage", "Recorded host-gene coverage")):
        metric = metrics[key]
        value = _display(metric["value"])
        if metric["unit"] == "fraction" and metric["value"] is not None:
            value = f'{metric["value"]:.1%}'
        detail = metric["reason"] or metric["definition"]
        if metric["unit"] == "fraction" and metric["value"] is not None:
            detail = f'{metric["numerator"]} / {metric["denominator"]}. {detail}'
        cards.append(f'<article class="metric"><div class="label">{_escape(label)}</div><div class="value">{_escape(value)}</div><p>{_escape(detail)}</p></article>')
    reasons = "".join(f"<li>{_escape(reason)}</li>" for reason in workflow["reasons"])
    status_note = f"<ul>{reasons}</ul>" if reasons else "<p>Recorded completion and available QC/status evidence agree. This is a report-level check, not independent validation of the analysis.</p>"
    warnings = "".join(f"<li>{_escape(warning)}</li>" for warning in payload["warnings"])
    warnings = f"<ul>{warnings}</ul>" if warnings else "<p>No missing-value or consistency warnings detected in the report inputs.</p>"
    notice = (f'<p class="notice"><a href="#warnings">{len(payload["warnings"])} report warning(s)</a> · '
              f'{_escape(payload["warnings"][0][:400])}</p>') if payload["warnings"] else ""
    evidence_rows = [[row["kind"], status_name, count, row["source"]] for row in payload["detector_evidence"] for status_name, count in row["status_counts"].items()]
    cells = payload["per_cell"]
    candidates = payload["per_candidate"]
    cell_table = _table(["Cell", "Candidates", "Support", "Detector", "Alignment"],
        [[row["cell_id"], row["candidate_count"], row["total_support"], row["detector_status"], row["alignment_status"]] for row in cells])
    candidate_table = _table(["circRNA candidate", "Cells detected", "Prevalence", "Support", "Host gene"],
        [[row["circ_id"], row["n_cells_detected"], f'{row["prevalence"]:.2%}' if row["prevalence"] is not None else None,
          row["total_support"], row["host_gene"] or "Not available"] for row in candidates])
    output_table = _table(["Path relative to results directory", "Contents", "File present"],
        [[row["path"], row["description"], "Yes" if row["available"] else "No"] for row in payload["outputs"]])
    source_table = _table(["Source", "SHA-256 / availability"], [[row["path"], row.get("sha256", row.get("reason", "Not available"))] for row in payload["sources"]])
    provenance_table = _table(["Field", "Recorded value"], [[key, value] for key, value in provenance.items()])
    checks = _table(["Consistency check", "Left value", "Right value", "Result"],
        [[row["check"], row.get("left"), row.get("right"), "Match" if row["ok"] is True else "Mismatch" if row["ok"] is False else "Not available"] for row in payload["consistency_checks"]])
    stage_table = _table(["Recorded stage", "Recorded status"], [[row["stage"], row["status"]] for row in workflow["stages"]])
    optional = payload["optional_summaries"]["rna"]
    rna_table = _table(["RNA summary field", "Recorded value"], [[key, value] for key, value in optional.items()])
    definition_table = _table(["Metric", "Source", "Definition / missing-value reason"],
        [[key, value["source"], value["definition"] + (" " + value["reason"] if value["reason"] else "")] for key, value in metrics.items()])
    prevalence_note = ("Prevalence: Not available. " + payload["candidate_prevalence_reason"]
                       if payload["candidate_prevalence_reason"] else
                       f'Prevalence = n_cells_detected / all {payload["candidate_prevalence_denominator"]} evaluated QC cells.')
    return f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<meta http-equiv="Content-Security-Policy" content="default-src 'none'; style-src 'unsafe-inline'; img-src data:; base-uri 'none'; form-action 'none'">
<title>CIRCYTO · circRNA candidate QC</title><style>{STYLE}</style></head><body>
<header><div class="brand">CIRCYTO <span>WORKFLOW REPORT</span></div><nav aria-label="Report sections"><a href="#overview">Overview</a><a href="#evidence">Evidence</a><a href="#outputs">Outputs</a><a href="#provenance">Provenance</a></nav></header>
<main><section id="overview"><div class="eyebrow">OFFLINE QUALITY CONTROL</div><h1>circRNA candidate summary</h1>
<p class="subtitle">{_escape(provenance.get("workflow_type", provenance.get("workflow", "Workflow not recorded")))} · {_escape(provenance.get("protocol", "Protocol not recorded"))} · {_escape(provenance.get("detector_backend", "Detector not recorded"))}</p>
{notice}
<div class="status {status}"><h2>Workflow: {_escape(status.replace("_", " ").capitalize())}</h2>{status_note}</div>
<p class="interpretation">Detector-supported circRNA candidates require orthogonal confirmation. Workflow completion, detector evidence, and biological validation are separate questions.</p>
<div class="metrics">{"".join(cards)}</div></section>
<section><div class="section-heading"><h2>Candidate detection across cells</h2><span>Source QC values · no filtering</span></div><div class="charts"><article><h3>Per-cell candidate counts</h3>{_histogram([row["candidate_count"] for row in cells], "Distribution of circRNA candidate counts; labels above bars count cells", "circRNA candidates per cell")}</article><article><h3>Candidate prevalence</h3>{_histogram([row["n_cells_detected"] for row in candidates] if all(row["prevalence"] is not None for row in candidates) else [], "Distribution of detected-cell counts; labels above bars count candidates", "QC cells detecting each candidate")}</article></div>
<p class="muted">Bar labels give frequencies. {_escape(prevalence_note)} Percentages are rounded for display; exact counts and denominators are retained in metrics.json. Raw cell counts include every QC row; failed or unprocessed cells must not be interpreted as negative detections.</p></section>
<section id="evidence"><h2>Detector evidence and workflow stages</h2><p>Statuses are reported separately for each source; counts across sources must not be added together. “Empty” means the detector recorded no candidates; “skipped_existing” means existing output was reused.</p>{_table(["Stage", "Status", "Cells", "Source"], evidence_rows)}<details><summary>Recorded workflow stages</summary>{stage_table}</details></section>
<section id="warnings" class="warnings"><h2>Warnings and availability</h2>{warnings}<p>Host-gene coverage measures recorded annotation fields only. A blank field does not establish that no host gene exists. Missing metrics are never replaced with zero.</p></section>
<section><h2>Per-cell QC</h2>{cell_table}</section><section><h2>Per-candidate QC</h2>{candidate_table}</section>
<section id="outputs"><h2>Find and reuse your outputs</h2>{output_table}<p>The Matrix Market files retain their original candidate and cell order. Selected zero-count cells can appear in QC/AnnData but be absent from the collected matrix; use the matrix index files for its identities.</p>
<p>Open <code>qc/report.html</code> directly in a browser. To regenerate this report from existing outputs:</p><pre>circyto report --workdir RESULTS_DIRECTORY</pre><p><code>qc/metrics.json</code> contains all report values, definitions, checks, and source checksums. This command only writes these two report files.</p></section>
<section id="provenance"><h2>Provenance and reproducibility</h2><p>Report schema: {_escape(payload["schema_version"])} · Reporter package: {_escape(payload["reporter_version"])}. Feature: {_escape(payload["reporter_feature"])}. Source checksums identify the exact input bytes. Original run version and timestamps are shown below when recorded.</p>{provenance_table}
<details><summary>Source checksums</summary>{source_table}</details><details><summary>QC consistency checks</summary>{checks}</details><details><summary>Metric definitions</summary>{definition_table}</details><details><summary>Optional existing RNA summary</summary>{rna_table}</details>
<p class="muted">Machine hostname, full paths, command lines, and environment variables are omitted. Reference basenames are labels, not reference checksums. Cell and candidate identifiers are retained; review them before sharing. Original summaries remain in the results directory.</p></section>
</main><footer>CIRCYTO · Standalone HTML · No network, scripts, or external assets required</footer></body></html>\n'''


def generate_report(workdir: Path) -> dict[str, Any]:
    """Generate both artifacts after validating destinations; preserve sources."""
    root = Path(workdir).resolve()
    if not root.is_dir():
        raise ValueError("--workdir must be an existing workflow directory")
    qc = root / "qc"
    if qc.is_symlink() or (qc.exists() and not qc.is_dir()):
        raise ValueError("qc must be a real directory, not a symlink or file")
    targets = [qc / "metrics.json", qc / "report.html"]
    for path in targets:
        if path.is_symlink() or (path.exists() and not path.is_file()):
            raise ValueError(f"Refusing unsafe report destination: qc/{path.name}")
    if any((root / name).resolve() in targets for name in SOURCES):
        raise ValueError("A source artifact resolves to a report destination; refusing to replace source data")
    payload = collect_report_metrics(root)
    contents = [json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n", render_report(payload)]
    qc.mkdir(exist_ok=True)
    staged: list[Path] = []
    try:
        for content in contents:
            with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=qc, prefix=".report-", delete=False) as handle:
                staged.append(Path(handle.name))
                handle.write(content)
        for temp, target in zip(staged, targets):
            os.replace(temp, target)
    finally:
        for temp in staged:
            temp.unlink(missing_ok=True)
    return payload


def report_after_workflow(workdir: Path, *, progress: Callable[[str], None] = print) -> None:
    """Surface reporting errors after scientific outputs have been persisted."""
    try:
        payload = generate_report(workdir)
    except Exception as exc:
        message = ("HTML report generation failed after workflow outputs were written; scientific results were retained. "
                   "Fix the reporting error and rerun `circyto report --workdir RESULTS_DIRECTORY`. "
                   f"Cause: {type(exc).__name__}: {exc}")
        progress(message)
        raise RuntimeError(message) from exc
    progress(f"HTML report: {workdir / 'qc' / 'report.html'} (workflow status: {payload['workflow']['status']})")
    if payload["workflow"]["status"] != "completed":
        progress("WARNING: Report evidence does not establish successful completion; inspect warnings and detector statuses.")
