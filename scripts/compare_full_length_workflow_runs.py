#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import Any


def _load_summary(script_dir: Path, workflow_dir: Path, *, top_files: int, top_errors: int) -> dict[str, Any]:
    script_path = script_dir / "check_full_length_workflow_outputs.py"
    result = subprocess.run(
        [
            sys.executable,
            str(script_path),
            str(workflow_dir),
            "--json",
            "--top-files",
            str(max(1, top_files)),
            "--top-errors",
            str(max(1, top_errors)),
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        stderr = result.stderr.strip() or result.stdout.strip() or f"failed to inspect {workflow_dir}"
        raise RuntimeError(stderr)
    return json.loads(result.stdout)


def _h5ad_size_bytes(summary: dict[str, Any]) -> int:
    h5ad_path = summary.get("h5ad_path")
    if not h5ad_path:
        return 0
    path = Path(str(h5ad_path))
    if not path.exists():
        return 0
    try:
        return int(path.stat().st_size)
    except OSError:
        return 0


def _top_recurrent_ids(summary: dict[str, Any], *, limit: int = 5) -> list[str]:
    items = summary.get("top_recurrent_circRNAs", [])
    if not isinstance(items, list):
        return []
    return [str(item.get("circ_id", "")) for item in items[:limit] if str(item.get("circ_id", "")).strip()]


def _warning_lines(summary: dict[str, Any]) -> list[str]:
    warnings: list[str] = []
    summary_warnings = summary.get("warnings", [])
    if isinstance(summary_warnings, list):
        warnings.extend(str(item) for item in summary_warnings if str(item).strip())
    obvious_error_lines = summary.get("obvious_error_lines", [])
    if isinstance(obvious_error_lines, list):
        warnings.extend(str(item) for item in obvious_error_lines if str(item).strip())
    return warnings


def _extract_metrics(summary: dict[str, Any]) -> dict[str, Any]:
    matrix = summary.get("matrix", {})
    detector = summary.get("detector", {})
    circ_summary = detector.get("circRNAs_per_cell_summary", {})
    return {
        "workflow_dir": summary.get("workflow_dir"),
        "protocol": summary.get("protocol"),
        "planned_cells": int(summary.get("planned_cells", 0) or 0),
        "completed_cells": int(summary.get("completed_cells", 0) or 0),
        "failed_cells": int(summary.get("failed_cells", 0) or 0),
        "circRNAs": int(matrix.get("n_circRNAs", matrix.get("n_rows", 0) or 0) or 0),
        "cells": int(matrix.get("n_cells", matrix.get("n_cols", 0) or 0) or 0),
        "nnz": int(matrix.get("nnz", 0) or 0),
        "circRNAs_per_cell_summary": circ_summary,
        "h5ad_size_bytes": _h5ad_size_bytes(summary),
        "workdir_size_bytes": int(summary.get("workdir_size_bytes", 0) or 0),
        "align_size_bytes": int(summary.get("align_size_bytes", 0) or 0),
        "top_recurrent_circRNAs": summary.get("top_recurrent_circRNAs", []),
        "top_recurrent_circRNA_ids": _top_recurrent_ids(summary),
        "warnings_and_errors": _warning_lines(summary),
    }


def build_comparison(
    left_summary: dict[str, Any],
    right_summary: dict[str, Any],
    *,
    left_label: str,
    right_label: str,
) -> dict[str, Any]:
    left_metrics = _extract_metrics(left_summary)
    right_metrics = _extract_metrics(right_summary)
    return {
        "left_label": left_label,
        "right_label": right_label,
        "left": left_metrics,
        "right": right_metrics,
        "delta": {
            "cells": right_metrics["cells"] - left_metrics["cells"],
            "circRNAs": right_metrics["circRNAs"] - left_metrics["circRNAs"],
            "nnz": right_metrics["nnz"] - left_metrics["nnz"],
            "completed_cells": right_metrics["completed_cells"] - left_metrics["completed_cells"],
            "failed_cells": right_metrics["failed_cells"] - left_metrics["failed_cells"],
            "h5ad_size_bytes": right_metrics["h5ad_size_bytes"] - left_metrics["h5ad_size_bytes"],
            "workdir_size_bytes": right_metrics["workdir_size_bytes"] - left_metrics["workdir_size_bytes"],
            "align_size_bytes": right_metrics["align_size_bytes"] - left_metrics["align_size_bytes"],
        },
        "shared_top_recurrent_circRNAs": sorted(
            set(left_metrics["top_recurrent_circRNA_ids"]) & set(right_metrics["top_recurrent_circRNA_ids"])
        ),
    }


def _print_side(label: str, metrics: dict[str, Any]) -> None:
    circ_summary = metrics.get("circRNAs_per_cell_summary", {})
    print(f"{label}: {metrics.get('workflow_dir')}")
    print(
        f"  cells={metrics.get('cells')} completed={metrics.get('completed_cells')} "
        f"failed={metrics.get('failed_cells')} circRNAs={metrics.get('circRNAs')} nnz={metrics.get('nnz')}"
    )
    print(
        f"  circRNAs_per_cell={circ_summary if circ_summary else 'n/a'}"
    )
    print(
        f"  h5ad_size_bytes={metrics.get('h5ad_size_bytes')} "
        f"workdir_size_bytes={metrics.get('workdir_size_bytes')} "
        f"align_size_bytes={metrics.get('align_size_bytes')}"
    )
    print(
        f"  top_recurrent_circRNAs={metrics.get('top_recurrent_circRNA_ids', [])}"
    )
    warnings = metrics.get("warnings_and_errors", [])
    print("  warnings_errors:")
    if warnings:
        for line in warnings[:8]:
            print(f"  - {line}")
    else:
        print("  - none")


def _print_text(comparison: dict[str, Any]) -> None:
    _print_side(comparison["left_label"], comparison["left"])
    print("")
    _print_side(comparison["right_label"], comparison["right"])
    print("")
    delta = comparison["delta"]
    print("Delta (right - left):")
    print(
        f"  cells={delta['cells']} circRNAs={delta['circRNAs']} nnz={delta['nnz']} "
        f"completed_cells={delta['completed_cells']} failed_cells={delta['failed_cells']}"
    )
    print(
        f"  h5ad_size_bytes={delta['h5ad_size_bytes']} "
        f"workdir_size_bytes={delta['workdir_size_bytes']} align_size_bytes={delta['align_size_bytes']}"
    )
    print(f"Shared top recurrent circRNAs: {comparison['shared_top_recurrent_circRNAs']}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Read-only comparison of two completed full-length workflow directories.")
    parser.add_argument("left_workflow_dir", type=Path, help="First workflow output directory")
    parser.add_argument("right_workflow_dir", type=Path, help="Second workflow output directory")
    parser.add_argument("--left-label", default="left", help="Display label for the first workflow")
    parser.add_argument("--right-label", default="right", help="Display label for the second workflow")
    parser.add_argument("--json", action="store_true", help="Emit JSON instead of human-readable text.")
    parser.add_argument("--top-files", type=int, default=8, help="Largest-file limit passed to the inspector.")
    parser.add_argument("--top-errors", type=int, default=12, help="Error-line limit passed to the inspector.")
    args = parser.parse_args(argv)

    script_dir = Path(__file__).resolve().parent
    try:
        left_summary = _load_summary(script_dir, args.left_workflow_dir, top_files=args.top_files, top_errors=args.top_errors)
        right_summary = _load_summary(script_dir, args.right_workflow_dir, top_files=args.top_files, top_errors=args.top_errors)
        comparison = build_comparison(
            left_summary,
            right_summary,
            left_label=args.left_label,
            right_label=args.right_label,
        )
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    if args.json:
        print(json.dumps(comparison, indent=2, sort_keys=True))
    else:
        _print_text(comparison)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
