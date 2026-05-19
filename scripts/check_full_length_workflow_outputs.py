#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _format_bytes(size: int) -> str:
    units = ["B", "KB", "MB", "GB", "TB"]
    value = float(size)
    for unit in units:
        if value < 1024.0 or unit == units[-1]:
            if unit == "B":
                return f"{int(value)} {unit}"
            return f"{value:.1f} {unit}"
        value /= 1024.0
    return f"{size} B"


def _read_matrix_header(path: Path) -> dict[str, int] | None:
    if not path.exists():
        return None
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            text = line.strip()
            if not text or text.startswith("%"):
                continue
            parts = text.split()
            if len(parts) != 3:
                return None
            n_rows, n_cols, nnz = (int(part) for part in parts)
            return {"n_rows": n_rows, "n_cols": n_cols, "nnz": nnz}
    return None


def _largest_files(root: Path, *, limit: int) -> list[tuple[Path, int]]:
    ranked: list[tuple[Path, int]] = []
    for path in root.rglob("*"):
        if not path.is_file():
            continue
        try:
            size = path.stat().st_size
        except OSError:
            continue
        ranked.append((path, size))
    ranked.sort(key=lambda item: (-item[1], str(item[0])))
    return ranked[:limit]


def _error_lines(root: Path, *, limit: int) -> list[str]:
    patterns = ("error", "exception", "traceback", "failed", "fatal")
    lines: list[str] = []
    for path in sorted(root.rglob("*")):
        if not path.is_file():
            continue
        if path.suffix.lower() not in {".log", ".txt", ".out", ".err", ".json"} and "log" not in path.name.lower():
            continue
        try:
            text = path.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        for raw_line in text.splitlines():
            line = raw_line.strip()
            lower = line.lower()
            if any(pattern in lower for pattern in patterns):
                rel = path.relative_to(root)
                lines.append(f"{rel}: {line}")
                if len(lines) >= limit:
                    return lines
    return lines


def _top_cells(detector_summary: dict[str, Any], workflow_summary: dict[str, Any]) -> list[dict[str, Any]]:
    cells = detector_summary.get("cells", [])
    if isinstance(cells, list) and cells:
        ranked = sorted(
            (
                {
                    "cell_id": str(cell.get("cell_id", "")),
                    "normalized_row_count": int(cell.get("normalized_row_count", 0) or 0),
                    "raw_row_count": int(cell.get("raw_row_count", 0) or 0),
                    "status": str(cell.get("status", "")),
                }
                for cell in cells
                if str(cell.get("cell_id", "")).strip()
            ),
            key=lambda item: (-item["normalized_row_count"], -item["raw_row_count"], item["cell_id"]),
        )
        return ranked[:10]
    detector_block = workflow_summary.get("detector", {})
    top = detector_block.get("top_cells_by_circRNA_count", [])
    if isinstance(top, list):
        return list(top)[:10]
    return []


def _top_recurrent(workflow_summary: dict[str, Any]) -> list[dict[str, Any]]:
    matrix = workflow_summary.get("matrix", {})
    top = matrix.get("top_recurrent_circRNAs", [])
    if isinstance(top, list):
        return list(top)[:10]
    return []


def _matrix_dims(workflow_summary: dict[str, Any], root: Path) -> dict[str, int | None]:
    matrix = workflow_summary.get("matrix", {})
    dims = {
        "n_rows": matrix.get("n_rows"),
        "n_cols": matrix.get("n_cols"),
        "nnz": matrix.get("nnz"),
        "n_cells": matrix.get("n_cells"),
        "n_circRNAs": matrix.get("n_circRNAs"),
    }
    if all(value is not None for value in dims.values()):
        return dims
    mm_path = root / "matrix" / "circ_counts.mtx"
    header = _read_matrix_header(mm_path)
    if header is None:
        return dims
    dims["n_rows"] = header["n_rows"]
    dims["n_cols"] = header["n_cols"]
    dims["nnz"] = header["nnz"]
    return dims


def build_summary(root: Path, *, top_files: int, top_errors: int) -> dict[str, Any]:
    workflow_path = root / "workflow_summary.json"
    detector_path = root / "ciri3" / "detector_run_summary.json"
    if not workflow_path.exists():
        raise FileNotFoundError(f"Missing workflow summary: {workflow_path}")
    workflow_summary = _load_json(workflow_path)
    detector_summary = _load_json(detector_path) if detector_path.exists() else {}

    h5ad_path = workflow_summary.get("paths", {}).get("h5ad")
    h5ad_exists = bool(h5ad_path) and Path(str(h5ad_path)).exists()
    matrix_dims = _matrix_dims(workflow_summary, root)

    return {
        "workflow_dir": str(root.resolve()),
        "workflow": workflow_summary.get("workflow"),
        "protocol": workflow_summary.get("protocol"),
        "dry_run": bool(workflow_summary.get("dry_run", False)),
        "warnings": workflow_summary.get("warnings", []),
        "stage_graph": workflow_summary.get("stage_graph", []),
        "planned_cells": int(workflow_summary.get("planned_cells", 0) or 0),
        "completed_cells": int(workflow_summary.get("completed_cells", 0) or 0),
        "failed_cells": int(detector_summary.get("failed_cells", 0) or 0),
        "detector_status_counts": detector_summary.get("status_counts", workflow_summary.get("detector_status_counts", {})),
        "circRNAs_per_cell": _top_cells(detector_summary, workflow_summary),
        "matrix": matrix_dims,
        "top_recurrent_circRNAs": _top_recurrent(workflow_summary),
        "h5ad_exists": h5ad_exists,
        "h5ad_path": h5ad_path,
        "workdir_size_bytes": int(workflow_summary.get("workdir_size_bytes", 0) or 0),
        "align_size_bytes": int(workflow_summary.get("align_size_bytes", 0) or 0),
        "largest_files": [
            {"path": str(path.relative_to(root)), "bytes": size, "size_human": _format_bytes(size)}
            for path, size in _largest_files(root, limit=top_files)
        ],
        "obvious_error_lines": _error_lines(root, limit=top_errors),
    }


def _print_text(summary: dict[str, Any]) -> None:
    print(f"Workflow dir: {summary['workflow_dir']}")
    print(f"Workflow: {summary.get('workflow')}  protocol={summary.get('protocol')}  dry_run={summary.get('dry_run')}")
    print("")
    print("Stage graph:")
    for item in summary.get("stage_graph", []):
        stage = item.get("stage", "?")
        status = item.get("status", "?")
        detail = item.get("detail")
        if detail:
            print(f"- {stage}: {status} ({detail})")
        else:
            print(f"- {stage}: {status}")
    print("")
    print(
        "Cells: "
        f"planned={summary.get('planned_cells', 0)} "
        f"completed={summary.get('completed_cells', 0)} "
        f"failed={summary.get('failed_cells', 0)}"
    )
    print(f"Detector status counts: {summary.get('detector_status_counts', {})}")
    print("")
    print("CircRNAs per cell:")
    circ_per_cell = summary.get("circRNAs_per_cell", [])
    if circ_per_cell:
        for item in circ_per_cell:
            if "circRNA_count" in item:
                print(f"- {item.get('cell_id')}: circRNA_count={item.get('circRNA_count')}")
            else:
                print(
                    f"- {item.get('cell_id')}: normalized_row_count={item.get('normalized_row_count', 0)} "
                    f"raw_row_count={item.get('raw_row_count', 0)} status={item.get('status', '')}"
                )
    else:
        print("- none")
    print("")
    matrix = summary.get("matrix", {})
    print(
        "Matrix: "
        f"n_rows={matrix.get('n_rows')} n_cols={matrix.get('n_cols')} nnz={matrix.get('nnz')} "
        f"n_cells={matrix.get('n_cells')} n_circRNAs={matrix.get('n_circRNAs')}"
    )
    print("")
    print("Top recurrent circRNAs:")
    top = summary.get("top_recurrent_circRNAs", [])
    if top:
        for item in top:
            print(
                f"- {item.get('circ_id')}: n_cells_detected={item.get('n_cells_detected')} "
                f"total_support={item.get('total_support')}"
            )
    else:
        print("- none")
    print("")
    print(f"h5ad exists: {summary.get('h5ad_exists')}  path={summary.get('h5ad_path')}")
    print("")
    print("Largest files:")
    for item in summary.get("largest_files", []):
        print(f"- {item['path']}: {item['size_human']}")
    print("")
    print("Obvious error lines:")
    errors = summary.get("obvious_error_lines", [])
    if errors:
        for line in errors:
            print(f"- {line}")
    else:
        print("- none")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Read-only inspection of full-length-circrna workflow outputs.")
    parser.add_argument("workflow_dir", type=Path, help="Workflow output directory containing workflow_summary.json")
    parser.add_argument("--json", action="store_true", help="Emit JSON instead of human-readable text.")
    parser.add_argument("--top-files", type=int, default=8, help="Number of largest files to report.")
    parser.add_argument("--top-errors", type=int, default=12, help="Maximum number of obvious error lines to report.")
    args = parser.parse_args(argv)

    try:
        summary = build_summary(args.workflow_dir, top_files=max(1, args.top_files), top_errors=max(1, args.top_errors))
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    if args.json:
        print(json.dumps(summary, indent=2, sort_keys=True))
    else:
        _print_text(summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
