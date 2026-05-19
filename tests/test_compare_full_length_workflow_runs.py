from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


def _write_workflow_fixture(
    root: Path,
    *,
    cells: int,
    circRNAs: int,
    nnz: int,
    circ_counts: list[int],
    top_circs: list[tuple[str, int, int]],
    workdir_size_bytes: int,
    align_size_bytes: int,
    h5ad_bytes: int,
    warning: str | None = None,
    error_line: str | None = None,
) -> Path:
    (root / "ciri3").mkdir(parents=True)
    (root / "matrix").mkdir(parents=True)
    (root / "logs").mkdir(parents=True)
    (root / "anndata").mkdir(parents=True)
    (root / "align").mkdir(parents=True)

    top_recurrent = [
        {"circ_id": circ_id, "n_cells_detected": n_cells_detected, "total_support": total_support}
        for circ_id, n_cells_detected, total_support in top_circs
    ]
    detector_cells = [
        {
            "cell_id": f"cell{i+1}",
            "status": "success",
            "raw_row_count": count,
            "normalized_row_count": count,
        }
        for i, count in enumerate(circ_counts)
    ]
    payload = {
        "workflow": "full-length-circrna",
        "protocol": "ramda",
        "dry_run": False,
        "planned_cells": cells,
        "completed_cells": cells,
        "warnings": [warning] if warning else [],
        "detector": {
            "circRNAs_per_cell_summary": {
                "count": len(circ_counts),
                "min": min(circ_counts) if circ_counts else 0,
                "max": max(circ_counts) if circ_counts else 0,
                "mean": (sum(circ_counts) / len(circ_counts)) if circ_counts else 0,
            }
        },
        "matrix": {
            "n_rows": circRNAs,
            "n_cols": cells,
            "nnz": nnz,
            "n_cells": cells,
            "n_circRNAs": circRNAs,
            "top_recurrent_circRNAs": top_recurrent,
        },
        "workdir_size_bytes": workdir_size_bytes,
        "align_size_bytes": align_size_bytes,
        "paths": {"h5ad": str((root / "anndata" / "circ_counts.h5ad").resolve())},
    }
    (root / "workflow_summary.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    (root / "ciri3" / "detector_run_summary.json").write_text(
        json.dumps({"failed_cells": 0, "status_counts": {"success": cells}, "cells": detector_cells}, indent=2) + "\n",
        encoding="utf-8",
    )
    (root / "matrix" / "circ_counts.mtx").write_text(
        f"%%MatrixMarket matrix coordinate integer general\n%\n{circRNAs} {cells} {nnz}\n",
        encoding="utf-8",
    )
    (root / "anndata" / "circ_counts.h5ad").write_bytes(b"x" * h5ad_bytes)
    (root / "align" / "Aligned.out.sam").write_bytes(b"a" * max(1, align_size_bytes // 2))
    log_text = "INFO complete\n"
    if error_line:
        log_text += f"{error_line}\n"
    (root / "logs" / "workflow.log").write_text(log_text, encoding="utf-8")
    return root


def test_compare_full_length_workflow_runs_text(tmp_path: Path) -> None:
    left = _write_workflow_fixture(
        tmp_path / "left",
        cells=2,
        circRNAs=221,
        nnz=250,
        circ_counts=[110, 111],
        top_circs=[("circA", 2, 9), ("circB", 2, 5)],
        workdir_size_bytes=5000,
        align_size_bytes=2000,
        h5ad_bytes=512,
        warning="left warning",
        error_line="ERROR left log issue",
    )
    right = _write_workflow_fixture(
        tmp_path / "right",
        cells=2,
        circRNAs=240,
        nnz=280,
        circ_counts=[120, 118],
        top_circs=[("circA", 2, 10), ("circC", 2, 6)],
        workdir_size_bytes=7000,
        align_size_bytes=1500,
        h5ad_bytes=1024,
        warning="right warning",
        error_line="ERROR right log issue",
    )

    result = subprocess.run(
        [
            sys.executable,
            "scripts/compare_full_length_workflow_runs.py",
            str(left),
            str(right),
            "--left-label",
            "imr90",
            "--right-label",
            "hap1",
        ],
        cwd=Path(__file__).resolve().parents[1],
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert "imr90:" in result.stdout
    assert "hap1:" in result.stdout
    assert "cells=2 completed=2 failed=0 circRNAs=221 nnz=250" in result.stdout
    assert "cells=2 completed=2 failed=0 circRNAs=240 nnz=280" in result.stdout
    assert "Shared top recurrent circRNAs: ['circA']" in result.stdout
    assert "right warning" in result.stdout
    assert "ERROR right log issue" in result.stdout


def test_compare_full_length_workflow_runs_json(tmp_path: Path) -> None:
    left = _write_workflow_fixture(
        tmp_path / "left",
        cells=1,
        circRNAs=10,
        nnz=12,
        circ_counts=[10],
        top_circs=[("circX", 1, 4)],
        workdir_size_bytes=100,
        align_size_bytes=40,
        h5ad_bytes=20,
    )
    right = _write_workflow_fixture(
        tmp_path / "right",
        cells=2,
        circRNAs=12,
        nnz=16,
        circ_counts=[7, 9],
        top_circs=[("circX", 2, 8), ("circY", 1, 3)],
        workdir_size_bytes=150,
        align_size_bytes=30,
        h5ad_bytes=30,
    )

    result = subprocess.run(
        [
            sys.executable,
            "scripts/compare_full_length_workflow_runs.py",
            str(left),
            str(right),
            "--json",
        ],
        cwd=Path(__file__).resolve().parents[1],
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    payload = json.loads(result.stdout)
    assert payload["left"]["cells"] == 1
    assert payload["right"]["cells"] == 2
    assert payload["delta"]["cells"] == 1
    assert payload["delta"]["circRNAs"] == 2
    assert payload["delta"]["nnz"] == 4
    assert payload["shared_top_recurrent_circRNAs"] == ["circX"]
