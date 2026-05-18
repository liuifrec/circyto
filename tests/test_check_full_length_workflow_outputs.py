from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


def _write_fixture_workflow(tmp_path: Path) -> Path:
    root = tmp_path / "workflow"
    (root / "ciri3").mkdir(parents=True)
    (root / "matrix").mkdir(parents=True)
    (root / "logs").mkdir(parents=True)
    (root / "anndata").mkdir(parents=True)

    (root / "workflow_summary.json").write_text(
        json.dumps(
            {
                "workflow": "full-length-circrna",
                "protocol": "ramda",
                "dry_run": False,
                "planned_cells": 2,
                "completed_cells": 2,
                "stage_graph": [
                    {"stage": "alignment", "status": "completed"},
                    {"stage": "detector", "status": "completed"},
                    {"stage": "matrix", "status": "completed"},
                    {"stage": "h5ad_export", "status": "completed"},
                ],
                "matrix": {
                    "n_rows": 221,
                    "n_cols": 2,
                    "nnz": 250,
                    "n_cells": 2,
                    "n_circRNAs": 221,
                    "top_recurrent_circRNAs": [
                        {"circ_id": "circA", "n_cells_detected": 2, "total_support": 9},
                        {"circ_id": "circB", "n_cells_detected": 2, "total_support": 5},
                    ],
                },
                "paths": {"h5ad": str((root / "anndata" / "circ_counts.h5ad").resolve())},
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    (root / "ciri3" / "detector_run_summary.json").write_text(
        json.dumps(
            {
                "failed_cells": 0,
                "status_counts": {"success": 2},
                "cells": [
                    {
                        "cell_id": "IMR90_A_100",
                        "status": "success",
                        "raw_row_count": 120,
                        "normalized_row_count": 110,
                    },
                    {
                        "cell_id": "IMR90_A_101",
                        "status": "success",
                        "raw_row_count": 118,
                        "normalized_row_count": 111,
                    },
                ],
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    (root / "matrix" / "circ_counts.mtx").write_text(
        "%%MatrixMarket matrix coordinate integer general\n%\n221 2 250\n1 1 4\n",
        encoding="utf-8",
    )
    (root / "anndata" / "circ_counts.h5ad").write_text("placeholder", encoding="utf-8")
    (root / "logs" / "workflow.log").write_text(
        "INFO completed\nERROR sample one had a transient warning path\n",
        encoding="utf-8",
    )
    (root / "ciri3" / "big_output.tsv").write_text("x" * 4096, encoding="utf-8")
    return root


def test_check_full_length_workflow_outputs_text(tmp_path: Path) -> None:
    root = _write_fixture_workflow(tmp_path)
    result = subprocess.run(
        [sys.executable, "scripts/check_full_length_workflow_outputs.py", str(root)],
        cwd=Path(__file__).resolve().parents[1],
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert "Stage graph:" in result.stdout
    assert "Cells: planned=2 completed=2 failed=0" in result.stdout
    assert "Matrix: n_rows=221 n_cols=2 nnz=250 n_cells=2 n_circRNAs=221" in result.stdout
    assert "h5ad exists: True" in result.stdout
    assert "circA" in result.stdout
    assert "workflow.log: ERROR sample one had a transient warning path" in result.stdout


def test_check_full_length_workflow_outputs_json_fallbacks_to_matrix_header(tmp_path: Path) -> None:
    root = _write_fixture_workflow(tmp_path)
    workflow_summary = json.loads((root / "workflow_summary.json").read_text(encoding="utf-8"))
    workflow_summary["matrix"] = {}
    (root / "workflow_summary.json").write_text(json.dumps(workflow_summary, indent=2) + "\n", encoding="utf-8")

    result = subprocess.run(
        [sys.executable, "scripts/check_full_length_workflow_outputs.py", str(root), "--json"],
        cwd=Path(__file__).resolve().parents[1],
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    payload = json.loads(result.stdout)
    assert payload["matrix"]["n_rows"] == 221
    assert payload["matrix"]["n_cols"] == 2
    assert payload["matrix"]["nnz"] == 250
    assert payload["largest_files"][0]["path"] == "ciri3/big_output.tsv"
