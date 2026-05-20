from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


def _write_fixture_workdir(tmp_path: Path) -> Path:
    root = tmp_path / "workflow"
    (root / "rna").mkdir(parents=True)
    (root / "matrix").mkdir(parents=True)

    (root / "rna" / "gene_counts.tsv").write_text(
        "gene_id\tgene_name\tcellA\tcellB\n"
        "G1\tGENE1\t10\t1\n"
        "G2\tGENE2\t0\t5\n"
        "G3\tGENE3\t3\t2\n",
        encoding="utf-8",
    )
    (root / "rna" / "gene_feature_table.tsv").write_text(
        "gene_id\tgene_name\tchrom\tstart\tend\tstrand\n"
        "G1\tGENE1\tchr21\t1\t10\t+\n"
        "G2\tGENE2\tchr21\t20\t30\t+\n"
        "G3\tGENE3\tchr21\t40\t50\t-\n",
        encoding="utf-8",
    )
    (root / "rna" / "rna_import_summary.json").write_text(
        json.dumps(
            {
                "n_genes": 3,
                "n_cells": 2,
                "total_counts_sum": 21,
                "assigned_templates": 20,
                "ambiguous_templates_excluded": 2,
                "unassigned_templates": 5,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    (root / "matrix" / "cell_index.txt").write_text("cellA\ncellB\n", encoding="utf-8")
    return root


def test_check_rna_profile_outputs_text(tmp_path: Path) -> None:
    root = _write_fixture_workdir(tmp_path)
    result = subprocess.run(
        [sys.executable, "scripts/check_rna_profile_outputs.py", "--workdir", str(root)],
        cwd=Path(__file__).resolve().parents[1],
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert "RNA outputs: gene_counts=True gene_feature_table=True rna_import_summary=True" in result.stdout
    assert "Counts: n_genes=3 n_cells=2 total_counts_sum=21" in result.stdout
    assert "G1 GENE1: total_count=11" in result.stdout
    assert "cellA: total_count=13" in result.stdout
    assert "cellB: total_count=8" in result.stdout


def test_check_rna_profile_outputs_json_and_cell_id_match(tmp_path: Path) -> None:
    root = _write_fixture_workdir(tmp_path)
    result = subprocess.run(
        [sys.executable, "scripts/check_rna_profile_outputs.py", "--workdir", str(root), "--json"],
        cwd=Path(__file__).resolve().parents[1],
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    payload = json.loads(result.stdout)
    assert payload["n_genes"] == 3
    assert payload["n_cells"] == 2
    assert payload["assigned_templates"] == 20
    assert payload["cell_ids_match_matrix"]["matches"] is True
    assert payload["top_expressed_genes"][0]["gene_id"] == "G1"
