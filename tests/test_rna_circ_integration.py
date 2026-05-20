from __future__ import annotations

import json
from pathlib import Path

from typer.testing import CliRunner

from circyto.cli.circyto import app


runner = CliRunner()


def _write_fixture_workdir(tmp_path: Path, *, include_rna_qc: bool) -> Path:
    root = tmp_path / "workflow"
    (root / "rna").mkdir(parents=True)
    (root / "matrix").mkdir(parents=True)
    (root / "qc").mkdir(parents=True)

    (root / "rna" / "gene_counts.tsv").write_text(
        "gene_id\tgene_name\tcellA\tDIYHEK_192\n"
        "G1\tMT-CO1\t10\t1\n"
        "G2\tRPLP0\t0\t5\n"
        "G3\tGENE3\t3\t2\n",
        encoding="utf-8",
    )
    (root / "rna" / "gene_feature_table.tsv").write_text(
        "gene_id\tgene_name\tgene_biotype\n"
        "G1\tMT-CO1\tprotein_coding\n"
        "G2\tRPLP0\tprotein_coding\n"
        "G3\tGENE3\tlncRNA\n",
        encoding="utf-8",
    )
    (root / "rna" / "rna_import_summary.json").write_text(
        json.dumps({"n_genes": 3, "n_cells": 2, "total_counts_sum": 21}, indent=2) + "\n",
        encoding="utf-8",
    )
    if include_rna_qc:
        (root / "qc" / "rna_qc.tsv").write_text(
            "cell_id\ttotal_counts\tdetected_genes\tmitochondrial_fraction\tribosomal_fraction\tcircRNA_count\n"
            "cellA\t13\t2\t0.7692307692\t0.0\t2\n"
            "DIYHEK_192\t8\t3\t0.125\t0.625\t0\n",
            encoding="utf-8",
        )
    (root / "matrix" / "circ_counts.mtx").write_text(
        "%%MatrixMarket matrix coordinate integer general\n%\n3 1 2\n1 1 2\n3 1 1\n",
        encoding="utf-8",
    )
    (root / "matrix" / "circ_index.txt").write_text("circ1\ncirc2\ncirc3\n", encoding="utf-8")
    (root / "matrix" / "cell_index.txt").write_text("cellA\n", encoding="utf-8")
    return root


def test_summarize_rna_circ_shared_and_rna_only_cells(tmp_path: Path) -> None:
    root = _write_fixture_workdir(tmp_path, include_rna_qc=True)
    result = runner.invoke(app, ["summarize-rna-circ", "--workdir", str(root), "--json"])
    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["n_rna_cells"] == 2
    assert payload["n_circ_cells"] == 1
    assert payload["n_shared_cells"] == 1
    assert payload["n_rna_only_cells"] == 1
    assert payload["rna_only_cells"] == ["DIYHEK_192"]
    preview = {item["cell_id"]: item for item in payload["joined_preview"]}
    assert preview["cellA"]["circRNA_count"] == 2
    assert preview["cellA"]["circRNA_total_support"] == 3
    assert preview["DIYHEK_192"]["circRNA_count"] == 0


def test_summarize_rna_circ_falls_back_when_rna_qc_missing(tmp_path: Path) -> None:
    root = _write_fixture_workdir(tmp_path, include_rna_qc=False)
    result = runner.invoke(app, ["summarize-rna-circ", "--workdir", str(root), "--json"])
    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["n_rna_cells"] == 2
    preview = {item["cell_id"]: item for item in payload["joined_preview"]}
    assert preview["cellA"]["total_rna_counts"] == 13
    assert preview["DIYHEK_192"]["detected_genes"] == 3


def test_summarize_rna_circ_write_summary_outputs(tmp_path: Path) -> None:
    root = _write_fixture_workdir(tmp_path, include_rna_qc=False)
    result = runner.invoke(app, ["summarize-rna-circ", "--workdir", str(root), "--write-summary"])
    assert result.exit_code == 0, result.output
    assert (root / "qc" / "rna_circ_cell_summary.tsv").exists()
    assert (root / "qc" / "rna_circ_summary.json").exists()
    table_text = (root / "qc" / "rna_circ_cell_summary.tsv").read_text(encoding="utf-8")
    assert "DIYHEK_192" in table_text
    assert "circRNA_total_support" in table_text


def test_summarize_rna_circ_text_output(tmp_path: Path) -> None:
    root = _write_fixture_workdir(tmp_path, include_rna_qc=True)
    result = runner.invoke(app, ["summarize-rna-circ", "--workdir", str(root)])
    assert result.exit_code == 0, result.output
    assert "RNA cells=2 circ cells=1 shared=1" in result.output
    assert "RNA-only=1 circ-only=0" in result.output
