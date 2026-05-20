from __future__ import annotations

import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

from circyto.cli.circyto import app
from circyto.pipeline.gene_expression_velocity import HAS_MUDATA


runner = CliRunner()


def _write_fixture_workdir(tmp_path: Path) -> Path:
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
        json.dumps({"method": "simple-overlap", "n_genes": 3, "n_cells": 2}, indent=2) + "\n",
        encoding="utf-8",
    )
    (root / "qc" / "rna_qc.tsv").write_text(
        "cell_id\ttotal_counts\tdetected_genes\tmitochondrial_fraction\tribosomal_fraction\tcircRNA_count\n"
        "cellA\t13\t2\t0.7692307692\t0.0\t2\n"
        "DIYHEK_192\t8\t3\t0.125\t0.625\t0\n",
        encoding="utf-8",
    )
    (root / "qc" / "rna_circ_cell_summary.tsv").write_text(
        "cell_id\tmembership\ttotal_rna_counts\tdetected_genes\tmitochondrial_fraction\tribosomal_fraction\tcircRNA_count\tcircRNA_total_support\n"
        "cellA\tshared\t13\t2\t0.7692307692\t0.0\t2\t3\n"
        "DIYHEK_192\trna_only\t8\t3\t0.125\t0.625\t0\t0\n",
        encoding="utf-8",
    )
    (root / "qc" / "rna_circ_summary.json").write_text(
        json.dumps({"n_shared_cells": 1, "n_rna_only_cells": 1, "n_circ_only_cells": 0}, indent=2) + "\n",
        encoding="utf-8",
    )
    (root / "matrix" / "circ_counts.mtx").write_text(
        "%%MatrixMarket matrix coordinate integer general\n%\n3 1 2\n1 1 2\n3 1 1\n",
        encoding="utf-8",
    )
    (root / "matrix" / "circ_index.txt").write_text("circ1\ncirc2\ncirc3\n", encoding="utf-8")
    (root / "matrix" / "cell_index.txt").write_text("cellA\n", encoding="utf-8")
    (root / "matrix" / "circ_feature_table.tsv").write_text(
        "circ_id\tchrom\tstart\tend\tstrand\thost_gene\n"
        "circ1\tchr1\t1\t10\t+\tGENE1\n"
        "circ2\tchr1\t11\t20\t+\tGENE2\n"
        "circ3\tchr1\t21\t30\t-\tGENE3\n",
        encoding="utf-8",
    )
    (root / "workflow_summary.json").write_text(
        json.dumps({"workflow": "full-length-circrna", "protocol": "ramda"}, indent=2) + "\n",
        encoding="utf-8",
    )
    return root


def test_export_mudata_missing_rna_profile_fails_clearly(tmp_path: Path) -> None:
    root = tmp_path / "workflow"
    root.mkdir()
    result = runner.invoke(app, ["export-mudata", "--workdir", str(root)])
    assert result.exit_code != 0
    assert "Missing RNA gene-count table" in result.output


def test_export_mudata_missing_circ_matrix_fails_clearly(tmp_path: Path) -> None:
    root = tmp_path / "workflow"
    (root / "rna").mkdir(parents=True)
    (root / "rna" / "gene_counts.tsv").write_text("gene_id\tgene_name\tcellA\nG1\tGENE1\t1\n", encoding="utf-8")
    (root / "rna" / "gene_feature_table.tsv").write_text("gene_id\tgene_name\nG1\tGENE1\n", encoding="utf-8")
    result = runner.invoke(app, ["export-mudata", "--workdir", str(root)])
    assert result.exit_code != 0
    assert "Missing circ matrix" in result.output


def test_export_mudata_creates_h5mu_and_preserves_metadata(tmp_path: Path) -> None:
    if not HAS_MUDATA:
        pytest.skip("mudata not installed")
    import mudata as mu
    import scipy.sparse as sp

    root = _write_fixture_workdir(tmp_path)
    out_path = tmp_path / "full_length.h5mu"
    result = runner.invoke(app, ["export-mudata", "--workdir", str(root), "--output", str(out_path)])
    assert result.exit_code == 0, result.output
    assert out_path.exists()
    payload = json.loads(result.output)
    assert payload["n_rna_only_cells"] == 1
    assert "DIYHEK_192" in payload["rna_only_cells"]

    mdata = mu.read_h5mu(out_path)
    assert "rna" in mdata.mod
    assert "circ" in mdata.mod
    assert list(mdata.obs_names) == ["cellA", "DIYHEK_192"]
    assert mdata.mod["rna"].shape == (2, 3)
    assert mdata.mod["circ"].shape == (2, 3)
    assert sp.issparse(mdata.mod["rna"].X)
    assert sp.issparse(mdata.mod["circ"].X)
    assert int(mdata.mod["circ"][1, :].X.sum()) == 0
    assert mdata.obs.loc["DIYHEK_192", "membership"] == "rna_only"
    assert int(mdata.obs.loc["DIYHEK_192", "circRNA_count"]) == 0
    assert "circyto" in mdata.uns
    assert mdata.uns["circyto"]["command_name"] == "circyto export-mudata"
    assert "workflow_summary" in mdata.uns["circyto"]


def test_export_mudata_transposes_circ_by_cell_matrix_to_cells_by_circ(tmp_path: Path) -> None:
    if not HAS_MUDATA:
        pytest.skip("mudata not installed")
    import mudata as mu

    root = tmp_path / "workflow"
    (root / "rna").mkdir(parents=True)
    (root / "matrix").mkdir(parents=True)
    (root / "qc").mkdir(parents=True)

    rna_cells = [f"cell{i:03d}" for i in range(191)] + ["DIYHEK_192"]
    (root / "rna" / "gene_counts.tsv").write_text(
        "gene_id\tgene_name\t" + "\t".join(rna_cells) + "\n"
        + "G1\tGENE1\t" + "\t".join(["1"] * 191 + ["2"]) + "\n",
        encoding="utf-8",
    )
    (root / "rna" / "gene_feature_table.tsv").write_text(
        "gene_id\tgene_name\tgene_biotype\nG1\tGENE1\tprotein_coding\n",
        encoding="utf-8",
    )
    qc_lines = ["cell_id\ttotal_counts\tdetected_genes\tmitochondrial_fraction\tribosomal_fraction\tcircRNA_count"]
    qc_lines.extend(f"{cell_id}\t1\t1\t0.0\t0.0\t1" for cell_id in rna_cells[:-1])
    qc_lines.append("DIYHEK_192\t2\t1\t0.0\t0.0\t0")
    (root / "qc" / "rna_qc.tsv").write_text("\n".join(qc_lines) + "\n", encoding="utf-8")
    (root / "matrix" / "circ_counts.mtx").write_text(
        "%%MatrixMarket matrix coordinate integer general\n%\n2503 191 2\n1 1 1\n2503 191 2\n",
        encoding="utf-8",
    )
    circ_ids = [f"circ{i:04d}" for i in range(2503)]
    circ_cells = [f"cell{i:03d}" for i in range(191)]
    (root / "matrix" / "circ_index.txt").write_text("\n".join(circ_ids) + "\n", encoding="utf-8")
    (root / "matrix" / "cell_index.txt").write_text("\n".join(circ_cells) + "\n", encoding="utf-8")
    circ_feature_lines = ["circ_id\tchrom\tstart\tend\tstrand\thost_gene"]
    circ_feature_lines.extend(
        f"{circ_id}\tchr1\t1\t2\t+\tGENE{i:04d}" for i, circ_id in enumerate(circ_ids, start=1)
    )
    (root / "matrix" / "circ_feature_table.tsv").write_text("\n".join(circ_feature_lines) + "\n", encoding="utf-8")
    out_path = tmp_path / "full_length_large_oriented.h5mu"

    result = runner.invoke(app, ["export-mudata", "--workdir", str(root), "--output", str(out_path)])
    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["n_obs"] == 192
    assert payload["n_circ_features"] == 2503
    assert payload["n_rna_only_cells"] == 1

    mdata = mu.read_h5mu(out_path)
    assert mdata.mod["circ"].shape == (192, 2503)
    assert mdata.mod["rna"].shape == (192, 1)
    assert int(mdata.mod["circ"][191, :].X.sum()) == 0


def test_export_mudata_no_overwrite_without_flag(tmp_path: Path) -> None:
    if not HAS_MUDATA:
        pytest.skip("mudata not installed")
    root = _write_fixture_workdir(tmp_path)
    out_path = tmp_path / "full_length.h5mu"
    out_path.write_text("existing", encoding="utf-8")
    result = runner.invoke(app, ["export-mudata", "--workdir", str(root), "--output", str(out_path)])
    assert result.exit_code != 0
    assert "Use --overwrite to replace it" in result.output


def test_export_mudata_missing_package_gives_clear_error(tmp_path: Path) -> None:
    if HAS_MUDATA:
        pytest.skip("mudata installed")
    root = _write_fixture_workdir(tmp_path)
    result = runner.invoke(app, ["export-mudata", "--workdir", str(root)])
    assert result.exit_code != 0
    assert "mudata is not installed" in result.output
    assert "circyto[mudata]" in result.output
