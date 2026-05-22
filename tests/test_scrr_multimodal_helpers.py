from __future__ import annotations

import json
from pathlib import Path

from typer.testing import CliRunner

from circyto.cli.circyto import app


runner = CliRunner()


def _write_multimodal_fixture(tmp_path: Path) -> Path:
    root = tmp_path / "workflow"
    (root / "matrix").mkdir(parents=True)
    (root / "rna").mkdir()
    (root / "qc").mkdir()
    (root / "mudata").mkdir()

    (root / "workflow_summary.json").write_text(
        json.dumps(
            {
                "workflow_type": "full-length-circrna",
                "protocol": "ramda",
                "read_layout": "single",
                "workflow_uuid": "uuid-1",
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    (root / "rna" / "gene_counts.tsv").write_text(
        "gene_id\tgene_name\tcellA\nG1\tGENE1\t5\n",
        encoding="utf-8",
    )
    (root / "rna" / "gene_feature_table.tsv").write_text(
        "gene_id\tgene_name\nG1\tGENE1\n",
        encoding="utf-8",
    )
    (root / "qc" / "rna_qc.tsv").write_text(
        "cell_id\ttotal_counts\tdetected_genes\tmitochondrial_fraction\tribosomal_fraction\tcircRNA_count\ncellA\t5\t1\t0.0\t0.0\t2\n",
        encoding="utf-8",
    )
    (root / "matrix" / "circ_counts.mtx").write_text(
        "%%MatrixMarket matrix coordinate integer general\n%\n2 1 2\n1 1 3\n2 1 1\n",
        encoding="utf-8",
    )
    (root / "matrix" / "circ_index.txt").write_text("circA\ncircB\n", encoding="utf-8")
    (root / "matrix" / "cell_index.txt").write_text("cellA\n", encoding="utf-8")
    (root / "matrix" / "circ_feature_table.tsv").write_text(
        "circ_id\tchrom\tstart\tend\tstrand\thost_gene\n"
        "circA\tchr1\t100\t200\t+\tGENE1\n"
        "circB\tchr1\t300\t400\t-\tGENE1\n",
        encoding="utf-8",
    )
    (root / "mudata" / "full_length.h5mu").write_text("", encoding="utf-8")
    return root


def test_inspect_workdir_reports_modalities(tmp_path: Path) -> None:
    root = _write_multimodal_fixture(tmp_path)
    result = runner.invoke(app, ["inspect-workdir", "--workdir", str(root), "--json"])
    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["available_modalities"] == ["rna", "circ"]
    assert payload["source_workflow_type"] == "full-length-circrna"
    assert payload["source_protocol"] == "ramda"
    assert payload["mudata_present"] is True
    assert payload["qc_present"] is True


def test_summarize_circ_host_genes_writes_summary(tmp_path: Path) -> None:
    root = _write_multimodal_fixture(tmp_path)
    result = runner.invoke(app, ["summarize-circ-host-genes", "--workdir", str(root), "--json"])
    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["n_host_genes"] == 1
    assert payload["n_circ_with_host_gene"] == 2
    table_path = Path(payload["output_path"])
    assert table_path.exists()
    table_text = table_path.read_text(encoding="utf-8")
    assert "GENE1" in table_text
    assert "n_circRNAs" in table_text


def test_export_circ_bed_writes_bed_like_rows(tmp_path: Path) -> None:
    root = _write_multimodal_fixture(tmp_path)
    out = root / "qc" / "custom_circs.bed"
    result = runner.invoke(app, ["export-circ-bed", "--workdir", str(root), "--output", str(out)])
    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["n_rows"] == 2
    bed_lines = out.read_text(encoding="utf-8").splitlines()
    assert bed_lines[0] == "chr1\t100\t200\tcircA\t3"
    assert bed_lines[1] == "chr1\t300\t400\tcircB\t1"
