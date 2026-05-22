from __future__ import annotations

import json
from pathlib import Path

from typer.testing import CliRunner

from circyto.cli.circyto import app


runner = CliRunner()


def _write_base_workdir(tmp_path: Path) -> Path:
    root = tmp_path / "workflow"
    (root / "rna").mkdir(parents=True)
    (root / "matrix").mkdir()
    (root / "qc").mkdir()
    (root / "rna" / "gene_counts.tsv").write_text(
        "gene_id\tgene_name\tcellA\tcellB\nG1\tGENE1\t5\t2\n",
        encoding="utf-8",
    )
    (root / "rna" / "gene_feature_table.tsv").write_text(
        "gene_id\tgene_name\nG1\tGENE1\n",
        encoding="utf-8",
    )
    (root / "qc" / "rna_qc.tsv").write_text(
        "cell_id\ttotal_counts\tdetected_genes\tmitochondrial_fraction\tribosomal_fraction\tcircRNA_count\ncellA\t5\t1\t0.0\t0.0\t2\ncellB\t2\t1\t0.0\t0.0\t0\n",
        encoding="utf-8",
    )
    (root / "matrix" / "circ_counts.mtx").write_text(
        "%%MatrixMarket matrix coordinate integer general\n%\n2 1 1\n1 1 2\n",
        encoding="utf-8",
    )
    (root / "matrix" / "circ_index.txt").write_text("circ1\ncirc2\n", encoding="utf-8")
    (root / "matrix" / "cell_index.txt").write_text("cellA\n", encoding="utf-8")
    return root


def test_import_dna_snv_summary_and_join(tmp_path: Path) -> None:
    root = _write_base_workdir(tmp_path)
    dna_cell = tmp_path / "dna_cell_summary.tsv"
    dna_cell.write_text(
        "cell_id\tdna_library_id\tcnv_burden\treplication_score\tcell_cycle_phase\tdna_variant_count\tnotes\n"
        "cellA\tlibA\t1.5\t0.7\tS\t3\tok\n"
        "cellC\tlibC\t2.0\t0.2\tG1\t1\tDNA-only\n",
        encoding="utf-8",
    )
    scomatic = tmp_path / "scomatic_candidate_summary.tsv"
    scomatic.write_text(
        "variant_id\tcell_id\tchrom\tpos\tref\talt\tgene\tfilter_status\tcandidate_variant_class\tread_support\tvaf\tcaller\n"
        "var1\tcellA\tchr1\t10\tA\tG\tGENE1\tPASS\tRNA-derived candidate variant signal\t4\t0.2\tSComatic\n",
        encoding="utf-8",
    )
    result = runner.invoke(
        app,
        [
            "import-dna-snv-summary",
            "--workdir",
            str(root),
            "--dna-cell-summary",
            str(dna_cell),
            "--scomatic-candidate-summary",
            str(scomatic),
        ],
    )
    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["n_shared_cells"] == 1
    assert payload["n_dna_only_cells"] == 1
    assert "RNA-derived candidate variant signals" in payload["terminology_note"]

    join_result = runner.invoke(app, ["summarize-dna-rna-circ", "--workdir", str(root), "--write-summary"])
    assert join_result.exit_code == 0, join_result.output
    join_payload = json.loads(join_result.output)
    assert join_payload["n_shared_cells"] == 1
    assert join_payload["n_dna_only_cells"] == 1
    table_text = (root / "qc" / "dna_rna_circ_cell_summary.tsv").read_text(encoding="utf-8")
    assert "scomatic_candidate_count" in table_text
    assert "cellC" in table_text


def test_import_dna_snv_summary_fails_on_missing_required_columns(tmp_path: Path) -> None:
    root = _write_base_workdir(tmp_path)
    dna_cell = tmp_path / "dna_cell_summary.tsv"
    dna_cell.write_text("cell_id\tcnv_burden\ncellA\t1.0\n", encoding="utf-8")
    result = runner.invoke(
        app,
        ["import-dna-snv-summary", "--workdir", str(root), "--dna-cell-summary", str(dna_cell)],
    )
    assert result.exit_code != 0
    assert "missing required columns" in result.output.lower()


def test_dna_snv_docs_preserve_candidate_signal_terminology() -> None:
    text = Path("docs/dna_snv_integration.md").read_text(encoding="utf-8")
    assert "RNA-derived candidate variant signals" in text
