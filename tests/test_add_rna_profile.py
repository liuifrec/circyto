from __future__ import annotations

import json
from pathlib import Path

from typer.testing import CliRunner

from circyto.cli.circyto import app


runner = CliRunner()


def _write_gtf(path: Path) -> None:
    path.write_text(
        'chr21\ttest\texon\t10\t20\t.\t+\t.\tgene_id "GENE1"; gene_name "GENE1";\n'
        'chr21\ttest\texon\t40\t50\t.\t+\t.\tgene_id "GENE2"; gene_name "GENE2";\n',
        encoding="utf-8",
    )


def _write_sam(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "@SQ\tSN:chr21\tLN:1000\n"
        "read1\t0\tchr21\t10\t255\t8M\t*\t0\t0\tACGTACGT\t!!!!!!!!\n"
        "read2\t0\tchr21\t40\t255\t8M\t*\t0\t0\tACGTACGT\t!!!!!!!!\n",
        encoding="utf-8",
    )


def _write_alignment_manifest(path: Path, sam_path: Path, *, cell_id: str = "cellA") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "cell_id\tbam\tsam\tgroup_id\tread_layout\taligner\treference\tcache_key\tsource_manifest\tmapper_mode\tartifact_bucket\tsortedness\n"
        f"{cell_id}\t\t{sam_path}\tgrp\tsingle-end\tbwa-mem\t/tmp/ref.fa\tk1\t/tmp/source.tsv\t0\tbwa_mem\tunsorted\n",
        encoding="utf-8",
    )


def test_add_rna_profile_detects_manifest_in_align_root(tmp_path: Path) -> None:
    workdir = tmp_path / "completed"
    sam = workdir / "align" / "cache" / "alignment.sam"
    _write_sam(sam)
    _write_alignment_manifest(workdir / "align" / "alignment_manifest.tsv", sam)
    gtf = tmp_path / "genes.gtf"
    _write_gtf(gtf)

    result = runner.invoke(
        app,
        ["add-rna-profile", "--workdir", str(workdir), "--gtf", str(gtf), "--dry-run"],
    )
    assert result.exit_code == 0, result.stdout
    payload = json.loads(result.stdout)
    assert payload["command_name"] == "circyto add-rna-profile"
    assert payload["dry_run"] is True
    assert payload["alignment_manifest"].endswith("align/alignment_manifest.tsv")


def test_add_rna_profile_detects_manifest_in_align_star(tmp_path: Path) -> None:
    workdir = tmp_path / "completed"
    sam = workdir / "align" / "star" / "cache" / "alignment.sam"
    _write_sam(sam)
    _write_alignment_manifest(workdir / "align" / "star" / "alignment_manifest.tsv", sam)
    gtf = tmp_path / "genes.gtf"
    _write_gtf(gtf)

    result = runner.invoke(
        app,
        ["add-rna-profile", "--workdir", str(workdir), "--gtf", str(gtf), "--dry-run"],
    )
    assert result.exit_code == 0, result.stdout
    payload = json.loads(result.stdout)
    assert payload["alignment_manifest"].endswith("align/star/alignment_manifest.tsv")


def test_add_rna_profile_dry_run_does_not_write_outputs(tmp_path: Path) -> None:
    workdir = tmp_path / "completed"
    sam = workdir / "align" / "cache" / "alignment.sam"
    _write_sam(sam)
    _write_alignment_manifest(workdir / "align" / "alignment_manifest.tsv", sam)
    gtf = tmp_path / "genes.gtf"
    _write_gtf(gtf)

    result = runner.invoke(
        app,
        ["add-rna-profile", "--workdir", str(workdir), "--gtf", str(gtf), "--dry-run"],
    )
    assert result.exit_code == 0, result.stdout
    assert not (workdir / "rna" / "gene_counts.tsv").exists()
    assert not (workdir / "rna" / "rna_import_summary.json").exists()


def test_add_rna_profile_real_run_writes_rna_outputs_and_updates_summary(tmp_path: Path) -> None:
    workdir = tmp_path / "completed"
    sam = workdir / "align" / "cache" / "alignment.sam"
    _write_sam(sam)
    _write_alignment_manifest(workdir / "align" / "alignment_manifest.tsv", sam, cell_id="cellA")
    (workdir / "matrix").mkdir(parents=True, exist_ok=True)
    (workdir / "matrix" / "circ_counts.mtx").write_text(
        "%%MatrixMarket matrix coordinate integer general\n%\n2 1 1\n1 1 1\n",
        encoding="utf-8",
    )
    (workdir / "matrix" / "circ_index.txt").write_text("circ1\ncirc2\n", encoding="utf-8")
    (workdir / "matrix" / "cell_index.txt").write_text("cellA\n", encoding="utf-8")
    (workdir / "workflow_summary.json").write_text(
        json.dumps({"workflow": "full-length-circrna", "existing": {"keep": 1}}, indent=2) + "\n",
        encoding="utf-8",
    )
    gtf = tmp_path / "genes.gtf"
    _write_gtf(gtf)

    result = runner.invoke(
        app,
        ["add-rna-profile", "--workdir", str(workdir), "--gtf", str(gtf)],
    )
    assert result.exit_code == 0, result.stdout
    assert (workdir / "rna" / "gene_counts.tsv").exists()
    assert (workdir / "rna" / "gene_feature_table.tsv").exists()
    assert (workdir / "rna" / "rna_import_summary.json").exists()
    assert (workdir / "qc" / "rna_qc.tsv").exists()
    assert (workdir / "qc" / "rna_gene_qc.tsv").exists()

    payload = json.loads(result.stdout)
    assert payload["dry_run"] is False
    assert payload["rna_import"]["method"] == "simple-overlap"
    assert payload["rna_import"]["assigned_templates"] == 2
    assert payload["rna_import"]["genes_detected_ge_1_cell"] == 2
    assert payload["rna_import"]["median_genes_per_cell"] == 2.0

    summary = json.loads((workdir / "workflow_summary.json").read_text(encoding="utf-8"))
    assert summary["command_name"] == "circyto add-rna-profile"
    assert summary["existing"] == {"keep": 1}
    assert summary["rna_import"]["n_cells"] == 1
    assert summary["rna_import"]["n_genes"] == 2
    qc_text = (workdir / "qc" / "rna_qc.tsv").read_text(encoding="utf-8")
    assert "circRNA_count" in qc_text


def test_add_rna_profile_missing_manifest_gives_clear_error(tmp_path: Path) -> None:
    workdir = tmp_path / "completed"
    workdir.mkdir()
    gtf = tmp_path / "genes.gtf"
    _write_gtf(gtf)

    result = runner.invoke(
        app,
        ["add-rna-profile", "--workdir", str(workdir), "--gtf", str(gtf)],
    )
    assert result.exit_code != 0
    assert "Could not find alignment manifest" in result.output


def test_refresh_rna_qc_succeeds_from_existing_rna_outputs_without_alignments(tmp_path: Path) -> None:
    workdir = tmp_path / "completed"
    (workdir / "rna").mkdir(parents=True, exist_ok=True)
    (workdir / "rna" / "gene_counts.tsv").write_text(
        "gene_id\tgene_name\tcellA\tDIYHEK_192\n"
        "GENE1\tMT-CO1\t2\t0\n"
        "GENE2\tRPLP0\t1\t3\n",
        encoding="utf-8",
    )
    (workdir / "rna" / "gene_feature_table.tsv").write_text(
        "gene_id\tgene_name\tgene_biotype\n"
        "GENE1\tMT-CO1\tprotein_coding\n"
        "GENE2\tRPLP0\tprotein_coding\n",
        encoding="utf-8",
    )
    (workdir / "rna" / "rna_import_summary.json").write_text(
        json.dumps({"method": "simple-overlap", "assigned_templates": 3}, indent=2) + "\n",
        encoding="utf-8",
    )
    (workdir / "matrix").mkdir(parents=True, exist_ok=True)
    (workdir / "matrix" / "circ_counts.mtx").write_text(
        "%%MatrixMarket matrix coordinate integer general\n%\n2 1 1\n1 1 4\n",
        encoding="utf-8",
    )
    (workdir / "matrix" / "circ_index.txt").write_text("circ1\ncirc2\n", encoding="utf-8")
    (workdir / "matrix" / "cell_index.txt").write_text("cellA\n", encoding="utf-8")
    (workdir / "workflow_summary.json").write_text(
        json.dumps({"workflow": "full-length-circrna", "rna_import": {"method": "simple-overlap"}}, indent=2) + "\n",
        encoding="utf-8",
    )

    result = runner.invoke(app, ["refresh-rna-qc", "--workdir", str(workdir)])
    assert result.exit_code == 0, result.stdout
    payload = json.loads(result.stdout)
    assert payload["command_name"] == "circyto refresh-rna-qc"
    assert (workdir / "qc" / "rna_qc.tsv").exists()
    assert (workdir / "qc" / "rna_gene_qc.tsv").exists()
    qc_text = (workdir / "qc" / "rna_qc.tsv").read_text(encoding="utf-8")
    assert "DIYHEK_192" in qc_text
    assert "\t0\n" in qc_text
    refreshed = json.loads((workdir / "rna" / "rna_import_summary.json").read_text(encoding="utf-8"))
    assert refreshed["genes_detected_ge_1_cell"] == 2
    assert refreshed["genes_detected_ge_3_cells"] == 0
    assert refreshed["median_genes_per_cell"] == 1.5
    assert refreshed["matrix_rna_only_cell_count"] == 1
    summary = json.loads((workdir / "workflow_summary.json").read_text(encoding="utf-8"))
    assert summary["command_name"] == "circyto refresh-rna-qc"
    assert summary["rna_import"]["output_rna_qc"].endswith("qc/rna_qc.tsv")


def test_refresh_rna_qc_missing_rna_outputs_fails_clearly(tmp_path: Path) -> None:
    workdir = tmp_path / "completed"
    workdir.mkdir()
    result = runner.invoke(app, ["refresh-rna-qc", "--workdir", str(workdir)])
    assert result.exit_code != 0
    assert "Missing RNA gene-count table" in result.output
