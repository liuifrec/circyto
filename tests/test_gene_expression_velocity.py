from __future__ import annotations

from pathlib import Path

import pytest

from circyto.pipeline.gene_expression_velocity import (
    build_cleanup_plan,
    count_gene_expression_from_alignments,
    execute_cleanup_plan,
    import_gene_counts_table,
    validate_cell_id_consistency,
    validate_feature_id_uniqueness,
    validate_gene_expression_table_schema,
    validate_velocity_layers_schema,
)


def test_validate_gene_expression_table_schema_passes_for_valid_tsv(tmp_path: Path) -> None:
    path = tmp_path / "gene_counts.tsv"
    path.write_text(
        "gene_id\tgene_name\tcellA\tcellB\n"
        "ENSG1\tGENE1\t1\t0\n"
        "ENSG2\tGENE2\t0\t2\n",
        encoding="utf-8",
    )
    summary = validate_gene_expression_table_schema(path)
    assert summary["n_genes"] == 2
    assert summary["n_cells"] == 2
    assert summary["cell_ids"] == ["cellA", "cellB"]


def test_validate_gene_expression_table_schema_fails_without_gene_name(tmp_path: Path) -> None:
    path = tmp_path / "gene_counts.tsv"
    path.write_text(
        "gene_id\tcellA\tcellB\n"
        "ENSG1\t1\t0\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="gene_name"):
        validate_gene_expression_table_schema(path)


def test_validate_cell_id_consistency_fails_on_mismatch() -> None:
    with pytest.raises(ValueError, match="Cell ID mismatch"):
        validate_cell_id_consistency(["cellA", "cellB"], ["cellA", "cellC"])


def test_validate_feature_id_uniqueness_fails_on_duplicate_ids() -> None:
    with pytest.raises(ValueError, match="Duplicate gene IDs detected"):
        validate_feature_id_uniqueness(["ENSG1", "ENSG1"], label="gene")


def test_validate_velocity_layers_schema_passes_for_minimal_layout(tmp_path: Path) -> None:
    velocity_dir = tmp_path / "velocity_layers"
    velocity_dir.mkdir()
    (velocity_dir / "barcodes.tsv").write_text("cellA\ncellB\n", encoding="utf-8")
    (velocity_dir / "features.tsv").write_text("ENSG1\tGENE1\nENSG2\tGENE2\n", encoding="utf-8")
    (velocity_dir / "spliced.mtx").write_text("%%MatrixMarket matrix coordinate integer general\n", encoding="utf-8")
    (velocity_dir / "unspliced.mtx").write_text("%%MatrixMarket matrix coordinate integer general\n", encoding="utf-8")
    summary = validate_velocity_layers_schema(velocity_dir)
    assert summary["n_cells"] == 2
    assert summary["n_genes"] == 2
    assert summary["has_ambiguous"] is False


def test_validate_velocity_layers_schema_fails_on_duplicate_feature_ids(tmp_path: Path) -> None:
    velocity_dir = tmp_path / "velocity_layers"
    velocity_dir.mkdir()
    (velocity_dir / "barcodes.tsv").write_text("cellA\ncellB\n", encoding="utf-8")
    (velocity_dir / "features.tsv").write_text("ENSG1\tGENE1\nENSG1\tGENE1_dup\n", encoding="utf-8")
    (velocity_dir / "spliced.mtx").write_text("%%MatrixMarket matrix coordinate integer general\n", encoding="utf-8")
    (velocity_dir / "unspliced.mtx").write_text("%%MatrixMarket matrix coordinate integer general\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Duplicate gene IDs detected"):
        validate_velocity_layers_schema(velocity_dir)


def test_import_gene_counts_table_writes_normalized_snapshot(tmp_path: Path) -> None:
    path = tmp_path / "gene_counts.tsv"
    path.write_text(
        "gene_id\tgene_name\tcellA\tcellB\n"
        "ENSG1\tGENE1\t1\t0\n"
        "ENSG2\tGENE2\t0\t2\n",
        encoding="utf-8",
    )
    outdir = tmp_path / "rna"
    summary = import_gene_counts_table(path=path, expected_cell_ids=["cellA", "cellB"], outdir=outdir)
    assert (outdir / "gene_counts.tsv").exists()
    assert (outdir / "gene_feature_table.tsv").exists()
    assert (outdir / "rna_import_summary.json").exists()
    assert summary["n_genes"] == 2
    assert summary["n_cells"] == 2


def test_import_gene_counts_table_fails_on_cell_id_mismatch(tmp_path: Path) -> None:
    path = tmp_path / "gene_counts.tsv"
    path.write_text(
        "gene_id\tgene_name\tcellA\tcellC\n"
        "ENSG1\tGENE1\t1\t0\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="Cell ID mismatch"):
        import_gene_counts_table(path=path, expected_cell_ids=["cellA", "cellB"], outdir=tmp_path / "rna")


def _write_gtf(path: Path) -> None:
    path.write_text(
        'chr21\ttest\tgene\t10\t30\t.\t+\t.\tgene_id "ENSG1"; gene_name "GENE1";\n'
        'chr21\ttest\tgene\t40\t60\t.\t+\t.\tgene_id "ENSG2"; gene_name "GENE2";\n'
        'chr21\ttest\tgene\t50\t70\t.\t+\t.\tgene_id "ENSG3"; gene_name "GENE3";\n',
        encoding="utf-8",
    )


def _write_alignment_manifest(path: Path, sam_path: Path, *, cell_id: str = "cellA", read_layout: str = "single-end") -> None:
    path.write_text(
        "cell_id\tbam\tsam\tgroup_id\tread_layout\taligner\treference\tcache_key\tsource_manifest\tmapper_mode\tartifact_bucket\tsortedness\n"
        f"{cell_id}\t\t{sam_path}\tgrp\t{read_layout}\tbwa-mem\t/tmp/ref.fa\tk1\t/tmp/source.tsv\t0\tbwa_mem\tunsorted\n",
        encoding="utf-8",
    )


def test_count_gene_expression_from_alignments_counts_unique_gene_overlaps(tmp_path: Path) -> None:
    gtf = tmp_path / "genes.gtf"
    _write_gtf(gtf)
    sam = tmp_path / "cellA.sam"
    sam.write_text(
        "@SQ\tSN:chr21\tLN:1000\n"
        "read1\t0\tchr21\t12\t255\t10M\t*\t0\t0\tACGTACGTAC\t!!!!!!!!!!\n"
        "read2\t0\tchr21\t40\t255\t10M\t*\t0\t0\tACGTACGTAC\t!!!!!!!!!!\n",
        encoding="utf-8",
    )
    manifest = tmp_path / "alignment_manifest.tsv"
    _write_alignment_manifest(manifest, sam)

    summary = count_gene_expression_from_alignments(
        alignment_manifest_path=manifest,
        gtf_path=gtf,
        expected_cell_ids=["cellA"],
        outdir=tmp_path / "rna",
    )

    counts = (tmp_path / "rna" / "gene_counts.tsv").read_text(encoding="utf-8")
    assert "ENSG1\tGENE1\t1" in counts
    assert "ENSG2\tGENE2\t1" in counts
    assert summary["assigned_templates"] == 2
    assert summary["ambiguous_templates_excluded"] == 0


def test_count_gene_expression_from_alignments_excludes_ambiguous_overlaps(tmp_path: Path) -> None:
    gtf = tmp_path / "genes.gtf"
    _write_gtf(gtf)
    sam = tmp_path / "cellA.sam"
    sam.write_text(
        "@SQ\tSN:chr21\tLN:1000\n"
        "ambig1\t0\tchr21\t52\t255\t10M\t*\t0\t0\tACGTACGTAC\t!!!!!!!!!!\n",
        encoding="utf-8",
    )
    manifest = tmp_path / "alignment_manifest.tsv"
    _write_alignment_manifest(manifest, sam)

    summary = count_gene_expression_from_alignments(
        alignment_manifest_path=manifest,
        gtf_path=gtf,
        expected_cell_ids=["cellA"],
        outdir=tmp_path / "rna",
    )

    counts = (tmp_path / "rna" / "gene_counts.tsv").read_text(encoding="utf-8")
    assert "ENSG2\tGENE2\t0" in counts
    assert "ENSG3\tGENE3\t0" in counts
    assert summary["assigned_templates"] == 0
    assert summary["ambiguous_templates_excluded"] == 1
    assert (tmp_path / "rna" / "rna_import_summary.json").exists()


def test_build_cleanup_plan_distinguishes_user_raw_fastq_from_generated_demux_fastq(tmp_path: Path) -> None:
    outdir = tmp_path / "work"
    demux_fastq = outdir / "demux" / "fastq" / "cellA_R1.fastq.gz"
    demux_fastq.parent.mkdir(parents=True, exist_ok=True)
    demux_fastq.write_text("demux", encoding="utf-8")
    raw_fastq = tmp_path / "raw_inputs" / "pooled_R1.fastq.gz"
    raw_fastq.parent.mkdir(parents=True, exist_ok=True)
    raw_fastq.write_text("raw", encoding="utf-8")

    plan = build_cleanup_plan(outdir=outdir, scope="all", workflow_succeeded=True)
    delete_paths = [item["path"] for item in plan["delete_candidates"]]

    assert str(demux_fastq.resolve()) in delete_paths
    assert str(raw_fastq.resolve()) not in delete_paths
    assert "original raw FASTQs" in plan["user_inputs_never_delete"]


def test_build_cleanup_plan_excludes_matrix_anndata_qc_and_logs(tmp_path: Path) -> None:
    outdir = tmp_path / "work"
    (outdir / "align").mkdir(parents=True, exist_ok=True)
    (outdir / "matrix").mkdir(parents=True, exist_ok=True)
    (outdir / "anndata").mkdir(parents=True, exist_ok=True)
    (outdir / "qc").mkdir(parents=True, exist_ok=True)
    (outdir / "logs").mkdir(parents=True, exist_ok=True)

    (outdir / "align" / "cache.sam").write_text("sam", encoding="utf-8")
    (outdir / "matrix" / "circ_counts.mtx").write_text("matrix", encoding="utf-8")
    (outdir / "anndata" / "circ_counts.h5ad").write_text("h5ad", encoding="utf-8")
    (outdir / "qc" / "cell_qc.tsv").write_text("qc", encoding="utf-8")
    (outdir / "logs" / "run.log").write_text("log", encoding="utf-8")

    plan = build_cleanup_plan(outdir=outdir, scope="all", workflow_succeeded=True)
    delete_paths = [item["path"] for item in plan["delete_candidates"]]

    assert any(path.endswith("cache.sam") for path in delete_paths)
    assert not any(path.endswith("circ_counts.mtx") for path in delete_paths)
    assert not any(path.endswith("circ_counts.h5ad") for path in delete_paths)
    assert not any(path.endswith("cell_qc.tsv") for path in delete_paths)
    assert not any(path.endswith("run.log") for path in delete_paths)


def test_build_cleanup_plan_for_failed_workflow_does_not_plan_deletion(tmp_path: Path) -> None:
    outdir = tmp_path / "work"
    align_sam = outdir / "align" / "failed.sam"
    align_sam.parent.mkdir(parents=True, exist_ok=True)
    align_sam.write_text("sam", encoding="utf-8")

    plan = build_cleanup_plan(outdir=outdir, scope="all", workflow_succeeded=False)

    assert plan["workflow_succeeded"] is False
    assert plan["candidate_count"] == 0
    assert plan["candidate_bytes"] == 0
    assert plan["delete_candidates"] == []


def test_execute_cleanup_plan_deletes_only_alignment_intermediates(tmp_path: Path) -> None:
    outdir = tmp_path / "work"
    align_sam = outdir / "align" / "cache.sam"
    matrix = outdir / "matrix" / "circ_counts.mtx"
    align_sam.parent.mkdir(parents=True, exist_ok=True)
    matrix.parent.mkdir(parents=True, exist_ok=True)
    align_sam.write_text("sam", encoding="utf-8")
    matrix.write_text("matrix", encoding="utf-8")

    result = execute_cleanup_plan(outdir=outdir, scope="alignments", workflow_succeeded=True)

    assert result["cleanup_performed"] is True
    assert not align_sam.exists()
    assert matrix.exists()
    assert result["cleanup_reclaimed_bytes"] >= 3
    assert any(path.endswith("cache.sam") for path in result["cleanup_deleted_paths"])


def test_execute_cleanup_plan_skips_failed_workflow(tmp_path: Path) -> None:
    outdir = tmp_path / "work"
    align_sam = outdir / "align" / "failed.sam"
    align_sam.parent.mkdir(parents=True, exist_ok=True)
    align_sam.write_text("sam", encoding="utf-8")

    result = execute_cleanup_plan(outdir=outdir, scope="alignments", workflow_succeeded=False)

    assert result["cleanup_performed"] is False
    assert align_sam.exists()
    assert result["cleanup_deleted_paths"] == []
