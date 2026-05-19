from __future__ import annotations

from pathlib import Path

import pytest

from circyto.pipeline.gene_expression_velocity import (
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
