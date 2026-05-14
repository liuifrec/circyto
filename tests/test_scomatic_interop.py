from __future__ import annotations

from pathlib import Path

import pandas as pd
from typer.testing import CliRunner

from circyto.cli.circyto import app


runner = CliRunner()


def test_export_scomatic_inputs_cli_writes_expected_files(tmp_path: Path) -> None:
    ref_fa = tmp_path / "reference.fa"
    ref_fa.write_text(">chr1\nACGT\n", encoding="utf-8")
    outdir = tmp_path / "scomatic_export"

    result = runner.invoke(
        app,
        [
            "export-scomatic-inputs",
            "--bam-manifest",
            "testdata/scomatic_interop/bam_manifest.tsv",
            "--cell-metadata",
            "testdata/scomatic_interop/cell_metadata.tsv",
            "--outdir",
            str(outdir),
            "--reference-fasta",
            str(ref_fa),
            "--protocol",
            "ramda",
        ],
    )

    assert result.exit_code == 0, result.stdout
    bam_list = outdir / "scomatic_bam_list.tsv"
    celltypes = outdir / "scomatic_celltypes.tsv"
    readme = outdir / "README_scomatic_next_steps.md"
    assert bam_list.exists()
    assert celltypes.exists()
    assert readme.exists()

    bam_df = pd.read_csv(bam_list, sep="\t")
    cell_df = pd.read_csv(celltypes, sep="\t")
    assert {"cell_id", "bam", "protocol", "reference_fasta"}.issubset(bam_df.columns)
    assert {"cell_id", "cell_type"}.issubset(cell_df.columns)
    assert bam_df["protocol"].tolist() == ["ramda", "ramda", "ramda"]

    readme_text = readme.read_text(encoding="utf-8")
    assert "does not install, invoke, or validate SComatic" not in readme_text
    assert "exploratory candidate SNVs" in readme_text


def test_join_circ_snv_summary_cli_writes_expected_files_and_columns(tmp_path: Path) -> None:
    outdir = tmp_path / "joined"

    result = runner.invoke(
        app,
        [
            "join-circ-snv-summary",
            "--circ-matrix",
            "testdata/scomatic_interop/circ_matrix.tsv",
            "--circ-feature-table",
            "testdata/scomatic_interop/circ_feature_table.tsv",
            "--scomatic-candidates",
            "testdata/scomatic_interop/scomatic_candidates.tsv",
            "--cell-metadata",
            "testdata/scomatic_interop/cell_metadata.tsv",
            "--outdir",
            str(outdir),
        ],
    )

    assert result.exit_code == 0, result.stdout
    cell_summary_path = outdir / "circ_snv_cell_summary.tsv"
    gene_summary_path = outdir / "circ_snv_host_gene_summary.tsv"
    assert cell_summary_path.exists()
    assert gene_summary_path.exists()

    cell_df = pd.read_csv(cell_summary_path, sep="\t")
    gene_df = pd.read_csv(gene_summary_path, sep="\t")

    assert {"cell_id", "circRNA_count", "circRNA_read_support_sum", "candidate_snv_count"}.issubset(cell_df.columns)
    assert {"gene", "circRNA_feature_count", "circRNA_cell_count", "circRNA_read_support_sum", "candidate_snv_count"}.issubset(gene_df.columns)

    cell_df = cell_df.set_index("cell_id")
    assert int(cell_df.loc["cellA", "circRNA_count"]) == 1
    assert float(cell_df.loc["cellA", "circRNA_read_support_sum"]) == 3.0
    assert int(cell_df.loc["cellA", "candidate_snv_count"]) == 2
    assert int(cell_df.loc["cellC", "candidate_snv_count"]) == 0

    gene_df = gene_df.set_index("gene")
    assert int(gene_df.loc["GENE1", "circRNA_feature_count"]) == 2
    assert int(gene_df.loc["GENE1", "candidate_snv_count"]) == 1
    assert int(gene_df.loc["GENE2", "candidate_snv_count"]) == 3
