from __future__ import annotations

import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

from circyto.cli.circyto import app
from circyto.pipeline.scrr_cnv_import import HAS_ANNDATA, canonical_scrr_cell_id


runner = CliRunner()


def _write_cnv_tables(tmp_path: Path) -> tuple[Path, Path]:
    states = tmp_path / "summary_CNV_states_bin50kb.txt"
    states.write_text(
        "seqname\tstart\tend\tDNA_IMR90_A_100\tDNA_IMR90_A_101\n"
        "chr1\t0\t50000\t2-somy\t1-somy\n"
        "chr1\t50000\t100000\t3-somy\t2-somy\n",
        encoding="utf-8",
    )
    mappability = tmp_path / "summary_CNV_mappabilitynorm_bin50kb.txt"
    mappability.write_text(
        "seqname\tstart\tend\tDNA_IMR90_A_100\tDNA_IMR90_A_101\n"
        "chr1\t0\t50000\t10\t5\n"
        "chr1\t50000\t100000\t12\t8\n",
        encoding="utf-8",
    )
    return states, mappability


def test_canonical_scrr_cell_id_strips_one_modality_prefix() -> None:
    assert canonical_scrr_cell_id("DNA_IMR90_A_100") == "IMR90_A_100"
    assert canonical_scrr_cell_id("RNA_IMR90_A_100") == "IMR90_A_100"
    assert canonical_scrr_cell_id("IMR90_A_100") == "IMR90_A_100"


@pytest.mark.skipif(not HAS_ANNDATA, reason="anndata not installed")
def test_import_scrr_cnv_writes_h5ad_and_metadata(tmp_path: Path) -> None:
    import anndata as ad

    states, mappability = _write_cnv_tables(tmp_path)
    outdir = tmp_path / "cnv"

    result = runner.invoke(
        app,
        [
            "import-scrr-cnv",
            "--cnv-states",
            str(states),
            "--cnv-mappabilitynorm",
            str(mappability),
            "--outdir",
            str(outdir),
        ],
    )

    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["n_cells"] == 2
    assert payload["n_bins"] == 2
    assert payload["bin_sizes"] == [50000]
    assert payload["copy_number_state_counts"] == {"1": 1, "2": 2, "3": 1}
    assert payload["modality_name"] == "cnv"

    adata = ad.read_h5ad(outdir / "cnv.h5ad")
    assert list(adata.obs_names) == ["IMR90_A_100", "IMR90_A_101"]
    assert list(adata.var_names) == ["chr1:0-50000", "chr1:50000-100000"]
    assert adata.X.tolist() == [[2, 3], [1, 2]]
    assert "mappabilitynorm" in adata.layers
    assert adata.obs.loc["IMR90_A_100", "dna_cell_id"] == "DNA_IMR90_A_100"
    assert adata.obs.loc["IMR90_A_100", "rna_cell_id"] == "RNA_IMR90_A_100"
    assert adata.uns["circyto"]["X_semantics"].startswith("integer copy number")

    cells_text = (outdir / "cnv_cells.tsv").read_text(encoding="utf-8")
    bins_text = (outdir / "cnv_bins.tsv").read_text(encoding="utf-8")
    assert "canonical_cell_id" in cells_text
    assert "bin_size" in bins_text


def test_import_scrr_cnv_no_h5ad_allows_metadata_only(tmp_path: Path) -> None:
    states, _ = _write_cnv_tables(tmp_path)
    outdir = tmp_path / "cnv"

    result = runner.invoke(
        app,
        [
            "import-scrr-cnv",
            "--cnv-states",
            str(states),
            "--outdir",
            str(outdir),
            "--no-h5ad",
        ],
    )

    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["output_cnv_h5ad"] is None
    assert (outdir / "cnv_cells.tsv").exists()
    assert (outdir / "cnv_bins.tsv").exists()


def test_import_scrr_cnv_fails_on_mismatched_mappability_bins(tmp_path: Path) -> None:
    states, mappability = _write_cnv_tables(tmp_path)
    mappability.write_text(
        "seqname\tstart\tend\tDNA_IMR90_A_100\tDNA_IMR90_A_101\n"
        "chr1\t0\t50000\t10\t5\n"
        "chr1\t60000\t110000\t12\t8\n",
        encoding="utf-8",
    )

    result = runner.invoke(
        app,
        [
            "import-scrr-cnv",
            "--cnv-states",
            str(states),
            "--cnv-mappabilitynorm",
            str(mappability),
            "--outdir",
            str(tmp_path / "cnv"),
        ],
    )

    assert result.exit_code != 0
    assert "genomic bins do not match" in result.output
