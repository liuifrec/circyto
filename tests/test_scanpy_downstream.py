from __future__ import annotations

import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

from circyto.cli.circyto import app
from circyto.pipeline.gene_expression_velocity import HAS_MUDATA
from circyto.pipeline.scanpy_downstream import HAS_SCANPY


runner = CliRunner()


def _write_mudata_fixture(tmp_path: Path) -> Path:
    if not HAS_MUDATA or not HAS_SCANPY:
        pytest.skip("mudata or scanpy not installed")
    import anndata as ad
    import mudata as mu
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp

    cells = [f"cell{i:02d}" for i in range(8)]
    obs = pd.DataFrame(
        {
            "membership": ["shared"] * 8,
            "total_rna_counts": [10, 12, 9, 15, 20, 8, 14, 11],
            "detected_genes": [3, 4, 3, 5, 5, 2, 4, 3],
            "mitochondrial_fraction": [0.1, 0.05, 0.2, 0.12, 0.03, 0.25, 0.08, 0.11],
            "ribosomal_fraction": [0.2, 0.15, 0.1, 0.05, 0.12, 0.08, 0.09, 0.07],
            "circRNA_count": [2, 1, 3, 2, 4, 1, 3, 2],
            "circRNA_total_support": [7, 3, 8, 5, 11, 2, 9, 6],
        },
        index=cells,
    )
    X_rna = sp.csr_matrix(
        np.array(
            [
                [5, 2, 1, 0, 2],
                [6, 3, 1, 1, 1],
                [3, 1, 4, 0, 1],
                [7, 4, 2, 1, 1],
                [8, 5, 3, 2, 2],
                [2, 1, 1, 0, 0],
                [6, 3, 2, 1, 2],
                [5, 2, 2, 1, 1],
            ],
            dtype=int,
        )
    )
    rna_var = pd.DataFrame(
        {
            "gene_id": [f"G{i}" for i in range(5)],
            "gene_name": ["MT-CO1", "GENE2", "GENE3", "GENE4", "RPLP0"],
            "gene_biotype": ["protein_coding"] * 5,
        },
        index=[f"G{i}" for i in range(5)],
    )
    X_circ = sp.csr_matrix(
        np.array(
            [
                [3, 1, 0],
                [1, 0, 0],
                [4, 2, 1],
                [2, 1, 0],
                [5, 3, 1],
                [1, 0, 0],
                [4, 2, 1],
                [3, 1, 1],
            ],
            dtype=int,
        )
    )
    circ_var = pd.DataFrame(
        {"circ_id": ["circ1", "circ2", "circ3"], "host_gene": ["GENE1", "GENE2", "GENE3"]},
        index=["circ1", "circ2", "circ3"],
    )
    rna = ad.AnnData(X=X_rna, obs=obs.copy(), var=rna_var)
    circ = ad.AnnData(X=X_circ, obs=obs.copy(), var=circ_var)
    mdata = mu.MuData({"rna": rna, "circ": circ})
    mdata.obs = obs.copy()
    mdata.uns["circyto"] = {
        "command_name": "circyto export-mudata",
        "circyto_version": "0.9.0",
        "source_workdir": "/tmp/workflow",
    }
    out_path = tmp_path / "full_length.h5mu"
    mdata.write_h5mu(out_path)
    return out_path


def test_scanpy_qc_report_outputs(tmp_path: Path) -> None:
    path = _write_mudata_fixture(tmp_path)
    outdir = tmp_path / "scanpy_qc"
    result = runner.invoke(app, ["scanpy-qc-report", "--input", str(path), "--output-dir", str(outdir)])
    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert (outdir / "qc_violin.png").exists()
    assert (outdir / "qc_scatter_counts_vs_genes.png").exists()
    assert (outdir / "qc_mt_vs_counts.png").exists()
    assert (outdir / "scanpy_qc_summary.json").exists()
    assert payload["exploratory_only"] is True


def test_scanpy_pca_outputs(tmp_path: Path) -> None:
    path = _write_mudata_fixture(tmp_path)
    outdir = tmp_path / "scanpy_pca"
    result = runner.invoke(app, ["scanpy-pca", "--input", str(path), "--output-dir", str(outdir)])
    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert (outdir / "rna_umap.png").exists()
    assert (outdir / "rna_leiden.tsv").exists()
    assert (outdir / "exploratory_rna_processed.h5ad").exists()
    assert (outdir / "scanpy_pca_summary.json").exists()
    assert payload["exploratory_only"] is True
