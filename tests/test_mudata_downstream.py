from __future__ import annotations

import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

from circyto.cli.circyto import app
from circyto.pipeline.gene_expression_velocity import HAS_MUDATA


runner = CliRunner()


def _write_mudata_fixture(tmp_path: Path) -> Path:
    if not HAS_MUDATA:
        pytest.skip("mudata not installed")
    import anndata as ad
    import mudata as mu
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp

    obs = pd.DataFrame(
        {
            "membership": ["shared", "rna_only", "circ_only"],
            "total_rna_counts": [10, 5, 0],
            "detected_genes": [3, 2, 0],
            "mitochondrial_fraction": [0.1, 0.2, 0.0],
            "ribosomal_fraction": [0.3, 0.1, 0.0],
            "circRNA_count": [2, 0, 4],
            "circRNA_total_support": [7, 0, 9],
        },
        index=["cellA", "DIYHEK_192", "circ_only_1"],
    )
    rna = ad.AnnData(
        X=sp.csr_matrix(np.array([[5, 1], [2, 3], [0, 0]], dtype=int)),
        obs=obs.copy(),
        var=pd.DataFrame(
            {"gene_id": ["G1", "G2"], "gene_name": ["GENE1", "GENE2"], "gene_biotype": ["protein_coding", "lncRNA"]},
            index=["G1", "G2"],
        ),
    )
    circ = ad.AnnData(
        X=sp.csr_matrix(np.array([[4, 3], [0, 0], [7, 2]], dtype=int)),
        obs=obs.copy(),
        var=pd.DataFrame(
            {"circ_id": ["circ1", "circ2"], "chrom": ["chr1", "chr1"], "host_gene": ["GENE1", "GENE2"]},
            index=["circ1", "circ2"],
        ),
    )
    mdata = mu.MuData({"rna": rna, "circ": circ})
    mdata.obs = obs.copy()
    mdata.uns["circyto"] = {
        "command_name": "circyto export-mudata",
        "circyto_version": "0.9.0",
        "source_workdir": "/tmp/workflow",
        "workflow_summary_json": "{}",
    }
    out_path = tmp_path / "full_length.h5mu"
    mdata.write_h5mu(out_path)
    return out_path


def test_inspect_mudata_json(tmp_path: Path) -> None:
    path = _write_mudata_fixture(tmp_path)
    result = runner.invoke(app, ["inspect-mudata", "--input", str(path), "--json"])
    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["modalities"] == ["rna", "circ"]
    assert payload["n_obs"] == 3
    assert payload["rna_shape"] == [3, 2]
    assert payload["circ_shape"] == [3, 2]
    assert "membership" in payload["obs_columns"]
    assert "gene_biotype" in payload["rna_var_columns"]
    assert "host_gene" in payload["circ_var_columns"]
    assert "command_name" in payload["circyto_uns_keys"]
    assert payload["n_shared_cells"] == 1
    assert payload["n_rna_only_cells"] == 1
    assert payload["n_circ_only_cells"] == 1


def test_summarize_mudata_qc_json(tmp_path: Path) -> None:
    path = _write_mudata_fixture(tmp_path)
    result = runner.invoke(app, ["summarize-mudata-qc", "--input", str(path), "--json"])
    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["n_obs"] == 3
    assert payload["total_rna_counts"]["max"] == 10.0
    assert payload["detected_genes"]["median"] == 2.0
    assert payload["circRNA_count"]["max"] == 4.0
    assert payload["circRNA_total_support"]["max"] == 9.0
    assert payload["pearson_total_rna_vs_circRNA_count"] is not None
