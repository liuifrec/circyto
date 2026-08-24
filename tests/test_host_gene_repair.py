from __future__ import annotations

from pathlib import Path
import json

import numpy as np
import pandas as pd
import pytest
from scipy import sparse as sp
from typer.testing import CliRunner

from circyto.cli.circyto import app
from circyto.multimodal.sync import mudata_from_modalities, read_h5mu, write_h5mu


runner = CliRunner()


HOST_COLUMNS = {
    "host_gene",
    "host_gene_source",
    "host_gene_from_gtf",
    "host_gene_from_circatlas",
    "host_gene_from_circatlas_id",
}


def _write_h5ad(path: Path, var: pd.DataFrame) -> None:
    ad = pytest.importorskip("anndata")
    obs = pd.DataFrame({"batch": ["b1", "b1"]}, index=["cellA", "cellB"])
    if var.shape[0] == 2:
        values = np.array([[1, 0], [0, 2]], dtype=np.int32)
    else:
        values = np.arange(1, (2 * max(1, var.shape[0])) + 1, dtype=np.int32).reshape(2, max(1, var.shape[0]))
    adata = ad.AnnData(X=sp.csr_matrix(values), obs=obs, var=var)
    adata.uns["keep"] = {"value": "metadata"}
    adata.write_h5ad(path)


def test_repair_h5ad_recovers_host_gene_from_circatlas_id(tmp_path: Path) -> None:
    ad = pytest.importorskip("anndata")
    input_path = tmp_path / "circ_counts.h5ad"
    output_path = tmp_path / "circ_counts.hostgene_fixed.h5ad"
    var = pd.DataFrame(
        {
            "chrom": ["chr1", "chr2"],
            "start": [100, 300],
            "end": [200, 400],
            "strand": ["+", "-"],
            "host_gene": ["", ""],
            "circatlas_v3_ids": ["hsa-UBAP2_0052", ""],
        },
        index=["circA", "circB"],
    )
    _write_h5ad(input_path, var)

    result = runner.invoke(
        app,
        ["repair-host-genes", "--input", str(input_path), "--output", str(output_path)],
    )

    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["n_host_gene_added"] == 1
    fixed = ad.read_h5ad(output_path)
    assert HOST_COLUMNS <= set(fixed.var.columns)
    assert fixed.var.loc["circA", "host_gene_from_circatlas_id"] == "UBAP2"
    assert fixed.var.loc["circA", "host_gene"] == "UBAP2"
    assert fixed.var.loc["circA", "host_gene_source"] == "circatlas_id"
    assert fixed.obs.equals(ad.read_h5ad(input_path).obs)
    assert fixed.X.toarray().tolist() == [[1, 0], [0, 2]]


def test_repair_h5ad_gtf_priority_over_circatlas_id(tmp_path: Path) -> None:
    ad = pytest.importorskip("anndata")
    input_path = tmp_path / "circ_counts.h5ad"
    output_path = tmp_path / "circ_counts.gtf_fixed.h5ad"
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(
        'chr1\ttest\tgene\t50\t250\t.\t+\t.\tgene_id "ENSG1"; gene_name "GTFGENE";\n',
        encoding="utf-8",
    )
    var = pd.DataFrame(
        {
            "chrom": ["chr1", "chr2"],
            "start": [100, 300],
            "end": [200, 400],
            "strand": ["+", "-"],
            "host_gene": ["", ""],
            "circatlas_v3_ids": ["hsa-UBAP2_0052", ""],
        },
        index=["circA", "circB"],
    )
    _write_h5ad(input_path, var)

    result = runner.invoke(
        app,
        [
            "repair-host-genes",
            "--input",
            str(input_path),
            "--output",
            str(output_path),
            "--gtf",
            str(gtf),
        ],
    )

    assert result.exit_code == 0, result.output
    fixed = ad.read_h5ad(output_path)
    assert fixed.var.loc["circA", "host_gene_from_gtf"] == "GTFGENE"
    assert fixed.var.loc["circA", "host_gene_from_circatlas_id"] == "UBAP2"
    assert fixed.var.loc["circA", "host_gene"] == "GTFGENE"
    assert fixed.var.loc["circA", "host_gene_source"] == "gtf"


def test_repair_h5ad_circatlas_table_priority_over_id_parsing(tmp_path: Path) -> None:
    ad = pytest.importorskip("anndata")
    input_path = tmp_path / "circ_counts.h5ad"
    output_path = tmp_path / "circ_counts.circatlas_fixed.h5ad"
    var = pd.DataFrame(
        {
            "chrom": ["chr1"],
            "start": [100],
            "end": [200],
            "strand": ["+"],
            "host_gene": [""],
            "circatlas_v3_ids": ["hsa-UBAP2_0052"],
            "circatlas_v3_host_genes": ["TABLEGENE"],
        },
        index=["circA"],
    )
    _write_h5ad(input_path, var)

    result = runner.invoke(
        app,
        ["repair-host-genes", "--input", str(input_path), "--output", str(output_path)],
    )

    assert result.exit_code == 0, result.output
    fixed = ad.read_h5ad(output_path)
    assert fixed.var.loc["circA", "host_gene_from_circatlas"] == "TABLEGENE"
    assert fixed.var.loc["circA", "host_gene_from_circatlas_id"] == "UBAP2"
    assert fixed.var.loc["circA", "host_gene"] == "TABLEGENE"
    assert fixed.var.loc["circA", "host_gene_source"] == "circatlas"


def test_repair_h5mu_patches_circ_modality_var(tmp_path: Path) -> None:
    ad = pytest.importorskip("anndata")
    mu = pytest.importorskip("mudata")
    input_path = tmp_path / "full_length.h5mu"
    output_path = tmp_path / "full_length.hostgene_fixed.h5mu"
    obs = pd.DataFrame(index=["cellA", "cellB"])
    rna = ad.AnnData(X=np.array([[1], [2]], dtype=np.int32), obs=obs.copy(), var=pd.DataFrame(index=["G1"]))
    circ_var = pd.DataFrame(
        {
            "chrom": ["chr1"],
            "start": [100],
            "end": [200],
            "strand": ["+"],
            "host_gene": [""],
            "circatlas_v3_ids": ["hsa-UBAP2_0052"],
        },
        index=["circA"],
    )
    circ = ad.AnnData(X=np.array([[0], [3]], dtype=np.int32), obs=obs.copy(), var=circ_var)
    write_h5mu(mudata_from_modalities({"rna": rna, "circ": circ}), input_path)

    result = runner.invoke(
        app,
        [
            "repair-host-genes",
            "--input",
            str(input_path),
            "--output",
            str(output_path),
            "--circ-mod",
            "circ",
        ],
    )

    assert result.exit_code == 0, result.output
    fixed = read_h5mu(output_path)
    circ_fixed = fixed.mod["circ"]
    assert HOST_COLUMNS <= set(circ_fixed.var.columns)
    assert circ_fixed.var.loc["circA", "host_gene_from_circatlas_id"] == "UBAP2"
    assert circ_fixed.var.loc["circA", "host_gene"] == "UBAP2"
    assert circ_fixed.var.loc["circA", "host_gene_source"] == "circatlas_id"
    assert fixed.var.loc["circA", "circ:host_gene"] == "UBAP2"
    circ_x = circ_fixed.X.toarray() if sp.issparse(circ_fixed.X) else np.asarray(circ_fixed.X)
    assert circ_x.tolist() == [[0], [3]]
