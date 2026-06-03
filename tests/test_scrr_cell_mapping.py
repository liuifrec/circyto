from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from typer.testing import CliRunner

from circyto.cli.circyto import app
from circyto.pipeline.scrr_cell_mapping import HAS_ANNDATA, HAS_MUDATA


runner = CliRunner()


def _write_soft(path: Path) -> None:
    text = """^SAMPLE = GSM8558852
!Sample_geo_accession = GSM8558852
!Sample_title = RNA_IMR90_A_100
!Sample_source_name_ch1 = IMR90 cell
!Sample_characteristics_ch1 = treatment: none
^SAMPLE = GSM8558875
!Sample_geo_accession = GSM8558875
!Sample_title = DNA_IMR90_A_100
!Sample_source_name_ch1 = IMR90 cell
!Sample_characteristics_ch1 = treatment: none
^SAMPLE = GSM8558853
!Sample_geo_accession = GSM8558853
!Sample_title = RNA_IMR90_A_101
!Sample_source_name_ch1 = IMR90 cell
!Sample_characteristics_ch1 = treatment: aphidicoline
"""
    if path.name.endswith(".gz"):
        with gzip.open(path, "wt", encoding="utf-8") as handle:
            handle.write(text)
    else:
        path.write_text(text, encoding="utf-8")


def test_build_scrr_cell_map_parses_soft_gz_and_reports_unpaired(tmp_path: Path) -> None:
    soft = tmp_path / "GSE278958_family.soft.gz"
    out = tmp_path / "scrr_cell_map.tsv"
    _write_soft(soft)

    result = runner.invoke(app, ["build-scrr-cell-map", "--soft", str(soft), "--out", str(out)])

    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["n_rows"] == 3
    assert payload["n_paired_canonical_cells"] == 1
    assert payload["rna_only_canonical_cells"] == ["IMR90_A_101"]

    df = pd.read_csv(out, sep="\t", keep_default_na=False)
    assert list(df.columns) == [
        "gsm_id",
        "rna_cell_id",
        "dna_cell_id",
        "canonical_cell_id",
        "sample_title",
        "molecule",
        "treatment",
        "source_name",
    ]
    rna = df[df["gsm_id"] == "GSM8558852"].iloc[0].to_dict()
    assert rna["rna_cell_id"] == "RNA_IMR90_A_100"
    assert rna["dna_cell_id"] == "DNA_IMR90_A_100"
    assert rna["canonical_cell_id"] == "IMR90_A_100"
    assert rna["molecule"] == "RNA"


def test_build_scrr_cell_map_fails_on_duplicate_rna_canonical_id(tmp_path: Path) -> None:
    soft = tmp_path / "dup.soft"
    soft.write_text(
        "^SAMPLE = GSM1\n"
        "!Sample_geo_accession = GSM1\n"
        "!Sample_title = RNA_IMR90_A_100\n"
        "^SAMPLE = GSM2\n"
        "!Sample_geo_accession = GSM2\n"
        "!Sample_title = RNA_IMR90_A_100\n",
        encoding="utf-8",
    )

    result = runner.invoke(app, ["build-scrr-cell-map", "--soft", str(soft), "--out", str(tmp_path / "map.tsv")])

    assert result.exit_code != 0
    assert "Duplicate Sample_title" in result.output


@pytest.mark.skipif(not (HAS_ANNDATA and HAS_MUDATA), reason="anndata or mudata not installed")
def test_remap_scrr_mudata_obs_from_gsm_to_canonical(tmp_path: Path) -> None:
    import anndata as ad
    import mudata as mu

    cell_map = tmp_path / "scrr_cell_map.tsv"
    cell_map.write_text(
        "gsm_id\trna_cell_id\tdna_cell_id\tcanonical_cell_id\tsample_title\tmolecule\ttreatment\tsource_name\n"
        "GSM8558852\tRNA_IMR90_A_100\tDNA_IMR90_A_100\tIMR90_A_100\tRNA_IMR90_A_100\tRNA\tnone\tIMR90\n"
        "GSM8558853\tRNA_IMR90_A_101\tDNA_IMR90_A_101\tIMR90_A_101\tRNA_IMR90_A_101\tRNA\tnone\tIMR90\n",
        encoding="utf-8",
    )
    obs = pd.DataFrame(index=["GSM8558852", "GSM8558853"])
    rna = ad.AnnData(X=np.array([[1, 2], [3, 4]], dtype=np.int32), obs=obs.copy(), var=pd.DataFrame(index=["G1", "G2"]))
    circ = ad.AnnData(X=np.array([[0, 1], [1, 0]], dtype=np.int32), obs=obs.copy(), var=pd.DataFrame(index=["C1", "C2"]))
    input_h5mu = tmp_path / "rna_circ.h5mu"
    mu.MuData({"rna": rna, "circ": circ}).write_h5mu(input_h5mu)
    output_h5mu = tmp_path / "rna_circ_remapped.h5mu"

    result = runner.invoke(
        app,
        [
            "remap-scrr-mudata-obs",
            "--input",
            str(input_h5mu),
            "--cell-map",
            str(cell_map),
            "--output",
            str(output_h5mu),
        ],
    )

    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["target_id"] == "canonical_cell_id"

    remapped = mu.read_h5mu(output_h5mu)
    assert list(remapped.obs_names) == ["IMR90_A_100", "IMR90_A_101"]
    assert list(remapped.mod["rna"].obs_names) == ["IMR90_A_100", "IMR90_A_101"]
    assert remapped.mod["rna"].obs.loc["IMR90_A_100", "gsm_id"] == "GSM8558852"
    assert remapped.mod["circ"].obs.loc["IMR90_A_100", "original_circ_obs_id"] == "GSM8558852"


@pytest.mark.skipif(not (HAS_ANNDATA and HAS_MUDATA), reason="anndata or mudata not installed")
def test_merge_scrr_cnv_creates_trimodal_mudata(tmp_path: Path) -> None:
    import anndata as ad
    import mudata as mu

    obs = pd.DataFrame(
        {
            "gsm_id": ["GSM8558852", "GSM8558853"],
            "canonical_cell_id": ["IMR90_A_100", "IMR90_A_101"],
        },
        index=["IMR90_A_100", "IMR90_A_101"],
    )
    rna = ad.AnnData(X=np.array([[1], [2]], dtype=np.int32), obs=obs.copy(), var=pd.DataFrame(index=["G1"]))
    circ = ad.AnnData(X=np.array([[0], [3]], dtype=np.int32), obs=obs.copy(), var=pd.DataFrame(index=["C1"]))
    rna_circ_h5mu = tmp_path / "rna_circ_remapped.h5mu"
    mu.MuData({"rna": rna, "circ": circ}).write_h5mu(rna_circ_h5mu)

    cnv_obs = pd.DataFrame(
        {
            "dna_cell_id": ["DNA_IMR90_A_100", "DNA_IMR90_A_101"],
            "rna_cell_id": ["RNA_IMR90_A_100", "RNA_IMR90_A_101"],
            "canonical_cell_id": ["IMR90_A_100", "IMR90_A_101"],
        },
        index=["IMR90_A_100", "IMR90_A_101"],
    )
    cnv_var = pd.DataFrame(
        {"seqname": ["chr1", "chr1"], "start": [0, 50000], "end": [50000, 100000], "bin_size": [50000, 50000]},
        index=["chr1:0-50000", "chr1:50000-100000"],
    )
    cnv = ad.AnnData(X=np.array([[2, 3], [1, 2]], dtype=np.int16), obs=cnv_obs, var=cnv_var)
    cnv.layers["mappabilitynorm"] = np.array([[10, 12], [5, 8]], dtype=np.float32)
    cnv_h5ad = tmp_path / "cnv.h5ad"
    cnv.write_h5ad(cnv_h5ad)

    output = tmp_path / "trimodal.h5mu"
    result = runner.invoke(
        app,
        [
            "merge-scrr-cnv",
            "--input",
            str(rna_circ_h5mu),
            "--cnv",
            str(cnv_h5ad),
            "--output",
            str(output),
        ],
    )

    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["modalities"] == ["rna", "circ", "cnv"]
    assert payload["overlap_counts"]["n_trimodal_overlap"] == 2
    assert Path(payload["output_summary_json"]).exists()

    merged = mu.read_h5mu(output)
    assert list(merged.mod.keys()) == ["rna", "circ", "cnv"]
    assert list(merged.obs_names) == ["IMR90_A_100", "IMR90_A_101"]
    assert "mappabilitynorm" in merged.mod["cnv"].layers
    assert np.issubdtype(merged.mod["cnv"].X.dtype, np.integer)
