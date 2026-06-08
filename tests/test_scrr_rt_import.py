from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from typer.testing import CliRunner

from circyto.cli.circyto import app
from circyto.pipeline.scrr_rt_import import HAS_ANNDATA, HAS_MUDATA, infer_paired_dna_cell_id


runner = CliRunner()


def _write_rt_table(tmp_path: Path) -> Path:
    table = tmp_path / "GSE278952_05_scRR-seq-DNA_HAP1_human_binarized_selectedsamples_all_sorted_hg38.txt"
    table.write_text(
        "seqname\tstart\tend\tDNA_HAP1_A_001\tDNA_HAP1_A_002\n"
        "chr1\t0\t50000\t0\t1\n"
        "chr1\t50000\t100000\t1\t1\n",
        encoding="utf-8",
    )
    return table


def test_infer_paired_dna_cell_id() -> None:
    assert infer_paired_dna_cell_id("RNA_HAP1_A_001") == "DNA_HAP1_A_001"
    assert infer_paired_dna_cell_id("DNA_HAP1_A_001") == "DNA_HAP1_A_001"
    assert infer_paired_dna_cell_id("HAP1_A_001") == "HAP1_A_001"


def test_import_scrr_rt_no_h5ad_writes_metadata_outputs(tmp_path: Path) -> None:
    table = _write_rt_table(tmp_path)
    outdir = tmp_path / "rt"

    result = runner.invoke(
        app,
        [
            "import-scrr-rt",
            "--rt-table",
            str(table),
            "--outdir",
            str(outdir),
            "--no-h5ad",
        ],
    )

    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["modality_name"] == "rt"
    assert payload["n_cells"] == 2
    assert payload["n_features"] == 2
    assert payload["value_semantics"] == "binary replication state"
    assert payload["output_rt_h5ad"] is None
    assert Path(payload["output_rt_cells"]).exists()
    assert Path(payload["output_rt_features"]).exists()

    cells = pd.read_csv(outdir / "rt_cells.tsv", sep="\t")
    assert cells["cell_id"].tolist() == ["HAP1_A_001", "HAP1_A_002"]
    assert cells["dna_cell_id"].tolist() == ["DNA_HAP1_A_001", "DNA_HAP1_A_002"]
    assert cells["rna_cell_id"].tolist() == ["RNA_HAP1_A_001", "RNA_HAP1_A_002"]

    features = pd.read_csv(outdir / "rt_features.tsv", sep="\t")
    assert features["feature_id"].tolist() == ["chr1:0-50000", "chr1:50000-100000"]
    assert features["feature_type"].tolist() == ["genomic_bin", "genomic_bin"]


def test_import_scrr_rt_detects_hap1_mids_sorted_bin_columns(tmp_path: Path) -> None:
    table = tmp_path / "hap1_sorted_bins.txt"
    table.write_text(
        "seqname\tstart\tend\tHAP1_scRR1_MidS_053\tHAP1_scRR1_MidS_036\n"
        "chr1\t0\t50000\t-1\t0\n"
        "chr1\t50000\t100000\t1\t-1\n",
        encoding="utf-8",
    )
    outdir = tmp_path / "rt"

    result = runner.invoke(
        app,
        [
            "import-scrr-rt",
            "--rt-table",
            str(table),
            "--outdir",
            str(outdir),
            "--no-h5ad",
        ],
    )

    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["n_cells"] == 2
    assert payload["n_features"] == 2
    assert payload["observed_values"] == [-1, 0, 1]

    cells = pd.read_csv(outdir / "rt_cells.tsv", sep="\t", keep_default_na=False)
    assert cells["cell_id"].tolist() == ["HAP1_scRR1_MidS_053", "HAP1_scRR1_MidS_036"]
    assert cells["canonical_cell_id"].tolist() == ["HAP1_scRR1_MidS_053", "HAP1_scRR1_MidS_036"]
    assert cells["dna_cell_id"].tolist() == ["HAP1_scRR1_MidS_053", "HAP1_scRR1_MidS_036"]
    assert cells["rna_cell_id"].tolist() == ["", ""]

    features = pd.read_csv(outdir / "rt_features.tsv", sep="\t")
    assert features["feature_id"].tolist() == ["chr1:0-50000", "chr1:50000-100000"]
    assert features["feature_type"].tolist() == ["genomic_bin", "genomic_bin"]


def test_import_scrr_rt_detects_hap1_mids_geneintersect_columns(tmp_path: Path) -> None:
    table = tmp_path / "hap1_geneintersect.txt"
    table.write_text(
        "gene_id\tchr\tstart\tend\tgene_name\tHAP1_scRR1_MidS_053\tHAP1_scRR1_MidS_074\n"
        "ENSG00000000419.14\tchr20\t50934867\t50958555\tDPM1\t-1\t1\n"
        "ENSG00000000457.14\tchr1\t169849631\t169894267\tSCYL3\t0\t-1\n",
        encoding="utf-8",
    )
    outdir = tmp_path / "rt"

    result = runner.invoke(
        app,
        [
            "import-scrr-rt",
            "--rt-table",
            str(table),
            "--outdir",
            str(outdir),
            "--no-h5ad",
        ],
    )

    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["n_cells"] == 2
    assert payload["n_features"] == 2

    cells = pd.read_csv(outdir / "rt_cells.tsv", sep="\t", keep_default_na=False)
    assert cells["cell_id"].tolist() == ["HAP1_scRR1_MidS_053", "HAP1_scRR1_MidS_074"]
    assert cells["rna_cell_id"].tolist() == ["", ""]

    features = pd.read_csv(outdir / "rt_features.tsv", sep="\t")
    assert features["feature_id"].tolist() == ["ENSG00000000419.14", "ENSG00000000457.14"]
    assert features["feature_type"].tolist() == [
        "gene_intersect_genomic_feature",
        "gene_intersect_genomic_feature",
    ]
    assert features["gene_id"].tolist() == ["ENSG00000000419.14", "ENSG00000000457.14"]
    assert features["gene_name"].tolist() == ["DPM1", "SCYL3"]
    assert features["seqname"].tolist() == ["chr20", "chr1"]


def test_import_scrr_rt_no_sample_columns_error_includes_header_context(tmp_path: Path) -> None:
    table = tmp_path / "not_rt_samples.txt"
    table.write_text(
        "seqname\tstart\tend\tunexpected_sample\n"
        "chr1\t0\t50000\t0\n",
        encoding="utf-8",
    )
    outdir = tmp_path / "rt"

    result = runner.invoke(
        app,
        [
            "import-scrr-rt",
            "--rt-table",
            str(table),
            "--outdir",
            str(outdir),
            "--no-h5ad",
        ],
    )

    assert result.exit_code != 0
    assert "Detected columns count: 4" in result.output
    assert "First 15 columns:" in result.output
    assert "unexpected_sample" in result.output
    assert "Accepted sample prefixes/patterns:" in result.output
    assert "HAP1_scRR1_MidS_*" in result.output


@pytest.mark.skipif(not HAS_ANNDATA, reason="anndata not installed")
def test_import_scrr_rt_writes_h5ad_and_avg_rt_var(tmp_path: Path) -> None:
    import anndata as ad

    table = _write_rt_table(tmp_path)
    bedgraph = tmp_path / "GSE278952_05_scRR-seq-DNA_HAP1_human_midS_Avg_RT_hg38.bedGraph"
    bedgraph.write_text(
        "chr1\t0\t50000\t0.25\n"
        "chr1\t50000\t100000\t0.75\n",
        encoding="utf-8",
    )
    outdir = tmp_path / "rt"

    result = runner.invoke(
        app,
        [
            "import-scrr-rt",
            "--rt-table",
            str(table),
            "--avg-rt-bedgraph",
            str(bedgraph),
            "--outdir",
            str(outdir),
        ],
    )

    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["avg_rt_bedgraph"]["avg_rt_storage"] == "var['avg_rt']"
    assert payload["output_rt_h5ad"].endswith("rt.h5ad")

    adata = ad.read_h5ad(outdir / "rt.h5ad")
    assert list(adata.obs_names) == ["HAP1_A_001", "HAP1_A_002"]
    assert list(adata.var_names) == ["chr1:0-50000", "chr1:50000-100000"]
    assert adata.X.tolist() == [[0, 1], [1, 1]]
    assert adata.var["avg_rt"].tolist() == [0.25, 0.75]
    assert adata.uns["circyto"]["modality_name"] == "rt"
    assert adata.uns["circyto"]["X_semantics"] == "binary replication state"


def test_import_scrr_rt_uses_cell_mapping_and_rna_obs_strategy(tmp_path: Path) -> None:
    table = _write_rt_table(tmp_path)
    mapping = tmp_path / "scrr_cell_map.tsv"
    mapping.write_text(
        "gsm_id\trna_cell_id\tdna_cell_id\tcanonical_cell_id\tsample_title\tmolecule\ttreatment\tsource_name\n"
        "GSM1\tRNA_HAP1_A_001\tDNA_HAP1_A_001\tHAP1_A_001\tDNA_HAP1_A_001\tDNA\tnone\tHAP1\n"
        "GSM2\tRNA_HAP1_A_002\tDNA_HAP1_A_002\tHAP1_A_002\tDNA_HAP1_A_002\tDNA\tnone\tHAP1\n",
        encoding="utf-8",
    )
    outdir = tmp_path / "rt"

    result = runner.invoke(
        app,
        [
            "import-scrr-rt",
            "--rt-table",
            str(table),
            "--cell-mapping",
            str(mapping),
            "--obs-id-strategy",
            "rna",
            "--outdir",
            str(outdir),
            "--no-h5ad",
        ],
    )

    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["n_cells_with_mapping"] == 2
    cells = pd.read_csv(outdir / "rt_cells.tsv", sep="\t")
    assert cells["cell_id"].tolist() == ["RNA_HAP1_A_001", "RNA_HAP1_A_002"]
    assert cells["canonical_cell_id"].tolist() == ["HAP1_A_001", "HAP1_A_002"]
    assert cells["pairing_strategy"].tolist() == ["cell_mapping", "cell_mapping"]


def test_import_scrr_rt_categorical_values_are_encoded(tmp_path: Path) -> None:
    table = tmp_path / "rt_states.txt"
    table.write_text(
        "gene_id\tDNA_HAP1_A_001\tDNA_HAP1_A_002\n"
        "GENE1\tlate\tearly\n"
        "GENE2\tearly\tearly\n",
        encoding="utf-8",
    )
    outdir = tmp_path / "rt"

    result = runner.invoke(
        app,
        [
            "import-scrr-rt",
            "--rt-table",
            str(table),
            "--outdir",
            str(outdir),
            "--no-h5ad",
        ],
    )

    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["value_semantics"] == "categorical replication state encoded as integer labels"
    assert payload["category_encoding"] == {"early": 0, "late": 1}
    features = pd.read_csv(outdir / "rt_features.tsv", sep="\t")
    assert features["feature_id"].tolist() == ["GENE1", "GENE2"]
    assert features["feature_type"].tolist() == ["gene_intersect_feature", "gene_intersect_feature"]


@pytest.mark.skipif(not (HAS_ANNDATA and HAS_MUDATA), reason="anndata or mudata not installed")
def test_merge_scrr_rt_creates_rna_circ_rt_mudata(tmp_path: Path) -> None:
    import anndata as ad
    import mudata as mu

    obs = pd.DataFrame(
        {
            "gsm_id": ["GSM1", "GSM2"],
            "canonical_cell_id": ["HAP1_A_001", "HAP1_A_002"],
        },
        index=["HAP1_A_001", "HAP1_A_002"],
    )
    rna = ad.AnnData(X=np.array([[1], [2]], dtype=np.int32), obs=obs.copy(), var=pd.DataFrame(index=["G1"]))
    circ = ad.AnnData(X=np.array([[0], [3]], dtype=np.int32), obs=obs.copy(), var=pd.DataFrame(index=["C1"]))
    input_h5mu = tmp_path / "rna_circ_remapped.h5mu"
    mu.MuData({"rna": rna, "circ": circ}).write_h5mu(input_h5mu)

    rt_obs = pd.DataFrame(
        {
            "dna_cell_id": ["DNA_HAP1_A_001", "DNA_HAP1_A_002"],
            "rna_cell_id": ["RNA_HAP1_A_001", "RNA_HAP1_A_002"],
            "canonical_cell_id": ["HAP1_A_001", "HAP1_A_002"],
        },
        index=["HAP1_A_001", "HAP1_A_002"],
    )
    rt_var = pd.DataFrame(
        {"seqname": ["chr1", "chr1"], "start": [0, 50000], "end": [50000, 100000]},
        index=["chr1:0-50000", "chr1:50000-100000"],
    )
    rt = ad.AnnData(X=np.array([[0, 1], [1, 1]], dtype=np.int8), obs=rt_obs, var=rt_var)
    rt_h5ad = tmp_path / "rt.h5ad"
    rt.write_h5ad(rt_h5ad)

    output = tmp_path / "rna_circ_rt.h5mu"
    result = runner.invoke(
        app,
        [
            "merge-scrr-rt",
            "--input",
            str(input_h5mu),
            "--rt",
            str(rt_h5ad),
            "--output",
            str(output),
        ],
    )

    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["modalities"] == ["rna", "circ", "rt"]
    assert payload["overlap_counts"]["n_trimodal_overlap"] == 2
    assert Path(payload["output_summary_json"]).exists()

    merged = mu.read_h5mu(output)
    assert list(merged.mod.keys()) == ["rna", "circ", "rt"]
    assert list(merged.obs_names) == ["HAP1_A_001", "HAP1_A_002"]
    assert np.issubdtype(merged.mod["rt"].X.dtype, np.integer)
