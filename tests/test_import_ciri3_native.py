"""Native acceptance fixtures are unchanged official output, not H5MU extractions."""

import csv
import hashlib
import json
import shutil
from pathlib import Path

import anndata as ad
import pandas as pd
import pytest
from scipy import sparse
from typer.testing import CliRunner

from circyto.cli.circyto import app
from circyto.pipeline.import_circ import import_circ
from circyto.pipeline.import_ciri3 import NATIVE_HEADER

FIXTURE = Path(__file__).resolve().parents[1] / "examples/native_ciri3"


@pytest.fixture
def native(tmp_path):
    for name in ("result.txt", "manifest.tsv", "cells.tsv"):
        shutil.copyfile(FIXTURE / name, tmp_path / name)
    return {
        "ciri3_manifest": tmp_path / "manifest.tsv",
        "cells": tmp_path / "cells.tsv",
        "outdir": tmp_path / "output",
        "genome_build": "upstream_demo_chr1_unverified",
        "coordinate_convention": "1-based-closed",
        "support_definition": "CIRI3_junction_reads",
        "detector": "CIRI3",
    }


def test_official_native_values_independently(native):
    manifest = json.loads((FIXTURE / "provenance.json").read_text())
    for name, entry in manifest["files"].items():
        assert (
            hashlib.sha256((FIXTURE / name).read_bytes()).hexdigest() == entry["sha256"]
        )
    metrics = import_circ(**native)
    result = ad.read_h5ad(native["outdir"] / "anndata/circ_counts.h5ad")
    with (FIXTURE / "result.txt").open() as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    assert result.var_names.tolist() == [r["circRNA_ID"] for r in rows]
    assert (
        result.X.data.tolist()
        == [int(r["#junction_reads"]) for r in rows]
        == [9, 9, 5, 6, 14, 7]
    )
    for output, source in [
        ("chrom", "chr"),
        ("strand", "strand"),
        ("ciri3_gene_id", "gene_id"),
        ("ciri3_circRNA_type", "circRNA_type"),
    ]:
        assert result.var[output].tolist() == [r[source] for r in rows]
    assert result.var.start.tolist() == [int(r["circRNA_start"]) for r in rows]
    assert result.var.end.tolist() == [int(r["circRNA_end"]) for r in rows]
    assert result.obs_names.tolist() == ["upstream_single"]
    assert result.obs.circ_observation_status.tolist() == ["unknown"]
    assert metrics["total_support"] == 50
    assert metrics["detection_fraction_completed"]["value"] is None
    assert result.var.host_gene.eq("").all() and result.var.host_gene_id.eq("").all()
    assert metrics["native_ciri3"]["samples"][0]["detector_version"] == "unknown"
    evidence = native["outdir"] / "evidence/ciri3_native_000001.tsv"
    assert evidence.read_bytes() == (FIXTURE / "result.txt").read_bytes()
    assert (
        metrics["source_sha256"][evidence.name]
        == manifest["files"]["result.txt"]["sha256"]
    )
    assert (
        json.loads(result.uns["circyto"]["native_ciri3_json"])
        == metrics["native_ciri3"]
    )
    report = (native["outdir"] / "qc/report.html").read_text()
    assert "unavailable" in report and "unknown" in report


def test_official_simulated_annotation_contract(native):
    shutil.copyfile(
        FIXTURE / "annotated_simulation.txt", native["cells"].parent / "result.txt"
    )
    import_circ(**native)
    a = ad.read_h5ad(native["outdir"] / "anndata/circ_counts.h5ad")
    source = pd.read_csv(
        FIXTURE / "annotated_simulation.txt", sep="\t", dtype=str, keep_default_na=False
    )
    assert a.var.ciri3_gene_id.tolist() == source.gene_id.tolist()
    assert a.var.ciri3_circRNA_type.tolist() == source.circRNA_type.tolist()
    assert a.var.host_gene.eq("").all()  # native gene IDs are not gene symbols
    assert (
        a.var.loc[a.var.host_gene_id != "", "host_gene_id_source"]
        .eq("CIRI3_gene_id")
        .all()
    )
    assert a.X.data.tolist() == source["#junction_reads"].astype(int).tolist()


@pytest.mark.parametrize(
    "mode",
    [
        "normalized",
        "missing_score",
        "extra_column",
        "reordered",
        "duplicate",
        "bad_id",
        "strand",
        "fraction",
        "negative",
        "row_width",
    ],
)
def test_reject_unsupported_native(native, mode):
    p = native["cells"].parent / "result.txt"
    lines = p.read_text().splitlines()
    if mode == "normalized":
        lines = ["circ_id\tchr\tstart\tend\tstrand\tsupport", "x\tchr1\t1\t2\t+\t2"]
    elif mode == "missing_score":
        lines = ["\t".join(x.split("\t")[:-1]) for x in lines]
    elif mode == "extra_column":
        lines = [x + "\textra" for x in lines]
    elif mode == "reordered":
        lines = [
            "\t".join([x.split("\t")[1], x.split("\t")[0], *x.split("\t")[2:]])
            for x in lines
        ]
    elif mode == "duplicate":
        lines.append(lines[1])
    elif mode == "bad_id":
        lines[1] = lines[1].replace("chr1:1379083|1414681", "chr1:1379082|1414681")
    elif mode == "strand":
        lines[1] = lines[1].replace("\t-\t", "\t.\t")
    elif mode in ["fraction", "negative"]:
        fields = lines[1].split("\t")
        fields[4] = "1.5" if mode == "fraction" else "-1"
        lines[1] = "\t".join(fields)
    elif mode == "row_width":
        lines[1] += "\textra"
    p.write_text("\n".join(lines) + "\n")
    with pytest.raises(ValueError):
        import_circ(**native)
    assert not native["outdir"].exists()


@pytest.mark.parametrize(
    "key,value",
    [
        ("coordinate_convention", "0-based-half-open"),
        ("support_definition", "Score"),
        ("detector", "CIRI2"),
    ],
)
def test_native_requires_explicit_correct_declarations(native, key, value):
    with pytest.raises(ValueError):
        import_circ(**{**native, key: value})


@pytest.mark.parametrize(
    "mode",
    [
        "cell_duplicate",
        "sample_duplicate",
        "file_duplicate",
        "version",
        "blank_version",
        "pooled",
        "unknown_cell",
        "extra_column",
    ],
)
def test_invalid_mapping(native, mode):
    p = native["ciri3_manifest"]
    lines = p.read_text().splitlines()
    if mode == "cell_duplicate":
        lines.append(lines[1])
    elif mode == "sample_duplicate":
        lines.append(lines[1].replace("upstream_single\t", "other\t", 1))
    elif mode == "file_duplicate":
        lines.append("other\tother_sample\tresult.txt\tunknown\tsingle-sample")
    elif mode == "version":
        lines[1] = lines[1].replace("\tunknown\t", "\t4.0\t")
    elif mode == "blank_version":
        lines[1] = lines[1].replace("\tunknown\t", "\t\t")
    elif mode == "pooled":
        lines[1] = lines[1].replace("single-sample", "pooled")
    elif mode == "unknown_cell":
        lines[1] = lines[1].replace("upstream_single\t", "other\t", 1)
    elif mode == "extra_column":
        lines = [x + "\tother" for x in lines]
    p.write_text("\n".join(lines) + "\n")
    with pytest.raises(ValueError):
        import_circ(**native)


def test_header_only_zero_requires_explicit_status(native):
    (native["cells"].parent / "result.txt").write_text("\t".join(NATIVE_HEADER) + "\n")
    native["cells"].write_text(
        "cell_id\tdetector_status\nupstream_single\tcompleted\nfailed\tfailed\nskipped\tskipped\nunprocessed\tunprocessed\nmissing\tunknown\n"
    )
    m = import_circ(**native)
    a = ad.read_h5ad(native["outdir"] / "anndata/circ_counts.h5ad")
    assert a.shape == (5, 0)
    assert a.obs.circ_observation_status.tolist() == [
        "completed_zero",
        "failed",
        "skipped",
        "unprocessed",
        "unknown",
    ]
    assert m["zero_fraction_completed"]["value"] == 1


def test_missing_completed_file_is_not_zero(native):
    native["cells"].write_text(
        "cell_id\tdetector_status\nupstream_single\tunknown\nmissing\tcompleted\n"
    )
    with pytest.raises(ValueError, match="Completed cells"):
        import_circ(**native)


def test_no_dense_candidate_matrix_and_cell_order(native, monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("Dense cell-by-candidate conversion")

    for cls in (sparse.csr_matrix, sparse.csc_matrix, sparse.coo_matrix):
        monkeypatch.setattr(cls, "toarray", forbidden)
        monkeypatch.setattr(cls, "todense", forbidden)
    native["cells"].write_text("cell_id\nmissing\nupstream_single\n")
    import_circ(**native)
    a = ad.read_h5ad(native["outdir"] / "anndata/circ_counts.h5ad")
    assert a.obs_names.tolist() == ["missing", "upstream_single"]
    assert a.X[0].nnz == 0 and a.X[1].sum() == 50


def test_cli_and_existing_rna_technical_integration(native):
    mu = pytest.importorskip("mudata")
    rna = ad.AnnData(
        sparse.csr_matrix([[2], [3]]),
        obs=pd.DataFrame(index=["rna_only", "upstream_single"]),
        var=pd.DataFrame(index=["test_gene"]),
    )
    rna.layers["original"] = rna.X.copy()
    path = native["cells"].parent / "rna.h5ad"
    rna.write_h5ad(path)
    args = ["import-circ"]
    for key, value in {**native, "rna": path}.items():
        args.extend(["--" + key.replace("_", "-"), str(value)])
    run = CliRunner().invoke(app, args)
    assert run.exit_code == 0, run.output
    m = mu.read_h5mu(native["outdir"] / "mudata/circyto_multimodal.h5mu")
    assert m["rna"].obs_names.tolist() == rna.obs_names.tolist()
    assert (m["rna"].X != rna.X).nnz == 0 and (
        m["rna"].layers["original"] != rna.X
    ).nnz == 0
    assert m.obs.loc["rna_only", "circ_observation_status"] == "unknown"
    assert m["circ"].X.sum() == 50


def test_input_modes_exclusive(native):
    with pytest.raises(ValueError, match="exactly one"):
        import_circ(**native, calls=FIXTURE / "result.txt")


def test_multiple_mapped_files_preserve_cells_not_aggregate(native):
    root = native["cells"].parent
    (root / "second.txt").write_text(
        (root / "result.txt").read_text().replace("\t9\t1_7_8", "\t11\t1_7_8")
    )
    with native["ciri3_manifest"].open("a") as f:
        f.write("second\tsecond_sample\tsecond.txt\tunknown\tsingle-sample\n")
    native["cells"].write_text(
        "cell_id\tdetector_status\nsecond\tfailed\nupstream_single\tcompleted\n"
    )
    m = import_circ(**native)
    a = ad.read_h5ad(native["outdir"] / "anndata/circ_counts.h5ad")
    assert a.obs_names.tolist() == ["second", "upstream_single"]
    assert a.X[0, 0] == 11 and a.X[1, 0] == 9
    assert a.obs.circ_observation_status.tolist() == ["failed", "completed_detected"]
    assert m["total_support"] == 102
    assert m["detection_fraction_completed"] == {
        "numerator": 1,
        "denominator": 1,
        "value": 1.0,
    }
    assert set(m["source_sha256"]) == {
        "cells.tsv",
        "ciri3_manifest.tsv",
        "ciri3_native_000001.tsv",
        "ciri3_native_000002.tsv",
    }


@pytest.mark.parametrize("conflict", ["strand", "gene_id", "version"])
def test_mixed_native_evidence_rejected(native, conflict):
    root = native["cells"].parent
    lines = (root / "result.txt").read_text().splitlines()
    fields = lines[1].split("\t")
    if conflict == "strand":
        fields[10] = "+"
    if conflict == "gene_id":
        fields[9] = "GENE_CONFLICT"
    lines[1] = "\t".join(fields)
    (root / "second.txt").write_text("\n".join(lines) + "\n")
    version = "3.0.1" if conflict == "version" else "unknown"
    with native["ciri3_manifest"].open("a") as f:
        f.write(f"second\tsecond_sample\tsecond.txt\t{version}\tsingle-sample\n")
    native["cells"].write_text("cell_id\nupstream_single\nsecond\n")
    with pytest.raises(ValueError):
        import_circ(**native)
    assert not native["outdir"].exists()


def test_native_chromosome_not_silently_stripped(native):
    path = native["cells"].parent / "result.txt"
    path.write_text(
        path.read_text().replace("chr1:", "chr1 :").replace("\tchr1\t", "\tchr1 \t")
    )
    with pytest.raises(ValueError, match="native CIRI3 chromosome"):
        import_circ(**native)
