import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from scipy.io import mmwrite
from typer.testing import CliRunner

from circyto.cli.circyto import app
from circyto.pipeline.import_circ import import_circ


@pytest.fixture
def inputs(tmp_path):
    calls = tmp_path / "calls.tsv"
    calls.write_text(
        "cell_id\tcirc_id\tchrom\tstart\tend\tstrand\tsupport\thost_gene\tevidence\n"
        "001\tNA\tchr1\t0\t20\t+\t7\tGENE1\t<script>alert(1)</script>\n"
        "failed\tNA\tchr1\t0\t20\t+\t2\tGENE1\tpartial\n"
    )
    cells = tmp_path / "cells.tsv"
    cells.write_text(
        "cell_id\tdetector_status\tbatch\nzero\tcompleted\t001\n001\tcompleted\t001\n"
        "failed\tfailed\t002\nskip\tskipped\t002\nnew\tunprocessed\t003\nNA\tunknown\t003\n"
    )
    return {
        "calls": calls,
        "cells": cells,
        "outdir": tmp_path / "result",
        "genome_build": "GRCh38",
        "coordinate_convention": "0-based-half-open",
        "support_definition": "BSJ_reads",
        "detector": "example",
    }


def test_counts_states_evidence_and_reproducibility(inputs):
    metrics = import_circ(**inputs)
    root = inputs["outdir"]
    a = ad.read_h5ad(root / "anndata/circ_counts.h5ad")
    assert a.obs_names.tolist() == ["zero", "001", "failed", "skip", "new", "NA"]
    assert a.var_names.tolist() == ["NA"]
    assert a.X.dtype == np.int64 and sparse.isspmatrix_csr(a.X)
    assert a.X.toarray().ravel().tolist() == [0, 7, 2, 0, 0, 0]
    assert a.obs.batch.tolist() == ["001", "001", "002", "002", "003", "003"]
    assert a.obs.circ_observation_status.tolist() == [
        "completed_zero",
        "completed_detected",
        "failed",
        "skipped",
        "unprocessed",
        "unknown",
    ]
    assert a.var.loc["NA", "host_gene"] == "GENE1"
    assert a.var.loc["NA", "host_gene_source"] == ""  # no fabricated GTF evidence
    assert a.var.loc["NA", "start"] == 0
    assert a.var.loc["NA", "n_cells_detected"] == 2
    assert a.var.loc["NA", "n_completed_cells_detected"] == 1
    assert a.var.loc["NA", "prevalence_completed"] == 0.5
    assert metrics["detection_fraction_completed"] == {
        "numerator": 1,
        "denominator": 2,
        "value": 0.5,
    }
    assert metrics["zero_fraction_completed"] == {
        "numerator": 1,
        "denominator": 2,
        "value": 0.5,
    }
    assert (root / "evidence/calls.tsv").read_bytes() == inputs["calls"].read_bytes()
    report = (root / "qc/report.html").read_text()
    assert "<script>" not in report and "biological absence" in report
    assert "1 / 2" in report
    second = root.with_name("second")
    assert import_circ(**{**inputs, "outdir": second}) == metrics
    for path in ["qc/metrics.json", "qc/report.html", "matrix/circ_counts.mtx"]:
        assert (root / path).read_bytes() == (second / path).read_bytes()
    with pytest.raises(ValueError, match="already exists"):
        import_circ(**inputs)


def test_matrix_roundtrip_and_missing_completed(inputs):
    import_circ(**inputs)
    kwargs = {
        **inputs,
        "calls": None,
        "matrix_dir": inputs["outdir"] / "matrix",
        "outdir": inputs["outdir"].with_name("second"),
    }
    import_circ(**kwargs)
    a = ad.read_h5ad(inputs["outdir"] / "anndata/circ_counts.h5ad")
    b = ad.read_h5ad(kwargs["outdir"] / "anndata/circ_counts.h5ad")
    assert (a.X != b.X).nnz == 0
    pd.testing.assert_frame_equal(a.obs, b.obs)
    pd.testing.assert_frame_equal(a.var, b.var)
    inputs["cells"].write_text(inputs["cells"].read_text() + "absent\tcompleted\t004\n")
    with pytest.raises(ValueError, match="missing from matrix"):
        import_circ(**{**kwargs, "outdir": kwargs["outdir"].with_name("third")})


def test_rna_native_axes_metadata_layers_and_overlap(inputs):
    pytest.importorskip("mudata")
    from circyto.multimodal.sync import read_h5mu

    rna = ad.AnnData(
        sparse.csr_matrix([[1.5], [3.0], [9.0]]),
        obs=pd.DataFrame({"sample": ["a", "b", "c"]}, index=["001", "zero", "RNAonly"]),
        var=pd.DataFrame(index=["g"]),
    )
    rna.layers["original"] = rna.X.copy()
    rna.obsm["embedding"] = np.array([[1, 2], [3, 4], [5, 6]])
    rna.uns["genome_build"] = "GRCh38"
    path = inputs["cells"].with_name("rna.h5ad")
    rna.write_h5ad(path)
    metrics = import_circ(**inputs, rna=path)
    m = read_h5mu(inputs["outdir"] / "mudata/circyto_multimodal.h5mu")
    assert metrics["rna_overlap"] == {
        "n_shared": 2,
        "n_circ_only": 4,
        "n_rna_only": 1,
        "rna_genome_build_verified": True,
    }
    assert m["rna"].obs_names.tolist() == rna.obs_names.tolist()
    assert m["circ"].obs_names.tolist() == [
        "zero",
        "001",
        "failed",
        "skip",
        "new",
        "NA",
    ]
    assert (m["rna"].X != rna.X).nnz == 0
    assert (m["rna"].layers["original"] != rna.X).nnz == 0
    np.testing.assert_array_equal(m["rna"].obsm["embedding"], rna.obsm["embedding"])
    assert m.obs.loc["RNAonly", "circ_observation_status"] == "unknown"
    assert not m.obs.loc["RNAonly", "circ_present"]
    assert pd.isna(m.obs.loc["RNAonly", "circRNA_count"])
    assert m.obs.loc["001", "membership"] == "shared"


@pytest.mark.parametrize("support", ["-1", "1.5", "NaN", "inf", "", str(2**63)])
def test_reject_bad_support(inputs, support):
    inputs["calls"].write_text(
        inputs["calls"].read_text().replace("\t7\t", "\t" + support + "\t")
    )
    with pytest.raises(ValueError):
        import_circ(**inputs)
    assert not inputs["outdir"].exists()


@pytest.mark.parametrize(
    "change,match",
    [
        (lambda s: s + s.splitlines()[1] + "\n", "Duplicate cell/candidate"),
        (
            lambda s: s.replace("failed\tNA\tchr1", "failed\tNA\tchr2"),
            "Conflicting candidate",
        ),
        (lambda s: s.replace("001\tNA", "other\tNA"), "absent from cell metadata"),
        (lambda s: s.replace("\t0\t20\t", "\t20\t20\t"), "Invalid interval"),
        (lambda s: s.replace("\t+\t", "\t?\t"), "Invalid interval"),
        (lambda s: s.replace("support\t", "support\tsupport\t"), "duplicate column"),
    ],
)
def test_reject_ambiguous_calls(inputs, change, match):
    inputs["calls"].write_text(change(inputs["calls"].read_text()))
    with pytest.raises(ValueError, match=match):
        import_circ(**inputs)


@pytest.mark.parametrize(
    "column,value",
    [
        ("genome_build", "hg19"),
        ("coordinate_convention", "1-based-closed"),
        ("support_definition", "UMIs"),
        ("detector", "other"),
    ],
)
def test_reject_mixed_declarations(inputs, column, value):
    lines = inputs["calls"].read_text().splitlines()
    inputs["calls"].write_text(
        "\n".join(
            line + "\t" + (column if i == 0 else value) for i, line in enumerate(lines)
        )
        + "\n"
    )
    with pytest.raises(ValueError, match=column):
        import_circ(**inputs)


def test_empty_calls_and_unknown_status(inputs):
    inputs["calls"].write_text(inputs["calls"].read_text().splitlines()[0] + "\n")
    inputs["cells"].write_text("cell_id\n001\nNA\n")
    m = import_circ(**inputs)
    a = ad.read_h5ad(inputs["outdir"] / "anndata/circ_counts.h5ad")
    assert a.shape == (2, 0)
    assert m["detection_fraction_completed"]["value"] is None
    assert m["observation_status_counts"] == {"unknown": 2}


def test_large_support_preserved_exactly(inputs):
    inputs["calls"].write_text(
        inputs["calls"].read_text().replace("\t7\t", "\t9007199254740993\t")
    )
    m = import_circ(**inputs)
    a = ad.read_h5ad(inputs["outdir"] / "anndata/circ_counts.h5ad")
    assert a.X[1, 0] == 9007199254740993
    assert m["total_support"] == 9007199254740995


def test_cli_help_and_error(inputs):
    runner = CliRunner()
    assert runner.invoke(app, ["import-circ", "--help"]).exit_code == 0
    args = ["import-circ"]
    for key, value in inputs.items():
        args.extend(["--" + key.replace("_", "-"), str(value)])
    result = runner.invoke(app, args)
    assert result.exit_code == 0, result.output
    assert runner.invoke(app, args).exit_code == 2


@pytest.mark.parametrize("status", ["success", "", "missing"])
def test_no_implicit_status_mapping(inputs, status):
    inputs["cells"].write_text(
        "cell_id\tdetector_status\n001\t" + status + "\nfailed\tfailed\n"
    )
    with pytest.raises(ValueError, match="detector_status"):
        import_circ(**inputs)


def test_reject_duplicate_cells(inputs):
    inputs["cells"].write_text(inputs["cells"].read_text() + "001\tunknown\tx\n")
    with pytest.raises(ValueError, match="duplicate IDs"):
        import_circ(**inputs)


@pytest.mark.parametrize(
    "problem", ["float", "negative", "duplicate", "shape", "duplicate_index"]
)
def test_invalid_matrix_contract(inputs, problem):
    matrix = inputs["outdir"].with_name("matrix")
    matrix.mkdir()
    (matrix / "cell_index.txt").write_text("001\nfailed\n")
    (matrix / "circ_index.txt").write_text("NA\n")
    (matrix / "circ_feature_table.tsv").write_text(
        "circ_id\tchrom\tstart\tend\tstrand\nNA\tchr1\t0\t20\t+\n"
    )
    raw = sparse.coo_matrix(([7, 2], ([0, 0], [0, 1])), shape=(1, 2))
    if problem == "float":
        raw = raw.astype(float)
    elif problem == "negative":
        raw.data[0] = -1
    elif problem == "shape":
        raw = sparse.coo_matrix((2, 2), dtype=np.int64)
    elif problem == "duplicate_index":
        (matrix / "cell_index.txt").write_text("001\n001\n")
    if problem == "duplicate":
        (matrix / "circ_counts.mtx").write_text(
            "%%MatrixMarket matrix coordinate integer general\n1 2 2\n1 1 7\n1 1 2\n"
        )
    else:
        mmwrite(matrix / "circ_counts.mtx", raw)
    with pytest.raises(ValueError):
        import_circ(**{**inputs, "calls": None, "matrix_dir": matrix})
    assert not inputs["outdir"].exists()


def test_metadata_order_and_square_matrix_orientation(inputs):
    matrix = inputs["outdir"].with_name("matrix")
    matrix.mkdir()
    (matrix / "cell_index.txt").write_text("001\nzero\n")
    (matrix / "circ_index.txt").write_text("a\nb\n")
    (matrix / "circ_feature_table.tsv").write_text(
        "circ_id\tchrom\tstart\tend\tstrand\nb\tchr2\t10\t30\t-\na\tchr1\t0\t20\t+\n"
    )
    mmwrite(
        matrix / "circ_counts.mtx", sparse.csr_matrix([[7, 0], [3, 2]], dtype=np.int64)
    )
    import_circ(**{**inputs, "calls": None, "matrix_dir": matrix})
    a = ad.read_h5ad(inputs["outdir"] / "anndata/circ_counts.h5ad")
    assert a.X.toarray().tolist() == [[0, 2], [7, 3], [0, 0], [0, 0], [0, 0], [0, 0]]
    assert a.var_names.tolist() == ["a", "b"]
    assert a.var.chrom.tolist() == ["chr1", "chr2"]
    assert a.obs.loc["failed", "circ_observation_status"] == "failed"


@pytest.mark.parametrize("problem", ["build", "no_overlap", "duplicate"])
def test_invalid_rna(inputs, problem):
    pytest.importorskip("mudata")
    ids = ["001", "zero"]
    if problem == "no_overlap":
        ids = ["x", "y"]
    if problem == "duplicate":
        ids = ["001", "001"]
    rna = ad.AnnData(
        sparse.csr_matrix([[1], [2]]),
        obs=pd.DataFrame(index=ids),
        var=pd.DataFrame(index=["g"]),
    )
    if problem == "build":
        rna.uns["genome_build"] = "hg19"
    path = inputs["cells"].with_name("rna.h5ad")
    rna.write_h5ad(path)
    with pytest.raises(ValueError):
        import_circ(**inputs, rna=path)
    assert not inputs["outdir"].exists()


def test_html_escapes_metadata_and_no_false_completion(inputs):
    inputs["cells"].write_text(
        inputs["cells"]
        .read_text()
        .replace("001\tcompleted\t001", "001\tunknown\t<script>bad</script>")
    )
    import_circ(**inputs)
    report = (inputs["outdir"] / "qc/report.html").read_text()
    assert "<script>" not in report
    assert "&lt;script&gt;bad&lt;/script&gt;" in report
    assert "0 / 1" in report


def test_total_support_overflow(inputs):
    inputs["calls"].write_text(
        inputs["calls"].read_text().replace("\t7\t", f"\t{2**63 - 1}\t")
    )
    with pytest.raises(ValueError, match="Total support exceeds"):
        import_circ(**inputs)


def test_preserve_explicit_host_provenance(inputs):
    calls = pd.read_csv(inputs["calls"], sep="\t", dtype=str, keep_default_na=False)
    calls["host_gene_source"] = "external_annotation"
    calls["host_gene_from_gtf"] = "different_gene"
    calls.to_csv(inputs["calls"], sep="\t", index=False)
    import_circ(**inputs)
    a = ad.read_h5ad(inputs["outdir"] / "anndata/circ_counts.h5ad")
    assert a.var.loc["NA", "host_gene"] == "GENE1"
    assert a.var.loc["NA", "host_gene_source"] == "external_annotation"
    assert a.var.loc["NA", "host_gene_from_gtf"] == "different_gene"


def test_changed_input_aborts_without_partial_bundle(inputs, monkeypatch):
    import importlib

    module = importlib.import_module("circyto.pipeline.import_circ")
    original = module._load_calls

    def changing_source(path, *args):
        result = original(path, *args)
        path.write_text(path.read_text().replace("\t7\t", "\t8\t"))
        return result

    monkeypatch.setattr(module, "_load_calls", changing_source)
    with pytest.raises(ValueError, match="Source changed"):
        import_circ(**inputs)
    assert not inputs["outdir"].exists()
    assert not list(inputs["outdir"].parent.glob(".circyto-import-*"))


def test_workflow_report_cannot_overwrite_import_qc(inputs):
    import hashlib

    import_circ(**inputs)
    root = inputs["outdir"]
    before = {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in root.rglob("*") if p.is_file()}
    result = CliRunner().invoke(app, ["report", "--workdir", str(root)])
    assert result.exit_code == 1
    assert "import-circ bundle" in result.output
    assert "new output directory" in result.output
    after = {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
             for p in root.rglob("*") if p.is_file()}
    assert after == before
