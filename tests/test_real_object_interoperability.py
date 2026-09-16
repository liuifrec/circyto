"""Small regressions for the audit; actual-object evidence comes from its CLI."""
from __future__ import annotations

import importlib.util
from pathlib import Path
import shutil

import anndata as ad
import h5py
import mudata as mu
import numpy as np
import pandas as pd
import pytest
from scipy import sparse


@pytest.fixture
def audit():
    path = Path("scripts/manuscript/verify_real_object_interoperability.py")
    spec = importlib.util.spec_from_file_location("real_object_audit", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def named_source(tmp_path):
    index = pd.Index(["cell_b", "cell_a"], name="cell_id")
    rna = ad.AnnData(
        X=sparse.csr_matrix([[2, 0], [0, 3]], dtype=np.int32),
        obs=pd.DataFrame({"batch": pd.Categorical(["b", None])}, index=index),
        var=pd.DataFrame({"gene": ["G2", "G1"]}, index=["g2", "g1"]),
    )
    rna.layers["counts"] = rna.X.copy()
    rna.obsm["embedding"] = np.array([[0, 1], [2, 3]], dtype=np.float32)
    rna.varm["loadings"] = np.array([[3], [4]], dtype=np.float64)
    rna.obsp["graph"] = sparse.eye(2, format="csr", dtype=np.float64)
    rna.varp["graph"] = sparse.eye(2, format="csr", dtype=np.float64)
    rna.uns["nested"] = {"array": np.array([1, 2]), "label": "retain"}
    rna.raw = rna.copy()
    circ = ad.AnnData(
        X=sparse.csr_matrix([[0, 5], [7, 0]], dtype=np.int32),
        obs=rna.obs.copy(),
        var=pd.DataFrame({"host_gene": pd.Categorical(["G1", None])}, index=["c2", "c1"]),
    )
    with mu.set_options(pull_on_update=True):
        obj = mu.MuData({"rna": rna, "circ": circ})
    obj.obs = pd.DataFrame({"global_quality": [0.5, np.nan]}, index=index.rename(None))
    obj.uns["note"] = {"scope": "global"}
    obj.obsm["global_embedding"] = np.array([[1.0], [np.nan]])
    obj.varm["global_loading"] = np.arange(4, dtype=np.float32).reshape(-1, 1)
    obj.obsp["global_graph"] = sparse.eye(2, format="csr")
    obj.varp["global_graph"] = sparse.eye(4, format="csr")
    path = tmp_path / "named.h5mu"
    with pytest.warns(FutureWarning):
        obj.write_h5mu(path)
    return path


def test_native_roundtrip_and_separate_named_index_helper_failure(audit, named_source, tmp_path):
    original_hash = audit.sha256(named_source)
    with pytest.warns(FutureWarning) as caught:
        result = audit.audit_object(named_source, tmp_path / "roundtrip.h5mu")
    assert result["native_read"]["status"] == "PASS"
    assert result["native_reread"]["status"] == "PASS"
    assert result["circyto_helper_read"]["status"] == "FAIL"
    assert "cannot join with no overlapping index names" in result["circyto_helper_read"]["traceback"]
    assert str(Path.cwd()) not in result["circyto_helper_read"]["traceback"]
    assert len(result["native_read"]["warnings"]) == 2
    assert len(caught) == 6  # Native read, write and reread each report two.
    assert result["hdf5_roundtrip"]["datasets_checked"] > 20
    assert audit.sha256(named_source) == original_hash
    diagnostic = result["name_only_diagnostic"]
    assert diagnostic["helper_read"]["status"] == "PASS"
    assert diagnostic["change"] == {
        "global_obs_index_name_before": None, "global_obs_index_name_after": "cell_id",
    }
    assert audit.sha256(named_source) == original_hash


@pytest.mark.parametrize("mismatched", ["smartseq3", "imr90"])
def test_both_checksums_precede_reads_and_writes(audit, monkeypatch, tmp_path, mismatched):
    paths = {label: tmp_path / f"{label}.h5mu" for label in audit.INPUTS}
    for path in paths.values():
        path.write_bytes(b"not HDF5; must never reach a reader")
    monkeypatch.setattr(audit, "INPUTS", {
        label: (path.name, "0" * 64 if label == mismatched else audit.sha256(path))
        for label, path in paths.items()
    })
    monkeypatch.setattr(audit, "audit_object", lambda *_: pytest.fail("Reader reached"))
    output = tmp_path / "output"
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        audit.main(["--smartseq3", str(paths["smartseq3"]), "--imr90", str(paths["imr90"]),
                    "--output", str(output)])
    assert not output.exists()


@pytest.mark.parametrize("target", ["existing", "checkout", "symlink"])
def test_output_guard_prevents_overwrite_and_checkout_artifacts(audit, monkeypatch, tmp_path, target):
    checkout = tmp_path / "repo"
    checkout.mkdir()
    link = tmp_path / "link"
    link.symlink_to(checkout, target_is_directory=True)
    monkeypatch.setattr(audit, "REPO_ROOT", checkout)
    monkeypatch.setattr(audit, "verify_inputs", lambda _: None)
    output = {"existing": tmp_path, "checkout": checkout / "new", "symlink": link / "new"}[target]
    with pytest.raises(ValueError, match="NEW directory outside"):
        audit.main(["--smartseq3", "unused", "--imr90", "unused", "--output", str(output)])


@pytest.mark.parametrize("change", ["matrix", "cell_order", "missingness", "layer", "uns", "map"])
def test_equivalence_detects_semantic_changes(audit, named_source, change):
    with pytest.warns(FutureWarning):
        original = mu.read_h5mu(named_source)
    with pytest.warns(FutureWarning):
        changed = mu.read_h5mu(named_source)
    if change == "matrix":
        # Same shape, nnz, support total; different per-cell values.
        changed.mod["rna"].X.data[:] = [3, 2]
    elif change == "cell_order":
        changed.mod["rna"].obs_names = changed.mod["rna"].obs_names[::-1]
    elif change == "missingness":
        changed.mod["rna"].obs.iloc[1, 0] = "b"
    elif change == "layer":
        changed.mod["rna"].layers["counts"].data[0] += 1
    elif change == "uns":
        changed.uns["note"]["scope"] = "lost"
    elif change == "map":
        changed.obsmap["rna"] = np.array([2, 1], dtype=np.uint32)
    with pytest.raises(AssertionError):
        audit.assert_native_equal(original, changed)


def test_hdf5_audit_rejects_unlisted_attribute_loss(audit, named_source, tmp_path):
    changed = tmp_path / "attribute_loss.h5mu"
    shutil.copyfile(named_source, changed)
    with h5py.File(changed, "r+") as handle:
        del handle["mod/rna"].attrs["encoding-version"]
    with pytest.raises(AssertionError, match="attributes changed"):
        audit.compare_hdf5(named_source, changed)
