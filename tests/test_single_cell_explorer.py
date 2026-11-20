from __future__ import annotations

import hashlib
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from typer.testing import CliRunner

from circyto.cli.circyto import app
from circyto.pipeline.single_cell_explorer import ExplorerOptions, build_explorer, load_explorer
from circyto.pipeline.qc_report import generate_report
from circyto.pipeline.explorer_assets import render_explorer_document, SCRIPT, script_policy


@pytest.fixture
def objects():
    circ = ad.AnnData(sparse.csr_matrix([[4, 0], [0, 0], [9, 1], [0, 0], [0, 0], [0, 0]]),
        obs=pd.DataFrame({"detector_status": ["completed", "completed", "failed", "unknown", "unprocessed", "completed"],
                          "condition": ["control", "exposed", "exposed", "control", "exposed", "control"]},
                         index=["A", "B", "F", "U", "P", "C-only"]),
        var=pd.DataFrame({"chrom": ["chr1", "chr2"], "start": [10, 30], "end": [20, 40], "strand": ["+", "-"],
                          "host_gene": ["DUP", "G3"], "circAtlas_id": ["hsa-DEMO_0001", "hsa-DEMO_0002"],
                          "support_definition": ["junction reads", "junction reads"]}, index=["chr1:10|20", "chr2:30|40"]))
    rna = ad.AnnData(sparse.csc_matrix([[0, 2, 3], [5, 0, 4], [7, 8, 0], [1, 1, 1]]),
        obs=pd.DataFrame({"condition": ["exposed", "control", "control", "exposed"],
                          "cluster": ["c2", "c1", "c1", "c2"], "total_counts": [5, 9, 15, 3]}, index=["B", "A", "R-only", "F"]),
        var=pd.DataFrame({"gene_name": ["DUP", "DUP", "G3"], "gene_id": ["ENSG000001.1", "ENSG000002.2", "ENSG000003.1"]},
                         index=["ENSG000001.1", "ENSG000002.2", "ENSG000003.1"]))
    rna.obsm["X_umap"] = np.array([[3., 4.], [1., 2.], [5., 6.], [7., 8.]])
    return circ, rna


def by_cell(payload, feature):
    return dict(zip(payload["cells"], zip(feature["values"], feature["states"])))


def test_shuffled_partial_join_units_and_processing_states(objects):
    circ, rna = objects
    c_before, r_before = circ.copy(), rna.copy()
    payload = build_explorer(circ, rna, options=ExplorerOptions(genes=("DUP",), expression_unit="log1p CP10K"))
    assert payload["join"]["overlap"] == 3
    assert payload["cells"] == ["B", "A", "R-only", "F", "U", "P", "C-only"]
    assert by_cell(payload, payload["circs"][0]) == {
        "A": (4, "completed_detected"), "B": (0, "completed_zero"), "F": (9, "failed"),
        "U": (0, "unknown"), "P": (0, "unprocessed"), "R-only": (None, "missing"), "C-only": (0, "completed_zero")}
    assert payload["coordinates"][payload["cells"].index("A")] == [1., 2.]
    assert payload["coordinates"][-1] is None
    assert len(payload["genes"]) == 2
    assert by_cell(payload, payload["genes"][0])["A"] == (5, "available")
    assert payload["genes"][0]["unit"] == "log1p CP10K"
    assert payload["circs"][0]["unit"] == "junction reads"
    assert (circ.X != c_before.X).nnz == (rna.X != r_before.X).nnz == 0
    pd.testing.assert_frame_equal(circ.obs, c_before.obs)
    np.testing.assert_array_equal(rna.obsm["X_umap"], r_before.obsm["X_umap"])


def test_shuffle_invariance_and_mudata_global_obs(objects):
    from circyto.multimodal.sync import mudata_from_modalities
    circ, rna = objects
    original = build_explorer(circ, rna)
    m = mudata_from_modalities({"rna": rna[::-1].copy(), "circ": circ[::-1].copy()})
    m.obs["cohort"] = "test cohort"
    reordered = build_explorer(m)
    for key in ("genes", "circs"):
        for a, b in zip(original[key], reordered[key]):
            assert by_cell(original, a) == by_cell(reordered, b)
    assert any(m["key"] == "global:cohort" for m in reordered["metadata"])
    # A supplied RNA object supersedes MuData's RNA expression and coordinates as a unit.
    supplied = rna.copy(); supplied.X = supplied.X * 10
    external = build_explorer(m, supplied)
    assert external["genes"][0]["values"][0] == supplied.X[0, 0]


@pytest.mark.parametrize("kind", ["rna", "circ", "rna_feature"])
def test_ambiguous_identities_fail(objects, kind):
    c, r = objects
    if kind == "rna": r.obs_names = ["same"] * r.n_obs
    if kind == "circ": c.obs_names = ["same"] * c.n_obs
    if kind == "rna_feature":
        r.var_names = ["same"] * r.n_vars
        r.var["gene_id"] = "duplicate"
    with pytest.raises(ValueError, match="duplicate identities"):
        build_explorer(c, r)


def test_duplicate_symbol_index_resolves_unique_gene_ids(objects):
    c, r = objects
    r.var_names = ["DUP", "DUP", "G3"]
    p = build_explorer(c, r, options=ExplorerOptions(genes=("DUP",)))
    assert [g["id"] for g in p["genes"]] == ["ENSG000001.1", "ENSG000002.2"]


def test_explicit_biological_ids_and_no_overlap(objects):
    c, r = objects
    c.obs["biological_id"] = c.obs_names
    r.obs["biological_id"] = r.obs_names
    c.obs_names = ["circ_"+str(i) for i in range(c.n_obs)]
    r.obs_names = ["rna_"+str(i) for i in range(r.n_obs)]
    with pytest.raises(ValueError, match="no exact"):
        build_explorer(c, r)
    p = build_explorer(c, r, options=ExplorerOptions(cell_id_key="biological_id"))
    assert p["join"]["overlap"] == 3
    r.obs.loc[r.obs_names[0], "biological_id"] = ""
    with pytest.raises(ValueError, match="nonempty strings"):
        build_explorer(c, r, options=ExplorerOptions(cell_id_key="biological_id"))


def test_missing_umap_rna_empty_circ_and_invalid_coordinates(objects):
    c, r = objects
    del r.obsm["X_umap"]
    p = build_explorer(c, r)
    assert p["embedding"]["plotted_cells"] == 0 and len(p["genes"]) == 3
    assert "UMAP panels omitted" in " ".join(p["warnings"])
    p = build_explorer(c, options=ExplorerOptions(genes=("DUP",)))
    assert p["genes"] == [] and p["join"]["rna_cells"] == 0
    p = build_explorer(c[:, :0].copy(), r)
    assert p["circs"] == []
    r.obsm["X_umap"] = np.zeros((r.n_obs, 3))
    with pytest.raises(ValueError, match="two columns"):
        build_explorer(c, r)
    r.obsm["X_umap"] = np.full((r.n_obs, 2), np.nan)
    assert build_explorer(c, r)["embedding"]["plotted_cells"] == 0


@pytest.mark.parametrize("state", ["unknown", "failed", "unprocessed", "pending", "missing", "new_unrecognized_state"])
def test_zero_with_incomplete_state_is_not_completed_zero(objects, state):
    c, r = objects
    c.obs.loc["B", "detector_status"] = state
    p = build_explorer(c, r)
    assert by_cell(p, p["circs"][0])["B"][1] != "completed_zero"


def test_failure_overrides_completion_nan_is_missing_and_raw_values_retained(objects):
    c, r = objects
    c.obs["circ_observation_status"] = "completed_zero"
    c.obs["alignment_status"] = "aligned"
    c.obs.loc["B", "alignment_status"] = "failed"
    c.X = c.X.astype(float)
    c.X[0, 0] = np.nan
    p = build_explorer(c, r)
    assert by_cell(p, p["circs"][0])["A"] == (None, "missing")
    assert by_cell(p, p["circs"][0])["B"] == (0, "failed")
    assert by_cell(p, p["circs"][0])["F"] == (9, "failed")


def test_conflicting_unknown_status_is_not_promoted_to_completed(objects):
    c, r = objects
    c.obs["circ_observation_status"] = "completed_detected"
    c.obs.loc["A", "detector_status"] = None
    r.obs["rna_status"] = "completed"
    r.obs.loc["A", "rna_status"] = None
    p = build_explorer(c, r)
    assert by_cell(p, p["circs"][0])["A"][1] == "unknown"
    assert by_cell(p, p["circs"][0])["U"][1] == "unknown"
    assert by_cell(p, p["genes"][0])["A"][1] == "unknown"


@pytest.mark.parametrize("selector", ["chr1:10|20", "chr1:10-20", "chr1:10|20:+", "hsa-DEMO_0001"])
def test_exact_coordinate_and_catalogue_alias_search(objects, selector):
    p = build_explorer(*objects, options=ExplorerOptions(circs=(selector,), genes=("ENSG000001",)))
    assert len(p["circs"]) == len(p["genes"]) == 1
    assert p["circs"][0]["id"] == "chr1:10|20"


def test_uniform_alias_exon_commas_and_explicit_cohort(objects):
    c, r = objects
    c.var["uniform_id"] = ["circGENE(2,3,4).1", "circOTHER(1,2).1"]
    c.var["host_gene_source"] = "CIRI3_gene_id"
    p = build_explorer(c, r, options=ExplorerOptions(circs=("circGENE(2,3,4).1",), cells=("F", "A", "C-only")))
    assert p["cells"] == ["F", "A", "C-only"] and p["join"]["overlap"] == 2
    assert p["circs"][0]["values"] == [9, 4, 0]
    assert "CIRI3_gene_id" not in p["circs"][0]["aliases"]
    assert "3" not in p["circs"][0]["aliases"]
    assert p["selection"]["cells_total"] == 7
    with pytest.raises(ValueError, match="absent"):
        build_explorer(c, r, options=ExplorerOptions(cells=("absent",)))


@pytest.mark.parametrize("options,match", [
    (ExplorerOptions(max_genes=1, genes=("DUP",)), "duplicate symbols"),
    (ExplorerOptions(genes=("absent",)), "no exact"),
    (ExplorerOptions(max_cells=2), "max-cells"),
    (ExplorerOptions(max_values=2), "max-values"),
    (ExplorerOptions(max_bytes=100), "max-bytes"),
    (ExplorerOptions(rna_layer="absent"), "layer"),
    (ExplorerOptions(obs_columns=("rna:absent",)), "observation columns")])
def test_budget_and_option_guards(objects, options, match):
    with pytest.raises(ValueError, match=match):
        build_explorer(*objects, options=options)


@pytest.mark.parametrize("bad", [-1, 0.5, 2**53])
def test_invalid_or_inexact_native_support_rejected(objects, bad):
    c, r = objects
    c.X = c.X.astype(float); c.X[0, 0] = bad
    with pytest.raises(ValueError, match="native junction|exact numeric"):
        build_explorer(c, r)


def test_complex_measurements_are_rejected(objects):
    c, r = objects
    r.X = r.X.astype(complex)
    with pytest.raises(ValueError, match="real numeric"):
        build_explorer(c, r)


def test_large_sparse_features_only_selected_columns_materialized(monkeypatch):
    n = 100_000
    rna = ad.AnnData(sparse.csr_matrix(([3., 7.], ([0, 39], [n-1, 0])), shape=(40, n)),
                     obs=pd.DataFrame(index=[f"c{i}" for i in range(40)]),
                     var=pd.DataFrame(index=[f"G{i:06d}" for i in range(n)]))
    original = sparse.csr_matrix.toarray
    def bounded(self, *args, **kwargs):
        assert self.shape[1] <= 1, "Full feature matrix densified"
        return original(self, *args, **kwargs)
    monkeypatch.setattr(sparse.csr_matrix, "toarray", bounded)
    p = build_explorer(None, rna, options=ExplorerOptions(genes=("G099999",), max_genes=1))
    assert p["genes"][0]["values"][0] == 3
    assert p["selection"]["genes_total"] == n
    assert len(json.dumps(p)) < 20000


def test_metadata_conflicts_links_and_script_injection(objects):
    c, r = objects
    evil = '</script><script>window.PWNED=1</script>'
    c.var["host_gene"] = [evil, "G3"]
    c.var["evidence_url"] = ["javascript:alert(1)", "https://example.org/evidence"]
    p = build_explorer(c, r)
    document = render_explorer_document(p)
    assert evil not in document and "\\u003c/script" in document
    assert p["circs"][0]["evidence_links"] == []
    assert p["circs"][1]["evidence_links"][0]["url"] == "https://example.org/evidence"
    assert "connect-src 'none'" in script_policy()
    assert "fetch(" not in SCRIPT and "XMLHttpRequest" not in SCRIPT
    sidecar = pd.DataFrame({"strand": ["-"]}, index=[c.var_names[0]])
    with pytest.raises(ValueError, match="conflicts"):
        build_explorer(c, r, candidate_metadata=sidecar)


def test_backed_cli_replay_preserves_inputs_and_default_report(objects, tmp_path, monkeypatch):
    c, r = objects
    cpath, rpath = tmp_path/'circ.h5ad', tmp_path/'rna.h5ad'
    c.write_h5ad(cpath); r.write_h5ad(rpath)
    before = [hashlib.sha256(p.read_bytes()).hexdigest() for p in (cpath, rpath)]
    import subprocess
    monkeypatch.setattr(subprocess, "Popen", lambda *a, **k: pytest.fail("reporting launched a pipeline"))
    runner = CliRunner()
    args = ["report", "--explorer", str(cpath), "--rna", str(rpath), "--output-dir", str(tmp_path/'report'), "--gene", "DUP"]
    result = runner.invoke(app, args)
    assert result.exit_code == 0, result.output
    html = (tmp_path/'report/report.html').read_bytes()
    p = json.loads((tmp_path/'report/explorer.json').read_text())
    assert len(p["genes"]) == 2 and p["sources"][0]["sha256"] == before[0]
    assert runner.invoke(app, args).exit_code == 0
    assert (tmp_path/'report/report.html').read_bytes() == html
    assert [hashlib.sha256(p.read_bytes()).hexdigest() for p in (cpath, rpath)] == before
    assert not (tmp_path/'report/metrics.json').exists()
    default = tmp_path/'workflow'; default.mkdir()
    original = generate_report(default)
    assert 'explorer' not in original and '<script' not in (default/'qc/report.html').read_text()
    generate_report(default, explorer=p, output_dir=tmp_path/'combined')
    assert 'id="explorer"' in (tmp_path/'combined/report.html').read_text()
    assert (tmp_path/'combined/metrics.json').exists()


def test_backed_mudata_layers_and_file_handles(objects, tmp_path):
    from circyto.multimodal.sync import mudata_from_modalities, write_h5mu
    c, r = objects
    r.layers["log"] = r.X.copy() * 0.1
    m = mudata_from_modalities({"circ": c, "rna": r})
    path = tmp_path/'multi.h5mu'; write_h5mu(m, path)
    p = load_explorer(path, options=ExplorerOptions(rna_layer="log", expression_unit="log values"))
    assert p["genes"][0]["values"][1] == 0.5
    assert p["genes"][0]["unit"] == "log values"


def test_output_guards_no_partial_writes(objects, tmp_path):
    p = build_explorer(*objects)
    target = tmp_path/'source'; target.write_text('native evidence')
    out = tmp_path/'out'; out.mkdir()
    (out/'report.html').symlink_to(target)
    with pytest.raises(ValueError, match="unsafe report"):
        generate_report(None, explorer=p, output_dir=out)
    assert not (out/'explorer.json').exists() and target.read_text() == 'native evidence'
    (out/'report.html').unlink()
    (out/'metrics.json').write_text('{"schema_version":"circyto.import_circ.v1"}')
    with pytest.raises(ValueError, match="import-circ bundle"):
        generate_report(None, explorer=p, output_dir=out)
    assert not (out/'report.html').exists()
    (out/'metrics.json').unlink()
    (out/'native').mkdir()
    (out/'explorer.json').write_text('source cell identities')
    with pytest.raises(ValueError, match="source artifact resolves"):
        generate_report(None, explorer=p, output_dir=out/'native'/'..', protected_sources=(out/'explorer.json',))
    assert (out/'explorer.json').read_text() == 'source cell identities'
