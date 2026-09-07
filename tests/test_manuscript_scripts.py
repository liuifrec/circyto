from __future__ import annotations

import shutil
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from circyto.multimodal.sync import mudata_from_modalities, read_h5mu, write_h5mu


SCRIPT_DIR = Path("scripts/manuscript")


@pytest.fixture
def regeneration():
    spec = importlib.util.spec_from_file_location(
        "manuscript_regeneration", SCRIPT_DIR / "regenerate_application_note_results.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("mismatched", ["smartseq3", "imr90"])
def test_regeneration_checksum_failure_precedes_reads_and_writes(tmp_path, monkeypatch, regeneration, mismatched):
    paths = {key: tmp_path / f"{key}.h5mu" for key in ("smartseq3", "imr90")}
    for key, path in paths.items():
        path.write_bytes(b"checksum fixture")
        monkeypatch.setattr(regeneration, "SMARTSEQ3_SHA256" if key == "smartseq3" else "IMR90_SHA256",
                            "0" * 64 if key == mismatched else regeneration.sha256_file(path))
    monkeypatch.setattr(regeneration, "read_modalities", lambda *_: pytest.fail("Must checksum both inputs before reading"))
    output = tmp_path / "manuscript" / "results"
    with pytest.raises(SystemExit, match="SHA-256 mismatch"):
        regeneration.main(["--smartseq3", str(paths["smartseq3"]), "--imr90", str(paths["imr90"]), "--output", str(output)])
    assert not output.exists()
    assert all(path.read_bytes() == b"checksum fixture" for path in paths.values())


def test_regeneration_output_scope_and_existing_directory(tmp_path, monkeypatch, regeneration):
    monkeypatch.setattr(regeneration, "REPO_ROOT", tmp_path)
    root = tmp_path / "manuscript"
    root.mkdir()
    outside = tmp_path / "outside"
    outside.mkdir()
    (root / "escape").symlink_to(outside, target_is_directory=True)
    for path in (root, outside / "results", root / "escape" / "results", root / ".." / "results"):
        with pytest.raises(SystemExit, match="Output"):
            regeneration.validate_output_directory(path, [])
    assert regeneration.validate_output_directory(root / "results", []) == root / "results"


def test_regeneration_sparse_stored_zeros_and_categorical_hosts(regeneration):
    from scipy import sparse
    from types import SimpleNamespace

    matrix = sparse.csr_matrix(([0, 2, 3], ([0, 0, 1], [0, 1, 0])), shape=(2, 2))
    assert regeneration.nonzero_per_row(matrix).tolist() == [1, 1]
    assert matrix.nnz == 3
    circ = SimpleNamespace(var=pd.DataFrame({"host_gene": pd.Categorical(["MAN1A2", None, "", "NA"])}))
    assert regeneration.nonempty_host_gene_count(circ) == 1


def test_regeneration_aligns_cells_and_publishes_disagreements(tmp_path, monkeypatch, regeneration):
    script_copy = tmp_path / "scripts" / "manuscript" / "regenerate_application_note_results.py"
    script_copy.parent.mkdir(parents=True)
    shutil.copyfile(regeneration.__file__, script_copy)
    monkeypatch.setattr(regeneration, "__file__", str(script_copy))
    monkeypatch.setattr(regeneration, "REPO_ROOT", tmp_path)
    source = tmp_path / "smart.h5mu"
    other = tmp_path / "imr.h5mu"
    _write_tiny_mudata(source, include_rt=False)
    data = read_h5mu(source)
    data.mod["circ"].var_names = [regeneration.MAN1A2_CIRC_ID, "circ2", "circ3", "circ4"]
    data.mod["circ"].var["host_gene"] = ["MAN1A2", "FN1", "VIM", ""]
    data.mod["circ"] = data.mod["circ"][list(reversed(data.mod["circ"].obs_names)), :].copy()
    write_h5mu(mudata_from_modalities(dict(data.mod)), source)
    shutil.copyfile(source, other)
    before = [regeneration.sha256_file(path) for path in (source, other)]
    monkeypatch.setattr(regeneration, "SMARTSEQ3_SHA256", before[0])
    monkeypatch.setattr(regeneration, "IMR90_SHA256", before[1])
    monkeypatch.setattr(regeneration, "repository_provenance", lambda: {"production_diff_from_baseline": ""})
    monkeypatch.setattr(regeneration, "compute_rna_umap", lambda rna, cells: (np.arange(10).reshape(5, 2), {"source_modality": "rna"}))
    monkeypatch.setattr(regeneration, "write_figures", lambda *_: [])
    output = tmp_path / "manuscript" / "results"
    assert regeneration.main(["--smartseq3", str(source), "--imr90", str(other), "--output", str(output)]) == 0
    cells = pd.read_csv(output / "smartseq3_umap_cells.tsv", sep="\t")
    candidate = pd.read_csv(output / "smartseq3_selected_candidate.tsv", sep="\t")
    assert cells["cell_id"].tolist() == [f"cell{i}" for i in range(1, 6)]
    assert candidate["man1a2_candidate_support"].tolist() == [1, 0, 1, 2, 3]
    pd.testing.assert_frame_equal(cells[["cell_id", "UMAP1", "UMAP2"]], candidate[["cell_id", "UMAP1", "UMAP2"]])
    comparison = pd.read_csv(output / "evidence_comparison.tsv", sep="\t")
    assert not comparison["match"].all()
    assert {"quantity", "historical_value", "regenerated_value", "match", "source_checksum", "evidence_class"} <= set(comparison)
    assert json.loads((output / "provenance.json").read_text())["all_evidence_values_agree"] is False
    assert [regeneration.sha256_file(path) for path in (source, other)] == before


def test_regeneration_rejects_duplicate_cell_ids(regeneration):
    from types import SimpleNamespace

    with pytest.raises(SystemExit, match="Duplicate cell IDs"):
        regeneration.shared_obs_in_left_order(SimpleNamespace(obs_names=pd.Index(["cell1", "cell1"])))


def _run_script(args: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, *args],
        cwd=Path.cwd(),
        check=False,
        capture_output=True,
        text=True,
    )


def _write_tiny_mudata(path: Path, *, include_rt: bool = True, include_cnv: bool = True) -> None:
    mu = pytest.importorskip("mudata")
    ad = pytest.importorskip("anndata")

    cells = [f"cell{i}" for i in range(1, 6)]
    rna = ad.AnnData(
        X=np.array(
            [
                [1, 0, 0, 0, 0],
                [1, 1, 0, 0, 0],
                [1, 1, 1, 0, 0],
                [1, 1, 1, 1, 0],
                [1, 1, 1, 1, 1],
            ],
            dtype=float,
        ),
        obs=pd.DataFrame({"detected_genes": [1, 2, 3, 4, 5]}, index=cells),
        var=pd.DataFrame({"gene_name": ["G1", "G2", "G3", "G4", "G5"]}, index=[f"G{i}" for i in range(1, 6)]),
    )
    circ = ad.AnnData(
        X=np.array(
            [
                [1, 0, 0, 0],
                [0, 2, 0, 0],
                [1, 1, 0, 0],
                [2, 0, 1, 0],
                [3, 1, 1, 0],
            ],
            dtype=float,
        ),
        obs=pd.DataFrame(index=cells),
        var=pd.DataFrame(
            {
                "host_gene": ["COL1A1", "FN1", "VIM", ""],
                "host_gene_source": ["gtf", "gtf", "circatlas", ""],
                "host_gene_from_gtf": ["COL1A1", "FN1", "", ""],
                "host_gene_from_circatlas": ["", "", "VIM", ""],
                "host_gene_from_circatlas_id": ["", "", "", ""],
                "host_gene_id": ["ENSG1", "ENSG2", "ENSG3", ""],
                "host_genes_multi": ["COL1A1", "FN1", "VIM", ""],
                "host_gene_ids_multi": ["ENSG1", "ENSG2", "ENSG3", ""],
                "host_gene_n": [1, 1, 1, 0],
                "chrom": ["chr1", "chr1", "chr2", "chr3"],
                "start": [10, 150, 20, 1],
                "end": [20, 180, 40, 10],
                "annotation_status": ["known", "novel", "known", "novel"],
            },
            index=["circ1", "circ2", "circ3", "circ4"],
        ),
    )
    modalities = {"rna": rna, "circ": circ}
    if include_rt:
        modalities["rt"] = ad.AnnData(
            X=np.array(
                [
                    [-1, -1, 0],
                    [-1, 0, 1],
                    [0, 1, 1],
                    [1, 1, 1],
                    [1, 0, 1],
                ],
                dtype=float,
            ),
            obs=pd.DataFrame(index=cells),
            var=pd.DataFrame(
                {"seqname": ["chr1", "chr1", "chr2"], "start": [0, 100, 0], "end": [100, 200, 100]},
                index=["chr1:0-100", "chr1:100-200", "chr2:0-100"],
            ),
        )
    if include_cnv:
        modalities["cnv"] = ad.AnnData(
            X=np.array(
                [
                    [2, 2, 2],
                    [1, 2, 2],
                    [2, 3, 2],
                    [3, 3, 2],
                    [1, 3, 3],
                ],
                dtype=float,
            ),
            obs=pd.DataFrame(index=cells),
            var=pd.DataFrame(
                {"seqname": ["chr1", "chr1", "chr2"], "start": [0, 100, 0], "end": [100, 200, 100]},
                index=["chr1:0-100", "chr1:100-200", "chr2:0-100"],
            ),
        )
    write_h5mu(mudata_from_modalities(modalities), path)


def test_manuscript_script_help_runs_without_mudata() -> None:
    for script_name in [
        "summarize_mudata_inventory.py",
        "hap1_rt_circ_analysis.py",
        "imr90_cnv_circ_analysis.py",
        "cross_dataset_host_overlap.py",
        "known_novel_circ_summary.py",
        "regenerate_application_note_results.py",
    ]:
        result = _run_script([str(SCRIPT_DIR / script_name), "--help"])
        assert result.returncode == 0, result.stderr + result.stdout
        assert "usage:" in result.stdout.lower()


def test_manuscript_inventory_hap1_imr90_and_known_novel_scripts(tmp_path: Path) -> None:
    h5mu = tmp_path / "tiny.hostgene_fixed.h5mu"
    _write_tiny_mudata(h5mu)

    inventory = tmp_path / "inventory.tsv"
    result = _run_script(
        [
            str(SCRIPT_DIR / "summarize_mudata_inventory.py"),
            str(h5mu),
            "--dataset-name",
            "tiny",
            "--out",
            str(inventory),
        ]
    )
    assert result.returncode == 0, result.stderr + result.stdout
    inventory_df = pd.read_csv(inventory, sep="\t")
    assert {"dataset", "modalities", "circRNA_count", "host_gene_recovery"}.issubset(inventory_df.columns)
    assert inventory_df.loc[0, "circRNA_count"] == 4
    assert inventory_df.loc[0, "host_gene_annotated"] == 3

    hap1_corr = tmp_path / "hap1_correlations.tsv"
    hap1_ols = tmp_path / "hap1_ols.tsv"
    hap1_cells = tmp_path / "hap1_cells.tsv"
    result = _run_script(
        [
            str(SCRIPT_DIR / "hap1_rt_circ_analysis.py"),
            "--input",
            str(h5mu),
            "--correlations-out",
            str(hap1_corr),
            "--ols-out",
            str(hap1_ols),
            "--scatter-out",
            str(hap1_cells),
        ]
    )
    assert result.returncode == 0, result.stderr + result.stdout
    assert {"x", "y", "pearson_r", "spearman_r"}.issubset(pd.read_csv(hap1_corr, sep="\t").columns)
    assert "frac_rt_pos" in set(pd.read_csv(hap1_ols, sep="\t")["term"])
    assert {"cell_id", "circRNA_count", "detected_genes", "frac_rt_pos"}.issubset(
        pd.read_csv(hap1_cells, sep="\t").columns
    )

    imr_corr = tmp_path / "imr90_correlations.tsv"
    imr_enrich = tmp_path / "imr90_enrichment.tsv"
    imr_local = tmp_path / "imr90_local.tsv"
    result = _run_script(
        [
            str(SCRIPT_DIR / "imr90_cnv_circ_analysis.py"),
            "--input",
            str(h5mu),
            "--correlations-out",
            str(imr_corr),
            "--enrichment-out",
            str(imr_enrich),
            "--local-cnv-out",
            str(imr_local),
        ]
    )
    assert result.returncode == 0, result.stderr + result.stdout
    assert {"frac_non_diploid", "frac_loss", "frac_gain"} & set(pd.read_csv(imr_corr, sep="\t")["x"])
    assert "COL1A1" in set(pd.read_csv(imr_enrich, sep="\t")["host_gene"])
    local_df = pd.read_csv(imr_local, sep="\t")
    assert "status" in local_df.columns
    assert "ok" in set(local_df["status"])

    known_out = tmp_path / "known_novel.tsv"
    result = _run_script(
        [
            str(SCRIPT_DIR / "known_novel_circ_summary.py"),
            str(h5mu),
            "--dataset-name",
            "tiny",
            "--out",
            str(known_out),
        ]
    )
    assert result.returncode == 0, result.stderr + result.stdout
    known_df = pd.read_csv(known_out, sep="\t")
    assert {"known", "novel", "all"}.issubset(set(known_df["known_status"]))


def test_known_novel_script_handles_categorical_status_with_missing_values(tmp_path: Path) -> None:
    mu = pytest.importorskip("mudata")
    h5mu = tmp_path / "categorical_status.h5mu"
    _write_tiny_mudata(h5mu)
    mdata = read_h5mu(h5mu)
    mdata.mod["circ"].var["annotation_status"] = pd.Categorical(
        ["known", "novel", None, "database_match"]
    )
    write_h5mu(mdata, h5mu)

    known_out = tmp_path / "known_novel.tsv"
    result = _run_script(
        [
            str(SCRIPT_DIR / "known_novel_circ_summary.py"),
            str(h5mu),
            "--out",
            str(known_out),
        ]
    )

    assert result.returncode == 0, result.stderr + result.stdout
    known_df = pd.read_csv(known_out, sep="\t").set_index("known_status")
    assert known_df.loc["known", "n_circRNAs"] == 2
    assert known_df.loc["novel", "n_circRNAs"] == 1
    assert known_df.loc["unknown", "n_circRNAs"] == 1


def test_cross_dataset_host_overlap_script(tmp_path: Path) -> None:
    h5mu = tmp_path / "tiny.hostgene_fixed.h5mu"
    _write_tiny_mudata(h5mu)
    copy_a = tmp_path / "tiny_a.h5mu"
    copy_b = tmp_path / "tiny_b.h5mu"
    shutil.copyfile(h5mu, copy_a)
    shutil.copyfile(h5mu, copy_b)
    enrichment = tmp_path / "enrichment.tsv"
    enrichment.write_text(
        "host_gene\tlog2_fc_pseudocount1\nCOL1A1\t1.5\nFN1\t0.2\nVIM\t-0.1\n",
        encoding="utf-8",
    )

    outdir = tmp_path / "overlap"
    result = _run_script(
        [
            str(SCRIPT_DIR / "cross_dataset_host_overlap.py"),
            str(h5mu),
            str(copy_a),
            str(copy_b),
            "--dataset-name",
            "smart",
            "--dataset-name",
            "hap1",
            "--dataset-name",
            "imr90",
            "--outdir",
            str(outdir),
            "--enrichment-table",
            f"HAP1={enrichment}",
            "--enrichment-table",
            f"IMR90={enrichment}",
        ]
    )
    assert result.returncode == 0, result.stderr + result.stdout
    pairwise = pd.read_csv(outdir / "pairwise_host_gene_overlap.tsv", sep="\t")
    three_way = pd.read_csv(outdir / "three_way_host_gene_overlap.tsv", sep="\t")
    program = pd.read_csv(outdir / "shared_positive_host_gene_program.tsv", sep="\t")
    assert pairwise["n_overlap"].min() == 3
    assert {"COL1A1", "FN1", "VIM"}.issubset(set(three_way["host_gene"]))
    assert {"COL1A1", "FN1"}.issubset(set(program["host_gene"]))


def test_hap1_script_missing_rt_modality_error_is_informative(tmp_path: Path) -> None:
    h5mu = tmp_path / "missing_rt.h5mu"
    _write_tiny_mudata(h5mu, include_rt=False)
    result = _run_script(
        [
            str(SCRIPT_DIR / "hap1_rt_circ_analysis.py"),
            "--input",
            str(h5mu),
            "--correlations-out",
            str(tmp_path / "corr.tsv"),
            "--ols-out",
            str(tmp_path / "ols.tsv"),
        ]
    )
    assert result.returncode != 0
    assert "Missing required modality" in (result.stderr + result.stdout)
    assert "rt" in (result.stderr + result.stdout)


def test_candidate_signal_terminology_does_not_use_banned_phrase() -> None:
    forbidden_variant_phrase = " ".join(["validated", "somatic", "mutation"])
    forbidden_exposure_phrase = " ".join(["radiation", "dataset"])
    roots = [Path("README.md"), Path("docs"), Path("manuscript"), Path("scripts"), Path("circyto"), Path("tests")]
    checked = 0
    for root in roots:
        files = [root] if root.is_file() else [path for path in root.rglob("*") if path.is_file()]
        for path in files:
            if path.suffix.lower() not in {".md", ".py", ".sh", ".txt", ".toml"}:
                continue
            text = path.read_text(encoding="utf-8", errors="ignore").lower()
            assert forbidden_variant_phrase not in text, f"Forbidden variant phrase found in {path}"
            assert forbidden_exposure_phrase not in text, f"Overclaim-prone exposure phrase found in {path}"
            checked += 1
    assert checked > 0
