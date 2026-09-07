#!/usr/bin/env python
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

# Always resolve the frozen package from this checkout, even when an unrelated
# circyto distribution happens to be installed in the analysis environment.
REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
sys.dont_write_bytecode = True
for _thread_variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[_thread_variable] = "1"

import numpy as np
import pandas as pd
from scipy import sparse

SMARTSEQ3_SHA256 = "0ecd36bb0a74455db7f0affb9ade5023c1934c1dd234aca975365c0b69d8b339"
IMR90_SHA256 = "bb2e12f7c3b36f9fa72d98cd71e8bea905a67f50e22af1d6b713550ee92b60c8"
MAN1A2_CIRC_ID = "chr1:117402186|117420649"
UMAP_SEED = 17


def read_modalities(path: Path) -> dict:
    """Read stored AnnData groups without synchronizing legacy global metadata.

    MuData 0.3.10 cannot synchronize these historical global obs indices with
    pandas 2.3.3. No global attributes are used for manuscript quantities.
    h5py mode 'r' and AnnData's public reader preserve the archived bytes.
    """
    import h5py
    from anndata.io import read_elem

    with h5py.File(path, "r") as handle:
        return {name: read_elem(group) for name, group in handle["mod"].items()}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_input(path: Path, expected: str, label: str) -> dict[str, object]:
    if not path.is_file():
        raise SystemExit(f"[ERROR] {label} object does not exist: {path}")
    observed = sha256_file(path)
    if observed != expected:
        raise SystemExit(
            f"[ERROR] {label} SHA-256 mismatch; expected {expected}, observed {observed}. "
            "No manuscript outputs were generated."
        )
    return {
        "dataset": label,
        "path": str(path),
        "expected_sha256": expected,
        "observed_sha256": observed,
        "checksum_match": True,
        "size_bytes": path.stat().st_size,
    }


def nonzero_per_row(matrix) -> np.ndarray:
    if sparse.issparse(matrix):
        # Stored sparse zeros are not detections; do not mutate the input.
        return np.asarray((matrix != 0).sum(axis=1)).ravel().astype(int)
    return np.count_nonzero(np.asarray(matrix), axis=1).astype(int)


def sum_per_row(matrix) -> np.ndarray:
    return np.asarray(matrix.sum(axis=1)).ravel().astype(float)


def column_values(matrix, index: int) -> np.ndarray:
    column = matrix[:, index]
    if sparse.issparse(column):
        return column.toarray().ravel().astype(float)
    return np.asarray(column).ravel().astype(float)


def nonempty_host_gene_count(circ) -> int:
    if "host_gene" not in circ.var.columns:
        return 0
    values = circ.var["host_gene"].astype(object).fillna("").astype(str).str.strip()
    valid = (values != "") & (~values.str.lower().isin({"nan", "none", "na"}))
    return int(valid.sum())


def shared_obs_in_left_order(*modalities) -> list[str]:
    if any(not modality.obs_names.is_unique for modality in modalities):
        raise SystemExit("[ERROR] Duplicate cell IDs; cannot align modalities safely.")
    shared = set(map(str, modalities[0].obs_names))
    for modality in modalities[1:]:
        shared &= set(map(str, modality.obs_names))
    return [name for name in map(str, modalities[0].obs_names) if name in shared]


def compute_rna_umap(rna, shared_cells: list[str]) -> tuple[np.ndarray, dict[str, object]]:
    try:
        import scanpy as sc
    except ImportError as exc:
        raise SystemExit("[ERROR] scanpy is required to regenerate the RNA-derived UMAP.") from exc

    adata = rna[shared_cells, :].copy()
    adata.var_names_make_unique()
    sc.pp.normalize_total(adata, target_sum=10_000)
    sc.pp.log1p(adata)
    n_top_genes = min(2_000, adata.n_vars)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor="seurat")
    n_comps = min(50, max(2, adata.n_obs - 1), max(2, adata.n_vars - 1))
    sc.tl.pca(
        adata,
        n_comps=n_comps,
        use_highly_variable=True,
        svd_solver="arpack",
        random_state=UMAP_SEED,
    )
    n_neighbors = min(15, max(2, adata.n_obs - 1))
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_comps,
                    method="umap", metric="euclidean", random_state=UMAP_SEED)
    sc.tl.umap(adata, min_dist=0.5, spread=1.0, n_components=2,
               alpha=1.0, gamma=1.0, negative_sample_rate=5,
               init_pos="spectral", random_state=UMAP_SEED)
    method = {
        "source_modality": "rna",
        "source_matrix": "rna.X; no circRNA features enter preprocessing or embedding",
        "embedding_decision": "Regenerate using c99cdda historical exporter method; source RNA has no stored embedding",
        "source_embedding_keys": list(rna.obsm.keys()),
        "scanpy_version": importlib.metadata.version("scanpy"),
        "shared_cells": len(shared_cells),
        "normalize_total_target_sum": 10_000,
        "transform": "scanpy.pp.log1p",
        "highly_variable_genes_flavor": "seurat",
        "n_top_genes": n_top_genes,
        "selected_variable_genes": int(adata.var["highly_variable"].sum()),
        "pca_components": n_comps,
        "pca_solver": "arpack",
        "neighbors": n_neighbors,
        "neighbors_method": "umap",
        "neighbors_metric": "euclidean",
        "neighbors_representation": "X_pca",
        "umap": {"min_dist": 0.5, "spread": 1.0, "n_components": 2,
                 "alpha": 1.0, "gamma": 1.0, "negative_sample_rate": 5,
                 "init_pos": "spectral", "maxiter": None,
                 "effective_epochs": 500 if adata.n_obs <= 10000 else 200,
                 "a": float(adata.uns["umap"]["params"]["a"]),
                 "b": float(adata.uns["umap"]["params"]["b"])},
        "threads": 1,
        "random_seed": UMAP_SEED,
    }
    return np.asarray(adata.obsm["X_umap"]), method


def write_json(path: Path, payload: object) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_figures(cells: pd.DataFrame, outdir: Path) -> list[Path]:
    import matplotlib

    matplotlib.use("Agg")
    matplotlib.rcParams.update({"svg.hashsalt": "circyto-application-note-v2",
                                "font.family": "DejaVu Sans", "font.size": 10})
    import matplotlib.pyplot as plt

    outputs: list[Path] = []
    plot_specs = [
        (
            "figure1b_smartseq3_circrna_burden",
            "circRNA_count",
            "Detected circRNA candidates per cell",
            "viridis",
        ),
    ]
    for stem, column, colorbar_label, cmap in plot_specs:
        fig, ax = plt.subplots(figsize=(5.2, 4.3), constrained_layout=True)
        points = ax.scatter(
            cells["UMAP1"],
            cells["UMAP2"],
            c=cells[column],
            cmap=cmap,
            s=28,
            linewidths=0.25,
            edgecolors="white",
        )
        colorbar = fig.colorbar(points, ax=ax, fraction=0.046, pad=0.04)
        colorbar.set_label(colorbar_label)
        ax.set_xlabel("RNA UMAP 1")
        ax.set_ylabel("RNA UMAP 2")
        ax.set_title("Smart-seq3 circRNA burden")
        ax.set_xticks([])
        ax.set_yticks([])
        for suffix in ("png", "svg"):
            path = outdir / f"{stem}.{suffix}"
            save_kwargs = {"dpi": 300} if suffix == "png" else {"metadata": {"Date": None}}
            fig.savefig(path, **save_kwargs)
            outputs.append(path)
        plt.close(fig)

    detected = cells["man1a2_candidate_detected"].astype(bool)
    fig, ax = plt.subplots(figsize=(5.2, 4.3), constrained_layout=True)
    ax.scatter(
        cells.loc[~detected, "UMAP1"],
        cells.loc[~detected, "UMAP2"],
        c="#d3d3d3",
        s=25,
        linewidths=0,
        label=f"Not detected (n={int((~detected).sum())})",
    )
    ax.scatter(
        cells.loc[detected, "UMAP1"],
        cells.loc[detected, "UMAP2"],
        c="#c44e52",
        s=40,
        linewidths=0.4,
        edgecolors="white",
        label=f"Detected (n={int(detected.sum())})",
    )
    ax.set_xlabel("RNA UMAP 1")
    ax.set_ylabel("RNA UMAP 2")
    ax.set_title("MAN1A2-associated circRNA candidate")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.legend(frameon=False, loc="best")
    for suffix in ("png", "svg"):
        path = outdir / f"figure1c_smartseq3_man1a2_detection.{suffix}"
        save_kwargs = {"dpi": 300} if suffix == "png" else {"metadata": {"Date": None}}
        fig.savefig(path, **save_kwargs)
        outputs.append(path)
    plt.close(fig)
    return outputs


def package_versions() -> dict[str, str]:
    import circyto

    names = ["anndata", "mudata", "numpy", "pandas", "scipy", "scanpy", "matplotlib", "umap-learn", "h5py", "numba", "llvmlite", "scikit-learn", "pynndescent"]
    versions: dict[str, str] = {"circyto": f"{circyto.__version__} (source tree)"}
    for name in names:
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            versions[name] = "not-installed"
    return versions


def agreement_rows(smart: dict[str, object], imr90: dict[str, object]) -> pd.DataFrame:
    rows = [
        ("Smart-seq3", "cells", "192", str(smart["cell_count"]), smart["cell_count"] == 192),
        ("Smart-seq3", "RNA features", "63187", str(smart["rna_feature_count"]), smart["rna_feature_count"] == 63187),
        ("Smart-seq3", "circRNA candidates", "2503", str(smart["circRNA_candidate_count"]), smart["circRNA_candidate_count"] == 2503),
        ("Smart-seq3", "nonzero circRNA entries", "2659", str(smart["circRNA_nonzero_entries"]), smart["circRNA_nonzero_entries"] == 2659),
        ("Smart-seq3", "host-gene annotated", "2379", str(smart["host_gene_annotated_count"]), smart["host_gene_annotated_count"] == 2379),
        ("Smart-seq3", "host-gene fraction", "0.950459448661606", repr(smart["host_gene_annotation_fraction"]), np.isclose(smart["host_gene_annotation_fraction"], 2379 / 2503)),
        ("Smart-seq3", "host-gene percent (1 dp)", "95.0", f"{100 * smart['host_gene_annotation_fraction']:.1f}", round(100 * smart["host_gene_annotation_fraction"], 1) == 95.0),
        ("Smart-seq3", "median detected circRNAs/cell", "12", repr(smart["median_detected_circRNAs_per_cell"]), smart["median_detected_circRNAs_per_cell"] == 12.0),
        ("Smart-seq3", "median total support/cell", "22.5", repr(smart["median_total_circRNA_support_per_cell"]), smart["median_total_circRNA_support_per_cell"] == 22.5),
        ("Smart-seq3", "MAN1A2 candidate detecting cells", "6", str(smart["man1a2_candidate_detecting_cells"]), smart["man1a2_candidate_detecting_cells"] == 6),
        ("Smart-seq3", "MAN1A2 candidate total support", "15", repr(smart["man1a2_candidate_total_support"]), smart["man1a2_candidate_total_support"] == 15.0),
        ("IMR90", "RNA shape", "23x63187", "x".join(map(str, imr90["rna_shape"])), imr90["rna_shape"] == [23, 63187]),
        ("IMR90", "circRNA shape", "23x2443", "x".join(map(str, imr90["circRNA_shape"])), imr90["circRNA_shape"] == [23, 2443]),
        ("IMR90", "CNV shape", "23x60607", "x".join(map(str, imr90["cnv_shape"])), imr90["cnv_shape"] == [23, 60607]),
        ("IMR90", "trimodal overlap", "23", str(imr90["trimodal_cell_overlap"]), imr90["trimodal_cell_overlap"] == 23),
        ("IMR90", "host-gene annotated", "2429", str(imr90["host_gene_annotated_count"]), imr90["host_gene_annotated_count"] == 2429),
        ("IMR90", "host-gene fraction", "0.994269340974212", repr(imr90["host_gene_annotation_fraction"]), np.isclose(imr90["host_gene_annotation_fraction"], 2429 / 2443)),
        ("IMR90", "host-gene percent (1 dp)", "99.4", f"{100 * imr90['host_gene_annotation_fraction']:.1f}", round(100 * imr90["host_gene_annotation_fraction"], 1) == 99.4),
    ]
    frame = pd.DataFrame(rows, columns=["dataset", "quantity", "historical_value", "regenerated_value", "match"])
    frame["source_checksum"] = frame["dataset"].map({"Smart-seq3": SMARTSEQ3_SHA256, "IMR90": IMR90_SHA256})
    frame["evidence_class"] = "checksum-verified processed-object regeneration"
    return frame


def validate_output_directory(path: Path, inputs: list[Path]) -> Path:
    resolved = path.resolve()
    manuscript_root = (REPO_ROOT / "manuscript").resolve()
    if manuscript_root not in resolved.parents:
        raise SystemExit("[ERROR] Output must be a new directory under this checkout's manuscript/.")
    if resolved.exists():
        raise SystemExit("[ERROR] Output directory already exists; choose a new manuscript directory.")
    if any(resolved == source.resolve() or resolved in source.resolve().parents for source in inputs):
        raise SystemExit("[ERROR] Output directory overlaps an input.")
    return resolved


def repository_provenance() -> dict[str, object]:
    def git(*args: str) -> str:
        return subprocess.check_output(["git", "-C", str(REPO_ROOT), *args], text=True).strip()

    return {"commit": git("rev-parse", "HEAD"),
            "worktree_status": git("status", "--short"),
            "frozen_baseline": "44697355bcab1c525ca7ef9b130e2ad0094d9e1b",
            "production_diff_from_baseline": git("diff", "44697355bcab1c525ca7ef9b130e2ad0094d9e1b", "--", "circyto", "pyproject.toml")}


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Regenerate checksum-gated Application Note numbers and Figure 1B/1C.")
    parser.add_argument("--smartseq3", required=True, type=Path, help="Checksum-matched Smart-seq3 H5MU object.")
    parser.add_argument("--imr90", required=True, type=Path, help="Checksum-matched IMR90 H5MU object.")
    parser.add_argument("--output", "--outdir", dest="outdir", required=True, type=Path, help="New output directory under this checkout's manuscript/.")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)

    # Check both complete byte streams before creating the output directory.
    checksums = [
        verify_input(args.smartseq3, SMARTSEQ3_SHA256, "Smart-seq3"),
        verify_input(args.imr90, IMR90_SHA256, "IMR90"),
    ]
    args.outdir = validate_output_directory(args.outdir, [args.smartseq3, args.imr90])
    repository = repository_provenance()
    if repository["production_diff_from_baseline"]:
        raise SystemExit("[ERROR] Production source differs from the frozen manuscript baseline.")
    args.outdir.mkdir(parents=True, exist_ok=False)
    os.environ["MPLCONFIGDIR"] = str(args.outdir / ".cache" / "matplotlib")
    os.environ["NUMBA_CACHE_DIR"] = str(args.outdir / ".cache" / "numba")

    smart_mdata = read_modalities(args.smartseq3)
    imr90_mdata = read_modalities(args.imr90)
    smart_rna = smart_mdata["rna"]
    smart_circ = smart_mdata["circ"]
    imr90_rna = imr90_mdata["rna"]
    imr90_circ = imr90_mdata["circ"]
    imr90_cnv = imr90_mdata["cnv"]

    shared_smart = shared_obs_in_left_order(smart_rna, smart_circ)
    if len(shared_smart) != smart_rna.n_obs or len(shared_smart) != smart_circ.n_obs:
        raise SystemExit("[ERROR] Smart-seq3 RNA/circ axes are not identical; no outputs generated.")
    smart_circ = smart_circ[shared_smart, :].copy()

    circ_count = nonzero_per_row(smart_circ.X)
    circ_support = sum_per_row(smart_circ.X)
    candidate_matches = np.flatnonzero(np.asarray(smart_circ.var_names == MAN1A2_CIRC_ID))
    if len(candidate_matches) != 1:
        raise SystemExit(f"[ERROR] Expected exactly one {MAN1A2_CIRC_ID} feature; found {len(candidate_matches)}.")
    candidate_index = int(candidate_matches[0])
    candidate_support = column_values(smart_circ.X, candidate_index)
    candidate_host_gene = str(smart_circ.var.iloc[candidate_index].get("host_gene", ""))
    if "MAN1A2" not in {part.strip() for part in candidate_host_gene.split(";")}:
        raise SystemExit(
            f"[ERROR] Candidate {MAN1A2_CIRC_ID} is not annotated to MAN1A2; observed host_gene={candidate_host_gene!r}."
        )

    coordinates, umap_method = compute_rna_umap(smart_rna, shared_smart)
    cells = pd.DataFrame(
        {
            "cell_id": shared_smart,
            "UMAP1": coordinates[:, 0],
            "UMAP2": coordinates[:, 1],
            "circRNA_count": circ_count,
            "circRNA_total_support": circ_support,
            "man1a2_candidate_id": MAN1A2_CIRC_ID,
            "man1a2_candidate_detected": candidate_support > 0,
            "man1a2_candidate_support": candidate_support,
        }
    )

    smart_annotated = nonempty_host_gene_count(smart_circ)
    smart_summary: dict[str, object] = {
        "dataset": "Smart-seq3 / E-MTAB-8735",
        "cell_count": int(smart_rna.n_obs),
        "rna_feature_count": int(smart_rna.n_vars),
        "circRNA_candidate_count": int(smart_circ.n_vars),
        "circRNA_nonzero_entries": int(circ_count.sum()),
        "host_gene_annotated_count": smart_annotated,
        "host_gene_annotation_fraction": smart_annotated / smart_circ.n_vars,
        "median_detected_circRNAs_per_cell": float(np.median(circ_count)),
        "median_total_circRNA_support_per_cell": float(np.median(circ_support)),
        "man1a2_candidate_id": MAN1A2_CIRC_ID,
        "man1a2_candidate_host_gene": candidate_host_gene,
        "man1a2_candidate_detecting_cells": int(np.count_nonzero(candidate_support > 0)),
        "man1a2_candidate_total_support": float(candidate_support.sum()),
        "umap_method": umap_method,
    }

    imr90_annotated = nonempty_host_gene_count(imr90_circ)
    imr90_summary: dict[str, object] = {
        "dataset": "IMR90 scRR / GSE278958",
        "rna_shape": [int(imr90_rna.n_obs), int(imr90_rna.n_vars)],
        "circRNA_shape": [int(imr90_circ.n_obs), int(imr90_circ.n_vars)],
        "cnv_shape": [int(imr90_cnv.n_obs), int(imr90_cnv.n_vars)],
        "trimodal_cell_overlap": len(shared_obs_in_left_order(imr90_rna, imr90_circ, imr90_cnv)),
        "host_gene_annotated_count": imr90_annotated,
        "host_gene_annotation_fraction": imr90_annotated / imr90_circ.n_vars,
    }

    comparison = agreement_rows(smart_summary, imr90_summary)
    all_agree = bool(comparison["match"].all())
    if not all_agree:
        disagreements = comparison.loc[~comparison["match"]]
        print(f"[WARNING] Regenerated values supersede historical expectations; investigate:\n{disagreements.to_string(index=False)}", file=sys.stderr)

    # Check again before publishing summaries; the script never writes an H5MU.
    verify_input(args.smartseq3, SMARTSEQ3_SHA256, "Smart-seq3")
    verify_input(args.imr90, IMR90_SHA256, "IMR90")
    write_json(args.outdir / "regeneration_summary.json", {"smartseq3": smart_summary, "imr90": imr90_summary, "all_evidence_values_agree": all_agree})
    for label, summary in (("smartseq3", smart_summary), ("imr90", imr90_summary)):
        pd.DataFrame([{"quantity": key, "value": json.dumps(value) if isinstance(value, (list, dict)) else value}
                      for key, value in summary.items() if key != "umap_method"]).to_csv(
                          args.outdir / f"{label}_summary.tsv", sep="\t", index=False)
    cells[["cell_id", "UMAP1", "UMAP2", "circRNA_count", "circRNA_total_support"]].to_csv(args.outdir / "smartseq3_umap_cells.tsv", sep="\t", index=False)
    cells[["cell_id", "UMAP1", "UMAP2", "man1a2_candidate_id", "man1a2_candidate_detected", "man1a2_candidate_support"]].to_csv(args.outdir / "smartseq3_selected_candidate.tsv", sep="\t", index=False)
    comparison.to_csv(args.outdir / "evidence_comparison.tsv", sep="\t", index=False)
    figure_paths = write_figures(cells, args.outdir)

    generated = sorted(path for path in args.outdir.iterdir() if path.is_file())
    provenance = {
        "schema": "circyto.application_note_results.v2",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "read_only_inputs": True,
        "reader": "h5py.File(mode='r') + anndata.io.read_elem on mod/*; no global MuData synchronization",
        "input_checksums": checksums,
        "candidate_id": MAN1A2_CIRC_ID,
        "umap_method": umap_method,
        "python": sys.version,
        "python_executable": Path(sys.executable).name,
        "platform": platform.platform(),
        "package_versions": package_versions(),
        "script": str(Path(__file__).resolve().relative_to(REPO_ROOT)),
        "repository": repository,
        "command": ["python", "scripts/manuscript/regenerate_application_note_results.py", "--smartseq3", str(args.smartseq3), "--imr90", str(args.imr90), "--output", str(args.outdir.relative_to(REPO_ROOT))],
        "input_path_convention": "Exact supplied paths; relative input paths resolve from the invocation directory. Materialization root is not needed for content identity.",
        "installed_distributions": sorted({f"{dist.metadata['Name']}=={dist.version}" for dist in importlib.metadata.distributions()}),
        "script_sha256": sha256_file(Path(__file__).resolve()),
        "figures": [path.name for path in figure_paths],
        "generated_files": {path.name: sha256_file(path) for path in generated},
        "all_evidence_values_agree": all_agree,
    }
    write_json(args.outdir / "provenance.json", provenance)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
