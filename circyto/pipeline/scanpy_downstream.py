from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from circyto.multimodal.sync import read_h5mu

from circyto.pipeline.workflow_reporting import numeric_summary, write_json

try:
    import anndata as ad  # type: ignore

    HAS_ANNDATA = True
except Exception:
    HAS_ANNDATA = False

try:
    import mudata as mu  # type: ignore

    HAS_MUDATA = True
except Exception:
    HAS_MUDATA = False

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import scanpy as sc  # type: ignore

    HAS_SCANPY = True
except Exception:
    HAS_SCANPY = False


def _require_scanpy_stack() -> None:
    if not HAS_ANNDATA:
        raise RuntimeError("anndata is required for Scanpy downstream analysis")
    if not HAS_MUDATA:
        raise RuntimeError("mudata is required for Scanpy downstream analysis; install circyto[scanpy]")
    if not HAS_SCANPY:
        raise RuntimeError("scanpy is not installed; install circyto[scanpy] or pip install scanpy")


def _load_rna_from_h5mu(input_path: Path):
    _require_scanpy_stack()
    if not input_path.exists():
        raise FileNotFoundError(f"MuData input not found: {input_path}")
    mdata = read_h5mu(str(input_path))
    if "rna" not in mdata.mod:
        raise ValueError(f"{input_path} does not contain an 'rna' modality")
    adata = mdata.mod["rna"].copy()
    if "gene_name" in adata.var.columns:
        gene_names = adata.var["gene_name"].astype(str)
    else:
        gene_names = pd.Series(adata.var_names.astype(str), index=adata.var_names)
    adata.var["mt"] = gene_names.str.upper().str.startswith("MT-")
    adata.var_names_make_unique()
    return adata


def _save_current_figure(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.gcf().savefig(path, dpi=150, bbox_inches="tight")
    plt.close(plt.gcf())


def scanpy_qc_report(*, input_path: Path, output_dir: Path) -> dict[str, Any]:
    adata = _load_rna_from_h5mu(input_path)
    output_dir.mkdir(parents=True, exist_ok=True)

    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, log1p=False, percent_top=None)

    sc.pl.violin(adata, ["total_counts", "n_genes_by_counts", "pct_counts_mt"], show=False)
    violin_path = output_dir / "qc_violin.png"
    _save_current_figure(violin_path)

    sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", show=False)
    counts_vs_genes_path = output_dir / "qc_scatter_counts_vs_genes.png"
    _save_current_figure(counts_vs_genes_path)

    sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt", show=False)
    mt_vs_counts_path = output_dir / "qc_mt_vs_counts.png"
    _save_current_figure(mt_vs_counts_path)

    summary = {
        "input": str(input_path.resolve()),
        "output_dir": str(output_dir.resolve()),
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "total_counts": numeric_summary(pd.to_numeric(adata.obs["total_counts"], errors="coerce").dropna().tolist()),
        "n_genes_by_counts": numeric_summary(pd.to_numeric(adata.obs["n_genes_by_counts"], errors="coerce").dropna().tolist()),
        "pct_counts_mt": numeric_summary(pd.to_numeric(adata.obs["pct_counts_mt"], errors="coerce").dropna().tolist()),
        "outputs": {
            "qc_violin_png": str(violin_path.resolve()),
            "qc_scatter_counts_vs_genes_png": str(counts_vs_genes_path.resolve()),
            "qc_mt_vs_counts_png": str(mt_vs_counts_path.resolve()),
        },
        "exploratory_only": True,
    }
    write_json(output_dir / "scanpy_qc_summary.json", summary)
    return summary


def scanpy_pca_workflow(*, input_path: Path, output_dir: Path) -> dict[str, Any]:
    adata = _load_rna_from_h5mu(input_path)
    output_dir.mkdir(parents=True, exist_ok=True)

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    n_top_genes = min(2000, adata.n_vars)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor="seurat")

    pca_comps = min(50, max(2, min(adata.n_obs - 1, adata.n_vars - 1)))
    sc.tl.pca(adata, n_comps=pca_comps, use_highly_variable=True if "highly_variable" in adata.var.columns else None)

    n_neighbors = min(15, max(2, adata.n_obs - 1))
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=min(pca_comps, max(2, pca_comps)))
    sc.tl.umap(adata)
    try:
        sc.tl.leiden(adata)
    except Exception as exc:
        raise RuntimeError(
            "Leiden clustering failed. Ensure scanpy clustering dependencies are installed."
        ) from exc

    sc.pl.umap(adata, color="leiden", show=False)
    umap_path = output_dir / "rna_umap.png"
    _save_current_figure(umap_path)

    leiden_path = output_dir / "rna_leiden.tsv"
    adata.obs[["leiden"]].reset_index(names="cell_id").to_csv(leiden_path, sep="\t", index=False)

    processed_path = output_dir / "exploratory_rna_processed.h5ad"
    adata.write_h5ad(processed_path)

    summary = {
        "input": str(input_path.resolve()),
        "output_dir": str(output_dir.resolve()),
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "n_highly_variable_genes": int(adata.var["highly_variable"].sum()) if "highly_variable" in adata.var.columns else 0,
        "pca_components": int(pca_comps),
        "n_neighbors": int(n_neighbors),
        "outputs": {
            "rna_umap_png": str(umap_path.resolve()),
            "rna_leiden_tsv": str(leiden_path.resolve()),
            "exploratory_rna_processed_h5ad": str(processed_path.resolve()),
        },
        "exploratory_only": True,
    }
    write_json(output_dir / "scanpy_pca_summary.json", summary)
    return summary
