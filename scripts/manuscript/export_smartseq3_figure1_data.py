#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from _mudata_utils import (
    CIRC_MODALITY_CANDIDATES,
    RNA_MODALITY_CANDIDATES,
    common_obs_names,
    dense_matrix,
    find_modality,
    get_or_compute_circ_metric,
    get_or_compute_rna_metric,
    load_mudata,
    nonempty_string_mask,
    split_host_genes,
    write_table,
)
from manuscript_object_specs import MANUSCRIPT_OBJECTS


DEFAULT_SELECTED_HOSTS = ["DYNC2I1", "MAN1A2", "PTP4A2"]


def normalize_label(value: object) -> str:
    text = str(value)
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", text)
    return text.strip("_") or "feature"


def host_gene_membership(circ) -> dict[str, list[int]]:
    membership: dict[str, list[int]] = {}
    if "host_gene" not in circ.var.columns:
        return membership
    for idx, value in enumerate(circ.var["host_gene"]):
        for gene in split_host_genes(value):
            membership.setdefault(gene, []).append(idx)
    return membership


def circ_candidate_table(circ, *, top_n: int) -> pd.DataFrame:
    matrix = dense_matrix(circ).astype(float)
    detected = matrix > 0
    var = circ.var.copy()
    host_gene = var["host_gene"].astype(str) if "host_gene" in var.columns else pd.Series("", index=var.index)
    host_source = (
        var["host_gene_source"].astype(str) if "host_gene_source" in var.columns else pd.Series("", index=var.index)
    )
    rows = []
    for idx, circ_id in enumerate(map(str, circ.var_names)):
        total_support = float(matrix[:, idx].sum())
        detected_cells = int(detected[:, idx].sum())
        rows.append(
            {
                "rank_metric": "prevalence_then_total_support",
                "circ_id": circ_id,
                "host_gene": host_gene.iloc[idx],
                "host_gene_source": host_source.iloc[idx],
                "host_gene_annotated": bool(nonempty_string_mask(pd.Series([host_gene.iloc[idx]])).iloc[0]),
                "detected_cells": detected_cells,
                "prevalence": detected_cells / circ.n_obs if circ.n_obs else 0.0,
                "total_support": total_support,
                "mean_support": float(matrix[:, idx].mean()) if circ.n_obs else 0.0,
                "max_cell_support": float(matrix[:, idx].max()) if circ.n_obs else 0.0,
                "sparsity": 1.0 - (detected_cells / circ.n_obs if circ.n_obs else 0.0),
            }
        )
    return (
        pd.DataFrame(rows)
        .sort_values(["prevalence", "total_support", "host_gene_annotated", "circ_id"], ascending=[False, False, False, True])
        .head(top_n)
        .reset_index(drop=True)
        .assign(rank=lambda frame: np.arange(1, len(frame) + 1))
    )


def host_candidate_table(circ, *, top_n: int) -> pd.DataFrame:
    matrix = dense_matrix(circ).astype(float)
    membership = host_gene_membership(circ)
    rows = []
    for host_gene, cols in sorted(membership.items()):
        values = matrix[:, cols].sum(axis=1)
        detected_cells = int(np.count_nonzero(values > 0))
        rows.append(
            {
                "rank_metric": "prevalence_then_total_support",
                "host_gene": host_gene,
                "n_circRNAs": int(len(cols)),
                "detected_cells": detected_cells,
                "prevalence": detected_cells / circ.n_obs if circ.n_obs else 0.0,
                "total_support": float(values.sum()),
                "mean_support": float(values.mean()) if circ.n_obs else 0.0,
                "max_cell_support": float(values.max()) if circ.n_obs else 0.0,
                "sparsity": 1.0 - (detected_cells / circ.n_obs if circ.n_obs else 0.0),
                "circ_ids": ";".join(map(str, circ.var_names[cols])),
            }
        )
    return (
        pd.DataFrame(rows)
        .sort_values(["prevalence", "total_support", "n_circRNAs", "host_gene"], ascending=[False, False, False, True])
        .head(top_n)
        .reset_index(drop=True)
        .assign(rank=lambda frame: np.arange(1, len(frame) + 1))
    )


def select_feature_values(circ, hosts: list[str]) -> tuple[pd.DataFrame, dict[str, np.ndarray]]:
    matrix = dense_matrix(circ).astype(float)
    membership = host_gene_membership(circ)
    cell_ids = list(map(str, circ.obs_names))
    selected_rows: list[dict[str, object]] = []
    cell_columns: dict[str, np.ndarray] = {}
    var = circ.var.copy()
    host_source = (
        var["host_gene_source"].astype(str) if "host_gene_source" in var.columns else pd.Series("", index=var.index)
    )

    def append_values(
        *,
        label: str,
        role: str,
        scope: str,
        feature_id: str,
        host_gene: str,
        source: str,
        values: np.ndarray,
    ) -> None:
        detected_cells = int(np.count_nonzero(values > 0))
        prevalence = detected_cells / len(values) if len(values) else 0.0
        total_support = float(values.sum())
        column_name = normalize_label(f"{scope}_{label}_support")
        cell_columns[column_name] = values
        for cell_id, support in zip(cell_ids, values):
            selected_rows.append(
                {
                    "candidate_label": label,
                    "selection_role": role,
                    "feature_scope": scope,
                    "feature_id": feature_id,
                    "host_gene": host_gene,
                    "host_gene_source": source,
                    "cell_id": cell_id,
                    "support": float(support),
                    "detected_cells": detected_cells,
                    "prevalence": prevalence,
                    "total_support": total_support,
                    "sparsity": 1.0 - prevalence,
                }
            )

    for host in hosts:
        cols = membership.get(host, [])
        if not cols:
            continue
        aggregate = matrix[:, cols].sum(axis=1)
        role = "primary_host_gene_candidate" if host == "MAN1A2" else "backup_host_gene_candidate"
        append_values(
            label=f"{host}_host",
            role=role,
            scope="host_gene",
            feature_id=host,
            host_gene=host,
            source="aggregate",
            values=aggregate,
        )
        circ_rank = sorted(
            cols,
            key=lambda idx: (
                int(np.count_nonzero(matrix[:, idx] > 0)),
                float(matrix[:, idx].sum()),
                str(circ.var_names[idx]),
            ),
            reverse=True,
        )
        top_idx = circ_rank[0]
        role = "primary_individual_circRNA_candidate" if host == "DYNC2I1" else "backup_individual_circRNA_candidate"
        append_values(
            label=f"{host}_top_circ",
            role=role,
            scope="circRNA",
            feature_id=str(circ.var_names[top_idx]),
            host_gene=host,
            source=str(host_source.iloc[top_idx]),
            values=matrix[:, top_idx],
        )
    return pd.DataFrame(selected_rows), cell_columns


def compute_umap_frame(mdata, rna, circ, *, seed: int, n_top_genes: int) -> pd.DataFrame:
    try:
        import scanpy as sc  # type: ignore
    except Exception as exc:
        raise SystemExit("scanpy is required to export RNA-derived UMAP coordinates; install circyto[scanpy].") from exc

    shared = common_obs_names(rna, circ)
    if not shared:
        raise SystemExit("No shared Smart-seq3 cells across RNA and circRNA modalities.")
    adata = rna[shared, :].copy()
    adata.var_names_make_unique()
    sc.pp.normalize_total(adata, target_sum=10_000)
    sc.pp.log1p(adata)
    hvg_count = min(n_top_genes, adata.n_vars)
    sc.pp.highly_variable_genes(adata, n_top_genes=hvg_count, flavor="seurat")
    n_comps = min(50, max(2, adata.n_obs - 1), max(2, adata.n_vars - 1))
    sc.tl.pca(adata, n_comps=n_comps, use_highly_variable=True, svd_solver="arpack", random_state=seed)
    n_neighbors = min(15, max(2, adata.n_obs - 1))
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_comps, random_state=seed)
    sc.tl.umap(adata, random_state=seed)
    coordinates = pd.DataFrame(
        {
            "cell_id": list(map(str, adata.obs_names)),
            "rna_umap1": adata.obsm["X_umap"][:, 0],
            "rna_umap2": adata.obsm["X_umap"][:, 1],
            "detected_genes": get_or_compute_rna_metric(mdata, rna, "detected_genes").reindex(shared).values,
            "circRNA_count": get_or_compute_circ_metric(mdata, circ, "circRNA_count").reindex(shared).values,
            "circRNA_total_support": get_or_compute_circ_metric(mdata, circ, "circRNA_total_support").reindex(shared).values,
        }
    )
    coordinates.attrs["preprocessing"] = {
        "normalization": "scanpy.pp.normalize_total target_sum=10000",
        "transform": "scanpy.pp.log1p",
        "highly_variable_genes": int(hvg_count),
        "pca_components": int(n_comps),
        "neighbors": int(n_neighbors),
        "umap_random_state": int(seed),
    }
    return coordinates


def method_text(input_path: Path, seed: int, n_top_genes: int) -> str:
    return f"""# Smart-seq3 Figure 1 Data Export

Input object: `{input_path}`

This script exports deterministic, notebook-ready data only. It does not render the final manuscript figure.

RNA preprocessing for `smartseq3_umap_cells.tsv`:

- subset to cells shared by RNA and circRNA modalities;
- `scanpy.pp.normalize_total(target_sum=10000)`;
- `scanpy.pp.log1p`;
- highly variable genes with `flavor="seurat"` and `n_top_genes={n_top_genes}`;
- PCA with at most 50 components and fixed `random_state={seed}`;
- neighbors with at most 15 neighbors;
- UMAP with fixed `random_state={seed}`.

Candidate tables are ranked by prevalence first, then total support. The selected-candidate export reports prevalence, total support, annotation source, sparsity, and per-cell support values so final notebook choices are not forced by total support alone.
"""


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Export deterministic Smart-seq3 Figure 1 data tables.")
    parser.add_argument("--input", type=Path, default=MANUSCRIPT_OBJECTS["Smart-seq3"]["default_path"])
    parser.add_argument("--outdir", type=Path, default=Path("manuscript/results"))
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument("--n-top-genes", type=int, default=2000)
    parser.add_argument("--top-n-circs", type=int, default=100)
    parser.add_argument("--top-n-hosts", type=int, default=100)
    parser.add_argument("--selected-host", action="append", default=None, help="Host gene to include in selected values.")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    mdata = load_mudata(args.input)
    rna = mdata.mod[find_modality(mdata, RNA_MODALITY_CANDIDATES)]
    circ = mdata.mod[find_modality(mdata, CIRC_MODALITY_CANDIDATES)]

    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    selected_hosts = args.selected_host or DEFAULT_SELECTED_HOSTS

    selected_values, selected_cell_columns = select_feature_values(circ, selected_hosts)
    umap = compute_umap_frame(mdata, rna, circ, seed=args.seed, n_top_genes=args.n_top_genes)
    for column_name, values in selected_cell_columns.items():
        umap[column_name] = pd.Series(values, index=list(map(str, circ.obs_names))).reindex(umap["cell_id"]).values

    write_table(umap, outdir / "smartseq3_umap_cells.tsv")
    write_table(circ_candidate_table(circ, top_n=args.top_n_circs), outdir / "smartseq3_top_circRNA_candidates.tsv")
    write_table(host_candidate_table(circ, top_n=args.top_n_hosts), outdir / "smartseq3_top_hostgene_features.tsv")
    write_table(selected_values, outdir / "smartseq3_selected_feature_candidates.tsv")
    (outdir / "smartseq3_figure1_data_method.md").write_text(
        method_text(args.input, args.seed, args.n_top_genes),
        encoding="utf-8",
    )
    (outdir / "smartseq3_figure1_data_summary.json").write_text(
        json.dumps(
            {
                "input": str(args.input),
                "seed": args.seed,
                "n_top_genes": args.n_top_genes,
                "selected_hosts": selected_hosts,
                "outputs": [
                    "smartseq3_umap_cells.tsv",
                    "smartseq3_top_circRNA_candidates.tsv",
                    "smartseq3_top_hostgene_features.tsv",
                    "smartseq3_selected_feature_candidates.tsv",
                    "smartseq3_figure1_data_method.md",
                ],
            },
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
