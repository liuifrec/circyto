#!/usr/bin/env python
from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy import sparse

from _mudata_utils import (
    CIRC_MODALITY_CANDIDATES,
    CNV_MODALITY_CANDIDATES,
    RNA_MODALITY_CANDIDATES,
    RT_MODALITY_CANDIDATES,
    common_obs_names,
    dense_matrix,
    find_modality,
    get_or_compute_circ_metric,
    host_gene_recovery,
    host_gene_source_counts,
    json_dumps,
    load_mudata,
    nonempty_string_mask,
    split_host_genes,
    write_table,
)


@dataclass(frozen=True)
class DatasetSpec:
    label: str
    source_accession: str
    protocol: str
    path: Path
    expected_aux_modality: str | None = None


DEFAULT_SPECS = (
    DatasetSpec(
        label="Smart-seq3",
        source_accession="E-MTAB-8735",
        protocol="Smart-seq3",
        path=Path("/mnt/d/circyto/load_work/emtab8735_smartseq3/full_length.hostgene_fixed.h5mu"),
    ),
    DatasetSpec(
        label="HAP1",
        source_accession="GSE278952",
        protocol="scRR-seq RNA+circ+RT",
        path=Path("/mnt/d/circyto/load_work/scrr_hap1/full_length_rna_circ_rt.hostgene_fixed.h5mu"),
        expected_aux_modality="rt",
    ),
    DatasetSpec(
        label="IMR90",
        source_accession="GSE278958",
        protocol="scRR-seq RNA+circ+CNV",
        path=Path("/mnt/d/circyto/load_work/scrr_imr90/full_length_rna_circ_cnv.hostgene_fixed.h5mu"),
        expected_aux_modality="cnv",
    ),
)


EXPECTED_MANUSCRIPT_NUMBERS = {
    "Smart-seq3": {
        "cells": 192,
        "rna_features": 63187,
        "circRNAs": 2503,
        "host_gene_annotated": 2379,
        "host_gene_total": 2503,
        "host_annotation_rate_pct_rounded": 95.0,
    },
    "HAP1": {
        "cells": 63,
        "circRNAs": 3209,
        "rt_overlap": 56,
        "host_gene_annotated": 3117,
        "host_gene_total": 3209,
        "host_annotation_rate_pct_rounded": 97.1,
    },
    "IMR90": {
        "cells": 23,
        "circRNAs": 2443,
        "cnv_overlap": 23,
        "host_gene_annotated": 2429,
        "host_gene_total": 2443,
        "host_annotation_rate_pct_rounded": 99.4,
    },
}


HOST_GENE_COLUMNS = [
    "host_gene",
    "host_gene_source",
    "host_gene_from_gtf",
    "host_gene_from_circatlas",
    "host_gene_from_circatlas_id",
    "host_gene_id",
    "host_genes_multi",
    "host_gene_ids_multi",
    "host_gene_n",
]


def _json_list(values: Iterable[object]) -> str:
    return json_dumps([str(value) for value in values])


def _matrix_as_array(adata) -> np.ndarray:
    return np.asarray(dense_matrix(adata), dtype=float)


def _rank_pct(values: Iterable[float]) -> np.ndarray:
    series = pd.Series(list(values), dtype=float)
    if len(series) == 0:
        return np.asarray([], dtype=float)
    if series.nunique(dropna=True) <= 1:
        return np.full(len(series), 0.5, dtype=float)
    return series.rank(method="average", pct=True).to_numpy(dtype=float)


def _pairwise_overlap(mdata, modalities: list[str]) -> dict[str, int]:
    overlap: dict[str, int] = {}
    for i, left in enumerate(modalities):
        for right in modalities[i + 1 :]:
            overlap[f"{left}&{right}"] = len(common_obs_names(mdata.mod[left], mdata.mod[right]))
    return overlap


def _all_modality_overlap(mdata, modalities: list[str]) -> int:
    if not modalities:
        return 0
    return len(common_obs_names(*(mdata.mod[name] for name in modalities)))


def _all_modality_union(mdata, modalities: list[str]) -> int:
    if not modalities:
        return 0
    return len(set().union(*(set(map(str, mdata.mod[name].obs_names)) for name in modalities)))


def _source_counts_for_columns(circ, cols: list[int]) -> dict[str, int]:
    if "host_gene_source" not in circ.var.columns or not cols:
        return {}
    values = circ.var.iloc[cols]["host_gene_source"].fillna("").astype(str).str.strip()
    values = values.where(values != "", other="missing")
    return {str(k): int(v) for k, v in values.value_counts().sort_index().items()}


def _dominant_source(counts: dict[str, int]) -> tuple[str, int]:
    if not counts:
        return "", 0
    key = max(counts, key=lambda name: (counts[name], name))
    return key, int(counts[key])


def _modality_shape(mdata, modality: str | None) -> tuple[int, int]:
    if modality is None or modality not in mdata.mod:
        return 0, 0
    adata = mdata.mod[modality]
    return int(adata.n_obs), int(adata.n_vars)


def build_inventory_rows(specs: list[DatasetSpec], loaded: dict[str, object]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for spec in specs:
        mdata = loaded[spec.label]
        modalities = list(mdata.mod.keys())
        rna_name = find_modality(mdata, RNA_MODALITY_CANDIDATES, required=False)
        circ_name = find_modality(mdata, CIRC_MODALITY_CANDIDATES, required=False)
        rt_name = find_modality(mdata, RT_MODALITY_CANDIDATES, required=False)
        cnv_name = find_modality(mdata, CNV_MODALITY_CANDIDATES, required=False)
        aux_name = spec.expected_aux_modality or rt_name or cnv_name
        rna_cells, rna_features = _modality_shape(mdata, rna_name)
        circ_cells, circ_features = _modality_shape(mdata, circ_name)
        aux_cells, aux_features = _modality_shape(mdata, aux_name)
        annotated = total = 0
        recovery = 0.0
        source_counts: dict[str, int] = {}
        if circ_name is not None:
            circ = mdata.mod[circ_name]
            annotated, total, recovery = host_gene_recovery(circ)
            source_counts = host_gene_source_counts(circ)
        rows.append(
            {
                "dataset": spec.label,
                "source_accession": spec.source_accession,
                "protocol": spec.protocol,
                "object_path": str(spec.path),
                "modalities": ";".join(map(str, modalities)),
                "modality_n_obs": json_dumps({name: int(mdata.mod[name].n_obs) for name in modalities}),
                "modality_n_vars": json_dumps({name: int(mdata.mod[name].n_vars) for name in modalities}),
                "all_modality_cell_overlap": _all_modality_overlap(mdata, modalities),
                "all_modality_cell_union": _all_modality_union(mdata, modalities),
                "pairwise_cell_overlap": json_dumps(_pairwise_overlap(mdata, modalities)),
                "obs_columns_by_modality": json_dumps(
                    {name: list(map(str, mdata.mod[name].obs.columns)) for name in modalities}
                ),
                "var_columns_by_modality": json_dumps(
                    {name: list(map(str, mdata.mod[name].var.columns)) for name in modalities}
                ),
                "obsm_keys_by_modality": json_dumps(
                    {name: list(map(str, getattr(mdata.mod[name], "obsm", {}).keys())) for name in modalities}
                ),
                "layers_by_modality": json_dumps(
                    {name: list(map(str, getattr(mdata.mod[name], "layers", {}).keys())) for name in modalities}
                ),
                "uns_keys_by_modality": json_dumps(
                    {name: list(map(str, getattr(mdata.mod[name], "uns", {}).keys())) for name in modalities}
                ),
                "rna_cells": rna_cells,
                "rna_features": rna_features,
                "circ_cells": circ_cells,
                "circRNAs": circ_features,
                "aux_modality": aux_name or "",
                "aux_cells": aux_cells,
                "aux_features": aux_features,
                "circ_host_gene_annotated": annotated,
                "circ_host_gene_total": total,
                "circ_host_gene_annotation_rate": recovery,
                "circ_host_gene_source_counts": json_dumps(source_counts),
            }
        )
    return pd.DataFrame(rows)


def build_host_gene_annotation_summary(specs: list[DatasetSpec], loaded: dict[str, object]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for spec in specs:
        mdata = loaded[spec.label]
        circ_name = find_modality(mdata, CIRC_MODALITY_CANDIDATES)
        circ = mdata.mod[circ_name]
        annotated, total, recovery = host_gene_recovery(circ)
        source_counts = host_gene_source_counts(circ)
        dominant_source, dominant_source_count = _dominant_source(source_counts)
        host_values = circ.var["host_gene"] if "host_gene" in circ.var.columns else pd.Series([], dtype=object)
        host_genes: set[str] = set()
        for value in host_values:
            host_genes.update(split_host_genes(value))
        multi_host = 0
        if "host_gene_n" in circ.var.columns:
            multi_host = int((pd.to_numeric(circ.var["host_gene_n"], errors="coerce").fillna(0) > 1).sum())
        rows.append(
            {
                "dataset": spec.label,
                "object_path": str(spec.path),
                "circ_modality": circ_name,
                "n_circRNAs": int(total),
                "host_gene_annotated": annotated,
                "host_gene_unannotated": int(total - annotated),
                "host_gene_annotation_rate": recovery,
                "n_unique_host_genes": int(len(host_genes)),
                "n_multi_host_circRNAs": multi_host,
                "host_gene_source_counts": json_dumps(source_counts),
                "dominant_host_gene_source": dominant_source,
                "dominant_source_count": dominant_source_count,
                "required_host_gene_columns_present": json_dumps(
                    {column: column in circ.var.columns for column in HOST_GENE_COLUMNS}
                ),
            }
        )
    return pd.DataFrame(rows)


def _feature_coordinate_columns(circ) -> pd.DataFrame:
    var = circ.var.copy()
    out = pd.DataFrame({"feature_id": list(map(str, circ.var_names))})
    for column in ["circ_id", "chrom", "start", "end", "strand", "host_gene", "host_gene_source"]:
        out[column] = var[column].astype(str).values if column in var.columns else ""
    return out


def _spread_metrics(values: np.ndarray, coords: pd.DataFrame | None) -> dict[str, float]:
    if coords is None or coords.empty:
        return {
            "umap_detected_spread_ratio": np.nan,
            "umap_detected_area_fraction": np.nan,
            "umap_detected_centroid_x": np.nan,
            "umap_detected_centroid_y": np.nan,
        }
    detected = values > 0
    coord_array = coords[["umap_1", "umap_2"]].to_numpy(dtype=float)
    if int(detected.sum()) == 0:
        return {
            "umap_detected_spread_ratio": 0.0,
            "umap_detected_area_fraction": 0.0,
            "umap_detected_centroid_x": np.nan,
            "umap_detected_centroid_y": np.nan,
        }
    det_coords = coord_array[detected]
    centroid = det_coords.mean(axis=0)
    global_centroid = coord_array.mean(axis=0)
    detected_spread = float(np.mean(np.linalg.norm(det_coords - centroid, axis=1))) if len(det_coords) > 1 else 0.0
    global_spread = float(np.mean(np.linalg.norm(coord_array - global_centroid, axis=1)))
    det_range = np.ptp(det_coords, axis=0) if len(det_coords) > 1 else np.asarray([0.0, 0.0])
    all_range = np.ptp(coord_array, axis=0)
    det_area = float(det_range[0] * det_range[1])
    all_area = float(all_range[0] * all_range[1])
    return {
        "umap_detected_spread_ratio": detected_spread / global_spread if global_spread > 0 else np.nan,
        "umap_detected_area_fraction": det_area / all_area if all_area > 0 else np.nan,
        "umap_detected_centroid_x": float(centroid[0]),
        "umap_detected_centroid_y": float(centroid[1]),
    }


def rank_circ_features(circ, coords: pd.DataFrame | None = None) -> pd.DataFrame:
    matrix = _matrix_as_array(circ)
    n_cells = (matrix > 0).sum(axis=0).astype(int)
    total_support = matrix.sum(axis=0).astype(float)
    max_support = matrix.max(axis=0).astype(float) if matrix.size else np.zeros(circ.n_vars, dtype=float)
    mean_detected = np.divide(
        total_support,
        n_cells,
        out=np.zeros_like(total_support, dtype=float),
        where=n_cells > 0,
    )
    rows = _feature_coordinate_columns(circ)
    rows["n_cells_detected"] = n_cells
    rows["total_support"] = total_support
    rows["max_cell_support"] = max_support
    rows["mean_support_detected"] = mean_detected
    if coords is not None:
        for idx in range(circ.n_vars):
            for key, value in _spread_metrics(matrix[:, idx], coords).items():
                rows.loc[idx, key] = value
    else:
        rows["umap_detected_spread_ratio"] = np.nan
        rows["umap_detected_area_fraction"] = np.nan
        rows["umap_detected_centroid_x"] = np.nan
        rows["umap_detected_centroid_y"] = np.nan
    rows["score_n_cells_pct"] = _rank_pct(rows["n_cells_detected"])
    rows["score_total_support_pct"] = _rank_pct(rows["total_support"])
    rows["score_spread_pct"] = _rank_pct(rows["umap_detected_spread_ratio"].fillna(0.0))
    rows["combined_score"] = (
        0.50 * rows["score_n_cells_pct"]
        + 0.45 * rows["score_total_support_pct"]
        + 0.05 * rows["score_spread_pct"]
    )
    return rows.sort_values(
        ["combined_score", "n_cells_detected", "total_support", "feature_id"],
        ascending=[False, False, False, True],
    ).reset_index(drop=True)


def rank_host_gene_features(circ, coords: pd.DataFrame | None = None) -> pd.DataFrame:
    matrix = _matrix_as_array(circ)
    host_to_cols: dict[str, list[int]] = {}
    if "host_gene" not in circ.var.columns:
        return pd.DataFrame()
    for col_idx, value in enumerate(circ.var["host_gene"]):
        for host_gene in split_host_genes(value):
            host_to_cols.setdefault(host_gene, []).append(col_idx)
    rows: list[dict[str, object]] = []
    for host_gene, cols in sorted(host_to_cols.items()):
        host_matrix = matrix[:, cols]
        per_cell_support = host_matrix.sum(axis=1).astype(float)
        n_detected = int(np.count_nonzero(per_cell_support > 0))
        total_support = float(per_cell_support.sum())
        row: dict[str, object] = {
            "host_gene": host_gene,
            "n_circRNAs": int(len(cols)),
            "n_cells_detected": n_detected,
            "total_support": total_support,
            "max_cell_support": float(per_cell_support.max()) if len(per_cell_support) else 0.0,
            "mean_support_detected": float(total_support / n_detected) if n_detected else 0.0,
            "host_gene_source_counts": json_dumps(_source_counts_for_columns(circ, cols)),
            "circRNA_ids": ";".join(map(str, circ.var_names[cols])),
        }
        row.update(_spread_metrics(per_cell_support, coords))
        rows.append(row)
    frame = pd.DataFrame(rows)
    if frame.empty:
        return frame
    frame["score_n_cells_pct"] = _rank_pct(frame["n_cells_detected"])
    frame["score_total_support_pct"] = _rank_pct(frame["total_support"])
    frame["score_n_circRNAs_pct"] = _rank_pct(frame["n_circRNAs"])
    frame["score_spread_pct"] = _rank_pct(frame["umap_detected_spread_ratio"].fillna(0.0))
    frame["combined_score"] = (
        0.30 * frame["score_n_cells_pct"]
        + 0.55 * frame["score_total_support_pct"]
        + 0.10 * frame["score_n_circRNAs_pct"]
        + 0.05 * frame["score_spread_pct"]
    )
    return frame.sort_values(
        ["combined_score", "n_cells_detected", "total_support", "host_gene"],
        ascending=[False, False, False, True],
    ).reset_index(drop=True)


def _normalize_log_rna_matrix(rna) -> np.ndarray:
    x = rna.X
    if sparse.issparse(x):
        x = x.tocsr().astype(float)
        library = np.asarray(x.sum(axis=1)).ravel()
        scale = np.divide(1e4, library, out=np.zeros_like(library, dtype=float), where=library > 0)
        x = sparse.diags(scale) @ x
        x = x.tocsr()
        x.data = np.log1p(x.data)
        mean = np.asarray(x.mean(axis=0)).ravel()
        mean_sq = np.asarray(x.power(2).mean(axis=0)).ravel()
        variance = np.maximum(mean_sq - mean**2, 0.0)
        top_n = min(2000, x.shape[1])
        top_idx = np.argsort(variance)[-top_n:]
        arr = x[:, top_idx].toarray()
    else:
        arr = np.asarray(x, dtype=float)
        library = arr.sum(axis=1)
        scale = np.divide(1e4, library, out=np.zeros_like(library, dtype=float), where=library > 0)
        arr = np.log1p(arr * scale[:, None])
        variance = arr.var(axis=0)
        top_n = min(2000, arr.shape[1])
        top_idx = np.argsort(variance)[-top_n:]
        arr = arr[:, top_idx]
    arr = np.asarray(arr, dtype=float)
    arr -= arr.mean(axis=0)
    std = arr.std(axis=0)
    arr[:, std > 0] /= std[std > 0]
    return arr


def _pca_scores(arr: np.ndarray, n_components: int) -> np.ndarray:
    n_components = max(2, min(n_components, arr.shape[0] - 1, arr.shape[1]))
    u, s, _ = np.linalg.svd(arr, full_matrices=False)
    return u[:, :n_components] * s[:n_components]


def _pairwise_distances(arr: np.ndarray) -> np.ndarray:
    squared = np.sum(arr * arr, axis=1)
    d2 = squared[:, None] + squared[None, :] - 2.0 * (arr @ arr.T)
    return np.sqrt(np.maximum(d2, 0.0))


def _smooth_knn_dist(knn_dists: np.ndarray, n_neighbors: int) -> tuple[np.ndarray, np.ndarray]:
    target = math.log2(n_neighbors)
    rhos = np.zeros(knn_dists.shape[0], dtype=float)
    sigmas = np.ones(knn_dists.shape[0], dtype=float)
    for i, distances in enumerate(knn_dists):
        positive = distances[distances > 0]
        rhos[i] = float(positive.min()) if len(positive) else 0.0
        lo = 0.0
        hi = np.inf
        mid = 1.0
        for _ in range(64):
            psum = 0.0
            for distance in distances:
                shifted = distance - rhos[i]
                psum += 1.0 if shifted <= 0 else math.exp(-(shifted / mid))
            if abs(psum - target) < 1e-5:
                break
            if psum > target:
                hi = mid
                mid = (lo + hi) / 2.0
            else:
                lo = mid
                mid = mid * 2.0 if not np.isfinite(hi) else (lo + hi) / 2.0
        sigmas[i] = max(mid, 1e-6)
    return sigmas, rhos


def _fuzzy_graph(pcs: np.ndarray, n_neighbors: int) -> np.ndarray:
    n_obs = pcs.shape[0]
    n_neighbors = max(2, min(n_neighbors, n_obs - 1))
    distances = _pairwise_distances(pcs)
    order = np.argsort(distances, axis=1)[:, 1 : n_neighbors + 1]
    knn_dists = np.take_along_axis(distances, order, axis=1)
    sigmas, rhos = _smooth_knn_dist(knn_dists, n_neighbors)
    graph = np.zeros((n_obs, n_obs), dtype=float)
    for i in range(n_obs):
        for j, distance in zip(order[i], knn_dists[i]):
            shifted = distance - rhos[i]
            weight = 1.0 if shifted <= 0 else math.exp(-(shifted / sigmas[i]))
            graph[i, j] = max(graph[i, j], weight)
    return graph + graph.T - graph * graph.T


def _spectral_init(graph: np.ndarray, pcs: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    try:
        from scipy import linalg

        degree = graph.sum(axis=1)
        laplacian = np.diag(degree) - graph
        mass = np.diag(degree + 1e-6)
        _, vectors = linalg.eigh(laplacian, mass)
        coords = vectors[:, 1:3]
    except Exception:
        coords = pcs[:, :2].copy()
    coords = np.asarray(coords, dtype=float)
    coords -= coords.mean(axis=0)
    scale = np.std(coords, axis=0)
    coords[:, scale > 0] /= scale[scale > 0]
    coords += rng.normal(0.0, 1e-4, size=coords.shape)
    return coords


def _internal_umap_embedding(pcs: np.ndarray, *, random_state: int = 13) -> np.ndarray:
    rng = np.random.default_rng(random_state)
    n_obs = pcs.shape[0]
    graph = _fuzzy_graph(pcs, n_neighbors=min(15, n_obs - 1))
    coords = _spectral_init(graph, pcs, rng)
    edge_i, edge_j = np.where(np.triu(graph, k=1) > 0)
    weights = graph[edge_i, edge_j]
    if len(weights) == 0:
        return coords
    probabilities = weights / weights.sum()
    a = 1.5769434603113077
    b = 0.8950608781227859
    epochs = 350
    negative_samples = 3
    initial_alpha = 1.0
    for epoch in range(epochs):
        alpha = initial_alpha * (1.0 - epoch / epochs)
        sample_size = min(len(weights), max(200, len(weights) // 2))
        sampled = rng.choice(len(weights), size=sample_size, replace=True, p=probabilities)
        for edge_idx in sampled:
            i = int(edge_i[edge_idx])
            j = int(edge_j[edge_idx])
            delta = coords[i] - coords[j]
            dist2 = float(delta @ delta) + 1e-6
            grad_coeff = -2.0 * a * b * (dist2 ** (b - 1.0)) / (a * (dist2**b) + 1.0)
            grad = np.clip(grad_coeff * delta, -4.0, 4.0)
            coords[i] += alpha * grad * 0.01
            coords[j] -= alpha * grad * 0.01
            for _ in range(negative_samples):
                k = int(rng.integers(0, n_obs))
                if k == i:
                    continue
                delta_neg = coords[i] - coords[k]
                dist2_neg = float(delta_neg @ delta_neg) + 1e-6
                repulse = 2.0 * b / ((0.001 + dist2_neg) * (a * (dist2_neg**b) + 1.0))
                grad_neg = np.clip(repulse * delta_neg, -4.0, 4.0)
                coords[i] += alpha * grad_neg * 0.01
        coords -= coords.mean(axis=0)
    return coords


def smartseq3_embedding_table(mdata) -> tuple[pd.DataFrame, str]:
    rna = mdata.mod[find_modality(mdata, RNA_MODALITY_CANDIDATES)]
    circ = mdata.mod[find_modality(mdata, CIRC_MODALITY_CANDIDATES)]
    embedding_key = None
    embedding = None
    for key in ["X_umap", "umap", "X_rna_umap", "rna_umap"]:
        if key in getattr(rna, "obsm", {}):
            values = np.asarray(rna.obsm[key])
            if values.ndim == 2 and values.shape[1] >= 2:
                embedding_key = key
                embedding = values[:, :2].astype(float)
                break
    if embedding is None:
        normalized = _normalize_log_rna_matrix(rna)
        pcs = _pca_scores(normalized, n_components=30)
        embedding = _internal_umap_embedding(pcs)
        embedding_key = "computed_rna_umap"
        embedding_source = "computed_internal_exact_knn_umap"
    else:
        embedding_source = "object_rna_obsm"
    cells = list(map(str, rna.obs_names))
    frame = pd.DataFrame(
        {
            "cell_id": cells,
            "umap_1": embedding[:, 0],
            "umap_2": embedding[:, 1],
            "embedding_key": embedding_key,
            "embedding_source": embedding_source,
            "circRNA_count": get_or_compute_circ_metric(mdata, circ, "circRNA_count").reindex(cells).to_numpy(),
            "circRNA_total_support": get_or_compute_circ_metric(
                mdata, circ, "circRNA_total_support"
            ).reindex(cells).to_numpy(),
        }
    )
    for column in [
        "membership",
        "total_rna_counts",
        "detected_genes",
        "mitochondrial_fraction",
        "ribosomal_fraction",
    ]:
        if column in rna.obs.columns:
            frame[column] = rna.obs[column].reindex(cells).to_numpy()
    return frame, embedding_source


def build_top_host_genes_by_dataset(specs: list[DatasetSpec], loaded: dict[str, object], top_n: int) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for spec in specs:
        circ = loaded[spec.label].mod[find_modality(loaded[spec.label], CIRC_MODALITY_CANDIDATES)]
        frame = rank_host_gene_features(circ).head(top_n).copy()
        frame.insert(0, "rank", np.arange(1, len(frame) + 1))
        frame.insert(0, "dataset", spec.label)
        frames.append(frame)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def select_smartseq3_features(
    circ_candidates: pd.DataFrame,
    host_candidates: pd.DataFrame,
    n_cells: int,
) -> pd.DataFrame:
    min_individual_cells = max(5, int(math.ceil(0.03 * n_cells)))
    min_host_cells = max(5, int(math.ceil(0.03 * n_cells)))
    individual_pool = circ_candidates[
        (circ_candidates["n_cells_detected"] >= min_individual_cells)
        & (circ_candidates["host_gene"].fillna("").astype(str).str.len() > 0)
    ]
    individual = individual_pool.iloc[0] if not individual_pool.empty else circ_candidates.iloc[0]
    host_pool = host_candidates[host_candidates["n_cells_detected"] >= min_host_cells]
    host = host_pool.iloc[0] if not host_pool.empty else host_candidates.iloc[0]
    backup_pool = circ_candidates[circ_candidates["feature_id"] != individual["feature_id"]]
    backup = backup_pool.iloc[0] if not backup_pool.empty else individual
    rows = [
        {
            "role": "individual_circRNA",
            "feature_type": "circRNA",
            "feature_id": individual["feature_id"],
            "host_gene": individual.get("host_gene", ""),
            "n_cells_detected": int(individual["n_cells_detected"]),
            "total_support": float(individual["total_support"]),
            "combined_score": float(individual["combined_score"]),
            "umap_detected_spread_ratio": float(individual["umap_detected_spread_ratio"]),
            "umap_detected_area_fraction": float(individual["umap_detected_area_fraction"]),
            "selection_note": (
                f"Top individual circRNA by combined recurrence/support score with at least "
                f"{min_individual_cells} detected cells."
            ),
        },
        {
            "role": "host_gene_signal",
            "feature_type": "host_gene",
            "feature_id": host["host_gene"],
            "host_gene": host["host_gene"],
            "n_cells_detected": int(host["n_cells_detected"]),
            "total_support": float(host["total_support"]),
            "combined_score": float(host["combined_score"]),
            "umap_detected_spread_ratio": float(host["umap_detected_spread_ratio"]),
            "umap_detected_area_fraction": float(host["umap_detected_area_fraction"]),
            "selection_note": (
                f"Top host-gene-level circRNA signal by recurrence/support/feature-count score "
                f"with at least {min_host_cells} detected cells."
            ),
        },
        {
            "role": "backup_candidate",
            "feature_type": "circRNA",
            "feature_id": backup["feature_id"],
            "host_gene": backup.get("host_gene", ""),
            "n_cells_detected": int(backup["n_cells_detected"]),
            "total_support": float(backup["total_support"]),
            "combined_score": float(backup["combined_score"]),
            "umap_detected_spread_ratio": float(backup["umap_detected_spread_ratio"]),
            "umap_detected_area_fraction": float(backup["umap_detected_area_fraction"]),
            "selection_note": "Second-ranked individual circRNA candidate by the same data-driven score.",
        },
    ]
    return pd.DataFrame(rows)


def build_detector_status_table() -> pd.DataFrame:
    try:
        from circyto.cli.doctor_meta import DETECTOR_SPECS, detector_runtime_status
        from circyto.detectors import build_default_engines, get_detector_capabilities
    except Exception as exc:
        return pd.DataFrame(
            [
                {
                    "detector": "unavailable",
                    "status": "error",
                    "available": False,
                    "input_modes": "",
                    "recommended_execution_mode": "",
                    "validation_status": "",
                    "missing_dependencies": "",
                    "reason": f"{type(exc).__name__}: {exc}",
                }
            ]
        )
    engines = build_default_engines()
    rows: list[dict[str, object]] = []
    for spec in DETECTOR_SPECS:
        runtime = detector_runtime_status(spec.name)
        engine = engines.get(spec.name)
        caps = get_detector_capabilities(engine) if engine is not None else None
        input_modes: list[str] = []
        if caps is not None:
            if caps.accepts_fastq:
                input_modes.append("fastq")
            if caps.accepts_alignment:
                input_modes.append("alignment")
        missing = list(runtime.get("details", {}).get("missing_cmds", []))
        missing.extend(runtime.get("details", {}).get("missing_assets", []))
        rows.append(
            {
                "detector": spec.name,
                "status": runtime["status"],
                "available": runtime["status"] == "READY",
                "type": spec.det_type,
                "input_modes": ";".join(input_modes),
                "recommended_execution_mode": caps.recommended_execution_mode if caps else "",
                "validation_status": spec.validation_status,
                "missing_dependencies": ";".join(map(str, sorted(set(missing)))),
                "reason": runtime["reason"],
            }
        )
    return pd.DataFrame(rows)


def build_consistency_checks(inventory: pd.DataFrame, smart_cells: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []

    def add(dataset: str, metric: str, manuscript_value: object, object_value: object, status: str) -> None:
        rows.append(
            {
                "dataset": dataset,
                "metric": metric,
                "manuscript_value": manuscript_value,
                "object_value": object_value,
                "status": status,
            }
        )

    for _, row in inventory.iterrows():
        dataset = str(row["dataset"])
        expected = EXPECTED_MANUSCRIPT_NUMBERS[dataset]
        if dataset == "Smart-seq3":
            actual_cells = int(row["all_modality_cell_overlap"])
            add(dataset, "cells", expected["cells"], actual_cells, "match" if actual_cells == expected["cells"] else "mismatch")
            actual_rna = int(row["rna_features"])
            add(
                dataset,
                "RNA features",
                expected["rna_features"],
                actual_rna,
                "match" if actual_rna == expected["rna_features"] else "mismatch",
            )
            actual_circ = int(row["circRNAs"])
            add(
                dataset,
                "circRNAs",
                expected["circRNAs"],
                actual_circ,
                "match" if actual_circ == expected["circRNAs"] else "mismatch",
            )
            host_pct = 100.0 * float(row["circ_host_gene_annotation_rate"])
            add(
                dataset,
                "host annotation rate",
                f"{expected['host_annotation_rate_pct_rounded']:.1f}%",
                f"{host_pct:.1f}% ({int(row['circ_host_gene_annotated'])}/{int(row['circ_host_gene_total'])})",
                "match" if round(host_pct, 1) == expected["host_annotation_rate_pct_rounded"] else "mismatch",
            )
            add(
                dataset,
                "median circRNAs/cell",
                "not documented",
                float(smart_cells["circRNA_count"].median()),
                "needs_manuscript_value",
            )
            add(
                dataset,
                "median support/cell",
                "not documented",
                float(smart_cells["circRNA_total_support"].median()),
                "needs_manuscript_value",
            )
        elif dataset == "HAP1":
            actual_cells = int(row["rna_cells"])
            add(dataset, "cells", expected["cells"], actual_cells, "match" if actual_cells == expected["cells"] else "mismatch")
            actual_circ = int(row["circRNAs"])
            add(
                dataset,
                "circRNAs",
                expected["circRNAs"],
                actual_circ,
                "match" if actual_circ == expected["circRNAs"] else "mismatch",
            )
            actual_overlap = int(row["all_modality_cell_overlap"])
            add(
                dataset,
                "RT overlap",
                expected["rt_overlap"],
                actual_overlap,
                "match" if actual_overlap == expected["rt_overlap"] else "mismatch",
            )
            host_pct = 100.0 * float(row["circ_host_gene_annotation_rate"])
            add(
                dataset,
                "host annotation rate",
                f"{expected['host_annotation_rate_pct_rounded']:.1f}%",
                f"{host_pct:.1f}% ({int(row['circ_host_gene_annotated'])}/{int(row['circ_host_gene_total'])})",
                "match" if round(host_pct, 1) == expected["host_annotation_rate_pct_rounded"] else "mismatch",
            )
        elif dataset == "IMR90":
            actual_cells = int(row["rna_cells"])
            add(dataset, "cells", expected["cells"], actual_cells, "match" if actual_cells == expected["cells"] else "mismatch")
            actual_circ = int(row["circRNAs"])
            add(
                dataset,
                "circRNAs",
                expected["circRNAs"],
                actual_circ,
                "match" if actual_circ == expected["circRNAs"] else "mismatch",
            )
            actual_overlap = int(row["all_modality_cell_overlap"])
            add(
                dataset,
                "CNV overlap",
                expected["cnv_overlap"],
                actual_overlap,
                "match" if actual_overlap == expected["cnv_overlap"] else "mismatch",
            )
            host_pct = 100.0 * float(row["circ_host_gene_annotation_rate"])
            add(
                dataset,
                "host annotation rate",
                f"{expected['host_annotation_rate_pct_rounded']:.1f}%",
                f"{host_pct:.1f}% ({int(row['circ_host_gene_annotated'])}/{int(row['circ_host_gene_total'])})",
                "match" if round(host_pct, 1) == expected["host_annotation_rate_pct_rounded"] else "mismatch",
            )
    return pd.DataFrame(rows)


def _markdown_table(frame: pd.DataFrame) -> str:
    if frame.empty:
        return ""
    columns = [str(column) for column in frame.columns]

    def cell(value: object) -> str:
        text = "" if pd.isna(value) else str(value)
        return text.replace("|", "\\|").replace("\n", "<br>")

    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join(["---"] * len(columns)) + " |",
    ]
    for _, row in frame.iterrows():
        lines.append("| " + " | ".join(cell(row[column]) for column in frame.columns) + " |")
    return "\n".join(lines)


def write_consistency_report(
    path: Path,
    specs: list[DatasetSpec],
    inventory: pd.DataFrame,
    host_summary: pd.DataFrame,
    checks: pd.DataFrame,
    selected: pd.DataFrame,
    embedding_source: str,
) -> None:
    mismatch_count = int((checks["status"] == "mismatch").sum())
    needs_value_count = int((checks["status"] == "needs_manuscript_value").sum())
    lines: list[str] = [
        "# Manuscript Consistency Report",
        "",
        "## Objects inspected",
        "",
    ]
    for spec in specs:
        lines.append(f"- {spec.label}: `{spec.path}`")
    lines.extend(
        [
            "",
            "## Object inventory",
            "",
            _markdown_table(
                inventory[
                    [
                        "dataset",
                        "modalities",
                        "modality_n_obs",
                        "modality_n_vars",
                        "all_modality_cell_overlap",
                        "pairwise_cell_overlap",
                    ]
                ]
            ),
            "",
            "## Host-gene annotation",
            "",
            _markdown_table(
                host_summary[
                    [
                        "dataset",
                        "n_circRNAs",
                        "host_gene_annotated",
                        "host_gene_unannotated",
                        "host_gene_annotation_rate",
                        "n_unique_host_genes",
                        "host_gene_source_counts",
                    ]
                ]
            ),
            "",
            "## Manuscript number checks",
            "",
            _markdown_table(checks),
            "",
            "## Smart-seq3 Figure 1 data",
            "",
            f"- RNA UMAP in object: `{'no' if embedding_source != 'object_rna_obsm' else 'yes'}`",
            f"- UMAP coordinate source: `{embedding_source}`",
            "- Final figure design and plotting were not performed.",
            "",
            "## Candidate Figure 1 panel C features",
            "",
            _markdown_table(selected),
            "",
            "## Metadata columns inspected",
            "",
        ]
    )
    for _, row in inventory.iterrows():
        lines.append(f"### {row['dataset']}")
        lines.append("")
        lines.append(f"- obs columns by modality: `{row['obs_columns_by_modality']}`")
        lines.append(f"- var columns by modality: `{row['var_columns_by_modality']}`")
        lines.append("")
    if mismatch_count == 0:
        lines.append("## Summary")
        lines.append("")
        lines.append("No documented count mismatches were found against the inspected objects.")
    else:
        lines.append("## Summary")
        lines.append("")
        lines.append(f"Found {mismatch_count} documented count mismatch(es).")
    if needs_value_count:
        lines.append(
            f"{needs_value_count} Smart-seq3 median value(s) were not documented in the current manuscript markdown and should be inserted from this report/table if needed."
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Prepare manuscript inventory, Smart-seq3 figure-ready tables, and consistency outputs."
    )
    parser.add_argument("--smartseq3", type=Path, default=DEFAULT_SPECS[0].path, help="Smart-seq3 RNA+circ .h5mu.")
    parser.add_argument("--hap1", type=Path, default=DEFAULT_SPECS[1].path, help="HAP1 RNA+circ+RT .h5mu.")
    parser.add_argument("--imr90", type=Path, default=DEFAULT_SPECS[2].path, help="IMR90 RNA+circ+CNV .h5mu.")
    parser.add_argument("--outdir", type=Path, default=Path("manuscript/results"), help="Output directory.")
    parser.add_argument("--top-n", type=int, default=100, help="Top feature rows to keep per dataset.")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    specs = [
        DatasetSpec("Smart-seq3", "E-MTAB-8735", "Smart-seq3", args.smartseq3),
        DatasetSpec("HAP1", "GSE278952", "scRR-seq RNA+circ+RT", args.hap1, "rt"),
        DatasetSpec("IMR90", "GSE278958", "scRR-seq RNA+circ+CNV", args.imr90, "cnv"),
    ]
    loaded = {spec.label: load_mudata(spec.path) for spec in specs}
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    inventory = build_inventory_rows(specs, loaded)
    host_summary = build_host_gene_annotation_summary(specs, loaded)
    top_hosts = build_top_host_genes_by_dataset(specs, loaded, top_n=args.top_n)

    smart_mdata = loaded["Smart-seq3"]
    smart_cells, embedding_source = smartseq3_embedding_table(smart_mdata)
    smart_circ = smart_mdata.mod[find_modality(smart_mdata, CIRC_MODALITY_CANDIDATES)]
    smart_circ_candidates = rank_circ_features(smart_circ, smart_cells).head(args.top_n)
    smart_host_candidates = rank_host_gene_features(smart_circ, smart_cells).head(args.top_n)
    selected = select_smartseq3_features(smart_circ_candidates, smart_host_candidates, len(smart_cells))
    detector_status = build_detector_status_table()
    checks = build_consistency_checks(inventory, smart_cells)

    write_table(inventory, outdir / "mudata_inventory.tsv")
    write_table(host_summary, outdir / "host_gene_annotation_summary.tsv")
    write_table(top_hosts, outdir / "top_host_genes_by_dataset.tsv")
    write_table(smart_cells, outdir / "smartseq3_umap_cells.tsv")
    write_table(smart_circ_candidates, outdir / "smartseq3_top_circRNA_candidates.tsv")
    write_table(smart_host_candidates, outdir / "smartseq3_top_hostgene_features.tsv")
    write_table(selected, outdir / "smartseq3_selected_feature_candidates.tsv")
    write_table(inventory, outdir / "supplement_table_1_dataset_inventory.tsv")
    write_table(detector_status, outdir / "supplement_table_2_detector_status.tsv")
    write_table(top_hosts, outdir / "supplement_table_3_top_host_genes.tsv")
    write_table(selected, outdir / "supplement_table_4_selected_smartseq3_features.tsv")
    write_consistency_report(
        outdir / "manuscript_consistency_report.md",
        specs,
        inventory,
        host_summary,
        checks,
        selected,
        embedding_source,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
