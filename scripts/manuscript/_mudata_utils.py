from __future__ import annotations

import json
import re
import sys
from pathlib import Path
from typing import Callable, Iterable, Sequence

import numpy as np
import pandas as pd

from circyto.multimodal.sync import read_h5mu
from scipy import sparse


CIRC_MODALITY_CANDIDATES = ("circ", "circrna", "circRNA", "circ_rna")
RNA_MODALITY_CANDIDATES = ("rna", "gene", "gene_expression")
RT_MODALITY_CANDIDATES = ("rt", "replication_timing", "replication_state")
CNV_MODALITY_CANDIDATES = ("cnv", "copy_number")


def fail(message: str) -> None:
    raise SystemExit(f"[ERROR] {message}")


def load_mudata(path: Path):
    try:
        import mudata as mu
    except ImportError as exc:
        fail(
            "mudata is required to read .h5mu files. Install circyto[mudata] "
            "or run `python -m pip install mudata`."
        )
        raise exc
    if not path.exists():
        fail(f"MuData file does not exist: {path}")
    return read_h5mu(path)


def find_modality(mdata, candidates: Sequence[str], *, required: bool = True) -> str | None:
    keys = list(mdata.mod.keys())
    for candidate in candidates:
        if candidate in keys:
            return candidate
    lowered = {str(key).lower(): key for key in keys}
    for candidate in candidates:
        match = lowered.get(candidate.lower())
        if match is not None:
            return match
    if required:
        fail(
            "Missing required modality. Expected one of "
            f"{', '.join(candidates)}; found {', '.join(map(str, keys)) or 'none'}."
        )
    return None


def dataset_name_from_path(path: Path) -> str:
    name = path.name
    for suffix in (".hostgene_fixed.h5mu", ".h5mu"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def matrix_sum_per_obs(adata) -> np.ndarray:
    values = adata.X.sum(axis=1)
    return np.asarray(values).ravel().astype(float)


def matrix_nnz_per_obs(adata) -> np.ndarray:
    if sparse.issparse(adata.X):
        return np.asarray(adata.X.getnnz(axis=1)).ravel().astype(float)
    return np.count_nonzero(np.asarray(adata.X), axis=1).astype(float)


def dense_matrix(adata) -> np.ndarray:
    if sparse.issparse(adata.X):
        return adata.X.toarray()
    return np.asarray(adata.X)


def metric_from_obs(mdata, adata, field: str) -> pd.Series | None:
    if field in mdata.obs.columns:
        return pd.to_numeric(mdata.obs[field], errors="coerce")
    if field in adata.obs.columns:
        return pd.to_numeric(adata.obs[field], errors="coerce")
    return None


def get_or_compute_rna_metric(mdata, rna, field: str) -> pd.Series:
    obs_metric = metric_from_obs(mdata, rna, field)
    if obs_metric is not None:
        return obs_metric
    if field == "detected_genes":
        return pd.Series(matrix_nnz_per_obs(rna), index=rna.obs_names, name=field)
    fail(
        f"Missing required RNA/cell metric `{field}` in mdata.obs or rna.obs, "
        "and no supported fallback is available."
    )


def get_or_compute_circ_metric(mdata, circ, field: str) -> pd.Series:
    obs_metric = metric_from_obs(mdata, circ, field)
    if obs_metric is not None:
        return obs_metric
    if field == "circRNA_count":
        return pd.Series(matrix_nnz_per_obs(circ), index=circ.obs_names, name=field)
    if field == "circRNA_total_support":
        return pd.Series(matrix_sum_per_obs(circ), index=circ.obs_names, name=field)
    fail(
        f"Missing required circRNA/cell metric `{field}` in mdata.obs or circ.obs, "
        "and no supported fallback is available."
    )


def fraction_by_predicate(adata, predicate: Callable[[np.ndarray], np.ndarray]) -> pd.Series:
    arr = dense_matrix(adata)
    if arr.ndim != 2:
        fail("Expected a two-dimensional modality matrix.")
    mask = predicate(arr)
    return pd.Series(mask.mean(axis=1).astype(float), index=adata.obs_names)


def get_or_compute_rt_metric(mdata, rt, field: str) -> pd.Series:
    obs_metric = metric_from_obs(mdata, rt, field)
    if obs_metric is not None:
        return obs_metric
    arr = dense_matrix(rt)
    values = set(np.unique(arr[~pd.isna(arr)]).tolist())
    if field == "frac_rt_pos":
        return fraction_by_predicate(rt, lambda x: x > 0).rename(field)
    if field == "frac_rt_neg":
        if values and values.issubset({0, 1, 0.0, 1.0}):
            return fraction_by_predicate(rt, lambda x: x == 0).rename(field)
        return fraction_by_predicate(rt, lambda x: x < 0).rename(field)
    fail(
        f"Missing required RT/cell metric `{field}` in mdata.obs or rt.obs, "
        "and no supported fallback is available."
    )


def get_or_compute_cnv_metric(mdata, cnv, field: str) -> pd.Series:
    obs_metric = metric_from_obs(mdata, cnv, field)
    if obs_metric is not None:
        return obs_metric
    if field == "frac_non_diploid":
        return fraction_by_predicate(cnv, lambda x: x != 2).rename(field)
    if field == "frac_loss":
        return fraction_by_predicate(cnv, lambda x: x < 2).rename(field)
    if field == "frac_gain":
        return fraction_by_predicate(cnv, lambda x: x > 2).rename(field)
    fail(
        f"Missing required CNV/cell metric `{field}` in mdata.obs or cnv.obs, "
        "and no supported fallback is available."
    )


def common_obs_names(*adatas) -> list[str]:
    if not adatas:
        return []
    shared = set(map(str, adatas[0].obs_names))
    for adata in adatas[1:]:
        shared &= set(map(str, adata.obs_names))
    return [cell for cell in map(str, adatas[0].obs_names) if cell in shared]


def nonempty_string_mask(series: pd.Series) -> pd.Series:
    text = series.fillna("").astype(str).str.strip()
    return (text != "") & (~text.str.lower().isin({"nan", "none", "na"}))


def split_host_genes(value: object) -> set[str]:
    if value is None or pd.isna(value):
        return set()
    text = str(value).strip()
    if not text or text.lower() in {"nan", "none", "na"}:
        return set()
    parts = re.split(r"[;,|]", text)
    return {part.strip() for part in parts if part.strip()}


def host_gene_set(circ) -> set[str]:
    if "host_gene" not in circ.var.columns:
        return set()
    genes: set[str] = set()
    for value in circ.var["host_gene"]:
        genes.update(split_host_genes(value))
    return genes


def host_gene_recovery(circ) -> tuple[int, int, float]:
    total = int(circ.n_vars)
    if "host_gene" not in circ.var.columns or total == 0:
        return 0, total, 0.0
    annotated = int(nonempty_string_mask(circ.var["host_gene"]).sum())
    return annotated, total, annotated / total if total else 0.0


def host_gene_source_counts(circ) -> dict[str, int]:
    if "host_gene_source" not in circ.var.columns:
        return {}
    text = circ.var["host_gene_source"].fillna("").astype(str).str.strip()
    text = text.where(text != "", other="missing")
    return {str(k): int(v) for k, v in text.value_counts().sort_index().items()}


def json_dumps(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def write_table(df: pd.DataFrame, path: Path | None, *, sep: str = "\t") -> None:
    if path is None or str(path) == "-":
        df.to_csv(sys.stdout, sep=sep, index=False)
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix.lower() == ".csv":
        df.to_csv(path, index=False)
    else:
        df.to_csv(path, sep=sep, index=False)


def infer_table_sep(path: Path) -> str:
    return "," if path.suffix.lower() == ".csv" else "\t"


def read_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep=infer_table_sep(path), keep_default_na=False)


def first_existing_column(frame: pd.DataFrame, candidates: Iterable[str]) -> str | None:
    lowered = {str(col).lower(): col for col in frame.columns}
    for candidate in candidates:
        if candidate in frame.columns:
            return candidate
        match = lowered.get(candidate.lower())
        if match is not None:
            return match
    return None


def parse_interval_label(label: object) -> tuple[str, int, int] | None:
    match = re.match(r"^([^:]+):(\d+)-(\d+)$", str(label))
    if not match:
        return None
    chrom, start, end = match.groups()
    return chrom, int(start), int(end)


def feature_coordinates(adata, *, kind: str) -> pd.DataFrame | None:
    var = adata.var.copy()
    chrom_col = first_existing_column(var, ["chrom", "chr", "chromosome", "seqname", "seqnames"])
    start_col = first_existing_column(var, ["start", "start_0based", "start_1based", "donor", "bin_start"])
    end_col = first_existing_column(var, ["end", "stop", "acceptor", "bin_end"])
    if chrom_col and start_col and end_col:
        out = pd.DataFrame(
            {
                "feature_position": np.arange(adata.n_vars, dtype=int),
                "feature_id": list(map(str, adata.var_names)),
                "chrom": var[chrom_col].astype(str).values,
                "start": pd.to_numeric(var[start_col], errors="coerce").values,
                "end": pd.to_numeric(var[end_col], errors="coerce").values,
            }
        )
    else:
        parsed = [parse_interval_label(name) for name in adata.var_names]
        if any(value is None for value in parsed):
            return None
        out = pd.DataFrame(
            {
                "feature_position": np.arange(adata.n_vars, dtype=int),
                "feature_id": list(map(str, adata.var_names)),
                "chrom": [value[0] for value in parsed if value is not None],
                "start": [value[1] for value in parsed if value is not None],
                "end": [value[2] for value in parsed if value is not None],
            }
        )
    out = out.dropna(subset=["chrom", "start", "end"]).copy()
    out["start"] = out["start"].astype(int)
    out["end"] = out["end"].astype(int)
    if kind == "circ":
        out["midpoint"] = ((out["start"] + out["end"]) / 2.0).astype(float)
    return out


def correlation_rows(frame: pd.DataFrame, pairs: Sequence[tuple[str, str]]) -> list[dict[str, object]]:
    from scipy import stats

    rows: list[dict[str, object]] = []
    for x_name, y_name in pairs:
        subset = frame[[x_name, y_name]].dropna()
        n = int(len(subset))
        row: dict[str, object] = {"x": x_name, "y": y_name, "n": n}
        if n < 3 or subset[x_name].nunique() < 2 or subset[y_name].nunique() < 2:
            row.update(
                {
                    "pearson_r": np.nan,
                    "pearson_p": np.nan,
                    "spearman_r": np.nan,
                    "spearman_p": np.nan,
                    "status": "insufficient_variation",
                }
            )
        else:
            pearson = stats.pearsonr(subset[x_name], subset[y_name])
            spearman = stats.spearmanr(subset[x_name], subset[y_name])
            row.update(
                {
                    "pearson_r": float(pearson.statistic),
                    "pearson_p": float(pearson.pvalue),
                    "spearman_r": float(spearman.statistic),
                    "spearman_p": float(spearman.pvalue),
                    "status": "ok",
                }
            )
        rows.append(row)
    return rows
