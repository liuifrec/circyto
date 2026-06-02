from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from circyto import __version__
from circyto.pipeline.workflow_reporting import (
    apply_standard_provenance,
    sanitize_frame_for_anndata,
    utc_now_iso,
    write_json,
)

try:
    import anndata as ad  # type: ignore

    HAS_ANNDATA = True
except Exception:
    HAS_ANNDATA = False


GENOMIC_COLUMNS = ("seqname", "start", "end")
STATE_RE = re.compile(r"^(?P<copy_number>\d+)-somy$")

CELL_MAPPING_DNA_COLUMNS = ("dna_cell_id", "dna_sample_id", "dna_id", "DNA")
CELL_MAPPING_RNA_COLUMNS = ("rna_cell_id", "rna_sample_id", "rna_id", "RNA")
CELL_MAPPING_CANONICAL_COLUMNS = ("canonical_cell_id", "scrr_cell_id", "cell_id")


def canonical_scrr_cell_id(sample_id: str) -> str:
    value = str(sample_id).strip()
    for prefix in ("DNA_", "RNA_"):
        if value.startswith(prefix):
            return value[len(prefix) :]
    return value


def infer_paired_rna_cell_id(dna_cell_id: str) -> str:
    value = str(dna_cell_id).strip()
    if value.startswith("DNA_"):
        return "RNA_" + value[len("DNA_") :]
    return value


def _find_column(df: pd.DataFrame, candidates: tuple[str, ...]) -> str | None:
    lookup = {str(column).strip().lower(): str(column) for column in df.columns}
    for candidate in candidates:
        match = lookup.get(candidate.lower())
        if match is not None:
            return match
    return None


def _validate_unique(values: list[str], *, label: str) -> None:
    seen: set[str] = set()
    duplicates: list[str] = []
    for value in values:
        if value in seen:
            duplicates.append(value)
        seen.add(value)
    if duplicates:
        preview = ", ".join(sorted(set(duplicates))[:5])
        raise ValueError(f"Duplicate {label} identifiers: {preview}")


def _read_cnv_wide_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", keep_default_na=False, compression="infer")
    df = df.rename(columns={column: str(column).strip() for column in df.columns})
    missing = [column for column in GENOMIC_COLUMNS if column not in df.columns[:3].tolist()]
    if missing:
        raise ValueError(
            f"{path} must start with genomic columns seqname, start, end; missing: {', '.join(missing)}"
        )
    if df.shape[1] <= len(GENOMIC_COLUMNS):
        raise ValueError(f"{path} has no DNA cell columns after seqname/start/end.")
    for column in ("start", "end"):
        try:
            df[column] = pd.to_numeric(df[column], errors="raise").astype(np.int64)
        except Exception as exc:
            raise ValueError(f"{path} column {column} must be integer genomic coordinates.") from exc
    return df


def _make_bins(df: pd.DataFrame) -> pd.DataFrame:
    bins = df.loc[:, list(GENOMIC_COLUMNS)].copy()
    bins["seqname"] = bins["seqname"].astype(str)
    bins["start"] = bins["start"].astype(np.int64)
    bins["end"] = bins["end"].astype(np.int64)
    bins["bin_size"] = bins["end"] - bins["start"]
    if (bins["bin_size"] <= 0).any():
        raise ValueError("CNV bins must have end > start for every row.")
    bins["bin_id"] = bins["seqname"] + ":" + bins["start"].astype(str) + "-" + bins["end"].astype(str)
    _validate_unique(bins["bin_id"].astype(str).tolist(), label="CNV bin")
    return bins.loc[:, ["bin_id", "seqname", "start", "end", "bin_size"]]


def _copy_number_from_state(value: object) -> int:
    text = str(value).strip()
    match = STATE_RE.match(text)
    if match is None:
        raise ValueError(f"CNV state value must match '<integer>-somy'; observed {text!r}")
    return int(match.group("copy_number"))


def _state_matrix_to_copy_number(df: pd.DataFrame, cell_columns: list[str]) -> np.ndarray:
    parsed = df[cell_columns].apply(lambda column: column.map(_copy_number_from_state))
    return parsed.to_numpy(dtype=np.int16).T


def _numeric_matrix(df: pd.DataFrame, cell_columns: list[str], *, source: Path) -> np.ndarray:
    numeric = df[cell_columns].apply(pd.to_numeric, errors="raise")
    if numeric.isna().any().any():
        raise ValueError(f"{source} contains missing numeric CNV values.")
    return numeric.to_numpy(dtype=np.float32).T


def _validate_matching_bins(states_bins: pd.DataFrame, other_bins: pd.DataFrame, *, other_path: Path) -> None:
    if not states_bins.loc[:, ["seqname", "start", "end"]].equals(other_bins.loc[:, ["seqname", "start", "end"]]):
        raise ValueError(f"{other_path} genomic bins do not match the CNV states table.")


def _load_cell_mapping(path: Path | None) -> dict[str, dict[str, str]]:
    if path is None:
        return {}
    df = pd.read_csv(path, sep="\t", keep_default_na=False)
    dna_col = _find_column(df, CELL_MAPPING_DNA_COLUMNS)
    if dna_col is None:
        raise ValueError(
            f"{path} must include a DNA cell column. Accepted names: {', '.join(CELL_MAPPING_DNA_COLUMNS)}"
        )
    rna_col = _find_column(df, CELL_MAPPING_RNA_COLUMNS)
    canonical_col = _find_column(df, CELL_MAPPING_CANONICAL_COLUMNS)
    _validate_unique(df[dna_col].astype(str).tolist(), label="DNA mapping")
    mapping: dict[str, dict[str, str]] = {}
    for _, row in df.iterrows():
        dna_cell_id = str(row[dna_col]).strip()
        rna_cell_id = str(row[rna_col]).strip() if rna_col is not None else infer_paired_rna_cell_id(dna_cell_id)
        canonical_cell_id = (
            str(row[canonical_col]).strip()
            if canonical_col is not None
            else canonical_scrr_cell_id(rna_cell_id or dna_cell_id)
        )
        mapping[dna_cell_id] = {
            "rna_cell_id": rna_cell_id,
            "canonical_cell_id": canonical_cell_id,
        }
    return mapping


def _build_obs(
    dna_cell_ids: list[str],
    *,
    cell_mapping_path: Path | None,
    obs_id_strategy: str,
) -> pd.DataFrame:
    mapping = _load_cell_mapping(cell_mapping_path)
    rows: list[dict[str, str]] = []
    missing_mapped = sorted(set(mapping) - set(dna_cell_ids))
    if missing_mapped:
        preview = ", ".join(missing_mapped[:5])
        raise ValueError(f"Cell mapping contains DNA IDs absent from CNV states table: {preview}")
    for dna_cell_id in dna_cell_ids:
        mapped = mapping.get(dna_cell_id, {})
        canonical_cell_id = mapped.get("canonical_cell_id", canonical_scrr_cell_id(dna_cell_id))
        rna_cell_id = mapped.get("rna_cell_id", infer_paired_rna_cell_id(dna_cell_id))
        if obs_id_strategy == "canonical":
            cell_id = canonical_cell_id
        elif obs_id_strategy == "dna":
            cell_id = dna_cell_id
        elif obs_id_strategy == "rna":
            cell_id = rna_cell_id
        else:
            raise ValueError(f"Unsupported obs_id_strategy: {obs_id_strategy}")
        rows.append(
            {
                "cell_id": cell_id,
                "canonical_cell_id": canonical_cell_id,
                "dna_cell_id": dna_cell_id,
                "rna_cell_id": rna_cell_id,
                "pairing_strategy": "cell_mapping" if dna_cell_id in mapping else "strip_modality_prefix",
            }
        )
    obs = pd.DataFrame(rows)
    _validate_unique(obs["cell_id"].astype(str).tolist(), label="CNV obs")
    obs = obs.set_index("cell_id", drop=False)
    obs.index.name = None
    return obs


def _safe_uns(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _safe_uns(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_safe_uns(item) for item in value]
    if isinstance(value, tuple):
        return [_safe_uns(item) for item in value]
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return str(value)


def import_scrr_cnv(
    *,
    cnv_states_path: Path,
    outdir: Path,
    cnv_mappabilitynorm_path: Path | None = None,
    cell_mapping_path: Path | None = None,
    obs_id_strategy: str = "canonical",
    write_h5ad: bool = True,
) -> dict[str, Any]:
    states_df = _read_cnv_wide_table(cnv_states_path)
    bins_df = _make_bins(states_df)
    dna_cell_ids = [str(column) for column in states_df.columns[len(GENOMIC_COLUMNS) :]]
    _validate_unique(dna_cell_ids, label="DNA cell")
    copy_number = _state_matrix_to_copy_number(states_df, dna_cell_ids)

    obs_df = _build_obs(
        dna_cell_ids,
        cell_mapping_path=cell_mapping_path,
        obs_id_strategy=obs_id_strategy,
    )
    layer_paths: dict[str, str] = {}
    mappabilitynorm = None
    if cnv_mappabilitynorm_path is not None:
        mappability_df = _read_cnv_wide_table(cnv_mappabilitynorm_path)
        mappability_bins = _make_bins(mappability_df)
        _validate_matching_bins(bins_df, mappability_bins, other_path=cnv_mappabilitynorm_path)
        mappability_cells = [str(column) for column in mappability_df.columns[len(GENOMIC_COLUMNS) :]]
        if set(mappability_cells) != set(dna_cell_ids):
            raise ValueError(f"{cnv_mappabilitynorm_path} cell columns do not match the CNV states table.")
        mappability_df = mappability_df.loc[:, list(GENOMIC_COLUMNS) + dna_cell_ids]
        mappabilitynorm = _numeric_matrix(mappability_df, dna_cell_ids, source=cnv_mappabilitynorm_path)

    outdir.mkdir(parents=True, exist_ok=True)
    cells_out = outdir / "cnv_cells.tsv"
    bins_out = outdir / "cnv_bins.tsv"
    obs_df.reset_index(drop=True).to_csv(cells_out, sep="\t", index=False)
    bins_df.to_csv(bins_out, sep="\t", index=False)

    h5ad_out = outdir / "cnv.h5ad"
    h5ad_written = False
    if write_h5ad:
        if not HAS_ANNDATA:
            raise RuntimeError("anndata is required for import-scrr-cnv h5ad output")
        var_df = bins_df.set_index("bin_id", drop=False)
        adata = ad.AnnData(
            X=copy_number,
            obs=sanitize_frame_for_anndata(obs_df),
            var=sanitize_frame_for_anndata(var_df),
        )
        if mappabilitynorm is not None:
            adata.layers["mappabilitynorm"] = mappabilitynorm
            layer_paths["mappabilitynorm"] = str(cnv_mappabilitynorm_path.resolve())
        adata.uns["circyto"] = _safe_uns(
            {
                "command_name": "circyto import-scrr-cnv",
                "circyto_version": __version__,
                "source_cnv_states": str(cnv_states_path.resolve()),
                "source_cnv_mappabilitynorm": str(cnv_mappabilitynorm_path.resolve())
                if cnv_mappabilitynorm_path is not None
                else "",
                "cell_mapping": str(cell_mapping_path.resolve()) if cell_mapping_path is not None else "",
                "obs_id_strategy": obs_id_strategy,
                "dna_modality_interpretation": "scRR-seq DNA CNV/copy-number state by genomic bin",
                "X_semantics": "integer copy number parsed from '<N>-somy' CNV states",
                "layers": layer_paths,
            }
        )
        adata.write_h5ad(str(h5ad_out))
        h5ad_written = True

    unique_bin_sizes = sorted({int(value) for value in bins_df["bin_size"].astype(int).tolist()})
    states, counts = np.unique(copy_number, return_counts=True)
    state_counts = {str(int(state)): int(count) for state, count in zip(states, counts)}
    shared_prefix_pairs = int(sum(1 for dna in dna_cell_ids if infer_paired_rna_cell_id(dna) != dna))

    summary: dict[str, Any] = {
        "description": "Imported processed scRR-seq DNA CNV summaries into a bin-level CNV modality.",
        "source_cnv_states": str(cnv_states_path.resolve()),
        "source_cnv_mappabilitynorm": str(cnv_mappabilitynorm_path.resolve())
        if cnv_mappabilitynorm_path is not None
        else None,
        "cell_mapping": str(cell_mapping_path.resolve()) if cell_mapping_path is not None else None,
        "obs_id_strategy": obs_id_strategy,
        "n_cells": int(len(dna_cell_ids)),
        "n_bins": int(bins_df.shape[0]),
        "bin_sizes": unique_bin_sizes,
        "copy_number_state_counts": state_counts,
        "n_dna_ids_with_inferred_rna_prefix_pair": shared_prefix_pairs,
        "output_cnv_cells": str(cells_out.resolve()),
        "output_cnv_bins": str(bins_out.resolve()),
        "output_cnv_h5ad": str(h5ad_out.resolve()) if h5ad_written else None,
        "modality_name": "cnv",
        "X_semantics": "cell x genomic-bin integer copy number parsed from CNV state labels",
        "layer_semantics": {
            "mappabilitynorm": "cell x genomic-bin numeric mappability-normalized DNA signal"
            if cnv_mappabilitynorm_path is not None
            else None
        },
    }
    summary_path = outdir / "scrr_cnv_import_summary.json"
    summary["output_summary"] = str(summary_path.resolve())
    write_json(summary_path, summary)
    return apply_standard_provenance(
        summary,
        command_name="circyto import-scrr-cnv",
        workflow_type="scrr-cnv-import",
        protocol="scRR-seq",
        read_layout="processed-geo-summary",
        genome_fasta=None,
        gtf=None,
        detector_backend="scRepli-seq-CNV-summary",
        started_at=utc_now_iso(),
        completed_at=utc_now_iso(),
        elapsed_seconds=0.0,
    )
