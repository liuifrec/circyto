"""Strict detector-free import; processing status never changes count values."""

from __future__ import annotations

import csv
import hashlib
import html
import json
import re
import shutil
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread, mmwrite

from circyto import __version__
from circyto.multimodal.sync import (
    MUDATA_SYNC_POLICY,
    mudata_from_modalities,
    write_h5mu,
)
from circyto.pipeline.annotate_host_gene import (
    HOST_GENE_COLUMNS,
    normalize_host_gene_annotations,
)
from circyto.pipeline.qc_report_style import STYLE

SCHEMA = "circyto.import_circ.v1"
STATUSES = {"completed", "failed", "skipped", "unprocessed", "unknown"}
COORDINATES = {"0-based-half-open", "1-based-closed"}
FEATURES = ["circ_id", "chrom", "start", "end", "strand"]
OBS_FIELDS = {
    "circ_observation_status",
    "circRNA_count",
    "circRNA_total_support",
    "membership",
    "circ_present",
    "rna_present",
}
NOTE = (
    "Completion is user-supplied processing evidence, not locus-level callability. "
    "Zero support does not establish a true negative or biological absence. "
    "Failed, skipped, unprocessed and unknown cells are excluded from detection-rate denominators; "
    "their supplied counts remain unchanged."
)


def _table(path: Path, required: list[str]) -> pd.DataFrame:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = csv.reader(handle, delimiter="\t")
        header = next(rows, [])
        if not header or len(set(header)) != len(header) or any(not x for x in header):
            raise ValueError(f"{path.name}: empty or duplicate column names")
        if not set(required) <= set(header):
            raise ValueError(f"{path.name}: required columns: {', '.join(required)}")
        for number, row in enumerate(rows, 2):
            if len(row) != len(header):
                raise ValueError(f"{path.name}: malformed TSV row {number}")
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def _ids(values, label: str, *, unique: bool = True) -> list[str]:
    ids = list(values)
    if any(
        not isinstance(x, str)
        or not x
        or x != x.strip()
        or any(c in x for c in "\r\n\t")
        for x in ids
    ):
        raise ValueError(
            f"{label}: IDs must be nonempty strings without surrounding whitespace"
        )
    if unique and len(set(ids)) != len(ids):
        raise ValueError(f"{label}: duplicate IDs")
    return ids


def _integers(values, label: str) -> np.ndarray:
    # Never pass integers through float: preserve support above 2**53 exactly.
    result = []
    for value in values:
        if not re.fullmatch(r"[0-9]+", str(value)):
            raise ValueError(f"{label}: expected nonnegative integers, found {value!r}")
        number = int(value)
        if number > np.iinfo(np.int64).max:
            raise ValueError(f"{label}: int64 overflow")
        result.append(number)
    return np.array(result, dtype=np.int64)


def _declarations(frame: pd.DataFrame, declarations: dict[str, str]) -> None:
    for column, expected in declarations.items():
        if column in frame and not frame[column].eq(expected).all():
            raise ValueError(
                f"Conflicting or missing {column}; expected only {expected!r}"
            )


def _features(frame: pd.DataFrame, declarations: dict[str, str]) -> pd.DataFrame:
    _ids(frame.circ_id, "candidate")
    _ids(frame.chrom, "chromosome", unique=False)
    frame = frame.copy()
    for col in ("start", "end"):
        frame[col] = _integers(frame[col], col)
    minimum = 0 if declarations["coordinate_convention"] == "0-based-half-open" else 1
    invalid = (frame.start < minimum) | (frame.end < frame.start)
    if minimum == 0:
        invalid |= frame.end == frame.start
    if invalid.any() or not frame.strand.isin(["+", "-", "."]).all():
        raise ValueError(
            "Invalid interval or strand for declared coordinate convention"
        )
    _declarations(frame, declarations)
    for key, value in declarations.items():
        frame[key] = value
    frame = frame.set_index("circ_id", drop=False)
    # Reuse schema normalization, but never overwrite supplied annotation/provenance.
    normalized = normalize_host_gene_annotations(frame)
    for col in HOST_GENE_COLUMNS:
        if col in frame:
            normalized[col] = frame[col]
    normalized.index.name = None
    return normalized


def _load_calls(path: Path, cells: list[str], declarations: dict[str, str]):
    calls = _table(path, ["cell_id", *FEATURES, "support"])
    _ids(calls.cell_id, "call cell", unique=False)
    _ids(calls.circ_id, "call candidate", unique=False)
    _declarations(calls, declarations)
    if not set(calls.cell_id) <= set(cells):
        raise ValueError("Calls contain cell IDs absent from cell metadata")
    if calls.duplicated(["cell_id", "circ_id"]).any():
        raise ValueError(
            "Duplicate cell/candidate pair; aggregate explicitly before import"
        )
    counts = _integers(calls.support, "support")
    # Extra evidence can vary by cell and remains in the byte-for-byte source copy.
    columns = FEATURES + [
        x for x in HOST_GENE_COLUMNS if x in calls and x not in FEATURES
    ]
    if calls.groupby("circ_id", sort=False)[columns[1:]].nunique().gt(1).any().any():
        raise ValueError("Conflicting candidate coordinates or host-gene annotations")
    features = _features(calls[columns].drop_duplicates("circ_id"), declarations)
    row = pd.Index(cells).get_indexer(calls.cell_id)
    col = features.index.get_indexer(calls.circ_id)
    X = sparse.csr_matrix(
        (counts, (row, col)), shape=(len(cells), len(features)), dtype=np.int64
    )
    return X, features, {"calls.tsv": path}


def _load_matrix(directory: Path, cells: list[str], declarations: dict[str, str]):
    paths = {
        name: directory / name
        for name in (
            "circ_counts.mtx",
            "cell_index.txt",
            "circ_index.txt",
            "circ_feature_table.tsv",
        )
    }
    matrix_cells = _ids(paths["cell_index.txt"].read_text().splitlines(), "matrix cell")
    candidates = _ids(
        paths["circ_index.txt"].read_text().splitlines(), "matrix candidate"
    )
    if not set(matrix_cells) <= set(cells):
        raise ValueError("Matrix contains cell IDs absent from cell metadata")
    features = _table(paths["circ_feature_table.tsv"], FEATURES)
    features = _features(features, declarations)
    if set(features.index) != set(candidates):
        raise ValueError(
            "Candidate index and feature table must contain exactly the same IDs"
        )
    features = features.loc[candidates]
    raw = mmread(paths["circ_counts.mtx"])
    if not sparse.issparse(raw) or raw.shape != (len(candidates), len(matrix_cells)):
        raise ValueError(
            "Expected sparse circ-by-cell Matrix Market matching the index files"
        )
    coo = raw.tocoo()
    if coo.dtype.kind not in "iu" or (coo.data < 0).any():
        raise ValueError(
            "Matrix support must use nonnegative integer Matrix Market values"
        )
    data = _integers(coo.data, "support")
    pairs = pd.MultiIndex.from_arrays([coo.row, coo.col])
    if pairs.has_duplicates:
        raise ValueError("Duplicate matrix entries; aggregate explicitly before import")
    # A completed zero cell must be represented in the source matrix, not inferred
    # from an absent column. Long tables use the explicit metadata cohort instead.
    missing = set(cells) - set(matrix_cells)
    row = pd.Index(cells).get_indexer(matrix_cells)[coo.col]
    X = sparse.csr_matrix((data, (row, coo.row)), shape=(len(cells), len(candidates)))
    return X, features, paths, missing


def _report(metrics: dict, obs: pd.DataFrame, var: pd.DataFrame) -> str:
    def table(frame):
        return frame.head(100).to_html(index=False, escape=True, border=0)

    rate = metrics["detection_fraction_completed"]
    value = (
        "unavailable"
        if rate["value"] is None
        else f"{rate['numerator']} / {rate['denominator']}"
    )
    declaration = html.escape(json.dumps(metrics["declarations"], sort_keys=True))
    return f"""<!doctype html><html lang="en"><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1"><title>circyto import QC</title>
<style>{STYLE}</style><body><header><div class="brand">CIRCYTO</div>Detector-free import QC</header>
<main><section><h1>Imported candidate support</h1><p>{html.escape(NOTE)}</p>
<p>Import succeeded; this is not evidence that a detector workflow completed.</p>
<p>{metrics["n_cells"]} metadata cells; {metrics["n_candidates"]} candidates; {metrics["total_support"]} total support.</p>
<p>Cells with candidates / explicitly completed cells: <strong>{value}</strong>.</p>
<p>Genome and support declarations (user supplied): <code>{declaration}</code></p>
<p>RNA overlap: {html.escape(json.dumps(metrics["rna_overlap"], sort_keys=True))}</p></section>
<section><h2>Observation states</h2>{table(pd.DataFrame(list(metrics["observation_status_counts"].items()), columns=["state", "cells"]))}</section>
<section><h2>Cells</h2><p>First 100 rows; full values in cell_qc.tsv.</p><div class="table-scroll">{table(obs.reset_index())}</div></section>
<section><h2>Candidates</h2><p>First 100 rows; completed-cell prevalence uses the same denominator above. Full values in circ_qc.tsv.</p><div class="table-scroll">{table(var.reset_index(drop=True))}</div></section>
<section><h2>Provenance</h2><p>Original supplied evidence is copied under evidence/. File hashes and exact rates are in metrics.json.</p></section></main></body></html>"""


def import_circ(
    *,
    cells: Path,
    outdir: Path,
    genome_build: str,
    coordinate_convention: str,
    support_definition: str,
    detector: str,
    calls: Path | None = None,
    matrix_dir: Path | None = None,
    rna: Path | None = None,
) -> dict:
    """Import one homogeneous support assay, preserving counts and declared states."""
    if (calls is None) == (matrix_dir is None):
        raise ValueError("Supply exactly one of --calls or --matrix-dir")
    if coordinate_convention not in COORDINATES:
        raise ValueError(f"coordinate-convention must be one of {sorted(COORDINATES)}")
    declarations = {
        "genome_build": genome_build,
        "coordinate_convention": coordinate_convention,
        "support_definition": support_definition,
        "detector": detector,
    }
    if any(not v.strip() or v != v.strip() for v in declarations.values()):
        raise ValueError(
            "Genome, coordinate, support and detector declarations must be nonempty"
        )
    if outdir.exists() or outdir.is_symlink():
        raise ValueError("Output directory already exists; choose a new directory")
    input_paths = [cells]
    if calls is not None:
        input_paths.append(calls)
    else:
        input_paths.extend(
            matrix_dir / name
            for name in (
                "circ_counts.mtx",
                "cell_index.txt",
                "circ_index.txt",
                "circ_feature_table.tsv",
            )
        )
    if rna is not None:
        input_paths.append(rna)
    before_hashes = {path: _sha256(path) for path in input_paths}
    obs = _table(cells, ["cell_id"])
    cell_ids = _ids(obs.cell_id, "metadata cell")
    if not cell_ids:
        raise ValueError("Cell metadata must contain at least one cell")
    if OBS_FIELDS & set(obs.columns):
        raise ValueError(
            f"Cell metadata uses reserved derived fields: {sorted(OBS_FIELDS & set(obs.columns))}"
        )
    _declarations(obs, declarations)
    if "detector_status" not in obs:
        obs["detector_status"] = "unknown"
    if not obs.detector_status.isin(STATUSES).all():
        raise ValueError(f"detector_status must be one of {sorted(STATUSES)}")
    obs = obs.set_index("cell_id")
    if calls is not None:
        X, var, sources = _load_calls(calls, cell_ids, declarations)
    else:
        X, var, sources, missing = _load_matrix(matrix_dir, cell_ids, declarations)
        if any(obs.loc[cell, "detector_status"] == "completed" for cell in missing):
            raise ValueError(
                "Completed metadata cells missing from matrix; supply explicit zero columns"
            )
    X.eliminate_zeros()
    # Sparse sums stay exact and cannot wrap int64.
    if sum(map(int, X.data)) > np.iinfo(np.int64).max:
        raise ValueError("Total support exceeds int64; cannot compute safe summaries")
    obs["circRNA_count"] = X.getnnz(axis=1)
    obs["circRNA_total_support"] = np.asarray(X.sum(axis=1)).ravel()
    completed = obs.detector_status.eq("completed").to_numpy()
    detected = obs.circRNA_count.to_numpy() > 0
    obs["circ_observation_status"] = np.where(
        completed,
        np.where(detected, "completed_detected", "completed_zero"),
        obs.detector_status,
    )
    denominator = int(completed.sum())
    numerator = int((completed & detected).sum())
    var["n_cells_detected"] = X.getnnz(axis=0)
    var["total_support"] = np.asarray(X.sum(axis=0)).ravel()
    var["n_completed_cells_detected"] = X[completed].getnnz(axis=0)
    var["completed_cell_denominator"] = denominator
    var["prevalence_completed"] = (
        var.n_completed_cells_detected / denominator if denominator else np.nan
    )
    sources["cells.tsv"] = cells
    rna_adata = None
    overlap = {}
    if rna is not None:
        try:
            import mudata  # noqa: F401
        except ImportError as exc:
            raise ValueError("Optional RNA export requires circyto[mudata]") from exc
        rna_adata = ad.read_h5ad(rna)
        rna_ids = _ids(rna_adata.obs_names, "RNA cell")
        _ids(rna_adata.var_names, "RNA feature")
        if not set(cell_ids) & set(rna_ids):
            raise ValueError(
                "RNA and circ cell IDs have no overlap; provide a shared cell namespace"
            )
        # Reject build conflicts wherever the usual schema records them; absence is
        # explicitly reported as unverified, never inferred from cell overlap.
        builds = []
        for frame in (rna_adata.obs, rna_adata.var):
            if "genome_build" in frame:
                builds.extend(
                    str(x) for x in frame.genome_build.dropna().unique() if str(x)
                )
        for metadata in (rna_adata.uns, rna_adata.uns.get("circyto", {})):
            if "genome_build" in metadata:
                builds.append(str(metadata["genome_build"]))
        if any(value != genome_build for value in builds):
            raise ValueError("RNA genome_build conflicts with circ genome_build")
        overlap = {
            "n_shared": len(set(cell_ids) & set(rna_ids)),
            "n_circ_only": len(set(cell_ids) - set(rna_ids)),
            "n_rna_only": len(set(rna_ids) - set(cell_ids)),
            "rna_genome_build_verified": bool(builds),
        }
        sources["rna.h5ad"] = rna
    hashes = {name: before_hashes[path] for name, path in sources.items()}
    metrics = {
        "schema_version": SCHEMA,
        "command_name": "circyto import-circ",
        "circyto_version": __version__,
        "declarations": declarations,
        "source_sha256": hashes,
        "n_cells": X.shape[0],
        "n_candidates": X.shape[1],
        "nnz": X.nnz,
        "total_support": int(X.sum()),
        "rna_overlap": overlap,
        "observation_status_counts": {
            str(k): int(v)
            for k, v in obs.circ_observation_status.value_counts().items()
        },
        "detection_fraction_completed": {
            "numerator": numerator,
            "denominator": denominator,
            "value": numerator / denominator if denominator else None,
        },
        "zero_fraction_completed": {
            "numerator": denominator - numerator,
            "denominator": denominator,
            "value": (denominator - numerator) / denominator if denominator else None,
        },
        "interpretation": NOTE,
    }
    provenance = {
        "schema_version": SCHEMA,
        "command_name": "circyto import-circ",
        "circyto_version": __version__,
        **declarations,
        "source_sha256": hashes,
        "metrics_json": json.dumps(metrics, sort_keys=True),
        "observation_status_definition": NOTE,
    }
    circ = ad.AnnData(X=X, obs=obs, var=var, uns={"circyto": provenance})
    # Stage a complete bundle before publishing it; no partial import on validation failure.
    outdir.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix=".circyto-import-", dir=outdir.parent
    ) as tmp:
        stage = Path(tmp) / "bundle"
        for name in ("evidence", "matrix", "anndata", "qc"):
            (stage / name).mkdir(parents=True)
        for name, source in sources.items():
            # RNA can be large: provenance hashes it, the H5MU preserves the object.
            if name != "rna.h5ad":
                shutil.copyfile(source, stage / "evidence" / name)
                if _sha256(stage / "evidence" / name) != hashes[name]:
                    raise ValueError("Source changed during import")
        circ.write_h5ad(stage / "anndata" / "circ_counts.h5ad")
        mmwrite(stage / "matrix" / "circ_counts.mtx", X.T)
        (stage / "matrix" / "cell_index.txt").write_text("\n".join(cell_ids) + "\n")
        (stage / "matrix" / "circ_index.txt").write_text(
            "\n".join(var.index) + ("\n" if len(var) else "")
        )
        var.to_csv(stage / "matrix" / "circ_feature_table.tsv", sep="\t", index=False)
        obs.to_csv(stage / "qc" / "cell_qc.tsv", sep="\t")
        var.to_csv(stage / "qc" / "circ_qc.tsv", sep="\t", index=False)
        (stage / "qc" / "metrics.json").write_text(
            json.dumps(metrics, indent=2, sort_keys=True) + "\n"
        )
        (stage / "qc" / "report.html").write_text(_report(metrics, obs, var))
        if rna_adata is not None:
            bundle = mudata_from_modalities({"circ": circ, "rna": rna_adata})
            global_obs = obs.reindex(bundle.obs_names).copy()
            for col in ("detector_status", "circ_observation_status"):
                global_obs[col] = global_obs[col].fillna("unknown")
            global_obs["circ_present"] = global_obs.index.isin(circ.obs_names)
            global_obs["rna_present"] = global_obs.index.isin(rna_adata.obs_names)
            global_obs["membership"] = np.where(
                global_obs.circ_present & global_obs.rna_present,
                "shared",
                np.where(global_obs.circ_present, "circ_only", "rna_only"),
            )
            bundle.obs = global_obs
            bundle.uns["circyto"] = {
                **provenance,
                "mudata_sync_policy": MUDATA_SYNC_POLICY,
            }
            (stage / "mudata").mkdir()
            write_h5mu(bundle, stage / "mudata" / "circyto_multimodal.h5mu")
        if any(_sha256(path) != before_hashes[path] for path in input_paths):
            raise ValueError("Source changed during import")
        stage.rename(outdir)
    return metrics


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()
