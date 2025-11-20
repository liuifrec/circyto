# circyto/pipeline/compare_detectors.py
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List

import json
import numpy as np
import pandas as pd


@dataclass
class CompareResult:
    """
    Simple container for compare-detectors outputs.
    """

    indir: Path
    outdir: Path
    detectors: List[str]
    jaccard_path: Path
    summary_path: Path
    metadata_path: Path


def _infer_detectors_from_union(df: pd.DataFrame) -> List[str]:
    """
    Infer detector names from columns like '<detector>_total_support'.
    """
    detectors: List[str] = []
    for col in df.columns:
        if col.endswith("_total_support"):
            det = col[: -len("_total_support")]
            detectors.append(det)
    return sorted(set(detectors))


def _compute_jaccard(df: pd.DataFrame, detectors: List[str]) -> pd.DataFrame:
    """
    Compute Jaccard index on *sets of circ_ids* per detector.

    Presence is defined as `total_support > 0` for that detector.
    """
    circ_id_col = "circ_id"
    if circ_id_col not in df.columns:
        raise ValueError("circ_union.tsv must contain a 'circ_id' column")

    # Build sets of circ_ids for each detector
    det_sets: Dict[str, set] = {}
    for det in detectors:
        col = f"{det}_total_support"
        if col not in df.columns:
            det_sets[det] = set()
            continue

        s = df[col].fillna(0).astype(float)
        ids = set(df.loc[s > 0, circ_id_col].tolist())
        det_sets[det] = ids

    # Jaccard matrix
    jaccard = pd.DataFrame(
        np.zeros((len(detectors), len(detectors)), dtype=float),
        index=detectors,
        columns=detectors,
    )

    for di in detectors:
        Si = det_sets[di]
        for dj in detectors:
            Sj = det_sets[dj]
            if not Si and not Sj:
                val = 0.0
            else:
                inter = len(Si & Sj)
                union = len(Si | Sj)
                val = float(inter) / float(union) if union > 0 else 0.0
            jaccard.loc[di, dj] = val

    return jaccard


def _build_detector_summary(df: pd.DataFrame, detectors: List[str]) -> pd.DataFrame:
    """
    Build a detector-level summary matching the tests' expectations:

      - n_circ: # circ_ids with total_support > 0
      - total_support: sum of total_support across all circ_ids
      - total_cells: sum of n_cells across all circ_ids
    """
    rows: List[Dict[str, object]] = []

    for det in detectors:
        support_col = f"{det}_total_support"
        cells_col = f"{det}_n_cells"

        if support_col in df.columns:
            s = df[support_col].fillna(0).astype(float)
            n_circ = int((s > 0).sum())
            total_support_sum = float(s.sum())
        else:
            n_circ = 0
            total_support_sum = 0.0

        if cells_col in df.columns:
            c = df[cells_col].fillna(0).astype(float)
            n_cells_sum = float(c.sum())
        else:
            n_cells_sum = 0.0

        rows.append(
            {
                "detector": det,
                "n_circ": n_circ,
                "total_support": total_support_sum,
                "total_cells": n_cells_sum,
            }
        )

    return pd.DataFrame(rows).set_index("detector")


def compare_detectors(indir: Path, outdir: Path) -> CompareResult:
    """
    Compare detectors based on a merged circ_union.tsv produced by
    `merge-detectors`.

    Expects:

      indir/
        circ_union.tsv

    Writes:

      outdir/
        jaccard.tsv
        detector_summary.tsv
        compare_metadata.json
    """
    indir = Path(indir)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    union_path = indir / "circ_union.tsv"
    if not union_path.exists():
        raise FileNotFoundError(f"Expected {union_path} (from merge-detectors)")

    df = pd.read_csv(union_path, sep="\t")
    if "circ_id" not in df.columns:
        raise ValueError("circ_union.tsv must contain a 'circ_id' column")

    detectors = _infer_detectors_from_union(df)
    if not detectors:
        raise ValueError("No detectors inferred from circ_union.tsv columns")

    # --- Jaccard ---
    jaccard = _compute_jaccard(df, detectors)
    jaccard_path = outdir / "jaccard.tsv"
    jaccard.to_csv(jaccard_path, sep="\t")

    # --- Summary ---
    summary = _build_detector_summary(df, detectors)
    summary_path = outdir / "detector_summary.tsv"
    summary.to_csv(summary_path, sep="\t")  # keep index=detector

    # --- Metadata ---
    meta = {
        "n_circ": int(df.shape[0]),
        "detectors": detectors,
        "jaccard_path": str(jaccard_path),
        "summary_path": str(summary_path),
    }
    metadata_path = outdir / "compare_metadata.json"
    with metadata_path.open("w") as f:
        json.dump(meta, f, indent=2)

    return CompareResult(
        indir=indir,
        outdir=outdir,
        detectors=detectors,
        jaccard_path=jaccard_path,
        summary_path=summary_path,
        metadata_path=metadata_path,
    )
