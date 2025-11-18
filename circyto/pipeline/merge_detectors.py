# circyto/pipeline/merge_detectors.py
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List

import pandas as pd

from circyto.utils import ensure_dir


@dataclass
class DetectorSummary:
    name: str
    outdir: Path
    n_cells: int


def _load_summary(indir: Path) -> Dict[str, DetectorSummary]:
    """
    Load summary.json from a multi-detector run.

    Expected format (from `run-multidetector`):

    {
      "ciri-full": {"n_cells": 2, "detector": "ciri-full", "outdir": "work/multi/ciri-full"},
      "ciri2":     {"n_cells": 2, "detector": "ciri2",     "outdir": "work/multi/ciri2"}
    }
    """
    summary_path = indir / "summary.json"
    if not summary_path.exists():
        raise FileNotFoundError(f"summary.json not found in {indir}")

    with summary_path.open() as f:
        raw = json.load(f)

    out: Dict[str, DetectorSummary] = {}
    for det_name, info in raw.items():
        outdir = Path(info.get("outdir", ""))  # stored as string
        if not outdir.is_absolute():
            outdir = indir / Path(outdir).name  # normalize to subdir
        out[det_name] = DetectorSummary(
            name=det_name,
            outdir=outdir,
            n_cells=int(info.get("n_cells", 0)),
        )
    return out


def _read_detector_tsvs(det: DetectorSummary) -> pd.DataFrame:
    """
    Read all per-cell TSVs for a single detector and aggregate per circ_id.

    Each TSV must have columns:
      circ_id, chr, start, end, strand, support
    """
    tsv_files = sorted(det.outdir.glob("*.tsv"))
    if not tsv_files:
        # Return empty frame with expected columns
        return pd.DataFrame(
            columns=["circ_id", "chr", "start", "end", "strand", "support", "cell_id"]
        )

    dfs: List[pd.DataFrame] = []
    for tsv in tsv_files:
        cell_id = tsv.stem
        df = pd.read_csv(tsv, sep="\t")
        if df.empty:
            continue
        # tolerate extra columns
        required = ["circ_id", "chr", "start", "end", "strand", "support"]
        missing = [c for c in required if c not in df.columns]
        if missing:
            raise ValueError(f"Missing columns {missing} in {tsv}")
        df = df[required].copy()
        df["cell_id"] = cell_id
        dfs.append(df)

    if not dfs:
        return pd.DataFrame(
            columns=["circ_id", "chr", "start", "end", "strand", "support", "cell_id"]
        )

    return pd.concat(dfs, axis=0, ignore_index=True)


def merge_detectors(indir: Path, outdir: Path) -> None:
    """
    Merge outputs from `run-multidetector` into detector-comparison tables.

    Produces:
      - circ_union.tsv
          circ_id, chr, start, end, strand,
          <det>_total_support, <det>_n_cells, <det>_present
      - circ_by_detector.tsv (long format)
          circ_id, chr, start, end, strand, detector,
          total_support, n_cells
      - metadata.json
          detectors, n_circs, per-detector stats, input summary
    """
    indir = Path(indir)
    outdir = Path(outdir)
    ensure_dir(outdir)

    summaries = _load_summary(indir)

    per_detector_frames: Dict[str, pd.DataFrame] = {}
    union_keys = set()

    # 1) Load and aggregate per detector
    det_stats = {}
    for det_name, det_sum in summaries.items():
        df = _read_detector_tsvs(det_sum)
        if df.empty:
            per_detector_frames[det_name] = df
            det_stats[det_name] = {
                "n_cells": det_sum.n_cells,
                "n_rows": 0,
                "n_circs": 0,
            }
            continue

        agg = (
            df.groupby(["circ_id", "chr", "start", "end", "strand"])
            .agg(total_support=("support", "sum"), n_cells=("cell_id", "nunique"))
            .reset_index()
        )

        per_detector_frames[det_name] = agg
        union_keys.update(agg["circ_id"].tolist())

        det_stats[det_name] = {
            "n_cells": det_sum.n_cells,
            "n_rows": int(df.shape[0]),
            "n_circs": int(agg.shape[0]),
        }

    # 2) Build union table with canonical coords
    if not union_keys:
        union_df = pd.DataFrame(
            columns=["circ_id", "chr", "start", "end", "strand"]
        )
    else:
        rows = []
        seen = set()
        for det_name, agg in per_detector_frames.items():
            if agg.empty:
                continue
            for _, r in agg.iterrows():
                cid = r["circ_id"]
                if cid in seen:
                    continue
                seen.add(cid)
                rows.append(
                    {
                        "circ_id": cid,
                        "chr": r["chr"],
                        "start": r["start"],
                        "end": r["end"],
                        "strand": r["strand"],
                    }
                )
        union_df = pd.DataFrame(rows)

    # 3) Add per-detector columns
    for det_name, agg in per_detector_frames.items():
        if union_df.empty:
            # no circRNAs at all; still add columns
            union_df[f"{det_name}_total_support"] = 0
            union_df[f"{det_name}_n_cells"] = 0
            union_df[f"{det_name}_present"] = 0
            continue

        if agg.empty:
            # detector saw nothing; all zeros
            union_df[f"{det_name}_total_support"] = 0
            union_df[f"{det_name}_n_cells"] = 0
            union_df[f"{det_name}_present"] = 0
            continue

        merged = union_df.merge(
            agg[["circ_id", "total_support", "n_cells"]],
            on="circ_id",
            how="left",
        )

        # Fill NaNs with 0 for this detector and rename
        union_df = merged
        union_df[f"{det_name}_total_support"] = (
            union_df["total_support"].fillna(0).astype(int)
        )
        union_df[f"{det_name}_n_cells"] = (
            union_df["n_cells"].fillna(0).astype(int)
        )
        union_df[f"{det_name}_present"] = (
            union_df[f"{det_name}_n_cells"] > 0
        ).astype(int)

        # Drop the generic columns so the next detector can reuse these names
        union_df = union_df.drop(columns=["total_support", "n_cells"])

    union_path = outdir / "circ_union.tsv"
    union_df.to_csv(union_path, sep="\t", index=False)

    # 4) Long-format table
    long_rows = []
    for det_name, agg in per_detector_frames.items():
        if agg.empty:
            continue
        for _, r in agg.iterrows():
            long_rows.append(
                {
                    "circ_id": r["circ_id"],
                    "chr": r["chr"],
                    "start": r["start"],
                    "end": r["end"],
                    "strand": r["strand"],
                    "detector": det_name,
                    "total_support": int(r["total_support"]),
                    "n_cells": int(r["n_cells"]),
                }
            )

    long_df = pd.DataFrame(long_rows)
    long_path = outdir / "circ_by_detector.tsv"
    long_df.to_csv(long_path, sep="\t", index=False)

    # 5) Metadata
    meta = {
        "input_dir": str(indir),
        "output_dir": str(outdir),
        "n_circs_union": int(union_df.shape[0]),
        "detectors": det_stats,
        "raw_summary": {
            k: {"outdir": str(v.outdir), "n_cells": v.n_cells}
            for k, v in summaries.items()
        },
    }

    with (outdir / "metadata.json").open("w") as f:
        json.dump(meta, f, indent=2)
