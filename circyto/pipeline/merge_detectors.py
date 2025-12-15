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

    Supported formats:

    1) Old (flat) schema:

        {
          "ciri-full": {"n_cells": 2, "detector": "ciri-full", "outdir": "work/multi/ciri-full"},
          "ciri2":     {"n_cells": 2, "detector": "ciri2",     "outdir": "work/multi/ciri2"}
        }

    2) New schema from `run-multidetector` (nested):

        {
          "summary": {
            "ciri-full":  {"n_cells": 16, "detector": "ciri-full",  "outdir": "work/multi/ciri-full"},
            "find-circ3": {"n_cells": 16, "detector": "find-circ3", "outdir": "work/multi/find-circ3"}
          },
          "results": {
            "ciri-full":  [...],
            "find-circ3": [...]
          }
        }
    """
    summary_path = indir / "summary.json"
    if not summary_path.exists():
        raise FileNotFoundError(f"summary.json not found in {indir}")

    with summary_path.open() as f:
        raw = json.load(f)

    # Detect new vs old schema
    if isinstance(raw, dict) and "summary" in raw and isinstance(raw["summary"], dict):
        summary_block = raw["summary"]
    else:
        summary_block = raw

    out: Dict[str, DetectorSummary] = {}
    for det_name, info in summary_block.items():
        if not isinstance(info, dict):
            continue
        outdir = Path(info.get("outdir", ""))

        # Normalize to subdir under indir if relative
        if not outdir.is_absolute():
            # Old runs sometimes stored a full path; if it's just a name, join with indir
            outdir = indir / Path(outdir).name

        out[det_name] = DetectorSummary(
            name=det_name,
            outdir=outdir,
            n_cells=int(info.get("n_cells", 0)),
        )
    return out


def _parse_circ_id(circ_id: str) -> tuple[str, int, int, str]:
    """
    Parse circ_id of form: 'chr:start|end|strand'

    Returns (chr, start, end, strand).

    This matches the IDs produced by CIRI-full and our find_circ3 collector.
    """
    chrom, rest = circ_id.split(":", 1)
    start_str, rest2 = rest.split("|", 1)
    end_str, strand = rest2.split("|", 1)
    return chrom, int(start_str), int(end_str), strand


def _read_detector_matrix_from_triplet(det: DetectorSummary) -> pd.DataFrame:
    """
    Fallback reader: reconstruct (circ, cell, support) from matrix triplet.

    Expects files in:

        det.outdir.parent / "matrices" / f"{det.name}.mtx"
        det.outdir.parent / "matrices" / f"{det.name}.circ.txt"
        det.outdir.parent / "matrices" / f"{det.name}.cell.txt"

    Returns DataFrame with columns:
        circ_id, chr, start, end, strand, support, cell_id
    """
    base = det.outdir.parent  # e.g. work/multidetector_chr21_16cells
    matrix_dir = base / "matrices"

    mtx_path = matrix_dir / f"{det.name}.mtx"
    circ_index_path = matrix_dir / f"{det.name}.circ.txt"
    cell_index_path = matrix_dir / f"{det.name}.cell.txt"

    if not (mtx_path.exists() and circ_index_path.exists() and cell_index_path.exists()):
        return pd.DataFrame(
            columns=["circ_id", "chr", "start", "end", "strand", "support", "cell_id"]
        )

    # Load circ and cell indices
    circ_ids = [line.strip() for line in circ_index_path.read_text().splitlines() if line.strip()]
    cell_ids = [line.strip() for line in cell_index_path.read_text().splitlines() if line.strip()]

    if not circ_ids or not cell_ids:
        return pd.DataFrame(
            columns=["circ_id", "chr", "start", "end", "strand", "support", "cell_id"]
        )

    rows: List[dict] = []
    header_seen = False

    with mtx_path.open() as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("%"):
                continue

            parts = line.split()

            if not header_seen:
                # First non-comment line is the header: n_rows n_cols nnz
                header_seen = True
                # We could sanity-check here, but we don't strictly need it
                continue

            if len(parts) != 3:
                continue

            i_str, j_str, v_str = parts
            val = int(v_str)
            if val == 0:
                continue

            r_idx = int(i_str) - 1  # MatrixMarket is 1-based
            c_idx = int(j_str) - 1

            if r_idx < 0 or r_idx >= len(circ_ids):
                continue
            if c_idx < 0 or c_idx >= len(cell_ids):
                continue

            circ_id = circ_ids[r_idx]
            cell_id = cell_ids[c_idx]

            chrom, start, end, strand = _parse_circ_id(circ_id)

            rows.append(
                dict(
                    circ_id=circ_id,
                    chr=chrom,
                    start=start,
                    end=end,
                    strand=strand,
                    support=val,
                    cell_id=cell_id,
                )
            )

    if not rows:
        return pd.DataFrame(
            columns=["circ_id", "chr", "start", "end", "strand", "support", "cell_id"]
        )

    return pd.DataFrame(rows)


def _read_detector_tsvs_or_matrix(det: DetectorSummary) -> pd.DataFrame:
    """
    Read per-cell TSVs for a single detector, or fall back to matrix triplet.

    TSV layout (per-cell):
        circ_id, chr, start, end, strand, support

    Matrix triplet layout (fallback):
        <matrices>/<det>.mtx + <det>.circ.txt + <det>.cell.txt
    """
    # First try per-cell TSVs
    tsv_files = sorted(det.outdir.glob("*.tsv"))
    if tsv_files:
        dfs: List[pd.DataFrame] = []
        for tsv in tsv_files:
            cell_id = tsv.stem
            df = pd.read_csv(tsv, sep="\t")
            if df.empty:
                continue

            required = ["circ_id", "chr", "start", "end", "strand", "support"]
            missing = [c for c in required if c not in df.columns]
            if missing:
                raise ValueError(f"Missing columns {missing} in {tsv}")

            df = df[required].copy()
            df["cell_id"] = cell_id
            dfs.append(df)

        if dfs:
            return pd.concat(dfs, axis=0, ignore_index=True)

    # Fallback: matrix triplet (for detectors like find_circ3)
    return _read_detector_matrix_from_triplet(det)


def merge_detectors(indir: Path, outdir: Path) -> None:
    """
    Merge outputs from `run-multidetector` into detector-comparison tables.

    Produces:
      - circ_union.tsv
          circ_id, chr, start, end, strand,
          <det>_total_support, <det>_n_cells, <det>_present

      - circ_by_detector.tsv (long format)
          circ_id, chr, start, end, strand,
          detector, total_support, n_cells
    """
    ensure_dir(outdir)

    summaries = _load_summary(indir)
    if not summaries:
        raise ValueError(f"No detector summaries found in {indir}")

    union_keys = set()
    per_detector_frames: Dict[str, pd.DataFrame] = {}
    det_stats: Dict[str, dict] = {}

    # 1) Aggregate each detector separately
    for det_name, det_sum in summaries.items():
        df = _read_detector_tsvs_or_matrix(det_sum)
        if df.empty:
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

    # If absolutely no circRNAs anywhere, write clean empties
    union_path = outdir / "circ_union.tsv"
    long_path = outdir / "circ_by_detector.tsv"

    if not union_keys:
        pd.DataFrame(
            columns=["circ_id", "chr", "start", "end", "strand"]
        ).to_csv(union_path, sep="\t", index=False)

        pd.DataFrame(
            columns=[
                "circ_id",
                "chr",
                "start",
                "end",
                "strand",
                "detector",
                "total_support",
                "n_cells",
            ]
        ).to_csv(long_path, sep="\t", index=False)

        meta = {
            "input_dir": str(indir),
            "output_dir": str(outdir),
            "n_circs_union": 0,
            "detectors": det_stats,
            "raw_summary": {
                k: {"outdir": str(v.outdir), "n_cells": v.n_cells}
                for k, v in summaries.items()
            },
        }
        with (outdir / "metadata.json").open("w") as f:
            json.dump(meta, f, indent=2)
        return

    # 2) Build union frame over all circ_ids
    union_keys = sorted(union_keys)
    union_df = pd.DataFrame({"circ_id": union_keys})
    base_cols = ["chr", "start", "end", "strand"]

    for det_name, agg in per_detector_frames.items():
        # On the first detector, attach coordinates
        if not all(c in union_df.columns for c in base_cols):
            union_df = union_df.merge(
                agg[["circ_id"] + base_cols].drop_duplicates("circ_id"),
                on="circ_id",
                how="left",
            )

        merged = union_df.merge(
            agg[["circ_id", "total_support", "n_cells"]],
            on="circ_id",
            how="left",
        )

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

        # Drop temporary columns before next detector
        union_df = union_df.drop(columns=["total_support", "n_cells"])

    # 3) Long-format frame
    long_rows: List[pd.DataFrame] = []
    for det_name, agg in per_detector_frames.items():
        tmp = agg.copy()
        tmp["detector"] = det_name
        long_rows.append(tmp)

    if long_rows:
        long_df = pd.concat(long_rows, axis=0, ignore_index=True)
    else:
        long_df = pd.DataFrame(
            columns=[
                "circ_id",
                "chr",
                "start",
                "end",
                "strand",
                "detector",
                "total_support",
                "n_cells",
            ]
        )

    # 4) Write outputs
    union_df.to_csv(union_path, sep="\t", index=False)
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
