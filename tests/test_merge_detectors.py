# tests/test_merge_detectors.py
from pathlib import Path

import pandas as pd

from circyto.pipeline.merge_detectors import merge_detectors


def _write_detector_tsv(outdir: Path, detector: str, cell: str, rows):
    det_dir = outdir / detector
    det_dir.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(
        rows, columns=["circ_id", "chr", "start", "end", "strand", "support"]
    )
    df.to_csv(det_dir / f"{cell}.tsv", sep="\t", index=False)


def test_merge_detectors_simple(tmp_path: Path):
    # layout: indir/summary.json + detector subdirs
    indir = tmp_path / "multi"
    indir.mkdir()

    # write summary
    summary = {
        "ciri-full": {
            "detector": "ciri-full",
            "n_cells": 2,
            "outdir": str(indir / "ciri-full"),
        },
        "ciri2": {
            "detector": "ciri2",
            "n_cells": 2,
            "outdir": str(indir / "ciri2"),
        },
    }
    import json

    with (indir / "summary.json").open("w") as f:
        json.dump(summary, f)

    # detector 1: two cells, same circ + one unique
    _write_detector_tsv(
        indir,
        "ciri-full",
        "cell1",
        [
            ["circA", "chr1", 100, 200, "+", 3],
            ["circB", "chr1", 300, 400, "+", 1],
        ],
    )
    _write_detector_tsv(
        indir,
        "ciri-full",
        "cell2",
        [
            ["circA", "chr1", 100, 200, "+", 2],
        ],
    )

    # detector 2: one cell, circA only
    _write_detector_tsv(
        indir,
        "ciri2",
        "cell1",
        [
            ["circA", "chr1", 100, 200, "+", 5],
        ],
    )

    outdir = tmp_path / "merged"
    merge_detectors(indir=indir, outdir=outdir)

    union_path = outdir / "circ_union.tsv"
    long_path = outdir / "circ_by_detector.tsv"
    meta_path = outdir / "metadata.json"

    assert union_path.exists()
    assert long_path.exists()
    assert meta_path.exists()

    union = pd.read_csv(union_path, sep="\t")
    # circA + circB
    assert set(union["circ_id"]) == {"circA", "circB"}

    # circA aggregated
    circA = union[union["circ_id"] == "circA"].iloc[0]
    assert circA["ciri-full_total_support"] == 5  # 3 + 2
    assert circA["ciri-full_n_cells"] == 2
    assert circA["ciri2_total_support"] == 5
    assert circA["ciri2_n_cells"] == 1

    # circB only in ciri-full
    circB = union[union["circ_id"] == "circB"].iloc[0]
    assert circB["ciri-full_total_support"] == 1
    assert circB["ciri2_total_support"] == 0

    long_df = pd.read_csv(long_path, sep="\t")
    # 3 rows: circA/ciri-full, circB/ciri-full, circA/ciri2
    assert long_df.shape[0] == 3
