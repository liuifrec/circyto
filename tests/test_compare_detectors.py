# tests/test_compare_detectors.py
from pathlib import Path

import pandas as pd

from circyto.pipeline.compare_detectors import compare_detectors


def test_merge_detectors_simple(tmp_path: Path):
    # --- 1) Create a fake merged circ_union.tsv with 2 detectors ---
    indir = tmp_path / "multi"
    indir.mkdir()

    union_df = pd.DataFrame(
        [
            # circA detected by both
            {
                "circ_id": "circA",
                "chr": "chr1",
                "start": 100,
                "end": 200,
                "strand": "+",
                "ciri-full_total_support": 5,
                "ciri-full_n_cells": 2,
                "ciri2_total_support": 5,
                "ciri2_n_cells": 1,
            },
            # circB only detected by ciri-full
            {
                "circ_id": "circB",
                "chr": "chr1",
                "start": 300,
                "end": 400,
                "strand": "+",
                "ciri-full_total_support": 1,
                "ciri-full_n_cells": 1,
                "ciri2_total_support": 0,
                "ciri2_n_cells": 0,
            },
        ]
    )
    union_path = indir / "circ_union.tsv"
    union_df.to_csv(union_path, sep="\t", index=False)

    outdir = tmp_path / "compare"
    result = compare_detectors(indir=indir, outdir=outdir)

    # --- 2) Files exist ---
    jaccard_path = outdir / "jaccard.tsv"
    summary_path = outdir / "detector_summary.tsv"
    meta_path = outdir / "compare_metadata.json"

    assert jaccard_path.exists()
    assert summary_path.exists()
    assert meta_path.exists()

    # --- 3) Jaccard matrix correctness ---
    jac = pd.read_csv(jaccard_path, sep="\t", index_col=0)
    assert set(jac.index) == {"ciri-full", "ciri2"}
    assert set(jac.columns) == {"ciri-full", "ciri2"}

    # circA is in both, circB only in ciri-full → intersection=1, union=2 → 0.5
    assert jac.loc["ciri-full", "ciri2"] == 0.5
    assert jac.loc["ciri2", "ciri-full"] == 0.5
    assert jac.loc["ciri-full", "ciri-full"] == 1.0
    assert jac.loc["ciri2", "ciri2"] == 1.0

    # --- 4) Detector summary correctness ---
    summary = pd.read_csv(summary_path, sep="\t", index_col=0)

    # ciri-full: circA + circB → 2 circRNAs, total_support = 5+1 = 6, cells = 2+1 = 3
    sf = summary.loc["ciri-full"]
    assert sf["n_circ"] == 2
    assert sf["total_support"] == 6
    assert sf["total_cells"] == 3

    # ciri2: circA only
    s2 = summary.loc["ciri2"]
    assert s2["n_circ"] == 1
    assert s2["total_support"] == 5
    assert s2["total_cells"] == 1

    # Sanity: detectors list in result matches
    assert result.detectors == ["ciri-full", "ciri2"]
