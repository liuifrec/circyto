# tests/test_collect_circexplorer2_matrix.py

from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy.io import mmread

from circyto.pipeline.collect_circexplorer2_matrix import (
    collect_circexplorer2_matrix,
)


def _write_dummy_cell(indir: Path, cell_id: str, rows: list[str]) -> None:
    """
    Helper: create <indir>/<cell_id>/<cell_id>_CIRCexplorer2_circ.txt
    with a simple header and provided data rows.
    """
    cell_dir = indir / cell_id
    cell_dir.mkdir(parents=True, exist_ok=True)

    path = cell_dir / f"{cell_id}_CIRCexplorer2_circ.txt"
    with path.open("w") as f:
        f.write("chrom\tstart\tend\tstrand\n")
        for row in rows:
            f.write(row + "\n")


def test_collect_circexplorer2_matrix_basic(tmp_path: Path) -> None:
    """
    Basic collector test:

      - Two cells, overlapping and unique circRNAs.
      - Check that outputs exist and matrix shape/indices look reasonable.
    """
    indir = tmp_path / "circexplorer2_percell"
    indir.mkdir(parents=True, exist_ok=True)

    # CELL_A has two circRNAs
    _write_dummy_cell(
        indir,
        "CELL_A",
        [
            "chr21\t100\t200\t+",
            "chr21\t300\t400\t-",
        ],
    )

    # CELL_B shares one circRNA and has one unique
    _write_dummy_cell(
        indir,
        "CELL_B",
        [
            "chr21\t100\t200\t+",
            "chr21\t500\t600\t+",
        ],
    )

    outdir = tmp_path / "collected"
    matrix_path, circ_index_path, cell_index_path = collect_circexplorer2_matrix(
        indir=indir,
        outdir=outdir,
        min_counts=1,
    )

    # Files exist
    assert matrix_path.exists()
    assert circ_index_path.exists()
    assert cell_index_path.exists()

    # Read matrix
    mat = mmread(matrix_path).tocsr()
    # 3 unique circRNAs: (100-200), (300-400), (500-600)
    assert mat.shape == (3, 2)

    # total non-zeros == 4 events
    assert mat.nnz == 4

    # Check circ index
    circ_ids = [line.strip() for line in circ_index_path.read_text().splitlines() if line.strip()]
    assert len(circ_ids) == 3
    assert "chr21:100|200|+" in circ_ids
    assert "chr21:300|400|-" in circ_ids
    assert "chr21:500|600|+" in circ_ids

    # Check cell index
    cell_ids = [line.strip() for line in cell_index_path.read_text().splitlines() if line.strip()]
    assert cell_ids == ["CELL_A", "CELL_B"]

    # Sanity: at least one shared circ has non-zero in both columns
    row_idx = circ_ids.index("chr21:100|200|+")
    row = mat.getrow(row_idx).toarray().ravel()
    # one event in each cell
    assert np.all(row == np.array([1, 1]))
