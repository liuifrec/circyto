# circyto/pipeline/collect_circexplorer2_matrix.py

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
from scipy import sparse
from scipy.io import mmwrite


def _parse_circexplorer2_file(path: Path) -> List[Tuple[str, int]]:
    """
    Parse a single CIRCexplorer2 per-cell output file.

    Assumes a tab-separated file with a header line containing at least:
      - chrom
      - start
      - end
      - strand

    Each row is treated as one circRNA event with count = 1.

    Returns
    -------
    events : list of (circ_id, count)
        circ_id = "chrom:start|end|strand"
        count  = 1 (can be extended to use a dedicated count column).
    """
    events: List[Tuple[str, int]] = []
    with path.open("r") as f:
        header = None
        for line in f:
            line = line.strip()
            if not line:
                continue
            # detect header
            if header is None:
                header = line.split("\t")
                # normalize names
                header = [h.strip().lower() for h in header]
                # basic sanity check
                required = {"chrom", "start", "end", "strand"}
                missing = required - set(header)
                if missing:
                    raise ValueError(
                        f"[collect_circexplorer2] {path} missing columns: {missing}. "
                        f"Found columns: {header}"
                    )
                chrom_idx = header.index("chrom")
                start_idx = header.index("start")
                end_idx = header.index("end")
                strand_idx = header.index("strand")
                continue

            # data line
            parts = line.split("\t")
            if len(parts) < 4:
                # skip malformed rows
                continue

            chrom = parts[chrom_idx]
            start = parts[start_idx]
            end = parts[end_idx]
            strand = parts[strand_idx]

            circ_id = f"{chrom}:{start}|{end}|{strand}"
            # For now: each row is one event (count = 1)
            events.append((circ_id, 1))

    return events


def collect_circexplorer2_matrix(
    indir: str | Path,
    outdir: str | Path,
    min_counts: int = 1,
) -> tuple[Path, Path, Path]:
    """
    Collect CIRCexplorer2 per-cell circRNA calls into a MatrixMarket matrix.

    Parameters
    ----------
    indir : str or Path
        Directory containing per-cell subdirectories. Each subdirectory
        is expected to be named after the cell_id and contain a file
        <cell_id>_CIRCexplorer2_circ.txt.
    outdir : str or Path
        Output directory. Will be created if needed.
    min_counts : int, default 1
        Minimum total count across all cells required to keep a circRNA.

    Returns
    -------
    matrix_path, circ_index_path, cell_index_path : tuple[Path, Path, Path]
        Paths to the generated MatrixMarket and index files.

    Output files
    ------------
    - circ_matrix.mtx        : MatrixMarket (rows=circ, cols=cells)
    - circ_index.tsv         : 1 column, circ_id
    - cell_index.tsv         : 1 column, cell_id
    """
    indir = Path(indir)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Discover cell directories
    cell_dirs = sorted(
        [d for d in indir.iterdir() if d.is_dir()],
        key=lambda p: p.name,
    )
    if not cell_dirs:
        raise ValueError(f"[collect_circexplorer2] No cell directories found in {indir}")

    cell_ids = [d.name for d in cell_dirs]

    # Maps
    circ_to_row: Dict[str, int] = {}
    row_counter = 0

    data: List[int] = []
    rows: List[int] = []
    cols: List[int] = []

    # First pass: collect events
    for col_idx, (cell_id, cell_dir) in enumerate(zip(cell_ids, cell_dirs)):
        circ_file = cell_dir / f"{cell_id}_CIRCexplorer2_circ.txt"
        if not circ_file.exists():
            # Be strict: fail loudly; easier to debug
            raise FileNotFoundError(
                f"[collect_circexplorer2] Expected {circ_file} for cell {cell_id}"
            )

        events = _parse_circexplorer2_file(circ_file)

        for circ_id, count in events:
            if circ_id not in circ_to_row:
                circ_to_row[circ_id] = row_counter
                row_counter += 1
            row_idx = circ_to_row[circ_id]

            rows.append(row_idx)
            cols.append(col_idx)
            data.append(count)

    if not data:
        raise ValueError(
            "[collect_circexplorer2] No circRNA events found; all matrices would be empty."
        )

    n_circ = row_counter
    n_cells = len(cell_ids)

    # Build sparse matrix (circ x cell)
    mat = sparse.coo_matrix(
        (np.array(data, dtype=np.int64), (np.array(rows), np.array(cols))),
        shape=(n_circ, n_cells),
    ).tocsr()

    # Optional filtering by min_counts
    if min_counts > 1:
        circ_totals = np.asarray(mat.sum(axis=1)).ravel()
        keep_mask = circ_totals >= min_counts
        if keep_mask.sum() == 0:
            raise ValueError(
                f"[collect_circexplorer2] No circRNAs remain after filtering with "
                f"min_counts={min_counts}"
            )
        mat = mat[keep_mask, :]
        # rebuild circ index order
        inv_map = {row: circ_id for circ_id, row in circ_to_row.items()}
        new_circ_ids = [
            inv_map[row] for row in range(n_circ) if keep_mask[row]
        ]
    else:
        inv_map = {row: circ_id for circ_id, row in circ_to_row.items()}
        new_circ_ids = [inv_map[row] for row in range(n_circ)]

    # Write outputs
    matrix_path = outdir / "circ_matrix.mtx"
    circ_index_path = outdir / "circ_index.tsv"
    cell_index_path = outdir / "cell_index.tsv"

    mmwrite(matrix_path, mat)

    with circ_index_path.open("w") as f:
        for cid in new_circ_ids:
            f.write(f"{cid}\n")

    with cell_index_path.open("w") as f:
        for cell_id in cell_ids:
            f.write(f"{cell_id}\n")

    return matrix_path, circ_index_path, cell_index_path
