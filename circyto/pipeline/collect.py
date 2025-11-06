from pathlib import Path
from typing import Dict, List, Tuple
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite

from ..parsers.cirifull import read_cirifull_tsv


def _gather_cell_tsvs(in_dir: str) -> List[Tuple[str, Path]]:
    """
    Return list of (cell_id, path_to_tsv) for all *.tsv under in_dir.
    """
    base = Path(in_dir)
    if not base.exists():
        return []
    items: List[Tuple[str, Path]] = []
    for tsv in sorted(base.glob("*.tsv")):
        cell_id = tsv.stem
        items.append((cell_id, tsv))
    return items


def _build_matrix(cirifull_dir: str):
    """
    Internal: build csr_matrix + circ_ids + cell_ids in memory.
    """
    pairs = _gather_cell_tsvs(cirifull_dir)
    if not pairs:
        return sparse.csr_matrix((0, 0), dtype=int), [], []

    cell_ids: List[str] = [cid for cid, _ in pairs]
    circ_index: Dict[str, int] = {}
    rows: List[int] = []
    cols: List[int] = []
    data: List[int] = []

    for col, (cid, path) in enumerate(pairs):
        df = read_cirifull_tsv(str(path))
        if "circ_id" not in df.columns:
            continue
        # normalize support column if present
        supports: Dict[str, int] = {}
        if "support" in df.columns:
            for _, r in df.iterrows():
                try:
                    supports[str(r["circ)id"])] = int(r["support"])
                except Exception:
                    pass
        for circ in df["circ_id"].astype(str).tolist():
            if circ not in circ_index:
                circ_index[circ] = len(circ_index)
            rows.append(circ_index[circ])
            cols.append(col)
            val = supports.get(circ, 1)
            data.append(val)

    if not data:
        return sparse.csr_matrix((0, len(cell_ids)), dtype=int), [], cell_ids

    n_rows = len(circ_index)
    n_cols = len(cell_ids)
    X = sparse.csr_matrix((data, (rows, cols)), shape=(n_rows, n_cols), dtype=int)
    circ_ids = [None] * n_rows
    for k, i in circ_index.items():
        circ_ids[i] = k
    return X, circ_ids, cell_ids


def collect_matrix(
    cirifull_dir: str,
    matrix_path: str,
    circ_index_path: str,
    cell_index_path: str,
    min_count_per_cell: int = 1,
):
    """
    Build circ-by-cell count matrix from per-cell CIRI-full TSVs and write:
      - Matrix Market sparse matrix (.mtx)
      - circ index (rows) as one ID per line
      - cell index (cols) as one ID per line

    Applies a simple per-cell filter: keep columns with sum >= min_count_per_cell.
    """
    X, circ_ids, cell_ids = _build_matrix(cirifull_dir)

    # Filter by per-cell totals
    if X.shape[1] > 0 and min_count_per_cell > 1:
        keep_cols = X.sum(axis=0).A1 >= min_count_per_cell
        if keep_cols.any():
            X = X[:, keep_cols]
            cell_ids = [cid for cid, keep in zip(cell_ids, keep_cols) if keep]
        else:
            # nothing passes; write empty consistent shapes
            X = sparse.csr_matrix((X.shape[0], 0), dtype=int)
            cell_ids = []

    # Write outputs
    out_mtx = Path(matrix_path)
    out_mtx.parent.mkdir(parents=True, exist_ok=True)
    mmwrite(out_mtx, X)

    Path(circ_index_path).write_text("\n".join(circ_ids) + ("\n" if circ_ids else ""))
    Path(cell_index_path).write_text("\n".join(cell_ids) + ("\n" if cell_ids else ""))
