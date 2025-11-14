from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
import scipy.io
import scipy.sparse as sp
import anndata as ad


def _read_lines(path: Path) -> list[str]:
    return [line.strip() for line in path.open() if line.strip()]


def _load_circ_matrix(
    matrix: Path,
    circ_index: Path,
    cell_index: Path,
):
    mtx = scipy.io.mmread(matrix)  # COO
    circ_ids = _read_lines(circ_index)
    cell_ids = _read_lines(cell_index)

    # assume matrix is circ × cell (from circyto.collect)
    if mtx.shape[0] == len(circ_ids) and mtx.shape[1] == len(cell_ids):
        X_circ = mtx.T.tocsr()  # cells × circ
    elif mtx.shape[0] == len(cell_ids) and mtx.shape[1] == len(circ_ids):
        X_circ = mtx.tocsr()     # cells × circ
    else:
        raise ValueError(
            f"Shape mismatch: matrix={mtx.shape}, "
            f"circ_ids={len(circ_ids)}, cell_ids={len(cell_ids)}"
        )

    return X_circ, circ_ids, cell_ids


def export_multimodal(
    genes_h5ad: Path,
    circ_matrix: Path,
    circ_index: Path,
    cell_index: Path,
    out: Path,
    circ_layer_name: str = "X_circ",
) -> None:
    """
    Load an existing gene-expression AnnData and attach circRNA counts as a
    separate modality in `obsm[circ_layer_name]`. Cells are aligned by name.
    """
    # 1) load base gene expression AnnData
    adata = ad.read_h5ad(genes_h5ad)

    # 2) load circRNA matrix and IDs
    X_circ, circ_ids, circ_cell_ids = _load_circ_matrix(
        circ_matrix,
        circ_index,
        cell_index,
    )

    # 3) align circRNA cells to adata.obs_names
    # Build mapping from circ cell name -> row index in X_circ
    circ_pos = {cid: i for i, cid in enumerate(circ_cell_ids)}

    # indices in circ matrix for each cell in adata
    idx = []
    missing = []
    for cell in adata.obs_names:
        if cell in circ_pos:
            idx.append(circ_pos[cell])
        else:
            missing.append(cell)

    if missing:
        # You can decide: error or partial.
        # For now, just warn by raising a clear exception.
        raise ValueError(
            f"{len(missing)} cells in genes_h5ad are missing from circ matrix, "
            f"example: {missing[:5]}"
        )

    # reindex circ matrix to match adata cell order
    X_circ_aligned = X_circ[idx, :]

    # 4) store circ counts in obsm
    adata.obsm[circ_layer_name] = X_circ_aligned

    # 5) store circ feature names + metadata
    adata.uns["circ"] = {
        "circ_ids": circ_ids,
        "circ_source": "circyto.collect",
        # Use a list instead of a tuple for HDF5/anndata friendliness
        "circ_matrix_shape": list(X_circ_aligned.shape),
    }


    out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out)
