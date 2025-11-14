from __future__ import annotations

from pathlib import Path
from typing import Optional, Dict, List, Tuple

import numpy as np
import scipy.io
import scipy.sparse as sp
import anndata as ad
import pandas as pd


def _read_lines(path: Path) -> List[str]:
    return [line.strip() for line in path.open() if line.strip()]


def _load_circ_matrix(
    matrix: Path,
    circ_index: Path,
    cell_index: Path,
) -> Tuple[sp.csr_matrix, List[str], List[str]]:
    mtx = scipy.io.mmread(matrix)  # COO
    circ_ids = _read_lines(circ_index)
    cell_ids = _read_lines(cell_index)

    # assume matrix is circ × cell (from circyto.collect)
    if mtx.shape[0] == len(circ_ids) and mtx.shape[1] == len(cell_ids):
        X_circ = mtx.T.tocsr()  # cells × circ
    elif mtx.shape[0] == len(cell_ids) and mtx.shape[1] == len(circ_ids):
        X_circ = mtx.tocsr()  # cells × circ
    else:
        raise ValueError(
            f"Shape mismatch: matrix={mtx.shape}, "
            f"circ_ids={len(circ_ids)}, cell_ids={len(cell_ids)}"
        )

    return X_circ, circ_ids, cell_ids


def _load_circ_features(
    circ_feature_table: Optional[Path],
    circ_ids: List[str],
) -> Optional[pd.DataFrame]:
    if circ_feature_table is None:
        return None
    df = pd.read_csv(circ_feature_table, sep="\t")

    if "circ_id" not in df.columns:
        return None

    # reorder to match circ_ids; drop any circ not in matrix
    df = df.set_index("circ_id")
    missing = [cid for cid in circ_ids if cid not in df.index]
    if missing:
        # For now, we don't error; just warn via a column.
        # But we ensure we only keep those circ_ids we know.
        pass

    # align rows in exactly circ_ids order, filling missing with NaN
    df_aligned = df.reindex(circ_ids)
    df_aligned.reset_index(inplace=True)
    return df_aligned


def export_multimodal(
    genes_h5ad: Path,
    circ_matrix: Path,
    circ_index: Path,
    cell_index: Path,
    out: Path,
    circ_layer_name: str = "X_circ",
    circ_feature_table: Optional[Path] = None,
) -> None:
    """
    Load an existing gene-expression AnnData and attach circRNA counts as a
    separate modality in `obsm[circ_layer_name]`. Cells are aligned by name.

    If circ_feature_table is provided, store it in:
      - adata.uns["circ"]["feature_table"]
      - adata.uns["circ_host_map"] (circ_id -> host_gene)
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
    circ_pos = {cid: i for i, cid in enumerate(circ_cell_ids)}

    idx = []
    missing_cells = []
    for cell in adata.obs_names:
        if cell in circ_pos:
            idx.append(circ_pos[cell])
        else:
            missing_cells.append(cell)

    if missing_cells:
        raise ValueError(
            f"{len(missing_cells)} cells in genes_h5ad are missing from circ matrix, "
            f"example: {missing_cells[:5]}"
        )

    X_circ_aligned = X_circ[idx, :]

    # 4) store circ counts in obsm
    adata.obsm[circ_layer_name] = X_circ_aligned

    # 5) base circ metadata
    circ_meta = {
        "circ_ids": circ_ids,
        "circ_source": "circyto.collect",
        "circ_matrix_shape": list(X_circ_aligned.shape),
    }

    # 6) optionally attach feature table + host map
    circ_host_map: Dict[str, str] = {}
    if circ_feature_table is not None:
        df_feat = _load_circ_features(circ_feature_table, circ_ids)
        if df_feat is not None:
            circ_meta["feature_table"] = df_feat.to_dict("list")

            if "host_gene" in df_feat.columns:
                for cid, gene in zip(df_feat["circ_id"], df_feat["host_gene"]):
                    if isinstance(gene, str) and gene:
                        circ_host_map[cid] = gene

    adata.uns["circ"] = circ_meta
    if circ_host_map:
        adata.uns["circ_host_map"] = circ_host_map

    out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out)
