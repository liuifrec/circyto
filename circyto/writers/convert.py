from __future__ import annotations
from pathlib import Path
from typing import List, Optional
import numpy as np
import pandas as pd
from scipy import sparse as sp
from scipy import io as scio

# Optional deps
try:
    import loompy

    HAS_LOOM = True
except Exception:
    HAS_LOOM = False

try:
    import anndata as ad

    HAS_ANNDATA = True
except Exception:
    HAS_ANNDATA = False


def _read_index_lines(p: Path) -> List[str]:
    if not p.exists():
        return []
    with p.open() as f:
        return [ln.strip() for ln in f if ln.strip()]


def to_loom(
    matrix_csr: sp.csr_matrix, circ_ids: List[str], cell_ids: List[str], out_path: str
) -> Optional[str]:
    """
    Write a loom file. loom expects (rows, cols) = (features, cells).
    Our internal matrix is (circs, cells) already, so no transpose here.
    Handle empty matrices by writing a small TSV instead.
    """
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    n_rows, n_cols = matrix_csr.shape
    if n_rows == 0 or n_cols == 0 or not HAS_LOOM:
        # Fallback: TSV snapshot for debugging — keep shapes consistent
        tsv = out.with_suffix(".empty.tsv")
        idx = circ_ids if n_rows > 0 else []
        cols = cell_ids if n_cols > 0 else []
        df = pd.DataFrame(matrix_csr.toarray(), index=idx, columns=cols)
        df.to_csv(tsv, sep="\t")
        return str(tsv)

    row_attrs = {"circ_id": np.array(circ_ids, dtype=object)}
    col_attrs = {"cell_id": np.array(cell_ids, dtype=object)}
    with loompy.new(str(out)) as ds:
        ds.add_layer("matrix", matrix_csr.astype(np.int32))
        for k, v in row_attrs.items():
            ds.ra[k] = v
        for k, v in col_attrs.items():
            ds.ca[k] = v
    return str(out)


def _to_h5ad(
    matrix_csr: sp.csr_matrix, circ_ids: List[str], cell_ids: List[str], out_path: str
) -> Optional[str]:
    """
    Write H5AD using AnnData. AnnData expects X shape (cells, features).
    We therefore pass X = matrix_csr.T (cells × circs).
    Handle empty matrices by constructing an empty CSR of shape (len(cells), len(circs)).
    """
    if not HAS_ANNDATA:
        return None

    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    n_rows, n_cols = matrix_csr.shape  # (circs, cells)
    # Construct X for AnnData: (cells, circs)
    if n_rows == 0 or n_cols == 0:
        X = sp.csr_matrix((len(cell_ids), len(circ_ids)), dtype=np.int32)
    else:
        X = matrix_csr.T.tocsr().astype(np.int32)

    obs = pd.DataFrame(index=cell_ids)
    var = pd.DataFrame(index=circ_ids)
    ad.AnnData(X=X, obs=obs, var=var).write_h5ad(str(out))
    return str(out)


def convert_matrix_files(
    matrix_path: Path,
    circ_index_path: Path,
    cell_index_path: Path,
    loom: Optional[Path] = None,
    h5ad: Optional[Path] = None,
) -> None:
    """
    Load a Matrix Market matrix (.mtx) and matching index lists,
    then write loom and/or h5ad. The .mtx is (circs, cells).
    """
    X = scio.mmread(str(matrix_path)).tocsr()
    circ_ids = _read_index_lines(circ_index_path)
    cell_ids = _read_index_lines(cell_index_path)

    # If indices were not written (e.g., all-empty), synthesize placeholders
    n_rows, n_cols = X.shape
    if not circ_ids:
        circ_ids = [f"circ_{i}" for i in range(n_rows)]
    if not cell_ids:
        cell_ids = [f"cell_{j}" for j in range(n_cols)]

    if loom:
        to_loom(X, circ_ids, cell_ids, str(loom))

    if h5ad:
        result = _to_h5ad(X, circ_ids, cell_ids, str(h5ad))
        if result is None:
            Path(str(h5ad) + ".need_anndata.txt").write_text(
                "Install anndata to enable H5AD export: pip install anndata\n"
            )
