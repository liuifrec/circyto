from pathlib import Path
from typing import List, Optional
import numpy as np
from scipy import sparse
from scipy.io import mmread

# Optional dependencies
try:
    import loompy

    HAS_LOOM = True
except Exception:
    HAS_LOOM = False

try:
    import anndata as ad
    import pandas as pd

    HAS_ANNDATA = True
except Exception:
    HAS_ANNDATA = False
    import pandas as pd  # fallback for dense TSV export if needed


def _read_index(path: str) -> List[str]:
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return []
    return [line.strip() for line in p.read_text().splitlines() if line.strip()]


def to_loom(
    matrix_csr: sparse.csr_matrix,
    circ_ids: List[str],
    cell_ids: List[str],
    out_path: str,
) -> str:
    """
    Write a minimal Loom file. If matrix is empty or loompy is missing,
    write a friendly .empty.txt (or TSV) and return its path.
    """
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    # If empty (no rows or no cols), don't try to create a loom file.
    if matrix_csr.shape[0] == 0 or matrix_csr.shape[1] == 0:
        note = out.with_suffix(".empty.txt")
        note.write_text(
            "Matrix is empty (no circRNAs or no cells). Loom not created.\n"
        )
        return str(note)

    if HAS_LOOM:
        # loompy >=3: use create() with row/col attrs
        loompy.create(
            str(out),
            matrix_csr.astype(np.int32),
            row_attrs={"circ_id": np.array(circ_ids, dtype=object)},
            col_attrs={"cell_id": np.array(cell_ids, dtype=object)},
        )
        return str(out)

    # Fallback: TSV dense export (only for small matrices)
    df = pd.DataFrame(matrix_csr.toarray(), index=circ_ids, columns=cell_ids)
    tsv = out.with_suffix(".tsv")
    df.to_csv(tsv, sep="\t")
    return str(tsv)


def _to_h5ad(
    matrix_csr: sparse.csr_matrix,
    circ_ids: List[str],
    cell_ids: List[str],
    out_path: str,
) -> Optional[str]:
    if not HAS_ANNDATA:
        return None
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    var = pd.DataFrame(index=circ_ids)
    obs = pd.DataFrame(index=cell_ids)
    ad.AnnData(X=matrix_csr, var=var, obs=obs).write_h5ad(str(out))
    return str(out)


def convert_matrix_files(
    matrix_path: str,
    circ_index_path: str,
    cell_index_path: str,
    loom: Optional[str] = None,
    h5ad: Optional[str] = None,
) -> None:
    """
    Load MatrixMarket (.mtx) + index files and write loom/h5ad outputs.
    """
    X = mmread(str(matrix_path))
    if not sparse.issparse(X):
        X = sparse.csr_matrix(X)
    else:
        X = X.tocsr()

    circ_ids = _read_index(circ_index_path)
    cell_ids = _read_index(cell_index_path)

    n_rows, n_cols = X.shape
    if len(circ_ids) != n_rows:
        circ_ids = [f"circ_{i}" for i in range(n_rows)]
    if len(cell_ids) != n_cols:
        cell_ids = [f"cell_{j}" for j in range(n_cols)]

    if loom:
        to_loom(X, circ_ids, cell_ids, str(loom))

    if h5ad:
        result = _to_h5ad(X, circ_ids, cell_ids, str(h5ad))
        if result is None:
            Path(h5ad).with_suffix(".need_anndata.txt").write_text(
                "Install anndata to enable H5AD export: pip install anndata\n"
            )

    # If neither loom nor h5ad requested, emit a dense TSV as a generic fallback
    if not loom and not h5ad:
        out = Path(matrix_path).with_suffix(".dense.tsv")
        df = pd.DataFrame(X.toarray(), index=circ_ids, columns=cell_ids)
        df.to_csv(out, sep="\t")
