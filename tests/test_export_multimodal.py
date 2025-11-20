from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.io
import anndata as ad

from circyto.pipeline.export_multimodal import export_multimodal


def test_export_multimodal_simple(tmp_path: Path):
    # --- 1) Make a tiny gene-expression AnnData (3 cells × 2 genes) ---
    cells = ["c1", "c2", "c3"]
    genes = ["g1", "g2"]

    X_genes = sp.csr_matrix(
        [
            [1, 0],
            [0, 2],
            [3, 4],
        ]
    )

    adata = ad.AnnData(X_genes)
    adata.obs_names = cells
    adata.var_names = genes

    genes_h5ad = tmp_path / "genes.h5ad"
    adata.write_h5ad(genes_h5ad)

    # --- 2) Make a tiny circRNA matrix (3 cells × 2 circRNAs) ---
    circ_ids = ["circA", "circB"]

    # circ counts per cell (cells × circ)
    X_circ_cells = sp.csr_matrix(
        [
            [5, 0],
            [0, 6],
            [7, 8],
        ]
    )

    # circyto.collect writes circ × cell, so transpose before writing
    X_circ_circ_by_cell = X_circ_cells.T  # shape: (2 circ, 3 cells)

    circ_mtx = tmp_path / "circ.mtx"
    scipy.io.mmwrite(str(circ_mtx), X_circ_circ_by_cell)

    circ_index = tmp_path / "circ_ids.txt"
    circ_index.write_text("\n".join(circ_ids) + "\n")

    cell_index = tmp_path / "cell_ids.txt"
    cell_index.write_text("\n".join(cells) + "\n")

    out = tmp_path / "multimodal.h5ad"

    # --- 3) Run the exporter ---
    export_multimodal(
        genes_h5ad=genes_h5ad,
        circ_matrix=circ_mtx,
        circ_index=circ_index,
        cell_index=cell_index,
        out=out,
    )

    # --- 4) Validate the multimodal AnnData ---
    assert out.is_file(), "multimodal .h5ad was not created"

    adata_mm = ad.read_h5ad(out)

    # gene layer unchanged
    np.testing.assert_array_equal(adata_mm.X.toarray(), X_genes.toarray())
    assert list(adata_mm.obs_names) == cells
    assert list(adata_mm.var_names) == genes

    # circRNA modality present and aligned
    assert "X_circ" in adata_mm.obsm
    X_circ_from_mm = adata_mm.obsm["X_circ"]

    assert X_circ_from_mm.shape == X_circ_cells.shape  # (3 cells × 2 circRNAs)
    np.testing.assert_array_equal(X_circ_from_mm.toarray(), X_circ_cells.toarray())

    # circ metadata stored in .uns
    assert "circ" in adata_mm.uns
    circ_meta = adata_mm.uns["circ"]
    shape_meta = circ_meta["circ_matrix_shape"]
    assert tuple(shape_meta) == X_circ_cells.shape
    assert circ_meta.get("circ_source") == "circyto.collect" or "circ_source" in circ_meta
