from pathlib import Path
from typing import Optional
import pandas as pd
from scipy import io as scio
import loompy as lp
import anndata as ad
from rich.console import Console

console = Console()

def convert_matrix_files(
matrix_path: Path,
circ_index_path: Path,
cell_index_path: Path,
loom: Optional[Path] = None,
h5ad: Optional[Path] = None,
):
M = scio.mmread(matrix_path).tocsr()
circ = [line.strip() for line in open(circ_index_path)]
cells = [line.strip() for line in open(cell_index_path)]
assert M.shape == (len(circ), len(cells)), f"Matrix {M.shape} != ({len(circ)}, {len(cells)})"

if loom:
    loom = Path(loom); loom.parent.mkdir(parents=True, exist_ok=True)
    row_attrs = {"circ_id": circ}
    col_attrs = {"CellID": cells}
    lp.create(str(loom), layers={"circ_counts": M.tocsc()}, row_attrs=row_attrs, col_attrs=col_attrs)
    console.print(f"[bold]Wrote Loom[/bold] → {loom}")

if h5ad:
    h5ad = Path(h5ad); h5ad.parent.mkdir(parents=True, exist_ok=True)
    adata = ad.AnnData(X=M.T)  # cells x circ
    adata.var_names = pd.Index(circ, dtype="string")
    adata.obs_names = pd.Index(cells, dtype="string")
    adata.layers["circ_counts"] = adata.X.copy()
    adata.write_h5ad(h5ad)
    console.print(f"[bold]Wrote H5AD[/bold] → {h5ad}")


