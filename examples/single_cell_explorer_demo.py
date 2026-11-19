"""Generate a fully synthetic offline explorer, including incomplete observations."""
from pathlib import Path
import argparse
import tempfile

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from circyto.pipeline.single_cell_explorer import ExplorerOptions, load_explorer
from circyto.pipeline.qc_report import generate_report


def make_objects():
    rng = np.random.default_rng(1729)
    n = 96
    ids = [f"synthetic_cell_{i:03d}" for i in range(n)]
    groups = np.repeat(np.arange(3), 32)
    xy = rng.normal(size=(n, 2)) * [1.2, 0.8] + np.array([[-3, 2], [1, -2], [4, 3]])[groups]
    expression = rng.poisson(np.array([3, 6, 2, 4, 8, 2])[None, :] * (1+groups[:, None]*0.3))
    names = ["MAN1A2", "MKI67", "SLC8A1", "CDKN1A", "DUPLICATE", "DUPLICATE"]
    rna = ad.AnnData(sparse.csr_matrix(expression),
        obs=pd.DataFrame({"RNA_cluster": [f"RNA group {x+1}" for x in groups],
                          "condition": np.where(np.arange(n)%2, "exposed", "control"),
                          "detected_genes": (expression>0).sum(1), "total_counts": expression.sum(1),
                          "mitochondrial_fraction": rng.uniform(.01,.12,n)}, index=ids),
        var=pd.DataFrame({"gene_name": names, "gene_id": [f"ENSG{i+1:011d}" for i in range(6)]},
                         index=[f"ENSG{i+1:011d}" for i in range(6)]))
    rna.obsm["X_umap"] = xy
    rna.layers["log1p_demo"] = sparse.csr_matrix(np.log1p(expression))
    support = rng.poisson(np.array([1.2, .5, .8])[None,:]*(1+groups[:,None]))
    support[[0, 10, 20, 30]] = 0
    states = ["completed"]*n
    for index, state in [(1,"failed"),(2,"unknown"),(3,"unprocessed"),(4,"failed"),(5,"missing")]:
        states[index] = state
    circ = ad.AnnData(sparse.csr_matrix(support),
        obs=pd.DataFrame({"detector_status": states, "condition": rna.obs.condition}, index=ids),
        var=pd.DataFrame({"chrom":["chr1","chr1","chr2"],"start":[100,300,500],"end":[200,400,600],
                          "strand":["+","+","-"],"host_gene":["MAN1A2","MAN1A2","SLC8A1"],
                          "circAtlas_id":["synthetic-MAN1A2_001","synthetic-MAN1A2_002","synthetic-SLC8A1_001"],
                          "annotation_source":["Synthetic catalogue fixture"]*3,
                          "genome_build":["synthetic coordinates"]*3, "coordinate_convention":["1-based-closed"]*3,
                          "support_definition":["synthetic junction reads"]*3},
                         index=["chr1:100|200","chr1:300|400","chr2:500|600"]))
    # Independent row orders and partial overlap are deliberate.
    return circ[rng.permutation(n)[1:]].copy(), rna[rng.permutation(n)[1:]].copy()


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--out", type=Path, default=Path("docs/examples/explorer/synthetic"))
    args = parser.parse_args()
    c, r = make_objects()
    with tempfile.TemporaryDirectory(prefix="circyto-explorer-demo-") as temp:
        cp, rp = Path(temp)/"synthetic_circ.h5ad", Path(temp)/"synthetic_rna.h5ad"
        c.write_h5ad(cp); r.write_h5ad(rp)
        p = load_explorer(cp, rp, options=ExplorerOptions(title="Synthetic demonstration · illustrative data, not biological evidence",
            genes=("MAN1A2","MKI67","SLC8A1","CDKN1A","DUPLICATE"), rna_layer="log1p_demo", expression_unit="log1p synthetic gene counts"))
        generate_report(None, explorer=p, output_dir=args.out)
    print(args.out/'report.html')


if __name__ == "__main__":
    main()
