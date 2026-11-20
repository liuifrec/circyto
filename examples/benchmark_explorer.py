"""Measure bounded reporting on large sparse feature axes; no analysis pipeline."""
from pathlib import Path
import argparse
import json
import resource
import time
import tempfile

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from circyto.pipeline.single_cell_explorer import ExplorerOptions, load_explorer
from circyto.pipeline.qc_report import generate_report


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=Path('docs/examples/explorer/validation/scaling.json'))
    args = parser.parse_args()
    n, ng, nc = 1000, 200_000, 20_000
    rng = np.random.default_rng(1729)
    def matrix(features, nnz):
        return sparse.coo_matrix((rng.integers(1,20,nnz),(rng.integers(n,size=nnz),rng.integers(features,size=nnz))),shape=(n,features)).tocsr()
    obs = pd.DataFrame({'detector_status':'completed','condition':np.where(np.arange(n)%2,'A','B')},index=[f'cell{i}' for i in range(n)])
    rna = ad.AnnData(matrix(ng, 40_000), obs=obs.copy(), var=pd.DataFrame(index=[f'GENE{i:06d}' for i in range(ng)]))
    rna.obsm['X_umap'] = rng.normal(size=(n,2))
    circ = ad.AnnData(matrix(nc, 4000), obs=obs.copy(), var=pd.DataFrame(index=[f'junction{i:05d}' for i in range(nc)]))
    with tempfile.TemporaryDirectory(prefix='circyto-explorer-benchmark-') as temp:
        root = Path(temp);cp,rp=root/'circ.h5ad',root/'rna.h5ad'
        circ.write_h5ad(cp);rna.write_h5ad(rp)
        start = time.perf_counter()
        p = load_explorer(cp,rp,options=ExplorerOptions(max_genes=8,max_circs=8,title='Sparse dimensionality benchmark'))
        generate_report(None,explorer=p,output_dir=root/'report')
        duration = time.perf_counter()-start
        record=dict(status='pass',cells=n,genes_available=ng,circs_available=nc,genes_embedded=len(p['genes']),circs_embedded=len(p['circs']),
                    embedded_measurements=n*(len(p['genes'])+len(p['circs'])),wall_seconds=duration,
                    process_peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                    html_bytes=(root/'report/report.html').stat().st_size,json_bytes=(root/'report/explorer.json').stat().st_size,
                    matrix_format='backed CSR; only selected feature columns materialized',input_rna_nnz=rna.X.nnz,input_circ_nnz=circ.X.nnz)
        assert record['embedded_measurements']==16000 and record['html_bytes']<1_000_000
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2))


if __name__=='__main__':
    main()
