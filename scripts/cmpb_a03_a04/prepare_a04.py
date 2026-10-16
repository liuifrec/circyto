"""Read-only inspection of the fixed source and small derived base extraction."""
import argparse
import json
from pathlib import Path
import anndata as ad
import h5py
import numpy as np
from scipy import sparse
from scipy.io import mmwrite
from run_stage import dump, digest

def main(root):
    env=json.loads((root/'records/environment.json').read_text())
    source=Path(env['original_imr90_path'])
    if digest(source)!='bb2e12f7c3b36f9fa72d98cd71e8bea905a67f50e22af1d6b713550ee92b60c8':raise ValueError('Original identity changed')
    with h5py.File(source) as f:
        circ=ad.io.read_elem(f['mod/circ']);rna=ad.io.read_elem(f['mod/rna'])
    if circ.shape!=(23,2443) or rna.shape!=(23,63187) or set(circ.obs_names)!=set(rna.obs_names):raise ValueError('Registered source axes differ')
    if not sparse.issparse(circ.X) or circ.X.dtype.kind not in 'iu' or np.any(circ.X.data<0):raise ValueError('Nonnegative integer source required')
    totals=np.asarray(rna.X.sum(axis=1)).ravel()
    if not np.isfinite(totals).all() or np.any(totals<0):raise ValueError('Invalid stored overlap totals')
    out=root/'a04/base';out.mkdir(exist_ok=False)
    circ.write_h5ad(out/'circ_source_copy.h5ad',compression='gzip');rna.write_h5ad(out/'rna_source_copy.h5ad',compression='gzip')
    weights={str(c):float(1/(1+v)) for c,v in zip(rna.obs_names,totals)}
    dump(out/'base_record.json',dict(source_sha256=digest(source),label='DERIVED/SIMULATED; source historical statuses remain unknown',circ_cells=list(circ.obs_names),rna_cells=list(rna.obs_names),circ_shape=list(circ.shape),rna_shape=list(rna.shape),circ_nnz=circ.X.nnz,circ_support=int(circ.X.sum()),circ_var_columns=list(circ.var.columns),rna_obs_columns=list(rna.obs.columns),rna_uns_keys=list(rna.uns),stored_RNA_overlap_totals={str(c):float(v) for c,v in zip(rna.obs_names,totals)},sampling_weights=weights,source_observation_fields={k:circ.obs[k].astype(str).tolist() for k in circ.obs if 'status' in k},output_sha256={p.name:digest(p) for p in out.glob('*.h5ad')}))
    print(json.dumps(dict(circ_shape=list(circ.shape),rna_shape=list(rna.shape),circ_var_columns=list(circ.var.columns),circ_obs_columns=list(circ.obs.columns),rna_uns_keys=list(rna.uns),stored_totals_finite_nonnegative=True)))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);a=ap.parse_args();main(a.root)
