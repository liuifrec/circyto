"""Read-only replay of the preserved 181-cell KYSE180 RNA embedding and native calls.

Only reports, explicit cell lists and validation receipts are written to --out.
No scripts from the scientific analysis are imported or executed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import resource
import sys
import time

import anndata as ad
import numpy as np
import pandas as pd
from typer.testing import CliRunner

from circyto.cli.circyto import app


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024**2), b''):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--source', type=Path, required=True,
                        help='Existing accepted primary results directory; read-only.')
    parser.add_argument('--candidate-metadata', type=Path, required=True,
                        help='Existing seven-candidate identity TSV; read-only and not bundled with the software.')
    parser.add_argument('--out', type=Path, default=Path('docs/examples/explorer/kyse180'))
    args = parser.parse_args()
    cpath, rpath = args.source/'circ_annotated.h5ad', args.source/'rna_analysis.h5ad'
    annotation = args.candidate_metadata
    source_files = [cpath, rpath, annotation, args.source/'summary.json', args.source/'candidate_summary.tsv', args.source/'cell_qc.tsv']
    before = {str(p): sha(p) for p in source_files}
    r = ad.read_h5ad(rpath, backed='r')
    ids = r.obs_names[r.obs.library_layout.eq('SINGLE')].tolist()
    assert len(ids) == 181 and r.obs.loc[ids, 'condition'].astype(str).value_counts().to_dict() == {'12 Gy': 92, '30 Gy': 89}
    args.out.mkdir(parents=True, exist_ok=True)
    cell_list = args.out/'cells.txt'; cell_list.write_text('\n'.join(ids)+'\n')
    candidates = pd.read_csv(annotation, sep='\t').circ_id.tolist()
    assert len(candidates) == 7
    command = ['report', '--explorer', str(cpath), '--rna', str(rpath), '--output-dir', str(args.out),
               '--cells', str(cell_list), '--umap-key', 'X_umap_single_primary',
               '--expression-unit', 'gene-associated single reads (not UMIs)',
               '--title', 'KYSE180 · 181 irradiated cells · preserved accepted RNA embedding',
               '--candidate-metadata', str(annotation)]
    for gene in ['MAN1A2', 'SLC8A1', 'MKI67', 'CDKN1A', 'KRT5', 'CCNB1']:
        command += ['--gene', gene]
    for circ in candidates:
        command += ['--circ', circ]
    for field in ['rna:rna_cluster_single_primary', 'rna:condition', 'rna:rna_total_assigned',
                  'rna:rna_detected_genes', 'rna:rna_unique_mapped_percent', 'rna:source_spots', 'circ:circ_observation_status']:
        command += ['--obs-column', field]
    start = time.perf_counter()
    result = CliRunner().invoke(app, command)
    elapsed = time.perf_counter() - start
    assert result.exit_code == 0, result.output
    print(result.output)
    p = json.loads((args.out/'explorer.json').read_text())
    c = ad.read_h5ad(cpath, backed='r')
    assert p['cells'] == ids and p['join']['overlap'] == 181 and p['embedding']['plotted_cells'] == 181
    np.testing.assert_array_equal(p['coordinates'], r.obsm['X_umap_single_primary'][r.obs_names.get_indexer(ids)])
    assert (c.obs.loc[ids, 'circ_observation_status'] == 'completed_zero').sum() == 2
    for feature in p['genes']:
        values = r.X[:, r.var_names.get_loc(feature['id']):r.var_names.get_loc(feature['id'])+1].toarray().ravel()
        np.testing.assert_array_equal(feature['values'], values[r.obs_names.get_indexer(ids)])
    summaries = []
    for feature in p['circs']:
        values = c.X[:, c.var_names.get_loc(feature['id']):c.var_names.get_loc(feature['id'])+1].toarray().ravel()
        selected = values[c.obs_names.get_indexer(ids)]
        np.testing.assert_array_equal(feature['values'], selected)
        assert all(s.startswith('completed_') for s in feature['states'])
        for group in ['12 Gy', '30 Gy']:
            mask = r.obs.loc[ids, 'condition'].eq(group).to_numpy()
            for threshold in [1, 2, 3]:
                summaries.append(dict(circ_id=feature['id'], condition=group, threshold=threshold,
                                      completed_cells=int(mask.sum()), detected_cells=int((selected[mask]>=threshold).sum()),
                                      detection_frequency=float((selected[mask]>=threshold).mean()), native_support=int(selected[mask].sum())))
    pd.DataFrame(summaries).to_csv(args.out/'verified_summary.tsv', sep='\t', index=False)
    c.file.close(); r.file.close()
    after = {str(p): sha(p) for p in source_files}
    assert before == after
    receipt = dict(status='pass', source_checkpoint='8f729b3', cells=181, exposures={'12 Gy':92,'30 Gy':89},
        completed_zero_cells=2, genes=len(p['genes']), circs=len(p['circs']), exact_cell_join=True,
        coordinates_exact=True, selected_native_values_exact=True, selected_expression_values_exact=True,
        input_hashes_unchanged=True, no_alignment_counting_calling_or_embedding_fit=True,
        scope='Preserved accepted X_umap_single_primary. The newer personal_RNA_181_annotated.h5ad temporary artifact is unavailable; no attempt to recreate it.',
        generation_wall_seconds=elapsed, process_peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        python=sys.version.split()[0], artifacts={name:(args.out/name).stat().st_size for name in ['report.html','explorer.json']},
        sources=[dict(name=Path(name).name, sha256=digest, unchanged=after[name]==digest) for name,digest in before.items()],
        cli=command)
    (args.out/'validation.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k not in {'cli','sources'}}, indent=2))


if __name__ == '__main__':
    main()
