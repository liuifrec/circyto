"""Independent raw-CSR integer oracle and complete draft source-table checks."""
import argparse
import collections
import csv
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import anndata as ad
import h5py
import numpy as np


def read(p):
    with p.open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))


def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--original',type=Path,required=True);ap.add_argument('--tables',type=Path,required=True)
    a=ap.parse_args();out=a.tables
    provenance=json.loads((out/'analysis_provenance.json').read_text())
    cells=read(out/'descriptive_cell_metrics.tsv');candidates=read(out/'descriptive_candidate_metrics.tsv')
    threshold_cells=read(out/'threshold_cell_metrics.tsv');thresholds=read(out/'threshold_sensitivity.tsv')
    edges=read(out/'circ_nonzero_entries.tsv');matrix=read(out/'matrix_audit.tsv')
    results=[]
    assertions=0
    def check(name,observed,expected):
        nonlocal assertions
        assertions+=1
        if observed!=expected:raise AssertionError((name,observed,expected))
    for ds,member in [('Smart-seq3','emtab8735_smartseq3/full_length.hostgene_fixed.h5mu'),('IMR90','scrr_imr90/full_length_rna_circ_cnv.hostgene_fixed.h5mu')]:
        before=assertions;path=a.original/member;check('source hash',sha(path),provenance['source_hashes'][ds])
        cs=[r for r in cells if r['dataset']==ds];vs=[r for r in candidates if r['dataset']==ds]
        e={(r['cell_id'],r['candidate_id']):int(r['support']) for r in edges if r['dataset']==ds}
        with h5py.File(path,'r') as f:
            global_ids=ad.io.read_elem(f['obs']).index.tolist();check('global cell order',[r['cell_id'] for r in cs],global_ids)
            totals=collections.Counter();numbers=collections.Counter();maximum=collections.Counter();cell_entries={}
            for mod in ['circ','rna']:
                g=f['mod'][mod];obs=ad.io.read_elem(g['obs']);var=ad.io.read_elem(g['var']);data=g['X/data'][:];ix=g['X/indices'][:];ptr=g['X/indptr'][:]
                check('unique ids',len(set(obs.index)),len(obs));check('unique features',len(set(var.index)),len(var))
                check('CSR orientation',tuple(g['X'].attrs['shape']),(len(obs),len(var)))
                total,nnz=0,0
                cm={r['cell_id']:r for r in cs}
                for i,cid in enumerate(obs.index):
                    row={};rawcols=ix[ptr[i]:ptr[i+1]]
                    check('no duplicate sparse coordinates',len(set(map(int,rawcols))),len(rawcols))
                    for k in range(int(ptr[i]),int(ptr[i+1])):
                        val=data[k].item();check('finite integer',bool(np.isfinite(val) and val>=0 and int(val)==val),True)
                        if val>0:row[int(ix[k])]=int(val)
                    sm=sum(row.values());nz=len(row);total+=sm;nnz+=nz
                    check(mod+' row signal',int(cm[cid]['circ_support' if mod=='circ' else 'rna_overlap_signal']),sm)
                    check(mod+' row positives',int(cm[cid]['circ_positive_candidates' if mod=='circ' else 'rna_positive_features']),nz)
                    check('completion unknown',cm[cid]['historical_circ_processing_status'],'unknown')
                    if mod=='circ':
                        cell_entries[cid]=row
                        for j,v in row.items():
                            totals[j]+=v;numbers[j]+=1;maximum[j]=max(maximum[j],v)
                            check('all sparse edge values',e.pop((cid,var.index[j])),v)
                mr=next(r for r in matrix if r['dataset']==ds and r['modality']==mod)
                check('independent Python integer total',int(mr['total_signal']),total)
                check('independent positive entries',int(mr['numerical_nonzeros']),nnz)
                if mod=='circ':
                    check('all original candidates ordered',[r['candidate_id'] for r in vs],var.index.tolist())
                    for j,r in enumerate(vs):
                        check('candidate pooled',int(r['pooled_support']),totals[j]);check('candidate positive cells',int(r['positive_stored_cells']),numbers[j]);check('candidate max',int(r['maximum_per_cell_support']),maximum[j])
                        check('stored denominator',int(r['stored_cell_denominator']),len(obs))
                        check('stored fraction',float(r['stored_cell_fraction']),numbers[j]/len(obs))
                        check('undefined completed denominator',r['completed_cohort_denominator'],'');check('undefined completed prevalence',r['completed_cohort_prevalence'],'')
                        for field in var.columns:
                            check('all original annotation fields',r[field],str(var.iloc[j][field]))
            check('no spurious sparse entries',len(e),0)
            for mode in ['individual_entry','pooled_candidate']:
                for t in [1,2,3,5]:
                    source=[r for r in threshold_cells if r['dataset']==ds and r['threshold_mode']==mode and int(r['threshold'])==t]
                    check('threshold keeps all stored rows',len(source),len(cell_entries))
                    ns,ss,features=[],[],set()
                    for r in source:
                        keep={j:v for j,v in cell_entries[r['cell_id']].items() if (v>=t if mode=='individual_entry' else totals[j]>=t)}
                        check('independent threshold cell candidates',int(r['positive_candidates']),len(keep));check('independent threshold cell support',int(r['support']),sum(keep.values()))
                        ns.append(len(keep));ss.append(sum(keep.values()));features.update(keep)
                    row=next(r for r in thresholds if r['dataset']==ds and r['threshold_mode']==mode and int(r['threshold'])==t)
                    for field,value in [('retained_candidates',len(features)),('positive_cell_candidate_entries',sum(ns)),('retained_support',sum(ss)),('stored_zero_rows',ss.count(0))]:check('threshold aggregate '+field,int(row[field]),value)
                    check('threshold candidate median',float(row['median_positive_candidates_per_stored_cell']),float(np.median(ns)))
                    check('threshold support median',float(row['median_support_per_stored_cell']),float(np.median(ss)))
            if ds=='Smart-seq3':
                padded=next(r for r in cs if r['cell_id']=='DIYHEK_192')
                check('padded member','rna_only',padded['recorded_membership']);check('padded circ zero',int(padded['circ_support']),0)
                check('padded match excluded',padded['matched_rna_circ'],'False');check('padded physical row retained',padded['stored_circ_row_present'],'True')
        check('source after hash',sha(path),provenance['source_hashes'][ds])
        results.append(dict(dataset=ds,independent_assertions=assertions-before,method='raw CSR indptr/indices/data with Python integer dictionaries; no scipy matrix reductions',passed=True))
    # Every plotted observation has a complete deterministic source-table match.
    plotted=read(out/'figure_observations.tsv');manifest=json.loads((out/'figure_manifest.json').read_text())
    for fig in manifest:
        name=fig['figure'];ds=fig['dataset'];rr=[r for r in plotted if r['figure']==name]
        panels=collections.defaultdict(list)
        for r in rr:panels[r['panel']].append(r)
        cs=[r for r in cells if r['dataset']==ds];vs=[r for r in candidates if r['dataset']==ds]
        for p,rows in panels.items():
            panel=int(p)
            if 'scatter' in name:
                sources=[r for r in cs if r['matched_rna_circ']=='True'];xm,ym=[('rna_overlap_signal','circ_support'),('rna_positive_features','circ_positive_candidates')][panel]
                expected=[(r['cell_id'],float(r[xm]),float(r[ym])) for r in sources]
            elif 'distributions' in name:
                source,field=[(cs,'circ_positive_candidates'),(cs,'circ_support'),(vs,'positive_stored_cells'),(vs,'pooled_support')][panel]
                order=sorted(enumerate(source),key=lambda z:(int(z[1][field]),z[0]))
                expected=[(r['cell_id'] if panel<2 else r['candidate_id'],float(r[field]),(i+1)/len(source)) for i,(_,r) in enumerate(order)]
            elif 'threshold' in name:
                field=['retained_candidates','retained_support'][panel]
                source=[r for r in thresholds if r['dataset']==ds]
                expected=[(r['threshold_mode']+':'+r['threshold'],float(r['threshold']),float(r[field])) for r in source]
            else:
                field=['circ_positive_candidates','circ_support','rna_positive_features','rna_overlap_signal'][panel]
                expected=[(r['cell_id'],float(r[field]),float(i)) for i,r in enumerate(cs)]
            check('all plot observations '+name+p,[(r['observation_id'],float(r['x']),float(r['y'])) for r in rows],expected)
        tree=ET.parse(out/name);ns={'s':'http://www.w3.org/2000/svg'}
        if 'distributions' not in name:
            circles=tree.findall('.//s:circle',ns)
            check('SVG point count',len(circles),len(rr))
            check('SVG raw plotted values',[(c.attrib['data-id'],float(c.attrib['data-x']),float(c.attrib['data-y'])) for c in circles],[(r['observation_id'],float(r['x']),float(r['y'])) for r in rr])
        else:
            paths=tree.findall('.//s:path',ns);check('ECDF path count',len(paths),4)
            for i,path in enumerate(paths):check('ECDF retained step count',path.attrib['d'].count('L'),2*len(panels[str(i)]))
    check('no invalid CNV display table',(out/'cnv_display_bins.tsv').exists(),False)
    gate=json.loads((out/'cnv_gate.json').read_text());check('CNV gate',gate['decision'],'NOT VERIFIED')
    (out/'validation_record.json').write_text(json.dumps(dict(passed=True,assertions=assertions,scientific_oracles=results,
        plotted_observations_checked=len(plotted),figures_checked=len(manifest),source_hashes_before_after=provenance['source_hashes'],
        limitations='Assertions are technical checks, not independent biological replications.'),indent=2)+'\n')
    print(json.dumps(dict(passed=True,assertions=assertions,plotted_observations=len(plotted))))


if __name__=='__main__':main()
