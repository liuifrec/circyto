"""Check portable figure recreation, quantiles and all proposed handoff file types."""
import csv
import json
import math
from pathlib import Path
import subprocess
import sys

from inspect_sources import digest


def read(p):
    with p.open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))


def main():
    root=Path(sys.argv[1]);repo=Path(__file__).resolve().parents[2];out=repo/'manuscript/cmpb/a05_descriptive'
    target=root/'portable_reproduction';assert not target.exists()
    manifest=json.loads((out/'figure_manifest.json').read_text())
    command=[sys.executable,'-B',str(Path(__file__).with_name('plot_drafts.py')),'--tables',str(out),'--out',str(target),'--renderer-commit',manifest[0]['renderer_commit']]
    subprocess.run(command,check=True)
    names=[r['figure'] for r in manifest]+['figure_manifest.json','figure_observations.tsv']
    for name in names:assert digest(target/name)==digest(out/name),name
    copied=[]
    for name in ['descriptive_cell_metrics','descriptive_candidate_metrics','modality_membership_audit','threshold_sensitivity','threshold_cell_metrics','distribution_summary','registration_checks','circ_nonzero_entries','matrix_audit']:
        assert digest(out/(name+'.tsv'))==digest(root/'results'/(name+'.tsv'))
        copied.append(name+'.tsv')
    cells=read(out/'descriptive_cell_metrics.tsv');candidates=read(out/'descriptive_candidate_metrics.tsv');checks=0
    for row in read(out/'distribution_summary.tsv'):
        pool=cells if row['observation_unit']=='cell' else candidates
        pool=[x for x in pool if x['dataset']==row['dataset']]
        if row['subset']=='recorded_circ_members':pool=[x for x in pool if x['recorded_circ_member']=='True']
        values=sorted(float(x[row['metric']]) for x in pool)
        assert int(row['n_total'])==len(values) and int(row['n_missing'])==0
        for field,q in [('minimum',0),('q25',.25),('median',.5),('q75',.75),('maximum',1)]:
            pos=(len(values)-1)*q;i=math.floor(pos);fraction=pos-i
            expected=values[i]+fraction*(values[min(i+1,len(values)-1)]-values[i])
            assert math.isclose(float(row[field]),expected,rel_tol=1e-14,abs_tol=1e-14),(row,field)
            checks+=1
    for row in json.loads((out/'png_render_record.json').read_text()):
        assert digest(out/row['svg'])==row['svg_sha256'];assert digest(out/row['png'])==row['png_sha256']
    forbidden=[]
    allowed={'.tsv','.json','.md','.svg','.png','.log'}
    for p in out.rglob('*'):
        if p.is_file() and (p.suffix not in allowed or p.stat().st_size>8*1024**2 or p.is_symlink()):forbidden.append(str(p))
    assert not forbidden
    record=dict(passed=True,byte_identical_recreated_files=names,unchanged_scientific_tables=copied,
        independent_sorted_linear_quantile_checks=checks,forbidden_files=forbidden,
        no_h5mu_raw_reads_or_full_cnv_arrays=True,plot_reproduction_command=command)
    (root/'results_review/packet_validation.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(dict(passed=True,recreated_files=len(names),quantiles_checked=checks)))


if __name__=='__main__':main()
