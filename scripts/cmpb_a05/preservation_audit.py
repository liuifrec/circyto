"""Verify exact original bytes, accepted evidence and frozen production without execution."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

from inspect_sources import SOURCES, digest


def main():
    root=Path(sys.argv[1]);repo=Path(__file__).resolve().parents[2]
    initial=json.loads((root/'records/protected_before.json').read_text())
    changed=[p for p,sha in initial.items() if not Path(p).is_file() or digest(p)!=sha]
    refs=json.loads((root/'records/refs_before.json').read_text())
    now={line.split()[1]:line.split()[0] for line in subprocess.check_output(['git','show-ref'],cwd=repo,text=True).splitlines()}
    changed_refs={k:[v,now.get(k)] for k,v in refs.items() if now.get(k)!=v and k!='refs/heads/analysis/cmpb-a05-descriptive'}
    archive={}
    for ds,(member,sha) in SOURCES.items():
        p=subprocess.Popen(['git','cat-file','blob','c99cddaae2314d3af7d5232150c5b86bb6b52e25:load_work/'+member],cwd=repo,stdout=subprocess.PIPE)
        h=hashlib.sha256()
        for block in iter(lambda:p.stdout.read(1024*1024),b''):h.update(block)
        assert p.wait()==0
        assert h.hexdigest()==sha
        current=digest(Path('/tmp/circyto-import/load_work')/member);assert current==sha
        archive[ds]=dict(git_archive_sha256=h.hexdigest(),before_sha256=sha,after_sha256=current,unchanged=True)
    prior=Path('/home/liuyuchen/data/circyto_cmpb_a02_20261015/a03_a04_20261016/a04_combined/completed_cohort_prevalence.tsv.gz')
    full=digest(prior);assert full=='8163178301f702fb4217fc0e48d4390ed82aa39244a8d1efd2cc681ad9152b8c'
    accepted=json.loads((repo/'manuscript/cmpb/a03_a04/environment.json').read_text())
    installed=Path('/home/liuyuchen/data/circyto_cmpb_a02_20261015/status_revalidation_20261016/final_candidate/installed-env/lib/python3.10/site-packages')
    modified=[p for p,s in accepted['installed_python_members'].items() if digest(installed/p)!=s]
    assert not changed and not changed_refs and not modified,(changed,changed_refs,modified)
    result=dict(protected_files_checked=len(initial),protected_file_changes=changed,frozen_ref_changes=changed_refs,
        original_objects=archive,installed_production_source='e6843a38469d49f96d44758a4e37893ebbe2734b',
        installed_python_members_checked=len(accepted['installed_python_members']),installed_member_changes=modified,
        reused_A02_commit='307758f57a39c009393d03fc0f69e4638b2498c6',reused_A03_A04_commit='02ee3d4acfed63c5b2f02c4baef6d22c822d4dac',
        complete_prevalence_gzip=dict(sha256=full,bytes=prior.stat().st_size,rows_from_accepted_record=5159619,
            action='read-only SHA256; prior gzip CRC and row-count evidence reused without reanalysis or copying'),
        passed=True)
    (root/'results/preservation_record.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(passed=True,protected_files=len(initial),production_members=len(accepted['installed_python_members']))))


if __name__=='__main__':main()
