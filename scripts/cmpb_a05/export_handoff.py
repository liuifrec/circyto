"""Create one exclusive, verified ZIP of committed A05 files in Downloads."""
import argparse
import datetime as dt
import hashlib
import json
from pathlib import Path
import subprocess
import time
import zipfile

from accounting import storage
from run_stage import digest

BASE='02ee3d4acfed63c5b2f02c4baef6d22c822d4dac'
PREFIXES=('manuscript/cmpb/a05_descriptive/','scripts/cmpb_a05/')


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--destination-base',type=Path,required=True)
    a=ap.parse_args();repo=Path(__file__).resolve().parents[2];started=time.monotonic()
    def git(*args):return subprocess.check_output(['git',*args],cwd=repo)
    head=git('rev-parse','HEAD').decode().strip();branch=git('branch','--show-current').decode().strip()
    assert branch=='analysis/cmpb-a05-descriptive';assert not git('status','--porcelain').strip()
    subprocess.run(['git','merge-base','--is-ancestor',BASE,head],cwd=repo,check=True)
    changes=git('diff','--name-only',BASE,head).decode().splitlines();assert all(p.startswith(PREFIXES) for p in changes)
    names=git('ls-tree','-r','--name-only',head,'--',*PREFIXES).decode().splitlines()
    assert names and all(Path(n).suffix in ['.py','.md','.json','.tsv','.svg','.png','.log'] for n in names)
    assert not any('cnv_display_bins' in n or '__pycache__' in n for n in names)
    files={n:git('show',head+':'+n) for n in names}
    assert max(map(len,files.values()))<8*1024**2
    # Count introduced Git object disk bytes in addition to the entire new checkout.
    objects=git('rev-list','--objects',BASE+'..'+head).decode().splitlines()
    ids=[x.split()[0] for x in objects]
    sizes=subprocess.check_output(['git','cat-file','--batch-check=%(objectsize:disk)'],input=('\n'.join(ids)+'\n').encode(),cwd=repo)
    git_bytes=sum(int(v) for v in sizes.splitlines())
    meta_path=Path(git('rev-parse','--absolute-git-dir').decode().strip())
    meta_bytes=storage(meta_path)['allocated_or_logical_bytes']
    existing=storage(a.root)['allocated_or_logical_bytes']+storage(repo)['allocated_or_logical_bytes']+git_bytes+meta_bytes
    # Conservative allowance for the ZIP (uncompressed member bytes plus 1 MiB), records and directory metadata.
    reserve=sum(map(len,files.values()))+2*1024**2
    assert existing+reserve<2*1024**3,(existing,reserve)
    final_time=dt.datetime.now(dt.timezone.utc)
    directory=a.destination_base/(final_time.strftime('handoff_%Y%m%dT%H%M%SZ_')+head[:12])
    directory.mkdir(parents=True,exist_ok=False)
    filename='CIRCYTO_CMPB_A05_'+head[:12]+'.zip';path=directory/filename
    readme=f'''CIRCYTO CMPB A05 handoff
Final packet commit: {head}
Branch: {branch}
Accepted base: {BASE}
Created UTC: {final_time.isoformat()}

Read manuscript/cmpb/a05_descriptive/README.md and a05_results_for_review.md.
CNV eligibility is NOT VERIFIED. Review the RNA/circ-only Figure 2C alternative.
All files below the manuscript/ and scripts/ prefixes are exact committed Git blobs.
SHA256SUMS covers every other ZIP member, including this README.
The outer .sha256 and handoff_receipt.json record ZIP identity and final storage.
No raw reads, H5MU, full CNV arrays, environments or genome indices are included.
A02-A04 evidence was reused; A06 was not started.
'''
    files['HANDOFF_README.txt']=readme.encode()
    sums={n:hashlib.sha256(b).hexdigest() for n,b in files.items()}
    files['SHA256SUMS']=''.join(sums[n]+'  '+n+'\n' for n in sorted(sums)).encode()
    with zipfile.ZipFile(path,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as z:
        for n in sorted(files):
            info=zipfile.ZipInfo(n,date_time=(2026,10,19,0,0,0));info.compress_type=zipfile.ZIP_DEFLATED
            info.external_attr=0o100644<<16;z.writestr(info,files[n])
    with zipfile.ZipFile(path) as z:
        assert z.testzip() is None
        assert len(z.namelist())==len(set(z.namelist()))==len(files)
        for n,b in files.items():assert z.read(n)==b,n
        for line in z.read('SHA256SUMS').decode().splitlines():
            expected,n=line.split('  ',1);assert hashlib.sha256(z.read(n)).hexdigest()==expected
    sha=digest(path)
    (directory/(filename+'.sha256')).write_text(sha+'  '+filename+'\n')
    receipt=dict(zip_path=str(path),zip_bytes=path.stat().st_size,zip_sha256=sha,commit=head,branch=branch,
        archive_members=len(files),committed_files=len(names),zip_crc_and_all_member_sha256_verified=True,
        byte_for_byte_matches_git=True,export_wall_seconds=time.monotonic()-started,
        introduced_git_object_disk_bytes=git_bytes,new_worktree_git_metadata_bytes=meta_bytes,
        pre_export_storage_bytes_including_full_worktree_and_git=existing,
        storage_bytes_after_zip=existing+storage(directory)['allocated_or_logical_bytes'],
        storage_cap_bytes=2*1024**3,
        scope='A05 only; science execution accounting is committed; export service accounting is appended after service completion')
    assert receipt['storage_bytes_after_zip']<receipt['storage_cap_bytes']
    (directory/'handoff_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    (a.root/'records/handoff_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt))


if __name__=='__main__':main()
