"""Close the completed controlled-job ledger; administrative bookkeeping only."""
import datetime as dt
import json
from pathlib import Path
import sys
import time

from run_stage import digest


def storage(root):
    logical=allocated=0;files=0
    for p in root.rglob('*'):
        if p.is_file():
            st=p.lstat();logical+=st.st_size;allocated+=max(st.st_size,st.st_blocks*512);files+=1
    return dict(files=files,logical_bytes=logical,allocated_or_logical_bytes=allocated)


def main():
    root=Path(sys.argv[1]);repo=Path(__file__).resolve().parents[2];out=repo/'manuscript/cmpb/a05_descriptive'
    def sanitize(text):
        for a,b in [(str(root),'${A05_ROOT}'),(str(repo),'${A05_WORKTREE}'),('/home/liuyuchen/data/circyto_cmpb_a02_20261015','${A02_ROOT}'),('/tmp/circyto-import/load_work','${ORIGINAL_OBJECT_ROOT}'),('/sys/fs/cgroup/user.slice/user-1001.slice/user@1001.service','${USER_CGROUP}')]:text=text.replace(a,b)
        return text
    stages=[]
    for p in sorted((root/'records').glob('*.json')):
        r=json.loads(p.read_text())
        if not isinstance(r,dict) or r.get('schema')!='circyto.cmpb.a05.stage.v1':continue
        assert r['status'] in ['succeeded','failed'],p
        g=json.loads((root/'records'/(r['stage']+'.guard.json')).read_text())
        props=dict(x.split('=',1) for x in r['systemd_properties'] if '=' in x)
        stages.append(dict(stage=r['stage'],status=r['status'],exit_code=r['launcher_exit_code'],started_utc=r['started_utc'],finished_utc=r['finished_utc'],
            elapsed_seconds=r['elapsed_seconds'],kernel_cpu_usage_nsec=int(props['CPUUsageNSec']) if props.get('CPUUsageNSec','').isdigit() else None,
            sampled_cgroup_memory_peak_bytes=g['memory_current_peak_bytes'],sampled_process_tree_rss_peak_bytes=g['process_tree_rss_peak_bytes'],
            observed_disk_peak_bytes=g['disk_peak_bytes'],controls=g['controls'],affinity=g['affinity'],
            interpretation='expected destructive-to-own-process control' if r['category']=='probe' else ('retained inventory stdout serialization failure, corrected by inspect-fixed' if r['stage']=='inspect' else 'A05 controlled task'),
            argv=r['argv'],helper_sha256=r['helper_sha256'],invoked_python_source_sha256=r['invoked_python_source_sha256']))
    stages.sort(key=lambda r:r['started_utc'])
    budget=json.loads((root/'records/budget.json').read_text());assert time.time()<budget['deadline_epoch']
    results_storage=storage(root);worktree_storage=storage(repo)
    total=results_storage['allocated_or_logical_bytes']+worktree_storage['allocated_or_logical_bytes']
    assert total<2*1024**3
    record=dict(schema='circyto.cmpb.a05.execution.v1',recorded_utc=dt.datetime.now(dt.timezone.utc).isoformat(),
        scope='A05 only; no A06 or sequencing rerun',budget=budget,stages=stages,
        active_controlled_stage_seconds=sum(r['elapsed_seconds'] for r in stages),
        kernel_cpu_usage_nsec_total=(sum(r['kernel_cpu_usage_nsec'] for r in stages) if all(r['kernel_cpu_usage_nsec'] is not None for r in stages) else None),
        known_stage_kernel_cpu_usage_nsec_sum=sum(r['kernel_cpu_usage_nsec'] or 0 for r in stages),
        stages_missing_kernel_cpu_accounting=[r['stage'] for r in stages if r['kernel_cpu_usage_nsec'] is None],
        cpu_accounting_limit='systemd removes successful transient service accounting before the post-exit query on this host; unavailable exact CPU totals stay null. Rounded CPU-time strings, when emitted, remain in original stage logs. CPU affinity enforcement is separately verified.',
        job_window_elapsed_seconds_at_close=time.time()-budget['started_epoch'],
        sampled_cgroup_memory_peak_bytes=max(r['sampled_cgroup_memory_peak_bytes'] for r in stages),
        sampled_process_tree_rss_peak_bytes=max(r['sampled_process_tree_rss_peak_bytes'] for r in stages),
        memory_measurement='0.25-second samples; sampled maxima are not exact kernel lifetime high-water marks. All science had hard 4 GiB cgroup memory.max and zero memory.swap.max.',
        cpu_enforcement='Inherited CPUs 0-1 plus 0.25-second all-thread affinity watchdog; cpu.max unavailable. Six enforcement probes passed; failure/timeout/OOM cleanup confirmed.',
        disk_enforcement='1-second allocated-or-logical output watchdog, stop at 1.75 GiB; final pre-export accounting includes entire new materialized worktree, retained failed/reviewed attempts, inputs and caches. ZIP/Git overhead checked separately in export receipt against 2 GiB.',
        results_root_storage=results_storage,entire_new_worktree_storage=worktree_storage,
        accounted_pre_export_allocated_or_logical_bytes=total,
        bookkeeping_scope='Small interactive source-file editing, code review, Git operations and this final scalar ledger closure are not included in controlled-stage CPU/wall totals; their files are included in storage. Export has a separate timed receipt.',
        no_background_analysis_left_running=True)
    (out/'execution_record.json').write_text(sanitize(json.dumps(record,indent=2))+'\n')
    for p in (root/'records').glob('*.json'):
        if p.name not in ['protected_before.json','refs_before.json']:(out/'execution_details'/p.name).write_text(sanitize(p.read_text()))
    for p in (root/'logs').glob('*.log'):(out/'execution_details'/p.name).write_text(sanitize(p.read_text()))
    (out/'packet_validation.json').write_text(sanitize((root/'results_review/packet_validation.json').read_text()))
    inputs={}
    for p in (root/'inputs').iterdir():
        if p.is_file():inputs[p.name]=dict(sha256=digest(p),bytes=p.stat().st_size,included_in_packet=p.suffix=='.md')
    original=json.loads((out/'preservation_record.json').read_text())['original_objects']
    files={}
    for folder in [out,repo/'scripts/cmpb_a05']:
        for p in sorted(folder.rglob('*')):
            if p.is_file() and '__pycache__' not in p.parts and p.name!='inputs_outputs_sha256.json':files[str(p.relative_to(repo))]=dict(sha256=digest(p),bytes=p.stat().st_size)
    (out/'inputs_outputs_sha256.json').write_text(json.dumps(dict(original_objects_before_after=original,public_and_authorization_inputs=inputs,
        committed_outputs_and_helpers=files,excluded_from_self_manifest='inputs_outputs_sha256.json; its hash is in the ZIP SHA256SUMS; ZIP and final commit are in handoff receipt'),indent=2)+'\n')
    print(json.dumps({k:record[k] for k in ['active_controlled_stage_seconds','kernel_cpu_usage_nsec_total','sampled_cgroup_memory_peak_bytes','sampled_process_tree_rss_peak_bytes','accounted_pre_export_allocated_or_logical_bytes']}))


if __name__=='__main__':main()
