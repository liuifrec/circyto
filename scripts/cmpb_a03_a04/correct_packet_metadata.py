"""One recorded correction of a stale packet label; no scientific output edits."""
import argparse
import json
from pathlib import Path
import shutil
import xml.etree.ElementTree as ET

from run_stage import digest, dump, stamp


def correct(root, repo):
    packet = repo / 'manuscript/cmpb/a03_a04/input_output_manifest.json'
    original = root / 'records/input_output_manifest_before_metadata_correction.json'
    if original.exists():
        raise FileExistsError(original)
    data = json.loads(packet.read_text())
    state = data['prevalence_integrity']
    assert state['error'] == 'file absent' and state['gzip_integrity'] is True
    assert data['A04_complete'] is True and data['prevalence_decoded_rows'] == 5159619
    assert len(state['sources']) == 3
    assert all(x['gzip_integrity'] is True and x['error'] is None for x in state['sources'])
    table = data['completed_cohort_table']
    assert digest(root / table['path']) == table['sha256']
    test = root / 'records/a04-metadata-correction-tests.xml'
    suites = list(ET.parse(test).getroot().iter('testsuite'))
    counts = {k: sum(int(s.get(k, 0)) for s in suites) for k in ['tests', 'failures', 'errors', 'skipped']}
    assert counts['tests'] == 32 and not any(counts[k] for k in ['failures', 'errors', 'skipped'])
    shutil.copyfile(packet, original)
    record = dict(recorded_utc=stamp(), original_manifest='records/'+original.name,
                  original_manifest_sha256=digest(original), corrected_field='prevalence_integrity.error',
                  before='file absent', after=None, scientific_outcomes_changed=False,
                  test_record='records/'+test.name, test_record_sha256=digest(test),
                  correction_helper_sha256=digest(__file__))
    state['error'] = None
    data['metadata_correction'] = record
    data['focused_tests'] = counts
    data['new_helper_source_hashes'] = {p.name: digest(p) for p in Path(__file__).parent.glob('*.py')}
    dump(packet, data)
    record['corrected_manifest_sha256'] = digest(packet)
    dump(root / 'records/metadata_correction.json', record)
    print(json.dumps(record))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--root', required=True, type=Path)
    ap.add_argument('--repo', required=True, type=Path)
    args = ap.parse_args()
    correct(args.root.resolve(), args.repo.resolve())
