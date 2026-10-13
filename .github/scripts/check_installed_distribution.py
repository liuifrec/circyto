"""Run outside a checkout after wheel installation; no detector is needed."""
import json
import subprocess
from pathlib import Path

import circyto
from circyto.detectors.ciri3 import Ciri3Detector
from circyto.detectors.ciri_full import CiriFullDetector
from circyto.paths import (find_ciri2_adapter, find_ciri3_jar,
                           find_ciri_full_adapter, find_ciri_full_jar, get_repo_root)

root = Path(circyto.__file__).parent
assert 'site-packages' in str(root), root
assert get_repo_root() is None
assert not list(root.rglob('*.jar')) and not list(root.rglob('*.class'))
assert not list(root.rglob('CIRI-full_manual.txt'))
assert find_ciri3_jar().resolved_path is None
assert find_ciri_full_jar().resolved_path is None
assert find_ciri2_adapter().resolved_path.is_file()
assert find_ciri_full_adapter().resolved_path.is_file()
ready, errors, details = Ciri3Detector().validate_runtime()
assert not ready and any('CIRCYTO_CIRI3_JAR' in message for message in errors)
try:
    CiriFullDetector()._resolve_ciri_full_jar()
except FileNotFoundError as exc:
    assert 'CIRCYTO_CIRI_FULL_JAR' in str(exc)
else:
    raise AssertionError('Missing CIRI-full did not fail')
for command in (['doctor'], ['detectors', '--json'], ['--help'],
                ['import-circ', '--help'], ['report', '--help'],
                ['workflow', 'full-length-circrna', '--help']):
    result = subprocess.run(['circyto', *command], capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
print(json.dumps({'installed_path': str(root), 'java_detectors_bundled': False,
                  'adapters_present': True, 'missing_detector_errors': 'actionable',
                  'cli_and_core_doctor': 'passed'}))
