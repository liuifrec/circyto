"""Build real artifacts and verify the external-detector distribution boundary."""
import shutil
import subprocess
import sys
import tarfile
import zipfile
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def source_copy(target):
    target.mkdir()
    for name in ('circyto', 'docs', 'examples'):
        shutil.copytree(ROOT / name, target / name,
                        ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
    for name in ('pyproject.toml', 'setup.py', 'MANIFEST.in', 'README.md',
                 'LICENSE', 'THIRD_PARTY_NOTICES.md', 'CITATION.cff'):
        shutil.copyfile(ROOT / name, target / name)
    return target


def run_build(source, *args):
    return subprocess.run([sys.executable, '-m', 'build', '--no-isolation', *args],
                          cwd=source, capture_output=True, text=True)


def check_members(names):
    assert not any(Path(n).suffix.lower() in {'.jar', '.class'} for n in names)
    assert not any(n.endswith('CIRI-full_manual.txt') for n in names)
    for suffix in ('circyto/resources/tools/CIRI-full_v2.0/bin/ciri_full_adapter.sh',
                   'circyto/resources/tools/CIRI-full_v2.0/bin/CIRI_v2.0.6/CIRI2.pl',
                   'circyto/resources/tools/CIRI-full_v2.0/bin/CIRI_v2.0.6/LICENSE.txt',
                   'circyto/resources/tools/CIRI-full_v2.0/bin/CIRI_AS_v1.2/CIRI_AS_v1.2.pl',
                   'circyto/resources/tools/CIRI-full_v2.0/bin/CIRI_AS_v1.2/LICENSE.txt',
                   'circyto/pipeline/import_ciri3.py'):
        assert any(n.endswith(suffix) for n in names), suffix


def test_wheel_sdist_and_rebuilt_wheel_exclude_unresolved_assets(tmp_path):
    source = source_copy(tmp_path / 'source')
    result = run_build(source, '--sdist', '--wheel')
    assert result.returncode == 0, result.stdout + result.stderr
    with zipfile.ZipFile(next((source / 'dist').glob('*.whl'))) as archive:
        check_members(archive.namelist())
        assert archive.testzip() is None
    expanded = tmp_path / 'expanded'
    with tarfile.open(next((source / 'dist').glob('*.tar.gz'))) as archive:
        check_members(archive.getnames())
        # This archive was just built from our own test source tree.
        archive.extractall(expanded, filter='data')
    result = run_build(next(expanded.iterdir()), '--wheel')
    assert result.returncode == 0, result.stdout + result.stderr
    with zipfile.ZipFile(next(expanded.glob('*/dist/*.whl'))) as archive:
        check_members(archive.namelist())


@pytest.mark.parametrize('name', ['stale.jar', 'stale.class', 'CIRI-full_manual.txt'])
def test_stale_build_assets_fail_closed(tmp_path, name):
    source = source_copy(tmp_path / 'source')
    stale = source / 'build/lib/circyto/resources/tools' / name
    stale.parent.mkdir(parents=True)
    stale.write_bytes(b'stale excluded asset')
    result = run_build(source, '--wheel')
    assert result.returncode != 0
    assert 'Excluded detector assets found in build output' in result.stdout + result.stderr
    assert not list((source / 'dist').glob('*.whl'))
