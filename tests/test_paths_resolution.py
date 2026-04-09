from __future__ import annotations

import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

from circyto.cli.detectors import detectors_app
from circyto.cli.doctor import doctor_app
from circyto.cli.doctor_meta import detector_runtime_status
from circyto.detectors.base import DetectorRunInputs
from circyto.detectors.ciri_full import CiriFullDetector
from circyto.paths import (
    clear_path_resolution_caches,
    find_ciri3_jar,
    find_ciri_full_jar,
    get_bundled_smoke_testdata_dir,
    get_repo_root,
    resolve_ciri3_installation,
)
from circyto.pipeline.run_cirifull import run_cirifull_with_manifest


runner = CliRunner()


@pytest.fixture(autouse=True)
def _clear_path_caches() -> None:
    clear_path_resolution_caches()
    yield
    clear_path_resolution_caches()


def test_find_ciri_full_jar_is_cwd_independent(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    repo_root = get_repo_root()
    assert repo_root is not None

    monkeypatch.chdir(tmp_path)
    resolution = find_ciri_full_jar()

    assert resolution.resolved_path == repo_root / "tools" / "CIRI-full_v2.0" / "CIRI-full.jar"
    assert resolution.source in {"repo-tools", "repo-tools-scan"}


def test_get_bundled_smoke_testdata_dir_is_cwd_independent(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    repo_root = get_repo_root()
    assert repo_root is not None

    monkeypatch.chdir(tmp_path)
    smoke_dir = get_bundled_smoke_testdata_dir()

    assert smoke_dir == repo_root / "testdata" / "smartseq2_smoke"
    assert smoke_dir.exists()


def test_find_ciri_full_jar_env_override(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    fake_jar = tmp_path / "custom" / "CIRI-full.jar"
    fake_jar.parent.mkdir(parents=True, exist_ok=True)
    fake_jar.write_text("jar", encoding="utf-8")

    monkeypatch.setenv("CIRCYTO_CIRI_FULL_JAR", str(fake_jar))
    resolution = find_ciri_full_jar()

    assert resolution.resolved_path == fake_jar
    assert resolution.source == "env:CIRCYTO_CIRI_FULL_JAR"


def test_missing_ciri_full_jar_reports_checked_paths(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    fake_repo = tmp_path / "repo"
    (fake_repo / "circyto").mkdir(parents=True)
    (fake_repo / "pyproject.toml").write_text("[project]\nname='circyto'\n", encoding="utf-8")

    monkeypatch.setenv("CIRCYTO_REPO_ROOT", str(fake_repo))
    resolution = find_ciri_full_jar()

    assert resolution.resolved_path is None
    assert fake_repo / "tools" / "CIRI-full_v2.0" / "CIRI-full.jar" in resolution.checked_paths
    assert fake_repo / "tools" / "CIRI-full_v2.0" / "CIRI_Full.jar" in resolution.checked_paths


def test_find_ciri3_jar_prefers_env_override(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    fake_home = tmp_path / "ciri3"
    fake_home.mkdir(parents=True, exist_ok=True)
    fake_jar = fake_home / "CIRI3_Java_18.0.1.jar"
    fake_jar.write_text("jar", encoding="utf-8")

    monkeypatch.setenv("CIRCYTO_CIRI3_JAR", str(fake_jar))
    resolution = find_ciri3_jar()

    assert resolution.resolved_path == fake_jar
    assert resolution.source == "env:CIRCYTO_CIRI3_JAR"


def test_resolve_ciri3_installation_is_cwd_independent(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    fake_home = tmp_path / "tools" / "CIRI3"
    fake_home.mkdir(parents=True, exist_ok=True)
    fake_jar = fake_home / "CIRI3_Java_18.0.1.jar"
    fake_jar.write_text("jar", encoding="utf-8")
    fake_java = tmp_path / "bin" / "java"
    fake_java.parent.mkdir(parents=True, exist_ok=True)
    fake_java.write_text("java", encoding="utf-8")

    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir(parents=True, exist_ok=True)
    monkeypatch.chdir(elsewhere)
    monkeypatch.setenv("CIRCYTO_CIRI3_HOME", str(fake_home))
    monkeypatch.setenv("CIRCYTO_CIRI3_JAVA", str(fake_java))

    resolution = resolve_ciri3_installation()

    assert resolution.home.resolved_path == fake_home
    assert resolution.jar.resolved_path == fake_jar
    assert resolution.java.resolved_path == fake_java


def test_doctor_and_detectors_report_ciri3_runtime_consistently(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    fake_home = tmp_path / "ciri3"
    fake_home.mkdir(parents=True, exist_ok=True)
    fake_jar = fake_home / "CIRI3_Java_18.0.1.jar"
    fake_jar.write_text("jar", encoding="utf-8")
    fake_java = tmp_path / "java"
    fake_java.write_text("java", encoding="utf-8")

    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("CIRCYTO_CIRI3_HOME", str(fake_home))
    monkeypatch.setenv("CIRCYTO_CIRI3_JAVA", str(fake_java))

    status = detector_runtime_status("ciri3")
    assert status["status"] == "READY"
    assert status["details"]["jar"] == str(fake_jar)
    assert status["details"]["java"] == str(fake_java)

    doctor_res = runner.invoke(doctor_app, [])
    assert doctor_res.exit_code in {0, 1}
    assert str(fake_jar) in doctor_res.stdout

    detectors_res = runner.invoke(detectors_app, ["--json"])
    assert detectors_res.exit_code == 0
    payload = json.loads(detectors_res.stdout)
    ciri3 = next(row for row in payload["detectors"] if row["name"] == "ciri3")
    assert ciri3["status"] == "READY"
    assert ciri3["runtime"]["jar"] == str(fake_jar)
    assert ciri3["runtime"]["java"] == str(fake_java)


def test_doctor_and_detectors_share_resolution_from_non_repo_cwd(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    monkeypatch.chdir(tmp_path)
    jar_path = find_ciri_full_jar().resolved_path
    assert jar_path is not None

    doctor_res = runner.invoke(doctor_app, [])
    assert doctor_res.exit_code in {0, 1}
    assert str(jar_path) in doctor_res.stdout

    detectors_res = runner.invoke(detectors_app, ["--json"])
    assert detectors_res.exit_code == 0
    payload = json.loads(detectors_res.stdout)
    ciri_full = next(row for row in payload["detectors"] if row["name"] == "ciri-full")
    assert ciri_full["assets"]["CIRI-full-jar"]["path"] == str(jar_path)
    assert ciri_full["status"] in {"READY", "PARTIAL", "NOT READY"}


def test_ciri_full_runtime_resolves_paths_independent_of_cwd(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    detector = CiriFullDetector()
    captured: dict[str, object] = {}

    def fake_run(cmd, **kwargs):
        captured["cmd"] = cmd
        captured.update(kwargs)

        class _Result:
            returncode = 0

        return _Result()

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr("subprocess.run", fake_run)

    inputs = DetectorRunInputs(
        cell_id="sc01",
        outdir=tmp_path,
        r1=tmp_path / "sc01_R1.fastq.gz",
        r2=tmp_path / "sc01_R2.fastq.gz",
        ref_fa=tmp_path / "ref.fa",
        gtf=tmp_path / "genes.gtf",
        threads=1,
    )

    result = detector.run(inputs)

    jar_path = find_ciri_full_jar().resolved_path
    assert jar_path is not None
    assert captured["cmd"][0] == "bash"
    assert Path(captured["cmd"][1]).is_absolute()
    assert captured["env"]["CIRCYTO_CIRI_FULL_JAR"] == str(jar_path)
    assert result.meta["ciri_full_jar"] == str(jar_path)


def test_run_cirifull_manifest_resolves_relative_paths_from_manifest_dir(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    captured: list[str] = []
    manifest_dir = tmp_path / "manifest_dir"
    fastq_dir = manifest_dir / "reads"
    fastq_dir.mkdir(parents=True, exist_ok=True)
    (fastq_dir / "sc01_R1.fastq.gz").write_text("r1", encoding="utf-8")
    (fastq_dir / "sc01_R2.fastq.gz").write_text("r2", encoding="utf-8")

    manifest = manifest_dir / "manifest.tsv"
    manifest.write_text(
        "cell_id\tr1\tr2\n"
        "sc01\treads/sc01_R1.fastq.gz\treads/sc01_R2.fastq.gz\n",
        encoding="utf-8",
    )

    def fake_run_shell(cmd: str) -> None:
        captured.append(cmd)

    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir(parents=True, exist_ok=True)
    monkeypatch.chdir(elsewhere)
    monkeypatch.setattr("circyto.pipeline.run_cirifull._run_shell", fake_run_shell)

    run_cirifull_with_manifest(
        manifest=manifest,
        outdir=tmp_path / "out",
        cmd_template="echo {r1} {r2} > {out_tsv}",
        ref_fa=tmp_path / "ref.fa",
        gtf=tmp_path / "genes.gtf",
    )

    assert str((fastq_dir / "sc01_R1.fastq.gz").resolve()) in captured[0]
    assert str((fastq_dir / "sc01_R2.fastq.gz").resolve()) in captured[0]
