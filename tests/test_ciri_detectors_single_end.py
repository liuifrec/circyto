from __future__ import annotations

import gzip
from pathlib import Path
from typing import Any

import pytest

from circyto.detectors.base import DetectorRunInputs
from circyto.detectors.ciri2 import Ciri2Detector
from circyto.detectors.ciri_full import CiriFullDetector
from circyto.paths import find_ciri2_adapter


def test_cirifull_single_end_sets_layout_and_fallback_metadata(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    det = CiriFullDetector()
    captured: dict[str, Any] = {}

    def fake_run(cmd, **kwargs):
        captured["cmd"] = cmd
        captured.update(kwargs)
        out_tsv = Path(kwargs["env"]["OUT_TSV"])
        out_tsv.write_text(
            "circ_id\tchr\tstart\tend\tstrand\tsupport\n"
            "circ1\tchr21\t1\t10\t+\t1\n",
            encoding="utf-8",
        )

        class _Result:
            returncode = 0

        return _Result()

    monkeypatch.setattr("subprocess.run", fake_run)

    inputs = DetectorRunInputs(
        cell_id="sc01",
        outdir=tmp_path,
        r1=tmp_path / "sc01_R1.fastq.gz",
        r2=None,
        ref_fa=tmp_path / "ref.fa",
        gtf=tmp_path / "genes.gtf",
        threads=1,
    )

    result = det.run(inputs)

    assert result.meta["read_layout"] == "single-end"
    assert result.meta["pipeline_mode"] == "ciri2-single-end-fallback"
    assert captured["env"]["CIRCYTO_READ_LAYOUT"] == "single"
    assert captured["env"]["R2"] == ""


def test_cirifull_paired_end_sets_layout_and_pipeline_metadata(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    det = CiriFullDetector()
    captured: dict[str, Any] = {}

    def fake_run(cmd, **kwargs):
        captured["cmd"] = cmd
        captured.update(kwargs)
        out_tsv = Path(kwargs["env"]["OUT_TSV"])
        out_tsv.write_text(
            "circ_id\tchr\tstart\tend\tstrand\tsupport\n"
            "circ1\tchr21\t1\t10\t+\t1\n",
            encoding="utf-8",
        )

        class _Result:
            returncode = 0

        return _Result()

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

    result = det.run(inputs)

    assert result.meta["read_layout"] == "paired-end"
    assert result.meta["pipeline_mode"] == "ciri-full-pipeline"
    assert captured["env"]["CIRCYTO_READ_LAYOUT"] == "paired"
    assert captured["env"]["R2"] == str(tmp_path / "sc01_R2.fastq.gz")


def test_cirifull_failure_includes_layout_command_and_log_tail(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    det = CiriFullDetector()

    def fake_run(cmd, **kwargs):
        log_handle = kwargs["stdout"]
        log_handle.write("first line\n")
        log_handle.write(">>> CMD[bwa]: bwa mem -T 19 ref.fa reads.fq\n")
        log_handle.write("fatal problem\n")

        class _Result:
            returncode = 2

        return _Result()

    monkeypatch.setattr("subprocess.run", fake_run)

    inputs = DetectorRunInputs(
        cell_id="sc01",
        outdir=tmp_path,
        r1=tmp_path / "sc01_R1.fastq.gz",
        r2=None,
        ref_fa=tmp_path / "ref.fa",
        gtf=tmp_path / "genes.gtf",
        threads=1,
    )

    with pytest.raises(RuntimeError) as excinfo:
        det.run(inputs)

    message = str(excinfo.value)
    assert "read_layout=single-end" in message
    assert "Adapter command: bash" in message
    assert "fatal problem" in message


def test_ciri2_uses_se_recommended_flags(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    det = Ciri2Detector()
    captured: dict[str, object] = {}

    def fake_run(cmd, **kwargs):
        captured["cmd"] = cmd
        captured.update(kwargs)

        class _Result:
            returncode = 0

        return _Result()

    monkeypatch.setattr("subprocess.run", fake_run)
    monkeypatch.chdir(tmp_path)

    inputs = DetectorRunInputs(
        cell_id="sc01",
        outdir=tmp_path,
        r1=tmp_path / "sc01_R1.fastq.gz",
        r2=None,
        ref_fa=tmp_path / "ref.fa",
        gtf=tmp_path / "genes.gtf",
        threads=1,
    )

    result = det.run(inputs)

    assert result.tsv_path == tmp_path / "sc01.tsv"
    adapter_path = find_ciri2_adapter().resolved_path
    assert adapter_path is not None
    assert captured["cmd"] == ["bash", str(adapter_path)]
    assert isinstance(captured["env"], dict)
    assert captured["check"] is False
    assert captured["stderr"] is not None
    assert captured["env"]["CIRI2_FLAGS"] == "-0 -U 15"
    assert captured["env"]["CIRI2_BWA_MEM_FLAGS"] == "-T 19"


def test_ciri2_uses_short_read_bwa_flags_for_short_se_fastq(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    det = Ciri2Detector()
    captured: dict[str, object] = {}

    def fake_run(cmd, **kwargs):
        captured["cmd"] = cmd
        captured.update(kwargs)

        class _Result:
            returncode = 0

        return _Result()

    monkeypatch.setattr("subprocess.run", fake_run)
    monkeypatch.chdir(tmp_path)

    r1 = tmp_path / "sc01_R1.fastq.gz"
    with gzip.open(r1, "wt", encoding="utf-8") as handle:
        handle.write("@r1\n" + ("A" * 43) + "\n+\n" + ("I" * 43) + "\n")

    inputs = DetectorRunInputs(
        cell_id="sc01",
        outdir=tmp_path,
        r1=r1,
        r2=None,
        ref_fa=tmp_path / "ref.fa",
        gtf=tmp_path / "genes.gtf",
        threads=1,
    )

    det.run(inputs)

    adapter_path = find_ciri2_adapter().resolved_path
    assert adapter_path is not None
    assert captured["cmd"] == ["bash", str(adapter_path)]
    assert isinstance(captured["env"], dict)
    assert captured["env"]["CIRI2_BWA_MEM_FLAGS"] == "-k 15 -T 15"
