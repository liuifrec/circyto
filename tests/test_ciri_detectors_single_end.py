from __future__ import annotations

from pathlib import Path

import pytest

from circyto.detectors.base import DetectorRunInputs
from circyto.detectors.ciri2 import Ciri2Detector
from circyto.detectors.ciri_full import CiriFullDetector


def test_cirifull_rejects_single_end_input(tmp_path: Path) -> None:
    det = CiriFullDetector()
    inputs = DetectorRunInputs(
        cell_id="sc01",
        outdir=tmp_path,
        r1=tmp_path / "sc01_R1.fastq.gz",
        r2=None,
        ref_fa=tmp_path / "ref.fa",
        gtf=tmp_path / "genes.gtf",
        threads=1,
    )

    with pytest.raises(ValueError, match="paired-end FASTQ input"):
        det.run(inputs)


def test_ciri2_uses_se_recommended_flags(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    det = Ciri2Detector()
    captured: dict[str, object] = {}

    def fake_run(cmd, shell=False, env=None):
        captured["cmd"] = cmd
        captured["shell"] = shell
        captured["env"] = env

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

    assert result.tsv_path == tmp_path / "sc01.tsv"
    assert captured["shell"] is True
    assert isinstance(captured["env"], dict)
    assert captured["env"]["CIRI2_FLAGS"] == "-0 -U 15"
