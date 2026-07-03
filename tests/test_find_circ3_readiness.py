from __future__ import annotations

from circyto.detectors.find_circ3 import FindCirc3Detector


def test_find_circ3_availability_uses_shutil_which(monkeypatch) -> None:
    present = {"bowtie2": "/usr/bin/bowtie2", "samtools": "/usr/bin/samtools"}

    def fake_which(name: str) -> str | None:
        return present.get(name)

    monkeypatch.setattr("circyto.detectors.find_circ3.shutil.which", fake_which)

    detector = FindCirc3Detector()
    assert detector.is_available() is False
    assert detector.missing_dependencies() == ["find-circ3"]

    present["find-circ3"] = "/usr/bin/find-circ3"
    assert detector.is_available() is True
    assert detector.missing_dependencies() == []
