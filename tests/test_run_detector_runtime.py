from __future__ import annotations

import json
from pathlib import Path

import pytest

from circyto.detectors.base import DetectorResult
from circyto.pipeline.run_detector import run_detector_manifest


class _FakeDetector:
    name = "ciri2"
    max_parallel = 4

    def __init__(self) -> None:
        self.calls: list[str] = []

    def run(self, inputs):
        self.calls.append(inputs.cell_id)
        out = Path(inputs.outdir) / f"{inputs.cell_id}.tsv"
        out.write_text(
            "circ_id\tchr\tstart\tend\tstrand\tsupport\n"
            f"{inputs.cell_id}\tchr1\t1\t10\t+\t1\n",
            encoding="utf-8",
        )
        return DetectorResult(
            detector=self.name,
            cell_id=inputs.cell_id,
            outdir=Path(inputs.outdir),
            tsv_path=out,
            meta={},
        )


def _write_manifest(tmp_path: Path, cells: list[str]) -> Path:
    manifest = tmp_path / "manifest.tsv"
    lines = ["cell_id\tr1\tr2"]
    for cell_id in cells:
        r1 = tmp_path / f"{cell_id}_R1.fastq.gz"
        r2 = tmp_path / f"{cell_id}_R2.fastq.gz"
        r1.write_text("r1", encoding="utf-8")
        r2.write_text("r2", encoding="utf-8")
        lines.append(f"{cell_id}\t{r1}\t{r2}")
    manifest.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return manifest


def test_run_detector_manifest_skips_completed_outputs_and_writes_summary(tmp_path: Path) -> None:
    manifest = _write_manifest(tmp_path, ["c1", "c2"])
    outdir = tmp_path / "run"
    outdir.mkdir(parents=True, exist_ok=True)
    existing_output = outdir / "c1.tsv"
    existing_output.write_text(
        "circ_id\tchr\tstart\tend\tstrand\tsupport\n"
        "circ1\tchr1\t1\t10\t+\t1\n",
        encoding="utf-8",
    )
    (outdir / "c1.tsv.provenance.json").write_text(
        json.dumps(
                {
                    "detector": "ciri2",
                    "detector_class": "_FakeDetector",
                    "cell_id": "c1",
                    "read1": str((tmp_path / "c1_R1.fastq.gz").resolve()),
                    "read2": str((tmp_path / "c1_R2.fastq.gz").resolve()),
                    "read_layout": "paired-end",
                    "ref_fa": str((tmp_path / "ref.fa").resolve()),
                    "gtf": str((tmp_path / "genes.gtf").resolve()),
                    "threads": 2,
                    "extra": {},
                }
        )
        + "\n",
        encoding="utf-8",
    )

    detector = _FakeDetector()
    results = run_detector_manifest(
        detector=detector,
        manifest=manifest,
        outdir=outdir,
        ref_fa=tmp_path / "ref.fa",
        gtf=tmp_path / "genes.gtf",
        threads=2,
        parallel=2,
    )

    assert detector.calls == ["c2"]
    assert len(results) == 2
    assert sum(1 for r in results if r.meta.get("skipped_existing")) == 1

    summary = json.loads((outdir / "detector_run_summary.json").read_text(encoding="utf-8"))
    assert summary["status_counts"]["skipped_existing"] == 1
    assert summary["status_counts"]["success"] == 1


def test_run_detector_manifest_failure_writes_summary(tmp_path: Path) -> None:
    manifest = _write_manifest(tmp_path, ["c1", "c2"])
    outdir = tmp_path / "run"

    class _FailingDetector(_FakeDetector):
        def run(self, inputs):
            if inputs.cell_id == "c2":
                raise RuntimeError("boom")
            return super().run(inputs)

    with pytest.raises(RuntimeError) as excinfo:
        run_detector_manifest(
            detector=_FailingDetector(),
            manifest=manifest,
            outdir=outdir,
            ref_fa=tmp_path / "ref.fa",
            gtf=tmp_path / "genes.gtf",
            threads=2,
            parallel=2,
        )

    assert "failed for 1/2 cells" in str(excinfo.value)
    summary = json.loads((outdir / "detector_run_summary.json").read_text(encoding="utf-8"))
    assert summary["status_counts"]["failed"] == 1
    assert summary["status_counts"]["success"] == 1
    failed_record = next(cell for cell in summary["cells"] if cell["cell_id"] == "c2")
    assert failed_record["read_layout"] == "paired-end"


def test_run_detector_manifest_reruns_when_provenance_changes(tmp_path: Path) -> None:
    manifest = _write_manifest(tmp_path, ["c1"])
    outdir = tmp_path / "run"
    detector = _FakeDetector()
    ref_a = tmp_path / "ref_a.fa"
    ref_b = tmp_path / "ref_b.fa"
    gtf = tmp_path / "genes.gtf"
    ref_a.write_text(">chr1\nACGT\n", encoding="utf-8")
    ref_b.write_text(">chr1\nTGCA\n", encoding="utf-8")
    gtf.write_text('chr1\ttest\texon\t1\t10\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n', encoding="utf-8")

    run_detector_manifest(
        detector=detector,
        manifest=manifest,
        outdir=outdir,
        ref_fa=ref_a,
        gtf=gtf,
        threads=2,
        parallel=1,
    )
    assert detector.calls == ["c1"]

    run_detector_manifest(
        detector=detector,
        manifest=manifest,
        outdir=outdir,
        ref_fa=ref_b,
        gtf=gtf,
        threads=2,
        parallel=1,
    )
    assert detector.calls == ["c1", "c1"]


def test_run_detector_manifest_records_single_end_layout_in_summary(tmp_path: Path) -> None:
    manifest = tmp_path / "manifest.tsv"
    r1 = tmp_path / "c1_R1.fastq.gz"
    r1.write_text("r1", encoding="utf-8")
    manifest.write_text(
        "cell_id\tr1\n"
        f"c1\t{r1}\n",
        encoding="utf-8",
    )

    detector = _FakeDetector()
    outdir = tmp_path / "run"
    run_detector_manifest(
        detector=detector,
        manifest=manifest,
        outdir=outdir,
        ref_fa=tmp_path / "ref.fa",
        gtf=tmp_path / "genes.gtf",
        threads=2,
        parallel=1,
    )

    summary = json.loads((outdir / "detector_run_summary.json").read_text(encoding="utf-8"))
    assert summary["cells"][0]["read_layout"] == "single-end"

    provenance = json.loads((outdir / "c1.tsv.provenance.json").read_text(encoding="utf-8"))
    assert provenance["read_layout"] == "single-end"
