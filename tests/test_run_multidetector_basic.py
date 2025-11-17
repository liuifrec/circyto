# tests/test_run_multidetector_basic.py
from pathlib import Path
from dataclasses import dataclass
from typing import Optional, Dict, Any

from circyto.pipeline.run_detector import run_multidetector
from circyto.detectors.base import DetectorBase, DetectorRunInputs, DetectorResult


@dataclass
class FakeDetector:
    name: str = "fake-det"
    input_type: str = "fastq"
    supports_paired_end: bool = True

    calls: int = 0

    def is_available(self) -> bool:
        return True

    def version(self) -> Optional[str]:
        return "fake-1.0"

    def run(self, inputs: DetectorRunInputs) -> DetectorResult:
        self.calls += 1
        outdir = inputs.outdir
        outdir.mkdir(parents=True, exist_ok=True)
        tsv_path = outdir / f"{inputs.cell_id}.tsv"
        # write a tiny TSV
        tsv_path.write_text("circ_id\tsupport\nfake_circ\t1\n")
        return DetectorResult(
            detector=self.name,
            cell_id=inputs.cell_id,
            outdir=outdir,
            tsv_path=tsv_path,
            meta={"fake": True},
        )


def test_run_multidetector_creates_per_detector_dirs(tmp_path):
    # manifest with two cells
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        "cell_id\tr1\tr2\n"
        "c1\tR1_1.fq\tR1_2.fq\n"
        "c2\tR2_1.fq\t\n"
    )

    root_out = tmp_path / "multi"
    detA = FakeDetector(name="detA")
    detB = FakeDetector(name="detB")

    detectors = {"detA": detA, "detB": detB}

    res = run_multidetector(
        detectors=detectors,
        manifest=manifest,
        root_outdir=root_out,
        ref_fa=Path("ref.fa"),
        gtf=Path("genes.gtf"),
        threads=2,
        parallel=1,
    )

    # Ensure keys
    assert set(res.keys()) == {"detA", "detB"}
    # Each detector should have 2 results (two cells)
    assert len(res["detA"]) == 2
    assert len(res["detB"]) == 2

    # Directories and TSV files exist
    for det_name in ["detA", "detB"]:
        det_dir = root_out / det_name
        assert det_dir.exists()
        assert (det_dir / "c1.tsv").exists()
        assert (det_dir / "c2.tsv").exists()

    # The fake detectors' run methods were called
    assert detA.calls == 2
    assert detB.calls == 2
