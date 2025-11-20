# tests/test_detectors_base.py
from pathlib import Path

from circyto.detectors import (
    available_detectors,
    DetectorRunInputs,
    DetectorResult,
    CiriFullDetector,
)


def test_detector_runinputs_basic():
    inputs = DetectorRunInputs(
        cell_id="cell1",
        r1=Path("R1.fastq.gz"),
        r2=Path("R2.fastq.gz"),
        outdir=Path("outdir"),
        ref_fa=Path("ref.fa"),
        gtf=Path("genes.gtf"),
        threads=4,
    )

    assert inputs.cell_id == "cell1"
    assert inputs.r1 == Path("R1.fastq.gz")
    assert inputs.r2 == Path("R2.fastq.gz")
    assert inputs.outdir == Path("outdir")
    assert inputs.ref_fa == Path("ref.fa")
    assert inputs.gtf == Path("genes.gtf")
    assert inputs.threads == 4
    assert isinstance(inputs.extra, dict)


def test_detectorresult_basic():
    res = DetectorResult(
        detector="dummy-det",
        cell_id="cell1",
        outdir=Path("outdir"),
        tsv_path=Path("outdir/cell1.tsv"),
    )

    assert res.detector == "dummy-det"
    assert res.cell_id == "cell1"
    assert res.outdir == Path("outdir")
    assert res.tsv_path == Path("outdir/cell1.tsv")
    assert res.run_dir is None or isinstance(res.run_dir, Path) or res.run_dir is None
    assert isinstance(res.meta, dict)


def test_available_detectors_does_not_crash():
    engines = available_detectors()
    assert isinstance(engines, dict)

    # If any engines are available (e.g. CIRI-full in a configured env),
    # they should have the expected core attributes.
    for name, eng in engines.items():
        assert isinstance(name, str)
        assert hasattr(eng, "name")
        assert hasattr(eng, "input_type")
        assert hasattr(eng, "supports_paired_end")
        assert callable(eng.is_available)
        assert callable(eng.version)
        assert callable(eng.run)


def test_cirifull_detector_basic_interface():
    """
    This is a shallow test that doesn't actually run the detector.
    It just checks that the class can be instantiated and reports
    availability in a non-crashing way.
    """
    det = CiriFullDetector()
    assert det.name == "ciri-full"
    # is_available() may return True or False depending on env; only require no crash
    available = det.is_available()
    assert isinstance(available, bool)

    ver = det.version()
    assert ver is None or isinstance(ver, str)
