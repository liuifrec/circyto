# tests/test_compare_detectors_normalization.py
from pathlib import Path

from circyto.compare.detector_compare import _load_circ_ids_from_tsv


def test_load_circ_ids_prefers_circ_id_column(tmp_path):
    tsv = tmp_path / "det.tsv"
    tsv.write_text(
        "circ_id\tsupport\n"
        "circA\t1\n"
        "circB\t2\n"
    )

    circs = _load_circ_ids_from_tsv(tsv)
    assert circs == {"circA", "circB"}


def test_load_circ_ids_coord_mode(tmp_path):
    tsv = tmp_path / "det.tsv"
    tsv.write_text(
        "chr\tstart\tend\tstrand\tsupport\n"
        "chr1\t100\t200\t+\t1\n"
        "chr1\t300\t400\t-\t1\n"
    )

    circs = _load_circ_ids_from_tsv(tsv)
    assert "chr1:100|200|+" in circs
    assert "chr1:300|400|-" in circs


def test_load_circ_ids_fallback_first_column(tmp_path):
    tsv = tmp_path / "det.tsv"
    tsv.write_text(
        "weird_name\tsupport\n"
        "X1\t1\n"
        "X2\t1\n"
    )

    circs = _load_circ_ids_from_tsv(tsv)
    assert circs == {"X1", "X2"}
