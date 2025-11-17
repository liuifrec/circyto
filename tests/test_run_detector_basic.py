# tests/test_run_detector_basic.py
from pathlib import Path
from circyto.pipeline.run_detector import read_manifest


def test_read_manifest_basic(tmp_path):
    m = tmp_path / "manifest.tsv"
    m.write_text(
        "cell_id\tr1\tr2\n"
        "c1\tR1_1.fq\tR1_2.fq\n"
        "c2\tR2_1.fq\t\n"
    )

    rows = read_manifest(m)
    assert len(rows) == 2
    assert rows[0][0] == "c1"
    assert rows[0][1] == Path("R1_1.fq")
    assert rows[0][2] == Path("R1_2.fq")

    assert rows[1][0] == "c2"
    assert rows[1][2] is None
