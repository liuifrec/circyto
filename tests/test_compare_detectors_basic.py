# tests/test_compare_detectors_basic.py
from pathlib import Path

from circyto.compare import compare_detectors_from_root


def _write_tsv(path: Path, circ_ids):
    path.write_text(
        "circ_id\tsupport\n"
        + "\n".join(f"{cid}\t1" for cid in circ_ids)
        + "\n"
    )


def test_compare_detectors_union_intersection_and_matrix(tmp_path):
    root = tmp_path / "root"
    detA = root / "detA"
    detB = root / "detB"
    detA.mkdir(parents=True)
    detB.mkdir(parents=True)

    # detA:
    #  cell1: circ1, circ2
    #  cell2: circ2, circ3
    _write_tsv(detA / "cell1.tsv", ["circ1", "circ2"])
    _write_tsv(detA / "cell2.tsv", ["circ2", "circ3"])

    # detB:
    #  cell1: circ2
    #  cell2: circ3, circ4
    _write_tsv(detB / "cell1.tsv", ["circ2"])
    _write_tsv(detB / "cell2.tsv", ["circ3", "circ4"])

    out = tmp_path / "out"

    summary = compare_detectors_from_root(root, ["detA", "detB"], out)

    # union = {circ1, circ2, circ3, circ4} → 4
    assert summary["union_size"] == 4

    # intersection = {circ2, circ3} → 2
    assert summary["intersection_size"] == 2

    # per-detector counts
    dets = summary["detectors"]
    assert dets["detA"]["n_circs"] == 3  # circ1,2,3
    assert dets["detB"]["n_circs"] == 3  # circ2,3,4

    # Check files exist
    assert (out / "circ_detector.mtx").exists()
    assert (out / "circ_ids.txt").exists()
    assert (out / "detectors.txt").exists()
    assert (out / "union.tsv").exists()
    assert (out / "intersection.tsv").exists()
    assert (out / "summary.json").exists()
