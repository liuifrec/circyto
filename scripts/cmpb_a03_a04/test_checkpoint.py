import csv
import json

import pytest

from robustness import archive, checkpoint_completed_cases


def test_checkpoint_skips_only_complete_retained_cases(tmp_path):
    records = tmp_path / "records"
    records.mkdir()
    (records / "storage_amendment.json").write_text("{}")
    folder = tmp_path / "a04"
    directory = folder / "cases/complete"
    directory.mkdir(parents=True)
    (directory / "data.txt").write_text("generated fixture")
    retained = archive(directory)
    with (folder / "archives.tsv").open("w") as f:
        w = csv.DictWriter(f, fieldnames=["case", "archive", "sha256"], delimiter="\t")
        w.writeheader()
        w.writerow(dict(case="complete", archive="a04/cases/complete.tar.gz", sha256=retained["sha256"]))
    (folder / "robustness_outcomes.tsv").write_text(
        "case\tarm\ncomplete\tcircyto\ncomplete\treference\ncomplete\tunsafe\npartial\tcircyto\n"
    )
    assert checkpoint_completed_cases(tmp_path) == {"complete"}
    (records / "a04_stop.json").write_text("{}")
    with pytest.raises(RuntimeError, match="Scientific stop"):
        checkpoint_completed_cases(tmp_path)


def test_checkpoint_reads_each_segment_and_blocks_its_scientific_stop(tmp_path):
    records = tmp_path / "records"
    records.mkdir()
    (records / "storage_amendment.json").write_text("{}")
    folder = tmp_path / "a04_continued"
    directory = folder / "cases/later"
    directory.mkdir(parents=True)
    (directory / "data.txt").write_text("generated fixture")
    retained = archive(directory)
    (folder / "archives.tsv").write_text("case\tarchive\tsha256\nlater\ta04_continued/cases/later.tar.gz\t"+retained["sha256"]+"\n")
    (folder / "robustness_outcomes.tsv").write_text("case\tarm\nlater\tcircyto\nlater\treference\nlater\tunsafe\n")
    assert checkpoint_completed_cases(tmp_path, "a04_continued") == {"later"}
    (records / "a04_continued_stop.json").write_text("{}")
    with pytest.raises(RuntimeError, match="Scientific stop"):
        checkpoint_completed_cases(tmp_path, "a04_continued")
