from __future__ import annotations

import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def _run_script(manifest: Path, dataset_root: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["bash", str(ROOT / "scripts" / "download_fastqs_from_manifest.sh"), str(manifest), str(dataset_root)],
        capture_output=True,
        text=True,
        cwd=ROOT,
    )


def test_download_manifest_parser_accepts_single_end_empty_fastq2_and_extra_columns(tmp_path: Path) -> None:
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        "sample_id\tfastq_1\tfastq_2\tprotocol\tstrandedness\tread_layout\tsrr\tgsm\n"
        "GSM8558868\traw/SRR30918094.fastq.gz\t\tramda\tunstranded\tsingle\tSRR30918094\tGSM8558868\n",
        encoding="utf-8",
    )
    dataset_root = tmp_path / "dataset"
    raw_dir = dataset_root / "raw"
    raw_dir.mkdir(parents=True)
    (raw_dir / "SRR30918094.fastq.gz").write_text("", encoding="utf-8")

    proc = _run_script(manifest, dataset_root)

    assert proc.returncode == 0
    assert "parsed manifest row sample_id=GSM8558868 read_layout=single" in proc.stdout
    assert "single-end FASTQ already present, skipping" in proc.stdout
    assert "unsupported read_layout" not in proc.stdout + proc.stderr


def test_download_manifest_parser_accepts_paired_end_row(tmp_path: Path) -> None:
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        "sample_id\tfastq_1\tfastq_2\tprotocol\tstrandedness\tread_layout\textra\n"
        "sample1\traw/SRR123_1.fastq.gz\traw/SRR123_2.fastq.gz\tramda\tunstranded\tpaired\tmeta\n",
        encoding="utf-8",
    )
    dataset_root = tmp_path / "dataset"
    raw_dir = dataset_root / "raw"
    raw_dir.mkdir(parents=True)
    (raw_dir / "SRR123_1.fastq.gz").write_text("", encoding="utf-8")
    (raw_dir / "SRR123_2.fastq.gz").write_text("", encoding="utf-8")

    proc = _run_script(manifest, dataset_root)

    assert proc.returncode == 0
    assert "parsed manifest row sample_id=sample1 read_layout=paired" in proc.stdout
    assert "paired-end FASTQs already present, skipping" in proc.stdout


def test_download_manifest_parser_fails_on_missing_required_columns(tmp_path: Path) -> None:
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        "sample_id\tfastq_1\tread_layout\n"
        "sample1\traw/SRR1.fastq.gz\tsingle\n",
        encoding="utf-8",
    )
    dataset_root = tmp_path / "dataset"

    proc = _run_script(manifest, dataset_root)

    assert proc.returncode != 0
    assert "manifest must contain header columns" in proc.stdout + proc.stderr
