from __future__ import annotations

import json
from pathlib import Path

import pytest
from scipy.io import mmread

from circyto.detectors.ciri3 import Ciri3Detector
from circyto.manifest.alignment import AlignmentManifestRow, write_alignment_manifest_tsv
from circyto.paths import Ciri3Resolution, PathResolution
from circyto.pipeline.align_manifest import run_detector_alignment_manifest
from circyto.cli.circyto import app
from typer.testing import CliRunner


runner = CliRunner()


def _fake_ciri3_resolution(tmp_path: Path) -> Ciri3Resolution:
    fake_home = tmp_path / "ciri3_home"
    fake_home.mkdir(parents=True, exist_ok=True)
    fake_jar = fake_home / "CIRI3_Java.jar"
    fake_jar.write_text("jar", encoding="utf-8")
    fake_java = tmp_path / "java"
    fake_java.write_text("java", encoding="utf-8")
    return Ciri3Resolution(
        home=PathResolution("home", fake_home, (fake_home,), "test"),
        jar=PathResolution("jar", fake_jar, (fake_jar,), "test"),
        bin=PathResolution("bin", None, tuple(), None),
        java=PathResolution("java", fake_java, (fake_java,), "test"),
    )


def _write_alignment_manifest(tmp_path: Path, *, read_layout: str) -> tuple[Path, Path]:
    sam = tmp_path / f"{read_layout}.sam"
    sam.write_text("@SQ\tSN:chr1\tLN:4\n", encoding="utf-8")
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n", encoding="utf-8")
    manifest = tmp_path / "alignment_manifest.tsv"
    write_alignment_manifest_tsv(
        [
            AlignmentManifestRow(
                cell_id=f"cell_{read_layout}",
                sam=str(sam),
                group_id="grp",
                read_layout=read_layout,
                aligner="bwa-mem",
                reference=str(ref),
                cache_key="cache1",
                source_manifest=str(tmp_path / "source_manifest.tsv"),
                mapper_mode="0",
                artifact_bucket="bwa_mem",
                sortedness="unsorted",
            )
        ],
        manifest,
    )
    return manifest, ref


@pytest.mark.parametrize("read_layout", ["single-end", "paired-end"])
def test_ciri3_alignment_manifest_handles_layouts(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    read_layout: str,
) -> None:
    manifest, ref = _write_alignment_manifest(tmp_path, read_layout=read_layout)
    detector = Ciri3Detector()

    def fake_validate_runtime(self, *, template=None):
        return True, [], {
            "preferred_mode": "direct-jar",
            "direct_ready": True,
            "jar": str(tmp_path / "fake.jar"),
            "bin": None,
            "java": str(tmp_path / "java"),
            "home": str(tmp_path),
        }

    def fake_subprocess_run(cmd, **kwargs):
        raw_output = Path(cmd[cmd.index("-O") + 1])
        internal_log = Path(cmd[cmd.index("-G") + 1])
        raw_output.parent.mkdir(parents=True, exist_ok=True)
        raw_output.write_text(
            "circRNA_ID\tchr\tcircRNA_start\tcircRNA_end\tstrand\tjunction_reads\n"
            "circ1\tchr1\t1\t10\t+\t2\n",
            encoding="utf-8",
        )
        internal_log.write_text("ok\n", encoding="utf-8")

        class _Result:
            returncode = 0

        return _Result()

    monkeypatch.setattr(Ciri3Detector, "validate_runtime", fake_validate_runtime)
    monkeypatch.setattr("circyto.detectors.ciri3.resolve_ciri3_installation", lambda: _fake_ciri3_resolution(tmp_path))
    monkeypatch.setattr("subprocess.run", fake_subprocess_run)

    outdir = tmp_path / "out"
    results = run_detector_alignment_manifest(
        detector=detector,
        manifest=manifest,
        outdir=outdir,
        ref_fa=ref,
        threads=2,
        parallel=1,
    )

    assert len(results) == 1
    summary = json.loads((outdir / "detector_run_summary.json").read_text(encoding="utf-8"))
    assert summary["cells"][0]["read_layout"] == read_layout
    assert summary["cells"][0]["normalized_row_count"] == 1
    assert (outdir / f"cell_{read_layout}.tsv").exists()


def test_collect_matrix_cli_accepts_ciri3_outputs(tmp_path: Path) -> None:
    indir = tmp_path / "ciri3_out"
    indir.mkdir(parents=True, exist_ok=True)
    (indir / "cell_a.tsv").write_text(
        "circ_id\tchr\tstart\tend\tstrand\tsupport\n"
        "circ1\tchr1\t1\t10\t+\t1\n"
        "circ2\tchr1\t20\t30\t-\t2\n",
        encoding="utf-8",
    )
    (indir / "cell_b.tsv").write_text(
        "circ_id\tchr\tstart\tend\tstrand\tsupport\n"
        "circ1\tchr1\t1\t10\t+\t3\n",
        encoding="utf-8",
    )

    outdir = tmp_path / "matrix"
    result = runner.invoke(
        app,
        [
            "collect-matrix",
            "--detector",
            "ciri3",
            "--indir",
            str(indir),
            "--outdir",
            str(outdir),
        ],
    )

    assert result.exit_code == 0, result.output
    mat = mmread(outdir / "circ_counts.mtx").tocsr()
    assert mat.shape == (2, 2)
    assert mat.nnz == 3
