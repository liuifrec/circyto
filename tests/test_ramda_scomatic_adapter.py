from __future__ import annotations

import json
import os
import subprocess
import textwrap
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "prepare_ramda_scomatic_input.sh"


def _write_fake_samtools(tmp_path: Path) -> Path:
    bindir = tmp_path / "bin"
    bindir.mkdir()
    samtools = bindir / "samtools"
    samtools.write_text(
        textwrap.dedent(
            """\
            #!/usr/bin/env python3
            from pathlib import Path
            import sys

            args = sys.argv[1:]
            if not args:
                raise SystemExit("missing command")
            cmd = args[0]
            rest = args[1:]

            if cmd == "view":
                if "-b" in rest:
                    sys.stdout.write(sys.stdin.read())
                else:
                    path = Path(rest[-1])
                    sys.stdout.write(path.read_text(encoding="utf-8"))
            elif cmd == "sort":
                out = rest[rest.index("-o") + 1]
                Path(out).write_text(sys.stdin.read(), encoding="utf-8")
            elif cmd == "merge":
                out = rest[rest.index("-f") + 1]
                start = rest.index("-f") + 2
                inputs = [item for item in rest[start:] if not item.startswith("-")]
                Path(out).write_text(
                    "".join(Path(item).read_text(encoding="utf-8") for item in inputs),
                    encoding="utf-8",
                )
            elif cmd == "index":
                Path(rest[0] + ".bai").write_text("index\\n", encoding="utf-8")
            else:
                raise SystemExit(f"unsupported fake samtools command: {cmd}")
            """
        ),
        encoding="utf-8",
    )
    samtools.chmod(0o755)
    return bindir


def _write_sam(path: Path, *, read_name: str, cb_tag: str = "") -> None:
    tag = f"\t{cb_tag}" if cb_tag else ""
    path.write_text(
        "@HD\tVN:1.6\tSO:unsorted\n"
        "@SQ\tSN:chr21\tLN:1000\n"
        f"{read_name}\t0\tchr21\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\tnM:i:0\tNH:i:1{tag}\n",
        encoding="utf-8",
    )


def _run_adapter(tmp_path: Path, manifest: Path, outdir: Path) -> subprocess.CompletedProcess[str]:
    env = os.environ.copy()
    env["PATH"] = f"{_write_fake_samtools(tmp_path)}{os.pathsep}{env['PATH']}"
    return subprocess.run(
        [
            "bash",
            str(SCRIPT),
            "--alignment-manifest",
            str(manifest),
            "--outdir",
            str(outdir),
            "--sample-id",
            "sample1",
        ],
        cwd=ROOT,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )


def test_prepare_ramda_scomatic_input_handles_one_cell_per_bam_inputs(tmp_path: Path) -> None:
    cell_a = tmp_path / "cellA.bam"
    cell_b = tmp_path / "cellB.bam"
    _write_sam(cell_a, read_name="readA")
    _write_sam(cell_b, read_name="readB", cb_tag="CB:Z:cellB")
    manifest = tmp_path / "alignment_manifest.tsv"
    manifest.write_text(
        "cell_id\tbam\tsam\tgroup_id\tread_layout\taligner\treference\tcache_key\tsource_manifest\n"
        f"cellA\t{cell_a}\t\tgrp\tsingle-end\tbwa-mem\t/tmp/ref.fa\tk1\t/tmp/source.tsv\n"
        f"cellB\t{cell_b}\t\tgrp\tsingle-end\tbwa-mem\t/tmp/ref.fa\tk2\t/tmp/source.tsv\n",
        encoding="utf-8",
    )
    outdir = tmp_path / "adapter"

    proc = _run_adapter(tmp_path, manifest, outdir)

    assert proc.returncode == 0, proc.stderr
    merged = outdir / "merged" / "sample1.scomatic.bam"
    assert merged.exists()
    assert Path(str(merged) + ".bai").exists()
    assert (outdir / "per_cell" / "cellA.scomatic.sorted.bam").exists()
    assert (outdir / "per_cell" / "cellB.scomatic.sorted.bam").exists()
    summary = json.loads((outdir / "adapter_summary.json").read_text(encoding="utf-8"))
    assert summary["input_cells"] == 2
    assert summary["runs_scomatic"] is False
    assert summary["scomatic_result_terminology"] == "RNA-derived candidate variant signals"


def test_prepare_ramda_scomatic_input_injects_missing_cb_tags(tmp_path: Path) -> None:
    cell_a = tmp_path / "cellA.sam"
    _write_sam(cell_a, read_name="readA")
    manifest = tmp_path / "alignment_manifest.tsv"
    manifest.write_text(
        "cell_id\tbam\tsam\tgroup_id\tread_layout\taligner\treference\tcache_key\tsource_manifest\n"
        f"cellA\t\t{cell_a}\tgrp\tsingle-end\tbwa-mem\t/tmp/ref.fa\tk1\t/tmp/source.tsv\n",
        encoding="utf-8",
    )
    outdir = tmp_path / "adapter"

    proc = _run_adapter(tmp_path, manifest, outdir)

    assert proc.returncode == 0, proc.stderr
    prepared = (outdir / "per_cell" / "cellA.scomatic.sorted.bam").read_text(encoding="utf-8")
    assert "CB:Z:cellA" in prepared
    summary = json.loads((outdir / "adapter_summary.json").read_text(encoding="utf-8"))
    assert summary["tag_summary"]["missing_cb_before_injection"] == 1
    assert summary["tag_summary"]["missing_nM"] == 0
    assert summary["tag_summary"]["missing_NH"] == 0


def test_prepare_ramda_scomatic_input_normalizes_existing_cb_to_cell_id(tmp_path: Path) -> None:
    cell_a = tmp_path / "cellA.sam"
    _write_sam(cell_a, read_name="readA", cb_tag="CB:Z:oldBarcode")
    manifest = tmp_path / "alignment_manifest.tsv"
    manifest.write_text(
        "cell_id\tbam\tsam\tgroup_id\tread_layout\taligner\treference\tcache_key\tsource_manifest\n"
        f"cellA\t\t{cell_a}\tgrp\tsingle-end\tbwa-mem\t/tmp/ref.fa\tk1\t/tmp/source.tsv\n",
        encoding="utf-8",
    )
    outdir = tmp_path / "adapter"

    proc = _run_adapter(tmp_path, manifest, outdir)

    assert proc.returncode == 0, proc.stderr
    prepared = (outdir / "per_cell" / "cellA.scomatic.sorted.bam").read_text(encoding="utf-8")
    assert "CB:Z:cellA" in prepared
    assert "CB:Z:oldBarcode" not in prepared
    summary = json.loads((outdir / "adapter_summary.json").read_text(encoding="utf-8"))
    assert summary["tag_summary"]["missing_cb_before_injection"] == 0
    assert summary["tag_summary"]["replaced_cb_not_matching_cell_id"] == 1


def test_prepare_ramda_scomatic_input_writes_annotation_table(tmp_path: Path) -> None:
    cell_a = tmp_path / "cellA.sam"
    _write_sam(cell_a, read_name="readA")
    manifest = tmp_path / "alignment_manifest.tsv"
    manifest.write_text(
        "cell_id\tbam\tsam\tgroup_id\tread_layout\taligner\treference\tcache_key\tsource_manifest\n"
        f"cellA\t\t{cell_a}\tgrp\tsingle-end\tbwa-mem\t/tmp/ref.fa\tk1\t/tmp/source.tsv\n",
        encoding="utf-8",
    )
    outdir = tmp_path / "adapter"

    proc = _run_adapter(tmp_path, manifest, outdir)

    assert proc.returncode == 0, proc.stderr
    annotation = (outdir / "cell_annotations.tsv").read_text(encoding="utf-8").splitlines()
    assert annotation == ["Index\tCell_type", "cellA\tcellA"]


def test_prepare_ramda_scomatic_input_fails_on_missing_bam(tmp_path: Path) -> None:
    manifest = tmp_path / "alignment_manifest.tsv"
    missing = tmp_path / "missing.bam"
    manifest.write_text(
        "cell_id\tbam\tsam\tgroup_id\tread_layout\taligner\treference\tcache_key\tsource_manifest\n"
        f"cellA\t{missing}\t\tgrp\tsingle-end\tbwa-mem\t/tmp/ref.fa\tk1\t/tmp/source.tsv\n",
        encoding="utf-8",
    )
    outdir = tmp_path / "adapter"

    proc = _run_adapter(tmp_path, manifest, outdir)

    assert proc.returncode != 0
    assert "alignment file not found for cell_id=cellA" in proc.stderr
