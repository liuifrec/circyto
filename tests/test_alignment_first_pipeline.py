from __future__ import annotations

import json
import os
from pathlib import Path

from typer.testing import CliRunner

from circyto.cli.circyto import app
from circyto.detectors import build_default_engines, get_detector_capabilities
from circyto.detectors.ciri3 import Ciri3Detector, _normalize_ciri3_output, validate_ciri3_template
from circyto.detectors.base import DetectorCapabilities, DetectorResult, DetectorRunInputs
from circyto.manifest.alignment import read_alignment_manifest_tsv
from circyto.pipeline.align_manifest import (
    export_manifest_subset,
    plan_alignment_cache,
    prepare_alignment_cache,
    run_detector_alignment_manifest,
    summarize_run_state,
)


class _FakeAlignmentDetector:
    name = "fake-align-detector"
    input_type = "alignment"
    supports_paired_end = True
    max_parallel = 4
    capabilities = DetectorCapabilities(
        accepts_fastq=False,
        accepts_alignment=True,
        prefers_paired=True,
        supports_single_end=True,
        supports_multisample_alignment=True,
        max_parallel=4,
        recommended_execution_mode="alignment-first",
    )

    def run_from_alignment(self, inputs):
        assert inputs.input_mode == "alignment"
        assert inputs.bam is not None
        out = Path(inputs.outdir) / f"{inputs.cell_id}.tsv"
        out.write_text(
            "circ_id\tchr\tstart\tend\tstrand\tsupport\n"
            f"{inputs.cell_id}\tchr1\t10\t20\t+\t2\n",
            encoding="utf-8",
        )
        return DetectorResult(
            detector=self.name,
            cell_id=inputs.cell_id,
            outdir=Path(inputs.outdir),
            tsv_path=out,
            meta={"alignment_group": inputs.alignment_group},
        )

    def run(self, inputs):
        raise AssertionError("alignment-first test should not call generic run()")


def _write_source_manifest(tmp_path: Path) -> Path:
    bam1 = tmp_path / "c1.bam"
    bam2 = tmp_path / "c2.bam"
    bam1.write_text("bam1", encoding="utf-8")
    bam2.write_text("bam2", encoding="utf-8")
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        "\n".join(
            [
                "cell_id\tplatform\tread1\tread2\tbam\tlibrary_id\tn_input_reads\tgroup_id",
                f"c1\tsmartseq2\t\t\t{bam1}\tlib1\t100\tplateA",
                f"c2\tsmartseq2\t\t\t{bam2}\tlib1\t200\tplateA",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return manifest


def _write_fastq_manifest(tmp_path: Path, cells: list[str]) -> Path:
    lines = ["cell_id\tplatform\tread1\tread2\tbam\tlibrary_id\tn_input_reads\tgroup_id"]
    for cell in cells:
        r1 = tmp_path / f"{cell}_R1.fastq"
        r2 = tmp_path / f"{cell}_R2.fastq"
        r1.write_text("@r1\nACGT\n+\n!!!!\n", encoding="utf-8")
        r2.write_text("@r2\nTGCA\n+\n!!!!\n", encoding="utf-8")
        lines.append(f"{cell}\tsmartseq2\t{r1}\t{r2}\t\tlib1\t100\tplateA")
    manifest = tmp_path / "fastq_manifest.tsv"
    manifest.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return manifest


def _install_fake_bwa_and_samtools(bin_dir: Path, *, fail_on: str | None = None) -> None:
    bwa = bin_dir / "bwa"
    bwa.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        "if [[ -n \"${BWA_ARGS_LOG:-}\" ]]; then\n"
        "  printf '%s\\n' \"$*\" >> \"$BWA_ARGS_LOG\"\n"
        "fi\n"
        "args=(\"$@\")\n"
        "for a in \"${args[@]}\"; do\n"
        f"  if [[ -n \"{fail_on or ''}\" && \"$a\" == *\"{fail_on or ''}\"* ]]; then exit 19; fi\n"
        "done\n"
        "printf '@HD\\tVN:1.0\\n'\n"
        "printf 'r001\\t0\\tchr1\\t7\\t60\\t8M\\t*\\t0\\t0\\tACGTACGT\\tFFFFFFFF\\n'\n",
        encoding="utf-8",
    )
    bwa.chmod(0o755)

    samtools = bin_dir / "samtools"
    samtools.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        "cmd=\"$1\"; shift\n"
        "if [[ \"$cmd\" == \"sort\" ]]; then\n"
        "  out=''\n"
        "  while [[ $# -gt 0 ]]; do\n"
        "    case \"$1\" in\n"
        "      -o) out=\"$2\"; shift 2 ;;\n"
        "      -@) shift 2 ;;\n"
        "      *) shift ;;\n"
        "    esac\n"
        "  done\n"
        "  cat > \"$out\"\n"
        "elif [[ \"$cmd\" == \"index\" ]]; then\n"
        "  touch \"$1.bai\"\n"
        "else\n"
        "  exit 2\n"
        "fi\n",
        encoding="utf-8",
    )
    samtools.chmod(0o755)


def _install_fake_star(bin_dir: Path) -> None:
    star = bin_dir / "STAR"
    star.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        "prefix=''\n"
        "while [[ $# -gt 0 ]]; do\n"
        "  case \"$1\" in\n"
        "    --outFileNamePrefix) prefix=\"$2\"; shift 2 ;;\n"
        "    *) shift ;;\n"
        "  esac\n"
        "done\n"
        "mkdir -p \"$(dirname \"$prefix\")\"\n"
        "printf '@HD\\tVN:1.0\\n' > \"${prefix}Aligned.out.sam\"\n"
        "printf 'chr1\\t10\\t20\\n' > \"${prefix}Chimeric.out.junction\"\n"
        "printf '@mate1\\nACGT\\n+\\n!!!!\\n' > \"${prefix}Unmapped.out.mate1\"\n"
        "printf '@mate2\\nTGCA\\n+\\n!!!!\\n' > \"${prefix}Unmapped.out.mate2\"\n",
        encoding="utf-8",
    )
    star.chmod(0o755)


def test_prepare_alignment_cache_reuses_bam_inputs(tmp_path: Path) -> None:
    manifest = _write_source_manifest(tmp_path)
    outdir = tmp_path / "align"

    alignment_manifest = prepare_alignment_cache(
        manifest=manifest,
        outdir=outdir,
        aligner="reuse-existing",
        detector_hint="fake-align-detector",
        threads=2,
        parallel=2,
        sentinel_cells=1,
        chunk_size=1,
    )

    rows = read_alignment_manifest_tsv(alignment_manifest, validate_files=True)
    assert [row.cell_id for row in rows] == ["c1", "c2"]
    assert all(row.bam for row in rows)
    assert rows[0].group_id == "plateA"
    assert rows[0].extra["reused_alignment"] == "true"
    assert rows[0].sortedness == "sorted"

    summary = json.loads((outdir / "alignment_prepare_summary.json").read_text(encoding="utf-8"))
    assert summary["status_counts"]["reused_input"] == 2
    assert summary["sentinel_cells"] == 1
    assert summary["chunk_size"] == 1


def test_run_detector_alignment_manifest_records_alignment_provenance(tmp_path: Path) -> None:
    manifest = _write_source_manifest(tmp_path)
    alignment_dir = tmp_path / "align"
    alignment_manifest = prepare_alignment_cache(
        manifest=manifest,
        outdir=alignment_dir,
        aligner="reuse-existing",
        detector_hint="fake-align-detector",
        threads=2,
        parallel=2,
    )

    outdir = tmp_path / "detector"
    ref = tmp_path / "ref.fa"
    gtf = tmp_path / "genes.gtf"
    ref.write_text(">chr1\nACGT\n", encoding="utf-8")
    gtf.write_text('chr1\ttest\texon\t1\t10\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n', encoding="utf-8")
    results = run_detector_alignment_manifest(
        detector=_FakeAlignmentDetector(),
        manifest=alignment_manifest,
        outdir=outdir,
        ref_fa=ref,
        gtf=gtf,
        threads=2,
        parallel=2,
    )

    assert len(results) == 2
    summary = json.loads((outdir / "detector_run_summary.json").read_text(encoding="utf-8"))
    assert summary["input_mode"] == "alignment"
    assert summary["execution_mode"] == "alignment-first"
    assert summary["status_counts"]["success"] == 2
    assert summary["cells"][0]["reused_alignment"] is True
    assert summary["cells"][0]["detector_backend"] == "fake-align-detector"


def test_run_detector_alignment_manifest_prefers_detector_input_metadata(tmp_path: Path) -> None:
    class _MetaAwareDetector(_FakeAlignmentDetector):
        def run_from_alignment(self, inputs):
            result = super().run_from_alignment(inputs)
            result.meta.update(
                {
                    "input_file_type": "STAR tuple",
                    "input_sortedness": "unsorted",
                    "mapper_mode": "1",
                }
            )
            return result

    manifest = _write_source_manifest(tmp_path)
    alignment_dir = tmp_path / "align"
    alignment_manifest = prepare_alignment_cache(
        manifest=manifest,
        outdir=alignment_dir,
        aligner="reuse-existing",
        detector_hint="fake-align-detector",
        threads=2,
        parallel=2,
    )

    outdir = tmp_path / "detector"
    results = run_detector_alignment_manifest(
        detector=_MetaAwareDetector(),
        manifest=alignment_manifest,
        outdir=outdir,
        threads=2,
        parallel=2,
    )

    assert len(results) == 2
    summary = json.loads((outdir / "detector_run_summary.json").read_text(encoding="utf-8"))
    assert summary["cells"][0]["input_file_type"] == "STAR tuple"
    assert summary["cells"][0]["input_sortedness"] == "unsorted"
    assert summary["cells"][0]["mapper_mode"] == "1"


def test_ciri3_registered_with_alignment_capabilities() -> None:
    engines = build_default_engines()
    assert "ciri3" in engines
    caps = get_detector_capabilities(engines["ciri3"])
    assert caps.accepts_alignment is True
    assert caps.accepts_fastq is False
    assert caps.supports_multisample_alignment is True
    assert caps.requires_unsorted_sam is True
    assert caps.supports_bwa is True
    assert caps.supports_star is True
    assert caps.recommended_execution_mode == "alignment-first"


def test_ciri3_runs_from_alignment_with_command_template(tmp_path: Path) -> None:
    bam = tmp_path / "cell1.bam"
    bam.write_text("bam", encoding="utf-8")
    script = tmp_path / "fake_ciri3.sh"
    script.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        "cat > \"$1\" <<'EOF'\n"
        "circRNA_ID\tchr\tstart\tend\tstrand\tbsj_reads\n"
        "circA\tchr1\t11\t22\t+\t7\n"
        "EOF\n",
        encoding="utf-8",
    )
    script.chmod(0o755)

    detector = Ciri3Detector(
        command_template=(
            f"bash {script} {{raw_output}}"
            " # {alignment} {alignment_format} {cell_id} {threads} {outdir}"
        )
    )
    result = detector.run_from_alignment(
        DetectorRunInputs(
            cell_id="cell1",
            bam=bam,
            outdir=tmp_path / "out",
            threads=2,
            input_mode="alignment",
            read_layout="paired-end",
        )
    )
    text = result.tsv_path.read_text(encoding="utf-8")
    assert "circ_id\tchr\tstart\tend\tstrand\tsupport" in text
    assert "circA\tchr1\t11\t22\t+\t7" in text


def test_ciri3_validate_runtime_accepts_direct_jar_contract(tmp_path: Path, monkeypatch) -> None:
    ciri3_home = tmp_path / "CIRI3"
    ciri3_home.mkdir()
    jar = ciri3_home / "CIRI3_Java_18.0.1.jar"
    jar.write_text("jar", encoding="utf-8")
    java = tmp_path / "java"
    java.write_text("java", encoding="utf-8")

    monkeypatch.setenv("CIRCYTO_CIRI3_HOME", str(ciri3_home))
    monkeypatch.setenv("CIRCYTO_CIRI3_JAVA", str(java))

    detector = Ciri3Detector()
    ok, errors, details = detector.validate_runtime()

    assert ok
    assert errors == []
    assert details["preferred_mode"] == "direct-jar"
    assert details["jar"] == str(jar)
    assert details["java"] == str(java)


def test_normalize_ciri3_output_handles_upstream_column_names(tmp_path: Path) -> None:
    raw = tmp_path / "raw.tsv"
    raw.write_text(
        "circRNA_ID\tchr\tcircRNA_start\tcircRNA_end\t#junction_reads\tstrand\n"
        "chr1:11|22\tchr1\t11\t22\t5\t+\n",
        encoding="utf-8",
    )
    out = tmp_path / "out.tsv"

    _normalize_ciri3_output(raw, out)

    text = out.read_text(encoding="utf-8")
    assert "circ_id\tchr\tstart\tend\tstrand\tsupport" in text
    assert "chr1:11|22\tchr1\t11\t22\t+\t5" in text


def test_ciri3_runs_from_alignment_with_direct_jar_contract(tmp_path: Path, monkeypatch) -> None:
    sam = tmp_path / "cell1.sam"
    sam.write_text("@HD\tVN:1.0\n", encoding="utf-8")
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n", encoding="utf-8")
    ciri3_home = tmp_path / "CIRI3"
    ciri3_home.mkdir()
    jar = ciri3_home / "CIRI3_Java_18.0.1.jar"
    jar.write_text("jar", encoding="utf-8")
    java = tmp_path / "java"
    java.write_text("java", encoding="utf-8")

    monkeypatch.setenv("CIRCYTO_CIRI3_HOME", str(ciri3_home))
    monkeypatch.setenv("CIRCYTO_CIRI3_JAVA", str(java))

    seen_cmds = []

    def fake_run(cmd, **kwargs):
        seen_cmds.append(cmd)
        out = tmp_path / "out" / "cell1.ciri3_run" / "ciri3_raw.tsv"
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(
            "circRNA_ID\tchr\tstart\tend\tstrand\tbsj_reads\n"
            "circA\tchr1\t11\t22\t+\t7\n",
            encoding="utf-8",
        )

        class _Result:
            returncode = 0

        return _Result()

    monkeypatch.setattr("circyto.detectors.ciri3.subprocess.run", fake_run)

    detector = Ciri3Detector()
    result = detector.run_from_alignment(
        DetectorRunInputs(
            cell_id="cell1",
            sam=sam,
            outdir=tmp_path / "out",
            ref_fa=ref,
            threads=2,
            input_mode="alignment",
            read_layout="single-end",
        )
    )

    text = result.tsv_path.read_text(encoding="utf-8")
    assert "circA\tchr1\t11\t22\t+\t7" in text
    assert result.meta["execution_mode"] == "direct-jar"
    assert result.meta["ciri3_jar"] == str(jar)
    assert result.meta["java_executable"] == str(java)
    assert result.meta["input_file_type"] == "SAM"
    assert result.meta["mapper_mode"] == "0"
    assert result.meta["stringency"] == "0"
    assert any(isinstance(cmd, list) and "-S" in cmd and "-Ma" in cmd for cmd in seen_cmds)


def test_ciri3_explicit_command_template_overrides_direct_jar(tmp_path: Path, monkeypatch) -> None:
    bam = tmp_path / "cell1.bam"
    bam.write_text("bam", encoding="utf-8")
    ciri3_home = tmp_path / "CIRI3"
    ciri3_home.mkdir()
    jar = ciri3_home / "CIRI3_Java_18.0.1.jar"
    jar.write_text("jar", encoding="utf-8")
    java = tmp_path / "java"
    java.write_text("java", encoding="utf-8")

    monkeypatch.setenv("CIRCYTO_CIRI3_HOME", str(ciri3_home))
    monkeypatch.setenv("CIRCYTO_CIRI3_JAVA", str(java))

    seen = []

    def fake_run(cmd, **kwargs):
        seen.append((cmd, kwargs))
        out = tmp_path / "out" / "cell1.ciri3_run" / "ciri3_raw.tsv"
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(
            "circRNA_ID\tchr\tstart\tend\tstrand\tbsj_reads\n"
            "circA\tchr1\t11\t22\t+\t7\n",
            encoding="utf-8",
        )

        class _Result:
            returncode = 0

        return _Result()

    monkeypatch.setattr("circyto.detectors.ciri3.subprocess.run", fake_run)

    detector = Ciri3Detector(
        command_template="bash -lc 'printf \"ignored\" >/dev/null' # {alignment} {raw_output} {cell_id}"
    )
    result = detector.run_from_alignment(
        DetectorRunInputs(
            cell_id="cell1",
            bam=bam,
            outdir=tmp_path / "out",
            threads=2,
            input_mode="alignment",
            read_layout="paired-end",
        )
    )

    assert result.meta["execution_mode"] == "template"
    assert seen
    cmd, kwargs = seen[0]
    assert isinstance(cmd, str)
    assert kwargs["shell"] is True


def test_cli_help_lists_alignment_first_commands() -> None:
    runner = CliRunner()
    result = runner.invoke(app, ["--help"])
    assert result.exit_code == 0
    assert "prepare-alignment-cache" in result.stdout
    assert "run-detector-from-alignments" in result.stdout
    assert "align-manifest" in result.stdout


def test_prepare_alignment_cache_bwa_mem_from_fastq_manifest(tmp_path: Path, monkeypatch) -> None:
    manifest = _write_fastq_manifest(tmp_path, ["c1", "c2"])
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    _install_fake_bwa_and_samtools(fake_bin)
    monkeypatch.setenv("PATH", f"{fake_bin}{os.pathsep}{os.environ['PATH']}")
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n", encoding="utf-8")

    alignment_manifest = prepare_alignment_cache(
        manifest=manifest,
        outdir=tmp_path / "align",
        aligner="bwa-mem",
        ref_fa=ref,
        threads=2,
        parallel=2,
        chunk_size=1,
        index_bam=True,
        output_format="bam",
    )
    rows = read_alignment_manifest_tsv(alignment_manifest, validate_files=True)
    assert len(rows) == 2
    for row in rows:
        bam = Path(row.bam)
        assert bam.exists()
        assert bam.suffix == ".bam"
        assert Path(str(bam) + ".bai").exists()


def test_prepare_alignment_cache_bwa_mem_for_ciri3_stages_unsorted_sam(tmp_path: Path, monkeypatch) -> None:
    manifest = _write_fastq_manifest(tmp_path, ["c1"])
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    _install_fake_bwa_and_samtools(fake_bin)
    monkeypatch.setenv("PATH", f"{fake_bin}{os.pathsep}{os.environ['PATH']}")
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n", encoding="utf-8")

    alignment_manifest = prepare_alignment_cache(
        manifest=manifest,
        outdir=tmp_path / "align",
        aligner="bwa-mem",
        detector_hint="ciri3",
        ref_fa=ref,
        threads=2,
        parallel=1,
    )
    rows = read_alignment_manifest_tsv(alignment_manifest, validate_files=True)
    assert len(rows) == 1
    assert rows[0].sam
    assert rows[0].bam == ""
    assert Path(rows[0].sam).suffix == ".sam"
    assert rows[0].artifact_bucket == "bwa_mem"
    assert rows[0].sortedness == "unsorted"
    assert rows[0].mapper_mode == "0"


def test_prepare_alignment_cache_bwa_mem_single_end_uses_only_r1(tmp_path: Path, monkeypatch) -> None:
    r1 = tmp_path / "c1_R1.fastq"
    r2 = tmp_path / "c1_R2.fastq"
    r1.write_text("@r1\nACGT\n+\n!!!!\n", encoding="utf-8")
    r2.write_text("@r2\nTGCA\n+\n!!!!\n", encoding="utf-8")
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        "\n".join(
            [
                "cell_id\tplatform\tread1\tread2\tbam\tlibrary_id\tn_input_reads\tgroup_id\tread_layout",
                f"c1\tsmartseq2\t{r1}\t{r2}\t\tlib1\t100\tplateA\tsingle-end",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    _install_fake_bwa_and_samtools(fake_bin)
    monkeypatch.setenv("PATH", f"{fake_bin}{os.pathsep}{os.environ['PATH']}")
    args_log = tmp_path / "bwa_args.log"
    monkeypatch.setenv("BWA_ARGS_LOG", str(args_log))
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n", encoding="utf-8")

    alignment_manifest = prepare_alignment_cache(
        manifest=manifest,
        outdir=tmp_path / "align",
        aligner="bwa-mem",
        detector_hint="ciri3",
        ref_fa=ref,
        threads=2,
        parallel=1,
    )
    rows = read_alignment_manifest_tsv(alignment_manifest, validate_files=True)
    assert len(rows) == 1
    assert rows[0].read_layout == "single-end"
    logged = args_log.read_text(encoding="utf-8").strip()
    assert str(r1) in logged
    assert str(r2) not in logged


def test_prepare_alignment_cache_bwa_mem_single_end_ignores_whitespace_read2(tmp_path: Path, monkeypatch) -> None:
    r1 = tmp_path / "c1_R1.fastq"
    r1.write_text("@r1\nACGT\n+\n!!!!\n", encoding="utf-8")
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        "\n".join(
            [
                "cell_id\tplatform\tread1\tread2\tbam\tlibrary_id\tn_input_reads\tgroup_id",
                f"c1\tsmartseq2\t{r1}\t \t\tlib1\t100\tplateA",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    _install_fake_bwa_and_samtools(fake_bin)
    monkeypatch.setenv("PATH", f"{fake_bin}{os.pathsep}{os.environ['PATH']}")
    args_log = tmp_path / "bwa_args.log"
    monkeypatch.setenv("BWA_ARGS_LOG", str(args_log))
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n", encoding="utf-8")

    alignment_manifest = prepare_alignment_cache(
        manifest=manifest,
        outdir=tmp_path / "align",
        aligner="bwa-mem",
        detector_hint="ciri3",
        ref_fa=ref,
        threads=2,
        parallel=1,
    )
    rows = read_alignment_manifest_tsv(alignment_manifest, validate_files=True)
    assert rows[0].read_layout == "single-end"
    logged = args_log.read_text(encoding="utf-8").strip()
    assert str(r1) in logged
    assert logged.split()[-1] == str(r1)


def test_prepare_alignment_cache_star_for_ciri3_records_chimeric_junction(tmp_path: Path, monkeypatch) -> None:
    manifest = _write_fastq_manifest(tmp_path, ["c1"])
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    _install_fake_star(fake_bin)
    _install_fake_bwa_and_samtools(fake_bin)
    monkeypatch.setenv("PATH", f"{fake_bin}{os.pathsep}{os.environ['PATH']}")
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n", encoding="utf-8")

    alignment_manifest = prepare_alignment_cache(
        manifest=manifest,
        outdir=tmp_path / "align",
        aligner="star",
        detector_hint="ciri3",
        ref_fa=ref,
        threads=2,
        parallel=1,
        extra_flags="--genomeDir /tmp/star-index",
    )
    rows = read_alignment_manifest_tsv(alignment_manifest, validate_files=True)
    assert len(rows) == 1
    assert rows[0].sam
    assert rows[0].mapper_mode == "1"
    assert rows[0].sortedness == "unsorted"
    assert Path(rows[0].chimeric_junction).exists()
    assert Path(rows[0].unmapped_mate1).exists()
    assert Path(rows[0].unmapped_mate2).exists()
    assert Path(rows[0].bwa_sam).exists()
    assert rows[0].artifact_bucket == "star"
    assert Path(rows[0].bwa_sam).read_text(encoding="utf-8").startswith("@HD\tVN:1.0")
    summary = json.loads((tmp_path / "align" / "alignment_prepare_summary.json").read_text(encoding="utf-8"))
    assert summary["status_counts"]["aligned"] == 1
    log_text = (tmp_path / "align" / "cache" / rows[0].cache_key / "c1.align.log").read_text(encoding="utf-8")
    assert "stage=prepare-start" in log_text
    assert "stage=star-start" in log_text
    assert "stage=star-end returncode=0" in log_text
    assert "stage=star-copy-end" in log_text
    assert "stage=bwa-rescue-start" in log_text
    assert "stage=bwa-rescue-end" in log_text
    assert "stage=stage-artifacts-end" in log_text


def test_prepare_alignment_cache_star_without_ciri3_does_not_require_rescue_outputs(tmp_path: Path, monkeypatch) -> None:
    manifest = _write_fastq_manifest(tmp_path, ["c1"])
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    _install_fake_star(fake_bin)
    monkeypatch.setenv("PATH", f"{fake_bin}{os.pathsep}{os.environ['PATH']}")

    alignment_manifest = prepare_alignment_cache(
        manifest=manifest,
        outdir=tmp_path / "align",
        aligner="star",
        detector_hint="circexplorer2",
        threads=2,
        parallel=1,
        extra_flags="--genomeDir /tmp/star-index",
    )
    rows = read_alignment_manifest_tsv(alignment_manifest, validate_files=True)
    assert len(rows) == 1
    assert rows[0].bam
    assert rows[0].sam == ""
    assert rows[0].chimeric_junction
    assert rows[0].unmapped_mate1 == ""
    assert rows[0].unmapped_mate2 == ""
    assert rows[0].bwa_sam == ""


def test_prepare_alignment_cache_resume_reuses_cached_chunks(tmp_path: Path) -> None:
    manifest = _write_source_manifest(tmp_path)
    outdir = tmp_path / "align"
    prepare_alignment_cache(
        manifest=manifest,
        outdir=outdir,
        aligner="reuse-existing",
        detector_hint="fake-align-detector",
        threads=2,
        parallel=2,
        chunk_size=1,
    )
    prepare_alignment_cache(
        manifest=manifest,
        outdir=outdir,
        aligner="reuse-existing",
        detector_hint="fake-align-detector",
        threads=2,
        parallel=2,
        chunk_size=1,
    )
    summary = json.loads((outdir / "alignment_prepare_summary.json").read_text(encoding="utf-8"))
    assert summary["status_counts"]["reused_cached"] == 2
    assert (outdir / "chunks" / "chunk_00001.json").exists()
    assert (outdir / "chunks" / "chunk_00002.json").exists()


def test_plan_alignment_cache_and_dry_run(tmp_path: Path) -> None:
    manifest = _write_source_manifest(tmp_path)
    plan = plan_alignment_cache(
        manifest=manifest,
        outdir=tmp_path / "align",
        aligner="reuse-existing",
        detector_hint="fake-align-detector",
        chunk_size=1,
    )
    assert plan["n_rows"] == 2
    assert plan["chunk_count"] == 2
    assert plan["errors"] == []
    assert plan["command_preview"]

    prepare_alignment_cache(
        manifest=manifest,
        outdir=tmp_path / "align",
        aligner="reuse-existing",
        detector_hint="fake-align-detector",
        dry_run=True,
    )
    assert (tmp_path / "align" / "alignment_prepare_plan.json").exists()


def test_prepare_alignment_cache_partial_failure_writes_summary(tmp_path: Path, monkeypatch) -> None:
    manifest = _write_fastq_manifest(tmp_path, ["c1", "c2"])
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    _install_fake_bwa_and_samtools(fake_bin, fail_on="c2_R1.fastq")
    monkeypatch.setenv("PATH", f"{fake_bin}{os.pathsep}{os.environ['PATH']}")
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n", encoding="utf-8")

    try:
        prepare_alignment_cache(
            manifest=manifest,
            outdir=tmp_path / "align",
            aligner="bwa-mem",
            ref_fa=ref,
            threads=2,
            parallel=2,
            chunk_size=1,
        )
        raise AssertionError("expected RuntimeError")
    except RuntimeError as exc:
        assert "completed with failures" in str(exc)

    summary = json.loads((tmp_path / "align" / "alignment_prepare_summary.json").read_text(encoding="utf-8"))
    assert summary["status_counts"]["aligned"] == 1
    assert summary["status_counts"]["failed"] == 1
    failed = next(cell for cell in summary["cells"] if cell["status"] == "failed")
    assert failed["read_layout"] == "paired-end"
    assert "bwa mem" in failed["command"]
    assert "c2_R1.fastq" in failed["command"]
    assert failed["log_path"].endswith("c2.align.log")
    assert failed["stderr_tail"]
    failure_detail = summary["failures"][0]
    assert failure_detail["read_layout"] == "paired-end"
    assert "bwa mem" in failure_detail["command"]
    assert failure_detail["read2"] is not None
    assert failure_detail["stderr_tail"]


def test_validate_ciri3_template_reports_missing_placeholders() -> None:
    ok, errors, details = validate_ciri3_template("ciri3 --bam {alignment}")
    assert not ok
    assert any("Missing required placeholders" in err for err in errors)
    assert "alignment" in details["fields"]


def test_validate_ciri3_template_cli(tmp_path: Path, monkeypatch) -> None:
    runner = CliRunner()
    template = "ciri3 --bam {alignment} --out {raw_output} --threads {threads} --cell {cell_id} --tmp {outdir} --fmt {alignment_format}"
    result = runner.invoke(app, ["validate-ciri3-template", "--template", template])
    assert result.exit_code == 0
    assert "OK" in result.stdout


def test_ciri3_star_manifest_row_requires_bwa_rescue(tmp_path: Path, monkeypatch) -> None:
    sam = tmp_path / "cell1.sam"
    sam.write_text("@HD\tVN:1.0\n", encoding="utf-8")
    chimeric = tmp_path / "cell1.Chimeric.out.junction"
    chimeric.write_text("chr1\t10\t+\tchr1\t20\t+\t0\t0\t10\tread1\t10\t10M\t20\t10M\n", encoding="utf-8")
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n", encoding="utf-8")
    ciri3_home = tmp_path / "CIRI3"
    ciri3_home.mkdir()
    (ciri3_home / "CIRI3_Java_18.0.1.jar").write_text("jar", encoding="utf-8")
    java = tmp_path / "java"
    java.write_text("java", encoding="utf-8")
    monkeypatch.setenv("CIRCYTO_CIRI3_HOME", str(ciri3_home))
    monkeypatch.setenv("CIRCYTO_CIRI3_JAVA", str(java))

    detector = Ciri3Detector()
    try:
        detector.run_from_alignment(
            DetectorRunInputs(
                cell_id="cell1",
                sam=sam,
                outdir=tmp_path / "out",
                ref_fa=ref,
                threads=2,
                input_mode="alignment",
                read_layout="paired-end",
                extra={
                    "alignment_manifest_row": {
                        "cell_id": "cell1",
                        "sam": str(sam),
                        "aligner": "star",
                        "mapper_mode": "1",
                        "sortedness": "unsorted",
                        "chimeric_junction": str(chimeric),
                    }
                },
            )
        )
        raise AssertionError("expected RuntimeError")
    except RuntimeError as exc:
        assert "requires bwa_sam" in str(exc)


def test_ciri3_runs_from_star_hybrid_manifest_row_with_direct_jar_contract(tmp_path: Path, monkeypatch) -> None:
    sam = tmp_path / "cell1.sam"
    sam.write_text("@HD\tVN:1.0\n", encoding="utf-8")
    chimeric = tmp_path / "cell1.Chimeric.out.junction"
    chimeric.write_text("chr1\t10\t+\tchr1\t20\t+\t0\t0\t10\tread1\t10\t10M\t20\t10M\n", encoding="utf-8")
    bwa_sam = tmp_path / "cell1.bwa_rescue.sam"
    bwa_sam.write_text("@SQ\tSN:chr1\tLN:4\n", encoding="utf-8")
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n", encoding="utf-8")
    ciri3_home = tmp_path / "CIRI3"
    ciri3_home.mkdir()
    jar = ciri3_home / "CIRI3_Java_18.0.1.jar"
    jar.write_text("jar", encoding="utf-8")
    java = tmp_path / "java"
    java.write_text("java", encoding="utf-8")

    monkeypatch.setenv("CIRCYTO_CIRI3_HOME", str(ciri3_home))
    monkeypatch.setenv("CIRCYTO_CIRI3_JAVA", str(java))

    seen_cmds = []

    def fake_run(cmd, **kwargs):
        seen_cmds.append(cmd)
        out = tmp_path / "out" / "cell1.ciri3_run" / "ciri3_raw.tsv"
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(
            "circRNA_ID\tchr\tstart\tend\tstrand\tbsj_reads\n"
            "circA\tchr1\t11\t22\t+\t7\n",
            encoding="utf-8",
        )

        class _Result:
            returncode = 0

        return _Result()

    monkeypatch.setattr("circyto.detectors.ciri3.subprocess.run", fake_run)

    detector = Ciri3Detector()
    result = detector.run_from_alignment(
        DetectorRunInputs(
            cell_id="cell1",
            sam=sam,
            outdir=tmp_path / "out",
            ref_fa=ref,
            threads=2,
            input_mode="alignment",
            read_layout="paired-end",
            extra={
                "alignment_manifest_row": {
                    "cell_id": "cell1",
                    "sam": str(sam),
                    "aligner": "star",
                    "mapper_mode": "1",
                    "sortedness": "unsorted",
                    "artifact_bucket": "star",
                    "chimeric_junction": str(chimeric),
                    "bwa_sam": str(bwa_sam),
                }
            },
        )
    )

    text = result.tsv_path.read_text(encoding="utf-8")
    assert "circA\tchr1\t11\t22\t+\t7" in text
    assert result.meta["execution_mode"] == "direct-jar"
    assert result.meta["input_file_type"] == "STAR tuple"
    assert result.meta["mapper_mode"] == "1"
    assert result.meta["bwa_sam"] == str(bwa_sam)
    assert f"{chimeric},{sam},{bwa_sam}" in result.meta["command"]
    assert " -Ma 1" in f" {result.meta['command']}"


def test_run_detector_from_alignments_cli_accepts_ciri3_template(tmp_path: Path) -> None:
    manifest = _write_fastq_manifest(tmp_path, ["c1", "c2"])
    alignment_dir = tmp_path / "align"
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    _install_fake_bwa_and_samtools(fake_bin)
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n", encoding="utf-8")
    original_path = os.environ["PATH"]
    os.environ["PATH"] = f"{fake_bin}{os.pathsep}{original_path}"
    alignment_manifest = prepare_alignment_cache(
        manifest=manifest,
        outdir=alignment_dir,
        aligner="bwa-mem",
        detector_hint="ciri3",
        ref_fa=ref,
    )
    runner = CliRunner()
    template = f"bash /mnt/d/circyto/tools/ciri3_smoke_template.sh {{alignment}} {{raw_output}} {{cell_id}} # {{alignment_format}} {{threads}} {{outdir}}"
    result = runner.invoke(
        app,
        [
            "run-detector-from-alignments",
            "--detector",
            "ciri3",
            "--manifest",
            str(alignment_manifest),
            "--outdir",
            str(tmp_path / 'ciri3'),
            "--ref-fa",
            str(ref),
            "--command-template",
            template,
        ],
    )
    os.environ["PATH"] = original_path
    assert result.exit_code == 0


def test_run_detector_from_alignment_manifest_dry_run(tmp_path: Path) -> None:
    manifest = _write_source_manifest(tmp_path)
    alignment_dir = tmp_path / "align"
    alignment_manifest = prepare_alignment_cache(
        manifest=manifest,
        outdir=alignment_dir,
        aligner="reuse-existing",
        detector_hint="fake-align-detector",
    )
    outdir = tmp_path / "detector"
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n", encoding="utf-8")
    run_detector_alignment_manifest(
        detector=_FakeAlignmentDetector(),
        manifest=alignment_manifest,
        outdir=outdir,
        ref_fa=ref,
        dry_run=True,
    )
    assert (outdir / "detector_alignment_plan.json").exists()
    payload = json.loads((outdir / "detector_alignment_plan.json").read_text(encoding="utf-8"))
    assert payload["n_manifest_rows"] == 2


def test_export_run_subset_failed_and_chunk(tmp_path: Path) -> None:
    manifest = _write_source_manifest(tmp_path)
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    (run_dir / "alignment_prepare_summary.json").write_text(
        json.dumps(
            {
                "cells": [
                    {"cell_id": "c1", "status": "failed"},
                    {"cell_id": "c2", "status": "reused_input"},
                ]
            }
        ),
        encoding="utf-8",
    )
    chunks = run_dir / "chunks"
    chunks.mkdir()
    (chunks / "chunk_00001.json").write_text(
        json.dumps(
            {
                "chunk_index": 1,
                "status": "partial_failure",
                "chunk_size": 2,
                "cells": [
                    {"cell_id": "c1", "status": "failed"},
                    {"cell_id": "c2", "status": "reused_input"},
                ],
            }
        ),
        encoding="utf-8",
    )
    failed_manifest = export_manifest_subset(
        manifest=manifest,
        run_dir=run_dir,
        out_path=tmp_path / "failed.tsv",
        subset="failed",
    )
    assert "c1" in failed_manifest.read_text(encoding="utf-8")
    chunk_manifest = export_manifest_subset(
        manifest=manifest,
        run_dir=run_dir,
        out_path=tmp_path / "chunk1.tsv",
        subset="failed",
        chunk_index=1,
    )
    text = chunk_manifest.read_text(encoding="utf-8")
    assert "c1" in text and "c2" in text


def test_export_run_subset_incomplete_includes_stale_cells(tmp_path: Path) -> None:
    manifest = _write_source_manifest(tmp_path)
    alignment_dir = tmp_path / "align"
    alignment_manifest = prepare_alignment_cache(
        manifest=manifest,
        outdir=alignment_dir,
        aligner="reuse-existing",
        detector_hint="fake-align-detector",
    )
    run_dir = tmp_path / "detector"
    run_dir.mkdir()
    (run_dir / "c1.tsv").write_text("circ_id\tchr\tstart\tend\tstrand\tsupport\nx\tchr1\t1\t2\t+\t1\n", encoding="utf-8")
    (run_dir / "detector_run_summary.json").write_text(
        json.dumps(
            {
                "cells": [
                    {"cell_id": "c1", "status": "success"},
                    {"cell_id": "c2", "status": "success"},
                ]
            }
        ),
        encoding="utf-8",
    )

    incomplete_manifest = export_manifest_subset(
        manifest=alignment_manifest,
        run_dir=run_dir,
        out_path=tmp_path / "incomplete.tsv",
        subset="incomplete",
    )
    stale_manifest = export_manifest_subset(
        manifest=alignment_manifest,
        run_dir=run_dir,
        out_path=tmp_path / "stale.tsv",
        subset="stale",
    )

    incomplete_ids = [line.split("\t", 1)[0] for line in incomplete_manifest.read_text(encoding="utf-8").splitlines()[1:] if line]
    stale_ids = [line.split("\t", 1)[0] for line in stale_manifest.read_text(encoding="utf-8").splitlines()[1:] if line]
    assert incomplete_ids == ["c1", "c2"]
    assert stale_ids == ["c1", "c2"]


def test_summarize_run_state_detects_stale_detector_outputs(tmp_path: Path) -> None:
    manifest = _write_source_manifest(tmp_path)
    alignment_dir = tmp_path / "align"
    alignment_manifest = prepare_alignment_cache(
        manifest=manifest,
        outdir=alignment_dir,
        aligner="reuse-existing",
        detector_hint="fake-align-detector",
    )
    outdir = tmp_path / "detector"
    outdir.mkdir()
    (outdir / "c1.tsv").write_text("circ_id\tchr\tstart\tend\tstrand\tsupport\nx\tchr1\t1\t2\t+\t1\n", encoding="utf-8")
    (outdir / "detector_run_summary.json").write_text(
        json.dumps(
            {
                "cells": [
                    {"cell_id": "c1", "status": "success"},
                    {"cell_id": "c2", "status": "failed"},
                ]
            }
        ),
        encoding="utf-8",
    )
    state = summarize_run_state(manifest=alignment_manifest, run_dir=outdir, mode="detector")
    assert "c1" in state["stale_cells"]
    assert "c2" in state["failed_cell_ids"]


def test_detector_subset_rerun_preserves_whole_run_summary(tmp_path: Path) -> None:
    manifest = _write_source_manifest(tmp_path)
    alignment_dir = tmp_path / "align"
    alignment_manifest = prepare_alignment_cache(
        manifest=manifest,
        outdir=alignment_dir,
        aligner="reuse-existing",
        detector_hint="fake-align-detector",
    )
    outdir = tmp_path / "detector"
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n", encoding="utf-8")
    gtf = tmp_path / "genes.gtf"
    gtf.write_text('chr1\ttest\texon\t1\t10\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n', encoding="utf-8")

    run_detector_alignment_manifest(
        detector=_FakeAlignmentDetector(),
        manifest=alignment_manifest,
        outdir=outdir,
        ref_fa=ref,
        gtf=gtf,
        threads=2,
        parallel=2,
        chunk_size=1,
    )

    subset_manifest = tmp_path / "subset.tsv"
    subset_manifest.write_text(
        "\n".join(
            [
                "cell_id\tbam\tsam\tgroup_id\tread_layout\taligner\treference\tcache_key\tsource_manifest\treused_alignment",
                next(line for line in alignment_manifest.read_text(encoding="utf-8").splitlines()[1:] if line.startswith("c2\t")),
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    run_detector_alignment_manifest(
        detector=_FakeAlignmentDetector(),
        manifest=subset_manifest,
        outdir=outdir,
        ref_fa=ref,
        gtf=gtf,
        threads=2,
        parallel=1,
        chunk_size=1,
    )

    summary = json.loads((outdir / "detector_run_summary.json").read_text(encoding="utf-8"))
    assert summary["planned_cells"] == 2
    assert sorted(cell["cell_id"] for cell in summary["cells"]) == ["c1", "c2"]
    assert summary["n_chunks"] == 3


def test_prepare_alignment_cache_fail_fast_stops_after_first_failed_chunk(tmp_path: Path, monkeypatch) -> None:
    manifest = _write_fastq_manifest(tmp_path, ["c1", "c2"])
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    _install_fake_bwa_and_samtools(fake_bin, fail_on="c1_R1.fastq")
    monkeypatch.setenv("PATH", f"{fake_bin}{os.pathsep}{os.environ['PATH']}")
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n", encoding="utf-8")

    try:
        prepare_alignment_cache(
            manifest=manifest,
            outdir=tmp_path / "align",
            aligner="bwa-mem",
            ref_fa=ref,
            chunk_size=1,
            fail_fast=True,
        )
        raise AssertionError("expected RuntimeError")
    except RuntimeError:
        pass

    assert (tmp_path / "align" / "chunks" / "chunk_00001.json").exists()
    assert not (tmp_path / "align" / "chunks" / "chunk_00002.json").exists()
