from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

from typer.testing import CliRunner

from circyto.cli.circyto import app
from circyto.detectors import build_default_engines
from circyto.manifest.ciri_long import (
    CiriLongManifestRow,
    write_ciri_long_manifest_tsv,
)


runner = CliRunner()


def _fake_tools(tmp_path: Path) -> tuple[Path, Path]:
    ciri_long = tmp_path / "CIRI-long"
    ciri_long.write_text(
        """#!/usr/bin/env python3
import sys
if "--version" in sys.argv:
    print("CIRI-long 1.1.0")
    raise SystemExit(0)
raise SystemExit("execution was not expected")
""",
        encoding="utf-8",
    )
    ciri_long.chmod(0o755)
    bwa = tmp_path / "bwa"
    bwa.write_text(
        """#!/usr/bin/env python3
import sys
sys.stderr.write("Version: 0.7.17-r1188\\n")
raise SystemExit(1)
""",
        encoding="utf-8",
    )
    bwa.chmod(0o755)
    return ciri_long, bwa


def _inputs(tmp_path: Path) -> tuple[Path, Path, Path]:
    reads = tmp_path / "reads.fastq"
    reads.write_text("@r1\nACGT\n+\nIIII\n", encoding="utf-8")
    reference = tmp_path / "reference.fa"
    reference.write_text(">chr1\nACGT\n", encoding="utf-8")
    for suffix in (".amb", ".ann", ".bwt", ".pac", ".sa"):
        Path(f"{reference}{suffix}").write_text("index\n", encoding="utf-8")
    gtf = tmp_path / "reference.gtf"
    gtf.write_text(
        'chr1\ttest\texon\t1\t4\t.\t+\t.\tgene_id "G1";\n',
        encoding="utf-8",
    )
    manifest = tmp_path / "manifest.tsv"
    write_ciri_long_manifest_tsv(
        [
            CiriLongManifestRow(
                sample_id="sampleA",
                reads_path=str(reads),
                source_accession="SYNTHETIC",
                dataset_id="synthetic",
                reference_id="reference",
                reference_build="v1",
            )
        ],
        manifest,
    )
    return manifest, reference, gtf


def test_ciri_long_cli_registers_expected_commands() -> None:
    result = runner.invoke(app, ["ciri-long", "--help"])
    assert result.exit_code == 0, result.output
    for command in (
        "doctor",
        "validate-manifest",
        "plan",
        "call",
        "collapse",
        "import",
    ):
        assert command in result.output
    assert "ciri-long" not in build_default_engines()


def test_ciri_long_cli_validate_and_plan_are_dry_run(
    tmp_path: Path,
) -> None:
    manifest, reference, gtf = _inputs(tmp_path)
    ciri_long, bwa = _fake_tools(tmp_path)
    validate = runner.invoke(
        app,
        [
            "ciri-long",
            "validate-manifest",
            "--manifest",
            str(manifest),
            "--strict",
        ],
    )
    assert validate.exit_code == 0, validate.output
    assert json.loads(validate.output)["chemistry_gate"] == (
        "accepted_rcrt_circrna_enriched"
    )

    outdir = tmp_path / "work"
    plan = runner.invoke(
        app,
        [
            "ciri-long",
            "plan",
            "--manifest",
            str(manifest),
            "--reference",
            str(reference),
            "--gtf",
            str(gtf),
            "--outdir",
            str(outdir),
            "--ciri-long",
            str(ciri_long),
            "--bwa",
            str(bwa),
        ],
    )
    assert plan.exit_code == 0, plan.output
    payload = json.loads(plan.output)
    assert payload["status"] == "planned"
    provenance = json.loads(
        (
            outdir
            / "ciri_long"
            / "call"
            / "sampleA"
            / "provenance.json"
        ).read_text(encoding="utf-8")
    )
    assert provenance["executed"] is False
    assert provenance["argv"][1] == "call"


def test_official_demo_script_dry_run_is_network_free(tmp_path: Path) -> None:
    script = Path(__file__).parents[1] / "scripts" / "run_ciri_long_official_demo.py"
    result = subprocess.run(
        [
            sys.executable,
            str(script),
            "--outdir",
            str(tmp_path / "demo"),
            "--dry-run",
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    payload = json.loads(result.stdout)
    assert payload["network_accessed"] is False
    assert payload["asset_name"] == "CIRI-long_test_data.tar.gz"
    assert payload["scientific_boundary"].startswith("Official bulk")
    assert not (tmp_path / "demo").exists()
