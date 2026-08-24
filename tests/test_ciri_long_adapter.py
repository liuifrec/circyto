from __future__ import annotations

import csv
import json
import os
from pathlib import Path

import pytest

from circyto.manifest.ciri_long import (
    CiriLongManifestRow,
    write_ciri_long_manifest_tsv,
)
from circyto.pipeline.ciri_long_adapter import (
    COORDINATE_CONVERSION_RULE,
    NORMALIZED_COORDINATE_SYSTEM,
    OFFICIAL_COORDINATE_SYSTEM,
    check_ciri_long_readiness,
    normalize_ciri_long_outputs,
    run_ciri_long_call_stage,
    run_ciri_long_collapse_stage,
    sha256_file,
)


FIXTURES = Path(__file__).parents[1] / "testdata" / "ciri_long"


def _write_reads(path: Path, read_name: str) -> None:
    path.write_text(
        f"@{read_name}\nACGTACGTACGT\n+\nIIIIIIIIIIII\n",
        encoding="utf-8",
    )


def _write_reference_and_index(tmp_path: Path) -> tuple[Path, Path]:
    reference = tmp_path / "reference.fa"
    reference.write_text(">chr1\n" + "A" * 1000 + "\n", encoding="utf-8")
    for suffix in (".amb", ".ann", ".bwt", ".pac", ".sa"):
        Path(f"{reference}{suffix}").write_text(
            f"synthetic {suffix}\n", encoding="utf-8"
        )
    gtf = tmp_path / "reference.gtf"
    gtf.write_text(
        'chr1\ttest\texon\t1\t1000\t.\t+\t.\tgene_id "GENE1";\n',
        encoding="utf-8",
    )
    return reference, gtf


def _write_manifest(tmp_path: Path, sample_ids: list[str]) -> Path:
    rows: list[CiriLongManifestRow] = []
    for sample_id in sample_ids:
        reads = tmp_path / f"{sample_id}.fastq"
        _write_reads(reads, f"{sample_id}_read")
        rows.append(
            CiriLongManifestRow(
                sample_id=sample_id,
                reads_path=str(reads),
                source_accession=f"SYNTHETIC_{sample_id}",
                dataset_id="synthetic_rcrt",
                reference_id="synthetic_reference",
                reference_build="synthetic_v1",
                extra={"protocol_evidence": "synthetic RCRT fixture"},
            )
        )
    manifest = tmp_path / "ciri_long_manifest.tsv"
    write_ciri_long_manifest_tsv(rows, manifest)
    return manifest


def _install_fake_tools(tmp_path: Path) -> tuple[Path, Path, Path]:
    bindir = tmp_path / "bin"
    bindir.mkdir()
    command_log = tmp_path / "commands.jsonl"
    ciri_long = bindir / "CIRI-long"
    ciri_long.write_text(
        """#!/usr/bin/env python3
import json
import os
from pathlib import Path
import shutil
import sys

if "--version" in sys.argv:
    print("CIRI-long 1.1.0")
    raise SystemExit(0)

args = sys.argv[1:]
with Path(os.environ["FAKE_CIRI_COMMAND_LOG"]).open("a", encoding="utf-8") as handle:
    handle.write(json.dumps(args) + "\\n")
stage = args[0]
outdir = Path(args[args.index("-o") + 1])
prefix = args[args.index("-p") + 1]
outdir.mkdir(parents=True, exist_ok=True)
(outdir / "tmp").mkdir(exist_ok=True)
if stage == "call":
    (outdir / f"{prefix}.cand_circ.fa").write_text(">candidate\\nACGT\\n", encoding="utf-8")
    (outdir / f"{prefix}.low_confidence.fa").write_text(">low\\nACGT\\n", encoding="utf-8")
    (outdir / f"{prefix}.json").write_text("{}\\n", encoding="utf-8")
    (outdir / f"{prefix}.log").write_text("call complete\\n", encoding="utf-8")
    (outdir / "tmp" / f"{prefix}.ccs.fa").write_text(">ccs\\nACGT\\n", encoding="utf-8")
    (outdir / "tmp" / f"{prefix}.raw.fa").write_text(">raw\\nACGT\\n", encoding="utf-8")
    (outdir / "tmp" / "ss.idx").write_text("index\\n", encoding="utf-8")
elif stage == "collapse":
    source = Path(os.environ["FAKE_CIRI_COLLAPSE_FIXTURES"])
    for suffix in ("info", "expression", "isoforms", "reads"):
        shutil.copyfile(source / f"cohort.{suffix}", outdir / f"{prefix}.{suffix}")
    (outdir / f"{prefix}.log").write_text("collapse complete\\n", encoding="utf-8")
    (outdir / "tmp" / "ss.idx").write_text("index\\n", encoding="utf-8")
    (outdir / "tmp" / f"{prefix}.corrected.pkl").write_bytes(b"pickle")
else:
    raise SystemExit(f"unsupported fake CIRI-long stage: {stage}")
""",
        encoding="utf-8",
    )
    ciri_long.chmod(0o755)
    bwa = bindir / "bwa"
    bwa.write_text(
        """#!/usr/bin/env python3
import sys
sys.stderr.write("Program: bwa (alignment via Burrows-Wheeler transformation)\\n")
sys.stderr.write("Version: 0.7.17-r1188\\n")
raise SystemExit(1)
""",
        encoding="utf-8",
    )
    bwa.chmod(0o755)
    return ciri_long, bwa, command_log


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_readiness_checks_versions_reference_and_indexes(tmp_path: Path) -> None:
    reference, gtf = _write_reference_and_index(tmp_path)
    ciri_long, bwa, _ = _install_fake_tools(tmp_path)
    readiness = check_ciri_long_readiness(
        ciri_long=str(ciri_long),
        bwa=str(bwa),
        reference_fasta=reference,
        gtf=gtf,
    )
    assert readiness["ok"] is True
    assert readiness["tools"]["ciri_long"]["version"] == "1.1.0"
    assert readiness["tools"]["bwa"]["version"] == "0.7.17-r1188"
    assert readiness["reference"]["sha256"] == sha256_file(reference)
    assert len(readiness["bwa_index_assets"]) == 5
    assert readiness["gtf"]["sha256"] == sha256_file(gtf)

    Path(f"{reference}.sa").unlink()
    missing = check_ciri_long_readiness(
        ciri_long=str(ciri_long),
        bwa=str(bwa),
        reference_fasta=reference,
        gtf=gtf,
    )
    assert missing["ok"] is False
    assert any("Missing BWA reference index asset" in item for item in missing["errors"])


def test_call_and_collapse_are_distinct_shell_free_stages(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    reference, gtf = _write_reference_and_index(tmp_path)
    manifest = _write_manifest(tmp_path, ["sampleA", "sampleB"])
    ciri_long, bwa, command_log = _install_fake_tools(tmp_path)
    monkeypatch.setenv("FAKE_CIRI_COMMAND_LOG", str(command_log))
    monkeypatch.setenv("FAKE_CIRI_COLLAPSE_FIXTURES", str(FIXTURES))
    outdir = tmp_path / "work"

    call_summary = run_ciri_long_call_stage(
        manifest_path=manifest,
        reference_fasta=reference,
        outdir=outdir,
        gtf=gtf,
        threads=2,
        ciri_long=str(ciri_long),
        bwa=str(bwa),
        execute=True,
    )
    assert call_summary["status"] == "completed"
    assert call_summary["sample_count"] == 2
    call_manifest = Path(call_summary["call_manifest"])
    call_rows = _read_tsv(call_manifest)
    assert len(call_rows) == 2
    assert len({row["candidate_fasta"] for row in call_rows}) == 2
    for row in call_rows:
        provenance = json.loads(
            Path(row["provenance"]).read_text(encoding="utf-8")
        )
        assert provenance["shell"] is False
        assert provenance["argv"][1] == "call"
        assert provenance["chemistry_gate"]["decision"] == (
            "accepted_rcrt_circrna_enriched_cdna"
        )
        assert provenance["input_reads"]["sha256"]
        assert len(provenance["bwa_index_assets"]) == 5
        assert provenance["reference"]["reference_id"] == "synthetic_reference"
        assert provenance["reference"]["reference_build"] == "synthetic_v1"
        assert provenance["gtf"]["sha256"] == sha256_file(gtf)

    collapse_plan = run_ciri_long_collapse_stage(
        call_manifest_path=call_manifest,
        reference_fasta=reference,
        outdir=outdir,
        prefix="cohort",
        gtf=gtf,
        threads=2,
        ciri_long=str(ciri_long),
        bwa=str(bwa),
        execute=False,
    )
    assert collapse_plan["status"] == "planned"
    collapse_summary = run_ciri_long_collapse_stage(
        call_manifest_path=call_manifest,
        reference_fasta=reference,
        outdir=outdir,
        prefix="cohort",
        gtf=gtf,
        threads=2,
        ciri_long=str(ciri_long),
        bwa=str(bwa),
        execute=True,
    )
    assert collapse_summary["status"] == "completed"
    collapse_provenance = json.loads(
        Path(collapse_summary["provenance"]).read_text(encoding="utf-8")
    )
    assert collapse_provenance["shell"] is False
    assert collapse_provenance["argv"][1] == "collapse"
    assert len(collapse_provenance["source_samples"]) == 2
    assert all(
        item["candidate_fasta_sha256"]
        for item in collapse_provenance["source_samples"]
    )

    logged = [
        json.loads(line)
        for line in command_log.read_text(encoding="utf-8").splitlines()
    ]
    assert [argv[0] for argv in logged] == ["call", "call", "collapse"]
    assert all("-i" in argv and "-r" in argv for argv in logged)


def test_dry_run_writes_argv_without_detector_execution(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    reference, gtf = _write_reference_and_index(tmp_path)
    manifest = _write_manifest(tmp_path, ["sampleA"])
    ciri_long, bwa, command_log = _install_fake_tools(tmp_path)
    monkeypatch.setenv("FAKE_CIRI_COMMAND_LOG", str(command_log))

    summary = run_ciri_long_call_stage(
        manifest_path=manifest,
        reference_fasta=reference,
        outdir=tmp_path / "plan",
        gtf=gtf,
        ciri_long=str(ciri_long),
        bwa=str(bwa),
        execute=False,
    )
    assert summary["status"] == "planned"
    assert not command_log.exists()
    provenance = json.loads(
        (
            tmp_path
            / "plan"
            / "ciri_long"
            / "call"
            / "sampleA"
            / "provenance.json"
        ).read_text(encoding="utf-8")
    )
    assert provenance["executed"] is False
    assert provenance["argv"][1] == "call"
    assert provenance["shell"] is False

    executed = run_ciri_long_call_stage(
        manifest_path=manifest,
        reference_fasta=reference,
        outdir=tmp_path / "plan",
        gtf=gtf,
        ciri_long=str(ciri_long),
        bwa=str(bwa),
        execute=True,
    )
    assert executed["status"] == "completed"
    assert command_log.exists()


def test_normalize_official_outputs_preserves_entities_and_coordinates(
    tmp_path: Path,
) -> None:
    collapse_dir = tmp_path / "collapse"
    collapse_dir.mkdir()
    for suffix in ("info", "expression", "isoforms", "reads"):
        (collapse_dir / f"cohort.{suffix}").write_bytes(
            (FIXTURES / f"cohort.{suffix}").read_bytes()
        )

    summary = normalize_ciri_long_outputs(
        collapse_dir=collapse_dir,
        outdir=tmp_path / "work",
        prefix="cohort",
    )
    assert summary["record_counts"] == {
        "bsj_records": 2,
        "isoform_records": 3,
        "expression_records": 4,
        "isoform_usage_records": 6,
        "read_assignment_records": 2,
    }
    normalized = tmp_path / "work" / "ciri_long" / "normalized"
    bsj = _read_tsv(normalized / "circRNA_bsj.tsv")
    assert bsj[0]["circRNA_id"] == "chr1:101-250"
    assert bsj[0]["start"] == "100"
    assert bsj[0]["end"] == "250"
    assert bsj[0]["original_start"] == "101"
    assert bsj[0]["original_coordinate_system"] == OFFICIAL_COORDINATE_SYSTEM
    assert bsj[0]["normalized_coordinate_system"] == NORMALIZED_COORDINATE_SYSTEM
    assert bsj[0]["coordinate_conversion_rule"] == COORDINATE_CONVERSION_RULE
    assert bsj[0]["host_gene_name"] == "GeneOne"
    assert bsj[0]["original_info_line"].startswith("chr1\tCIRI-long")

    isoforms = _read_tsv(normalized / "circRNA_isoforms.tsv")
    assert isoforms[0]["parent_circRNA_id"] == "chr1:101-250"
    assert isoforms[0]["transcript_length"] == "100"
    assert json.loads(isoforms[0]["exon_block_structure"]) == [
        [100, 150],
        [200, 250],
    ]
    assert isoforms[0]["major_isoform_status"] == "reported_major"
    assert isoforms[1]["major_isoform_status"] == "other_reported_isoform"

    expression = _read_tsv(normalized / "circRNA_expression.tsv")
    assert {row["sample_id"] for row in expression} == {"sampleA", "sampleB"}
    usage = _read_tsv(normalized / "circRNA_isoform_usage.tsv")
    assert all(row["parent_circRNA_id"] for row in usage)
    reads = _read_tsv(normalized / "circRNA_read_assignments.tsv")
    assert reads[0]["read_id"] == "readA1"
    assert reads[0]["cirexons"] == "101-150|**,201-250|**"
    assert json.loads(reads[0]["original_row_json"])["sample"] == "sampleA"

    import_summary = json.loads(
        (normalized / "import_summary.json").read_text(encoding="utf-8")
    )
    assert import_summary["coordinate_policy"]["conversion_rule"] == (
        COORDINATE_CONVERSION_RULE
    )
    assert set(import_summary["normalized_outputs"]) == {
        "circRNA_bsj",
        "circRNA_isoforms",
        "circRNA_expression",
        "circRNA_isoform_usage",
        "circRNA_read_assignments",
    }


def test_normalizer_rejects_malformed_coordinates(tmp_path: Path) -> None:
    collapse_dir = tmp_path / "collapse"
    collapse_dir.mkdir()
    for suffix in ("expression", "isoforms", "reads"):
        (collapse_dir / f"cohort.{suffix}").write_bytes(
            (FIXTURES / f"cohort.{suffix}").read_bytes()
        )
    (collapse_dir / "cohort.info").write_text(
        'chr1\tCIRI-long\tcircRNA\t0\t10\t1\t+\t.\tcirc_id "bad";\n',
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="coordinates must be >= 1"):
        normalize_ciri_long_outputs(
            collapse_dir=collapse_dir,
            outdir=tmp_path / "work",
            prefix="cohort",
        )
