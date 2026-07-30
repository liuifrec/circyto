from __future__ import annotations

import csv
import json
import os
from pathlib import Path

import pytest

from circyto.manifest.alignment import read_alignment_manifest_tsv
from circyto.manifest.long_read import LongReadManifestRow, write_long_read_manifest_tsv
from circyto.pipeline.nanopore_archive import sha256_file
from circyto.pipeline.nanopore_interop import (
    build_minimap2_argv,
    prepare_nanopore_alignments,
    summarize_sam_records,
)


SAM_TEXT = """@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:1000
q_linear\t0\tchr1\t101\t60\t10M90N10M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIII
q_back\t0\tchr1\t501\t60\t10M10S\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIII\tSA:Z:chr1,101,+,10S10M,55,0;
q_back\t2048\tchr1\t101\t55\t10S10M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIII
q_forward\t0\tchr1\t101\t60\t10M10S\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIII
q_forward\t2048\tchr1\t501\t55\t10S10M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIII
q_secondary\t256\tchr1\t501\t50\t10M10S\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIII
q_secondary\t2304\tchr1\t101\t50\t10S10M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIII
q_unmapped\t4\t*\t0\t0\t*\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIII
"""


def _write_fastq(path: Path, names: list[str]) -> None:
    lines: list[str] = []
    for name in names:
        lines.extend([f"@{name}", "A" * 20, "+", "I" * 20])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _row(cell_id: str, fastq: Path, *, molecule_type: str = "cdna") -> LongReadManifestRow:
    return LongReadManifestRow(
        cell_id=cell_id,
        long_read_fastq=str(fastq),
        protocol="synthetic_nanopore",
        sequencing_platform="OXFORD_NANOPORE",
        archive_library_selection="RACE",
        library_preparation_summary="Synthetic cDNA interoperability fixture.",
        molecule_type=molecule_type,
        barcode_status="not_applicable_physical_single_cell",
        source_accession=f"SYNTHETIC_{cell_id}",
        dataset_id="synthetic_dataset",
    )


def _install_fake_tools(bin_dir: Path, sam_path: Path, command_log: Path) -> tuple[Path, Path]:
    bin_dir.mkdir(parents=True, exist_ok=True)
    minimap2 = bin_dir / "minimap2"
    minimap2.write_text(
        """#!/usr/bin/env python3
import json
import os
from pathlib import Path
import sys
if "--version" in sys.argv:
    print("2.28-r1209")
    raise SystemExit(0)
with Path(os.environ["FAKE_COMMAND_LOG"]).open("a", encoding="utf-8") as handle:
    handle.write(json.dumps({"tool": "minimap2", "argv": sys.argv[1:]}) + "\\n")
sys.stdout.write(Path(os.environ["FAKE_SAM_PATH"]).read_text(encoding="utf-8"))
""",
        encoding="utf-8",
    )
    minimap2.chmod(0o755)
    samtools = bin_dir / "samtools"
    samtools.write_text(
        """#!/usr/bin/env python3
import json
import os
from pathlib import Path
import sys
if "--version" in sys.argv:
    print("samtools 1.21")
    raise SystemExit(0)
with Path(os.environ["FAKE_COMMAND_LOG"]).open("a", encoding="utf-8") as handle:
    handle.write(json.dumps({"tool": "samtools", "argv": sys.argv[1:]}) + "\\n")
command = sys.argv[1]
if command == "sort":
    output = Path(sys.argv[sys.argv.index("-o") + 1])
    output.write_bytes(sys.stdin.buffer.read())
elif command == "index":
    Path(sys.argv[2] + ".bai").write_bytes(b"fake-index")
elif command == "view":
    sys.stdout.buffer.write(Path(sys.argv[-1]).read_bytes())
else:
    raise SystemExit(f"unsupported fake samtools command: {command}")
""",
        encoding="utf-8",
    )
    samtools.chmod(0o755)
    return minimap2, samtools


def _prepare_inputs(tmp_path: Path, cell_ids: list[str]) -> tuple[Path, Path, str]:
    rows: list[LongReadManifestRow] = []
    names = ["q_linear", "q_back", "q_forward", "q_secondary", "q_unmapped"]
    for cell_id in cell_ids:
        fastq = tmp_path / f"{cell_id}.fastq"
        _write_fastq(fastq, names)
        rows.append(_row(cell_id, fastq))
    manifest = tmp_path / "long_read_manifest.tsv"
    write_long_read_manifest_tsv(rows, manifest)
    reference = tmp_path / "reference.fa"
    reference.write_text(">chr1\n" + "A" * 1000 + "\n", encoding="utf-8")
    return manifest, reference, sha256_file(reference)


def test_exploratory_patterns_require_order_inversion_not_n_cigar() -> None:
    qc, evidence = summarize_sam_records(
        SAM_TEXT.splitlines(keepends=True),
        input_query_count=5,
        cell_id="cell_a",
    )

    assert qc["input_query_count"] == 5
    assert qc["mapped_primary_query_count"] == 3
    assert qc["unmapped_query_count"] == 2
    assert qc["secondary_alignment_record_count"] == 2
    assert qc["supplementary_alignment_record_count"] == 3
    assert qc["queries_with_supplementary_alignments"] == 3
    assert qc["queries_with_sa_tag"] == 1
    assert qc["spliced_primary_query_count"] == 1
    assert [row["read_name"] for row in evidence] == ["q_back"]
    assert evidence[0]["reason_code"] == (
        "query_reference_order_inversion_same_contig_strand"
    )
    assert evidence[0]["circRNA_call"] == "false"
    assert "q_linear" not in {row["read_name"] for row in evidence}
    assert "q_secondary" not in {row["read_name"] for row in evidence}


def test_reverse_strand_is_normalized_for_query_reference_order() -> None:
    reverse_sam = """@SQ\tSN:chr1\tLN:1000
q_reverse\t16\tchr1\t101\t60\t10S10M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIII
q_reverse\t2064\tchr1\t501\t55\t10M10S\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIII
"""
    _, evidence = summarize_sam_records(
        reverse_sam.splitlines(keepends=True),
        input_query_count=1,
        cell_id="cell_reverse",
    )
    assert len(evidence) == 1
    assert evidence[0]["strand"] == "-"


def test_cdna_and_direct_rna_minimap2_options_are_not_interchangeable(
    tmp_path: Path,
) -> None:
    fastq = tmp_path / "reads.fastq"
    _write_fastq(fastq, ["q1"])
    reference = tmp_path / "reference.fa"
    reference.write_text(">chr1\nA\n", encoding="utf-8")
    cdna = _row("cdna_cell", fastq)
    direct = _row("direct_cell", fastq, molecule_type="direct_rna")

    cdna_argv = build_minimap2_argv(
        cdna,
        minimap2="minimap2",
        reference_fasta=reference,
        threads=1,
    )
    direct_argv = build_minimap2_argv(
        direct,
        minimap2="minimap2",
        reference_fasta=reference,
        threads=1,
    )

    assert cdna_argv[:3] == ["minimap2", "-ax", "splice"]
    assert "-uf" not in cdna_argv
    assert "-k14" not in cdna_argv
    assert "-uf" in direct_argv
    assert "-k14" in direct_argv
    with pytest.raises(ValueError, match="controlled by molecule_type"):
        build_minimap2_argv(
            cdna,
            minimap2="minimap2",
            reference_fasta=reference,
            threads=1,
            extra_args=["-uf", "-k14"],
        )


@pytest.mark.parametrize("cell_ids", [["cell_a"], ["cell_a", "cell_b"]])
def test_alignment_streaming_and_cell_scoped_outputs(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    cell_ids: list[str],
) -> None:
    manifest, reference, reference_sha256 = _prepare_inputs(tmp_path, cell_ids)
    sam_fixture = tmp_path / "fixture.sam"
    sam_fixture.write_text(SAM_TEXT, encoding="utf-8")
    command_log = tmp_path / "commands.jsonl"
    minimap2, samtools = _install_fake_tools(tmp_path / "bin", sam_fixture, command_log)
    monkeypatch.setenv("FAKE_SAM_PATH", str(sam_fixture))
    monkeypatch.setenv("FAKE_COMMAND_LOG", str(command_log))
    outdir = tmp_path / "work"

    alignment_manifest = prepare_nanopore_alignments(
        manifest_path=manifest,
        reference_fasta=reference,
        reference_id="synthetic_reference",
        reference_build="synthetic_build_v1",
        reference_sha256=reference_sha256,
        outdir=outdir,
        threads=2,
        minimap2=str(minimap2),
        samtools=str(samtools),
    )

    root = outdir / "nanopore_alignment"
    assert not (root / "alignment.bam").exists()
    manifest_rows = read_alignment_manifest_tsv(
        alignment_manifest,
        validate_files=True,
    )
    assert len(manifest_rows) == len(cell_ids)
    assert len({row.bam for row in manifest_rows}) == len(cell_ids)
    for cell_id in cell_ids:
        cell_dir = root / cell_id
        assert (cell_dir / "alignment.bam").is_file()
        assert (cell_dir / "alignment.bam.bai").is_file()
        assert (cell_dir / "alignment_qc.json").is_file()
        assert (cell_dir / "exploratory_bsj_evidence.tsv").is_file()
        assert (cell_dir / "provenance.json").is_file()
        assert not (cell_dir / "alignment.sam").exists()
        qc = json.loads((cell_dir / "alignment_qc.json").read_text(encoding="utf-8"))
        assert qc["input_query_count"] == 5
        assert qc["mapped_primary_query_count"] == 3
        provenance = json.loads(
            (cell_dir / "provenance.json").read_text(encoding="utf-8")
        )
        assert provenance["commands"]["shell"] is False
        assert provenance["commands"]["streamed_minimap2_to_samtools_sort"] is True
        assert provenance["reference"]["reference_id"] == "synthetic_reference"
        assert provenance["reference"]["reference_build"] == "synthetic_build_v1"
        assert provenance["reference"]["fasta_sha256"] == reference_sha256
        assert provenance["detector_invoked"] is False
        assert provenance["circRNA_validation_status"] is False
        with (cell_dir / "exploratory_bsj_evidence.tsv").open(
            "r", encoding="utf-8", newline=""
        ) as handle:
            evidence = list(csv.DictReader(handle, delimiter="\t"))
        assert [row["read_name"] for row in evidence] == ["q_back"]
        assert all(row["circRNA_call"] == "false" for row in evidence)

    logged = [
        json.loads(line)
        for line in command_log.read_text(encoding="utf-8").splitlines()
    ]
    minimap_commands = [item for item in logged if item["tool"] == "minimap2"]
    sort_commands = [
        item
        for item in logged
        if item["tool"] == "samtools" and item["argv"][0] == "sort"
    ]
    assert len(minimap_commands) == len(cell_ids)
    assert len(sort_commands) == len(cell_ids)
    assert all(item["argv"][:2] == ["-ax", "splice"] for item in minimap_commands)
    assert all(item["argv"][-1].endswith(".fastq") for item in minimap_commands)
    assert all(item["argv"][-1] == "-" for item in sort_commands)
    if len(cell_ids) == 2:
        with pytest.raises(FileExistsError, match="output already exists"):
            prepare_nanopore_alignments(
                manifest_path=manifest,
                reference_fasta=reference,
                reference_id="synthetic_reference",
                reference_build="synthetic_build_v1",
                reference_sha256=reference_sha256,
                outdir=outdir,
                threads=2,
                minimap2=str(minimap2),
                samtools=str(samtools),
            )


def test_keep_sam_is_opt_in(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    manifest, reference, reference_sha256 = _prepare_inputs(tmp_path, ["cell_a"])
    sam_fixture = tmp_path / "fixture.sam"
    sam_fixture.write_text(SAM_TEXT, encoding="utf-8")
    command_log = tmp_path / "commands.jsonl"
    minimap2, samtools = _install_fake_tools(tmp_path / "bin", sam_fixture, command_log)
    monkeypatch.setenv("FAKE_SAM_PATH", str(sam_fixture))
    monkeypatch.setenv("FAKE_COMMAND_LOG", str(command_log))

    prepare_nanopore_alignments(
        manifest_path=manifest,
        reference_fasta=reference,
        reference_id="synthetic_reference",
        reference_build="synthetic_build_v1",
        reference_sha256=reference_sha256,
        outdir=tmp_path / "work",
        threads=1,
        minimap2=str(minimap2),
        samtools=str(samtools),
        keep_sam=True,
    )

    assert (
        tmp_path / "work" / "nanopore_alignment" / "cell_a" / "alignment.sam"
    ).read_text(encoding="utf-8") == SAM_TEXT


def test_reference_sha256_is_explicit_and_verified(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    manifest, reference, _ = _prepare_inputs(tmp_path, ["cell_a"])
    sam_fixture = tmp_path / "fixture.sam"
    sam_fixture.write_text(SAM_TEXT, encoding="utf-8")
    command_log = tmp_path / "commands.jsonl"
    minimap2, samtools = _install_fake_tools(tmp_path / "bin", sam_fixture, command_log)
    monkeypatch.setenv("FAKE_SAM_PATH", str(sam_fixture))
    monkeypatch.setenv("FAKE_COMMAND_LOG", str(command_log))

    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        prepare_nanopore_alignments(
            manifest_path=manifest,
            reference_fasta=reference,
            reference_id="explicit_id",
            reference_build="explicit_build",
            reference_sha256="0" * 64,
            outdir=tmp_path / "work",
            minimap2=str(minimap2),
            samtools=str(samtools),
        )
