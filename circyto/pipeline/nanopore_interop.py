from __future__ import annotations

import csv
import gzip
import json
import os
import re
import shlex
import shutil
import subprocess
import time
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Iterator, Sequence, TextIO

from circyto.manifest.alignment import AlignmentManifestRow, write_alignment_manifest_tsv
from circyto.manifest.long_read import LongReadManifestRow, read_long_read_manifest_tsv
from circyto.pipeline.nanopore_archive import sha256_file
from circyto.pipeline.workflow_reporting import (
    apply_standard_provenance,
    generate_workflow_uuid,
    utc_now_iso,
    write_json,
)


SCIENTIFIC_BOUNDARY_WARNING = (
    "EXPERIMENTAL LONG-READ INTEROPERABILITY ONLY. This workflow does not validate "
    "circRNA detection. SRR4048177 was generated using a poly(A)-oriented cDNA "
    "protocol, so conventional non-polyadenylated circRNAs are expected to be "
    "depleted. Non-detection must not be interpreted as biological absence."
)
EXPLORATORY_EVIDENCE_WARNING = (
    "exploratory_bsj_evidence contains alignment patterns compatible with a possible "
    "back-splice. Entries are not circRNA calls. PCR and template-switch artifacts, "
    "chimeric cDNA, basecalling errors, and long-read alignment errors can produce "
    "similar patterns."
)
NO_DETECTOR_WARNING = "No CIRI3 or other short-read circRNA detector was run."

EVIDENCE_COLUMNS = [
    "candidate_pattern_id",
    "cell_id",
    "read_name",
    "contig",
    "strand",
    "reason_code",
    "circRNA_call",
    "segment_a_flag",
    "segment_a_mapq",
    "segment_a_cigar",
    "segment_a_query_start",
    "segment_a_query_end",
    "segment_a_reference_start",
    "segment_a_reference_end",
    "segment_b_flag",
    "segment_b_mapq",
    "segment_b_cigar",
    "segment_b_query_start",
    "segment_b_query_end",
    "segment_b_reference_start",
    "segment_b_reference_end",
    "supporting_segments_json",
    "warning",
]

_CIGAR_TOKEN = re.compile(r"(\d+)([MIDNSHP=X])")
_SHA256 = re.compile(r"^[0-9a-f]{64}$")


@dataclass(frozen=True)
class SamSegment:
    read_name: str
    flag: int
    contig: str
    reference_start: int
    reference_end: int
    mapq: int
    cigar: str
    strand: str
    query_start: int
    query_end: int
    has_sa_tag: bool

    @property
    def is_unmapped(self) -> bool:
        return bool(self.flag & 0x4)

    @property
    def is_secondary(self) -> bool:
        return bool(self.flag & 0x100)

    @property
    def is_supplementary(self) -> bool:
        return bool(self.flag & 0x800)

    @property
    def is_primary(self) -> bool:
        return not self.is_secondary and not self.is_supplementary

    def to_dict(self) -> dict[str, Any]:
        return {
            "read_name": self.read_name,
            "flag": self.flag,
            "contig": self.contig,
            "reference_start": self.reference_start,
            "reference_end": self.reference_end,
            "mapq": self.mapq,
            "cigar": self.cigar,
            "strand": self.strand,
            "query_start": self.query_start,
            "query_end": self.query_end,
            "alignment_role": (
                "secondary"
                if self.is_secondary
                else "supplementary"
                if self.is_supplementary
                else "primary"
            ),
            "has_sa_tag": self.has_sa_tag,
        }


def _parse_cigar(cigar: str, *, reverse: bool) -> tuple[int, int, int]:
    if cigar == "*":
        return 0, 0, 0
    tokens = [(int(length), operation) for length, operation in _CIGAR_TOKEN.findall(cigar)]
    if not tokens or "".join(f"{length}{operation}" for length, operation in tokens) != cigar:
        raise ValueError(f"Unsupported or malformed CIGAR: {cigar!r}")
    left_clip = 0
    for length, operation in tokens:
        if operation not in {"S", "H"}:
            break
        left_clip += length
    right_clip = 0
    for length, operation in reversed(tokens):
        if operation not in {"S", "H"}:
            break
        right_clip += length
    query_aligned = sum(
        length for length, operation in tokens if operation in {"M", "I", "=", "X"}
    )
    query_length = sum(
        length for length, operation in tokens if operation in {"M", "I", "S", "=", "X"}
    ) + sum(length for length, operation in tokens if operation == "H")
    reference_length = sum(
        length for length, operation in tokens if operation in {"M", "D", "N", "=", "X"}
    )
    if reverse:
        query_start = right_clip
        query_end = right_clip + query_aligned
    else:
        query_start = left_clip
        query_end = left_clip + query_aligned
    if query_end > query_length:
        raise ValueError(f"CIGAR query interval exceeds query length: {cigar!r}")
    return query_start, query_end, reference_length


def parse_sam_segment(line: str) -> SamSegment:
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 11:
        raise ValueError(f"SAM record has fewer than 11 fields: {line[:120]!r}")
    flag = int(fields[1])
    reverse = bool(flag & 0x10)
    query_start, query_end, reference_length = _parse_cigar(fields[5], reverse=reverse)
    reference_start = max(0, int(fields[3]) - 1) if fields[3] != "0" else 0
    return SamSegment(
        read_name=fields[0],
        flag=flag,
        contig=fields[2],
        reference_start=reference_start,
        reference_end=reference_start + reference_length,
        mapq=int(fields[4]),
        cigar=fields[5],
        strand="-" if reverse else "+",
        query_start=query_start,
        query_end=query_end,
        has_sa_tag=any(field.startswith("SA:Z:") for field in fields[11:]),
    )


def _reference_order_coordinate(segment: SamSegment) -> int:
    if segment.strand == "+":
        return segment.reference_start
    return -segment.reference_end


def exploratory_bsj_patterns(
    segments_by_query: dict[str, list[SamSegment]],
    *,
    cell_id: str,
) -> list[dict[str, Any]]:
    evidence: list[dict[str, Any]] = []
    for read_name in sorted(segments_by_query):
        all_segments = segments_by_query[read_name]
        usable = [
            segment
            for segment in all_segments
            if not segment.is_unmapped and not segment.is_secondary
        ]
        primaries = [segment for segment in usable if segment.is_primary]
        supplementary = [segment for segment in usable if segment.is_supplementary]
        if not primaries or not supplementary:
            continue
        support_json = json.dumps(
            [segment.to_dict() for segment in sorted(usable, key=lambda item: item.query_start)],
            separators=(",", ":"),
            sort_keys=True,
        )
        pair_index = 0
        for primary in primaries:
            for supplement in supplementary:
                if primary.contig == "*" or primary.contig != supplement.contig:
                    continue
                if primary.strand != supplement.strand:
                    continue
                ordered = sorted(
                    (primary, supplement),
                    key=lambda item: (item.query_start, item.query_end, item.reference_start),
                )
                first, second = ordered
                if first.query_end > second.query_start:
                    continue
                first_reference_order = _reference_order_coordinate(first)
                second_reference_order = _reference_order_coordinate(second)
                if first_reference_order <= second_reference_order:
                    continue
                pair_index += 1
                evidence.append(
                    {
                        "candidate_pattern_id": f"{cell_id}:{read_name}:{pair_index}",
                        "cell_id": cell_id,
                        "read_name": read_name,
                        "contig": first.contig,
                        "strand": first.strand,
                        "reason_code": "query_reference_order_inversion_same_contig_strand",
                        "circRNA_call": "false",
                        "segment_a_flag": first.flag,
                        "segment_a_mapq": first.mapq,
                        "segment_a_cigar": first.cigar,
                        "segment_a_query_start": first.query_start,
                        "segment_a_query_end": first.query_end,
                        "segment_a_reference_start": first.reference_start,
                        "segment_a_reference_end": first.reference_end,
                        "segment_b_flag": second.flag,
                        "segment_b_mapq": second.mapq,
                        "segment_b_cigar": second.cigar,
                        "segment_b_query_start": second.query_start,
                        "segment_b_query_end": second.query_end,
                        "segment_b_reference_start": second.reference_start,
                        "segment_b_reference_end": second.reference_end,
                        "supporting_segments_json": support_json,
                        "warning": EXPLORATORY_EVIDENCE_WARNING,
                    }
                )
    return evidence


def summarize_sam_records(
    lines: Iterable[str],
    *,
    input_query_count: int,
    cell_id: str,
    retained_sam: TextIO | None = None,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    mapped_primary_queries: set[str] = set()
    supplementary_queries: set[str] = set()
    sa_queries: set[str] = set()
    spliced_primary_queries: set[str] = set()
    secondary_records = 0
    supplementary_records = 0
    alignment_records = 0
    segments_by_query: dict[str, list[SamSegment]] = defaultdict(list)

    for line in lines:
        if retained_sam is not None:
            retained_sam.write(line)
        if not line or line.startswith("@"):
            continue
        segment = parse_sam_segment(line)
        alignment_records += 1
        if segment.has_sa_tag:
            sa_queries.add(segment.read_name)
        if segment.is_secondary:
            secondary_records += 1
        if segment.is_supplementary:
            supplementary_records += 1
            supplementary_queries.add(segment.read_name)
        if segment.is_primary and not segment.is_unmapped:
            mapped_primary_queries.add(segment.read_name)
            if "N" in segment.cigar:
                spliced_primary_queries.add(segment.read_name)
        if not segment.is_unmapped:
            segments_by_query[segment.read_name].append(segment)

    mapped_count = len(mapped_primary_queries)
    if mapped_count > input_query_count:
        raise ValueError(
            f"Mapped primary query count exceeds input FASTQ query count for {cell_id}: "
            f"{mapped_count}>{input_query_count}"
        )
    unmapped_count = input_query_count - mapped_count
    evidence = exploratory_bsj_patterns(segments_by_query, cell_id=cell_id)
    qc = {
        "cell_id": cell_id,
        "input_query_count": input_query_count,
        "alignment_record_count": alignment_records,
        "mapped_primary_query_count": mapped_count,
        "unmapped_query_count": unmapped_count,
        "secondary_alignment_record_count": secondary_records,
        "supplementary_alignment_record_count": supplementary_records,
        "queries_with_supplementary_alignments": len(supplementary_queries),
        "queries_with_sa_tag": len(sa_queries),
        "spliced_primary_query_count": len(spliced_primary_queries),
        "mapped_primary_query_fraction": (
            mapped_count / input_query_count if input_query_count else 0.0
        ),
        "unmapped_query_fraction": (
            unmapped_count / input_query_count if input_query_count else 0.0
        ),
        "spliced_primary_query_fraction": (
            len(spliced_primary_queries) / mapped_count if mapped_count else 0.0
        ),
        "exploratory_bsj_pattern_count": len(evidence),
        "circRNA_validation_status": False,
        "warnings": [
            SCIENTIFIC_BOUNDARY_WARNING,
            EXPLORATORY_EVIDENCE_WARNING,
            NO_DETECTOR_WARNING,
        ],
    }
    return qc, evidence


def write_exploratory_evidence(path: Path, evidence: Iterable[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=EVIDENCE_COLUMNS, delimiter="\t")
        writer.writeheader()
        for row in evidence:
            writer.writerow(row)


def inspect_alignment_for_bsj(
    *,
    alignment_path: Path,
    cell_id: str,
    input_query_count: int,
    output_path: Path,
    qc_output_path: Path | None = None,
    samtools: str = "samtools",
) -> dict[str, Any]:
    if input_query_count < 0:
        raise ValueError("input_query_count must be >= 0")
    if not alignment_path.is_file():
        raise FileNotFoundError(f"Alignment not found: {alignment_path}")
    if alignment_path.suffix.lower() == ".bam":
        samtools_path = _resolve_executable(samtools)
        qc, evidence, _ = _stream_bam_as_sam(
            bam_path=alignment_path,
            samtools=samtools_path,
            input_query_count=input_query_count,
            cell_id=cell_id,
            retained_sam_path=None,
        )
    elif alignment_path.suffix.lower() == ".sam":
        with alignment_path.open("r", encoding="utf-8") as handle:
            qc, evidence = summarize_sam_records(
                handle,
                input_query_count=input_query_count,
                cell_id=cell_id,
            )
    else:
        raise ValueError("Exploratory BSJ inspection accepts only .bam or .sam input")
    write_exploratory_evidence(output_path, evidence)
    if qc_output_path is not None:
        write_json(qc_output_path, qc)
    return qc


def _open_fastq(path: Path) -> TextIO:
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def count_fastq_queries(path: Path) -> int:
    count = 0
    with _open_fastq(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline()
            plus = handle.readline()
            qualities = handle.readline()
            if not sequence or not plus or not qualities:
                raise ValueError(f"Incomplete FASTQ record at query {count + 1}: {path}")
            if not header.startswith("@") or not plus.startswith("+"):
                raise ValueError(f"Malformed FASTQ record at query {count + 1}: {path}")
            if len(sequence.rstrip("\r\n")) != len(qualities.rstrip("\r\n")):
                raise ValueError(
                    f"FASTQ sequence/quality length mismatch at query {count + 1}: {path}"
                )
            count += 1
    return count


def _resolve_executable(value: str) -> str:
    resolved = shutil.which(value)
    if resolved is None:
        candidate = Path(value)
        if candidate.is_file() and os.access(candidate, os.X_OK):
            resolved = str(candidate.resolve())
    if resolved is None:
        raise RuntimeError(f"Required executable not found or not executable: {value}")
    return resolved


def _tool_version(executable: str) -> str:
    result = subprocess.run(
        [executable, "--version"],
        check=False,
        capture_output=True,
        text=True,
        shell=False,
    )
    text = (result.stdout or result.stderr).strip().splitlines()
    if result.returncode != 0 or not text:
        raise RuntimeError(f"Could not determine tool version: {executable}")
    return text[0]


def _reserved_minimap2_option(tokens: Sequence[str]) -> str | None:
    for token in tokens:
        if token in {"-uf", "-u", "-k14", "-k", "-x", "-a"}:
            return token
        if token.startswith("-k") or token.startswith("-x"):
            return token
    return None


def build_minimap2_argv(
    row: LongReadManifestRow,
    *,
    minimap2: str,
    reference_fasta: Path,
    threads: int,
    extra_args: Sequence[str] = (),
) -> list[str]:
    reserved = _reserved_minimap2_option(extra_args)
    if reserved is not None:
        raise ValueError(
            f"minimap2 option {reserved!r} is controlled by molecule_type and cannot "
            "be supplied through extra arguments"
        )
    if threads < 1:
        raise ValueError("threads must be >= 1")
    argv = [minimap2, "-ax", "splice"]
    if row.molecule_type == "direct_rna":
        argv.extend(["-uf", "-k14"])
    elif row.molecule_type != "cdna":
        raise ValueError(f"Unsupported molecule_type: {row.molecule_type!r}")
    argv.extend(["-t", str(threads)])
    argv.extend(extra_args)
    argv.extend([str(reference_fasta), row.long_read_fastq])
    return argv


def _validate_reference(
    reference_fasta: Path,
    *,
    reference_id: str,
    reference_build: str,
    expected_sha256: str,
) -> dict[str, Any]:
    if not reference_id.strip():
        raise ValueError("reference_id is required and must not be inferred from a filename")
    if not reference_build.strip():
        raise ValueError("reference_build is required and must not be inferred from a filename")
    normalized_sha = expected_sha256.strip().lower()
    if not _SHA256.fullmatch(normalized_sha):
        raise ValueError("reference_sha256 must be an explicit 64-character hexadecimal SHA-256")
    if not reference_fasta.is_file():
        raise FileNotFoundError(f"Reference FASTA not found: {reference_fasta}")
    observed_sha = sha256_file(reference_fasta)
    if observed_sha != normalized_sha:
        raise ValueError(
            f"Reference FASTA SHA-256 mismatch: expected {normalized_sha}, observed {observed_sha}"
        )
    contigs: list[str] = []
    with reference_fasta.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                contig = line[1:].strip().split()[0]
                if not contig:
                    raise ValueError(f"Reference FASTA contains an empty header: {reference_fasta}")
                contigs.append(contig)
    if not contigs:
        raise ValueError(f"Reference FASTA contains no contigs: {reference_fasta}")
    return {
        "reference_id": reference_id.strip(),
        "reference_build": reference_build.strip(),
        "fasta_path": str(reference_fasta.resolve()),
        "fasta_sha256": observed_sha,
        "fasta_bytes": reference_fasta.stat().st_size,
        "contigs": contigs,
    }


def _run_alignment_stream(
    *,
    minimap2_argv: list[str],
    samtools_sort_argv: list[str],
    temporary_bam: Path,
    minimap2_log: Path,
    samtools_log: Path,
) -> None:
    temporary_bam.unlink(missing_ok=True)
    with minimap2_log.open("w", encoding="utf-8") as minimap2_stderr, samtools_log.open(
        "w", encoding="utf-8"
    ) as samtools_stderr:
        minimap2_process = subprocess.Popen(
            minimap2_argv,
            stdout=subprocess.PIPE,
            stderr=minimap2_stderr,
            shell=False,
        )
        assert minimap2_process.stdout is not None
        samtools_process = subprocess.Popen(
            samtools_sort_argv,
            stdin=minimap2_process.stdout,
            stdout=subprocess.DEVNULL,
            stderr=samtools_stderr,
            shell=False,
        )
        minimap2_process.stdout.close()
        minimap2_returncode = minimap2_process.wait()
        samtools_returncode = samtools_process.wait()
    if minimap2_returncode != 0 or samtools_returncode != 0:
        temporary_bam.unlink(missing_ok=True)
        raise RuntimeError(
            "Nanopore alignment pipeline failed: "
            f"minimap2={minimap2_returncode}, samtools_sort={samtools_returncode}; "
            f"logs: {minimap2_log}, {samtools_log}"
        )
    if not temporary_bam.is_file():
        raise RuntimeError(
            f"samtools sort exited successfully but did not create BAM: {temporary_bam}"
        )


def _stream_bam_as_sam(
    *,
    bam_path: Path,
    samtools: str,
    input_query_count: int,
    cell_id: str,
    retained_sam_path: Path | None,
) -> tuple[dict[str, Any], list[dict[str, Any]], list[str]]:
    argv = [samtools, "view", "-h", str(bam_path)]
    process = subprocess.Popen(
        argv,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        shell=False,
    )
    assert process.stdout is not None
    retained_handle: TextIO | None = None
    try:
        if retained_sam_path is not None:
            retained_handle = retained_sam_path.open("w", encoding="utf-8")
        qc, evidence = summarize_sam_records(
            process.stdout,
            input_query_count=input_query_count,
            cell_id=cell_id,
            retained_sam=retained_handle,
        )
    finally:
        if retained_handle is not None:
            retained_handle.close()
        process.stdout.close()
    stderr = process.stderr.read() if process.stderr is not None else ""
    returncode = process.wait()
    if returncode != 0:
        if retained_sam_path is not None:
            retained_sam_path.unlink(missing_ok=True)
        raise RuntimeError(
            f"samtools view failed for {cell_id} with exit code {returncode}: {stderr.strip()}"
        )
    return qc, evidence, argv


def _run_checked(argv: list[str], *, log_path: Path) -> None:
    with log_path.open("w", encoding="utf-8") as log:
        result = subprocess.run(
            argv,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
            text=True,
            shell=False,
        )
    if result.returncode != 0:
        raise RuntimeError(
            f"Command failed with exit code {result.returncode}: {shlex.join(argv)}; "
            f"log: {log_path}"
        )


def prepare_nanopore_alignments(
    *,
    manifest_path: Path,
    reference_fasta: Path,
    reference_id: str,
    reference_build: str,
    reference_sha256: str,
    outdir: Path,
    threads: int = 8,
    minimap2: str = "minimap2",
    samtools: str = "samtools",
    keep_sam: bool = False,
    minimap2_extra_args: Sequence[str] = (),
    archive_metadata_path: Path | None = None,
) -> Path:
    rows = read_long_read_manifest_tsv(manifest_path, validate_files=True)
    reference = _validate_reference(
        reference_fasta,
        reference_id=reference_id,
        reference_build=reference_build,
        expected_sha256=reference_sha256,
    )
    minimap2_path = _resolve_executable(minimap2)
    samtools_path = _resolve_executable(samtools)
    minimap2_version = _tool_version(minimap2_path)
    samtools_version = _tool_version(samtools_path)
    archive_metadata: dict[str, Any] | None = None
    if archive_metadata_path is not None:
        archive_metadata = json.loads(archive_metadata_path.read_text(encoding="utf-8"))

    root = outdir / "nanopore_alignment"
    root.mkdir(parents=True, exist_ok=True)
    alignment_manifest_path = root / "alignment_manifest.tsv"
    root_summary_path = root / "nanopore_alignment_summary.json"
    root_provenance_path = root / "provenance.json"
    if (
        alignment_manifest_path.exists()
        or root_summary_path.exists()
        or root_provenance_path.exists()
    ):
        raise FileExistsError(
            f"Nanopore alignment output already exists under {root}; choose a new outdir"
        )

    workflow_uuid = generate_workflow_uuid()
    workflow_started = utc_now_iso()
    workflow_start_clock = time.perf_counter()
    manifest_sha256 = sha256_file(manifest_path)
    alignment_rows: list[AlignmentManifestRow] = []
    cell_summaries: list[dict[str, Any]] = []

    for row in rows:
        cell_dir = root / row.cell_id
        if cell_dir.exists() and any(cell_dir.iterdir()):
            raise FileExistsError(
                f"Cell-scoped output directory is not empty for {row.cell_id}: {cell_dir}"
            )
        cell_dir.mkdir(parents=True, exist_ok=True)
        final_bam = cell_dir / "alignment.bam"
        temporary_bam = cell_dir / "alignment.partial.bam"
        bam_index = cell_dir / "alignment.bam.bai"
        qc_path = cell_dir / "alignment_qc.json"
        evidence_path = cell_dir / "exploratory_bsj_evidence.tsv"
        provenance_path = cell_dir / "provenance.json"
        retained_sam_path = cell_dir / "alignment.sam" if keep_sam else None
        minimap2_log = cell_dir / "minimap2.stderr.log"
        samtools_sort_log = cell_dir / "samtools_sort.stderr.log"
        samtools_index_log = cell_dir / "samtools_index.log"

        input_query_count = count_fastq_queries(Path(row.long_read_fastq))
        minimap2_argv = build_minimap2_argv(
            row,
            minimap2=minimap2_path,
            reference_fasta=reference_fasta.resolve(),
            threads=threads,
            extra_args=minimap2_extra_args,
        )
        samtools_sort_argv = [
            samtools_path,
            "sort",
            "-@",
            str(threads),
            "-o",
            str(temporary_bam),
            "-",
        ]
        cell_started = utc_now_iso()
        cell_start_clock = time.perf_counter()
        _run_alignment_stream(
            minimap2_argv=minimap2_argv,
            samtools_sort_argv=samtools_sort_argv,
            temporary_bam=temporary_bam,
            minimap2_log=minimap2_log,
            samtools_log=samtools_sort_log,
        )
        os.replace(temporary_bam, final_bam)
        samtools_index_argv = [samtools_path, "index", str(final_bam)]
        _run_checked(samtools_index_argv, log_path=samtools_index_log)
        if not bam_index.is_file():
            raise RuntimeError(f"samtools index did not create expected index: {bam_index}")

        qc, evidence, samtools_view_argv = _stream_bam_as_sam(
            bam_path=final_bam,
            samtools=samtools_path,
            input_query_count=input_query_count,
            cell_id=row.cell_id,
            retained_sam_path=retained_sam_path,
        )
        write_json(qc_path, qc)
        write_exploratory_evidence(evidence_path, evidence)
        cell_completed = utc_now_iso()
        cell_elapsed = time.perf_counter() - cell_start_clock
        paths = {
            "bam": str(final_bam.resolve()),
            "bam_index": str(bam_index.resolve()),
            "alignment_qc": str(qc_path.resolve()),
            "exploratory_bsj_evidence": str(evidence_path.resolve()),
            "retained_sam": str(retained_sam_path.resolve()) if retained_sam_path else None,
            "minimap2_log": str(minimap2_log.resolve()),
            "samtools_sort_log": str(samtools_sort_log.resolve()),
            "samtools_index_log": str(samtools_index_log.resolve()),
            "provenance": str(provenance_path.resolve()),
        }
        provenance = {
            "schema_version": "circyto.nanopore_alignment_provenance.v1",
            "cell_id": row.cell_id,
            "dataset_id": row.dataset_id,
            "source_accession": row.source_accession,
            "manifest": {
                "path": str(manifest_path.resolve()),
                "sha256": manifest_sha256,
            },
            "input": {
                "long_read_fastq": row.long_read_fastq,
                "fastq_sha256": sha256_file(Path(row.long_read_fastq)),
                "input_query_count": input_query_count,
            },
            "archive_metadata": archive_metadata,
            "sequencing_platform": row.sequencing_platform,
            "protocol": row.protocol,
            "archive_library_selection": row.archive_library_selection,
            "library_preparation_summary": row.library_preparation_summary,
            "molecule_type": row.molecule_type,
            "barcode_status": row.barcode_status,
            "biological_interpretation_boundary": row.biological_interpretation_boundary,
            "reference": reference,
            "tools": {
                "minimap2": {
                    "path": minimap2_path,
                    "version": minimap2_version,
                },
                "samtools": {
                    "path": samtools_path,
                    "version": samtools_version,
                },
            },
            "commands": {
                "minimap2_argv": minimap2_argv,
                "minimap2_display": shlex.join(minimap2_argv),
                "samtools_sort_argv": samtools_sort_argv,
                "samtools_sort_display": shlex.join(samtools_sort_argv),
                "samtools_index_argv": samtools_index_argv,
                "samtools_index_display": shlex.join(samtools_index_argv),
                "samtools_view_argv": samtools_view_argv,
                "samtools_view_display": shlex.join(samtools_view_argv),
                "shell": False,
                "streamed_minimap2_to_samtools_sort": True,
            },
            "alignment_qc": qc,
            "paths": paths,
            "detector_invoked": False,
            "detector_backend": None,
            "circRNA_validation_status": False,
            "scientific_disclaimers": [
                SCIENTIFIC_BOUNDARY_WARNING,
                EXPLORATORY_EVIDENCE_WARNING,
                NO_DETECTOR_WARNING,
            ],
            "stage_graph": [
                {"stage": "minimap2_alignment", "status": "completed"},
                {"stage": "coordinate_sort", "status": "completed"},
                {"stage": "bam_index", "status": "completed"},
                {"stage": "alignment_qc", "status": "completed"},
                {"stage": "exploratory_bsj_evidence", "status": "completed"},
                {"stage": "circRNA_detector", "status": "disabled"},
            ],
        }
        provenance = apply_standard_provenance(
            provenance,
            command_name="circyto nanopore align",
            workflow_type="experimental_nanopore_interoperability",
            protocol=row.protocol,
            read_layout="single-end",
            genome_fasta=str(reference_fasta.resolve()),
            gtf=None,
            detector_backend=None,
            started_at=cell_started,
            completed_at=cell_completed,
            elapsed_seconds=cell_elapsed,
            workflow_uuid=workflow_uuid,
        )
        write_json(provenance_path, provenance)
        alignment_rows.append(
            AlignmentManifestRow(
                cell_id=row.cell_id,
                bam=str(final_bam.resolve()),
                group_id=row.dataset_id,
                read_layout="single-end",
                aligner="minimap2",
                reference=str(reference_fasta.resolve()),
                cache_key="",
                source_manifest=str(manifest_path.resolve()),
                mapper_mode="splice",
                artifact_bucket="nanopore_alignment",
                sortedness="coordinate",
                extra={
                    "bam_index": str(bam_index.resolve()),
                    "alignment_qc": str(qc_path.resolve()),
                    "provenance_json": str(provenance_path.resolve()),
                    "exploratory_bsj_evidence": str(evidence_path.resolve()),
                    "retained_sam": (
                        str(retained_sam_path.resolve()) if retained_sam_path else ""
                    ),
                    "source_accession": row.source_accession,
                    "dataset_id": row.dataset_id,
                    "molecule_type": row.molecule_type,
                    "circRNA_validation_status": "false",
                    "biological_interpretation_boundary": (
                        row.biological_interpretation_boundary
                    ),
                },
            )
        )
        cell_summaries.append(
            {
                "cell_id": row.cell_id,
                "cell_directory": str(cell_dir.resolve()),
                "paths": paths,
                "alignment_qc": qc,
                "commands": provenance["commands"],
            }
        )

    write_alignment_manifest_tsv(alignment_rows, alignment_manifest_path)
    workflow_completed = utc_now_iso()
    summary = {
        "schema_version": "circyto.nanopore_alignment_summary.v1",
        "workflow_uuid": workflow_uuid,
        "workflow_type": "experimental_nanopore_interoperability",
        "started_at": workflow_started,
        "completed_at": workflow_completed,
        "elapsed_seconds": round(time.perf_counter() - workflow_start_clock, 3),
        "manifest": str(manifest_path.resolve()),
        "alignment_manifest": str(alignment_manifest_path.resolve()),
        "run_provenance": str(root_provenance_path.resolve()),
        "reference": reference,
        "cell_count": len(cell_summaries),
        "cells": cell_summaries,
        "detector_invoked": False,
        "circRNA_validation_status": False,
        "scientific_disclaimers": [
            SCIENTIFIC_BOUNDARY_WARNING,
            EXPLORATORY_EVIDENCE_WARNING,
            NO_DETECTOR_WARNING,
        ],
    }
    write_json(root_summary_path, summary)
    run_provenance = {
        "schema_version": "circyto.nanopore_alignment_run_provenance.v1",
        "workflow_uuid": workflow_uuid,
        "dataset_ids": sorted({row.dataset_id for row in rows}),
        "source_accessions": sorted({row.source_accession for row in rows}),
        "biological_interpretation_boundaries": sorted(
            {row.biological_interpretation_boundary for row in rows}
        ),
        "manifest": {
            "path": str(manifest_path.resolve()),
            "sha256": manifest_sha256,
        },
        "archive_metadata": archive_metadata,
        "reference": reference,
        "tools": {
            "minimap2": {
                "path": minimap2_path,
                "version": minimap2_version,
            },
            "samtools": {
                "path": samtools_path,
                "version": samtools_version,
            },
        },
        "cells": cell_summaries,
        "paths": {
            "alignment_manifest": str(alignment_manifest_path.resolve()),
            "alignment_summary": str(root_summary_path.resolve()),
            "run_provenance": str(root_provenance_path.resolve()),
        },
        "detector_invoked": False,
        "detector_backend": None,
        "circRNA_validation_status": False,
        "scientific_disclaimers": [
            SCIENTIFIC_BOUNDARY_WARNING,
            EXPLORATORY_EVIDENCE_WARNING,
            NO_DETECTOR_WARNING,
        ],
        "stage_graph": [
            {"stage": "minimap2_alignment", "status": "completed"},
            {"stage": "coordinate_sort", "status": "completed"},
            {"stage": "bam_index", "status": "completed"},
            {"stage": "alignment_qc", "status": "completed"},
            {"stage": "exploratory_bsj_evidence", "status": "completed"},
            {"stage": "circRNA_detector", "status": "disabled"},
        ],
    }
    run_provenance = apply_standard_provenance(
        run_provenance,
        command_name="circyto nanopore align",
        workflow_type="experimental_nanopore_interoperability",
        protocol=";".join(sorted({row.protocol for row in rows})),
        read_layout="single-end",
        genome_fasta=str(reference_fasta.resolve()),
        gtf=None,
        detector_backend=None,
        started_at=workflow_started,
        completed_at=workflow_completed,
        elapsed_seconds=time.perf_counter() - workflow_start_clock,
        workflow_uuid=workflow_uuid,
    )
    write_json(root_provenance_path, run_provenance)
    return alignment_manifest_path
