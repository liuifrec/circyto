from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import zip_longest
from pathlib import Path
from typing import Dict, Iterable, Optional, TextIO, Tuple
import csv
import gzip
import json

from circyto.manifest.v1 import ManifestRow, write_manifest_tsv


def _open_text_maybe_gz(path: Path, mode: str = "rt") -> TextIO:
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return path.open(mode, encoding="utf-8")


def _detect_delimiter(path: Path) -> str:
    if path.suffix.lower() == ".csv":
        return ","
    sample = path.read_text(encoding="utf-8", errors="replace")[:4096]
    if not sample.strip():
        return "\t"
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
        return dialect.delimiter
    except csv.Error:
        return "\t"


def _fastq_iter_strict(handle: TextIO, *, path: Path) -> Iterable[Tuple[str, str, str, str]]:
    record_no = 0
    while True:
        header = handle.readline()
        if not header:
            return
        record_no += 1
        seq = handle.readline()
        plus = handle.readline()
        qual = handle.readline()
        if not seq or not plus or not qual:
            raise ValueError(f"Malformed FASTQ in {path} near record {record_no}: truncated 4-line record")
        header = header.rstrip("\n\r")
        seq = seq.rstrip("\n\r")
        plus = plus.rstrip("\n\r")
        qual = qual.rstrip("\n\r")
        if not header.startswith("@"):
            raise ValueError(f"Malformed FASTQ in {path} near record {record_no}: header does not start with '@'")
        if not plus.startswith("+"):
            raise ValueError(f"Malformed FASTQ in {path} near record {record_no}: plus line does not start with '+'")
        if len(seq) != len(qual):
            raise ValueError(
                f"Malformed FASTQ in {path} near record {record_no}: sequence/quality lengths differ ({len(seq)} != {len(qual)})"
            )
        yield header, seq, plus, qual


def _normalize_read_id(header: str) -> str:
    token = header.strip().split()[0]
    token = token[1:] if token.startswith("@") else token
    if token.endswith("/1") or token.endswith("/2"):
        token = token[:-2]
    return token


def _sanitize_header(header: str) -> str:
    token = header.strip().split()[0]
    if not token.startswith("@"):
        token = "@" + token.lstrip("@")
    return token


def _write_fastq_record(handle: TextIO, record: Tuple[str, str, str, str]) -> None:
    header, seq, _, qual = record
    handle.write(f"{_sanitize_header(header)}\n{seq}\n+\n{qual}\n")


def _hamming(a: str, b: str) -> int:
    if len(a) != len(b):
        raise ValueError("Hamming distance requires equal-length strings")
    return sum(ch1 != ch2 for ch1, ch2 in zip(a, b))


def _match_index_pair(
    index1: str,
    index2: str,
    *,
    pair_to_cell: Dict[Tuple[str, str], str],
    max_mismatch: int,
) -> Optional[str]:
    exact = pair_to_cell.get((index1, index2))
    if exact is not None or max_mismatch <= 0:
        return exact

    best_cell: Optional[str] = None
    best_distance: Optional[int] = None
    tied = False
    for (candidate_i1, candidate_i2), cell_id in pair_to_cell.items():
        if len(candidate_i1) != len(index1) or len(candidate_i2) != len(index2):
            continue
        distance = _hamming(candidate_i1, index1) + _hamming(candidate_i2, index2)
        if distance > max_mismatch:
            continue
        if best_distance is None or distance < best_distance:
            best_cell = cell_id
            best_distance = distance
            tied = False
        elif distance == best_distance:
            tied = True
    if tied:
        return None
    return best_cell


def load_smartseq3_annotation(
    path: Path,
    *,
    cell_id_column: str,
    index1_column: str,
    index2_column: str,
) -> Dict[Tuple[str, str], str]:
    delimiter = _detect_delimiter(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        columns = list(reader.fieldnames or [])
        required = [cell_id_column, index1_column, index2_column]
        missing = [column for column in required if column not in columns]
        if missing:
            available = ", ".join(columns) if columns else "<none>"
            missing_text = ", ".join(missing)
            raise ValueError(
                f"Annotation file {path} is missing requested column(s): {missing_text}. Available columns: {available}"
            )

        mapping: Dict[Tuple[str, str], str] = {}
        for line_no, row in enumerate(reader, start=2):
            cell_id = (row.get(cell_id_column) or "").strip()
            i1 = (row.get(index1_column) or "").strip().upper()
            i2 = (row.get(index2_column) or "").strip().upper()
            if not cell_id:
                raise ValueError(f"Annotation file {path}:{line_no} has empty cell ID in column '{cell_id_column}'")
            if not i1 or not i2:
                raise ValueError(
                    f"Annotation file {path}:{line_no} has empty index value in columns '{index1_column}'/'{index2_column}'"
                )
            key = (i1, i2)
            existing = mapping.get(key)
            if existing is not None and existing != cell_id:
                raise ValueError(
                    f"Duplicate SMART-Seq3 index pair {i1}+{i2} maps to multiple cell IDs: {existing}, {cell_id}"
                )
            mapping[key] = cell_id

    if not mapping:
        raise ValueError(f"Annotation file {path} contains no data rows")
    return mapping


@dataclass(frozen=True)
class SmartSeq3DemuxParams:
    read1: Path
    read2: Path
    index1: Path
    index2: Path
    annotation: Path
    outdir: Path
    cell_id_column: str
    index1_column: str
    index2_column: str
    max_mismatch: int = 0
    max_records: Optional[int] = None
    write_sink: bool = True
    compress_output: bool = True
    emit_manifest: bool = True


def demux_smartseq3(params: SmartSeq3DemuxParams) -> dict[str, object]:
    if params.max_mismatch < 0:
        raise ValueError("max_mismatch must be >= 0")
    if params.max_records is not None and params.max_records <= 0:
        raise ValueError("max_records must be > 0 when provided")

    fastq_dir = params.outdir / "fastq"
    sink_dir = params.outdir / "sink"
    fastq_dir.mkdir(parents=True, exist_ok=True)
    if params.write_sink:
        sink_dir.mkdir(parents=True, exist_ok=True)

    suffix = ".fastq.gz" if params.compress_output else ".fastq"
    pair_to_cell = load_smartseq3_annotation(
        params.annotation,
        cell_id_column=params.cell_id_column,
        index1_column=params.index1_column,
        index2_column=params.index2_column,
    )

    writers: Dict[str, Tuple[TextIO, TextIO]] = {}
    reads_per_cell: Counter[str] = Counter()
    unmatched_pairs: Counter[str] = Counter()
    total_records = 0
    assigned_records = 0
    unmatched_records = 0
    malformed_records = 0

    sink_r1: Optional[TextIO] = None
    sink_r2: Optional[TextIO] = None

    def get_sink_writers() -> Tuple[TextIO, TextIO]:
        nonlocal sink_r1, sink_r2
        if sink_r1 is not None and sink_r2 is not None:
            return sink_r1, sink_r2
        if params.compress_output:
            sink_r1 = gzip.open(sink_dir / f"unmatched_R1{suffix}", "wt")
            sink_r2 = gzip.open(sink_dir / f"unmatched_R2{suffix}", "wt")
        else:
            sink_r1 = (sink_dir / f"unmatched_R1{suffix}").open("w", encoding="utf-8")
            sink_r2 = (sink_dir / f"unmatched_R2{suffix}").open("w", encoding="utf-8")
        return sink_r1, sink_r2

    def get_writers(cell_id: str) -> Tuple[TextIO, TextIO]:
        existing = writers.get(cell_id)
        if existing is not None:
            return existing
        r1_path = fastq_dir / f"{cell_id}_R1{suffix}"
        r2_path = fastq_dir / f"{cell_id}_R2{suffix}"
        if params.compress_output:
            handle1 = gzip.open(r1_path, "wt")
            handle2 = gzip.open(r2_path, "wt")
        else:
            handle1 = r1_path.open("w", encoding="utf-8")
            handle2 = r2_path.open("w", encoding="utf-8")
        writers[cell_id] = (handle1, handle2)
        return handle1, handle2

    try:
        with (
            _open_text_maybe_gz(params.read1, "rt") as read1_handle,
            _open_text_maybe_gz(params.read2, "rt") as read2_handle,
            _open_text_maybe_gz(params.index1, "rt") as index1_handle,
            _open_text_maybe_gz(params.index2, "rt") as index2_handle,
        ):
            read1_iter = _fastq_iter_strict(read1_handle, path=params.read1)
            read2_iter = _fastq_iter_strict(read2_handle, path=params.read2)
            index1_iter = _fastq_iter_strict(index1_handle, path=params.index1)
            index2_iter = _fastq_iter_strict(index2_handle, path=params.index2)

            for record_no, records in enumerate(
                zip_longest(read1_iter, read2_iter, index1_iter, index2_iter),
                start=1,
            ):
                if params.max_records is not None and total_records >= params.max_records:
                    break
                if any(record is None for record in records):
                    raise ValueError(
                        "FASTQ inputs have different record counts near record "
                        f"{record_no}: read1={params.read1} read2={params.read2} index1={params.index1} index2={params.index2}"
                    )

                read1_record, read2_record, index1_record, index2_record = records
                assert read1_record is not None
                assert read2_record is not None
                assert index1_record is not None
                assert index2_record is not None

                read_ids = {
                    _normalize_read_id(read1_record[0]),
                    _normalize_read_id(read2_record[0]),
                    _normalize_read_id(index1_record[0]),
                    _normalize_read_id(index2_record[0]),
                }
                if len(read_ids) != 1:
                    raise ValueError(
                        "FASTQ read ID mismatch near record "
                        f"{record_no}: R1={read1_record[0]} R2={read2_record[0]} I1={index1_record[0]} I2={index2_record[0]}"
                    )

                total_records += 1
                i1_seq = index1_record[1].upper()
                i2_seq = index2_record[1].upper()
                cell_id = _match_index_pair(
                    i1_seq,
                    i2_seq,
                    pair_to_cell=pair_to_cell,
                    max_mismatch=params.max_mismatch,
                )

                if cell_id is None:
                    unmatched_records += 1
                    unmatched_pairs[f"{i1_seq}+{i2_seq}"] += 1
                    if params.write_sink:
                        sink_r1, sink_r2 = get_sink_writers()
                        _write_fastq_record(sink_r1, read1_record)
                        _write_fastq_record(sink_r2, read2_record)
                    continue

                assigned_records += 1
                reads_per_cell[cell_id] += 1
                writer_r1, writer_r2 = get_writers(cell_id)
                _write_fastq_record(writer_r1, read1_record)
                _write_fastq_record(writer_r2, read2_record)
    except ValueError as exc:
        if "Malformed FASTQ" in str(exc):
            malformed_records += 1
        raise
    finally:
        for writer_r1, writer_r2 in writers.values():
            writer_r1.close()
            writer_r2.close()
        if sink_r1 is not None:
            sink_r1.close()
        if sink_r2 is not None:
            sink_r2.close()

    manifest_path = params.outdir / "manifest.tsv"
    if params.emit_manifest:
        rows = [
            ManifestRow(
                cell_id=cell_id,
                platform="smartseq3",
                read1=str((fastq_dir / f"{cell_id}_R1{suffix}").resolve()),
                read2=str((fastq_dir / f"{cell_id}_R2{suffix}").resolve()),
                bam="",
                library_id="smartseq3_demux",
                n_input_reads=count,
                extra={"read_layout": "paired-end"},
            )
            for cell_id, count in sorted(reads_per_cell.items())
            if count > 0
        ]
        write_manifest_tsv(rows, manifest_path)

    summary = {
        "total_records": total_records,
        "assigned_records": assigned_records,
        "unmatched_records": unmatched_records,
        "malformed_records": malformed_records,
        "assignment_rate": (assigned_records / total_records) if total_records else 0.0,
        "number_of_cells_detected": len(reads_per_cell),
        "reads_per_cell": {cell_id: reads_per_cell[cell_id] for cell_id in sorted(reads_per_cell)},
        "top_unmatched_index_pairs": [
            {"index_pair": key, "count": count} for key, count in unmatched_pairs.most_common(10)
        ],
        "input_paths": {
            "read1": str(params.read1),
            "read2": str(params.read2),
            "index1": str(params.index1),
            "index2": str(params.index2),
        },
        "annotation_path": str(params.annotation),
        "max_records": params.max_records,
        "max_mismatch": params.max_mismatch,
    }
    (params.outdir / "demux_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return summary
