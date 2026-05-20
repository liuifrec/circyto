from __future__ import annotations

import csv
import os
import json
import re
import time
from bisect import bisect_left, bisect_right
from collections import defaultdict
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterator

import numpy as np
import pandas as pd

from circyto.manifest.alignment import read_alignment_manifest_tsv
from circyto.pipeline.workflow_reporting import (
    apply_standard_provenance,
    cleanup_summary_block,
    load_circ_matrix,
    load_json,
    read_index_lines,
    summarize_read_layouts,
    utc_now_iso,
    write_json,
)
from circyto.pipeline.workflow_integrity import check_workflow_integrity


VALID_GENE_EXPRESSION_METHODS = {"none", "simple-overlap", "featurecounts", "velocyto"}
VALID_VELOCITY_LAYERS = {"none", "velocyto"}
SUPPORTED_CLEANUP_SCOPES = ("alignments", "demux", "all")
REGENERABLE_ALIGNMENT_SUFFIXES = (
    ".sam",
    ".bam",
    ".bai",
    ".junction",
)
REGENERABLE_ALIGNMENT_NAMES = (
    "Aligned.out.sam",
    "Chimeric.out.junction",
    "Unmapped.out.mate1",
    "Unmapped.out.mate2",
    "bwa_rescue.sam",
)
REGENERABLE_DEMUX_SUFFIXES = (
    ".fastq",
    ".fq",
    ".fastq.gz",
    ".fq.gz",
)
REGENERABLE_CHUNK_SUFFIXES = (
    ".chunk",
    ".chunks",
    ".part",
    ".tmp",
)
MUST_KEEP_OUTPUTS = [
    "workflow_summary.json",
    "detector summaries",
    "matrix/",
    "anndata/",
    "mudata/",
    "qc/",
    "manifests/",
    "logs/",
    "provenance JSON",
    "final detector TSVs",
]
USER_INPUTS_NEVER_DELETE = [
    "original raw FASTQs",
    "original pooled Smart-seq FASTQs",
    "user-supplied manifests",
    "user-supplied gene_counts.tsv",
]
SAFE_TO_DELETE_ALIGNMENTS = [
    "align/cache/",
    "large *.sam / *.bam / *.bai intermediates",
    "STAR temporary outputs",
    "BWA rescue SAM/BAM",
    "temporary chunk files if reproducible",
]
SAFE_TO_DELETE_DEMUX = [
    "demux per-cell FASTQs generated from pooled Smart-seq2/3 input",
    "temporary chunk files if reproducible",
]

VALID_VELOCITY_LAYOUT_FILES = {
    "barcodes.tsv",
    "features.tsv",
    "spliced.mtx",
    "unspliced.mtx",
}


def validate_full_length_future_options(
    *,
    export_mudata: bool,
    gene_counts: Path | None,
    gene_expression_method: str,
    velocity_layers: str,
    cleanup_intermediates: str | None,
    dry_run: bool,
) -> None:
    if gene_expression_method not in VALID_GENE_EXPRESSION_METHODS:
        raise ValueError(
            f"Unsupported --gene-expression-method: {gene_expression_method}. "
            f"Choose from: {', '.join(sorted(VALID_GENE_EXPRESSION_METHODS))}"
        )
    if velocity_layers not in VALID_VELOCITY_LAYERS:
        raise ValueError(
            f"Unsupported --velocity-layers: {velocity_layers}. "
            f"Choose from: {', '.join(sorted(VALID_VELOCITY_LAYERS))}"
        )
    if velocity_layers == "velocyto" and gene_expression_method != "velocyto":
        raise ValueError("--velocity-layers velocyto requires --gene-expression-method velocyto")

    if gene_counts is not None and gene_expression_method == "velocyto":
        raise NotImplementedError(
            "--gene-counts TSV import currently supports only normalized gene-count table staging, "
            "not velocyto-based RNA layer import."
        )

    if gene_counts is not None and gene_expression_method == "simple-overlap":
        raise ValueError(
            "--gene-counts and --gene-expression-method simple-overlap are mutually exclusive. "
            "Provide either an external gene_counts.tsv import or request internal simple-overlap counting."
        )

    if export_mudata:
        raise NotImplementedError(
            "Gene-expression and velocity-compatible outputs are not implemented yet for "
            "`circyto workflow full-length-circrna`. Current behavior remains circ-only h5ad. "
            "See docs/gene_expression_velocity_integration.md for the staged design."
        )

    if gene_expression_method in {"featurecounts", "velocyto"} or velocity_layers != "none":
        raise NotImplementedError(
            "featureCounts/velocyto-based gene-expression and velocity-compatible outputs are not implemented yet for "
            "`circyto workflow full-length-circrna`. Current production behavior remains circ-only h5ad plus optional "
            "external gene-count import or internal simple-overlap RNA sanity counting. "
            "See docs/gene_expression_velocity_integration.md for the staged design."
        )

    if cleanup_intermediates is not None and cleanup_intermediates not in SUPPORTED_CLEANUP_SCOPES:
        raise ValueError(
            f"Unsupported --cleanup-intermediates scope: {cleanup_intermediates}. "
            f"Choose from: {', '.join(SUPPORTED_CLEANUP_SCOPES)}"
        )


def validate_feature_id_uniqueness(feature_ids: list[str], *, label: str) -> None:
    normalized = [str(feature_id).strip() for feature_id in feature_ids]
    duplicates = sorted({feature_id for feature_id in normalized if feature_id and normalized.count(feature_id) > 1})
    if duplicates:
        raise ValueError(
            f"Duplicate {label} IDs detected: {', '.join(duplicates[:5])}. "
            f"{label} identifiers must be unique."
        )


def validate_cell_id_consistency(
    circ_cell_ids: list[str],
    rna_cell_ids: list[str],
    *,
    circ_label: str = "circ",
    rna_label: str = "rna",
) -> None:
    circ = [str(cell_id).strip() for cell_id in circ_cell_ids if str(cell_id).strip()]
    rna = [str(cell_id).strip() for cell_id in rna_cell_ids if str(cell_id).strip()]
    if not circ:
        raise ValueError(f"{circ_label} cell ID set is empty")
    if not rna:
        raise ValueError(f"{rna_label} cell ID set is empty")
    if circ != rna:
        only_circ = sorted(set(circ) - set(rna))
        only_rna = sorted(set(rna) - set(circ))
        raise ValueError(
            f"Cell ID mismatch between {circ_label} and {rna_label}. "
            f"Only in {circ_label}: {only_circ[:5]}. Only in {rna_label}: {only_rna[:5]}."
        )


def validate_gene_expression_table_schema(path: Path) -> dict[str, Any]:
    df = pd.read_csv(path, sep="\t", keep_default_na=False)
    required = ["gene_id", "gene_name"]
    missing = [column for column in required if column not in df.columns]
    if missing:
        raise ValueError(
            f"{path} is missing required gene-expression columns: {', '.join(missing)}. "
            "Expected at minimum: gene_id, gene_name, and one or more cell columns."
        )
    if df.empty:
        raise ValueError(f"{path} has no gene-expression rows")
    cell_columns = [column for column in df.columns if column not in {"gene_id", "gene_name"}]
    if not cell_columns:
        raise ValueError(f"{path} has no cell columns after gene_id and gene_name")
    validate_feature_id_uniqueness(df["gene_id"].astype(str).tolist(), label="gene")
    validate_feature_id_uniqueness([str(column) for column in cell_columns], label="cell")
    numeric = df[cell_columns].apply(pd.to_numeric, errors="coerce")
    if numeric.isna().any().any():
        raise ValueError(
            f"{path} contains non-numeric gene-count values in one or more cell columns."
        )
    return {
        "path": str(path.resolve()),
        "n_genes": int(df.shape[0]),
        "n_cells": int(len(cell_columns)),
        "cell_ids": [str(column) for column in cell_columns],
        "feature_id_column": "gene_id",
        "feature_name_column": "gene_name",
    }


def import_gene_counts_table(
    *,
    path: Path,
    expected_cell_ids: list[str],
    outdir: Path,
) -> dict[str, Any]:
    summary = validate_gene_expression_table_schema(path)
    df = pd.read_csv(path, sep="\t", keep_default_na=False)
    cell_columns = [column for column in df.columns if column not in {"gene_id", "gene_name"}]
    validate_cell_id_consistency(expected_cell_ids, [str(column) for column in cell_columns], circ_label="circ", rna_label="rna")

    numeric = df[cell_columns].apply(pd.to_numeric, errors="raise")
    normalized_counts = pd.concat([df[["gene_id", "gene_name"]].copy(), numeric], axis=1)
    feature_table = df[["gene_id", "gene_name"]].copy()

    outdir.mkdir(parents=True, exist_ok=True)
    gene_counts_out = outdir / "gene_counts.tsv"
    feature_out = outdir / "gene_feature_table.tsv"
    summary_out = outdir / "rna_import_summary.json"
    normalized_counts.to_csv(gene_counts_out, sep="\t", index=False)
    feature_table.to_csv(feature_out, sep="\t", index=False)

    payload = {
        **summary,
        "cell_ids": [str(column) for column in cell_columns],
        "output_gene_counts": str(gene_counts_out.resolve()),
        "output_gene_feature_table": str(feature_out.resolve()),
        "feature_table_columns": ["gene_id", "gene_name"],
        "count_table_columns": ["gene_id", "gene_name", *[str(column) for column in cell_columns]],
        "total_counts_sum": float(numeric.to_numpy().sum()),
    }
    payload.update(
        _write_rna_qc_outputs(
            counts_df=normalized_counts,
            feature_df=feature_table,
            work_root=outdir.parent,
        )
    )
    summary_out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    payload["output_rna_import_summary"] = str(summary_out.resolve())
    return payload


def _parse_gtf_attributes(text: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for match in re.finditer(r'([A-Za-z0-9_:-]+)\s+"([^"]*)"', text):
        attrs[match.group(1)] = match.group(2)
    return attrs


def _load_gene_features_from_gtf(path: Path) -> list[dict[str, Any]]:
    genes: list[dict[str, Any]] = []
    exon_spans: dict[str, dict[str, Any]] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            attrs = _parse_gtf_attributes(fields[8])
            gene_id = attrs.get("gene_id", "").strip()
            gene_biotype = (
                attrs.get("gene_biotype", "").strip()
                or attrs.get("gene_type", "").strip()
                or attrs.get("biotype", "").strip()
            )
            feature_type = fields[2]
            if feature_type == "gene":
                if not gene_id:
                    continue
                genes.append(
                    {
                        "gene_id": gene_id,
                        "gene_name": attrs.get("gene_name", gene_id).strip() or gene_id,
                        "chrom": fields[0].strip(),
                        "start": int(fields[3]),
                        "end": int(fields[4]),
                        "strand": fields[6].strip(),
                        "gene_biotype": gene_biotype,
                    }
                )
            elif feature_type == "exon" and gene_id:
                gene_name = attrs.get("gene_name", gene_id).strip() or gene_id
                start = int(fields[3])
                end = int(fields[4])
                existing = exon_spans.get(gene_id)
                if existing is None:
                    exon_spans[gene_id] = {
                        "gene_id": gene_id,
                        "gene_name": gene_name,
                        "chrom": fields[0].strip(),
                        "start": start,
                        "end": end,
                        "strand": fields[6].strip(),
                        "gene_biotype": gene_biotype,
                    }
                else:
                    existing["start"] = min(int(existing["start"]), start)
                    existing["end"] = max(int(existing["end"]), end)
                    if not str(existing.get("gene_biotype", "")).strip() and gene_biotype:
                        existing["gene_biotype"] = gene_biotype
    if not genes:
        genes = list(exon_spans.values())
    if not genes:
        raise ValueError(f"{path} contains no gene or exon features with gene_id")
    validate_feature_id_uniqueness([item["gene_id"] for item in genes], label="gene")
    return genes


def _gtf_cache_key(path: Path) -> tuple[str, int, int]:
    stat = path.stat()
    return (str(path.resolve()), int(stat.st_mtime_ns), int(stat.st_size))


@lru_cache(maxsize=8)
def _load_cached_gene_index(cache_key: tuple[str, int, int]) -> tuple[list[dict[str, Any]], dict[str, dict[str, Any]]]:
    path = Path(cache_key[0])
    genes = _load_gene_features_from_gtf(path)
    return genes, _gene_index_by_chrom(genes)


def _gene_index_by_chrom(genes: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    by_chrom: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for gene in genes:
        by_chrom[str(gene["chrom"])].append(gene)
    indexed: dict[str, dict[str, Any]] = {}
    for chrom, chrom_genes in by_chrom.items():
        ordered = sorted(chrom_genes, key=lambda item: (int(item["start"]), int(item["end"]), str(item["gene_id"])))
        starts = [int(item["start"]) for item in ordered]
        prefix_max_ends: list[int] = []
        running_max = -1
        for item in ordered:
            running_max = max(running_max, int(item["end"]))
            prefix_max_ends.append(running_max)
        indexed[chrom] = {
            "genes": ordered,
            "starts": starts,
            "prefix_max_ends": prefix_max_ends,
        }
    return indexed


def _cigar_reference_blocks(pos: int, cigar: str) -> list[tuple[int, int]]:
    if not cigar or cigar == "*":
        return []
    blocks: list[tuple[int, int]] = []
    ref_pos = int(pos)
    for length_text, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        length = int(length_text)
        if op in {"M", "=", "X"}:
            start = ref_pos
            end = ref_pos + length - 1
            blocks.append((start, end))
            ref_pos += length
        elif op in {"D", "N"}:
            ref_pos += length
        else:
            continue
    return blocks


def _overlapping_gene_ids(
    *,
    chrom: str,
    blocks: list[tuple[int, int]],
    genes_by_chrom: dict[str, dict[str, Any]],
) -> set[str]:
    overlaps: set[str] = set()
    chrom_index = genes_by_chrom.get(chrom)
    if chrom_index is None:
        return overlaps
    genes = chrom_index["genes"]
    starts = chrom_index["starts"]
    prefix_max_ends = chrom_index["prefix_max_ends"]
    for block_start, block_end in blocks:
        hi = bisect_right(starts, int(block_end))
        if hi <= 0:
            continue
        lo = bisect_left(prefix_max_ends, int(block_start), 0, hi)
        for gene in genes[lo:hi]:
            gene_start = int(gene["start"])
            gene_end = int(gene["end"])
            if gene_end < block_start or gene_start > block_end:
                continue
            overlaps.add(str(gene["gene_id"]))
    return overlaps


def _iter_sam_templates(path: Path) -> Iterator[set[tuple[str, int, str]]]:
    current_qname: str | None = None
    current_spans: set[tuple[str, int, str]] = set()
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("@"):
                continue
            if len(row) < 6:
                continue
            qname = row[0].strip()
            flag = int(row[1])
            rname = row[2].strip()
            pos = int(row[3])
            cigar = row[5].strip()
            if not qname or rname in {"", "*"} or cigar in {"", "*"}:
                continue
            if flag & 0x4 or flag & 0x100 or flag & 0x800 or flag & 0x200 or flag & 0x400:
                continue
            if current_qname is None:
                current_qname = qname
            if qname != current_qname:
                if current_spans:
                    yield current_spans
                current_qname = qname
                current_spans = set()
            current_spans.add((rname, pos, cigar))
    if current_spans:
        yield current_spans


def _iter_bam_templates(path: Path) -> Iterator[set[tuple[str, int, str]]]:
    try:
        import pysam  # type: ignore
    except Exception as exc:  # pragma: no cover - optional dependency path
        raise NotImplementedError(
            "simple-overlap gene counting on BAM inputs requires optional pysam, which is not installed."
        ) from exc
    current_qname: str | None = None
    current_spans: set[tuple[str, int, str]] = set()
    with pysam.AlignmentFile(str(path), "rb") as handle:  # pragma: no cover - optional dependency path
        for record in handle.fetch(until_eof=True):
            if (
                record.is_unmapped
                or record.is_secondary
                or record.is_supplementary
                or record.is_qcfail
                or record.is_duplicate
            ):
                continue
            qname = str(record.query_name or "").strip()
            if not qname:
                continue
            cigar = record.cigarstring or ""
            rname = str(handle.get_reference_name(record.reference_id) or "")
            pos = int(record.reference_start) + 1
            if not rname or not cigar:
                continue
            if current_qname is None:
                current_qname = qname
            if qname != current_qname:
                if current_spans:
                    yield current_spans
                current_qname = qname
                current_spans = set()
            current_spans.add((rname, pos, cigar))
    if current_spans:
        yield current_spans


def _alignment_templates(path: Path) -> Iterator[set[tuple[str, int, str]]]:
    if path.suffix.lower() == ".bam":
        return _iter_bam_templates(path)
    return _iter_sam_templates(path)


def count_gene_expression_from_alignments(
    *,
    alignment_manifest_path: Path,
    gtf_path: Path,
    expected_cell_ids: list[str],
    outdir: Path,
) -> dict[str, Any]:
    started = time.perf_counter()
    genes, genes_by_chrom = _load_cached_gene_index(_gtf_cache_key(gtf_path))
    rows = read_alignment_manifest_tsv(alignment_manifest_path, validate_files=True)
    row_by_cell = {row.cell_id: row for row in rows}
    validate_cell_id_consistency(expected_cell_ids, list(row_by_cell.keys()), circ_label="circ", rna_label="alignment")

    counts_by_gene: dict[str, dict[str, int]] = {
        str(gene["gene_id"]): {cell_id: 0 for cell_id in expected_cell_ids}
        for gene in genes
    }
    per_cell_assigned: dict[str, int] = {cell_id: 0 for cell_id in expected_cell_ids}
    per_cell_ambiguous: dict[str, int] = {cell_id: 0 for cell_id in expected_cell_ids}
    per_cell_unassigned: dict[str, int] = {cell_id: 0 for cell_id in expected_cell_ids}
    per_cell_templates_seen: dict[str, int] = {cell_id: 0 for cell_id in expected_cell_ids}
    per_cell_elapsed_seconds: dict[str, float] = {cell_id: 0.0 for cell_id in expected_cell_ids}

    for cell_id in expected_cell_ids:
        cell_started = time.perf_counter()
        row = row_by_cell[cell_id]
        alignment_path = Path(row.alignment_path)
        for spans in _alignment_templates(alignment_path):
            per_cell_templates_seen[cell_id] += 1
            overlapping_genes: set[str] = set()
            for span in spans:
                chrom, pos, cigar = span
                blocks = _cigar_reference_blocks(int(pos), cigar)
                overlapping_genes.update(
                    _overlapping_gene_ids(chrom=chrom, blocks=blocks, genes_by_chrom=genes_by_chrom)
                )
            if len(overlapping_genes) == 1:
                gene_id = next(iter(overlapping_genes))
                counts_by_gene[gene_id][cell_id] += 1
                per_cell_assigned[cell_id] += 1
            elif len(overlapping_genes) > 1:
                per_cell_ambiguous[cell_id] += 1
            else:
                per_cell_unassigned[cell_id] += 1
        per_cell_elapsed_seconds[cell_id] = round(time.perf_counter() - cell_started, 3)

    feature_df = pd.DataFrame(genes, columns=["gene_id", "gene_name", "chrom", "start", "end", "strand", "gene_biotype"])
    count_rows: list[dict[str, Any]] = []
    for gene in genes:
        row = {
            "gene_id": gene["gene_id"],
            "gene_name": gene["gene_name"],
        }
        row.update({cell_id: counts_by_gene[str(gene["gene_id"])][cell_id] for cell_id in expected_cell_ids})
        count_rows.append(row)
    counts_df = pd.DataFrame(count_rows, columns=["gene_id", "gene_name", *expected_cell_ids])

    outdir.mkdir(parents=True, exist_ok=True)
    gene_counts_out = outdir / "gene_counts.tsv"
    feature_out = outdir / "gene_feature_table.tsv"
    summary_out = outdir / "rna_import_summary.json"
    counts_df.to_csv(gene_counts_out, sep="\t", index=False)
    feature_df.to_csv(feature_out, sep="\t", index=False)

    payload = {
        "method": "simple-overlap",
        "counting_rule": "Count each read or read-pair template once if it overlaps exactly one GTF gene feature; exclude multi-gene ambiguous templates.",
        "paired_end_rule": "Templates are grouped by QNAME per cell to reduce mate double-counting when primary alignments share the same read name.",
        "path": str(alignment_manifest_path.resolve()),
        "gtf": str(gtf_path.resolve()),
        "n_genes": int(feature_df.shape[0]),
        "n_cells": int(len(expected_cell_ids)),
        "cell_ids": list(expected_cell_ids),
        "feature_id_column": "gene_id",
        "feature_name_column": "gene_name",
        "feature_table_columns": ["gene_id", "gene_name", "chrom", "start", "end", "strand", "gene_biotype"],
        "count_table_columns": ["gene_id", "gene_name", *expected_cell_ids],
        "output_gene_counts": str(gene_counts_out.resolve()),
        "output_gene_feature_table": str(feature_out.resolve()),
        "output_rna_import_summary": str(summary_out.resolve()),
        "total_counts_sum": int(counts_df[expected_cell_ids].to_numpy().sum()),
        "cells_processed": int(len(expected_cell_ids)),
        "templates_seen": int(sum(per_cell_templates_seen.values())),
        "assigned_templates": int(sum(per_cell_assigned.values())),
        "ambiguous_templates_excluded": int(sum(per_cell_ambiguous.values())),
        "unassigned_templates": int(sum(per_cell_unassigned.values())),
        "per_cell_templates_seen": per_cell_templates_seen,
        "per_cell_assigned_templates": per_cell_assigned,
        "per_cell_ambiguous_templates_excluded": per_cell_ambiguous,
        "per_cell_unassigned_templates": per_cell_unassigned,
        "per_cell_elapsed_seconds": per_cell_elapsed_seconds,
        "elapsed_seconds": round(time.perf_counter() - started, 3),
    }
    payload.update(
        _write_rna_qc_outputs(
            counts_df=counts_df,
            feature_df=feature_df,
            work_root=outdir.parent,
        )
    )
    summary_out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return payload


def _circ_counts_by_cell(work_root: Path) -> dict[str, int]:
    matrix_dir = work_root / "matrix"
    matrix_path = matrix_dir / "circ_counts.mtx"
    circ_index_path = matrix_dir / "circ_index.txt"
    cell_index_path = matrix_dir / "cell_index.txt"
    if not matrix_path.exists() or not cell_index_path.exists():
        return {}
    X, _, cell_ids = load_circ_matrix(
        matrix_path=matrix_path,
        circ_index_path=circ_index_path if circ_index_path.exists() else cell_index_path,
        cell_index_path=cell_index_path,
    )
    if not cell_ids:
        return {}
    if X.shape[1] == len(cell_ids):
        detected = np.asarray((X > 0).sum(axis=0)).ravel()
    elif X.shape[0] == len(cell_ids):
        detected = np.asarray((X > 0).sum(axis=1)).ravel()
    else:
        return {}
    return {str(cell_id): int(detected[idx]) for idx, cell_id in enumerate(cell_ids)}


def _matrix_cell_ids(work_root: Path) -> list[str]:
    path = work_root / "matrix" / "cell_index.txt"
    return read_index_lines(path)


def _mitochondrial_gene_mask(feature_df: pd.DataFrame) -> pd.Series:
    names = feature_df.get("gene_name", pd.Series([""] * len(feature_df))).astype(str).str.upper()
    ids = feature_df.get("gene_id", pd.Series([""] * len(feature_df))).astype(str).str.upper()
    return names.str.startswith("MT-") | ids.str.startswith("MT-")


def _ribosomal_gene_mask(feature_df: pd.DataFrame) -> pd.Series:
    names = feature_df.get("gene_name", pd.Series([""] * len(feature_df))).astype(str).str.upper()
    prefixes = ("RPL", "RPS", "MRPL", "MRPS")
    return names.str.startswith(prefixes)


def _write_rna_qc_outputs(
    *,
    counts_df: pd.DataFrame,
    feature_df: pd.DataFrame,
    work_root: Path,
) -> dict[str, Any]:
    qc_dir = work_root / "qc"
    qc_dir.mkdir(parents=True, exist_ok=True)
    cell_columns = [column for column in counts_df.columns if column not in {"gene_id", "gene_name"}]
    count_matrix = counts_df[cell_columns].apply(pd.to_numeric, errors="raise")
    feature_qc = feature_df.copy()
    if "gene_biotype" not in feature_qc.columns:
        feature_qc["gene_biotype"] = ""
    mt_mask = _mitochondrial_gene_mask(feature_qc)
    ribo_mask = _ribosomal_gene_mask(feature_qc)

    total_counts = count_matrix.sum(axis=0)
    detected_genes = (count_matrix > 0).sum(axis=0)
    mt_counts = count_matrix.loc[mt_mask, cell_columns].sum(axis=0) if mt_mask.any() else pd.Series(0, index=cell_columns)
    ribo_counts = count_matrix.loc[ribo_mask, cell_columns].sum(axis=0) if ribo_mask.any() else pd.Series(0, index=cell_columns)
    circ_counts = _circ_counts_by_cell(work_root)
    matrix_cell_ids = _matrix_cell_ids(work_root)
    shared_cells = sorted(set(cell_columns) & set(matrix_cell_ids))
    rna_only_cells = sorted(set(cell_columns) - set(matrix_cell_ids))
    circ_only_cells = sorted(set(matrix_cell_ids) - set(cell_columns))

    per_cell = pd.DataFrame(
        {
            "cell_id": cell_columns,
            "total_counts": [int(total_counts[cell_id]) for cell_id in cell_columns],
            "detected_genes": [int(detected_genes[cell_id]) for cell_id in cell_columns],
            "mitochondrial_fraction": [
                float(mt_counts[cell_id] / total_counts[cell_id]) if float(total_counts[cell_id]) > 0 else 0.0
                for cell_id in cell_columns
            ],
            "ribosomal_fraction": [
                float(ribo_counts[cell_id] / total_counts[cell_id]) if float(total_counts[cell_id]) > 0 else 0.0
                for cell_id in cell_columns
            ],
            "circRNA_count": [int(circ_counts.get(str(cell_id), 0)) for cell_id in cell_columns],
        }
    )

    gene_totals = count_matrix.sum(axis=1)
    gene_detected = (count_matrix > 0).sum(axis=1)
    nonzero_means: list[float] = []
    for row_idx in range(count_matrix.shape[0]):
        nonzero = count_matrix.iloc[row_idx]
        nonzero = nonzero[nonzero > 0]
        nonzero_means.append(float(nonzero.mean()) if not nonzero.empty else 0.0)
    per_gene = feature_qc[["gene_id", "gene_name", "gene_biotype"]].copy()
    per_gene["n_cells_detected"] = gene_detected.astype(int)
    per_gene["total_counts"] = gene_totals.astype(int)
    per_gene["mean_counts_nonzero"] = nonzero_means

    cell_qc_path = qc_dir / "rna_qc.tsv"
    gene_qc_path = qc_dir / "rna_gene_qc.tsv"
    per_cell.to_csv(cell_qc_path, sep="\t", index=False)
    per_gene.to_csv(gene_qc_path, sep="\t", index=False)

    genes_detected = gene_detected.astype(int)
    return {
        "output_rna_qc": str(cell_qc_path.resolve()),
        "output_rna_gene_qc": str(gene_qc_path.resolve()),
        "genes_detected_ge_1_cell": int((genes_detected >= 1).sum()),
        "genes_detected_ge_3_cells": int((genes_detected >= 3).sum()),
        "genes_detected_ge_10_cells": int((genes_detected >= 10).sum()),
        "median_genes_per_cell": float(detected_genes.median()) if len(detected_genes) else 0.0,
        "matrix_cell_id_match": bool(matrix_cell_ids) and not rna_only_cells and not circ_only_cells,
        "matrix_shared_cell_count": int(len(shared_cells)),
        "matrix_rna_only_cell_count": int(len(rna_only_cells)),
        "matrix_circ_only_cell_count": int(len(circ_only_cells)),
        "filtering_report": {
            "automatic_filtering_applied": False,
            "threshold_summary_only": True,
        },
    }


def refresh_rna_qc_from_existing_outputs(
    *,
    workdir: Path,
) -> dict[str, Any]:
    rna_dir = workdir / "rna"
    gene_counts_path = rna_dir / "gene_counts.tsv"
    feature_path = rna_dir / "gene_feature_table.tsv"
    summary_path = rna_dir / "rna_import_summary.json"
    if not gene_counts_path.exists():
        raise FileNotFoundError(f"Missing RNA gene-count table: {gene_counts_path}")
    if not feature_path.exists():
        raise FileNotFoundError(f"Missing RNA gene feature table: {feature_path}")

    counts_df = pd.read_csv(gene_counts_path, sep="\t", keep_default_na=False)
    feature_df = pd.read_csv(feature_path, sep="\t", keep_default_na=False)
    required_count_columns = {"gene_id", "gene_name"}
    missing_counts = sorted(required_count_columns - set(counts_df.columns))
    if missing_counts:
        raise ValueError(
            f"{gene_counts_path} is missing required columns: {', '.join(missing_counts)}"
        )
    if "gene_id" not in feature_df.columns or "gene_name" not in feature_df.columns:
        raise ValueError(
            f"{feature_path} must contain at least gene_id and gene_name columns."
        )

    count_cell_ids = [str(column) for column in counts_df.columns if column not in {"gene_id", "gene_name"}]
    if not count_cell_ids:
        raise ValueError(f"{gene_counts_path} contains no RNA cell columns.")
    validate_feature_id_uniqueness(counts_df["gene_id"].astype(str).tolist(), label="gene")
    validate_feature_id_uniqueness(count_cell_ids, label="cell")
    if summary_path.exists():
        existing_summary = load_json(summary_path)
    else:
        existing_summary = {}

    qc_summary = _write_rna_qc_outputs(
        counts_df=counts_df,
        feature_df=feature_df,
        work_root=workdir,
    )

    refreshed_summary = {
        **existing_summary,
        "method": existing_summary.get("method", "existing-rna-refresh"),
        "n_cells": int(len(count_cell_ids)),
        "n_genes": int(counts_df.shape[0]),
        "cell_ids": count_cell_ids,
        "output_gene_counts": str(gene_counts_path.resolve()),
        "output_gene_feature_table": str(feature_path.resolve()),
        "output_rna_import_summary": str(summary_path.resolve()),
        **qc_summary,
    }
    summary_path.write_text(json.dumps(refreshed_summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    workflow_summary_path = workdir / "workflow_summary.json"
    workflow_summary_updated = False
    if workflow_summary_path.exists():
        workflow_summary = load_json(workflow_summary_path)
        workflow_summary["rna_import"] = {
            **dict(workflow_summary.get("rna_import", {})),
            **refreshed_summary,
        }
        workflow_summary = apply_standard_provenance(
            workflow_summary,
            command_name="circyto refresh-rna-qc",
            workflow_type=str(workflow_summary.get("workflow_type", workflow_summary.get("workflow", "posthoc-rna-profile"))),
            protocol=str(workflow_summary.get("protocol", "unknown")),
            read_layout=str(workflow_summary.get("read_layout", "unknown")),
            genome_fasta=str(workflow_summary.get("genome_fasta")) if workflow_summary.get("genome_fasta") else None,
            gtf=str(workflow_summary.get("gtf")) if workflow_summary.get("gtf") else None,
            detector_backend=str(workflow_summary.get("detector_backend")) if workflow_summary.get("detector_backend") else None,
            started_at=str(workflow_summary.get("started_at", utc_now_iso())),
            completed_at=utc_now_iso(),
            elapsed_seconds=float(workflow_summary.get("elapsed_seconds", 0.0) or 0.0),
            cleanup=workflow_summary.get("cleanup"),
            cleanup_scope=workflow_summary.get("cleanup_scope"),
            cleanup_performed=bool(workflow_summary.get("cleanup_performed", False)),
            cleanup_deleted_paths=list(workflow_summary.get("cleanup_deleted_paths", [])),
            cleanup_reclaimed_bytes=int(workflow_summary.get("cleanup_reclaimed_bytes", 0) or 0),
            workflow_uuid=str(workflow_summary.get("workflow_uuid")) if workflow_summary.get("workflow_uuid") else None,
        )
        write_json(workflow_summary_path, workflow_summary)
        workflow_summary_updated = True

    return apply_standard_provenance(
        {
            "workdir": str(workdir.resolve()),
            "rna_import": refreshed_summary,
            "workflow_summary_updated": workflow_summary_updated,
        },
        command_name="circyto refresh-rna-qc",
        workflow_type="refresh-rna-qc",
        protocol="unknown",
        read_layout="unknown",
        genome_fasta=None,
        gtf=None,
        detector_backend=None,
        started_at=utc_now_iso(),
        completed_at=utc_now_iso(),
        elapsed_seconds=0.0,
    )


def find_completed_workflow_alignment_manifest(workdir: Path) -> Path:
    candidates = [
        workdir / "align" / "alignment_manifest.tsv",
        workdir / "align" / "star" / "alignment_manifest.tsv",
        workdir / "align" / "bwa_mem" / "alignment_manifest.tsv",
    ]
    for path in candidates:
        if path.exists():
            return path
    tried = ", ".join(str(path) for path in candidates)
    raise FileNotFoundError(
        f"Could not find alignment manifest under completed workflow directory {workdir}. "
        f"Tried: {tried}"
    )


def build_posthoc_rna_profile_plan(
    *,
    workdir: Path,
    gtf_path: Path,
    method: str,
) -> dict[str, Any]:
    if method != "simple-overlap":
        raise ValueError(f"Unsupported post-hoc RNA profiling method: {method}")
    alignment_manifest_path = find_completed_workflow_alignment_manifest(workdir)
    return {
        "command_name": "circyto add-rna-profile",
        "dry_run": True,
        "workdir": str(workdir.resolve()),
        "method": method,
        "alignment_manifest": str(alignment_manifest_path.resolve()),
        "gtf": str(gtf_path.resolve()),
        "output_gene_counts": str((workdir / "rna" / "gene_counts.tsv").resolve()),
        "output_gene_feature_table": str((workdir / "rna" / "gene_feature_table.tsv").resolve()),
        "output_rna_import_summary": str((workdir / "rna" / "rna_import_summary.json").resolve()),
    }


def add_posthoc_rna_profile(
    *,
    workdir: Path,
    gtf_path: Path,
    method: str,
    dry_run: bool,
) -> dict[str, Any]:
    plan = build_posthoc_rna_profile_plan(workdir=workdir, gtf_path=gtf_path, method=method)
    if dry_run:
        return apply_standard_provenance(
            dict(plan),
            command_name="circyto add-rna-profile",
            workflow_type="posthoc-rna-profile",
            protocol="unknown",
            read_layout="unknown",
            genome_fasta=None,
            gtf=str(gtf_path.resolve()),
            detector_backend=None,
            started_at=utc_now_iso(),
            completed_at=utc_now_iso(),
            elapsed_seconds=0.0,
        )

    alignment_manifest_path = Path(plan["alignment_manifest"])
    rows = read_alignment_manifest_tsv(alignment_manifest_path, validate_files=True)
    started_at = utc_now_iso()
    rna_import = count_gene_expression_from_alignments(
        alignment_manifest_path=alignment_manifest_path,
        gtf_path=gtf_path,
        expected_cell_ids=[row.cell_id for row in rows],
        outdir=workdir / "rna",
    )

    workflow_summary_path = workdir / "workflow_summary.json"
    workflow_uuid = None
    protocol = "unknown"
    read_layout = summarize_read_layouts([row.read_layout for row in rows])
    genome_fasta = None
    detector_backend = None
    if workflow_summary_path.exists():
        try:
            workflow_summary = load_json(workflow_summary_path)
        except json.JSONDecodeError as exc:
            raise ValueError(f"Could not parse existing workflow summary: {workflow_summary_path}") from exc
        workflow_uuid = workflow_summary.get("workflow_uuid")
        protocol = str(workflow_summary.get("protocol", protocol))
        read_layout = str(workflow_summary.get("read_layout", read_layout))
        genome_fasta = workflow_summary.get("genome_fasta")
        detector_backend = workflow_summary.get("detector_backend")
        workflow_summary["rna_import"] = rna_import
        workflow_summary = apply_standard_provenance(
            workflow_summary,
            command_name="circyto add-rna-profile",
            workflow_type=str(workflow_summary.get("workflow_type", workflow_summary.get("workflow", "posthoc-rna-profile"))),
            protocol=protocol,
            read_layout=read_layout,
            genome_fasta=str(genome_fasta) if genome_fasta else None,
            gtf=str(gtf_path.resolve()),
            detector_backend=str(detector_backend) if detector_backend else None,
            started_at=workflow_summary.get("started_at", started_at),
            completed_at=utc_now_iso(),
            elapsed_seconds=float(workflow_summary.get("elapsed_seconds", 0.0) or 0.0),
            cleanup=workflow_summary.get("cleanup"),
            cleanup_scope=workflow_summary.get("cleanup_scope"),
            cleanup_performed=bool(workflow_summary.get("cleanup_performed", False)),
            cleanup_deleted_paths=list(workflow_summary.get("cleanup_deleted_paths", [])),
            cleanup_reclaimed_bytes=int(workflow_summary.get("cleanup_reclaimed_bytes", 0) or 0),
            workflow_uuid=str(workflow_uuid) if workflow_uuid else None,
        )
        write_json(workflow_summary_path, workflow_summary)

    result = {
        "dry_run": False,
        "workdir": str(workdir.resolve()),
        "method": method,
        "alignment_manifest": str(alignment_manifest_path.resolve()),
        "rna_import": rna_import,
        "workflow_summary_updated": workflow_summary_path.exists(),
    }
    return apply_standard_provenance(
        result,
        command_name="circyto add-rna-profile",
        workflow_type="posthoc-rna-profile",
        protocol=protocol,
        read_layout=read_layout,
        genome_fasta=str(genome_fasta) if genome_fasta else None,
        gtf=str(gtf_path.resolve()),
        detector_backend=str(detector_backend) if detector_backend else None,
        started_at=started_at,
        completed_at=utc_now_iso(),
        elapsed_seconds=0.0,
        workflow_uuid=str(workflow_uuid) if workflow_uuid else None,
    )


def validate_velocity_layers_schema(path: Path) -> dict[str, Any]:
    if not path.exists() or not path.is_dir():
        raise FileNotFoundError(f"Velocity layer directory not found: {path}")
    names = {child.name for child in path.iterdir()}
    missing = sorted(VALID_VELOCITY_LAYOUT_FILES - names)
    if missing:
        raise ValueError(
            f"{path} is missing required velocity-layer files: {', '.join(missing)}. "
            "Expected at minimum: barcodes.tsv, features.tsv, spliced.mtx, unspliced.mtx."
        )
    barcodes = [
        line.strip()
        for line in (path / "barcodes.tsv").read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if not barcodes:
        raise ValueError(f"{path / 'barcodes.tsv'} has no cell IDs")
    feature_rows = [
        line.rstrip("\n").split("\t")
        for line in (path / "features.tsv").read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if not feature_rows:
        raise ValueError(f"{path / 'features.tsv'} has no feature rows")
    gene_ids = [row[0].strip() for row in feature_rows if row and row[0].strip()]
    if len(gene_ids) != len(feature_rows):
        raise ValueError(f"{path / 'features.tsv'} contains blank gene_id entries")
    validate_feature_id_uniqueness(gene_ids, label="gene")
    validate_feature_id_uniqueness(barcodes, label="cell")
    return {
        "path": str(path.resolve()),
        "n_cells": int(len(barcodes)),
        "n_genes": int(len(gene_ids)),
        "has_ambiguous": (path / "ambiguous.mtx").exists(),
        "cell_ids": barcodes,
        "feature_ids": gene_ids,
    }


def _collect_cleanup_candidates(outdir: Path, *, include_demux_fastq: bool) -> list[dict[str, Any]]:
    candidates: list[dict[str, Any]] = []
    for subdir_name in ("align", "ciri3", "demux"):
        subdir = outdir / subdir_name
        if not subdir.exists():
            continue
        for root, _, filenames in os.walk(subdir):
            for filename in filenames:
                path = Path(root) / filename
                suffixes = "".join(path.suffixes[-2:]).lower()
                include_alignment = (
                    path.suffix.lower() in REGENERABLE_ALIGNMENT_SUFFIXES
                    or filename in REGENERABLE_ALIGNMENT_NAMES
                    or path.suffix.lower() in REGENERABLE_CHUNK_SUFFIXES
                )
                include_demux = include_demux_fastq and (
                    suffixes in REGENERABLE_DEMUX_SUFFIXES
                    or path.suffix.lower() in REGENERABLE_CHUNK_SUFFIXES
                )
                if (
                    include_alignment
                    or (
                        subdir_name == "demux"
                        and include_demux
                    )
                ):
                    try:
                        size = int(path.stat().st_size)
                    except OSError:
                        size = 0
                    candidates.append({"path": str(path.resolve()), "bytes": size})
    return _deduplicate_cleanup_candidates(candidates)


def _deduplicate_cleanup_candidates(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    ordered = sorted(candidates, key=lambda item: (-int(item["bytes"]), str(item["path"])))
    seen: set[str] = set()
    unique: list[dict[str, Any]] = []
    for item in ordered:
        path = str(item.get("path", ""))
        if not path or path in seen:
            continue
        seen.add(path)
        unique.append({"path": path, "bytes": int(item.get("bytes", 0) or 0)})
    return unique


def _resolve_cleanup_candidates(
    *,
    outdir: Path,
    scope: str,
    workflow_succeeded: bool,
) -> list[dict[str, Any]]:
    if not workflow_succeeded:
        return []
    candidates = _collect_cleanup_candidates(outdir, include_demux_fastq=scope in {"demux", "all"})
    if scope == "alignments":
        return _deduplicate_cleanup_candidates(
            [item for item in candidates if "/demux/" not in str(item["path"]).replace("\\", "/")]
        )
    if scope == "demux":
        return _deduplicate_cleanup_candidates(
            [item for item in candidates if "/demux/" in str(item["path"]).replace("\\", "/")]
        )
    return _deduplicate_cleanup_candidates(candidates)


def build_cleanup_plan(
    *,
    outdir: Path,
    scope: str = "all",
    workflow_succeeded: bool = True,
) -> dict[str, Any]:
    if scope not in SUPPORTED_CLEANUP_SCOPES:
        raise ValueError(
            f"Unsupported cleanup scope: {scope}. "
            f"Choose from: {', '.join(SUPPORTED_CLEANUP_SCOPES)}"
        )

    safe_to_delete_after_success: list[str] = []
    include_demux_fastq = scope in {"demux", "all"}
    if scope in {"alignments", "all"}:
        safe_to_delete_after_success.extend(SAFE_TO_DELETE_ALIGNMENTS)
    if scope in {"demux", "all"}:
        safe_to_delete_after_success.extend(SAFE_TO_DELETE_DEMUX)

    if not workflow_succeeded:
        return {
            "enabled": True,
            "mode": "planned-only",
            "workflow_succeeded": False,
            "planned_scope": scope,
            "supported_scopes": list(SUPPORTED_CLEANUP_SCOPES),
            "must_keep": list(MUST_KEEP_OUTPUTS),
            "user_inputs_never_delete": list(USER_INPUTS_NEVER_DELETE),
            "safe_to_delete_after_success": safe_to_delete_after_success,
            "candidate_count": 0,
            "candidate_bytes": 0,
            "delete_candidates": [],
            "reason": "Cleanup planning is disabled for failed workflows.",
        }

    candidates = _resolve_cleanup_candidates(outdir=outdir, scope=scope, workflow_succeeded=workflow_succeeded)
    return {
        "enabled": True,
        "mode": "planned-only",
        "workflow_succeeded": True,
        "planned_scope": scope,
        "supported_scopes": list(SUPPORTED_CLEANUP_SCOPES),
        "would_delete_only_regenerable_workflow_intermediates": True,
        "must_keep": list(MUST_KEEP_OUTPUTS),
        "user_inputs_never_delete": list(USER_INPUTS_NEVER_DELETE),
        "safe_to_delete_after_success": safe_to_delete_after_success,
        "candidate_count": len(candidates),
        "candidate_bytes": int(sum(int(item["bytes"]) for item in candidates)),
        "delete_candidates": candidates[:20],
    }


def execute_cleanup_plan(
    *,
    outdir: Path,
    scope: str = "all",
    workflow_succeeded: bool = True,
) -> dict[str, Any]:
    plan = build_cleanup_plan(outdir=outdir, scope=scope, workflow_succeeded=workflow_succeeded)
    if not workflow_succeeded:
        return {
            **plan,
            "cleanup_performed": False,
            "cleanup_deleted_paths": [],
            "cleanup_reclaimed_bytes": 0,
            "cleanup_scope": scope,
        }

    deleted_paths: list[str] = []
    reclaimed_bytes = 0
    for item in _resolve_cleanup_candidates(outdir=outdir, scope=scope, workflow_succeeded=workflow_succeeded):
        path = Path(str(item["path"]))
        if not path.exists():
            continue
        try:
            size = int(path.stat().st_size)
        except OSError:
            size = int(item.get("bytes", 0) or 0)
        try:
            path.unlink()
        except OSError:
            continue
        deleted_paths.append(str(path.resolve()))
        reclaimed_bytes += size

    return {
        **plan,
        "mode": "executed",
        "cleanup_performed": True,
        "cleanup_deleted_paths": deleted_paths,
        "cleanup_reclaimed_bytes": int(reclaimed_bytes),
        "cleanup_scope": scope,
    }


def _workflow_is_marked_successful(summary: dict[str, Any] | None) -> bool:
    if not summary:
        return False
    failed_stages = [str(item).strip() for item in summary.get("failed_stages", []) if str(item).strip()]
    partial_outputs = [str(item).strip() for item in summary.get("partial_outputs_detected", []) if str(item).strip()]
    completed_at = str(summary.get("completed_at", "")).strip()
    return bool(completed_at) and not failed_stages and not partial_outputs


def cleanup_completed_workflow(
    *,
    workdir: Path,
    scope: str,
    dry_run: bool,
    force: bool,
) -> dict[str, Any]:
    integrity = check_workflow_integrity(workdir)
    workflow_summary_path = workdir / "workflow_summary.json"
    workflow_summary: dict[str, Any] | None = None
    if workflow_summary_path.exists():
        try:
            workflow_summary = load_json(workflow_summary_path)
        except json.JSONDecodeError as exc:
            if not force:
                raise ValueError(f"Could not parse workflow summary: {workflow_summary_path}") from exc

    if not integrity.get("ok", False) and not force:
        raise ValueError(
            "Workflow integrity check failed; refusing cleanup without --force. "
            "Run `circyto check-workflow --workdir ...` to inspect the problems first."
        )

    eligible_for_cleanup = _workflow_is_marked_successful(workflow_summary) or force
    plan = build_cleanup_plan(outdir=workdir, scope=scope, workflow_succeeded=eligible_for_cleanup)

    if dry_run:
        return {
            "command_name": "circyto cleanup-workflow",
            "dry_run": True,
            "force": force,
            "workdir": str(workdir.resolve()),
            "scope": scope,
            "workflow_integrity_ok": bool(integrity.get("ok", False)),
            "workflow_succeeded": bool(eligible_for_cleanup),
            "cleanup": plan,
            "estimated_reclaimed_bytes": int(plan.get("candidate_bytes", 0) or 0),
            "planned_deletions": list(plan.get("delete_candidates", [])),
            "protected_paths": {
                "must_keep": list(MUST_KEEP_OUTPUTS),
                "user_inputs_never_delete": list(USER_INPUTS_NEVER_DELETE),
            },
            "integrity": integrity,
        }

    result = execute_cleanup_plan(outdir=workdir, scope=scope, workflow_succeeded=eligible_for_cleanup)
    if workflow_summary is not None:
        workflow_summary["cleanup"] = result
        workflow_summary["cleanup_performed"] = bool(result.get("cleanup_performed", False))
        workflow_summary["cleanup_deleted_paths"] = list(result.get("cleanup_deleted_paths", []))
        workflow_summary["cleanup_reclaimed_bytes"] = int(result.get("cleanup_reclaimed_bytes", 0) or 0)
        workflow_summary["cleanup_scope"] = result.get("cleanup_scope")
        workflow_summary["cleanup_summary"] = cleanup_summary_block(
            cleanup=result,
            cleanup_scope=scope,
            cleanup_performed=bool(result.get("cleanup_performed", False)),
            cleanup_deleted_paths=list(result.get("cleanup_deleted_paths", [])),
            cleanup_reclaimed_bytes=int(result.get("cleanup_reclaimed_bytes", 0) or 0),
        )
        write_json(workflow_summary_path, workflow_summary)

    return {
        "command_name": "circyto cleanup-workflow",
        "dry_run": False,
        "force": force,
        "workdir": str(workdir.resolve()),
        "scope": scope,
        "workflow_integrity_ok": bool(integrity.get("ok", False)),
        "workflow_succeeded": bool(eligible_for_cleanup),
        "cleanup": result,
        "protected_paths": {
            "must_keep": list(MUST_KEEP_OUTPUTS),
            "user_inputs_never_delete": list(USER_INPUTS_NEVER_DELETE),
        },
        "workflow_summary_updated": workflow_summary is not None,
        "integrity": integrity,
    }
