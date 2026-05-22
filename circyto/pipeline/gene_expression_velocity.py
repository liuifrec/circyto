from __future__ import annotations

import csv
import os
import json
import platform
import re
import time
from bisect import bisect_left, bisect_right
from collections import defaultdict
from importlib import metadata
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterator

import numpy as np
import pandas as pd
import scipy.sparse as sp

from circyto import __version__
from circyto.manifest.alignment import read_alignment_manifest_tsv
from circyto.pipeline.workflow_reporting import (
    apply_standard_provenance,
    cleanup_summary_block,
    directory_size_bytes,
    load_circ_feature_table,
    load_circ_matrix,
    load_json,
    numeric_summary,
    read_index_lines,
    sanitize_frame_for_anndata,
    summarize_read_layouts,
    utc_now_iso,
    write_json,
)
from circyto.pipeline.workflow_integrity import check_workflow_integrity

try:
    import anndata as ad  # type: ignore

    HAS_ANNDATA = True
except Exception:
    HAS_ANNDATA = False

try:
    import mudata as mu  # type: ignore

    HAS_MUDATA = True
except Exception:
    HAS_MUDATA = False


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
    require_ordered_equality: bool = False,
) -> list[str]:
    circ = [str(cell_id).strip() for cell_id in circ_cell_ids if str(cell_id).strip()]
    rna = [str(cell_id).strip() for cell_id in rna_cell_ids if str(cell_id).strip()]
    if not circ:
        raise ValueError(f"{circ_label} cell ID set is empty")
    if not rna:
        raise ValueError(f"{rna_label} cell ID set is empty")
    if set(circ) != set(rna):
        only_circ = sorted(set(circ) - set(rna))
        only_rna = sorted(set(rna) - set(circ))
        raise ValueError(
            f"Cell ID mismatch between {circ_label} and {rna_label}. "
            f"Only in {circ_label}: {only_circ[:5]}. Only in {rna_label}: {only_rna[:5]}."
        )
    if require_ordered_equality and circ != rna:
        raise ValueError(
            f"Cell ID order mismatch between {circ_label} and {rna_label}. "
            f"Expected {circ_label} order starts with {circ[:5]}. "
            f"Observed {rna_label} order starts with {rna[:5]}."
        )
    return list(circ)


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
    ordered_cell_ids = validate_cell_id_consistency(
        expected_cell_ids,
        [str(column) for column in cell_columns],
        circ_label="circ",
        rna_label="rna",
    )

    numeric = df[cell_columns].apply(pd.to_numeric, errors="raise")
    numeric = numeric[ordered_cell_ids]
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
        "cell_ids": list(ordered_cell_ids),
        "output_gene_counts": str(gene_counts_out.resolve()),
        "output_gene_feature_table": str(feature_out.resolve()),
        "feature_table_columns": ["gene_id", "gene_name"],
        "count_table_columns": ["gene_id", "gene_name", *ordered_cell_ids],
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
    ordered_cell_ids = validate_cell_id_consistency(
        expected_cell_ids,
        list(row_by_cell.keys()),
        circ_label="circ",
        rna_label="alignment",
    )

    counts_by_gene: dict[str, dict[str, int]] = {
        str(gene["gene_id"]): {cell_id: 0 for cell_id in ordered_cell_ids}
        for gene in genes
    }
    per_cell_assigned: dict[str, int] = {cell_id: 0 for cell_id in ordered_cell_ids}
    per_cell_ambiguous: dict[str, int] = {cell_id: 0 for cell_id in ordered_cell_ids}
    per_cell_unassigned: dict[str, int] = {cell_id: 0 for cell_id in ordered_cell_ids}
    per_cell_templates_seen: dict[str, int] = {cell_id: 0 for cell_id in ordered_cell_ids}
    per_cell_elapsed_seconds: dict[str, float] = {cell_id: 0.0 for cell_id in ordered_cell_ids}

    for cell_id in ordered_cell_ids:
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
        row.update({cell_id: counts_by_gene[str(gene["gene_id"])][cell_id] for cell_id in ordered_cell_ids})
        count_rows.append(row)
    counts_df = pd.DataFrame(count_rows, columns=["gene_id", "gene_name", *ordered_cell_ids])

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
        "n_cells": int(len(ordered_cell_ids)),
        "cell_ids": list(ordered_cell_ids),
        "feature_id_column": "gene_id",
        "feature_name_column": "gene_name",
        "feature_table_columns": ["gene_id", "gene_name", "chrom", "start", "end", "strand", "gene_biotype"],
        "count_table_columns": ["gene_id", "gene_name", *ordered_cell_ids],
        "output_gene_counts": str(gene_counts_out.resolve()),
        "output_gene_feature_table": str(feature_out.resolve()),
        "output_rna_import_summary": str(summary_out.resolve()),
        "total_counts_sum": int(counts_df[ordered_cell_ids].to_numpy().sum()),
        "cells_processed": int(len(ordered_cell_ids)),
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
    return _circ_metrics_by_cell(work_root)["circRNA_count"]


def _circ_metrics_by_cell(work_root: Path) -> dict[str, dict[str, int]] | dict[str, Any]:
    matrix_dir = work_root / "matrix"
    matrix_path = matrix_dir / "circ_counts.mtx"
    circ_index_path = matrix_dir / "circ_index.txt"
    cell_index_path = matrix_dir / "cell_index.txt"
    if not matrix_path.exists() or not cell_index_path.exists():
        return {"circRNA_count": {}, "circRNA_total_support": {}}
    X, _, cell_ids = load_circ_matrix(
        matrix_path=matrix_path,
        circ_index_path=circ_index_path if circ_index_path.exists() else cell_index_path,
        cell_index_path=cell_index_path,
    )
    if not cell_ids:
        return {"circRNA_count": {}, "circRNA_total_support": {}}
    try:
        X_cells_by_circ = _orient_circ_matrix_cells_by_circ(
            X=X,
            circ_ids=read_index_lines(circ_index_path) if circ_index_path.exists() else [],
            cell_ids=cell_ids,
            matrix_path=matrix_path,
        )
    except ValueError:
        return {"circRNA_count": {}, "circRNA_total_support": {}}
    detected = np.asarray((X_cells_by_circ > 0).sum(axis=1)).ravel()
    total_support = np.asarray(X_cells_by_circ.sum(axis=1)).ravel()
    return {
        "circRNA_count": {str(cell_id): int(detected[idx]) for idx, cell_id in enumerate(cell_ids)},
        "circRNA_total_support": {str(cell_id): int(total_support[idx]) for idx, cell_id in enumerate(cell_ids)},
    }


def _matrix_cell_ids(work_root: Path) -> list[str]:
    path = work_root / "matrix" / "cell_index.txt"
    return read_index_lines(path)


def sanitize_for_h5ad_uns(obj: Any) -> Any:
    if obj is None:
        return ""
    if isinstance(obj, (str, int, float, bool, np.integer, np.floating, np.bool_)):
        return obj
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, dict):
        safe: dict[str, Any] = {}
        for key, value in obj.items():
            key_str = str(key)
            if isinstance(value, dict):
                safe[key_str] = json.dumps(value, sort_keys=True, default=str)
            elif isinstance(value, list):
                if all(isinstance(item, (str, int, float, bool, np.integer, np.floating, np.bool_)) or item is None for item in value):
                    safe[key_str] = ["" if item is None else item for item in value]
                else:
                    safe[key_str] = json.dumps(value, sort_keys=True, default=str)
            else:
                safe[key_str] = sanitize_for_h5ad_uns(value)
        return safe
    if isinstance(obj, list):
        if all(isinstance(item, (str, int, float, bool, np.integer, np.floating, np.bool_)) or item is None for item in obj):
            return ["" if item is None else item for item in obj]
        return json.dumps(obj, sort_keys=True, default=str)
    return json.dumps(obj, sort_keys=True, default=str)


def _require_mudata_for_readers() -> None:
    if not HAS_ANNDATA:
        raise RuntimeError("anndata is required for MuData inspection")
    if not HAS_MUDATA:
        raise RuntimeError("mudata is not installed; install circyto[mudata] or pip install mudata")


def inspect_mudata_file(input_path: Path) -> dict[str, Any]:
    _require_mudata_for_readers()
    if not input_path.exists():
        raise FileNotFoundError(f"MuData input not found: {input_path}")

    mdata = mu.read_h5mu(str(input_path))
    modalities = list(mdata.mod.keys())
    rna = mdata.mod["rna"] if "rna" in mdata.mod else None
    circ = mdata.mod["circ"] if "circ" in mdata.mod else None
    membership_counts: dict[str, int] = {}
    if "membership" in mdata.obs.columns:
        membership_counts = {
            str(key): int(value)
            for key, value in mdata.obs["membership"].astype(str).value_counts().to_dict().items()
        }

    circyto_uns = mdata.uns.get("circyto", {})
    if isinstance(circyto_uns, dict):
        uns_keys = sorted(str(key) for key in circyto_uns.keys())
    else:
        uns_keys = []

    return {
        "input": str(input_path.resolve()),
        "modalities": modalities,
        "n_obs": int(mdata.n_obs),
        "rna_shape": list(rna.shape) if rna is not None else None,
        "circ_shape": list(circ.shape) if circ is not None else None,
        "obs_columns": [str(column) for column in mdata.obs.columns.tolist()],
        "rna_var_columns": [str(column) for column in rna.var.columns.tolist()] if rna is not None else [],
        "circ_var_columns": [str(column) for column in circ.var.columns.tolist()] if circ is not None else [],
        "circyto_uns_keys": uns_keys,
        "membership_counts": membership_counts,
        "n_shared_cells": int(membership_counts.get("shared", 0)),
        "n_rna_only_cells": int(membership_counts.get("rna_only", 0)),
        "n_circ_only_cells": int(membership_counts.get("circ_only", 0)),
    }


def summarize_mudata_qc(input_path: Path) -> dict[str, Any]:
    _require_mudata_for_readers()
    if not input_path.exists():
        raise FileNotFoundError(f"MuData input not found: {input_path}")
    mdata = mu.read_h5mu(str(input_path))
    obs = mdata.obs.copy()

    summary: dict[str, Any] = {
        "input": str(input_path.resolve()),
        "n_obs": int(mdata.n_obs),
    }
    for column in (
        "total_rna_counts",
        "detected_genes",
        "mitochondrial_fraction",
        "ribosomal_fraction",
        "circRNA_count",
        "circRNA_total_support",
    ):
        if column in obs.columns:
            values = pd.to_numeric(obs[column], errors="coerce").dropna().tolist()
            summary[column] = numeric_summary(values)
        else:
            summary[column] = None

    correlation = None
    if "total_rna_counts" in obs.columns and "circRNA_count" in obs.columns:
        pair_df = obs[["total_rna_counts", "circRNA_count"]].apply(pd.to_numeric, errors="coerce").dropna()
        if len(pair_df) >= 2 and pair_df["total_rna_counts"].nunique() > 1 and pair_df["circRNA_count"].nunique() > 1:
            correlation = float(pair_df["total_rna_counts"].corr(pair_df["circRNA_count"]))
    summary["pearson_total_rna_vs_circRNA_count"] = correlation
    return summary


def _installed_version(package_name: str) -> str | None:
    try:
        return metadata.version(package_name)
    except metadata.PackageNotFoundError:
        return None


def get_environment_summary() -> dict[str, Any]:
    return {
        "circyto_version": __version__,
        "python_version": platform.python_version(),
        "anndata_version": _installed_version("anndata"),
        "mudata_version": _installed_version("mudata"),
        "scanpy_version": _installed_version("scanpy"),
        "platform": platform.platform(),
        "timestamp": utc_now_iso(),
    }


def inspect_completed_workdir(workdir: Path) -> dict[str, Any]:
    workflow_summary_path = workdir / "workflow_summary.json"
    workflow_summary = load_json(workflow_summary_path) if workflow_summary_path.exists() else {}

    matrix_dir = workdir / "matrix"
    rna_dir = workdir / "rna"
    qc_dir = workdir / "qc"
    mudata_dir = workdir / "mudata"
    dna_cnv_dir = workdir / "dna_cnv"
    dna_snv_dir = workdir / "dna_snv"

    matrices_present = {
        "circ_counts_mtx": (matrix_dir / "circ_counts.mtx").exists(),
        "circ_index": (matrix_dir / "circ_index.txt").exists(),
        "cell_index": (matrix_dir / "cell_index.txt").exists(),
        "circ_feature_table": (matrix_dir / "circ_feature_table.tsv").exists(),
        "gene_counts_tsv": (rna_dir / "gene_counts.tsv").exists(),
        "gene_feature_table": (rna_dir / "gene_feature_table.tsv").exists(),
    }
    modalities: list[str] = []
    if matrices_present["gene_counts_tsv"]:
        modalities.append("rna")
    if matrices_present["circ_counts_mtx"] and matrices_present["cell_index"]:
        modalities.append("circ")
    if dna_cnv_dir.exists():
        modalities.append("dna_cnv")
    if dna_snv_dir.exists():
        modalities.append("dna_snv")

    mudata_paths = sorted(str(path.relative_to(workdir)) for path in mudata_dir.glob("*.h5mu")) if mudata_dir.exists() else []
    qc_files = sorted(str(path.relative_to(workdir)) for path in qc_dir.glob("*")) if qc_dir.exists() else []

    return apply_standard_provenance(
        {
            "workdir": str(workdir.resolve()),
            "available_modalities": modalities,
            "source_workflow_type": str(workflow_summary.get("workflow_type", workflow_summary.get("workflow", "unknown"))),
            "source_protocol": str(workflow_summary.get("protocol", "unknown")),
            "source_read_layout": str(workflow_summary.get("read_layout", "unknown")),
            "matrices_present": matrices_present,
            "mudata_present": bool(mudata_paths),
            "mudata_files": mudata_paths,
            "qc_present": bool(qc_files),
            "qc_files": qc_files,
            "has_workflow_summary": workflow_summary_path.exists(),
        },
        command_name="circyto inspect-workdir",
        workflow_type="workdir-inspection",
        protocol=str(workflow_summary.get("protocol", "unknown")),
        read_layout=str(workflow_summary.get("read_layout", "unknown")),
        genome_fasta=str(workflow_summary.get("genome_fasta")) if workflow_summary.get("genome_fasta") else None,
        gtf=str(workflow_summary.get("gtf")) if workflow_summary.get("gtf") else None,
        detector_backend=str(workflow_summary.get("detector_backend")) if workflow_summary.get("detector_backend") else None,
        started_at=utc_now_iso(),
        completed_at=utc_now_iso(),
        elapsed_seconds=0.0,
        workflow_uuid=str(workflow_summary.get("workflow_uuid")) if workflow_summary.get("workflow_uuid") else None,
    )


def validate_completed_workdir(workdir: Path) -> dict[str, Any]:
    base = check_workflow_integrity(workdir)
    errors = list(base.get("errors", []))
    warnings = list(base.get("warnings", []))
    details = dict(base.get("details", {}))

    matrix_dir = workdir / "matrix"
    qc_dir = workdir / "qc"
    anndata_path = workdir / "anndata" / "circ_counts.h5ad"
    mudata_path = workdir / "mudata" / "full_length.h5mu"
    circ_feature_path = matrix_dir / "circ_feature_table.tsv"

    if not matrix_dir.exists():
        errors.append(f"Missing matrix directory: {matrix_dir}")

    if circ_feature_path.exists():
        circ_feature_df = pd.read_csv(circ_feature_path, sep="\t", keep_default_na=False)
        required_circ_columns = {"circ_id", "chrom", "start", "end", "strand"}
        missing_circ_columns = sorted(required_circ_columns - set(circ_feature_df.columns))
        if missing_circ_columns:
            errors.append(
                f"circ_feature_table.tsv missing required columns: {', '.join(missing_circ_columns)}"
            )
    else:
        warnings.append(f"Optional circ feature table missing: {circ_feature_path}")

    if (workdir / "rna").exists():
        for name in ("gene_counts.tsv", "gene_feature_table.tsv", "rna_import_summary.json"):
            path = workdir / "rna" / name
            if not path.exists():
                errors.append(f"RNA outputs enabled but missing {name}")
        if not (qc_dir / "rna_qc.tsv").exists():
            warnings.append(f"Optional RNA QC summary missing: {qc_dir / 'rna_qc.tsv'}")
        if not (qc_dir / "rna_gene_qc.tsv").exists():
            warnings.append(f"Optional RNA gene QC summary missing: {qc_dir / 'rna_gene_qc.tsv'}")

    if anndata_path.exists():
        if not HAS_ANNDATA:
            warnings.append("anndata not installed; could not validate circ_counts.h5ad readability")
        else:
            try:
                adata = ad.read_h5ad(str(anndata_path))
                details["h5ad_shape"] = [int(adata.n_obs), int(adata.n_vars)]
            except Exception as exc:
                errors.append(f"h5ad unreadable: {anndata_path} ({exc})")
    else:
        warnings.append(f"Optional h5ad missing: {anndata_path}")

    if mudata_path.exists():
        if not HAS_MUDATA:
            warnings.append("mudata not installed; could not validate full_length.h5mu readability")
        else:
            try:
                mdata = mu.read_h5mu(str(mudata_path))
                details["h5mu_modalities"] = list(mdata.mod.keys())
                details["h5mu_n_obs"] = int(mdata.n_obs)
                if "circ" in mdata.mod:
                    circ_var_columns = list(mdata.mod["circ"].var.columns)
                    details["h5mu_circ_var_columns"] = circ_var_columns
                    missing = [column for column in ("chrom", "start", "end", "strand") if column not in circ_var_columns]
                    if missing:
                        errors.append(f"h5mu circ modality missing expected var columns: {', '.join(missing)}")
                if "rna" in mdata.mod and "circ" in mdata.mod:
                    if list(mdata.mod["rna"].obs_names) != list(mdata.mod["circ"].obs_names):
                        errors.append("h5mu obs inconsistency: RNA and circ obs_names differ")
            except Exception as exc:
                errors.append(f"h5mu unreadable: {mudata_path} ({exc})")
    else:
        warnings.append(f"Optional h5mu missing: {mudata_path}")

    if (qc_dir / "rna_circ_summary.json").exists() and not (qc_dir / "rna_circ_cell_summary.tsv").exists():
        warnings.append("rna_circ_summary.json exists but rna_circ_cell_summary.tsv is missing")
    if (workdir / "dna").exists():
        if not (workdir / "dna" / "dna_snv_import_summary.json").exists():
            warnings.append(f"DNA directory present but missing dna_snv_import_summary.json")

    return {
        "ok": len(errors) == 0,
        "workdir": str(workdir.resolve()),
        "errors": errors,
        "warnings": warnings,
        "details": details,
    }


def summarize_circ_host_genes(
    *,
    workdir: Path,
    output_path: Path | None = None,
) -> dict[str, Any]:
    matrix_dir = workdir / "matrix"
    matrix_path = matrix_dir / "circ_counts.mtx"
    circ_index_path = matrix_dir / "circ_index.txt"
    cell_index_path = matrix_dir / "cell_index.txt"
    feature_path = matrix_dir / "circ_feature_table.tsv"
    if not matrix_path.exists():
        raise FileNotFoundError(f"Missing circ matrix: {matrix_path}")
    if not circ_index_path.exists():
        raise FileNotFoundError(f"Missing circ index: {circ_index_path}")
    if not cell_index_path.exists():
        raise FileNotFoundError(f"Missing cell index: {cell_index_path}")
    if not feature_path.exists():
        raise FileNotFoundError(f"Missing circ feature table: {feature_path}")

    X, circ_ids, cell_ids = load_circ_matrix(
        matrix_path=matrix_path,
        circ_index_path=circ_index_path,
        cell_index_path=cell_index_path,
    )
    X_cells_by_circ = _orient_circ_matrix_cells_by_circ(X=X, circ_ids=circ_ids, cell_ids=cell_ids, matrix_path=matrix_path)
    feature_df = load_circ_feature_table(circ_ids, feature_path).reset_index(names="circ_id")
    feature_df["host_gene"] = feature_df.get("host_gene", pd.Series([""] * len(feature_df))).fillna("").astype(str)

    total_support = np.asarray(X_cells_by_circ.sum(axis=0)).ravel()
    n_cells_detected = np.asarray(X_cells_by_circ.getnnz(axis=0)).ravel()
    feature_df["total_support"] = total_support.astype(int)
    feature_df["n_cells_detected"] = n_cells_detected.astype(int)
    nonempty = feature_df[feature_df["host_gene"].astype(str).str.strip() != ""].copy()
    if nonempty.empty:
        raise ValueError(f"{feature_path} contains no non-empty host_gene annotations.")

    host_summary = (
        nonempty.groupby("host_gene", dropna=False)
        .agg(
            n_circRNAs=("circ_id", "nunique"),
            total_support=("total_support", "sum"),
            n_cells_detected_sum=("n_cells_detected", "sum"),
        )
        .reset_index()
        .sort_values(["total_support", "n_circRNAs", "host_gene"], ascending=[False, False, True])
    )
    if output_path is None:
        output_path = workdir / "qc" / "circ_host_gene_summary.tsv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    host_summary.to_csv(output_path, sep="\t", index=False)

    top_host_genes = host_summary.head(10).to_dict(orient="records")
    return apply_standard_provenance(
        {
            "workdir": str(workdir.resolve()),
            "output_path": str(output_path.resolve()),
            "n_host_genes": int(host_summary.shape[0]),
            "n_circ_with_host_gene": int(nonempty.shape[0]),
            "top_host_genes": top_host_genes,
        },
        command_name="circyto summarize-circ-host-genes",
        workflow_type="circ-host-gene-summary",
        protocol="unknown",
        read_layout="unknown",
        genome_fasta=None,
        gtf=None,
        detector_backend=None,
        started_at=utc_now_iso(),
        completed_at=utc_now_iso(),
        elapsed_seconds=0.0,
    )


def export_circ_bed(
    *,
    workdir: Path,
    output_path: Path | None = None,
) -> dict[str, Any]:
    matrix_dir = workdir / "matrix"
    matrix_path = matrix_dir / "circ_counts.mtx"
    circ_index_path = matrix_dir / "circ_index.txt"
    cell_index_path = matrix_dir / "cell_index.txt"
    feature_path = matrix_dir / "circ_feature_table.tsv"
    if not feature_path.exists():
        raise FileNotFoundError(f"Missing circ feature table: {feature_path}")
    if not matrix_path.exists():
        raise FileNotFoundError(f"Missing circ matrix: {matrix_path}")
    if not circ_index_path.exists():
        raise FileNotFoundError(f"Missing circ index: {circ_index_path}")
    if not cell_index_path.exists():
        raise FileNotFoundError(f"Missing cell index: {cell_index_path}")

    X, circ_ids, cell_ids = load_circ_matrix(
        matrix_path=matrix_path,
        circ_index_path=circ_index_path,
        cell_index_path=cell_index_path,
    )
    X_cells_by_circ = _orient_circ_matrix_cells_by_circ(X=X, circ_ids=circ_ids, cell_ids=cell_ids, matrix_path=matrix_path)
    feature_df = load_circ_feature_table(circ_ids, feature_path).reset_index(names="circ_id")
    if output_path is None:
        output_path = workdir / "qc" / "circs.bed"

    missing_columns = [column for column in ("chrom", "start", "end") if column not in feature_df.columns]
    if missing_columns:
        raise ValueError(f"{feature_path} is missing required BED columns: {', '.join(missing_columns)}")
    feature_df["support"] = np.asarray(X_cells_by_circ.sum(axis=0)).ravel().astype(int)
    bed_df = feature_df[["chrom", "start", "end", "circ_id", "support"]].copy()
    bed_df["start"] = pd.to_numeric(bed_df["start"], errors="coerce")
    bed_df["end"] = pd.to_numeric(bed_df["end"], errors="coerce")
    if bed_df[["chrom", "start", "end"]].isna().any().any():
        raise ValueError(f"{feature_path} contains incomplete chrom/start/end values required for BED export.")
    bed_df["start"] = bed_df["start"].astype(int)
    bed_df["end"] = bed_df["end"].astype(int)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    bed_df.to_csv(output_path, sep="\t", index=False, header=False)
    return apply_standard_provenance(
        {
            "workdir": str(workdir.resolve()),
            "output_path": str(output_path.resolve()),
            "n_rows": int(bed_df.shape[0]),
            "columns": ["chrom", "start", "end", "circ_id", "support"],
        },
        command_name="circyto export-circ-bed",
        workflow_type="circ-bed-export",
        protocol="unknown",
        read_layout="unknown",
        genome_fasta=None,
        gtf=None,
        detector_backend=None,
        started_at=utc_now_iso(),
        completed_at=utc_now_iso(),
        elapsed_seconds=0.0,
    )


def _infer_dataset_name_from_workdir(workdir: Path) -> str:
    joined = "/".join(part.lower() for part in workdir.parts)
    if "emtab8735" in joined or "all192" in joined or "smartseq3" in joined:
        return "E-MTAB-8735 Smart-seq3"
    if "imr90" in joined or "gse278958" in joined or "scrr_imr90" in joined:
        return "GSE278958 IMR90 scRR"
    if "hap1" in joined or "gse278952" in joined or "scrr_hap1" in joined:
        return "GSE278952 HAP1 scRR"
    return workdir.name


def summarize_benchmark_workdirs(
    *,
    workdirs: list[Path],
    output_tsv: Path | None = None,
    output_json: Path | None = None,
) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    for workdir in workdirs:
        joined_df, _ = _build_rna_circ_joined_table(workdir)
        workflow_summary_path = workdir / "workflow_summary.json"
        workflow_summary = load_json(workflow_summary_path) if workflow_summary_path.exists() else {}
        mudata_path = workdir / "mudata" / "full_length.h5mu"
        gene_counts_path = workdir / "rna" / "gene_counts.tsv"
        n_rna_features = 0
        if gene_counts_path.exists():
            n_rna_features = int(pd.read_csv(gene_counts_path, sep="\t", keep_default_na=False).shape[0])
        circ_ids = read_index_lines(workdir / "matrix" / "circ_index.txt")
        rows.append(
            {
                "workdir": str(workdir.resolve()),
                "dataset_name": _infer_dataset_name_from_workdir(workdir),
                "workflow_type": str(workflow_summary.get("workflow_type", workflow_summary.get("workflow", "unknown"))),
                "protocol": str(workflow_summary.get("protocol", "unknown")),
                "read_layout": str(workflow_summary.get("read_layout", "unknown")),
                "n_cells": int(joined_df.shape[0]),
                "n_rna_features": int(n_rna_features),
                "n_circ_features": int(len(circ_ids)),
                "median_rna_counts": float(joined_df["total_rna_counts"].median()),
                "median_detected_genes": float(joined_df["detected_genes"].median()),
                "median_circRNA_count": float(joined_df["circRNA_count"].median()),
                "median_circRNA_total_support": float(joined_df["circRNA_total_support"].median()),
                "h5mu_exists": mudata_path.exists(),
                "h5mu_size_bytes": int(mudata_path.stat().st_size) if mudata_path.exists() else 0,
                "workdir_size_bytes": directory_size_bytes(workdir),
                "cleanup_status": str(workflow_summary.get("cleanup_summary", {}).get("performed", workflow_summary.get("cleanup_performed", False))),
                "workflow_succeeded": _workflow_is_marked_successful(workflow_summary),
            }
        )
    summary_df = pd.DataFrame(rows)
    if output_tsv is not None:
        output_tsv.parent.mkdir(parents=True, exist_ok=True)
        summary_df.to_csv(output_tsv, sep="\t", index=False)
    payload = {
        "n_workdirs": int(len(rows)),
        "columns": list(summary_df.columns),
        "rows": summary_df.to_dict(orient="records"),
    }
    if output_tsv is not None:
        payload["output_tsv"] = str(output_tsv.resolve())
    if output_json is not None:
        output_json.parent.mkdir(parents=True, exist_ok=True)
        write_json(output_json, payload)
        payload["output_json"] = str(output_json.resolve())
    return apply_standard_provenance(
        payload,
        command_name="circyto summarize-benchmark",
        workflow_type="benchmark-summary",
        protocol="mixed",
        read_layout="mixed",
        genome_fasta=None,
        gtf=None,
        detector_backend=None,
        started_at=utc_now_iso(),
        completed_at=utc_now_iso(),
        elapsed_seconds=0.0,
    )


DNA_CELL_SUMMARY_REQUIRED = ["cell_id", "dna_library_id", "cnv_burden", "replication_score", "cell_cycle_phase", "dna_variant_count", "notes"]
DNA_VARIANT_SUMMARY_REQUIRED = ["variant_id", "cell_id", "chrom", "pos", "ref", "alt", "gene", "consequence", "evidence_type", "caller", "filter_status"]
SCOMATIC_CANDIDATE_REQUIRED = ["variant_id", "cell_id", "chrom", "pos", "ref", "alt", "gene", "filter_status", "candidate_variant_class", "read_support", "vaf", "caller"]


def _validate_required_columns(path: Path, df: pd.DataFrame, required: list[str]) -> None:
    missing = [column for column in required if column not in df.columns]
    if missing:
        raise ValueError(f"{path} is missing required columns: {', '.join(missing)}")


def _read_required_tsv(path: Path, required: list[str]) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", keep_default_na=False)
    _validate_required_columns(path, df, required)
    return df


def import_dna_snv_summary(
    *,
    workdir: Path,
    dna_cell_summary_path: Path,
    dna_variant_summary_path: Path | None = None,
    scomatic_candidate_summary_path: Path | None = None,
) -> dict[str, Any]:
    dna_dir = workdir / "dna"
    dna_dir.mkdir(parents=True, exist_ok=True)
    dna_cell_df = _read_required_tsv(dna_cell_summary_path, DNA_CELL_SUMMARY_REQUIRED)
    validate_feature_id_uniqueness(dna_cell_df["cell_id"].astype(str).tolist(), label="cell")

    expected_cells: list[str] = []
    try:
        joined_df, _ = _build_rna_circ_joined_table(workdir)
        expected_cells = joined_df["cell_id"].astype(str).tolist()
    except Exception:
        expected_cells = []
    dna_cells = dna_cell_df["cell_id"].astype(str).tolist()
    shared = sorted(set(expected_cells) & set(dna_cells))
    rna_circ_only = sorted(set(expected_cells) - set(dna_cells))
    dna_only = sorted(set(dna_cells) - set(expected_cells))

    dna_variant_df = _read_required_tsv(dna_variant_summary_path, DNA_VARIANT_SUMMARY_REQUIRED) if dna_variant_summary_path is not None else None
    scomatic_df = _read_required_tsv(scomatic_candidate_summary_path, SCOMATIC_CANDIDATE_REQUIRED) if scomatic_candidate_summary_path is not None else None

    dna_cell_out = dna_dir / "dna_cell_summary.tsv"
    dna_cell_df.to_csv(dna_cell_out, sep="\t", index=False)
    dna_variant_out = None
    if dna_variant_df is not None:
        dna_variant_out = dna_dir / "dna_variant_summary.tsv"
        dna_variant_df.to_csv(dna_variant_out, sep="\t", index=False)
    scomatic_out = None
    if scomatic_df is not None:
        scomatic_out = dna_dir / "scomatic_candidate_summary.tsv"
        scomatic_df.to_csv(scomatic_out, sep="\t", index=False)

    summary = {
        "workdir": str(workdir.resolve()),
        "dna_cell_summary_rows": int(dna_cell_df.shape[0]),
        "dna_variant_summary_rows": int(dna_variant_df.shape[0]) if dna_variant_df is not None else 0,
        "scomatic_candidate_summary_rows": int(scomatic_df.shape[0]) if scomatic_df is not None else 0,
        "shared_cells": shared,
        "n_shared_cells": int(len(shared)),
        "rna_circ_only_cells": rna_circ_only,
        "n_rna_circ_only_cells": int(len(rna_circ_only)),
        "dna_only_cells": dna_only,
        "n_dna_only_cells": int(len(dna_only)),
        "terminology_note": "SComatic outputs are treated as RNA-derived candidate variant signals, not validated somatic mutations.",
        "output_dna_cell_summary": str(dna_cell_out.resolve()),
        "output_dna_variant_summary": str(dna_variant_out.resolve()) if dna_variant_out is not None else None,
        "output_scomatic_candidate_summary": str(scomatic_out.resolve()) if scomatic_out is not None else None,
    }
    summary_path = dna_dir / "dna_snv_import_summary.json"
    write_json(summary_path, summary)
    summary["output_dna_snv_import_summary"] = str(summary_path.resolve())
    return apply_standard_provenance(
        summary,
        command_name="circyto import-dna-snv-summary",
        workflow_type="dna-snv-import",
        protocol="unknown",
        read_layout="unknown",
        genome_fasta=None,
        gtf=None,
        detector_backend=None,
        started_at=utc_now_iso(),
        completed_at=utc_now_iso(),
        elapsed_seconds=0.0,
    )


def summarize_dna_rna_circ(
    *,
    workdir: Path,
    write_summary: bool = False,
) -> dict[str, Any]:
    joined_df, _ = _build_rna_circ_joined_table(workdir)
    dna_cell_path = workdir / "dna" / "dna_cell_summary.tsv"
    if not dna_cell_path.exists():
        raise FileNotFoundError(f"Missing DNA cell summary: {dna_cell_path}")
    dna_df = _read_required_tsv(dna_cell_path, DNA_CELL_SUMMARY_REQUIRED)
    dna_df["cell_id"] = dna_df["cell_id"].astype(str)
    scomatic_counts: dict[str, int] = {}
    scomatic_path = workdir / "dna" / "scomatic_candidate_summary.tsv"
    if scomatic_path.exists():
        scomatic_df = _read_required_tsv(scomatic_path, SCOMATIC_CANDIDATE_REQUIRED)
        scomatic_counts = scomatic_df.groupby("cell_id").size().to_dict()

    rna_circ_cells = list(joined_df["cell_id"].astype(str))
    dna_cells = list(dna_df["cell_id"].astype(str))
    all_cells = rna_circ_cells + [cell_id for cell_id in dna_cells if cell_id not in set(rna_circ_cells)]
    dna_indexed = dna_df.set_index("cell_id", drop=False)
    rows: list[dict[str, Any]] = []
    for cell_id in all_cells:
        rna_match = joined_df[joined_df["cell_id"] == cell_id]
        if not rna_match.empty:
            rna_row = rna_match.iloc[0]
            total_rna_counts = int(rna_row.get("total_rna_counts", 0))
            detected_genes = int(rna_row.get("detected_genes", 0))
            circ_count = int(rna_row.get("circRNA_count", 0))
            circ_support = int(rna_row.get("circRNA_total_support", 0))
        else:
            total_rna_counts = detected_genes = circ_count = circ_support = 0
        if cell_id in dna_indexed.index:
            dna_row = dna_indexed.loc[cell_id]
            if isinstance(dna_row, pd.DataFrame):
                dna_row = dna_row.iloc[0]
            cnv_burden = float(dna_row.get("cnv_burden", 0) or 0)
            replication_score = float(dna_row.get("replication_score", 0) or 0)
            cell_cycle_phase = str(dna_row.get("cell_cycle_phase", ""))
            dna_variant_count = int(float(dna_row.get("dna_variant_count", 0) or 0))
        else:
            cnv_burden = replication_score = 0.0
            cell_cycle_phase = ""
            dna_variant_count = 0
        membership = "shared" if cell_id in set(rna_circ_cells) and cell_id in set(dna_cells) else ("rna_circ_only" if cell_id in set(rna_circ_cells) else "dna_only")
        rows.append(
            {
                "cell_id": cell_id,
                "membership": membership,
                "total_rna_counts": total_rna_counts,
                "detected_genes": detected_genes,
                "circRNA_count": circ_count,
                "circRNA_total_support": circ_support,
                "cnv_burden": cnv_burden,
                "replication_score": replication_score,
                "cell_cycle_phase": cell_cycle_phase,
                "dna_variant_count": dna_variant_count,
                "scomatic_candidate_count": int(scomatic_counts.get(cell_id, 0)),
            }
        )
    summary_df = pd.DataFrame(rows)
    summary = {
        "workdir": str(workdir.resolve()),
        "n_joined_cells": int(summary_df.shape[0]),
        "n_shared_cells": int((summary_df["membership"] == "shared").sum()),
        "n_rna_circ_only_cells": int((summary_df["membership"] == "rna_circ_only").sum()),
        "n_dna_only_cells": int((summary_df["membership"] == "dna_only").sum()),
        "terminology_note": "SComatic counts are RNA-derived candidate variant signals, not validated somatic mutations.",
        "joined_preview": summary_df.head(10).to_dict(orient="records"),
    }
    if write_summary:
        qc_dir = workdir / "qc"
        qc_dir.mkdir(parents=True, exist_ok=True)
        cell_summary_path = qc_dir / "dna_rna_circ_cell_summary.tsv"
        summary_json_path = qc_dir / "dna_rna_circ_summary.json"
        summary_df.to_csv(cell_summary_path, sep="\t", index=False)
        write_json(summary_json_path, summary)
        summary["output_dna_rna_circ_cell_summary"] = str(cell_summary_path.resolve())
        summary["output_dna_rna_circ_summary"] = str(summary_json_path.resolve())
    return apply_standard_provenance(
        summary,
        command_name="circyto summarize-dna-rna-circ",
        workflow_type="dna-rna-circ-summary",
        protocol="unknown",
        read_layout="unknown",
        genome_fasta=None,
        gtf=None,
        detector_backend=None,
        started_at=utc_now_iso(),
        completed_at=utc_now_iso(),
        elapsed_seconds=0.0,
    )


def _orient_circ_matrix_cells_by_circ(
    *,
    X: sp.csr_matrix,
    circ_ids: list[str],
    cell_ids: list[str],
    matrix_path: Path,
) -> sp.csr_matrix:
    if X.shape == (len(circ_ids), len(cell_ids)):
        return X.T.tocsr()
    if X.shape == (len(cell_ids), len(circ_ids)):
        return X.tocsr()
    raise ValueError(
        f"Circ matrix shape mismatch for {matrix_path}: matrix={X.shape}, "
        f"circ_index={len(circ_ids)}, cell_index={len(cell_ids)}. "
        "Expected either circRNAs x cells or cells x circRNAs."
    )


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
    per_cell, per_gene, summary = _build_rna_qc_tables(
        counts_df=counts_df,
        feature_df=feature_df,
        work_root=work_root,
    )

    cell_qc_path = qc_dir / "rna_qc.tsv"
    gene_qc_path = qc_dir / "rna_gene_qc.tsv"
    per_cell.to_csv(cell_qc_path, sep="\t", index=False)
    per_gene.to_csv(gene_qc_path, sep="\t", index=False)

    return {
        "output_rna_qc": str(cell_qc_path.resolve()),
        "output_rna_gene_qc": str(gene_qc_path.resolve()),
        **summary,
    }


def _build_rna_qc_tables(
    *,
    counts_df: pd.DataFrame,
    feature_df: pd.DataFrame,
    work_root: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
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
    circ_metrics = _circ_metrics_by_cell(work_root)
    circ_counts = dict(circ_metrics.get("circRNA_count", {}))
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

    genes_detected = gene_detected.astype(int)
    summary = {
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
    return per_cell, per_gene, summary


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


def summarize_rna_circ_integration(
    *,
    workdir: Path,
    write_summary: bool,
) -> dict[str, Any]:
    rna_dir = workdir / "rna"
    qc_dir = workdir / "qc"
    gene_counts_path = rna_dir / "gene_counts.tsv"
    if not gene_counts_path.exists():
        raise FileNotFoundError(f"Missing RNA gene-count table: {gene_counts_path}")
    counts_df = pd.read_csv(gene_counts_path, sep="\t", keep_default_na=False)
    required_columns = {"gene_id", "gene_name"}
    missing = sorted(required_columns - set(counts_df.columns))
    if missing:
        raise ValueError(f"{gene_counts_path} is missing required columns: {', '.join(missing)}")

    rna_qc_path = qc_dir / "rna_qc.tsv"
    feature_path = rna_dir / "gene_feature_table.tsv"
    if rna_qc_path.exists():
        rna_qc_df = pd.read_csv(rna_qc_path, sep="\t", keep_default_na=False)
    else:
        if not feature_path.exists():
            raise FileNotFoundError(
                f"Missing RNA gene feature table for RNA QC fallback: {feature_path}"
            )
        feature_df = pd.read_csv(feature_path, sep="\t", keep_default_na=False)
        rna_qc_df, _, _ = _build_rna_qc_tables(
            counts_df=counts_df,
            feature_df=feature_df,
            work_root=workdir,
        )

    if "cell_id" not in rna_qc_df.columns:
        raise ValueError(f"{rna_qc_path if rna_qc_path.exists() else 'derived RNA QC'} must contain a cell_id column.")

    rna_cells = [str(cell_id) for cell_id in rna_qc_df["cell_id"].astype(str).tolist()]
    validate_feature_id_uniqueness(rna_cells, label="cell")

    circ_metrics = _circ_metrics_by_cell(workdir)
    circ_counts = dict(circ_metrics.get("circRNA_count", {}))
    circ_support = dict(circ_metrics.get("circRNA_total_support", {}))
    circ_cells = _matrix_cell_ids(workdir)

    shared_cells = sorted(set(rna_cells) & set(circ_cells))
    rna_only_cells = sorted(set(rna_cells) - set(circ_cells))
    circ_only_cells = sorted(set(circ_cells) - set(rna_cells))
    all_cells = list(rna_cells) + [cell_id for cell_id in circ_cells if cell_id not in set(rna_cells)]

    rna_qc_indexed = rna_qc_df.set_index("cell_id", drop=False)
    rows: list[dict[str, Any]] = []
    for cell_id in all_cells:
        if cell_id in rna_qc_indexed.index:
            row = rna_qc_indexed.loc[cell_id]
            if isinstance(row, pd.DataFrame):
                row = row.iloc[0]
            total_counts = int(float(row.get("total_counts", 0) or 0))
            detected_genes = int(float(row.get("detected_genes", 0) or 0))
            mitochondrial_fraction = float(row.get("mitochondrial_fraction", 0.0) or 0.0)
            ribosomal_fraction = float(row.get("ribosomal_fraction", 0.0) or 0.0)
        else:
            total_counts = 0
            detected_genes = 0
            mitochondrial_fraction = 0.0
            ribosomal_fraction = 0.0
        circ_count = int(circ_counts.get(cell_id, 0))
        circ_total_support = int(circ_support.get(cell_id, 0))
        if cell_id in shared_cells:
            membership = "shared"
        elif cell_id in rna_only_cells:
            membership = "rna_only"
        else:
            membership = "circ_only"
        rows.append(
            {
                "cell_id": cell_id,
                "membership": membership,
                "total_rna_counts": total_counts,
                "detected_genes": detected_genes,
                "mitochondrial_fraction": mitochondrial_fraction,
                "ribosomal_fraction": ribosomal_fraction,
                "circRNA_count": circ_count,
                "circRNA_total_support": circ_total_support,
            }
        )

    joined_df = pd.DataFrame(rows)
    relationship = {
        "shared_cells_considered": int(len(shared_cells)),
        "pearson_correlation_total_rna_vs_circ_count": None,
    }
    if len(shared_cells) >= 2:
        shared_df = joined_df[joined_df["membership"] == "shared"]
        if shared_df["total_rna_counts"].nunique() > 1 and shared_df["circRNA_count"].nunique() > 1:
            relationship["pearson_correlation_total_rna_vs_circ_count"] = float(
                shared_df["total_rna_counts"].corr(shared_df["circRNA_count"])
            )

    summary = {
        "workdir": str(workdir.resolve()),
        "n_rna_cells": int(len(rna_cells)),
        "n_circ_cells": int(len(circ_cells)),
        "n_shared_cells": int(len(shared_cells)),
        "n_rna_only_cells": int(len(rna_only_cells)),
        "n_circ_only_cells": int(len(circ_only_cells)),
        "rna_only_cells": rna_only_cells,
        "circ_only_cells": circ_only_cells,
        "rna_total_count_vs_circRNA_count_relationship": relationship,
        "joined_rows": int(joined_df.shape[0]),
        "joined_columns": list(joined_df.columns),
        "joined_preview": joined_df.head(10).to_dict(orient="records"),
    }

    if write_summary:
        qc_dir.mkdir(parents=True, exist_ok=True)
        cell_summary_path = qc_dir / "rna_circ_cell_summary.tsv"
        summary_json_path = qc_dir / "rna_circ_summary.json"
        joined_df.to_csv(cell_summary_path, sep="\t", index=False)
        write_json(summary_json_path, summary)
        summary["output_rna_circ_cell_summary"] = str(cell_summary_path.resolve())
        summary["output_rna_circ_summary"] = str(summary_json_path.resolve())

    return apply_standard_provenance(
        summary,
        command_name="circyto summarize-rna-circ",
        workflow_type="rna-circ-integration-summary",
        protocol="unknown",
        read_layout="unknown",
        genome_fasta=None,
        gtf=None,
        detector_backend=None,
        started_at=utc_now_iso(),
        completed_at=utc_now_iso(),
        elapsed_seconds=0.0,
    )


def _build_rna_circ_joined_table(workdir: Path) -> tuple[pd.DataFrame, dict[str, Any]]:
    rna_dir = workdir / "rna"
    qc_dir = workdir / "qc"
    gene_counts_path = rna_dir / "gene_counts.tsv"
    if not gene_counts_path.exists():
        raise FileNotFoundError(f"Missing RNA gene-count table: {gene_counts_path}")
    counts_df = pd.read_csv(gene_counts_path, sep="\t", keep_default_na=False)
    required_columns = {"gene_id", "gene_name"}
    missing = sorted(required_columns - set(counts_df.columns))
    if missing:
        raise ValueError(f"{gene_counts_path} is missing required columns: {', '.join(missing)}")

    rna_qc_path = qc_dir / "rna_qc.tsv"
    feature_path = rna_dir / "gene_feature_table.tsv"
    if rna_qc_path.exists():
        rna_qc_df = pd.read_csv(rna_qc_path, sep="\t", keep_default_na=False)
    else:
        if not feature_path.exists():
            raise FileNotFoundError(
                f"Missing RNA gene feature table for RNA QC fallback: {feature_path}"
            )
        feature_df = pd.read_csv(feature_path, sep="\t", keep_default_na=False)
        rna_qc_df, _, _ = _build_rna_qc_tables(
            counts_df=counts_df,
            feature_df=feature_df,
            work_root=workdir,
        )

    if "cell_id" not in rna_qc_df.columns:
        raise ValueError(f"{rna_qc_path if rna_qc_path.exists() else 'derived RNA QC'} must contain a cell_id column.")

    rna_cells = [str(cell_id) for cell_id in rna_qc_df["cell_id"].astype(str).tolist()]
    validate_feature_id_uniqueness(rna_cells, label="cell")

    circ_metrics = _circ_metrics_by_cell(workdir)
    circ_counts = dict(circ_metrics.get("circRNA_count", {}))
    circ_support = dict(circ_metrics.get("circRNA_total_support", {}))
    circ_cells = _matrix_cell_ids(workdir)

    shared_cells = sorted(set(rna_cells) & set(circ_cells))
    rna_only_cells = sorted(set(rna_cells) - set(circ_cells))
    circ_only_cells = sorted(set(circ_cells) - set(rna_cells))
    all_cells = list(rna_cells) + [cell_id for cell_id in circ_cells if cell_id not in set(rna_cells)]

    rna_qc_indexed = rna_qc_df.set_index("cell_id", drop=False)
    rows: list[dict[str, Any]] = []
    for cell_id in all_cells:
        if cell_id in rna_qc_indexed.index:
            row = rna_qc_indexed.loc[cell_id]
            if isinstance(row, pd.DataFrame):
                row = row.iloc[0]
            total_counts = int(float(row.get("total_counts", 0) or 0))
            detected_genes = int(float(row.get("detected_genes", 0) or 0))
            mitochondrial_fraction = float(row.get("mitochondrial_fraction", 0.0) or 0.0)
            ribosomal_fraction = float(row.get("ribosomal_fraction", 0.0) or 0.0)
        else:
            total_counts = 0
            detected_genes = 0
            mitochondrial_fraction = 0.0
            ribosomal_fraction = 0.0
        circ_count = int(circ_counts.get(cell_id, 0))
        circ_total_support = int(circ_support.get(cell_id, 0))
        if cell_id in shared_cells:
            membership = "shared"
        elif cell_id in rna_only_cells:
            membership = "rna_only"
        else:
            membership = "circ_only"
        rows.append(
            {
                "cell_id": cell_id,
                "membership": membership,
                "total_rna_counts": total_counts,
                "detected_genes": detected_genes,
                "mitochondrial_fraction": mitochondrial_fraction,
                "ribosomal_fraction": ribosomal_fraction,
                "circRNA_count": circ_count,
                "circRNA_total_support": circ_total_support,
            }
        )

    joined_df = pd.DataFrame(rows)
    relationship = {
        "shared_cells_considered": int(len(shared_cells)),
        "pearson_correlation_total_rna_vs_circ_count": None,
    }
    if len(shared_cells) >= 2:
        shared_df = joined_df[joined_df["membership"] == "shared"]
        if shared_df["total_rna_counts"].nunique() > 1 and shared_df["circRNA_count"].nunique() > 1:
            relationship["pearson_correlation_total_rna_vs_circ_count"] = float(
                shared_df["total_rna_counts"].corr(shared_df["circRNA_count"])
            )

    summary = {
        "workdir": str(workdir.resolve()),
        "n_rna_cells": int(len(rna_cells)),
        "n_circ_cells": int(len(circ_cells)),
        "n_shared_cells": int(len(shared_cells)),
        "n_rna_only_cells": int(len(rna_only_cells)),
        "n_circ_only_cells": int(len(circ_only_cells)),
        "rna_only_cells": rna_only_cells,
        "circ_only_cells": circ_only_cells,
        "rna_total_count_vs_circRNA_count_relationship": relationship,
        "joined_rows": int(joined_df.shape[0]),
        "joined_columns": list(joined_df.columns),
        "joined_preview": joined_df.head(10).to_dict(orient="records"),
    }
    return joined_df, summary


def export_completed_workflow_mudata(
    *,
    workdir: Path,
    out_path: Path | None = None,
    overwrite: bool = False,
) -> dict[str, Any]:
    if out_path is None:
        out_path = workdir / "mudata" / "full_length.h5mu"
    if out_path.exists() and not overwrite:
        raise ValueError(f"Output already exists: {out_path}. Use --overwrite to replace it.")

    rna_dir = workdir / "rna"
    matrix_dir = workdir / "matrix"
    gene_counts_path = rna_dir / "gene_counts.tsv"
    gene_feature_path = rna_dir / "gene_feature_table.tsv"
    circ_matrix_path = matrix_dir / "circ_counts.mtx"
    circ_index_path = matrix_dir / "circ_index.txt"
    cell_index_path = matrix_dir / "cell_index.txt"
    circ_feature_path = matrix_dir / "circ_feature_table.tsv"
    if not gene_counts_path.exists():
        raise FileNotFoundError(f"Missing RNA gene-count table: {gene_counts_path}")
    if not gene_feature_path.exists():
        raise FileNotFoundError(f"Missing RNA gene feature table: {gene_feature_path}")
    if not circ_matrix_path.exists():
        raise FileNotFoundError(f"Missing circ matrix: {circ_matrix_path}")
    if not circ_index_path.exists():
        raise FileNotFoundError(f"Missing circ index: {circ_index_path}")
    if not cell_index_path.exists():
        raise FileNotFoundError(f"Missing circ cell index: {cell_index_path}")
    if not HAS_ANNDATA:
        raise RuntimeError("anndata is required for MuData export")
    if not HAS_MUDATA:
        raise RuntimeError("mudata is not installed; install circyto[mudata] or pip install mudata")

    joined_df, rna_circ_summary = _build_rna_circ_joined_table(workdir)
    joined_df = joined_df.set_index("cell_id", drop=True)

    rna_counts_df = pd.read_csv(gene_counts_path, sep="\t", keep_default_na=False)
    rna_feature_df = pd.read_csv(gene_feature_path, sep="\t", keep_default_na=False)
    validate_feature_id_uniqueness(rna_counts_df["gene_id"].astype(str).tolist(), label="gene")
    rna_cell_ids = [column for column in rna_counts_df.columns if column not in {"gene_id", "gene_name"}]
    rna_count_matrix = sp.csr_matrix(rna_counts_df[rna_cell_ids].apply(pd.to_numeric, errors="raise").to_numpy().T)
    rna_var = rna_feature_df.copy()
    if "gene_biotype" not in rna_var.columns:
        rna_var["gene_biotype"] = ""
    if "gene_id" in rna_var.columns:
        rna_var = rna_var.set_index("gene_id", drop=False)
    else:
        raise ValueError(f"{gene_feature_path} must contain a gene_id column.")
    if list(rna_var["gene_id"].astype(str)) != list(rna_counts_df["gene_id"].astype(str)):
        rna_var = rna_var.reindex(rna_counts_df["gene_id"].astype(str).tolist())
    union_cell_ids = list(joined_df.index.astype(str))
    rna_pos = {cell_id: idx for idx, cell_id in enumerate(rna_cell_ids)}
    zero_rna_row = sp.csr_matrix((1, rna_count_matrix.shape[1]), dtype=rna_count_matrix.dtype)
    rna_rows = [rna_count_matrix[rna_pos[cell_id]] if cell_id in rna_pos else zero_rna_row for cell_id in union_cell_ids]
    rna_aligned = sp.vstack(rna_rows, format="csr")

    circ_X_by_cell, circ_ids, circ_cell_ids = load_circ_matrix(
        matrix_path=circ_matrix_path,
        circ_index_path=circ_index_path,
        cell_index_path=cell_index_path,
    )
    circ_X_by_cell = _orient_circ_matrix_cells_by_circ(
        X=circ_X_by_cell,
        circ_ids=circ_ids,
        cell_ids=circ_cell_ids,
        matrix_path=circ_matrix_path,
    )
    circ_pos = {cell_id: idx for idx, cell_id in enumerate(circ_cell_ids)}
    zero_circ_row = sp.csr_matrix((1, circ_X_by_cell.shape[1]), dtype=circ_X_by_cell.dtype)
    circ_rows = [circ_X_by_cell[circ_pos[cell_id]] if cell_id in circ_pos else zero_circ_row for cell_id in union_cell_ids]
    circ_aligned = sp.vstack(circ_rows, format="csr")

    if circ_feature_path.exists():
        circ_var = pd.read_csv(circ_feature_path, sep="\t", keep_default_na=False)
        circ_id_col = "circ_id" if "circ_id" in circ_var.columns else circ_var.columns[0]
        circ_var = circ_var.set_index(circ_id_col, drop=False).reindex(circ_ids)
    else:
        circ_var = pd.DataFrame({"circ_id": circ_ids}, index=circ_ids)
    if circ_aligned.shape[1] != circ_var.shape[0]:
        raise ValueError(
            f"Circ AnnData axis mismatch after orientation: X has {circ_aligned.shape[1]} circRNA columns, "
            f"but circ var has {circ_var.shape[0]} rows. "
            f"matrix={circ_matrix_path}, circ_index={circ_index_path}, circ_feature_table={circ_feature_path if circ_feature_path.exists() else 'missing'}"
        )

    shared_obs = joined_df.copy()
    shared_obs.index = union_cell_ids
    shared_obs.index.name = None
    shared_obs_clean = sanitize_frame_for_anndata(shared_obs.reset_index(names="cell_id").set_index("cell_id"))

    rna_adata = ad.AnnData(
        X=rna_aligned.tocsr().astype(np.int32),
        obs=shared_obs_clean.copy(),
        var=sanitize_frame_for_anndata(rna_var),
    )
    circ_adata = ad.AnnData(
        X=circ_aligned.tocsr().astype(np.int32),
        obs=shared_obs_clean.copy(),
        var=sanitize_frame_for_anndata(circ_var),
    )
    mdata = mu.MuData({"rna": rna_adata, "circ": circ_adata})
    mdata.obs = shared_obs_clean.copy()
    workflow_summary_payload = load_json(workdir / "workflow_summary.json") if (workdir / "workflow_summary.json").exists() else {}
    rna_import_summary_payload = load_json(rna_dir / "rna_import_summary.json") if (rna_dir / "rna_import_summary.json").exists() else {}
    rna_circ_summary_payload = load_json(workdir / "qc" / "rna_circ_summary.json") if (workdir / "qc" / "rna_circ_summary.json").exists() else {}

    provenance = {
        "command_name": "circyto export-mudata",
        "circyto_version": __version__,
        "source_workdir": str(workdir.resolve()),
        "workflow_uuid": str(workflow_summary_payload.get("workflow_uuid", "")),
        "workflow_type": str(workflow_summary_payload.get("workflow_type", workflow_summary_payload.get("workflow", ""))),
        "protocol": str(workflow_summary_payload.get("protocol", "")),
        "read_layout": str(workflow_summary_payload.get("read_layout", "")),
        "n_obs": int(len(union_cell_ids)),
        "n_rna_cells": int(rna_aligned.shape[0]),
        "n_rna_features": int(rna_aligned.shape[1]),
        "n_circ_cells": int(circ_aligned.shape[0]),
        "n_circ_features": int(circ_aligned.shape[1]),
        "n_shared_cells": int(rna_circ_summary.get("n_shared_cells", 0)),
        "n_rna_only_cells": int(rna_circ_summary.get("n_rna_only_cells", 0)),
        "n_circ_only_cells": int(rna_circ_summary.get("n_circ_only_cells", 0)),
        "environment": {
            **get_environment_summary(),
            "export_timestamp": utc_now_iso(),
        },
        "workflow_summary_json": json.dumps(workflow_summary_payload, sort_keys=True, default=str),
        "rna_import_summary_json": json.dumps(rna_import_summary_payload, sort_keys=True, default=str),
        "rna_circ_summary_json": json.dumps(rna_circ_summary_payload, sort_keys=True, default=str),
    }
    mdata.uns["circyto"] = sanitize_for_h5ad_uns(provenance)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    mdata.write_h5mu(str(out_path))
    return {
        "path": str(out_path.resolve()),
        "n_obs": int(mdata.n_obs),
        "n_rna_cells": int(rna_aligned.shape[0]),
        "n_rna_features": int(rna_aligned.shape[1]),
        "n_circ_cells": int(circ_aligned.shape[0]),
        "n_circ_features": int(circ_aligned.shape[1]),
        "n_shared_cells": int(rna_circ_summary["n_shared_cells"]),
        "n_rna_only_cells": int(rna_circ_summary["n_rna_only_cells"]),
        "n_circ_only_cells": int(rna_circ_summary["n_circ_only_cells"]),
        "rna_only_cells": list(rna_circ_summary["rna_only_cells"]),
        "circ_only_cells": list(rna_circ_summary["circ_only_cells"]),
        "rna_only_cell_handling": "Included in shared mdata.obs and RNA modality; circ modality is zero-filled for missing cells.",
    }


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
