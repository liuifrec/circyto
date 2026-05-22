from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

import pandas as pd

from circyto.pipeline.workflow_reporting import write_json


SCOMATIC_CIRCYTO_COLUMNS = [
    "variant_id",
    "cell_id",
    "chrom",
    "pos",
    "ref",
    "alt",
    "gene",
    "filter_status",
    "candidate_variant_class",
    "read_support",
    "vaf",
    "caller",
]


def _find_column(df: pd.DataFrame, candidates: tuple[str, ...]) -> Optional[str]:
    lowered = {str(column).strip().lower(): str(column) for column in df.columns}
    for candidate in candidates:
        match = lowered.get(candidate.lower())
        if match is not None:
            return match
    return None


def normalize_scomatic_candidate_table(raw_df: pd.DataFrame) -> pd.DataFrame:
    cell_col = _find_column(raw_df, ("cell_id", "cell", "barcode", "sample", "sample_id"))
    chrom_col = _find_column(raw_df, ("chrom", "chr", "chromosome"))
    pos_col = _find_column(raw_df, ("pos", "position", "start"))
    ref_col = _find_column(raw_df, ("ref", "reference", "ref_allele"))
    alt_col = _find_column(raw_df, ("alt", "alternative", "alt_allele"))
    if cell_col is None or chrom_col is None or pos_col is None or ref_col is None or alt_col is None:
        raise ValueError(
            "Raw SComatic candidate table must contain cell, chrom, pos, ref, and alt columns."
        )

    gene_col = _find_column(raw_df, ("gene", "gene_name", "host_gene"))
    filter_col = _find_column(raw_df, ("filter_status", "filter", "FILTER"))
    class_col = _find_column(raw_df, ("candidate_variant_class", "variant_class", "class"))
    read_support_col = _find_column(raw_df, ("read_support", "support", "alt_reads", "alt_count"))
    vaf_col = _find_column(raw_df, ("vaf", "VAF", "allele_fraction"))
    caller_col = _find_column(raw_df, ("caller",))
    variant_id_col = _find_column(raw_df, ("variant_id", "id"))

    normalized = pd.DataFrame(
        {
            "cell_id": raw_df[cell_col].astype(str),
            "chrom": raw_df[chrom_col].astype(str),
            "pos": pd.to_numeric(raw_df[pos_col], errors="raise").astype(int),
            "ref": raw_df[ref_col].astype(str),
            "alt": raw_df[alt_col].astype(str),
            "gene": raw_df[gene_col].astype(str) if gene_col is not None else "",
            "filter_status": raw_df[filter_col].astype(str) if filter_col is not None else "PASS",
            "candidate_variant_class": raw_df[class_col].astype(str) if class_col is not None else "RNA-derived candidate variant signal",
            "read_support": pd.to_numeric(raw_df[read_support_col], errors="coerce").fillna(0).astype(int) if read_support_col is not None else 0,
            "vaf": pd.to_numeric(raw_df[vaf_col], errors="coerce").fillna(0.0).astype(float) if vaf_col is not None else 0.0,
            "caller": raw_df[caller_col].astype(str) if caller_col is not None else "SComatic",
        }
    )
    if variant_id_col is not None:
        normalized["variant_id"] = raw_df[variant_id_col].astype(str)
    else:
        normalized["variant_id"] = (
            normalized["chrom"]
            + ":"
            + normalized["pos"].astype(str)
            + ":"
            + normalized["ref"]
            + ">"
            + normalized["alt"]
            + ":"
            + normalized["cell_id"]
        )
    return normalized[SCOMATIC_CIRCYTO_COLUMNS]


def _load_workdir_cells(workdir: Path) -> list[str]:
    gene_counts = workdir / "rna" / "gene_counts.tsv"
    if gene_counts.exists():
        df = pd.read_csv(gene_counts, sep="\t", keep_default_na=False)
        return [str(column) for column in df.columns if column not in {"gene_id", "gene_name"}]
    cell_index = workdir / "matrix" / "cell_index.txt"
    if cell_index.exists():
        return [line.strip() for line in cell_index.read_text(encoding="utf-8").splitlines() if line.strip()]
    return []


def write_synthetic_scomatic_poc(
    *,
    workdir: Path,
    outdir: Path,
    reference: Path,
    gtf: Path,
) -> dict[str, object]:
    outdir.mkdir(parents=True, exist_ok=True)
    cell_ids = _load_workdir_cells(workdir)
    if not cell_ids:
        raise ValueError(f"Could not infer cell IDs from {workdir}; expected rna/gene_counts.tsv or matrix/cell_index.txt")
    first_cell = cell_ids[0]
    raw_df = pd.DataFrame(
        [
            {
                "cell_id": first_cell,
                "chrom": "chr21",
                "pos": 1000,
                "ref": "A",
                "alt": "G",
                "gene": "TEST_GENE1",
                "filter_status": "PASS",
                "candidate_variant_class": "RNA-derived candidate variant signal",
                "read_support": 3,
                "vaf": 0.25,
                "caller": "SComatic-synthetic",
            }
        ]
    )
    normalized = normalize_scomatic_candidate_table(raw_df)
    out_tsv = outdir / "scomatic_candidate_summary.tsv"
    normalized.to_csv(out_tsv, sep="\t", index=False)
    summary = {
        "mode": "synthetic",
        "workdir": str(workdir.resolve()),
        "reference": str(reference.resolve()),
        "gtf": str(gtf.resolve()),
        "n_candidates": int(normalized.shape[0]),
        "cell_ids": cell_ids,
        "output_scomatic_candidate_summary": str(out_tsv.resolve()),
        "terminology_note": "Outputs are RNA-derived candidate variant signals, not validated somatic mutations.",
    }
    write_json(outdir / "scomatic_poc_summary.json", summary)
    return summary
