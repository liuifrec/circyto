from __future__ import annotations

import warnings
from pathlib import Path
from typing import Optional

import pandas as pd


CELL_ID_CANDIDATES = (
    "cell_id",
    "cell",
    "barcode",
    "sample",
    "sample_id",
)

BAM_CANDIDATES = (
    "bam",
    "bam_path",
    "alignment_bam",
    "bam_file",
)

CELLTYPE_CANDIDATES = (
    "cell_type",
    "celltype",
    "cell_label",
    "cluster",
    "annotation",
    "cell_type_or_system",
)

GENE_CANDIDATES = (
    "host_gene",
    "gene",
    "gene_name",
    "gene_symbol",
    "gene_id",
    "Gene",
    "GENE",
)


def _read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", keep_default_na=False)


def _find_column(df: pd.DataFrame, candidates: tuple[str, ...]) -> Optional[str]:
    lowered = {str(column).strip().lower(): str(column) for column in df.columns}
    for candidate in candidates:
        match = lowered.get(candidate.lower())
        if match is not None:
            return match
    return None


def _require_column(df: pd.DataFrame, candidates: tuple[str, ...], label: str, source: Path) -> str:
    column = _find_column(df, candidates)
    if column is None:
        expected = ", ".join(candidates)
        raise ValueError(f"{source} is missing a {label} column. Tried: {expected}")
    return column


def _coerce_numeric(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def _load_circ_matrix_summary(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    df = _read_tsv(path)
    if df.empty:
        return (
            pd.DataFrame(columns=["cell_id", "circRNA_count", "circRNA_read_support_sum"]),
            pd.DataFrame(columns=["circ_id", "cell_count", "support_sum"]),
        )

    circ_id_col = _find_column(df, ("circ_id", "circRNA_id", "feature_id"))
    if circ_id_col is None:
        circ_id_col = str(df.columns[0])
        warnings.warn(
            f"{path} has no explicit circ_id column; using first column '{circ_id_col}'.",
            stacklevel=2,
        )

    cell_columns = [column for column in df.columns if column != circ_id_col]
    if not cell_columns:
        warnings.warn(f"{path} has no cell columns; writing zero-valued circ summaries.", stacklevel=2)
        return (
            pd.DataFrame(columns=["cell_id", "circRNA_count", "circRNA_read_support_sum"]),
            pd.DataFrame(
                {
                    "circ_id": df[circ_id_col].astype(str),
                    "cell_count": 0,
                    "support_sum": 0.0,
                }
            ),
        )

    numeric = df[cell_columns].apply(lambda column: column.map(_coerce_numeric))
    per_cell = pd.DataFrame(
        {
            "cell_id": [str(column) for column in cell_columns],
            "circRNA_count": (numeric > 0).sum(axis=0).astype(int).tolist(),
            "circRNA_read_support_sum": numeric.sum(axis=0).tolist(),
        }
    )
    per_circ = pd.DataFrame(
        {
            "circ_id": df[circ_id_col].astype(str),
            "cell_count": (numeric > 0).sum(axis=1).astype(int),
            "support_sum": numeric.sum(axis=1),
        }
    )
    return per_cell, per_circ


def export_scomatic_inputs(
    *,
    bam_manifest: Path,
    cell_metadata: Path,
    outdir: Path,
    reference_fasta: Path,
    protocol: str,
) -> None:
    bam_df = _read_tsv(bam_manifest)
    meta_df = _read_tsv(cell_metadata)

    bam_cell_col = _require_column(bam_df, CELL_ID_CANDIDATES, "cell identifier", bam_manifest)
    bam_path_col = _require_column(bam_df, BAM_CANDIDATES, "BAM path", bam_manifest)
    meta_cell_col = _require_column(meta_df, CELL_ID_CANDIDATES, "cell identifier", cell_metadata)
    meta_celltype_col = _find_column(meta_df, CELLTYPE_CANDIDATES)

    if meta_celltype_col is None:
        warnings.warn(
            f"{cell_metadata} has no cell-type-like column; defaulting all exported labels to 'unknown'.",
            stacklevel=2,
        )
        meta_export = pd.DataFrame(
            {
                "cell_id": meta_df[meta_cell_col].astype(str),
                "cell_type": "unknown",
            }
        )
    else:
        meta_export = pd.DataFrame(
            {
                "cell_id": meta_df[meta_cell_col].astype(str),
                "cell_type": meta_df[meta_celltype_col].astype(str).replace("", "unknown"),
            }
        )

    bam_export = pd.DataFrame(
        {
            "cell_id": bam_df[bam_cell_col].astype(str),
            "bam": bam_df[bam_path_col].astype(str),
            "protocol": protocol,
            "reference_fasta": str(reference_fasta),
        }
    )

    missing_metadata = sorted(set(bam_export["cell_id"]) - set(meta_export["cell_id"]))
    if missing_metadata:
        warnings.warn(
            f"{len(missing_metadata)} BAM entries have no matching cell metadata; example: {missing_metadata[:3]}",
            stacklevel=2,
        )

    outdir.mkdir(parents=True, exist_ok=True)
    bam_export.to_csv(outdir / "scomatic_bam_list.tsv", sep="\t", index=False)
    meta_export.drop_duplicates(subset=["cell_id"]).to_csv(outdir / "scomatic_celltypes.tsv", sep="\t", index=False)

    readme = (
        "# SComatic interoperability scaffold\n\n"
        "This export is intentionally limited to metadata and input tables.\n\n"
        "Important scientific framing:\n"
        "RNA-derived SComatic calls should be treated as exploratory candidate SNVs.\n"
        "They are not validated DNA somatic mutations unless orthogonal DNA validation exists.\n\n"
        "Generated files:\n"
        "- `scomatic_bam_list.tsv`: per-cell BAM list prepared from the circyto manifest.\n"
        "- `scomatic_celltypes.tsv`: per-cell labels for downstream grouping.\n\n"
        "Next steps outside circyto:\n"
        "1. Run SComatic externally using these tables plus the same reference FASTA.\n"
        "2. Review RNA-specific artifacts and filtering assumptions carefully.\n"
        "3. Re-import the exploratory candidate SNV table with `circyto join-circ-snv-summary`.\n\n"
        f"Protocol hint: `{protocol}`\n"
        f"Reference FASTA: `{reference_fasta}`\n"
    )
    (outdir / "README_scomatic_next_steps.md").write_text(readme, encoding="utf-8")


def join_circ_snv_summary(
    *,
    circ_matrix: Path,
    circ_feature_table: Path,
    scomatic_candidates: Path,
    cell_metadata: Path,
    outdir: Path,
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    cell_meta_df = _read_tsv(cell_metadata)
    cell_meta_col = _require_column(cell_meta_df, CELL_ID_CANDIDATES, "cell identifier", cell_metadata)
    cell_meta_export = cell_meta_df.rename(columns={cell_meta_col: "cell_id"}).copy()
    cell_meta_export["cell_id"] = cell_meta_export["cell_id"].astype(str)

    per_cell_circ, per_circ = _load_circ_matrix_summary(circ_matrix)

    snv_df = _read_tsv(scomatic_candidates)
    snv_cell_col = _require_column(snv_df, CELL_ID_CANDIDATES, "cell identifier", scomatic_candidates)
    snv_gene_col = _find_column(snv_df, GENE_CANDIDATES)
    if snv_gene_col is None:
        warnings.warn(
            f"{scomatic_candidates} has no gene-like column; host-gene candidate SNV summary will contain zeros.",
            stacklevel=2,
        )

    per_cell_snv = (
        snv_df.assign(cell_id=snv_df[snv_cell_col].astype(str))
        .groupby("cell_id", as_index=False)
        .size()
        .rename(columns={"size": "candidate_snv_count"})
    )

    cell_summary = cell_meta_export.merge(per_cell_circ, on="cell_id", how="left")
    cell_summary = cell_summary.merge(per_cell_snv, on="cell_id", how="left")
    if "circRNA_count" not in cell_summary.columns:
        cell_summary["circRNA_count"] = 0
    if "circRNA_read_support_sum" not in cell_summary.columns:
        cell_summary["circRNA_read_support_sum"] = 0.0
    cell_summary["candidate_snv_count"] = cell_summary["candidate_snv_count"].fillna(0).astype(int)
    cell_summary["circRNA_count"] = cell_summary["circRNA_count"].fillna(0).astype(int)
    cell_summary["circRNA_read_support_sum"] = cell_summary["circRNA_read_support_sum"].fillna(0.0)
    cell_summary.to_csv(outdir / "circ_snv_cell_summary.tsv", sep="\t", index=False)

    feature_df = _read_tsv(circ_feature_table)
    circ_id_col = _find_column(feature_df, ("circ_id", "circRNA_id", "feature_id"))
    host_gene_col = _find_column(feature_df, ("host_gene", "gene", "gene_name", "gene_symbol"))

    circ_gene_summary = pd.DataFrame(columns=["gene", "circRNA_feature_count", "circRNA_cell_count", "circRNA_read_support_sum"])
    if circ_id_col is None:
        warnings.warn(
            f"{circ_feature_table} has no circ_id column; circ host-gene summary will omit circ-level aggregation.",
            stacklevel=2,
        )
    elif host_gene_col is None:
        warnings.warn(
            f"{circ_feature_table} has no host_gene-like column; circ host-gene summary will omit circ gene aggregation.",
            stacklevel=2,
        )
    else:
        merged = feature_df[[circ_id_col, host_gene_col]].copy()
        merged[circ_id_col] = merged[circ_id_col].astype(str)
        merged[host_gene_col] = merged[host_gene_col].astype(str)
        merged = merged[merged[host_gene_col] != ""]
        merged = merged.merge(per_circ, left_on=circ_id_col, right_on="circ_id", how="left")
        merged["cell_count"] = merged["cell_count"].fillna(0).astype(int)
        merged["support_sum"] = merged["support_sum"].fillna(0.0)
        circ_gene_summary = (
            merged.groupby(host_gene_col, as_index=False)
            .agg(
                circRNA_feature_count=(circ_id_col, "nunique"),
                circRNA_cell_count=("cell_count", "sum"),
                circRNA_read_support_sum=("support_sum", "sum"),
            )
            .rename(columns={host_gene_col: "gene"})
        )

    snv_gene_summary = pd.DataFrame(columns=["gene", "candidate_snv_count"])
    if snv_gene_col is not None:
        snv_gene_summary = (
            snv_df.assign(gene=snv_df[snv_gene_col].astype(str))
            .loc[lambda frame: frame["gene"] != "", ["gene"]]
            .groupby("gene", as_index=False)
            .size()
            .rename(columns={"size": "candidate_snv_count"})
        )

    host_gene_summary = circ_gene_summary.merge(snv_gene_summary, on="gene", how="outer")
    if host_gene_summary.empty:
        host_gene_summary = pd.DataFrame(
            columns=[
                "gene",
                "circRNA_feature_count",
                "circRNA_cell_count",
                "circRNA_read_support_sum",
                "candidate_snv_count",
            ]
        )
    else:
        for column in ("circRNA_feature_count", "circRNA_cell_count", "candidate_snv_count"):
            if column not in host_gene_summary.columns:
                host_gene_summary[column] = 0
            host_gene_summary[column] = host_gene_summary[column].fillna(0).astype(int)
        if "circRNA_read_support_sum" not in host_gene_summary.columns:
            host_gene_summary["circRNA_read_support_sum"] = 0.0
        host_gene_summary["circRNA_read_support_sum"] = host_gene_summary["circRNA_read_support_sum"].fillna(0.0)
        host_gene_summary = host_gene_summary.sort_values("gene").reset_index(drop=True)

    host_gene_summary.to_csv(outdir / "circ_snv_host_gene_summary.tsv", sep="\t", index=False)
