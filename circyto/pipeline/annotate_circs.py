from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional
import json

import pandas as pd

from circyto.pipeline.annotate_host_gene import normalize_host_gene_annotations

try:
    import anndata as ad

    HAS_ANNDATA = True
except Exception:
    HAS_ANNDATA = False


MATCH_PRIORITY = {
    "known_exact_bsj": 3,
    "known_near_bsj": 2,
    "known_host_gene_only": 1,
    "novel": 0,
}


@dataclass(frozen=True)
class AnnotationDBSpec:
    name: str
    path: Path
    chrom_column: str
    start_column: str
    end_column: str
    strand_column: Optional[str] = None
    id_column: Optional[str] = None
    host_gene_column: Optional[str] = None
    extra_columns: tuple[str, ...] = ()


def parse_annotation_db_spec(spec_text: str) -> AnnotationDBSpec:
    spec_text = spec_text.strip()
    if not spec_text:
        raise ValueError("Empty annotation DB specification")
    if "=" not in spec_text:
        parts = spec_text.split(":", 1)
        if len(parts) != 2 or not parts[0].strip() or not parts[1].strip():
            raise ValueError(
                "Shorthand annotation DB spec must look like 'name:path/to/db.tsv'"
            )
        return AnnotationDBSpec(
            name=parts[0].strip(),
            path=Path(parts[1].strip()),
            chrom_column="chrom",
            start_column="start",
            end_column="end",
            strand_column="strand",
            id_column="id",
            host_gene_column="host_gene",
            extra_columns=(),
        )

    fields: dict[str, str] = {}
    for item in spec_text.split(";"):
        item = item.strip()
        if not item:
            continue
        key, sep, value = item.partition("=")
        if not sep:
            raise ValueError(f"Invalid annotation DB field: {item}")
        fields[key.strip().lower()] = value.strip()

    for key in ("name", "path", "chrom", "start", "end"):
        if not fields.get(key):
            raise ValueError(f"Missing required annotation DB field '{key}'")

    extras = tuple(
        column.strip()
        for column in fields.get("extra", "").split(",")
        if column.strip()
    )
    return AnnotationDBSpec(
        name=fields["name"],
        path=Path(fields["path"]),
        chrom_column=fields["chrom"],
        start_column=fields["start"],
        end_column=fields["end"],
        strand_column=fields.get("strand") or None,
        id_column=fields.get("id") or None,
        host_gene_column=fields.get("host_gene") or None,
        extra_columns=extras,
    )


def _load_circ_table(circ_table: Path) -> pd.DataFrame:
    df = pd.read_csv(circ_table, sep="\t", keep_default_na=False)
    for column in ("circ_id", "chrom", "start", "end"):
        if column not in df.columns:
            raise ValueError(f"circ table is missing required column '{column}'")
    for column in ("strand", "host_gene"):
        if column not in df.columns:
            df[column] = ""
    df["circ_id"] = df["circ_id"].astype(str)
    df["chrom"] = df["chrom"].astype(str)
    df["strand"] = df["strand"].astype(str)
    df["host_gene"] = df["host_gene"].astype(str)
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    if df["start"].isna().any() or df["end"].isna().any():
        raise ValueError("circ table contains non-numeric start/end coordinates")
    return df


def _normalize_db_table(spec: AnnotationDBSpec) -> pd.DataFrame:
    if not spec.path.exists():
        raise FileNotFoundError(f"Annotation DB file not found: {spec.path}")
    raw = pd.read_csv(spec.path, sep="\t", keep_default_na=False)
    required = [spec.chrom_column, spec.start_column, spec.end_column]
    for column in required:
        if column not in raw.columns:
            raise ValueError(
                f"Annotation DB '{spec.name}' is missing required column '{column}'"
            )
    for column in spec.extra_columns:
        if column not in raw.columns:
            raise ValueError(
                f"Annotation DB '{spec.name}' is missing extra column '{column}'"
            )

    start = pd.to_numeric(raw[spec.start_column], errors="coerce")
    end = pd.to_numeric(raw[spec.end_column], errors="coerce")
    valid_mask = start.notna() & end.notna()
    host_gene_column = (
        spec.host_gene_column
        if spec.host_gene_column and spec.host_gene_column in raw.columns
        else None
    )
    if host_gene_column is None:
        for candidate in (
            "host_gene",
            "host_genes",
            "gene_symbol",
            "gene_name",
            "host_gene_symbol",
            "parent_gene",
            "gene",
            "symbol",
        ):
            if candidate in raw.columns:
                host_gene_column = candidate
                break

    db = pd.DataFrame(
        {
            "_db_chrom": raw[spec.chrom_column].astype(str),
            "_db_start": start,
            "_db_end": end,
            "_db_strand": raw[spec.strand_column].astype(str)
            if spec.strand_column and spec.strand_column in raw.columns
            else "",
            "_db_id": raw[spec.id_column].astype(str)
            if spec.id_column and spec.id_column in raw.columns
            else "",
            "_db_host_gene": raw[host_gene_column].astype(str)
            if host_gene_column
            else "",
        }
    )
    if spec.id_column is None or spec.id_column not in raw.columns:
        db["_db_id"] = (
            db["_db_chrom"]
            + ":"
            + db["_db_start"].astype("Int64").astype(str)
            + "|"
            + db["_db_end"].astype("Int64").astype(str)
        )
    for column in spec.extra_columns:
        db[column] = raw[column].astype(str)
    db = db.loc[valid_mask].reset_index(drop=True)
    return db


def _join_unique(values: list[str]) -> str:
    seen: set[str] = set()
    ordered: list[str] = []
    for value in values:
        if not value or value in seen:
            continue
        seen.add(value)
        ordered.append(value)
    return "|".join(ordered)


def _clean_tsv_cell_text(value: Any) -> str:
    text = str(value)
    if text == "":
        return ""
    return text.replace("\t", " ").replace("\r", " ").replace("\n", " ")


def _build_extra_summary(matches: pd.DataFrame, extra_columns: tuple[str, ...]) -> str:
    if matches.empty or not extra_columns:
        return ""
    parts: list[str] = []
    for _, row in matches.iterrows():
        values = [
            f"{column}={_clean_tsv_cell_text(row[column])}"
            for column in extra_columns
            if _clean_tsv_cell_text(row[column]) != ""
        ]
        if values:
            parts.append(",".join(values))
    return _join_unique(parts)


def _strand_mask(circ_row: pd.Series, db: pd.DataFrame) -> pd.Series:
    circ_strand = str(circ_row.get("strand", ""))
    if circ_strand == "" or "_db_strand" not in db.columns:
        return pd.Series(True, index=db.index)
    return (db["_db_strand"] == "") | (db["_db_strand"] == circ_strand)


def _match_database(
    circ_row: pd.Series,
    db: pd.DataFrame,
    *,
    max_bsj_distance: int,
    enable_host_gene_match: bool,
) -> tuple[str, pd.DataFrame]:
    candidates = db.loc[db["_db_chrom"] == str(circ_row["chrom"])].copy()
    if candidates.empty:
        return "novel", candidates

    strand_ok = _strand_mask(circ_row, candidates)
    exact = candidates.loc[
        strand_ok
        & (candidates["_db_start"] == circ_row["start"])
        & (candidates["_db_end"] == circ_row["end"])
    ]
    if not exact.empty:
        return "known_exact_bsj", exact

    if max_bsj_distance > 0:
        near = candidates.loc[
            strand_ok
            & (candidates["_db_start"].sub(circ_row["start"]).abs() <= max_bsj_distance)
            & (candidates["_db_end"].sub(circ_row["end"]).abs() <= max_bsj_distance)
        ]
        if not near.empty:
            return "known_near_bsj", near

    circ_host_gene = str(circ_row.get("host_gene", ""))
    if enable_host_gene_match and circ_host_gene:
        overlap = (
            (candidates["_db_start"] <= circ_row["end"])
            & (candidates["_db_end"] >= circ_row["start"])
        )
        host_gene = candidates["_db_host_gene"] == circ_host_gene
        host_only = candidates.loc[strand_ok & overlap & host_gene]
        if not host_only.empty:
            return "known_host_gene_only", host_only

    return "novel", candidates.iloc[0:0]


def annotate_circ_table(
    *,
    circ_table: Path,
    db_specs: list[AnnotationDBSpec],
    out_path: Path,
    summary_path: Optional[Path] = None,
    max_bsj_distance: int = 0,
    enable_host_gene_match: bool = False,
    update_h5ad: Optional[Path] = None,
) -> dict[str, Any]:
    circ_df = _load_circ_table(circ_table)
    db_tables = {spec.name: _normalize_db_table(spec) for spec in db_specs}

    annotated = circ_df.copy()
    for spec in db_specs:
        annotated[f"{spec.name}_known"] = False
        annotated[f"{spec.name}_match_type"] = "novel"
        annotated[f"{spec.name}_ids"] = ""
        annotated[f"{spec.name}_host_genes"] = ""
        annotated[f"{spec.name}_extra_summary"] = ""

    known_counts: list[int] = []
    known_dbs: list[str] = []
    best_statuses: list[str] = []

    for _, row in circ_df.iterrows():
        matched_dbs: list[str] = []
        best_status = "novel"
        values: dict[str, Any] = {}
        for spec in db_specs:
            match_type, matches = _match_database(
                row,
                db_tables[spec.name],
                max_bsj_distance=max_bsj_distance,
                enable_host_gene_match=enable_host_gene_match,
            )
            known = match_type != "novel"
            if known:
                matched_dbs.append(spec.name)
                if MATCH_PRIORITY[match_type] > MATCH_PRIORITY[best_status]:
                    best_status = match_type
            values[f"{spec.name}_known"] = bool(known)
            values[f"{spec.name}_match_type"] = match_type
            values[f"{spec.name}_ids"] = _join_unique(
                [_clean_tsv_cell_text(item) for item in matches["_db_id"].astype(str).tolist()]
            )
            values[f"{spec.name}_host_genes"] = _join_unique(
                [
                    _clean_tsv_cell_text(item)
                    for item in matches["_db_host_gene"].astype(str).tolist()
                ]
            )
            values[f"{spec.name}_extra_summary"] = _build_extra_summary(
                matches, spec.extra_columns
            )
        known_counts.append(len(matched_dbs))
        known_dbs.append("|".join(matched_dbs))
        best_statuses.append(best_status)
        row_mask = annotated["circ_id"] == str(row["circ_id"])
        for column, value in values.items():
            annotated.loc[row_mask, column] = value

    annotated["known_database_count"] = known_counts
    annotated["known_databases"] = known_dbs
    annotated["best_annotation_status"] = best_statuses
    annotated = normalize_host_gene_annotations(
        annotated,
        legacy_host_gene_source="gtf",
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    annotated.to_csv(out_path, sep="\t", index=False, lineterminator="\n")

    if summary_path is None:
        summary_path = out_path.with_name("annotation_summary.json")
    summary = build_annotation_summary(annotated, db_specs)
    summary["circ_table"] = str(circ_table.resolve())
    summary["annotated_circ_table"] = str(out_path.resolve())
    summary["max_bsj_distance"] = int(max_bsj_distance)
    summary["enable_host_gene_match"] = bool(enable_host_gene_match)
    summary["databases"] = [
        {
            "name": spec.name,
            "path": str(spec.path.resolve()),
            "chrom_column": spec.chrom_column,
            "start_column": spec.start_column,
            "end_column": spec.end_column,
            "strand_column": spec.strand_column,
            "id_column": spec.id_column,
            "host_gene_column": spec.host_gene_column,
            "extra_columns": list(spec.extra_columns),
        }
        for spec in db_specs
    ]
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    if update_h5ad is not None:
        update_h5ad_annotations(
            h5ad_path=update_h5ad,
            annotated_df=annotated,
            db_specs=db_specs,
            summary_path=summary_path,
            source_circ_table=circ_table,
        )

    return summary


def build_annotation_summary(
    annotated_df: pd.DataFrame,
    db_specs: list[AnnotationDBSpec],
) -> dict[str, Any]:
    payload: dict[str, Any] = {
        "total_circRNAs": int(len(annotated_df)),
        "novel_count": int((annotated_df["best_annotation_status"] == "novel").sum()),
        "exact_matches_per_database": {},
        "near_matches_per_database": {},
        "host_gene_only_matches_per_database": {},
    }
    for spec in db_specs:
        match_col = f"{spec.name}_match_type"
        payload["exact_matches_per_database"][spec.name] = int(
            (annotated_df[match_col] == "known_exact_bsj").sum()
        )
        payload["near_matches_per_database"][spec.name] = int(
            (annotated_df[match_col] == "known_near_bsj").sum()
        )
        payload["host_gene_only_matches_per_database"][spec.name] = int(
            (annotated_df[match_col] == "known_host_gene_only").sum()
        )

    if "n_cells_detected" in annotated_df.columns:
        recurrent = annotated_df.loc[
            pd.to_numeric(annotated_df["n_cells_detected"], errors="coerce").fillna(0) > 1
        ].copy()
        recurrent["n_cells_detected"] = pd.to_numeric(
            recurrent["n_cells_detected"], errors="coerce"
        ).fillna(0).astype(int)
        if "total_support" in recurrent.columns:
            recurrent["total_support"] = pd.to_numeric(
                recurrent["total_support"], errors="coerce"
            ).fillna(0).astype(int)
        else:
            recurrent["total_support"] = 0
        payload["recurrent_known_count"] = int(
            (recurrent["best_annotation_status"] != "novel").sum()
        )
        payload["recurrent_novel_count"] = int(
            (recurrent["best_annotation_status"] == "novel").sum()
        )
        novel_recurrent = recurrent.loc[
            recurrent["best_annotation_status"] == "novel"
        ].sort_values(
            by=["n_cells_detected", "total_support", "circ_id"],
            ascending=[False, False, True],
        )
        top_columns = [
            column
            for column in (
                "circ_id",
                "chrom",
                "start",
                "end",
                "strand",
                "host_gene",
                "n_cells_detected",
                "total_support",
            )
            if column in novel_recurrent.columns
        ]
        payload["top_novel_recurrent_circRNAs"] = novel_recurrent[top_columns].head(20).to_dict("records")
    return payload


def update_h5ad_annotations(
    *,
    h5ad_path: Path,
    annotated_df: pd.DataFrame,
    db_specs: list[AnnotationDBSpec],
    summary_path: Path,
    source_circ_table: Path,
) -> None:
    if not HAS_ANNDATA:
        raise RuntimeError("anndata is required for --update-h5ad")
    adata = ad.read_h5ad(str(h5ad_path))
    if "circ_id" not in annotated_df.columns:
        raise ValueError("annotated circ table is missing 'circ_id'")
    annotation = annotated_df.set_index("circ_id").reindex(adata.var_names)
    if annotation.index.tolist() != list(adata.var_names):
        raise AssertionError("annotated circ table could not be aligned to h5ad var_names")

    for column in annotation.columns:
        adata.var[column] = annotation[column].tolist()

    circyto_uns: dict[str, Any]
    if isinstance(adata.uns.get("circyto"), dict):
        circyto_uns = dict(adata.uns["circyto"])
    else:
        circyto_uns = {}
    circyto_uns["circ_annotation"] = {
        "source_circ_table": str(source_circ_table.resolve()),
        "summary_json": str(summary_path.resolve()),
        "databases": [spec.name for spec in db_specs],
    }
    adata.uns["circyto"] = circyto_uns
    adata.write_h5ad(str(h5ad_path))
