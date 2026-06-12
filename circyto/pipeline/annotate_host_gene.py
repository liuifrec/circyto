from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple
import re

import pandas as pd


HOST_GENE_COLUMNS = [
    "host_gene",
    "host_gene_source",
    "host_gene_from_gtf",
    "host_gene_from_circatlas",
    "host_gene_from_circatlas_id",
]
HOST_GENE_SOURCES = {"gtf", "circatlas", "circatlas_id"}
CIRCATLAS_ID_HOST_RE = re.compile(r"hsa-([^_]+)_")
EMPTY_TEXT_VALUES = {"", "nan", "none", "na", "n/a", "<na>"}


def _clean_annotation_text(value: Any) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    text = str(value).strip()
    if text.lower() in EMPTY_TEXT_VALUES:
        return ""
    return text


def _split_annotation_values(value: Any) -> list[str]:
    text = _clean_annotation_text(value)
    if not text:
        return []
    return [
        item
        for item in (_clean_annotation_text(part) for part in re.split(r"[;|,]", text))
        if item
    ]


def _join_unique_semicolon(values: Iterable[Any]) -> str:
    seen: set[str] = set()
    ordered: list[str] = []
    for value in values:
        for item in _split_annotation_values(value):
            if item in seen:
                continue
            seen.add(item)
            ordered.append(item)
    return ";".join(ordered)


def _parse_gtf_attrs(attr_str: str) -> Dict[str, str]:
    """Parse a GTF/GFF attribute column into a dict."""
    attrs: Dict[str, str] = {}
    for field in attr_str.strip().split(";"):
        field = field.strip()
        if not field:
            continue
        if "=" in field:
            key, value = field.split("=", 1)
        elif " " in field:
            key, value = field.split(" ", 1)
        else:
            continue
        key = key.strip()
        value = value.strip().strip('"')
        if key:
            attrs[key] = value
    return attrs


def _attr_first(attrs: dict[str, str], keys: Iterable[str]) -> str:
    for key in keys:
        value = _clean_annotation_text(attrs.get(key, ""))
        if value:
            return value
    return ""


def load_genes_from_gtf(gtf_path: Path) -> pd.DataFrame:
    """Load gene intervals from a GTF/GFF file.

    Prefer explicit gene feature rows. If none exist, fall back to exon spans
    grouped by gene_id where possible.
    """
    gene_rows: list[dict[str, Any]] = []
    exon_spans: dict[str, dict[str, Any]] = {}

    with gtf_path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, _source, feature, start, end, _score, strand, _frame, attrs_text = parts
            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                continue

            attrs = _parse_gtf_attrs(attrs_text)
            gene_id = _attr_first(attrs, ("gene_id", "gene", "ID", "Parent"))
            gene_name = _attr_first(
                attrs,
                ("gene_name", "gene_symbol", "gene", "Name", "symbol", "gene_id", "ID"),
            )
            if not gene_name:
                gene_name = gene_id

            if feature == "gene":
                if not (gene_id or gene_name):
                    continue
                gene_rows.append(
                    {
                        "chrom": chrom,
                        "start": start_i,
                        "end": end_i,
                        "strand": strand,
                        "gene_id": gene_id or gene_name,
                        "gene_name": gene_name or gene_id,
                    }
                )
            elif feature == "exon" and (gene_id or gene_name):
                key = gene_id or gene_name
                existing = exon_spans.get(key)
                if existing is None:
                    exon_spans[key] = {
                        "chrom": chrom,
                        "start": start_i,
                        "end": end_i,
                        "strand": strand,
                        "gene_id": gene_id or gene_name,
                        "gene_name": gene_name or gene_id,
                    }
                else:
                    existing["start"] = min(int(existing["start"]), start_i)
                    existing["end"] = max(int(existing["end"]), end_i)

    rows = gene_rows if gene_rows else list(exon_spans.values())
    return pd.DataFrame(
        rows,
        columns=["chrom", "start", "end", "strand", "gene_id", "gene_name"],
    )


def _to_int_or_none(value: Any) -> Optional[int]:
    try:
        if pd.isna(value):
            return None
        return int(value)
    except Exception:
        return None


def _find_host_genes_for_circ(
    chrom: str,
    start: Optional[int],
    end: Optional[int],
    strand: str,
    genes: pd.DataFrame,
    max_genes: int = 5,
) -> Tuple[str, str, str, str, int]:
    """Return (host_gene, host_gene_id, multi_names, multi_ids, n) for one circ."""
    chrom = _clean_annotation_text(chrom)
    strand = _clean_annotation_text(strand)
    if not chrom or start is None or end is None or genes.empty:
        return "", "", "", "", 0

    g_chr = genes[genes["chrom"].astype(str) == chrom]
    if g_chr.empty:
        return "", "", "", "", 0

    mask_overlap = ~((end < g_chr["start"]) | (start > g_chr["end"]))
    g_ov = g_chr[mask_overlap]
    if g_ov.empty:
        return "", "", "", "", 0

    if strand in {"+", "-"}:
        g_same = g_ov[g_ov["strand"].astype(str) == strand]
        candidates = g_same if not g_same.empty else g_ov
    else:
        candidates = g_ov

    overlap_len = candidates[["start", "end"]].apply(
        lambda row: min(end, int(row["end"])) - max(start, int(row["start"])) + 1,
        axis=1,
    )
    gene_span = candidates["end"] - candidates["start"]
    order = (
        pd.DataFrame({"overlap": overlap_len, "span": gene_span})
        .assign(idx=candidates.index)
        .sort_values(["overlap", "span", "idx"], ascending=[False, True, True])["idx"]
    )

    names: list[str] = []
    ids: list[str] = []
    seen_pairs: set[tuple[str, str]] = set()
    for _, row in candidates.loc[order].iterrows():
        gene_id = _clean_annotation_text(row.get("gene_id", ""))
        gene_name = _clean_annotation_text(row.get("gene_name", "")) or gene_id
        if not gene_name:
            continue
        pair = (gene_name, gene_id)
        if pair in seen_pairs:
            continue
        seen_pairs.add(pair)
        names.append(gene_name)
        ids.append(gene_id)
        if len(names) >= max_genes:
            break

    if not names:
        return "", "", "", "", 0

    joined_names = ";".join(names)
    joined_ids = ";".join(gene_id for gene_id in ids if gene_id)
    return joined_names, joined_ids, joined_names, joined_ids, len(names)


def _apply_gtf_host_gene_columns(
    df: pd.DataFrame,
    *,
    gtf_path: Path,
    max_genes_per_circ: int,
) -> pd.DataFrame:
    out = df.copy()
    for column in ("chrom", "start", "end", "strand"):
        if column not in out.columns:
            out[column] = pd.NA if column in {"start", "end"} else ""

    genes = load_genes_from_gtf(gtf_path)

    host_gene: list[str] = []
    host_gene_id: list[str] = []
    host_genes_multi: list[str] = []
    host_gene_ids_multi: list[str] = []
    host_gene_n: list[int] = []

    for _, row in out.iterrows():
        h, h_id, h_multi, h_ids_multi, n = _find_host_genes_for_circ(
            _clean_annotation_text(row.get("chrom", "")),
            _to_int_or_none(row.get("start")),
            _to_int_or_none(row.get("end")),
            _clean_annotation_text(row.get("strand", "")),
            genes,
            max_genes=max_genes_per_circ,
        )
        host_gene.append(h)
        host_gene_id.append(h_id)
        host_genes_multi.append(h_multi)
        host_gene_ids_multi.append(h_ids_multi)
        host_gene_n.append(n)

    out["host_gene_from_gtf"] = host_gene
    out["host_gene_id"] = host_gene_id
    out["host_genes_multi"] = host_genes_multi
    out["host_gene_ids_multi"] = host_gene_ids_multi
    out["host_gene_n"] = host_gene_n
    return out


def _circatlas_host_gene_columns(columns: Iterable[str]) -> list[str]:
    selected: list[str] = []
    for column in columns:
        lower = column.lower()
        if lower in {"host_gene_from_circatlas", "host_gene_from_circatlas_id"}:
            continue
        if "circatlas" in lower and (
            lower.endswith("_host_genes") or lower.endswith("_host_gene")
        ):
            selected.append(column)
    return sorted(
        selected,
        key=lambda item: (
            item.lower() != "circatlas_v3_host_genes",
            item.lower(),
        ),
    )


def _circatlas_id_columns(columns: Iterable[str]) -> list[str]:
    selected: list[str] = []
    for column in columns:
        lower = column.lower()
        if lower == "host_gene_from_circatlas_id":
            continue
        if "circatlas" in lower and lower.endswith("_ids"):
            selected.append(column)
    return sorted(
        selected,
        key=lambda item: (
            item.lower() != "circatlas_v3_ids",
            item.lower(),
        ),
    )


def parse_circatlas_id_host_genes(value: Any) -> str:
    text = _clean_annotation_text(value)
    if not text:
        return ""
    return _join_unique_semicolon(CIRCATLAS_ID_HOST_RE.findall(text))


def normalize_host_gene_annotations(
    df: pd.DataFrame,
    *,
    gtf_path: Optional[Path] = None,
    max_genes_per_circ: int = 5,
    legacy_host_gene_source: str = "",
) -> pd.DataFrame:
    """Ensure final host-gene and provenance columns using fixed priority.

    Priority: GTF/GFF overlap, circAtlas host-gene table fields, circAtlas ID
    parsing. Existing non-empty host_gene values are preserved only when no
    derived source is available.
    """
    out = df.copy()
    if gtf_path is not None:
        out = _apply_gtf_host_gene_columns(
            out,
            gtf_path=gtf_path,
            max_genes_per_circ=max_genes_per_circ,
        )

    if "host_gene" not in out.columns:
        out["host_gene"] = ""
    for column in HOST_GENE_COLUMNS:
        if column not in out.columns:
            out[column] = ""

    existing_host = out["host_gene"].map(_clean_annotation_text)
    legacy_host_gene_source = _clean_annotation_text(legacy_host_gene_source)
    if legacy_host_gene_source in {"gtf", "circatlas"}:
        target_column = f"host_gene_from_{legacy_host_gene_source}"
        target = out[target_column].map(_clean_annotation_text)
        fill_mask = (target == "") & (existing_host != "")
        out.loc[fill_mask, target_column] = existing_host[fill_mask]

    circatlas_host_columns = _circatlas_host_gene_columns(out.columns)
    circatlas_id_columns = _circatlas_id_columns(out.columns)

    host_gene_from_circatlas: list[str] = []
    host_gene_from_circatlas_id: list[str] = []
    for _, row in out.iterrows():
        table_values = _split_annotation_values(row.get("host_gene_from_circatlas", ""))
        for column in circatlas_host_columns:
            table_values.extend(_split_annotation_values(row.get(column, "")))
        host_gene_from_circatlas.append(_join_unique_semicolon(table_values))

        id_values = _split_annotation_values(row.get("host_gene_from_circatlas_id", ""))
        for column in circatlas_id_columns:
            parsed = parse_circatlas_id_host_genes(row.get(column, ""))
            id_values.extend(_split_annotation_values(parsed))
        host_gene_from_circatlas_id.append(_join_unique_semicolon(id_values))

    out["host_gene_from_gtf"] = out["host_gene_from_gtf"].map(_clean_annotation_text)
    out["host_gene_from_circatlas"] = host_gene_from_circatlas
    out["host_gene_from_circatlas_id"] = host_gene_from_circatlas_id

    final_host: list[str] = []
    final_source: list[str] = []
    current_source = out["host_gene_source"].map(_clean_annotation_text)
    for idx, row in out.iterrows():
        gtf_gene = _clean_annotation_text(row.get("host_gene_from_gtf", ""))
        circatlas_gene = _clean_annotation_text(row.get("host_gene_from_circatlas", ""))
        circatlas_id_gene = _clean_annotation_text(row.get("host_gene_from_circatlas_id", ""))
        existing_gene = _clean_annotation_text(existing_host.loc[idx])
        source = _clean_annotation_text(current_source.loc[idx])

        if gtf_gene:
            final_host.append(gtf_gene)
            final_source.append("gtf")
        elif circatlas_gene:
            final_host.append(circatlas_gene)
            final_source.append("circatlas")
        elif circatlas_id_gene:
            final_host.append(circatlas_id_gene)
            final_source.append("circatlas_id")
        elif existing_gene:
            final_host.append(existing_gene)
            final_source.append(source if source in HOST_GENE_SOURCES else "")
        else:
            final_host.append("")
            final_source.append("")

    out["host_gene"] = final_host
    out["host_gene_source"] = final_source
    for column in HOST_GENE_COLUMNS:
        out[column] = out[column].map(_clean_annotation_text)
    return out


def annotate_host_gene_frame(
    circ_df: pd.DataFrame,
    *,
    gtf_path: Path,
    max_genes_per_circ: int = 5,
) -> pd.DataFrame:
    annotated = _apply_gtf_host_gene_columns(
        circ_df,
        gtf_path=gtf_path,
        max_genes_per_circ=max_genes_per_circ,
    )
    return normalize_host_gene_annotations(annotated)


def annotate_host_genes(
    circ_feature_table: Path,
    gtf_path: Path,
    out: Optional[Path] = None,
    max_genes_per_circ: int = 5,
) -> None:
    """Annotate circ_feature_table.tsv with host-gene information from GTF/GFF."""
    df_circ = pd.read_csv(circ_feature_table, sep="\t", keep_default_na=False)
    annotated = annotate_host_gene_frame(
        df_circ,
        gtf_path=gtf_path,
        max_genes_per_circ=max_genes_per_circ,
    )

    out_path = out if out is not None else circ_feature_table
    out_path.parent.mkdir(parents=True, exist_ok=True)
    annotated.to_csv(out_path, sep="\t", index=False, lineterminator="\n")


def _host_gene_summary(before: pd.DataFrame, after: pd.DataFrame) -> dict[str, Any]:
    before_host = before["host_gene"].map(_clean_annotation_text) if "host_gene" in before.columns else pd.Series([""] * len(before), index=before.index)
    after_host = after["host_gene"].map(_clean_annotation_text)
    source_counts = {
        str(key): int(value)
        for key, value in after["host_gene_source"].value_counts(dropna=False).items()
        if _clean_annotation_text(key)
    }
    return {
        "n_circRNAs": int(after.shape[0]),
        "n_host_gene_before": int((before_host != "").sum()),
        "n_host_gene_after": int((after_host != "").sum()),
        "n_host_gene_added": int(((before_host == "") & (after_host != "")).sum()),
        "host_gene_source_counts": source_counts,
    }


def repair_host_gene_file(
    *,
    input_path: Path,
    output_path: Path,
    circ_mod: str = "circ",
    gtf_path: Optional[Path] = None,
    overwrite: bool = False,
) -> dict[str, Any]:
    """Patch host-gene provenance columns in an existing h5ad or h5mu file."""
    if output_path.exists() and not overwrite:
        raise ValueError(f"Output already exists: {output_path}. Use --overwrite to replace it.")
    if input_path.resolve() == output_path.resolve():
        raise ValueError("Use a distinct --output path; in-place repair is not supported.")
    suffix = input_path.suffix.lower()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if suffix == ".h5ad":
        import anndata as ad

        adata = ad.read_h5ad(str(input_path))
        before = adata.var.copy()
        adata.var = normalize_host_gene_annotations(
            adata.var,
            gtf_path=gtf_path,
        )
        adata.write_h5ad(str(output_path))
        return {
            "input": str(input_path.resolve()),
            "output": str(output_path.resolve()),
            "file_type": "h5ad",
            "gtf": str(gtf_path.resolve()) if gtf_path is not None else None,
            **_host_gene_summary(before, adata.var),
        }

    if suffix == ".h5mu":
        import mudata as mu

        mdata = mu.read_h5mu(str(input_path))
        if circ_mod not in mdata.mod:
            raise ValueError(f"{input_path} is missing circ modality '{circ_mod}'")
        before = mdata.mod[circ_mod].var.copy()
        mdata.mod[circ_mod].var = normalize_host_gene_annotations(
            mdata.mod[circ_mod].var,
            gtf_path=gtf_path,
        )
        mdata.write_h5mu(str(output_path))
        return {
            "input": str(input_path.resolve()),
            "output": str(output_path.resolve()),
            "file_type": "h5mu",
            "circ_mod": circ_mod,
            "gtf": str(gtf_path.resolve()) if gtf_path is not None else None,
            **_host_gene_summary(before, mdata.mod[circ_mod].var),
        }

    raise ValueError(f"Unsupported input extension '{input_path.suffix}'. Expected .h5ad or .h5mu.")
