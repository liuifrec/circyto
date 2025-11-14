from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd


def _parse_gtf_attrs(attr_str: str) -> Dict[str, str]:
    """Parse the GTF attribute column into a dict."""
    attrs: Dict[str, str] = {}
    for field in attr_str.strip().split(";"):
        field = field.strip()
        if not field:
            continue
        if " " not in field:
            continue
        key, value = field.split(" ", 1)
        value = value.strip().strip('"')
        attrs[key] = value
    return attrs


def load_genes_from_gtf(gtf_path: Path) -> pd.DataFrame:
    """Load gene intervals from a GTF file.

    Prefer 'gene' feature rows; fall back to exon-based span if needed.
    """
    chroms: List[str] = []
    starts: List[int] = []
    ends: List[int] = []
    strands: List[str] = []
    gene_ids: List[str] = []
    gene_names: List[str] = []

    with gtf_path.open() as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, source, feature, start, end, score, strand, frame, attrs = parts

            if feature != "gene":
                continue

            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                continue

            a = _parse_gtf_attrs(attrs)
            gid = a.get("gene_id", "")
            gname = a.get("gene_name", gid)

            chroms.append(chrom)
            starts.append(start_i)
            ends.append(end_i)
            strands.append(strand)
            gene_ids.append(gid)
            gene_names.append(gname)

    df = pd.DataFrame(
        {
            "chrom": chroms,
            "start": starts,
            "end": ends,
            "strand": strands,
            "gene_id": gene_ids,
            "gene_name": gene_names,
        }
    )
    return df


def _find_host_genes_for_circ(
    chrom: str,
    start: Optional[int],
    end: Optional[int],
    strand: str,
    genes: pd.DataFrame,
    max_genes: int = 5,
) -> Tuple[str, str, str, str, int]:
    """Return (host_gene, host_gene_id, multi_names, multi_ids, n) for one circ.

    Strategy:
      1) Same-chrom, same-strand, interval overlap with gene.
      2) If none, same-chrom, any-strand overlap.
      3) If still none, leave empty.
      4) If multiple, report all in multi_* and pick first as primary.
    """
    if not chrom or start is None or end is None:
        return "", "", "", "", 0

    # same chrom
    g_chr = genes[genes["chrom"] == chrom]
    if g_chr.empty:
        return "", "", "", "", 0

    # interval overlap: not (circ_end < gene_start or circ_start > gene_end)
    mask_overlap = ~((end < g_chr["start"]) | (start > g_chr["end"]))
    g_ov = g_chr[mask_overlap]
    if g_ov.empty:
        return "", "", "", "", 0

    # same-strand first
    g_same = g_ov[g_ov["strand"] == strand]
    if not g_same.empty:
        candidates = g_same
    else:
        candidates = g_ov

    # sort by overlap length (descending), then smaller gene span
    overlap_len = (
        candidates[["start", "end"]]
        .apply(lambda r: min(end, r["end"]) - max(start, r["start"]) + 1, axis=1)
    )
    gene_span = candidates["end"] - candidates["start"]
    order = (
        pd.DataFrame({"overlap": overlap_len, "span": gene_span})
        .assign(idx=candidates.index)
        .sort_values(["overlap", "span"], ascending=[False, True])["idx"]
    )

    candidates = candidates.loc[order]

    # take up to max_genes
    top = candidates.head(max_genes)
    names = list(top["gene_name"])
    ids = list(top["gene_id"])

    if not names:
        return "", "", "", "", 0

    primary_name = names[0]
    primary_id = ids[0]
    multi_names = ",".join(names)
    multi_ids = ",".join(ids)
    return primary_name, primary_id, multi_names, multi_ids, len(names)


def annotate_host_genes(
    circ_feature_table: Path,
    gtf_path: Path,
    out: Optional[Path] = None,
    max_genes_per_circ: int = 5,
) -> None:
    """Annotate circ_feature_table.tsv with host gene information from a GTF.

    Input:
      - circ_feature_table: TSV with columns at least [circ_id, chrom, start, end, strand]
      - gtf_path: GTF with 'gene' features

    Output:
      - out (if given) or overwrite circ_feature_table with added columns:
        [host_gene, host_gene_id, host_genes_multi, host_gene_ids_multi, host_gene_n]
    """
    df_circ = pd.read_csv(circ_feature_table, sep="\t")

    # Ensure required columns exist
    for col in ["chrom", "start", "end", "strand"]:
        if col not in df_circ.columns:
            df_circ[col] = "" if col != "start" and col != "end" else None

    # Normalize types
    def _to_int_or_none(x):
        try:
            if pd.isna(x):
                return None
            return int(x)
        except Exception:
            return None

    df_circ["start"] = df_circ["start"].apply(_to_int_or_none)
    df_circ["end"] = df_circ["end"].apply(_to_int_or_none)
    df_circ["chrom"] = df_circ["chrom"].astype(str)
    df_circ["strand"] = df_circ["strand"].astype(str)

    genes = load_genes_from_gtf(gtf_path)

    host_gene: List[str] = []
    host_gene_id: List[str] = []
    host_genes_multi: List[str] = []
    host_gene_ids_multi: List[str] = []
    host_gene_n: List[int] = []

    for _, row in df_circ.iterrows():
        chrom = row["chrom"]
        start = row["start"]
        end = row["end"]
        strand = row["strand"]

        h, h_id, h_multi, h_ids_multi, n = _find_host_genes_for_circ(
            chrom, start, end, strand, genes, max_genes=max_genes_per_circ
        )
        host_gene.append(h)
        host_gene_id.append(h_id)
        host_genes_multi.append(h_multi)
        host_gene_ids_multi.append(h_ids_multi)
        host_gene_n.append(n)

    df_circ["host_gene"] = host_gene
    df_circ["host_gene_id"] = host_gene_id
    df_circ["host_genes_multi"] = host_genes_multi
    df_circ["host_gene_ids_multi"] = host_gene_ids_multi
    df_circ["host_gene_n"] = host_gene_n

    out_path = out if out is not None else circ_feature_table
    df_circ.to_csv(out_path, sep="\t", index=False)
