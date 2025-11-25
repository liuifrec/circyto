#!/usr/bin/env python
from __future__ import annotations

from pathlib import Path
from typing import List

import pandas as pd


def parse_splice_sites(path: Path, cell_id: str) -> pd.DataFrame:
    """
    Parse a find_circ3 splice_sites.bed file and add circ_id + cell_id.

    Assumes 18-column legacy format:
      chrom, start, end, name, n_reads, strand,
      bridge_reads, uniq_bridges,
      mapq1, mapq2,
      sample,
      n_uniq_anchors,
      flag1, flag2, flag3, flag4,
      splice_signal,
      tags
    """
    cols = [
        "chrom",
        "start",
        "end",
        "name",
        "n_reads",
        "strand",
        "bridge_reads",
        "uniq_bridges",
        "mapq1",
        "mapq2",
        "sample",
        "n_uniq_anchors",
        "flag1",
        "flag2",
        "flag3",
        "flag4",
        "splice_signal",
        "tags",
    ]
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=cols,
        dtype={"chrom": str},
    )
    df["cell_id"] = cell_id

    # Standard circ_id format for circyto: chr:start|end|strand
    df["circ_id"] = (
        df["chrom"].astype(str)
        + ":"
        + df["start"].astype(str)
        + "|"
        + df["end"].astype(str)
        + "|"
        + df["strand"].astype(str)
    )
    return df


def collect_find_circ3(outdir_root: Path) -> None:
    """
    Collect all <cell_id>_splice_sites.bed under outdir_root
    and build:
      - find_circ3_circ_feature_table.tsv
      - find_circ3_circ_counts.tsv (circ x cell counts)
    """
    outdir_root = outdir_root.resolve()
    cell_dirs = [p for p in outdir_root.iterdir() if p.is_dir()]

    all_dfs: List[pd.DataFrame] = []
    for cdir in sorted(cell_dirs):
        cell_id = cdir.name
        bed = cdir / f"{cell_id}_splice_sites.bed"
        if not bed.exists():
            print(f"[collect_find_circ3] WARNING: missing {bed}, skipping")
            continue
        df = parse_splice_sites(bed, cell_id)
        all_dfs.append(df)

    if not all_dfs:
        raise SystemExit("[collect_find_circ3] No splice_sites.bed files found.")

    df_all = pd.concat(all_dfs, ignore_index=True)

    # Feature table: one row per circ_id
    feat_cols = [
        "circ_id",
        "chrom",
        "start",
        "end",
        "strand",
        "splice_signal",
        "tags",
    ]
    feature_df = (
        df_all[feat_cols]
        .drop_duplicates("circ_id")
        .sort_values(["chrom", "start", "end"])
        .reset_index(drop=True)
    )

    feature_path = outdir_root / "find_circ3_circ_feature_table.tsv"
    feature_df.to_csv(feature_path, sep="\t", index=False)
    print(f"[collect_find_circ3] Wrote feature table: {feature_path}")

    # Counts matrix: sum n_reads per circ_id x cell_id
    counts_df = (
        df_all.groupby(["circ_id", "cell_id"])["n_reads"]
        .sum()
        .reset_index()
        .pivot(index="circ_id", columns="cell_id", values="n_reads")
        .fillna(0)
        .astype(int)
    )

    counts_path = outdir_root / "find_circ3_circ_counts.tsv"
    counts_df.to_csv(counts_path, sep="\t")
    print(f"[collect_find_circ3] Wrote counts matrix: {counts_path}")


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(
        description=(
            "Collect find_circ3 outputs into a circ x cell matrix "
            "for a small manifest (e.g. manifest_2.tsv)."
        )
    )
    ap.add_argument(
        "outdir_root",
        type=Path,
        help=(
            "Output directory used for run-detector "
            "(e.g. tests/find_circ3_manifest2)"
        ),
    )
    args = ap.parse_args()
    collect_find_circ3(args.outdir_root)
