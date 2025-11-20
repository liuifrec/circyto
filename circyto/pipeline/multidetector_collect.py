# circyto/pipeline/multidetector_collect.py
from __future__ import annotations

import pandas as pd
from pathlib import Path
from typing import Dict, List
from scipy import sparse
import scipy.io


def parse_detector_tsv(tsv_path: Path):
    """
    Returns: (list_of_rows)
    Each row is (circ_id, chrom, start, end, strand, support)
    """
    if not tsv_path.exists() or tsv_path.stat().st_size == 0:
        return []

    df = pd.read_csv(tsv_path, sep="\t")

    fields = ["circ_id", "chr", "start", "end", "strand", "support"]
    for f in fields:
        if f not in df:
            raise ValueError(f"Missing column {f} in {tsv_path}")

    rows = []
    for _, r in df.iterrows():
        circ = str(r["circ_id"])
        chrom = str(r["chr"])
        start = int(r["start"])
        end = int(r["end"])
        strand = str(r["strand"])
        support = int(r["support"])
        rows.append((circ, chrom, start, end, strand, support))

    return rows


def build_detector_matrix(detector_tsv_dir: Path):
    """
    Build circRNA x cell sparse matrix from detector TSV directory.
    """
    tsvs = sorted(detector_tsv_dir.glob("*.tsv"))
    cell_ids = [f.stem for f in tsvs]

    circ_index = {}
    data = []
    rows = []
    cols = []

    for col, tsv in enumerate(tsvs):
        parsed = parse_detector_tsv(tsv)

        for circ_id, chrom, start, end, strand, support in parsed:
            if circ_id not in circ_index:
                circ_index[circ_id] = len(circ_index)

            row = circ_index[circ_id]
            rows.append(row)
            cols.append(col)
            data.append(support)

    if len(circ_index) == 0:
        X = sparse.csr_matrix((0, len(cell_ids)))
    else:
        X = sparse.csr_matrix((data, (rows, cols)),
                              shape=(len(circ_index), len(cell_ids)))

    circ_ids = list(circ_index.keys())

    return X, circ_ids, cell_ids


def write_matrix(X, circ_ids, cell_ids, out_prefix: Path):
    """
    Write:
      <prefix>.mtx
      <prefix>_circ_ids.txt
      <prefix>_cell_ids.txt
    """
    mtx_path = out_prefix.with_suffix(".mtx")
    circ_path = out_prefix.parent / (out_prefix.name + "_circ_ids.txt")
    cell_path = out_prefix.parent / (out_prefix.name + "_cell_ids.txt")

    scipy.io.mmwrite(str(mtx_path), X)
    circ_path.write_text("\n".join(circ_ids) + "\n")
    cell_path.write_text("\n".join(cell_ids) + "\n")

    return mtx_path, circ_path, cell_path
