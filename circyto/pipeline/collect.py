from pathlib import Path
from typing import Dict, List, Tuple, Optional

import pandas as pd
from scipy import sparse
from scipy.io import mmwrite

from ..parsers.cirifull import read_cirifull_tsv


def _gather_cell_tsvs(in_dir: str) -> List[Tuple[str, Path]]:
    """
    Return list of (cell_id, path_to_tsv) for all *.tsv under in_dir.
    """
    base = Path(in_dir)
    if not base.exists():
        return []

    items: List[Tuple[str, Path]] = []
    for tsv in sorted(base.glob("*.tsv")):
        cell_id = tsv.stem
        items.append((cell_id, tsv))
    return items


def _build_matrix_and_features(
    cirifull_dir: str,
) -> Tuple[sparse.csr_matrix, List[str], List[str], Dict[str, Dict[str, Optional[str]]]]:
    """
    Internal: build csr_matrix + circ_ids + cell_ids + feature dict in memory.

    Features dict: circ_id -> {
        "chrom": str | "",
        "start": int | None,
        "end": int | None,
        "strand": str | "",
        "host_gene": str | "",
    }
    """
    pairs = _gather_cell_tsvs(cirifull_dir)
    if not pairs:
        empty = sparse.csr_matrix((0, 0), dtype=int)
        return empty, [], [], {}

    cell_ids: List[str] = [cid for cid, _ in pairs]

    circ_index: Dict[str, int] = {}
    rows: List[int] = []
    cols: List[int] = []
    data: List[int] = []

    features: Dict[str, Dict[str, Optional[str]]] = {}

    for col, (cid, path) in enumerate(pairs):
        df = read_cirifull_tsv(str(path))
        if "circ_id" not in df.columns:
            # Skip files that don't look like CIRI-full outputs
            continue

        # Normalize support column if present: circ_id -> support count
        supports: Dict[str, int] = {}
        if "support" in df.columns:
            for _, r in df.iterrows():
                try:
                    circ = str(r["circ_id"])
                    supports[circ] = int(r["support"])
                except Exception:
                    # Be tolerant of weird/NaN entries
                    continue

        # Iterate rows, register circ_ids + features + matrix entries
        for _, r in df.iterrows():
            circ = str(r["circ_id"])

            if circ not in circ_index:
                circ_index[circ] = len(circ_index)

                chrom = str(r["chr"]) if "chr" in df.columns and pd.notna(r["chr"]) else ""

                # start
                start = None
                if "start" in df.columns and pd.notna(r["start"]):
                    try:
                        start = int(r["start"])
                    except (ValueError, TypeError):
                        start = None

    # end
                end = None
                if "end" in df.columns and pd.notna(r["end"]):
                    try:
                        end = int(r["end"])
                    except (ValueError, TypeError):
                        end = None

                strand = str(r["strand"]) if "strand" in df.columns and pd.notna(r["strand"]) else ""

                host_gene = (
                    str(r["host_gene"])
                    if "host_gene" in df.columns and pd.notna(r["host_gene"])
                    else ""
                )


                features[circ] = {
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "host_gene": host_gene,
                }

            rows.append(circ_index[circ])
            cols.append(col)
            val = supports.get(circ, 1)
            data.append(val)

    if not data:
        # No non-zero entries, but we may still have cells; keep consistent shapes
        X_empty = sparse.csr_matrix((0, len(cell_ids)), dtype=int)
        return X_empty, [], cell_ids, {}

    n_rows = len(circ_index)
    n_cols = len(cell_ids)
    X = sparse.csr_matrix((data, (rows, cols)), shape=(n_rows, n_cols), dtype=int)

    circ_ids: List[str] = [None] * n_rows  # type: ignore[assignment]
    for circ, i in circ_index.items():
        circ_ids[i] = circ

    return X, circ_ids, cell_ids, features


def collect_matrix(
    cirifull_dir: str,
    matrix_path: str,
    circ_index_path: str,
    cell_index_path: str,
    min_count_per_cell: int = 1,
) -> None:
    """
    Build circ-by-cell count matrix from per-cell CIRI-full TSVs and write:

      - Matrix Market sparse matrix (.mtx)
      - circ index (rows) as one ID per line
      - cell index (cols) as one ID per line
      - circ_feature_table.tsv (per-circ features; same dir as matrix)

    Applies a simple per-cell filter: keep columns with sum >= min_count_per_cell.
    """
    X, circ_ids, cell_ids, features = _build_matrix_and_features(cirifull_dir)

    # Filter by per-cell totals (on the matrix)
    if X.shape[1] > 0 and min_count_per_cell > 1:
        keep_cols = X.sum(axis=0).A1 >= min_count_per_cell
        if keep_cols.any():
            X = X[:, keep_cols]
            cell_ids = [cid for cid, keep in zip(cell_ids, keep_cols) if keep]
        else:
            # Nothing passes filter; keep an empty matrix with 0 cells
            X = sparse.csr_matrix((X.shape[0], 0), dtype=int)
            cell_ids = []

    # Write matrix
    out_mtx = Path(matrix_path)
    out_mtx.parent.mkdir(parents=True, exist_ok=True)
    mmwrite(out_mtx, X)

    # Write circ & cell indices
    Path(circ_index_path).write_text(
        "\n".join(circ_ids) + ("\n" if circ_ids else "")
    )
    Path(cell_index_path).write_text(
        "\n".join(cell_ids) + ("\n" if cell_ids else "")
    )

    # Write circ_feature_table.tsv (if we have any circ features)
    if circ_ids:
        feature_rows = []
        for circ in circ_ids:
            f = features.get(circ, {})
            feature_rows.append(
                {
                    "circ_id": circ,
                    "chrom": f.get("chrom", ""),
                    "start": f.get("start", ""),
                    "end": f.get("end", ""),
                    "strand": f.get("strand", ""),
                    "host_gene": f.get("host_gene", ""),
                }
            )

        df_feat = pd.DataFrame(feature_rows)

        # Fixed filename next to the matrix
        feat_path = out_mtx.with_name("circ_feature_table.tsv")
        df_feat.to_csv(feat_path, sep="\t", index=False)
