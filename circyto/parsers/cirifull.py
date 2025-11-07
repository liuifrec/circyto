# /workspaces/circyto/circyto/parsers/cirifull.py
from __future__ import annotations

import pandas as pd
from pathlib import Path

CIRI_COLUMNS = ["circ_id", "chr", "start", "end", "strand", "support"]


def read_cirifull_tsv(path: str | Path) -> pd.DataFrame:
    """
    Read a CIRI-full output TSV file and normalize column names.

    Returns a DataFrame with the canonical columns:
    circ_id, chr, start, end, strand, support
    """
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        # Return an empty DataFrame with canonical columns
        return pd.DataFrame(columns=CIRI_COLUMNS)

    try:
        df = pd.read_csv(path, sep="\t", comment="#")
    except Exception as e:
        print(f"[circyto] Failed to parse {path}: {e}")
        return pd.DataFrame(columns=CIRI_COLUMNS)

    # Normalize lowercase
    lower = {c.lower(): c for c in df.columns}
    out = pd.DataFrame(columns=CIRI_COLUMNS)
    for c in CIRI_COLUMNS:
        if c in df.columns:
            out[c] = df[c]
        elif c in lower:
            out[c] = df[lower[c]]
        else:
            out[c] = ""

    return out
