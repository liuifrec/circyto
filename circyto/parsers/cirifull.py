# /workspaces/circyto/circyto/parsers/cirifull.py
from __future__ import annotations

from pathlib import Path
import pandas as pd

REQUIRED_CIRI_COLUMNS = ["circ_id", "chr", "start", "end", "strand", "support"]
OPTIONAL_CIRI_COLUMNS = ["host_gene", "gene_name"]


def read_cirifull_tsv(path: str | Path) -> pd.DataFrame:
    """
    Read a CIRI-full output TSV file and normalize column names.

    Returns a DataFrame with the canonical columns:
    circ_id, chr, start, end, strand, support

    Preserves optional annotation columns when present:
    host_gene, gene_name
    """
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        # Return an empty DataFrame with canonical columns
        return pd.DataFrame(columns=REQUIRED_CIRI_COLUMNS + OPTIONAL_CIRI_COLUMNS)

    try:
        df = pd.read_csv(path, sep="\t", comment="#", keep_default_na=False)
    except Exception as e:
        print(f"[circyto] Failed to parse {path}: {e}")
        return pd.DataFrame(columns=REQUIRED_CIRI_COLUMNS + OPTIONAL_CIRI_COLUMNS)

    # Normalize lowercase
    lower = {c.lower(): c for c in df.columns}
    out = pd.DataFrame(columns=REQUIRED_CIRI_COLUMNS + OPTIONAL_CIRI_COLUMNS)
    for c in REQUIRED_CIRI_COLUMNS:
        if c in df.columns:
            out[c] = df[c]
        elif c in lower:
            out[c] = df[lower[c]]
        else:
            out[c] = ""
    for c in OPTIONAL_CIRI_COLUMNS:
        if c in df.columns:
            out[c] = df[c]
        elif c in lower:
            out[c] = df[lower[c]]
        else:
            out[c] = ""

    return out
