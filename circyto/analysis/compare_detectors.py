# circyto/analysis/compare_detectors.py

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
import csv

@dataclass(frozen=True, order=True)
class CircKey:
    chrom: str
    start: int
    end: int
    strand: str  # "+", "-", "?"

# Accept:
#   chr1:100-200:+
#   chr1:100:200:+
#   chr1:100|200|-
#   chr1:100-200
#   chr1:100:200
_CIRC_ID_RE = re.compile(
    r"^(?P<chrom>[^:]+):(?P<start>\d+)(?:[-:|])(?P<end>\d+)(?:[:|](?P<strand>[+\-?]))?$"
)

def parse_circ_id(s: str, default_strand: str = "?") -> CircKey:
    s = s.strip()
    m = _CIRC_ID_RE.match(s)
    if not m:
        raise ValueError(f"Unrecognized circ_id format: {s}")
    return CircKey(
        chrom=m.group("chrom"),
        start=int(m.group("start")),
        end=int(m.group("end")),
        strand=m.group("strand") or default_strand,
    )

def load_keys(path: str | Path, col: str = "circ_id") -> set[CircKey]:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)

    keys: set[CircKey] = set()

    # circ_index-style
    if path.suffix in {".txt", ".list"}:
        for line in path.read_text().splitlines():
            line = line.strip()
            if line:
                keys.add(parse_circ_id(line))
        return keys

    # TSV
    with path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not reader.fieldnames or col not in reader.fieldnames:
            raise ValueError(f"Missing column '{col}' in {path} (has {reader.fieldnames})")
        for row in reader:
            v = (row.get(col) or "").strip()
            if v:
                keys.add(parse_circ_id(v))
    return keys

def fuzzy_hits(a: set[CircKey], b: set[CircKey], window: int = 5) -> int:
    """
    Count how many elements in A have at least one fuzzy match in B (±window on start/end).

    Strand rules:
      - If either side is '?', treat as wildcard (strand-agnostic match).
      - Otherwise require strand equality.
    """
    # index B by chrom only (strand handled by wildcard logic below)
    idx: dict[str, list[CircKey]] = {}
    for k in b:
        idx.setdefault(k.chrom, []).append(k)

    hits = 0
    for ka in a:
        cand = idx.get(ka.chrom, [])
        ok = any(
            (ka.strand == "?" or kb.strand == "?" or ka.strand == kb.strand)
            and abs(ka.start - kb.start) <= window
            and abs(ka.end - kb.end) <= window
            for kb in cand
        )
        if ok:
            hits += 1
    return hits


def fuzzy_jaccard(a: set[CircKey], b: set[CircKey], window: int = 5) -> float:
    ha = fuzzy_hits(a, b, window)
    hb = fuzzy_hits(b, a, window)
    i = min(ha, hb)
    u = len(a) + len(b) - i
    return (i / u) if u else 0.0
