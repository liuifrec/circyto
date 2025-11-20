# circyto/utils.py
from __future__ import annotations
from pathlib import Path

def ensure_dir(p: str | Path) -> Path:
    p = Path(p)
    p.mkdir(parents=True, exist_ok=True)
    return p
