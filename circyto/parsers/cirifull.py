from pathlib import Path
from typing import List, Dict, Tuple
import pandas as pd

def discover_plate_tsvs(cirifull_dir: Path) -> List[Tuple[str, Path]]:
items = []
for tsv in sorted(cirifull_dir.glob("*.tsv")):
cell_id = tsv.stem
if not cell_id.startswith("batch_"):
items.append((cell_id, tsv))
return items

def discover_batch_results(cirifull_dir: Path) -> List[dict]:
items = []
for tsv in sorted(cirifull_dir.glob("batch_*.tsv")):
batch_name = tsv.stem
prep_dir = cirifull_dir.parent / "fastq_by_cell" / batch_name
cells = prep_dir / "cells.txt"
if cells.exists():
items.append({"tsv": tsv, "cells": cells})
return items

def parse_cell_tsv(tsv_path: Path) -> Dict[str, int]:
df = pd.read_csv(tsv_path, sep="\t", comment="#")
if df.shape[1] == 0:
return {}
cand_id = [c for c in df.columns if c.lower() in ("circ_id","circid","circ","id","bsj","name")]
circ_col = cand_id[0] if cand_id else df.columns[0]
cand_cnt = [c for c in df.columns if c.lower() in ("count","reads","support","num_supports","supporting_reads")]
count_col = cand_cnt[0] if cand_cnt else df.columns[-1]
df[circ_col] = df[circ_col].astype(str).str.replace(",", "|", regex=False)
grouped = df.groupby(circ_col)[count_col].sum()
return grouped.to_dict()
