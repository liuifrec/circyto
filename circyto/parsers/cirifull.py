from pathlib import Path
from typing import List, Dict
import pandas as pd
def discover_cell_results(cirifull_dir: Path) -> List[dict]:
    items = []
    for tsv in sorted(cirifull_dir.glob('batch_*.tsv')):
        batch_name = tsv.stem
        prep_dir = cirifull_dir.parent / 'fastq_by_cell' / batch_name
        cells = prep_dir / 'cells.txt'
        if cells.exists(): items.append({'tsv': tsv, 'cells': cells})
    if not items:
        for tsv in cirifull_dir.glob('*.tsv'):
            cells = tsv.with_suffix('').parent / 'cells.txt'
            if cells.exists(): items.append({'tsv': tsv, 'cells': cells})
    if not items: raise FileNotFoundError('No CIRI-full TSVs paired with cells.txt were found.')
    return items
def parse_cell_tsv(tsv_path: Path) -> Dict[str, int]:
    df = pd.read_csv(tsv_path, sep='\t', comment='#')
    if df.shape[1] == 0: return {}
    cand_id = [c for c in df.columns if c.lower() in ('circ_id','circid','circ','id')]
    circ_col = cand_id[0] if cand_id else df.columns[0]
    cand_cnt = [c for c in df.columns if c.lower() in ('count','reads','support')]
    count_col = cand_cnt[0] if cand_cnt else df.columns[-1]
    df[circ_col] = df[circ_col].astype(str).str.replace(',', '|', regex=False)
    grouped = df.groupby(circ_col)[count_col].sum()
    return grouped.to_dict()
