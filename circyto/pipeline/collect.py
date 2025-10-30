from pathlib import Path
from rich.console import Console
import numpy as np
from scipy import sparse
from scipy.io import mmwrite
from ..parsers.cirifull import discover_cell_results, parse_cell_tsv
console = Console()
def collect_matrix(cirifull_dir: Path, matrix_out: Path, circ_index_out: Path, cell_index_out: Path, min_count_per_cell: int = 1):
    discovered = discover_cell_results(cirifull_dir)
    circ_set = set()
    for item in discovered:
        circ_counts = parse_cell_tsv(item['tsv'])
        circ_set.update(circ_counts.keys())
    circ_list = sorted(circ_set)
    circ_to_row = {c:i for i,c in enumerate(circ_list)}
    cell_list = []
    for item in discovered:
        with open(item['cells'], 'r') as f:
            cell_list.extend([line.strip() for line in f if line.strip()])
    cell_list = sorted(set(cell_list))
    cell_to_col = {c:i for i,c in enumerate(cell_list)}
    rows, cols, data = [], [], []
    for item in discovered:
        circ_counts = parse_cell_tsv(item['tsv'])
        with open(item['cells'], 'r') as f:
            cells = [line.strip() for line in f if line.strip()]
        ncell = max(1, len(cells))
        for circ_id, count in circ_counts.items():
            r = circ_to_row[circ_id]
            share = float(count) / ncell
            for cb in cells:
                c = cell_to_col[cb]
                rows.append(r); cols.append(c); data.append(share)
    M = sparse.coo_matrix((data, (rows, cols)), shape=(len(circ_list), len(cell_list))).tocsr()
    keep = np.array(M.sum(axis=0)).ravel() >= min_count_per_cell
    M = M[:, keep]
    cell_list = [c for c,k in zip(cell_list, keep) if k]
    matrix_out.parent.mkdir(parents=True, exist_ok=True)
    mmwrite(matrix_out, M)
    Path(circ_index_out).write_text('\n'.join(circ_list))
    Path(cell_index_out).write_text('\n'.join(cell_list))
    console.print(f'[green]Collected[/green] matrix: {M.shape[0]} circ × {M.shape[1]} cells')
