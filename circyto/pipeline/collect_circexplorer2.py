# circyto/pipeline/collect_circexplorer2.py

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple, Optional


@dataclass
class CircRecord:
    chrom: str
    start: int
    end: int
    strand: str
    gene_name: Optional[str]
    # total support across all cells (optional summary)
    total_support: int = 0


def _parse_circularRNA_known(path: Path, min_support: int) -> Dict[str, Tuple[CircRecord, int]]:
    """
    Parse a CIRCexplorer2 circularRNA_known.txt file for a single cell.

    Returns a mapping:
      circ_id -> (CircRecord, support_in_this_cell)

    Format (CIRCexplorer2 annotate):
      chrom, start, end, name, score, strand, thickStart, thickEnd, itemRgb,
      exonCount, exonSizes, exonOffsets, readNumber, circType,
      geneName, isoformName, index, flankIntron
    """
    per_circ: Dict[str, Tuple[CircRecord, int]] = {}

    if not path.exists():
        return per_circ

    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 18:
                continue

            chrom = parts[0]
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue

            strand = parts[5]
            try:
                read_number = int(parts[12])
            except ValueError:
                # if readNumber is malformed, skip
                continue

            if read_number < min_support:
                continue

            gene_name = parts[14] if parts[14] not in {"", "NA", "None"} else None

            # We define a detector-agnostic circ_id: chr:start|end
            circ_id = f"{chrom}:{start}|{end}"

            rec = CircRecord(
                chrom=chrom,
                start=start,
                end=end,
                strand=strand,
                gene_name=gene_name,
                total_support=read_number,
            )
            per_circ[circ_id] = (rec, read_number)

    return per_circ


def collect_circexplorer2_matrix(
    indir: Path,
    outdir: Path,
    min_support: int = 1,
) -> None:
    """
    Collect per-cell CIRCexplorer2 outputs into a sparse MatrixMarket circ × cell matrix.

    Expects a directory layout like:

      indir/
        CELL_A/
          circularRNA_known.txt
        CELL_B/
          circularRNA_known.txt
        ...

    Produces in outdir:
      - circ_counts.mtx
      - circ_index.txt
      - cell_index.txt
      - circ_feature_table.tsv
    """
    indir = Path(indir)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # 1) Discover cells and parse their circularRNA_known.txt
    cell_ids: List[str] = []
    cell_circ_maps: Dict[str, Dict[str, Tuple[CircRecord, int]]] = {}

    for sub in sorted(indir.iterdir()):
        if not sub.is_dir():
            continue
        cell_id = sub.name
        circ_file = sub / "circularRNA_known.txt"
        if not circ_file.exists():
            continue

        circ_map = _parse_circularRNA_known(circ_file, min_support=min_support)
        if not circ_map:
            # allow empty cells; they just won't contribute rows
            cell_ids.append(cell_id)
            cell_circ_maps[cell_id] = {}
            continue

        cell_ids.append(cell_id)
        cell_circ_maps[cell_id] = circ_map

    if not cell_ids:
        raise RuntimeError(f"No per-cell CIRCexplorer2 outputs found in {indir}")

    # 2) Build a global circ_id index and aggregate metadata
    circ_meta: Dict[str, CircRecord] = {}
    for cell_id, circ_map in cell_circ_maps.items():
        for circ_id, (rec, support) in circ_map.items():
            if circ_id not in circ_meta:
                circ_meta[circ_id] = CircRecord(
                    chrom=rec.chrom,
                    start=rec.start,
                    end=rec.end,
                    strand=rec.strand,
                    gene_name=rec.gene_name,
                    total_support=support,
                )
            else:
                circ_meta[circ_id].total_support += support

    if not circ_meta:
        raise RuntimeError(
            f"No circRNAs passed min_support={min_support} across any cells in {indir}"
        )

    # Sort circ_ids for deterministic output
    circ_ids_sorted = sorted(
        circ_meta.keys(),
        key=lambda cid: (
            circ_meta[cid].chrom,
            circ_meta[cid].start,
            circ_meta[cid].end,
        ),
    )
    circ_index = {cid: idx for idx, cid in enumerate(circ_ids_sorted, start=1)}
    cell_index = {cid: idx for idx, cid in enumerate(cell_ids, start=1)}

    # 3) Build sparse triplets (i, j, value)
    triplets: List[Tuple[int, int, int]] = []
    for cell_id, circ_map in cell_circ_maps.items():
        j = cell_index[cell_id]
        for circ_id, (_rec, support) in circ_map.items():
            i = circ_index[circ_id]
            triplets.append((i, j, support))

    n_circ = len(circ_ids_sorted)
    n_cells = len(cell_ids)
    nnz = len(triplets)

    if nnz == 0:
        raise RuntimeError(
            f"No non-zero entries in circ × cell matrix from {indir} "
            f"(after min_support={min_support})."
        )

    # 4) Write MatrixMarket file
    mtx_path = outdir / "circ_counts.mtx"
    with mtx_path.open("w") as f:
        f.write("%%MatrixMarket matrix coordinate integer general\n")
        f.write("% Generated by circyto collect-circexplorer2\n")
        f.write(f"{n_circ} {n_cells} {nnz}\n")
        for i, j, v in sorted(triplets):
            f.write(f"{i} {j} {v}\n")

    # 5) Write circ_index and cell_index
    circ_index_path = outdir / "circ_index.txt"
    with circ_index_path.open("w") as f:
        for cid in circ_ids_sorted:
            f.write(f"{cid}\n")

    cell_index_path = outdir / "cell_index.txt"
    with cell_index_path.open("w") as f:
        for cid in cell_ids:
            f.write(f"{cid}\n")

    # 6) Write circ_feature_table.tsv
    feat_path = outdir / "circ_feature_table.tsv"
    with feat_path.open("w") as f:
        f.write(
            "circ_id\tchrom\tstart\tend\tstrand\tgene_name\ttotal_support\n"
        )
        for cid in circ_ids_sorted:
            rec = circ_meta[cid]
            gene_name = rec.gene_name if rec.gene_name is not None else ""
            f.write(
                f"{cid}\t{rec.chrom}\t{rec.start}\t{rec.end}\t"
                f"{rec.strand}\t{gene_name}\t{rec.total_support}\n"
            )
