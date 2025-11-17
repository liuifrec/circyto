# circyto/compare/detector_compare.py
from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import scipy.sparse as sp
from scipy.io import mmwrite


def _ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def _load_circ_ids_from_tsv(tsv_path: Path) -> Set[str]:
    """
    Load circRNA IDs from a detector TSV.
    We prefer a 'circ_id' column, but if missing, we fall back to the first column.
    """
    circs: Set[str] = set()
    with tsv_path.open() as f:
        rd = csv.DictReader(f, delimiter="\t")
        fieldnames = rd.fieldnames or []
        if not fieldnames:
            return circs

        if "circ_id" in fieldnames:
            key = "circ_id"
        else:
            # fall back to first column
            key = fieldnames[0]

        for row in rd:
            cid = row.get(key)
            if cid:
                circs.add(str(cid))
    return circs


def _collect_circ_sets(
    root_dir: Path,
    detector_names: List[str],
) -> Dict[str, Set[str]]:
    """
    For each detector name, look under root_dir / <detector> for *.tsv files
    and union all circ_ids across cells.
    """
    circ_sets: Dict[str, Set[str]] = {}

    for det in detector_names:
        det_dir = root_dir / det
        if not det_dir.exists():
            raise FileNotFoundError(f"Detector directory not found: {det_dir}")

        circs: Set[str] = set()
        for tsv in det_dir.glob("*.tsv"):
            circs |= _load_circ_ids_from_tsv(tsv)
        circ_sets[det] = circs

    return circ_sets


def _build_presence_matrix(
    circ_sets: Dict[str, Set[str]]
) -> Tuple[sp.csr_matrix, List[str], List[str]]:
    """
    Build a sparse circ × detector presence matrix (1 if detector detects circ, else 0).
    """
    detector_names = sorted(circ_sets.keys())
    all_circs: Set[str] = set()
    for s in circ_sets.values():
        all_circs |= s

    circ_ids = sorted(all_circs)

    n_circ = len(circ_ids)
    n_det = len(detector_names)

    if n_circ == 0 or n_det == 0:
        mat = sp.csr_matrix((n_circ, n_det), dtype=int)
        return mat, circ_ids, detector_names

    circ_index = {cid: i for i, cid in enumerate(circ_ids)}
    det_index = {d: j for j, d in enumerate(detector_names)}

    rows = []
    cols = []
    data = []

    for det, circs in circ_sets.items():
        j = det_index[det]
        for cid in circs:
            i = circ_index[cid]
            rows.append(i)
            cols.append(j)
            data.append(1)

    mat = sp.csr_matrix((data, (rows, cols)), shape=(n_circ, n_det), dtype=int)
    return mat, circ_ids, detector_names


def _compute_summary(
    circ_sets: Dict[str, Set[str]],
    union_set: Set[str],
    intersection_set: Set[str],
) -> Dict:
    det_names = sorted(circ_sets.keys())
    summary = {
        "n_detectors": len(det_names),
        "detectors": {},
        "union_size": len(union_set),
        "intersection_size": len(intersection_set),
        "jaccard": {},
    }

    for d in det_names:
        s = circ_sets[d]
        summary["detectors"][d] = {
            "n_circs": len(s),
            "fraction_of_union": (len(s) / len(union_set)) if union_set else 0.0,
        }

    # pairwise Jaccard
    for i, d1 in enumerate(det_names):
        for d2 in det_names[i + 1 :]:
            s1 = circ_sets[d1]
            s2 = circ_sets[d2]
            inter = len(s1 & s2)
            uni = len(s1 | s2)
            j = (inter / uni) if uni else 0.0
            summary["jaccard"][f"{d1}|{d2}"] = j

    return summary


def compare_detectors_from_root(
    root_dir: Path,
    detector_names: List[str],
    outdir: Path,
) -> Dict:
    """
    High-level entry point:

      - root_dir: directory containing per-detector subdirs:
          root_dir/
            ciri-full/
              cell1.tsv
              cell2.tsv
            ciri-long/
              cell1.tsv
              cell2.tsv
      - detector_names: list of detector subdir names
      - outdir: where to write comparison outputs

    Outputs:
      - circ_detector.mtx       (MatrixMarket, circ × detector)
      - circ_ids.txt            (one circ_id per line)
      - detectors.txt           (one detector name per line, in column order)
      - union.tsv               (circ_id column)
      - intersection.tsv        (circ_id column)
      - summary.json            (stats & Jaccard)
    """
    _ensure_dir(outdir)

    circ_sets = _collect_circ_sets(root_dir, detector_names)

    # union & intersection
    values = list(circ_sets.values())
    union_set: Set[str] = set().union(*values) if values else set()
    intersection_set: Set[str] = set.intersection(*values) if values else set()

    # presence matrix
    mat, circ_ids, det_names_sorted = _build_presence_matrix(circ_sets)

    # write matrix + indices
    mtx_path = outdir / "circ_detector.mtx"
    mmwrite(str(mtx_path), mat)

    circ_ids_path = outdir / "circ_ids.txt"
    circ_ids_path.write_text("\n".join(circ_ids) + ("\n" if circ_ids else ""))

    dets_path = outdir / "detectors.txt"
    dets_path.write_text("\n".join(det_names_sorted) + ("\n" if det_names_sorted else ""))

    # union / intersection TSVs
    union_path = outdir / "union.tsv"
    inter_path = outdir / "intersection.tsv"

    with union_path.open("w") as f:
        f.write("circ_id\n")
        for cid in sorted(union_set):
            f.write(f"{cid}\n")

    with inter_path.open("w") as f:
        f.write("circ_id\n")
        for cid in sorted(intersection_set):
            f.write(f"{cid}\n")

    # summary JSON
    summary = _compute_summary(circ_sets, union_set, intersection_set)
    summary_path = outdir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True))

    return summary
