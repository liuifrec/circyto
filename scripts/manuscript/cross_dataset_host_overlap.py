#!/usr/bin/env python
from __future__ import annotations

import argparse
from itertools import combinations
from pathlib import Path

import pandas as pd

from _mudata_utils import (
    CIRC_MODALITY_CANDIDATES,
    dataset_name_from_path,
    fail,
    find_modality,
    first_existing_column,
    host_gene_set,
    load_mudata,
    read_table,
    split_host_genes,
    write_table,
)


def parse_label_path(value: str) -> tuple[str, Path]:
    if "=" not in value:
        fail(f"Expected LABEL=PATH, got: {value}")
    label, path = value.split("=", 1)
    return label, Path(path)


def load_host_genes(path: Path) -> set[str]:
    mdata = load_mudata(path)
    circ = mdata.mod[find_modality(mdata, CIRC_MODALITY_CANDIDATES)]
    return host_gene_set(circ)


def positive_program_genes(path: Path, host_gene_col: str | None, effect_col: str | None) -> set[str]:
    frame = read_table(path)
    gene_col = host_gene_col or first_existing_column(frame, ["host_gene", "gene", "gene_symbol", "host_gene_name"])
    if gene_col is None:
        fail(f"Could not find a host-gene column in enrichment table: {path}")
    effect = effect_col or first_existing_column(frame, ["log2_fc_pseudocount1", "log2_fc", "difference", "effect"])
    if effect is None:
        selected = frame
    else:
        values = pd.to_numeric(frame[effect], errors="coerce")
        selected = frame.loc[values > 0]
    genes: set[str] = set()
    for value in selected[gene_col]:
        genes.update(split_host_genes(value))
    return genes


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Summarize cross-dataset circRNA host-gene overlap.")
    parser.add_argument("h5mu", nargs="+", type=Path, help="Input MuData .h5mu files.")
    parser.add_argument(
        "--dataset-name",
        action="append",
        default=None,
        help="Optional dataset label. Repeat once per input, in the same order.",
    )
    parser.add_argument("--outdir", required=True, type=Path, help="Directory for overlap outputs.")
    parser.add_argument(
        "--enrichment-table",
        action="append",
        default=[],
        help="Optional positive-program table as LABEL=PATH. Repeat for HAP1/IMR90 tables.",
    )
    parser.add_argument("--host-gene-column", default=None, help="Optional host-gene column for enrichment tables.")
    parser.add_argument("--effect-column", default=None, help="Optional numeric effect column for enrichment tables.")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    if args.dataset_name and len(args.dataset_name) != len(args.h5mu):
        fail("--dataset-name must be repeated once per input .h5mu file.")
    labels = [
        args.dataset_name[i] if args.dataset_name else dataset_name_from_path(path)
        for i, path in enumerate(args.h5mu)
    ]
    gene_sets = {label: load_host_genes(path) for label, path in zip(labels, args.h5mu)}
    args.outdir.mkdir(parents=True, exist_ok=True)

    pair_rows = []
    for left, right in combinations(labels, 2):
        intersection = gene_sets[left] & gene_sets[right]
        union = gene_sets[left] | gene_sets[right]
        pair_rows.append(
            {
                "dataset_a": left,
                "dataset_b": right,
                "n_host_genes_a": len(gene_sets[left]),
                "n_host_genes_b": len(gene_sets[right]),
                "n_overlap": len(intersection),
                "n_union": len(union),
                "jaccard": len(intersection) / len(union) if union else 0.0,
            }
        )
    write_table(pd.DataFrame(pair_rows), args.outdir / "pairwise_host_gene_overlap.tsv")

    if labels:
        shared = set.intersection(*(gene_sets[label] for label in labels))
    else:
        shared = set()
    write_table(
        pd.DataFrame({"host_gene": sorted(shared), "n_datasets": len(labels)}),
        args.outdir / "three_way_host_gene_overlap.tsv",
    )

    if args.enrichment_table:
        program_sets = {}
        for item in args.enrichment_table:
            label, table_path = parse_label_path(item)
            program_sets[label] = positive_program_genes(table_path, args.host_gene_column, args.effect_column)
        shared_program = set.intersection(*program_sets.values()) if program_sets else set()
        shared_program &= shared if shared else set.union(*(gene_sets[label] for label in labels))
        write_table(
            pd.DataFrame(
                {
                    "host_gene": sorted(shared_program),
                    "positive_in_tables": ";".join(program_sets.keys()),
                    "present_in_datasets": ";".join(labels),
                }
            ),
            args.outdir / "shared_positive_host_gene_program.tsv",
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
