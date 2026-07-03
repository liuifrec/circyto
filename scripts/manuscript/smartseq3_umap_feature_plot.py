#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Planned Smart-seq3 UMAP/feature-plot entry point. "
            "This stub documents the expected CLI while manuscript plotting "
            "style and feature selection remain TODO."
        )
    )
    parser.add_argument("--input", required=True, type=Path, help="Smart-seq3 RNA+circ .h5mu file.")
    parser.add_argument("--outdir", required=True, type=Path, help="Output directory for plots and plot-ready tables.")
    parser.add_argument(
        "--feature",
        action="append",
        default=None,
        help="Feature/circRNA/gene to plot. Repeat for multiple features.",
    )
    parser.add_argument("--basis", default="X_umap", help="Embedding key to use when implemented.")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_arg_parser()
    parser.parse_args(argv)
    parser.error(
        "Smart-seq3 UMAP/feature plotting is planned but not implemented in this repository snapshot."
    )
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
