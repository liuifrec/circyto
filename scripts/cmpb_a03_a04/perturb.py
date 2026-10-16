"""Prespecified simulation inputs, independent of either converter or its oracle."""

import csv
import json
import math
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite

SEEDS = list(range(2026101500, 2026101520))
FRACTIONS = [0, 0.05, 0.10, 0.25, 0.50]


def plan():
    cases = [
        dict(case="BASE", scenario="CONTROL", variant="matrix-clean", valid=True),
        dict(
            case="NATIVE-CLEAN", scenario="CONTROL", variant="native-clean", valid=True
        ),
    ]
    for scenario, variants, valid in [
        ("R03", ["duplicate-metadata", "duplicate-RNA", "duplicate-native-map"], False),
        ("R04", ["literal-IDs"], True),
        ("R07", ["mixed-build", "mixed-coordinate", "mixed-support"], False),
        ("R08", ["duplicate-long", "duplicate-native"], False),
        ("R09", ["explicit-zero"], True),
        ("R10", ["no-denominator"], True),
        ("R11", ["missing-matrix-column", "missing-native-map"], False),
    ]:
        for v in variants:
            cases.append(
                dict(case=scenario + "-" + v, scenario=scenario, variant=v, valid=valid)
            )
    for scenario in ["R01", "R02", "R05", "R06"]:
        for model in (
            ["uniform"] if scenario == "R01" else ["uniform", "profile-associated"]
        ):
            for fraction in FRACTIONS:
                for seed in SEEDS:
                    cases.append(
                        dict(
                            case=f"{scenario}-{model}-{int(fraction * 100):02d}-{seed}",
                            scenario=scenario,
                            variant="grid",
                            valid=True,
                            model=model,
                            nominal_fraction=fraction,
                            seed=seed,
                        )
                    )
    return cases


def table(path, columns, rows):
    with Path(path).open("x", newline="") as f:
        w = csv.DictWriter(f, fieldnames=columns, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def choose(ids, weights, fraction, seed, model):
    k = math.floor(fraction * len(ids) + 0.5)
    p = (
        np.array([weights[c] for c in ids], dtype=float)
        if model == "profile-associated"
        else None
    )
    if p is not None:
        if not np.isfinite(p).all() or (p <= 0).any():
            raise ValueError("Invalid registered weights")
        p = p / p.sum()
    rng = np.random.default_rng(seed)
    return rng, rng.choice(len(ids), size=k, replace=False, p=p).tolist()


def tiny_native(directory, variant):
    header = [
        "circRNA_ID",
        "chr",
        "circRNA_start",
        "circRNA_end",
        "#junction_reads",
        "SM_MS_SMS",
        "#non_junction_reads",
        "junction_reads_ratio",
        "circRNA_type",
        "gene_id",
        "strand",
        "junction_reads_ID",
        "Score",
    ]
    rows = [
        dict(
            zip(
                header,
                [
                    "chr1:2|9",
                    "chr1",
                    "2",
                    "9",
                    "5",
                    "1",
                    "0",
                    "1",
                    "exon",
                    "SYNTHETIC_GENE",
                    "+",
                    "synthetic_read",
                    "1",
                ],
            )
        )
    ]
    if variant == "duplicate-native":
        rows = rows + rows
    table(directory / "a.tsv", header, rows)
    table(directory / "zero.tsv", header, [])
    mapping = [
        dict(
            cell_id="001",
            sample_id="S1",
            path="a.tsv",
            detector_version="3.0.1",
            result_scope="single-sample",
        ),
        dict(
            cell_id="NA",
            sample_id="S2",
            path="zero.tsv",
            detector_version="3.0.1",
            result_scope="single-sample",
        ),
    ]
    if variant == "duplicate-native-map":
        mapping[1]["cell_id"] = "001"
    if variant == "missing-native-map":
        mapping = mapping[:1]
    table(directory / "native.tsv", list(mapping[0]), mapping)
    table(
        directory / "cells.tsv",
        ["cell_id", "detector_status"],
        [dict(cell_id=c, detector_status="completed") for c in ["001", "NA"]],
    )
    rna = ad.AnnData(
        sparse.csr_matrix(np.array([[9, 1], [2, 7]], dtype=np.int64)),
        obs=pd.DataFrame(index=["NA", "001"]),
        var=pd.DataFrame(index=["synthetic_g1", "synthetic_g2"]),
    )
    rna.uns["fixture"] = "fully synthetic native-edge RNA"
    rna.write_h5ad(directory / "rna.h5ad", compression="gzip")
    return dict(
        mode="native",
        input=str(directory / "native.tsv"),
        cells=str(directory / "cells.tsv"),
        rna=str(directory / "rna.h5ad"),
        genome_build="synthetic",
        coordinate_convention="1-based-closed",
        support_definition="CIRI3_junction_reads",
        detector="CIRI3",
    ), dict(
        label="SYNTHETIC native edge case; not recovered calls",
        n_cells=2,
        selected_ids=[],
        moved_identities=0,
    )


def generate(base, case, directory):
    directory.mkdir(parents=True, exist_ok=False)
    variant = case["variant"]
    if "native" in variant:
        config, details = tiny_native(directory, variant)
    else:
        circ = ad.read_h5ad(base / "circ_source_copy.h5ad")
        rna = ad.read_h5ad(base / "rna_source_copy.h5ad")
        ids = list(circ.obs_names)
        features = circ.var.copy()
        features["circ_id"] = list(circ.var_names)
        matrix = circ.X.astype(np.int64).copy()
        status = {c: "completed" for c in ids}
        selected = []
        moved = 0
        removed_rna = []
        removed_circ = []
        if variant == "grid":
            weights = json.loads((base / "base_record.json").read_text())[
                "sampling_weights"
            ]
            rng, positions = choose(
                ids, weights, case["nominal_fraction"], case["seed"], case["model"]
            )
            selected = [ids[i] for i in positions]
            if case["scenario"] == "R01":
                order = list(rna.obs_names)
                slots = [order.index(c) for c in selected]
                replacement = rng.permutation(selected).tolist()
                for slot, value in zip(slots, replacement):
                    order[slot] = value
                moved = sum(a != b for a, b in zip(order, rna.obs_names))
                rna = rna[order].copy()
            elif case["scenario"] == "R02":
                removed_rna = selected[::2]
                removed_circ = selected[1::2]
                rna = rna[[c for c in rna.obs_names if c not in removed_rna]].copy()
                keep = [i for i, c in enumerate(ids) if c not in removed_circ]
                matrix = matrix[keep]
                ids = [ids[i] for i in keep]
            elif case["scenario"] == "R05":
                for i, c in enumerate(selected):
                    status[c] = ["failed", "unknown", "unprocessed", "skipped"][i % 4]
            elif case["scenario"] == "R06":
                for c in selected:
                    status[c] = "failed"
                    i = ids.index(c)
                    lo, hi = matrix.indptr[i : i + 2]
                    matrix.data[lo:hi] = np.maximum(1, matrix.data[lo:hi] // 2)
        if variant == "literal-IDs":
            rename = {ids[0]: "001", ids[1]: "NA"}
            ids = [rename.get(c, c) for c in ids]
            rna.obs_names = [rename.get(c, c) for c in rna.obs_names]
            status = {rename.get(c, c): v for c, v in status.items()}
        if variant == "no-denominator":
            status = {c: "unknown" for c in ids}
        if variant == "duplicate-RNA":
            new = list(rna.obs_names)
            new[1] = new[0]
            rna.obs_names = new
        metadata = [dict(cell_id=c, detector_status=status[c]) for c in ids]
        if variant == "duplicate-metadata":
            metadata[1]["cell_id"] = metadata[0]["cell_id"]
        declared = dict(
            genome_build="historical_IMR90_assembly_unverified",
            coordinate_convention="1-based-closed",
            support_definition="CIRI3_junction_reads",
            detector="CIRI3",
        )
        if variant.startswith("mixed-"):
            key = {
                "mixed-build": "genome_build",
                "mixed-coordinate": "coordinate_convention",
                "mixed-support": "support_definition",
            }[variant]
            for row in metadata:
                row[key] = declared[key]
            metadata[0][key] = {
                "mixed-build": "GRCh37",
                "mixed-coordinate": "0-based-half-open",
                "mixed-support": "UMI",
            }[variant]
        table(directory / "cells.tsv", list(metadata[0]), metadata)
        rna.write_h5ad(directory / "rna.h5ad", compression="gzip")
        folder = directory / "matrix"
        folder.mkdir()
        source_ids = ids.copy()
        if variant == "missing-matrix-column":
            source_ids = source_ids[1:]
            matrix = matrix[1:]
        coo = matrix.T.tocoo()
        if variant == "explicit-zero":
            present = set(zip(coo.row, coo.col))
            row, col = next(
                (i, j)
                for i in range(matrix.shape[1])
                for j in range(matrix.shape[0])
                if (i, j) not in present
            )
            coo = sparse.coo_matrix(
                (
                    np.append(coo.data, 0),
                    (np.append(coo.row, row), np.append(coo.col, col)),
                ),
                shape=coo.shape,
            )
        mmwrite(folder / "circ_counts.mtx", coo)
        (folder / "cell_index.txt").write_text("\n".join(source_ids) + "\n")
        (folder / "circ_index.txt").write_text("\n".join(circ.var_names) + "\n")
        features.to_csv(folder / "circ_feature_table.tsv", sep="\t", index=False)
        config = dict(
            mode="matrix",
            input=str(folder),
            cells=str(directory / "cells.tsv"),
            rna=str(directory / "rna.h5ad"),
            **declared,
        )
        if variant == "duplicate-long":
            r, c, v = next(
                zip(matrix.tocoo().row, matrix.tocoo().col, matrix.tocoo().data)
            )
            row = {
                k: str(features.iloc[c][k])
                for k in ["circ_id", "chrom", "start", "end", "strand"]
            }
            row.update(cell_id=ids[r], support=str(v))
            table(directory / "calls.tsv", list(row), [row, row])
            config.update(mode="calls", input=str(directory / "calls.tsv"))
        details = dict(
            label="DERIVED/SIMULATED from historical IMR90; invented starting completion",
            n_source_cells=23,
            n_circ_cells=len(ids),
            n_rna_cells=rna.n_obs,
            selected_ids=selected,
            selected_count=len(selected),
            realized_fraction=len(selected) / 23,
            moved_identities=moved,
            removed_RNA_membership=removed_rna,
            removed_circ_membership=removed_circ,
            status_rule=case["scenario"],
            historical_completion="unknown",
            full_cohort_positive_counts=np.asarray((circ.X > 0).sum(axis=0))
            .ravel()
            .astype(int)
            .tolist(),
            full_cohort_candidates=list(circ.var_names),
        )
    (directory / "config.json").write_text(
        json.dumps(config, indent=2, sort_keys=True) + "\n"
    )
    (directory / "perturbation.json").write_text(
        json.dumps(dict(**case, **details), indent=2, sort_keys=True) + "\n"
    )
    return config, details
