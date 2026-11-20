"""Offline acceptance using original Git evidence, never simulated detector calls.

Run against the installed wheel outside the checkout:
  python validate_import_circ_evidence.py --repo /path/to/circyto --out /tmp/evidence

Processed-object matrix extraction is a round-trip test, NOT independent
validation of native CIRI3 files. All input recovery is read-only Git access.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from collections import Counter
from pathlib import Path

import anndata as ad
import mudata as mu
import numpy as np
import pandas as pd
from scipy.io import mmwrite

OBJECT_COMMIT = "c99cddaae2314d3af7d5232150c5b86bb6b52e25"
OBJECTS = {
    "smartseq3": (
        "load_work/emtab8735_smartseq3/full_length.hostgene_fixed.h5mu",
        "0ecd36bb0a74455db7f0affb9ade5023c1934c1dd234aca975365c0b69d8b339",
    ),
    "imr90": (
        "load_work/scrr_imr90/full_length_rna_circ_cnv.hostgene_fixed.h5mu",
        "bb2e12f7c3b36f9fa72d98cd71e8bea905a67f50e22af1d6b713550ee92b60c8",
    ),
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def recover(repo, revision, member, target):
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("wb") as handle:
        subprocess.run(
            ["git", "-C", str(repo), "show", f"{revision}:{member}"],
            stdout=handle,
            check=True,
        )
    return sha(target)


def cli(*args):
    subprocess.run(["circyto", *map(str, args)], check=True)


def run_import(matrix, cells, out, detector, convention, build, support, rna=None):
    args = [
        "import-circ",
        "--matrix-dir",
        matrix,
        "--cells",
        cells,
        "--outdir",
        out,
        "--genome-build",
        build,
        "--coordinate-convention",
        convention,
        "--support-definition",
        support,
        "--detector",
        detector,
    ]
    if rna is not None:
        args += ["--rna", rna]
    cli(*args)
    return ad.read_h5ad(out / "anndata/circ_counts.h5ad"), json.loads(
        (out / "qc/metrics.json").read_text()
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    results = {
        "scope": "Archived collector acceptance and checksum-verified processed-object round trips; no raw analysis",
        "native_CIRI3_and_CIRCexplorer_acceptance": "pending: original native call files not recovered",
        "processed_objects": {},
        "collectors": {},
    }
    for name, (member, digest) in OBJECTS.items():
        root = args.out / name
        source = root / "source.h5mu"
        assert recover(args.repo, OBJECT_COMMIT, member, source) == digest
        original = mu.read_h5mu(source)
        circ, rna = original["circ"], original["rna"]
        matrix = root / "processed_matrix_extraction"
        matrix.mkdir()
        mmwrite(matrix / "circ_counts.mtx", circ.X.T)
        (matrix / "circ_index.txt").write_text("\n".join(circ.var_names) + "\n")
        (matrix / "cell_index.txt").write_text("\n".join(circ.obs_names) + "\n")
        circ.var.to_csv(matrix / "circ_feature_table.tsv", sep="\t", index=False)
        obs = circ.obs.copy()
        # Keep legacy values as source evidence; do not promote aggregate success
        # summaries or padded rows to per-cell processing-status assertions.
        obs = obs.rename(
            columns={
                c: "source_" + c
                for c in ("circRNA_count", "circRNA_total_support", "membership")
                if c in obs
            }
        )
        obs.insert(0, "cell_id", circ.obs_names)
        obs["detector_status"] = "unknown"
        obs.to_csv(root / "cells.tsv", sep="\t", index=False)
        rna.write_h5ad(root / "rna.h5ad")
        restored, metrics = run_import(
            matrix,
            root / "cells.tsv",
            root / "imported",
            "CIRI3",
            "1-based-closed",
            "GRCh38",
            "archived_CIRI3_support",
            root / "rna.h5ad",
        )
        assert restored.obs_names.tolist() == circ.obs_names.tolist()
        assert restored.var_names.tolist() == circ.var_names.tolist()
        assert (restored.X != circ.X).nnz == 0
        for col in circ.var:
            assert (
                restored.var[col].astype(str).tolist()
                == circ.var[col].astype(str).tolist()
            ), col
        for col in obs.columns.drop("cell_id"):
            assert (
                restored.obs[col].astype(str).tolist() == obs[col].astype(str).tolist()
            ), col
        paired = mu.read_h5mu(root / "imported/mudata/circyto_multimodal.h5mu")
        assert paired["rna"].obs_names.tolist() == rna.obs_names.tolist()
        assert (paired["rna"].X != rna.X).nnz == 0
        pd.testing.assert_frame_equal(paired["rna"].var, rna.var)
        pd.testing.assert_frame_equal(paired["rna"].obs, rna.obs)
        assert paired["circ"].obs_names.tolist() == circ.obs_names.tolist()
        assert (paired["circ"].X != circ.X).nnz == 0
        assert metrics["detection_fraction_completed"]["value"] is None
        assert metrics["rna_overlap"]["n_shared"] == len(
            set(rna.obs_names) & set(circ.obs_names)
        )
        assert sha(source) == digest
        result = {
            "source_commit": OBJECT_COMMIT,
            "source_member": member,
            "source_sha256": digest,
            "n_cells": circ.n_obs,
            "n_candidates": circ.n_vars,
            "nnz": int(restored.X.nnz),
            "total_support": int(restored.X.sum()),
            "host_gene_annotated": int((circ.var.host_gene.astype(str) != "").sum()),
            "matched_cell_overlap": metrics["rna_overlap"]["n_shared"],
            "median_candidates_per_cell": float(np.median(restored.X.getnnz(axis=1))),
            "median_support_per_cell": float(
                np.median(np.asarray(restored.X.sum(axis=1)).ravel())
            ),
            "observation_status": "unknown; no per-cell source status asserted",
            "declaration_scope": "GRCh38 from archived workflow reference paths; CIRI 1-based convention declared, no reference revalidation",
            "source_annotation_columns_preserved": list(circ.var.columns),
        }
        if name == "smartseq3":
            candidate = "chr1:117402186|117420649"
            x = restored[:, candidate].X
            result["MAN1A2"] = {
                "candidate": candidate,
                "detecting_cells": int(x.getnnz()),
                "support": int(x.sum()),
                "host_gene": restored.var.loc[candidate, "host_gene"],
            }
            assert (
                result["MAN1A2"]["detecting_cells"] == 6
                and result["MAN1A2"]["support"] == 15
            )
            assert (
                circ.n_obs,
                circ.n_vars,
                restored.X.nnz,
                result["host_gene_annotated"],
            ) == (192, 2503, 2659, 2379)
            assert (
                restored.obs.loc["DIYHEK_192", "circ_observation_status"] == "unknown"
            )
            assert restored.obs.loc["DIYHEK_192", "source_membership"] == "rna_only"
            result["DIYHEK_192"] = (
                "source RNA-only padded circ row preserved; unknown status, not a completed zero"
            )
        else:
            assert (circ.n_obs, circ.n_vars, result["host_gene_annotated"]) == (
                23,
                2443,
                2429,
            )
        results["processed_objects"][name] = result

    # CIRI2: committed normalized detector output from a public accession run.
    root = args.out / "ciri2"
    source = root / "source/ERR1041421.tsv"
    member = "tests/ciri2_err1041421_integration/ERR1041421.tsv"
    digest = recover(
        args.repo, "a2c934ec1a629528aae69dfb61df10af0f1407cd", member, source
    )
    cli(
        "collect-matrix",
        "-d",
        "ciri2",
        "--indir",
        source.parent,
        "--outdir",
        root / "matrix",
        "--min-count-per-cell",
        "0",
    )
    (root / "cells.tsv").write_text("cell_id\tdetector_status\nERR1041421\tunknown\n")
    a, m = run_import(
        root / "matrix",
        root / "cells.tsv",
        root / "imported",
        "CIRI2",
        "1-based-closed",
        "historical_chr21_assembly_unverified",
        "BSJ_reads",
    )
    expected = pd.read_csv(source, sep="\t")
    assert a.var_names.tolist() == expected.circ_id.tolist()
    assert a.X.toarray().ravel().tolist() == expected.support.tolist()
    for output, source_col in [
        ("chrom", "chr"),
        ("start", "start"),
        ("end", "end"),
        ("strand", "strand"),
    ]:
        assert a.var[output].tolist() == expected[source_col].tolist()
    results["collectors"]["ciri2"] = {
        "source_commit": "a2c934ec1a629528aae69dfb61df10af0f1407cd",
        "source_member": member,
        "source_sha256": digest,
        "n_candidates": a.n_vars,
        "total_support": m["total_support"],
        "scope": "normalized historical detector output, not native CIRI2 parser validation",
    }

    # find-circ: original large archived output, not the subsequently edited fixture.
    root = args.out / "find_circ"
    source = root / "source/ERR2139486/ERR2139486_splice_sites.bed"
    member = (
        "tests/find_circ3_manifest2_integration/ERR2139486/ERR2139486_splice_sites.bed"
    )
    revision = "bfb52c5499f4c64ac0a360e180d7acd6bef29076"
    digest = recover(args.repo, revision, member, source)
    cli(
        "collect-matrix",
        "-d",
        "find-circ3",
        "--indir",
        source.parent.parent,
        "--outdir",
        root / "matrix",
        "--min-count-per-cell",
        "0",
    )
    # Independently apply the documented existing collector rule to original rows.
    expected = Counter()
    features = {}
    raw_rows = 0
    for line in source.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        raw_rows += 1
        tags = set(fields[17].split(","))
        if (
            {"CIRCULAR", "UNAMBIGUOUS_BP", "CANONICAL"} <= tags
            and {"PERFECT_EXT", "GOOD_EXT"} & tags
            and "NO_UNIQ_BRIDGES" not in tags
            and int(fields[7]) >= 1
        ):
            chrom, start, end, strand = (
                fields[0],
                int(fields[1]),
                int(fields[2]),
                fields[5],
            )
            candidate = f"{chrom}:{start}|{end}|{strand}"
            expected[candidate] += int(fields[7])
            features[candidate] = [candidate, chrom, start, end, strand]
    pd.DataFrame(
        list(features.values()), columns=["circ_id", "chrom", "start", "end", "strand"]
    ).to_csv(root / "matrix/circ_feature_table.tsv", sep="\t", index=False)
    (root / "cells.tsv").write_text("cell_id\tdetector_status\nERR2139486\tunknown\n")
    a, m = run_import(
        root / "matrix",
        root / "cells.tsv",
        root / "imported",
        "find-circ3",
        "0-based-half-open",
        "historical_chr21_assembly_unverified",
        "unique_bridges",
    )
    assert set(a.var_names) == set(expected)
    assert a.X.toarray().ravel().tolist() == [expected[x] for x in a.var_names]
    results["collectors"]["find_circ"] = {
        "source_commit": revision,
        "source_member": member,
        "source_sha256": digest,
        "raw_rows": raw_rows,
        "n_candidates": a.n_vars,
        "total_support": m["total_support"],
        "scope": "original archived splice-sites output; existing flag-filter and unique-bridge semantics",
    }
    (args.out / "validation.json").write_text(
        json.dumps(results, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
