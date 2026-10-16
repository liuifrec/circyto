#!/usr/bin/env python3
"""Bind the remaining-three real executions and prepare joint postprocessing."""

import argparse
import copy
import csv
from datetime import datetime
import json
from pathlib import Path
import shutil

import anndata as ad
import h5py

from continuation_inputs import ROWS, SOURCE
from execution_records import alignment_output_binding
from independent_checks import SOURCE_SHA, equal_adata, native
from run_stage import digest, dump, stamp

CELLS = [f"IMR90_A_{i}" for i in range(100, 104)]


def freeze_recheck(root, original, revision="v2"):
    directories = (
        ["acceptance"] if revision == "v2" else ["acceptance", "acceptance_v2"]
    )
    failures = {
        directory: json.loads(
            (root / "records" / directory / "acceptance_failure.json").read_text()
        )
        for directory in directories
    }
    if "1.000000" not in failures["acceptance"]["error"]:
        raise ValueError("Expected preserved padding-assumption failure")
    if (
        revision == "v3"
        and "FileNotFoundError" not in failures["acceptance_v2"]["error"]
    ):
        raise ValueError("Expected preserved RNA-copy assumption failure")
    record = root / f"records/continuation_helpers_{revision}.json"
    if record.exists():
        raise FileExistsError(record)
    dump(
        record,
        {
            "frozen_utc": stamp(),
            "sha256": {p.name: digest(p) for p in Path(__file__).parent.glob("*.py")},
            "retained_failures": failures,
            "retained_audit_files": {
                str(p.relative_to(root)): digest(p)
                for directory in directories
                for p in (root / "records" / directory).rglob("*")
                if p.is_file()
            },
            "original_checker_stage": "records/continuation-independent-acceptance.json",
            "new_audit_directory": f"records/acceptance_{revision}",
            "correction": "Exact rational equality for zero padding; RNA hash checked against command-bound original input; no numerical tolerance or changed scientific assertion",
        },
    )
    observation = Path("/tmp/continuation-memory-observation.json")
    if (
        observation.exists()
        and not (root / "records/memory_ceiling_observation.json").exists()
    ):
        target = root / "records/memory_ceiling_observation.json"
        if target.exists():
            raise FileExistsError(target)
        shutil.copyfile(observation, target)
    print(
        json.dumps(
            {
                "original_failure_retained": True,
                "new_audit_directory": f"records/acceptance_{revision}",
            }
        )
    )


def start(root, original):
    if (root / "workflow_remaining3").exists():
        raise FileExistsError("Remaining-three workflow must have a new directory")
    with (root / "inputs/manifest_remaining3.tsv").open(newline="") as h:
        rows = list(csv.DictReader(h, delimiter="\t"))
    if [r["cell_id"] for r in rows] != CELLS[1:]:
        raise ValueError("Manifest scope differs from the fixed remaining three")
    repo = Path(__file__).resolve().parents[2]
    reviewed = json.loads(
        (
            repo
            / "manuscript/cmpb/a02_firstcell/status_revalidation/output_artifacts.json"
        ).read_text()
    )
    final = original / "status_revalidation_20261016/final_candidate"
    for path, record in reviewed["files"].items():
        source = final / path
        if (
            digest(source) != record["sha256"]
            or source.stat().st_size != record["bytes"]
        ):
            raise ValueError("Accepted A100 derived artifact differs: " + path)
    dump(
        root / "records/repaired_A100_artifacts_verified.json",
        {
            "reviewed_manifest_sha256": digest(
                repo
                / "manuscript/cmpb/a02_firstcell/status_revalidation/output_artifacts.json"
            ),
            "files": reviewed["files"],
            "all_match": True,
            "checked_utc": stamp(),
        },
    )
    dump(
        root / "records/continuation_helpers.json",
        {
            "frozen_utc": stamp(),
            "sha256": {p.name: digest(p) for p in Path(__file__).parent.glob("*.py")},
        },
    )
    dump(
        root / "records/pilot_progress.json",
        {
            "updated_utc": stamp(),
            "stage": "remaining_three_workflow_requested",
            "four_cell_acceptance": "pending",
            "cells": [
                {
                    "cell_id": c,
                    "execution_status": "preserved_original"
                    if i == 0
                    else "not_started",
                    "detector_status": "completed" if i == 0 else "unprocessed",
                    "acceptance": "status_revalidation_passed" if i == 0 else "pending",
                    "scientific_source": "6c280426ebf10abcbcaeab72f87c9d34060d6ee8"
                    if i == 0
                    else SOURCE,
                    "derived_status_source": SOURCE if i == 0 else None,
                    "historical_detector_status": "unknown",
                }
                for i, c in enumerate(CELLS)
            ],
        },
    )
    print(json.dumps({"new_workflow_scope": CELLS[1:], "A100_rerun": False}))


def four_rna(original):
    source = original / "inputs/original/full_length_rna_circ_cnv.hostgene_fixed.h5mu"
    if digest(source) != SOURCE_SHA:
        raise ValueError("Original H5MU identity changed")
    with h5py.File(source, "r") as h:
        rna = ad.io.read_elem(h["mod/rna"])
    ids = [str(c) for c in rna.obs_names if c in CELLS]
    if ids != CELLS[::-1] or rna.n_vars != 63187 or not rna.obs_names.is_unique:
        raise ValueError("Original RNA identities/order/features disagree")
    for i, cell in enumerate(CELLS):
        for k, v in dict(
            gsm_id=f"GSM{8558852 + i}",
            canonical_cell_id=cell,
            rna_cell_id="RNA_" + cell,
            dna_cell_id="DNA_" + cell,
        ).items():
            if str(rna.obs.loc[cell, k]) != v:
                raise ValueError("Same-cell RNA identity mismatch: " + cell + "/" + k)
    out = rna[ids, :].copy()
    out.obs = rna.obs.loc[ids].copy()
    out.var = rna.var.copy()
    out.uns = copy.deepcopy(rna.uns)
    return out


def native_union(paths):
    """Independent first-seen union, rejecting incompatible shared annotations."""
    per_cell, union = {}, {}
    fields = [
        "chr",
        "circRNA_start",
        "circRNA_end",
        "strand",
        "gene_id",
        "circRNA_type",
    ]
    for cell, path in paths.items():
        records, _ = native(path)
        per_cell[cell] = records
        for r in records:
            key = r["circRNA_ID"]
            if key in union:
                for field in fields:
                    if r[field] != union[key][field]:
                        raise ValueError(
                            f"Native cross-cell conflict: {cell}/{key}/{field}: {union[key][field]!r} vs {r[field]!r}"
                        )
            else:
                union[key] = r
    return per_cell, union


def bind(root, original):
    env = json.loads((root / "records/preflight.json").read_text())
    records = [
        (p, json.loads(p.read_text()))
        for p in (root / "records/tool_commands").glob("*.json")
    ]
    align = [(p, r) for p, r in records if r["kind"] == "alignment"]
    detect = [(p, r) for p, r in records if r["kind"] == "detector"]
    if len(align) != 3 or len(detect) != 3:
        raise ValueError(
            "Exactly three alignment and detector outcomes required; retain partial evidence"
        )
    intervals = sorted(
        (
            datetime.fromisoformat(r["started_utc"]),
            datetime.fromisoformat(r["finished_utc"]),
        )
        for _, r in align + detect
    )
    if any(
        previous[1] > current[0] for previous, current in zip(intervals, intervals[1:])
    ):
        raise ValueError(
            "Biological-cell alignment/detection command intervals overlap"
        )
    work = root / "workflow_remaining3"
    child = json.loads(
        (work / "align/bwa_mem/alignment_prepare_summary.json").read_text()
    )
    combined = json.loads((work / "align/alignment_prepare_summary.json").read_text())
    detector = json.loads((work / "ciri3/detector_run_summary.json").read_text())
    ids = CELLS[1:]
    for summary in [child, combined, detector]:
        if [r["cell_id"] for r in summary["cells"]] != ids:
            raise ValueError(
                "Real child/combined/detector cell identities/order disagree"
            )
    if combined["cells"] != child["cells"] or combined["status_counts"] != {
        "aligned": 3
    }:
        raise ValueError(
            "Normal workflow status handoff disagrees with real child records"
        )
    bindings = []
    fa = original / "ref/GRCh38.primary_assembly.genome.fa"
    gtf = original / "ref/gencode.v45.primary_assembly.annotation.gtf"
    jar = original / "external/CIRI3_Java_18.0.1.jar"
    for i, row in enumerate(ROWS):
        cell = ids[i]
        fq = json.loads((root / f"records/fastq-{row[2]}.json").read_text())
        raw = work / "ciri3" / (cell + ".ciri3_run/ciri3_raw.tsv")
        aa = [(p, r) for p, r in align if r["argv"][-1] == fq["path"]]
        dd = [(p, r) for p, r in detect if r["output_path"] == str(raw)]
        if len(aa) != 1 or len(dd) != 1:
            raise ValueError("Ambiguous/missing cell-to-command identity: " + cell)
        ap, a = aa[0]
        dp, d = dd[0]
        if (
            d.get("exit_code") != 0
            or not d.get("completion_bound")
            or not d.get("inputs_unchanged")
        ):
            raise ValueError("Detector completion is unbound: " + cell)
        linked = alignment_output_binding(
            root,
            a,
            d,
            workflow_name="workflow_remaining3",
            cell_id=cell,
            group_id=row[3],
            reference=fa,
            cohort_ids=ids,
        )
        bwa = env["external_executables"]["bwa"]
        java = env["external_executables"]["java"]
        if a["argv"] != [
            bwa["path"],
            "mem",
            "-t",
            "8",
            "-k",
            "15",
            "-T",
            "15",
            str(fa),
            fq["path"],
        ]:
            raise ValueError("Actual BWA argv differs: " + cell)
        if d["argv"] != [
            java["path"],
            "-jar",
            str(jar),
            "-I",
            linked["output_path"],
            "-O",
            str(raw),
            "-F",
            str(fa),
            "-T",
            "8",
            "-G",
            str(raw.parent / "ciri3_internal.log"),
            "-A",
            str(gtf),
            "-S",
            "0",
            "-Ma",
            "0",
        ]:
            raise ValueError("Actual detector argv differs: " + cell)
        for actual, expected in [
            (a["executable_sha256"], bwa["sha256"]),
            (d["executable_sha256"], java["sha256"]),
            (a["input_sha256"][fq["path"]], fq["sha256"]),
            (d["input_sha256"][str(jar)], env["external_jar_sha256"]),
            (digest(raw), d["output_sha256"]),
        ]:
            if actual != expected:
                raise ValueError("Execution identity/hash link differs: " + cell)
        for rec in [a, d]:
            for path, sha in rec["input_sha256"].items():
                if digest(path) != sha:
                    raise ValueError("Execution input changed: " + path)
        if (
            child["cells"][i]["status"] != "aligned"
            or detector["cells"][i]["status"] != "success"
        ):
            raise ValueError("Child completion contradicts bound exits: " + cell)
        bindings.append(
            dict(
                cell_id=cell,
                gsm=row[1],
                run=row[2],
                alignment_status="aligned",
                detector_status="completed",
                source_commit=SOURCE,
                native_path=str(raw),
                native_sha256=digest(raw),
                alignment_record=str(ap),
                alignment_record_sha256=digest(ap),
                detector_record=str(dp),
                detector_record_sha256=digest(dp),
                emitted_bwa_argv=a["argv"],
                emitted_detector_argv=d["argv"],
                recorded_environment=d["environment"],
                executable_sha256={
                    "bwa": a["executable_sha256"],
                    "java": d["executable_sha256"],
                },
                input_sha256={**a["input_sha256"], **d["input_sha256"]},
                alignment_output_binding=linked,
                original_detector_status="unknown",
            )
        )
        dump(root / "records/native_completion_bindings.json", bindings)
    print(
        json.dumps(
            {
                "bound_new_cells": len(bindings),
                "normal_combined_status": combined["status_counts"],
            }
        )
    )


def joint_inputs(root, original):
    bindings = json.loads(
        (root / "records/native_completion_bindings.json").read_text()
    )
    if [b["cell_id"] for b in bindings] != CELLS[1:]:
        raise ValueError("Remaining three are not all bound")
    old = json.loads((original / "records/native_completion_binding.json").read_text())
    from replay_alignment_status import bound_evidence

    bound_evidence(original)  # Read-only reuse of the accepted original binding checks.
    repaired = json.loads(
        (
            original
            / "status_revalidation_20261016/final_candidate/workflow_firstcell/align/alignment_prepare_summary.json"
        ).read_text()
    )
    if [r["cell_id"] for r in repaired["cells"]] != [CELLS[0]]:
        raise ValueError("Repaired A100 alignment identity differs")
    old_alignment = repaired["cells"][0]["status"]
    if old_alignment != "aligned" or old["new_detector_status"] != "completed":
        raise ValueError("A100 accepted completion evidence differs")
    raw = original / "workflow_firstcell/ciri3/IMR90_A_100.ciri3_run/ciri3_raw.tsv"
    if digest(raw) != old["native_sha256"]:
        raise ValueError("Original A100 native changed")
    bindings.insert(
        0,
        dict(
            cell_id=CELLS[0],
            gsm="GSM8558852",
            run="SRR30918126",
            native_path=str(raw),
            native_sha256=old["native_sha256"],
            detector_status=old["new_detector_status"],
            alignment_status=old_alignment,
            source_commit="6c280426ebf10abcbcaeab72f87c9d34060d6ee8",
            derived_status_source=SOURCE,
            binding_record=str(original / "records/native_completion_binding.json"),
            original_detector_status="unknown",
        ),
    )
    paths = {b["cell_id"]: Path(b["native_path"]) for b in bindings}
    per_cell, union = native_union(paths)
    dump(
        root / "records/joint_native_sources.json",
        dict(
            cells=bindings,
            union_candidates=len(union),
            per_cell_candidates={c: len(r) for c, r in per_cell.items()},
            per_cell_support={
                c: sum(int(x["#junction_reads"]) for x in r)
                for c, r in per_cell.items()
            },
            mixed_execution_lineage=True,
        ),
    )
    for filename, header, rows in [
        (
            "native_four.tsv",
            ["cell_id", "sample_id", "path", "detector_version", "result_scope"],
            [
                [b["cell_id"], b["gsm"], b["native_path"], "3.0.1", "single-sample"]
                for b in bindings
            ],
        ),
        (
            "cells_four.tsv",
            ["cell_id", "detector_status"],
            [[b["cell_id"], b["detector_status"]] for b in bindings],
        ),
    ]:
        with (root / "inputs" / filename).open("x", newline="") as h:
            writer = csv.writer(h, delimiter="\t", lineterminator="\n")
            writer.writerow(header)
            writer.writerows(rows)
    rna = four_rna(original)
    target = root / "inputs/rna/imr90_original_overlap_profile.four.h5ad"
    target.parent.mkdir(parents=True)
    if target.exists():
        raise FileExistsError(target)
    rna.write_h5ad(target)
    equal_adata(rna, ad.read_h5ad(target), "four-cell RNA extraction")
    dump(
        root / "records/rna_extraction.json",
        dict(
            source_sha256=SOURCE_SHA,
            shape=list(rna.shape),
            cells=list(rna.obs_names),
            output=str(target),
            sha256=digest(target),
            all_available_fields_equal=True,
            historical_reference_hash="unknown",
            profile="historical same-cell lightweight simple-overlap",
            metadata_policy="Original RNA declarations unchanged; new provenance separate; parent uns retained in original source",
        ),
    )
    print(
        json.dumps(
            {
                "native_union": len(union),
                "rna_shape": list(rna.shape),
                "rna_order": list(rna.obs_names),
            }
        )
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "action", choices=["start", "bind", "joint-inputs", "freeze-recheck"]
    )
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--original", type=Path, required=True)
    ap.add_argument("--audit-revision", choices=["v2", "v3"], default="v2")
    args = ap.parse_args()
    if args.action == "freeze-recheck":
        freeze_recheck(
            args.root.resolve(), args.original.resolve(), args.audit_revision
        )
        return
    {
        "start": start,
        "bind": bind,
        "joint-inputs": joint_inputs,
        "freeze-recheck": freeze_recheck,
    }[args.action](args.root.resolve(), args.original.resolve())


if __name__ == "__main__":
    main()
