#!/usr/bin/env python3
"""Bounded A02 input preparation; original files are opened read-only."""

import argparse
import csv
import gzip
import hashlib
import json
import re
import shutil
import subprocess
from pathlib import Path

from run_stage import digest, dump, stamp

SOURCES = {
    "wheel": (
        "/tmp/circyto-distribution/reviewed-packages/circyto-0.10.0-py3-none-any.whl",
        "inputs/artifacts/circyto-0.10.0-py3-none-any.whl",
        "eca5f45f35b351ec44df4e5b0f73eb4e5e0564e2691631f8eb417b0e8a7ca03b",
    ),
    "ciri3_jar": (
        "/tmp/circyto-distribution/external/CIRI3_Java_18.0.1.jar",
        "external/CIRI3_Java_18.0.1.jar",
        "569428a4fe0d6573dbeeda6d3d4a0d457ac01d01bbf0b49ab1a42b84e24f1be0",
    ),
    "original_imr90": (
        "/tmp/circyto-import/load_work/scrr_imr90/full_length_rna_circ_cnv.hostgene_fixed.h5mu",
        "inputs/original/full_length_rna_circ_cnv.hostgene_fixed.h5mu",
        "bb2e12f7c3b36f9fa72d98cd71e8bea905a67f50e22af1d6b713550ee92b60c8",
    ),
}
COHORT = [
    ("IMR90_A_100", "GSM8558852", "SRR30918126"),
    ("IMR90_A_101", "GSM8558853", "SRR30918117"),
    ("IMR90_A_102", "GSM8558854", "SRR30918116"),
    ("IMR90_A_103", "GSM8558855", "SRR30918115"),
]
EXPECTED = dict(
    run_accession="SRR30918126",
    experiment_accession="SRX26321174",
    sample_accession="SAMN44108155",
    study_accession="PRJNA1169833",
    library_layout="SINGLE",
    read_count="3824354",
    base_count="577477454",
    fastq_ftp="ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/026/SRR30918126/SRR30918126.fastq.gz",
    fastq_bytes="225876952",
    fastq_md5="dfd377ce81eeeb23503f8e4756ef60ae",
)
REF_BASE = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/"
REF_NAMES = [
    "GRCh38.primary_assembly.genome.fa.gz",
    "gencode.v45.primary_assembly.annotation.gtf.gz",
]


def fetch(root, url, relative, required=True, max_seconds=300):
    dest = root / relative
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        raise FileExistsError(dest)
    argv = [
        "/usr/bin/curl",
        "--fail",
        "--location",
        "--show-error",
        "--silent",
        "--connect-timeout",
        "30",
        "--max-time",
        str(max_seconds),
        "--retry",
        "1",
        "--output",
        str(dest),
        url,
    ]
    result = subprocess.run(argv, check=False)
    record = {
        "argv": argv,
        "accessed_utc": stamp(),
        "exit_code": result.returncode,
        "path": relative,
        "sha256": digest(dest) if dest.is_file() else None,
        "bytes": dest.stat().st_size if dest.is_file() else None,
    }
    dump(root / "records" / ("acquire-" + dest.name + ".json"), record)
    if required and result.returncode:
        raise RuntimeError("Acquisition failed: " + url)
    return record


def initialize(root):
    records = {}
    for key, (source, target, expected) in SOURCES.items():
        if digest(source) != expected:
            raise ValueError("Source checksum mismatch: " + key)
        target = root / target
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists():
            raise FileExistsError(target)
        shutil.copyfile(source, target)
        if digest(target) != expected:
            raise ValueError("Copy checksum mismatch: " + key)
        records[key] = {
            "source": source,
            "copy": str(target.relative_to(root)),
            "bytes": target.stat().st_size,
            "sha256": expected,
        }
    dump(root / "records/input_artifacts.json", records)
    (root / "inputs/python-constraints.txt").write_text(
        "anndata==0.11.4\nmudata==0.3.10\nnumpy==2.2.6\nscipy==1.15.3\npandas==2.3.3\n"
    )
    with (root / "inputs/registered_cohort.tsv").open("x", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["cell_id", "gsm", "srr", "session_selection"])
        for index, row in enumerate(COHORT):
            writer.writerow(
                [*row, "first_cell_only" if index == 0 else "pending_separate_review"]
            )
    dump(
        root / "records/pilot_progress.json",
        {
            "schema": "circyto.cmpb.a02.progress.v1",
            "route": "NEW technical validation",
            "registration_head": "17b733c8f21bb128a035e7b100e704099e2bdb0e",
            "source_commit": "6c280426ebf10abcbcaeab72f87c9d34060d6ee8",
            "updated_utc": stamp(),
            "stage": "preflight",
            "four_cell_pilot_complete": False,
            "cells": [
                {
                    "cell_id": c,
                    "gsm": g,
                    "srr": s,
                    "execution_status": "not_started",
                    "original_detector_status": "unknown",
                    "new_detector_status": "unprocessed",
                    "acceptance": "pending",
                }
                for c, g, s in COHORT
            ],
        },
    )
    print(
        json.dumps(
            {
                "copied_hash_matched_artifacts": list(records),
                "selected_cell": COHORT[0][0],
            }
        )
    )


def metadata(root):
    url = (
        "https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=experiment_accession%3D%22SRX26321174%22"
        "&fields=" + ",".join(EXPECTED) + "&format=tsv&limit=20"
    )
    fetch(root, url, "inputs/metadata/ena-SRR30918126.tsv")
    with (root / "inputs/metadata/ena-SRR30918126.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 1:
        raise ValueError("ENA did not return exactly the registered run")
    differences = {
        key: {"registered": value, "live": rows[0].get(key)}
        for key, value in EXPECTED.items()
        if rows[0].get(key) != value
    }
    dump(
        root / "records/ena_comparison.json",
        {
            "expected": EXPECTED,
            "observed": rows[0],
            "differences": differences,
            "matches_registration": not differences,
            "layout_conflict": "Generic GEO prose says paired-end; run-specific ENA says SINGLE, one file. Acquired FASTQ still requires structural/layout validation.",
        },
    )
    if differences:
        raise ValueError("Live ENA differs from registration; review required")
    checksum = fetch(
        root,
        REF_BASE + "MD5SUMS",
        "inputs/metadata/gencode-v45-MD5SUMS",
        required=False,
    )
    evidence = {
        "provider": "GENCODE release 45",
        "source_url": REF_BASE + "MD5SUMS",
        "provider_checksum_response": checksum,
        "files": {},
    }
    text = (
        (root / "inputs/metadata/gencode-v45-MD5SUMS").read_text()
        if checksum["exit_code"] == 0
        else ""
    )
    for name in REF_NAMES:
        matches = [
            line
            for line in text.splitlines()
            if line.split() and line.split()[-1].lstrip("*./") == name
        ]
        evidence["files"][name] = {
            "url": REF_BASE + name,
            "provider_md5": matches[0].split()[0] if len(matches) == 1 else None,
            "checksum_status": "recorded" if len(matches) == 1 else "missing",
        }
    dump(root / "records/reference_provider_evidence.json", evidence)
    print(
        json.dumps(
            {"ena_matches": not differences, "provider_checksums": evidence["files"]}
        )
    )


def hashes(path):
    sha, md5 = hashlib.sha256(), hashlib.md5()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            sha.update(block)
            md5.update(block)
    return {
        "sha256": sha.hexdigest(),
        "md5": md5.hexdigest(),
        "bytes": Path(path).stat().st_size,
    }


def check_fastq(path, expected):
    identity = hashes(path)
    if identity["md5"] != expected["fastq_md5"] or identity["bytes"] != int(
        expected["fastq_bytes"]
    ):
        raise ValueError("FASTQ compressed identity differs from ENA")
    seen = set()
    names_hash = hashlib.sha256()
    reads = bases = 0
    lengths, mates = set(), set()
    mate_counts = {"1": 0, "2": 0, "unspecified": 0}
    first, last = [], []
    with gzip.open(path, "rt", encoding="ascii", newline="") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            seq, plus, quality = [handle.readline().rstrip("\r\n") for _ in range(3)]
            header = header.rstrip("\r\n")
            if (
                not header.startswith("@")
                or not plus.startswith("+")
                or not seq
                or len(seq) != len(quality)
            ):
                raise ValueError("Invalid FASTQ structure at record " + str(reads + 1))
            if any(ord(q) < 33 or ord(q) > 126 for q in quality):
                raise ValueError("Invalid quality encoding")
            identifier = header[1:].split()[0]
            if len(plus) > 1 and plus[1:].lstrip("@").split()[0] != identifier:
                raise ValueError("FASTQ repeated plus identifier disagrees")
            if identifier in seen:
                raise ValueError(
                    "Repeated read identifier / possible interleaved mates: "
                    + identifier
                )
            seen.add(identifier)
            if not identifier.startswith(expected["run_accession"] + "."):
                raise ValueError(
                    "Read identifier does not bind to fixed run: " + identifier
                )
            # Preserve and inspect explicit mate suffixes/Illumina second fields;
            # a run-level SINGLE declaration alone is not a FASTQ layout check.
            read_mates = set()
            if identifier.endswith(("/1", "/2")):
                read_mates.add(identifier[-1])
            fields = header.split()
            if len(fields) > 1:
                # A bare number here is often the SRA spot ordinal, not a mate.
                match = re.match(r"([12]):[YN]:", fields[1])
                if match is None:
                    match = re.search(r"/([12])$", fields[1])
                if match:
                    read_mates.add(match.group(1))
            mates.update(read_mates)
            for mate in read_mates or {"unspecified"}:
                mate_counts[mate] += 1
            names_hash.update(identifier.encode("ascii") + b"\n")
            reads += 1
            bases += len(seq)
            lengths.add(len(seq))
            if reads <= 3:
                first.append(header)
            last = (last + [header])[-3:]
    if "2" in mates:
        raise ValueError("Mate-2 identifier contradicts registered SINGLE layout")
    if reads != int(expected["read_count"]) or bases != int(expected["base_count"]):
        raise ValueError("Full FASTQ read/base totals differ from ENA")
    return {
        **identity,
        "gzip_crc_valid": True,
        "structure_valid": True,
        "reads": reads,
        "bases": bases,
        "read_lengths": sorted(lengths),
        "unique_read_identifiers": len(seen),
        "ordered_read_ids_sha256": names_hash.hexdigest(),
        "first_headers": first,
        "last_headers": last,
        "explicit_mates_observed": sorted(mates),
        "explicit_mate_header_counts": mate_counts,
        "layout": "single-end",
        "subsampling": False,
    }


def check_reference(fasta, gtf):
    lengths, contig, count = {}, None, 0
    with Path(fasta).open() as handle:
        for line in handle:
            if line.startswith(">"):
                if contig is not None:
                    lengths[contig] = count
                contig = line[1:].split()[0]
                if contig in lengths:
                    raise ValueError("Duplicate FASTA contig")
                count = 0
            else:
                if contig is None:
                    raise ValueError("Sequence before FASTA header")
                count += len(line.strip())
    if contig is not None:
        lengths[contig] = count
    if not lengths or any(n <= 0 for n in lengths.values()):
        raise ValueError("Empty FASTA/contig")
    annotations, headers, rows = set(), [], 0
    with Path(gtf).open() as handle:
        for line in handle:
            if line.startswith("#"):
                headers.append(line.rstrip())
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[0] not in lengths:
                raise ValueError("GTF contig/column mismatch")
            start, end = int(fields[3]), int(fields[4])
            if not 1 <= start <= end <= lengths[fields[0]]:
                raise ValueError("GTF coordinate outside FASTA")
            annotations.add(fields[0])
            rows += 1
    if not rows:
        raise ValueError("Empty annotation")
    return {
        "fasta_contig_lengths": lengths,
        "fasta_bases": sum(lengths.values()),
        "annotation_contigs": sorted(annotations),
        "annotation_records": rows,
        "annotation_headers": headers,
        "all_gtf_contigs_and_coordinates_compatible": True,
    }


def acquire(root):
    if not json.loads((root / "records/ena_comparison.json").read_text())[
        "matches_registration"
    ]:
        raise ValueError("Metadata gate not passed")
    progress = json.loads((root / "records/pilot_progress.json").read_text())
    progress.update(stage="acquiring_first_cell_inputs", updated_utc=stamp())
    progress["cells"][0]["execution_status"] = "preparing_inputs"
    dump(root / "records/pilot_progress.json", progress)
    fetch(
        root,
        "https://" + EXPECTED["fastq_ftp"],
        "inputs/fastq/SRR30918126.fastq.gz",
        max_seconds=3600,
    )
    evidence = check_fastq(root / "inputs/fastq/SRR30918126.fastq.gz", EXPECTED)
    dump(root / "records/fastq_validation.json", evidence)
    progress.update(stage="preparing_shared_reference", updated_utc=stamp())
    progress["cells"][0]["fastq_sha256"] = evidence["sha256"]
    dump(root / "records/pilot_progress.json", progress)
    provider = json.loads(
        (root / "records/reference_provider_evidence.json").read_text()
    )
    for name in REF_NAMES:
        fetch(root, REF_BASE + name, "ref/" + name, max_seconds=3600)
        compressed = root / "ref" / name
        identity = hashes(compressed)
        expected = provider["files"][name]["provider_md5"]
        if expected and identity["md5"] != expected:
            raise ValueError("Reference provider MD5 mismatch: " + name)
        unpacked = compressed.with_suffix("")
        with gzip.open(compressed, "rb") as source, unpacked.open("xb") as dest:
            shutil.copyfileobj(source, dest, length=1024 * 1024)
        dump(
            root / "records" / ("reference-" + name + ".json"),
            {
                "compressed": identity,
                "uncompressed": hashes(unpacked),
                "provider_md5": expected,
                "provider_md5_verified": bool(expected),
                "gzip_crc_valid": True,
                "genome_chain": "NEW validation only",
            },
        )
    compat = check_reference(
        root / "ref" / REF_NAMES[0][:-3], root / "ref" / REF_NAMES[1][:-3]
    )
    dump(root / "records/reference_compatibility.json", compat)
    with (root / "inputs/manifest_firstcell.tsv").open("x", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "cell_id",
                "platform",
                "read1",
                "read2",
                "bam",
                "library_id",
                "n_input_reads",
                "protocol",
                "strandedness",
                "read_layout",
            ]
        )
        writer.writerow(
            [
                "IMR90_A_100",
                "plate",
                str(root / "inputs/fastq/SRR30918126.fastq.gz"),
                "",
                "",
                "SRX26321174",
                evidence["reads"],
                "ramda",
                "unstranded",
                "single",
            ]
        )
    progress.update(stage="full_inputs_validated_index_pending", updated_utc=stamp())
    dump(root / "records/pilot_progress.json", progress)
    print(
        json.dumps(
            {
                "full_library_validated": evidence["reads"],
                "fasta_contigs": len(compat["fasta_contig_lengths"]),
                "gtf_records": compat["annotation_records"],
            }
        )
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("action", choices=["init", "metadata", "acquire", "recheck-fastq"])
    ap.add_argument("--root", type=Path, required=True)
    args = ap.parse_args()
    if args.action == "recheck-fastq":
        root = args.root.resolve()
        result = check_fastq(root / "inputs/fastq/SRR30918126.fastq.gz", EXPECTED)
        result["amendment"] = (
            "Recognize /1 or /2 on the original flow-cell identifier in the second header field; do not mistake bare SRA spot ordinals for mates. Initial validation retained separately."
        )
        dump(root / "records/fastq_validation_header_v2.json", result)
        print(
            json.dumps(
                {
                    "reads": result["reads"],
                    "mate_headers": result["explicit_mate_header_counts"],
                }
            )
        )
    else:
        {"init": initialize, "metadata": metadata, "acquire": acquire}[args.action](
            args.root.resolve()
        )


if __name__ == "__main__":
    main()
