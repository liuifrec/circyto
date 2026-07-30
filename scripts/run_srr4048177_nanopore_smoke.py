#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

from circyto.manifest.long_read import LongReadManifestRow, write_long_read_manifest_tsv
from circyto.paths import get_package_root
from circyto.pipeline.nanopore_archive import (
    download_ena_fastq,
    load_expected_run,
    query_ena_run,
    validate_ena_run_identity,
)
from circyto.pipeline.nanopore_interop import prepare_nanopore_alignments
from circyto.pipeline.workflow_reporting import utc_now_iso, write_json


ACCESSION = "SRR4048177"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Opt-in real-data smoke test for experimental Mandalorion Nanopore "
            "interoperability. This is not a circRNA validation workflow."
        )
    )
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--reference", type=Path)
    parser.add_argument("--reference-id")
    parser.add_argument("--reference-build")
    parser.add_argument("--reference-sha256")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--minimap2", default="minimap2")
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--timeout-seconds", type=float, default=60.0)
    parser.add_argument("--metadata-only", action="store_true")
    parser.add_argument("--download-only", action="store_true")
    parser.add_argument("--keep-sam", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.threads < 1:
        raise SystemExit("--threads must be >= 1")
    args.outdir.mkdir(parents=True, exist_ok=True)
    expectation_path = (
        get_package_root() / "resources" / "nanopore" / "srr4048177_expected.json"
    )
    expectation = load_expected_run(expectation_path)
    metadata = query_ena_run(ACCESSION, timeout_seconds=args.timeout_seconds)
    warnings = validate_ena_run_identity(metadata, expectation)
    metadata_path = args.outdir / "archive_metadata.json"
    metadata_payload = {
        "schema_version": "circyto.ena_metadata_snapshot.v1",
        "metadata_retrieved_at": utc_now_iso(),
        "hard_identity_check": "passed",
        "warning_level_comparisons": warnings,
        "run": metadata.to_dict(),
    }
    write_json(metadata_path, metadata_payload)
    if args.metadata_only:
        print(json.dumps(metadata_payload, indent=2, sort_keys=True))
        return 0

    if len(metadata.fastq_files) != 1:
        raise SystemExit(
            f"{ACCESSION} must resolve to exactly one FASTQ after identity validation; "
            f"observed {len(metadata.fastq_files)}"
        )
    remote = metadata.fastq_files[0]
    fastq_path = args.outdir / "downloads" / remote.filename
    download_summary = download_ena_fastq(
        remote,
        fastq_path,
        timeout_seconds=args.timeout_seconds,
    )

    defaults = expectation["manifest_defaults"]
    manifest_path = args.outdir / "long_read_manifest.tsv"
    write_long_read_manifest_tsv(
        [
            LongReadManifestRow(
                cell_id=defaults["cell_id"],
                long_read_fastq=str(fastq_path.resolve()),
                protocol=defaults["protocol"],
                sequencing_platform=metadata.instrument_platform,
                archive_library_selection=metadata.archive_library_selection,
                library_preparation_summary=defaults["library_preparation_summary"],
                molecule_type=defaults["molecule_type"],
                barcode_status=defaults["barcode_status"],
                source_accession=metadata.run_accession,
                dataset_id=metadata.study_accession,
                biological_interpretation_boundary=defaults[
                    "biological_interpretation_boundary"
                ],
            )
        ],
        manifest_path,
    )
    if args.download_only:
        print(
            json.dumps(
                {
                    "archive_metadata": str(metadata_path.resolve()),
                    "download": download_summary,
                    "manifest": str(manifest_path.resolve()),
                    "alignment_run": False,
                },
                indent=2,
                sort_keys=True,
            )
        )
        return 0

    missing_reference_args = [
        name
        for name, value in (
            ("--reference", args.reference),
            ("--reference-id", args.reference_id),
            ("--reference-build", args.reference_build),
            ("--reference-sha256", args.reference_sha256),
        )
        if not value
    ]
    if missing_reference_args:
        raise SystemExit(
            "Alignment requires explicit reference identity: "
            + ", ".join(missing_reference_args)
        )
    alignment_manifest = prepare_nanopore_alignments(
        manifest_path=manifest_path,
        reference_fasta=args.reference,
        reference_id=args.reference_id,
        reference_build=args.reference_build,
        reference_sha256=args.reference_sha256,
        outdir=args.outdir,
        threads=args.threads,
        minimap2=args.minimap2,
        samtools=args.samtools,
        keep_sam=args.keep_sam,
        archive_metadata_path=metadata_path,
    )
    print(
        json.dumps(
            {
                "archive_metadata": str(metadata_path.resolve()),
                "download": download_summary,
                "manifest": str(manifest_path.resolve()),
                "alignment_manifest": str(alignment_manifest.resolve()),
                "circRNA_validation_status": False,
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
