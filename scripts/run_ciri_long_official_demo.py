#!/usr/bin/env python3
"""Opt-in CIRI-long official-demo smoke runner.

The demo is fetched from its official GitHub release on demand and is never
redistributed with circyto. Default tests exercise only ``--dry-run``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import tarfile
import urllib.request
from pathlib import Path
from typing import Any

from circyto.manifest.ciri_long import (
    CiriLongManifestRow,
    write_ciri_long_manifest_tsv,
)
from circyto.pipeline.ciri_long_adapter import (
    normalize_ciri_long_outputs,
    run_ciri_long_call_stage,
    run_ciri_long_collapse_stage,
)
from circyto.pipeline.workflow_reporting import utc_now_iso


RELEASE_TAG = "v0.6-alpha"
ASSET_NAME = "CIRI-long_test_data.tar.gz"
RELEASE_API = (
    "https://api.github.com/repos/bioinfo-biols/CIRI-long/releases/tags/"
    + RELEASE_TAG
)
OFFICIAL_SOURCE = (
    "https://github.com/bioinfo-biols/CIRI-long/releases/tag/" + RELEASE_TAG
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def query_demo_metadata(*, timeout_seconds: float) -> dict[str, Any]:
    request = urllib.request.Request(
        RELEASE_API,
        headers={
            "Accept": "application/vnd.github+json",
            "User-Agent": "circyto-ciri-long-demo",
        },
    )
    try:
        with urllib.request.urlopen(request, timeout=timeout_seconds) as response:
            payload = json.loads(response.read().decode("utf-8"))
    except Exception as exc:
        raise RuntimeError(f"Unable to query official CIRI-long release: {exc}") from exc
    if payload.get("tag_name") != RELEASE_TAG:
        raise RuntimeError(
            f"Official release identity mismatch: expected tag {RELEASE_TAG!r}, "
            f"got {payload.get('tag_name')!r}"
        )
    matches = [
        asset for asset in payload.get("assets", []) if asset.get("name") == ASSET_NAME
    ]
    if len(matches) != 1:
        raise RuntimeError(
            f"Official release must contain exactly one {ASSET_NAME!r} asset; "
            f"found {len(matches)}"
        )
    asset = matches[0]
    size = asset.get("size")
    url = asset.get("browser_download_url")
    if not isinstance(size, int) or size <= 0:
        raise RuntimeError(f"Official demo asset has invalid size metadata: {size!r}")
    if not isinstance(url, str) or not url.startswith("https://github.com/"):
        raise RuntimeError(f"Official demo asset has invalid download URL: {url!r}")
    return {
        "schema_version": "circyto.ciri_long_demo_metadata.v1",
        "retrieved_at": utc_now_iso(),
        "release_tag": RELEASE_TAG,
        "asset_name": ASSET_NAME,
        "asset_id": asset.get("id"),
        "compressed_bytes": size,
        "download_url": url,
        "digest": asset.get("digest"),
        "official_source": OFFICIAL_SOURCE,
        "redistribution": "not_permitted_by_this_script",
    }


def _download(
    *,
    metadata: dict[str, Any],
    destination: Path,
    timeout_seconds: float,
) -> dict[str, Any]:
    destination.parent.mkdir(parents=True, exist_ok=True)
    part = Path(f"{destination}.part")
    if destination.exists():
        observed_size = destination.stat().st_size
        if observed_size != int(metadata["compressed_bytes"]):
            raise RuntimeError(
                f"Existing demo archive size mismatch: {observed_size} != "
                f"{metadata['compressed_bytes']}"
            )
    else:
        request = urllib.request.Request(
            str(metadata["download_url"]),
            headers={"User-Agent": "circyto-ciri-long-demo"},
        )
        try:
            with urllib.request.urlopen(request, timeout=timeout_seconds) as response, part.open(
                "wb"
            ) as handle:
                shutil.copyfileobj(response, handle)
        except Exception as exc:
            raise RuntimeError(f"Official demo download failed: {exc}") from exc
        observed_size = part.stat().st_size
        if observed_size != int(metadata["compressed_bytes"]):
            raise RuntimeError(
                f"Downloaded demo size mismatch: {observed_size} != "
                f"{metadata['compressed_bytes']}; incomplete content remains at {part}"
            )
        part.replace(destination)

    observed_sha256 = _sha256(destination)
    official_digest = metadata.get("digest")
    checksum_status = "locally_recorded_no_official_digest"
    if isinstance(official_digest, str) and official_digest.startswith("sha256:"):
        expected = official_digest.split(":", 1)[1].lower()
        if observed_sha256 != expected:
            raise RuntimeError(
                f"Official demo SHA-256 mismatch: {observed_sha256} != {expected}"
            )
        checksum_status = "verified_official_sha256"
    return {
        "path": str(destination.resolve()),
        "compressed_bytes": destination.stat().st_size,
        "sha256": observed_sha256,
        "checksum_status": checksum_status,
    }


def _safe_extract(archive: Path, destination: Path) -> Path:
    destination.mkdir(parents=True, exist_ok=True)
    root = destination.resolve()
    with tarfile.open(archive, mode="r:gz") as handle:
        for member in handle.getmembers():
            if member.issym() or member.islnk():
                raise RuntimeError(
                    f"Refusing link in official demo archive: {member.name}"
                )
            target = (destination / member.name).resolve()
            if target != root and root not in target.parents:
                raise RuntimeError(
                    f"Refusing unsafe path in official demo archive: {member.name}"
                )
        handle.extractall(destination)
    candidates = sorted(destination.rglob("test_reads.fa"))
    if len(candidates) != 1:
        raise RuntimeError(
            "Expected exactly one test_reads.fa after extraction; "
            f"found {len(candidates)}"
        )
    return candidates[0].parent


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Fetch and optionally run the official bulk CIRI-long demo. "
            "This is detector interoperability, not single-cell validation."
        )
    )
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the no-network plan and exit.",
    )
    parser.add_argument(
        "--metadata-only",
        action="store_true",
        help="Query official release metadata but do not download.",
    )
    parser.add_argument(
        "--execute",
        action="store_true",
        help="After download/extraction, build the BWA index and run call/collapse/import.",
    )
    parser.add_argument("--ciri-long", default="CIRI-long")
    parser.add_argument("--bwa", default="bwa")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--timeout-seconds", type=float, default=60.0)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.threads < 1:
        raise SystemExit("--threads must be >= 1")
    plan = {
        "schema_version": "circyto.ciri_long_demo_plan.v1",
        "official_source": OFFICIAL_SOURCE,
        "release_api": RELEASE_API,
        "release_tag": RELEASE_TAG,
        "asset_name": ASSET_NAME,
        "outdir": str(args.outdir.resolve()),
        "metadata_only": bool(args.metadata_only),
        "execute": bool(args.execute),
        "network_accessed": False,
        "scientific_boundary": (
            "Official bulk RCRT-library detector interoperability only; "
            "not single-cell validation."
        ),
        "redistribution": "The demo archive is fetched on demand and is not redistributed.",
    }
    if args.dry_run:
        print(json.dumps(plan, indent=2, sort_keys=True))
        return 0

    args.outdir.mkdir(parents=True, exist_ok=True)
    metadata = query_demo_metadata(timeout_seconds=args.timeout_seconds)
    metadata_path = args.outdir / "official_demo_metadata.json"
    metadata_path.write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    if args.metadata_only:
        print(json.dumps(metadata, indent=2, sort_keys=True))
        return 0

    archive = args.outdir / ASSET_NAME
    download = _download(
        metadata=metadata,
        destination=archive,
        timeout_seconds=args.timeout_seconds,
    )
    demo_dir = _safe_extract(archive, args.outdir / "extracted")
    result: dict[str, Any] = {
        **plan,
        "network_accessed": True,
        "metadata": metadata,
        "download": download,
        "demo_dir": str(demo_dir.resolve()),
        "executed": False,
    }
    if args.execute:
        reads = demo_dir / "test_reads.fa"
        reference = demo_dir / "mm10_chr12.fa"
        gtf = demo_dir / "mm10_chr12.gtf"
        for path in (reads, reference, gtf):
            if not path.is_file():
                raise RuntimeError(f"Official demo asset is missing: {path}")
        bwa_result = subprocess.run(
            [args.bwa, "index", "-a", "bwtsw", str(reference)],
            capture_output=True,
            text=True,
            check=False,
            shell=False,
        )
        if bwa_result.returncode != 0:
            raise RuntimeError(
                f"BWA indexing failed ({bwa_result.returncode}): {bwa_result.stderr}"
            )
        manifest = args.outdir / "ciri_long_demo_manifest.tsv"
        write_ciri_long_manifest_tsv(
            [
                CiriLongManifestRow(
                    sample_id="official_demo",
                    reads_path=str(reads),
                    source_accession="CIRI-long-v0.6-alpha-demo",
                    dataset_id="ciri_long_official_demo",
                    reference_id="mm10_chr12_demo",
                    reference_build="mm10",
                    extra={
                        "observation_unit": "bulk_demo",
                        "protocol_evidence": OFFICIAL_SOURCE,
                    },
                )
            ],
            manifest,
        )
        call_summary = run_ciri_long_call_stage(
            manifest_path=manifest,
            reference_fasta=reference,
            outdir=args.outdir,
            gtf=gtf,
            threads=args.threads,
            ciri_long=args.ciri_long,
            bwa=args.bwa,
            execute=True,
        )
        collapse_summary = run_ciri_long_collapse_stage(
            call_manifest_path=Path(call_summary["call_manifest"]),
            reference_fasta=reference,
            outdir=args.outdir,
            prefix="cohort",
            gtf=gtf,
            threads=args.threads,
            ciri_long=args.ciri_long,
            bwa=args.bwa,
            execute=True,
        )
        import_summary = normalize_ciri_long_outputs(
            collapse_dir=args.outdir / "ciri_long" / "collapse",
            outdir=args.outdir,
            prefix="cohort",
        )
        result.update(
            {
                "executed": True,
                "bwa_index_argv": [
                    args.bwa,
                    "index",
                    "-a",
                    "bwtsw",
                    str(reference),
                ],
                "call_summary": call_summary,
                "collapse_summary": collapse_summary,
                "import_summary": import_summary,
            }
        )
    result_path = args.outdir / "official_demo_smoke_summary.json"
    result_path.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
