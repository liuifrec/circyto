from __future__ import annotations

import hashlib
import io
import json
import socket
from pathlib import Path
from typing import Any

import pytest

from circyto.pipeline.nanopore_archive import (
    ArchiveIdentityError,
    DownloadIntegrityError,
    EnaFastqFile,
    download_ena_fastq,
    parse_ena_run_record,
    validate_ena_run_identity,
)


def _record(**updates: Any) -> dict[str, str]:
    record = {
        "run_accession": "SRR4048177",
        "study_accession": "PRJNA339767",
        "sample_accession": "SAMN05607017",
        "experiment_accession": "SRX2039007",
        "scientific_name": "Mus musculus",
        "instrument_platform": "OXFORD_NANOPORE",
        "instrument_model": "MinION",
        "library_strategy": "RNA-Seq",
        "library_source": "TRANSCRIPTOMIC",
        "library_selection": "RACE",
        "base_count": "70823600",
        "read_count": "52696",
        "fastq_bytes": "8",
        "fastq_ftp": (
            "ftp.sra.ebi.ac.uk/vol1/fastq/SRR404/007/SRR4048177/"
            "SRR4048177_1.fastq.gz"
        ),
        "fastq_md5": hashlib.md5(b"abcdefgh", usedforsecurity=False).hexdigest(),
    }
    record.update({key: str(value) for key, value in updates.items()})
    return record


def _expectation(metadata: Any) -> dict[str, Any]:
    fastq = metadata.fastq_files[0]
    return {
        "hard_identity": {
            "run_accession": metadata.run_accession,
            "study_accession": metadata.study_accession,
            "instrument_platform": metadata.instrument_platform,
            "instrument_model": metadata.instrument_model,
            "library_strategy": metadata.library_strategy,
            "library_source": metadata.library_source,
            "archive_library_selection": metadata.archive_library_selection,
            "fastq_files": [{"filename": fastq.filename, "md5": fastq.md5}],
        },
        "warning_baseline": {
            "read_count": metadata.read_count,
            "base_count": metadata.base_count,
            "fastq_files": [
                {
                    "filename": fastq.filename,
                    "url": fastq.url,
                    "compressed_bytes": fastq.compressed_bytes,
                }
            ],
        },
    }


def test_ena_identity_mismatch_is_hard_failure() -> None:
    metadata = parse_ena_run_record(_record())
    expectation = _expectation(metadata)
    expectation["hard_identity"]["instrument_platform"] = "ILLUMINA"

    with pytest.raises(ArchiveIdentityError, match="instrument_platform"):
        validate_ena_run_identity(metadata, expectation)


def test_ena_count_and_size_changes_are_warning_only() -> None:
    baseline = parse_ena_run_record(_record())
    expectation = _expectation(baseline)
    changed = parse_ena_run_record(
        _record(read_count=53000, base_count=71000000, fastq_bytes=9)
    )

    warnings = validate_ena_run_identity(changed, expectation)

    assert len(warnings) == 3
    assert any("read_count" in warning for warning in warnings)
    assert any("base_count" in warning for warning in warnings)
    assert any("compressed_bytes" in warning for warning in warnings)


def test_ena_checksum_change_is_hard_failure() -> None:
    metadata = parse_ena_run_record(_record())
    expectation = _expectation(metadata)
    expectation["hard_identity"]["fastq_files"][0]["md5"] = "0" * 32

    with pytest.raises(ArchiveIdentityError, match="checksum integrity violation"):
        validate_ena_run_identity(metadata, expectation)


def test_ena_file_lists_must_have_matching_lengths() -> None:
    with pytest.raises(ValueError, match="different lengths"):
        parse_ena_run_record(_record(fastq_md5="a" * 32 + ";" + "b" * 32))


class _Response(io.BytesIO):
    def __init__(
        self,
        content: bytes,
        *,
        status: int,
        headers: dict[str, str],
    ) -> None:
        super().__init__(content)
        self.status = status
        self.headers = headers

    def __enter__(self) -> "_Response":
        return self

    def __exit__(self, *args: Any) -> None:
        self.close()

    def getcode(self) -> int:
        return self.status


class _FakeUrlOpen:
    def __init__(
        self,
        content: bytes,
        *,
        ignore_range: bool = False,
        etag: str = '"fixture-v1"',
    ) -> None:
        self.content = content
        self.ignore_range = ignore_range
        self.etag = etag
        self.requests: list[Any] = []

    def __call__(self, request: Any, *, timeout: float) -> _Response:
        self.requests.append(request)
        common = {
            "ETag": self.etag,
            "Last-Modified": "Wed, 29 Jul 2026 00:00:00 GMT",
        }
        if request.get_method() == "HEAD":
            return _Response(
                b"",
                status=200,
                headers={**common, "Content-Length": str(len(self.content))},
            )
        range_header = request.headers.get("Range")
        if range_header and not self.ignore_range:
            offset = int(range_header.removeprefix("bytes=").removesuffix("-"))
            return _Response(
                self.content[offset:],
                status=206,
                headers={
                    **common,
                    "Content-Range": (
                        f"bytes {offset}-{len(self.content) - 1}/{len(self.content)}"
                    ),
                },
            )
        return _Response(self.content, status=200, headers=common)


def _remote(content: bytes, *, md5: str | None = None) -> EnaFastqFile:
    return EnaFastqFile(
        filename="SRR4048177_1.fastq.gz",
        url=(
            "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR404/007/SRR4048177/"
            "SRR4048177_1.fastq.gz"
        ),
        md5=md5 or hashlib.md5(content, usedforsecurity=False).hexdigest(),
        compressed_bytes=len(content),
    )


def test_download_restarts_when_server_ignores_range(tmp_path: Path) -> None:
    content = b"abcdefgh"
    remote = _remote(content)
    destination = tmp_path / remote.filename
    part = Path(str(destination) + ".part")
    state = Path(str(part) + ".json")
    part.write_bytes(content[:3])
    state.write_text(
        json.dumps(
            {
                "url": remote.url,
                "compressed_bytes": len(content),
                "md5": remote.md5,
                "etag": '"fixture-v1"',
                "last_modified": "Wed, 29 Jul 2026 00:00:00 GMT",
            }
        ),
        encoding="utf-8",
    )
    opener = _FakeUrlOpen(content, ignore_range=True)

    summary = download_ena_fastq(remote, destination, urlopen=opener)

    assert destination.read_bytes() == content
    assert summary["range_ignored"] is True
    assert summary["resumed_from_bytes"] == 0
    assert any(request.headers.get("Range") == "bytes=3-" for request in opener.requests)


def test_download_rejects_md5_failure_without_promoting_file(tmp_path: Path) -> None:
    content = b"abcdefgh"
    remote = _remote(content, md5="0" * 32)
    destination = tmp_path / remote.filename

    with pytest.raises(DownloadIntegrityError, match="MD5 verification failed"):
        download_ena_fastq(remote, destination, urlopen=_FakeUrlOpen(content))

    assert not destination.exists()
    assert Path(str(destination) + ".part.md5_failed").read_bytes() == content


def test_oversized_partial_is_discarded_before_download(tmp_path: Path) -> None:
    content = b"abcdefgh"
    remote = _remote(content)
    destination = tmp_path / remote.filename
    Path(str(destination) + ".part").write_bytes(content + b"too-large")

    summary = download_ena_fastq(remote, destination, urlopen=_FakeUrlOpen(content))

    assert summary["status"] == "downloaded_verified"
    assert destination.read_bytes() == content


def test_changed_remote_identity_discards_stale_partial(tmp_path: Path) -> None:
    content = b"abcdefgh"
    remote = _remote(content)
    destination = tmp_path / remote.filename
    part = Path(str(destination) + ".part")
    state = Path(str(part) + ".json")
    part.write_bytes(content[:3])
    state.write_text(
        json.dumps(
            {
                "url": remote.url,
                "compressed_bytes": len(content),
                "md5": remote.md5,
                "etag": '"old-object"',
                "last_modified": "Wed, 29 Jul 2026 00:00:00 GMT",
            }
        ),
        encoding="utf-8",
    )
    opener = _FakeUrlOpen(content, etag='"new-object"')

    summary = download_ena_fastq(remote, destination, urlopen=opener)

    assert summary["resumed_from_bytes"] == 0
    assert destination.read_bytes() == content
    assert not any(request.headers.get("Range") for request in opener.requests)


def test_timeout_keeps_partial_for_future_resume(tmp_path: Path) -> None:
    content = b"abcdefgh"
    remote = _remote(content)
    destination = tmp_path / remote.filename

    class _TimeoutAfterHead(_FakeUrlOpen):
        def __call__(self, request: Any, *, timeout: float) -> _Response:
            if request.get_method() == "HEAD":
                return super().__call__(request, timeout=timeout)
            raise socket.timeout("synthetic timeout")

    with pytest.raises(RuntimeError, match="partial file is resumable"):
        download_ena_fastq(remote, destination, urlopen=_TimeoutAfterHead(content))

    assert Path(str(destination) + ".part.json").is_file()
    assert not destination.exists()
