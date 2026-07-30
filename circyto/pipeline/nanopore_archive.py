from __future__ import annotations

import hashlib
import json
import os
import socket
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Mapping


ENA_REPORT_FIELDS = [
    "run_accession",
    "instrument_platform",
    "instrument_model",
    "library_strategy",
    "library_source",
    "library_selection",
    "base_count",
    "read_count",
    "fastq_bytes",
    "fastq_ftp",
    "fastq_md5",
    "study_accession",
    "sample_accession",
    "experiment_accession",
    "scientific_name",
]

UrlOpen = Callable[..., Any]


class ArchiveIdentityError(ValueError):
    """The archive record no longer identifies the documented run."""


class DownloadIntegrityError(RuntimeError):
    """A remote or downloaded archive object failed an integrity check."""


@dataclass(frozen=True)
class EnaFastqFile:
    filename: str
    url: str
    md5: str
    compressed_bytes: int

    def to_dict(self) -> dict[str, Any]:
        return {
            "filename": self.filename,
            "url": self.url,
            "md5": self.md5,
            "compressed_bytes": self.compressed_bytes,
        }


@dataclass(frozen=True)
class EnaRunMetadata:
    run_accession: str
    study_accession: str
    sample_accession: str
    experiment_accession: str
    scientific_name: str
    instrument_platform: str
    instrument_model: str
    library_strategy: str
    library_source: str
    archive_library_selection: str
    base_count: int
    read_count: int
    fastq_files: tuple[EnaFastqFile, ...]
    query_url: str = ""

    def to_dict(self) -> dict[str, Any]:
        return {
            "run_accession": self.run_accession,
            "study_accession": self.study_accession,
            "sample_accession": self.sample_accession,
            "experiment_accession": self.experiment_accession,
            "scientific_name": self.scientific_name,
            "instrument_platform": self.instrument_platform,
            "instrument_model": self.instrument_model,
            "library_strategy": self.library_strategy,
            "library_source": self.library_source,
            "archive_library_selection": self.archive_library_selection,
            "base_count": self.base_count,
            "read_count": self.read_count,
            "fastq_files": [item.to_dict() for item in self.fastq_files],
            "query_url": self.query_url,
        }


def _split_ena_file_list(value: Any, *, field: str) -> list[str]:
    text = str(value or "").strip()
    if not text:
        return []
    values = [part.strip() for part in text.split(";")]
    if any(not part for part in values):
        raise ValueError(f"ENA field {field!r} contains an empty file-list entry")
    return values


def _positive_int(value: Any, *, field: str) -> int:
    try:
        parsed = int(str(value).strip())
    except (TypeError, ValueError) as exc:
        raise ValueError(f"ENA field {field!r} is not an integer: {value!r}") from exc
    if parsed < 0:
        raise ValueError(f"ENA field {field!r} must be non-negative: {parsed}")
    return parsed


def _normalize_ena_fastq_url(value: str) -> str:
    text = value.strip()
    if "://" not in text:
        text = "https://" + text
    parsed = urllib.parse.urlsplit(text)
    if parsed.scheme == "ftp":
        parsed = parsed._replace(scheme="https")
    if parsed.scheme != "https":
        raise ValueError(f"ENA FASTQ URL must use HTTPS or archive FTP form: {value!r}")
    hostname = (parsed.hostname or "").lower()
    if hostname != "ftp.sra.ebi.ac.uk":
        raise ValueError(f"Unexpected ENA FASTQ host {hostname!r} in {value!r}")
    if parsed.username or parsed.password or parsed.query or parsed.fragment:
        raise ValueError(f"Unsafe ENA FASTQ URL components in {value!r}")
    path = urllib.parse.unquote(parsed.path)
    if not path.startswith("/vol1/fastq/") or "/../" in path or path.endswith("/.."):
        raise ValueError(f"Unexpected ENA FASTQ path in {value!r}")
    return urllib.parse.urlunsplit(parsed)


def parse_ena_run_record(record: Mapping[str, Any], *, query_url: str = "") -> EnaRunMetadata:
    urls = _split_ena_file_list(record.get("fastq_ftp"), field="fastq_ftp")
    md5s = _split_ena_file_list(record.get("fastq_md5"), field="fastq_md5")
    byte_values = _split_ena_file_list(record.get("fastq_bytes"), field="fastq_bytes")
    if not urls:
        raise ValueError("ENA run record contains no FASTQ URLs")
    if len({len(urls), len(md5s), len(byte_values)}) != 1:
        raise ValueError(
            "ENA FASTQ file lists have different lengths: "
            f"urls={len(urls)}, md5={len(md5s)}, bytes={len(byte_values)}"
        )

    files: list[EnaFastqFile] = []
    seen_names: set[str] = set()
    seen_urls: set[str] = set()
    for raw_url, raw_md5, raw_bytes in zip(urls, md5s, byte_values):
        url = _normalize_ena_fastq_url(raw_url)
        filename = Path(urllib.parse.urlsplit(url).path).name
        if not filename.endswith((".fastq.gz", ".fq.gz")):
            raise ValueError(f"ENA object is not a compressed FASTQ: {filename!r}")
        md5 = raw_md5.lower()
        if len(md5) != 32 or any(char not in "0123456789abcdef" for char in md5):
            raise ValueError(f"Invalid ENA MD5 for {filename}: {raw_md5!r}")
        compressed_bytes = _positive_int(raw_bytes, field=f"fastq_bytes[{filename}]")
        if compressed_bytes <= 0:
            raise ValueError(f"ENA compressed FASTQ size must be positive for {filename}")
        if filename in seen_names or url in seen_urls:
            raise ValueError(f"Duplicate ENA FASTQ file entry: {filename}")
        seen_names.add(filename)
        seen_urls.add(url)
        files.append(
            EnaFastqFile(
                filename=filename,
                url=url,
                md5=md5,
                compressed_bytes=compressed_bytes,
            )
        )

    required_text = (
        "run_accession",
        "study_accession",
        "sample_accession",
        "experiment_accession",
        "scientific_name",
        "instrument_platform",
        "instrument_model",
        "library_strategy",
        "library_source",
        "library_selection",
    )
    values: dict[str, str] = {}
    for field in required_text:
        value = str(record.get(field) or "").strip()
        if not value:
            raise ValueError(f"ENA run record has empty required field {field!r}")
        values[field] = value

    return EnaRunMetadata(
        run_accession=values["run_accession"],
        study_accession=values["study_accession"],
        sample_accession=values["sample_accession"],
        experiment_accession=values["experiment_accession"],
        scientific_name=values["scientific_name"],
        instrument_platform=values["instrument_platform"],
        instrument_model=values["instrument_model"],
        library_strategy=values["library_strategy"],
        library_source=values["library_source"],
        archive_library_selection=values["library_selection"],
        base_count=_positive_int(record.get("base_count"), field="base_count"),
        read_count=_positive_int(record.get("read_count"), field="read_count"),
        fastq_files=tuple(files),
        query_url=query_url,
    )


def ena_run_report_url(accession: str) -> str:
    query = urllib.parse.urlencode(
        {
            "accession": accession,
            "result": "read_run",
            "fields": ",".join(ENA_REPORT_FIELDS),
            "format": "json",
            "download": "true",
        }
    )
    return f"https://www.ebi.ac.uk/ena/portal/api/filereport?{query}"


def query_ena_run(
    accession: str,
    *,
    timeout_seconds: float = 30.0,
    urlopen: UrlOpen = urllib.request.urlopen,
) -> EnaRunMetadata:
    if not accession.strip():
        raise ValueError("ENA accession must not be empty")
    url = ena_run_report_url(accession.strip())
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "circyto-nanopore-interoperability/1"},
        method="GET",
    )
    try:
        with urlopen(request, timeout=timeout_seconds) as response:
            payload = json.loads(response.read().decode("utf-8"))
    except (TimeoutError, socket.timeout, urllib.error.URLError) as exc:
        raise RuntimeError(
            f"ENA metadata request timed out or failed for {accession}: {exc}"
        ) from exc
    except json.JSONDecodeError as exc:
        raise RuntimeError(f"ENA returned invalid JSON for {accession}: {exc}") from exc
    if not isinstance(payload, list) or len(payload) != 1 or not isinstance(payload[0], dict):
        count = len(payload) if isinstance(payload, list) else "non-list"
        raise RuntimeError(
            f"ENA metadata query for {accession} returned {count} records; expected exactly 1"
        )
    metadata = parse_ena_run_record(payload[0], query_url=url)
    if metadata.run_accession != accession:
        raise ArchiveIdentityError(
            f"ENA query requested {accession!r} but returned {metadata.run_accession!r}"
        )
    return metadata


def load_expected_run(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("schema_version") != "circyto.ena_run_expectation.v1":
        raise ValueError(f"Unsupported ENA expectation schema in {path}")
    return payload


def validate_ena_run_identity(
    metadata: EnaRunMetadata,
    expectation: Mapping[str, Any],
) -> list[str]:
    hard = expectation.get("hard_identity")
    if not isinstance(hard, Mapping):
        raise ValueError("ENA expectation is missing hard_identity")
    scalar_fields = (
        "run_accession",
        "study_accession",
        "instrument_platform",
        "instrument_model",
        "library_strategy",
        "library_source",
        "archive_library_selection",
    )
    errors: list[str] = []
    for field in scalar_fields:
        expected = str(hard.get(field) or "")
        observed = str(getattr(metadata, field))
        if not expected:
            errors.append(f"expectation missing hard identity field {field}")
        elif observed != expected:
            errors.append(f"{field}: expected {expected!r}, observed {observed!r}")

    expected_files_raw = hard.get("fastq_files")
    if not isinstance(expected_files_raw, list) or not expected_files_raw:
        errors.append("expectation missing hard fastq_files identity")
        expected_files_raw = []
    expected_files = {
        str(item.get("filename") or ""): str(item.get("md5") or "").lower()
        for item in expected_files_raw
        if isinstance(item, Mapping)
    }
    observed_files = {item.filename: item.md5 for item in metadata.fastq_files}
    if set(observed_files) != set(expected_files):
        errors.append(
            "FASTQ file identity changed: "
            f"expected {sorted(expected_files)}, observed {sorted(observed_files)}"
        )
    for filename in sorted(set(observed_files) & set(expected_files)):
        if observed_files[filename] != expected_files[filename]:
            errors.append(
                f"checksum integrity violation for {filename}: "
                f"expected {expected_files[filename]}, observed {observed_files[filename]}"
            )
    if errors:
        raise ArchiveIdentityError(
            "ENA hard identity check failed for "
            f"{metadata.run_accession}: " + "; ".join(errors)
        )

    warnings: list[str] = []
    baseline = expectation.get("warning_baseline")
    if not isinstance(baseline, Mapping):
        return warnings
    for field in ("read_count", "base_count"):
        expected = baseline.get(field)
        observed = getattr(metadata, field)
        if expected is not None and int(expected) != observed:
            warnings.append(
                f"ENA mutable metadata changed for {field}: expected {int(expected)}, "
                f"observed {observed}"
            )
    expected_by_name = {
        str(item.get("filename") or ""): item
        for item in baseline.get("fastq_files", [])
        if isinstance(item, Mapping)
    }
    for observed in metadata.fastq_files:
        expected = expected_by_name.get(observed.filename, {})
        if expected.get("url") and str(expected["url"]) != observed.url:
            warnings.append(
                f"ENA mutable metadata changed for {observed.filename} URL: "
                f"expected {expected['url']!r}, observed {observed.url!r}"
            )
        if expected.get("compressed_bytes") is not None and int(
            expected["compressed_bytes"]
        ) != observed.compressed_bytes:
            warnings.append(
                f"ENA mutable metadata changed for {observed.filename} compressed_bytes: "
                f"expected {int(expected['compressed_bytes'])}, "
                f"observed {observed.compressed_bytes}"
            )
    return warnings


def sha256_file(path: Path, *, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def md5_file(path: Path, *, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def _response_status(response: Any) -> int:
    status = getattr(response, "status", None)
    if status is None and hasattr(response, "getcode"):
        status = response.getcode()
    return int(status or 200)


def _header(response: Any, name: str) -> str:
    headers = getattr(response, "headers", {})
    return str(headers.get(name, "") or "").strip()


def _discard_partial(part_path: Path, state_path: Path) -> None:
    if part_path.exists():
        part_path.unlink()
    if state_path.exists():
        state_path.unlink()


def download_ena_fastq(
    remote: EnaFastqFile,
    destination: Path,
    *,
    timeout_seconds: float = 60.0,
    chunk_size: int = 1024 * 1024,
    urlopen: UrlOpen = urllib.request.urlopen,
) -> dict[str, Any]:
    destination.parent.mkdir(parents=True, exist_ok=True)
    part_path = Path(str(destination) + ".part")
    state_path = Path(str(part_path) + ".json")

    if destination.exists():
        if destination.stat().st_size != remote.compressed_bytes:
            raise DownloadIntegrityError(
                f"Existing destination has unexpected size: {destination}"
            )
        observed_md5 = md5_file(destination)
        if observed_md5 != remote.md5:
            raise DownloadIntegrityError(
                f"Existing destination MD5 mismatch for {destination}: "
                f"expected {remote.md5}, observed {observed_md5}"
            )
        return {
            "path": str(destination.resolve()),
            "status": "reused_verified",
            "compressed_bytes": destination.stat().st_size,
            "md5": observed_md5,
            "resumed_from_bytes": destination.stat().st_size,
            "range_ignored": False,
        }

    head_request = urllib.request.Request(
        remote.url,
        headers={"User-Agent": "circyto-nanopore-interoperability/1"},
        method="HEAD",
    )
    try:
        with urlopen(head_request, timeout=timeout_seconds) as head_response:
            head_size_text = _header(head_response, "Content-Length")
            remote_size = int(head_size_text) if head_size_text else remote.compressed_bytes
            etag = _header(head_response, "ETag")
            last_modified = _header(head_response, "Last-Modified")
    except (TimeoutError, socket.timeout, urllib.error.URLError) as exc:
        raise RuntimeError(f"ENA HEAD request timed out or failed for {remote.url}: {exc}") from exc
    except ValueError as exc:
        raise DownloadIntegrityError(
            f"ENA returned an invalid Content-Length for {remote.url}: {head_size_text!r}"
        ) from exc
    if remote_size != remote.compressed_bytes:
        _discard_partial(part_path, state_path)
        raise DownloadIntegrityError(
            f"Remote object size differs from current ENA metadata for {remote.filename}: "
            f"HEAD={remote_size}, ENA={remote.compressed_bytes}"
        )

    current_identity = {
        "url": remote.url,
        "compressed_bytes": remote.compressed_bytes,
        "md5": remote.md5,
        "etag": etag,
        "last_modified": last_modified,
    }
    previous_identity: dict[str, Any] | None = None
    if state_path.exists():
        try:
            previous_identity = json.loads(state_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            previous_identity = None
    if part_path.exists() and (
        part_path.stat().st_size > remote_size or previous_identity != current_identity
    ):
        _discard_partial(part_path, state_path)

    offset = part_path.stat().st_size if part_path.exists() else 0
    state_path.write_text(
        json.dumps(current_identity, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    if offset == remote_size:
        observed_md5 = md5_file(part_path)
        if observed_md5 != remote.md5:
            failed_path = Path(str(part_path) + ".md5_failed")
            if failed_path.exists():
                failed_path.unlink()
            os.replace(part_path, failed_path)
            state_path.unlink(missing_ok=True)
            raise DownloadIntegrityError(
                f"MD5 verification failed for completed partial {remote.filename}: "
                f"expected {remote.md5}, observed {observed_md5}; retained as {failed_path}"
            )
        os.replace(part_path, destination)
        state_path.unlink(missing_ok=True)
        return {
            "path": str(destination.resolve()),
            "status": "resumed_verified",
            "compressed_bytes": remote_size,
            "md5": observed_md5,
            "resumed_from_bytes": offset,
            "range_ignored": False,
        }

    headers = {"User-Agent": "circyto-nanopore-interoperability/1"}
    if offset:
        headers["Range"] = f"bytes={offset}-"
        if etag:
            headers["If-Range"] = etag
        elif last_modified:
            headers["If-Range"] = last_modified
    get_request = urllib.request.Request(remote.url, headers=headers, method="GET")
    range_ignored = False
    try:
        with urlopen(get_request, timeout=timeout_seconds) as response:
            status = _response_status(response)
            response_etag = _header(response, "ETag")
            response_modified = _header(response, "Last-Modified")
            if etag and response_etag and response_etag != etag:
                _discard_partial(part_path, state_path)
                raise DownloadIntegrityError(
                    f"Remote object ETag changed during download of {remote.filename}"
                )
            if last_modified and response_modified and response_modified != last_modified:
                _discard_partial(part_path, state_path)
                raise DownloadIntegrityError(
                    f"Remote object Last-Modified changed during download of {remote.filename}"
                )

            mode = "wb"
            if offset and status == 206:
                content_range = _header(response, "Content-Range")
                expected_prefix = f"bytes {offset}-"
                if not content_range.startswith(expected_prefix) or not content_range.endswith(
                    f"/{remote_size}"
                ):
                    _discard_partial(part_path, state_path)
                    raise DownloadIntegrityError(
                        f"Invalid Content-Range while resuming {remote.filename}: "
                        f"{content_range!r}"
                    )
                mode = "ab"
            elif offset and status == 200:
                range_ignored = True
                offset = 0
                mode = "wb"
            elif status not in {200, 206}:
                raise RuntimeError(
                    f"Unexpected HTTP status {status} while downloading {remote.url}"
                )

            with part_path.open(mode) as handle:
                while chunk := response.read(chunk_size):
                    handle.write(chunk)
    except DownloadIntegrityError:
        raise
    except (TimeoutError, socket.timeout, urllib.error.URLError) as exc:
        raise RuntimeError(
            f"ENA download timed out or failed for {remote.url}; partial file is resumable: {exc}"
        ) from exc

    observed_size = part_path.stat().st_size
    if observed_size < remote_size:
        raise RuntimeError(
            f"ENA download ended early for {remote.filename}: "
            f"{observed_size}/{remote_size} bytes; partial file is resumable"
        )
    if observed_size > remote_size:
        _discard_partial(part_path, state_path)
        raise DownloadIntegrityError(
            f"Downloaded object exceeds expected size for {remote.filename}: "
            f"{observed_size}>{remote_size}"
        )

    observed_md5 = md5_file(part_path)
    if observed_md5 != remote.md5:
        failed_path = Path(str(part_path) + ".md5_failed")
        if failed_path.exists():
            failed_path.unlink()
        os.replace(part_path, failed_path)
        state_path.unlink(missing_ok=True)
        raise DownloadIntegrityError(
            f"MD5 verification failed for {remote.filename}: expected {remote.md5}, "
            f"observed {observed_md5}; retained as {failed_path}"
        )
    os.replace(part_path, destination)
    state_path.unlink(missing_ok=True)
    return {
        "path": str(destination.resolve()),
        "status": "downloaded_verified",
        "compressed_bytes": observed_size,
        "md5": observed_md5,
        "resumed_from_bytes": offset,
        "range_ignored": range_ignored,
    }
