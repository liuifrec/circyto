from __future__ import annotations

import os
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class PathResolution:
    label: str
    resolved_path: Path | None
    checked_paths: tuple[Path, ...]
    source: str | None = None

    @property
    def found(self) -> bool:
        return self.resolved_path is not None


def _normalize_path(path: str | Path) -> Path:
    candidate = Path(path).expanduser()
    if candidate.is_absolute():
        return candidate.resolve(strict=False)
    return candidate.resolve(strict=False)


def _dedupe_paths(paths: Iterable[Path]) -> tuple[Path, ...]:
    seen: set[Path] = set()
    ordered: list[Path] = []
    for path in paths:
        if path in seen:
            continue
        seen.add(path)
        ordered.append(path)
    return tuple(ordered)


@lru_cache(maxsize=1)
def get_package_root() -> Path:
    return Path(__file__).resolve().parent


@lru_cache(maxsize=1)
def get_repo_root() -> Path | None:
    env_repo_root = os.environ.get("CIRCYTO_REPO_ROOT")
    if env_repo_root:
        repo_root = _normalize_path(env_repo_root)
        if repo_root.exists():
            return repo_root

    package_root = get_package_root()
    for candidate in (package_root, *package_root.parents):
        if (candidate / "pyproject.toml").is_file() and (candidate / "circyto").is_dir():
            return candidate
    return None


def clear_path_resolution_caches() -> None:
    get_package_root.cache_clear()
    get_repo_root.cache_clear()


def get_tools_dir() -> Path | None:
    env_tools_dir = os.environ.get("CIRCYTO_TOOLS_DIR")
    if env_tools_dir:
        return _normalize_path(env_tools_dir)

    repo_root = get_repo_root()
    if repo_root is None:
        return None
    return repo_root / "tools"


def get_testdata_dir() -> Path | None:
    env_testdata_dir = os.environ.get("CIRCYTO_TESTDATA_DIR")
    if env_testdata_dir:
        return _normalize_path(env_testdata_dir)

    repo_root = get_repo_root()
    if repo_root is None:
        return None
    return repo_root / "testdata"


def get_bundled_smoke_testdata_dir(name: str = "smartseq2_smoke") -> Path | None:
    testdata_dir = get_testdata_dir()
    if testdata_dir is None:
        return None
    return testdata_dir / name


def resolve_relative_to(base_dir: Path, value: str | Path) -> Path:
    path = Path(value).expanduser()
    if path.is_absolute():
        return path.resolve(strict=False)
    return (base_dir / path).resolve(strict=False)


def resolve_manifest_path(manifest_path: Path, value: str | Path) -> Path:
    return resolve_relative_to(manifest_path.resolve(strict=False).parent, value)


def _resolve_from_candidates(
    label: str,
    *,
    override: str | Path | None,
    override_source: str | None,
    candidates: Iterable[tuple[Path, str]],
) -> PathResolution:
    checked: list[Path] = []

    if override is not None:
        candidate = _normalize_path(override)
        checked.append(candidate)
        if candidate.exists():
            return PathResolution(label=label, resolved_path=candidate, checked_paths=(candidate,), source=override_source)

    for candidate, source in candidates:
        normalized = _normalize_path(candidate)
        checked.append(normalized)
        if normalized.exists():
            return PathResolution(
                label=label,
                resolved_path=normalized,
                checked_paths=_dedupe_paths(checked),
                source=source,
            )

    return PathResolution(label=label, resolved_path=None, checked_paths=_dedupe_paths(checked), source=None)


def find_ciri_full_jar(override: str | Path | None = None) -> PathResolution:
    tools_dir = get_tools_dir()
    candidates: list[tuple[Path, str]] = []

    if tools_dir is not None:
        candidates.extend(
            [
                (tools_dir / "CIRI-full_v2.0" / "CIRI-full.jar", "repo-tools"),
                (tools_dir / "CIRI-full_v2.0" / "CIRI_Full.jar", "repo-tools"),
                (tools_dir / "CIRI-full.jar", "repo-tools"),
                (tools_dir / "CIRI_Full.jar", "repo-tools"),
            ]
        )
        if tools_dir.exists():
            for match in sorted(tools_dir.rglob("CIRI-full*.jar")):
                candidates.append((match, "repo-tools-scan"))
            for match in sorted(tools_dir.rglob("CIRI_Full*.jar")):
                candidates.append((match, "repo-tools-scan"))

    return _resolve_from_candidates(
        "CIRI-full jar",
        override=override if override is not None else os.environ.get("CIRCYTO_CIRI_FULL_JAR"),
        override_source="override" if override is not None else "env:CIRCYTO_CIRI_FULL_JAR",
        candidates=candidates,
    )


def find_ciri_full_adapter(override: str | Path | None = None) -> PathResolution:
    tools_dir = get_tools_dir()
    candidates: list[tuple[Path, str]] = []
    if tools_dir is not None:
        candidates.append((tools_dir / "CIRI-full_v2.0" / "bin" / "ciri_full_adapter.sh", "repo-tools"))

    return _resolve_from_candidates(
        "CIRI-full adapter",
        override=override if override is not None else os.environ.get("CIRCYTO_CIRI_FULL_ADAPTER"),
        override_source="override" if override is not None else "env:CIRCYTO_CIRI_FULL_ADAPTER",
        candidates=candidates,
    )


def find_ciri2_adapter(override: str | Path | None = None) -> PathResolution:
    tools_dir = get_tools_dir()
    candidates: list[tuple[Path, str]] = []
    if tools_dir is not None:
        candidates.append((tools_dir / "CIRI-full_v2.0" / "bin" / "ciri2_adapter.sh", "repo-tools"))

    return _resolve_from_candidates(
        "CIRI2 adapter",
        override=override if override is not None else os.environ.get("CIRCYTO_CIRI2_ADAPTER"),
        override_source="override" if override is not None else "env:CIRCYTO_CIRI2_ADAPTER",
        candidates=candidates,
    )


def format_missing_resolution(label: str, resolution: PathResolution) -> str:
    checked = "\n".join(f"  - {path}" for path in resolution.checked_paths) or "  - <no candidates checked>"
    return f"{label} not found.\nChecked:\n{checked}"
