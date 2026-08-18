from __future__ import annotations

import os
import shutil
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


@dataclass(frozen=True)
class Ciri3Resolution:
    home: PathResolution
    jar: PathResolution
    bin: PathResolution
    java: PathResolution

    @property
    def has_home(self) -> bool:
        return self.home.found

    @property
    def has_jar(self) -> bool:
        return self.jar.found

    @property
    def has_bin(self) -> bool:
        return self.bin.found

    @property
    def has_java(self) -> bool:
        return self.java.found


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
    if repo_root is not None:
        return repo_root / "tools"

    packaged_tools = get_package_root() / "resources" / "tools"
    if packaged_tools.is_dir():
        return packaged_tools
    return None


def _tools_source_label(tools_dir: Path) -> str:
    packaged_tools = get_package_root() / "resources" / "tools"
    if tools_dir == packaged_tools:
        return "packaged-tools"
    if os.environ.get("CIRCYTO_TOOLS_DIR"):
        return "env:CIRCYTO_TOOLS_DIR"
    return "repo-tools"


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


def get_packaged_smoke_demo_dir(name: str = "smoke_demo") -> Path:
    return get_package_root() / "resources" / name


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


def _which_resolution(
    label: str,
    *,
    override: str | Path | None,
    override_source: str | None,
    names: Iterable[str],
) -> PathResolution:
    checked: list[Path] = []

    if override is not None:
        candidate = _normalize_path(override)
        checked.append(candidate)
        if candidate.exists():
            return PathResolution(label=label, resolved_path=candidate, checked_paths=(candidate,), source=override_source)

    for name in names:
        found = shutil.which(name)
        if not found:
            continue
        candidate = _normalize_path(found)
        checked.append(candidate)
        if candidate.exists():
            return PathResolution(label=label, resolved_path=candidate, checked_paths=_dedupe_paths(checked), source=f"path:{name}")

    return PathResolution(label=label, resolved_path=None, checked_paths=_dedupe_paths(checked), source=None)


def _iter_configured_search_roots() -> tuple[Path, ...]:
    roots: list[Path] = []

    env_value = os.environ.get("CIRCYTO_CIRI3_SEARCH_ROOTS", "")
    for raw in env_value.split(os.pathsep):
        raw = raw.strip()
        if raw:
            roots.append(_normalize_path(raw))

    tools_dir = get_tools_dir()
    if tools_dir is not None:
        roots.extend(
            [
                tools_dir / "CIRI3",
                tools_dir / "ciri3",
            ]
        )

    return _dedupe_paths(roots)


def _parse_version_tuple(path: Path) -> tuple[int, ...]:
    stem = path.stem
    digits = "".join(ch if ch.isdigit() or ch == "." else " " for ch in stem).split()
    if not digits:
        return tuple()
    parts = digits[-1].split(".")
    values: list[int] = []
    for part in parts:
        try:
            values.append(int(part))
        except ValueError:
            values.append(0)
    return tuple(values)


def find_java_executable(override: str | Path | None = None) -> PathResolution:
    return _which_resolution(
        "java",
        override=override if override is not None else os.environ.get("CIRCYTO_CIRI3_JAVA"),
        override_source="override" if override is not None else "env:CIRCYTO_CIRI3_JAVA",
        names=("java",),
    )


def find_ciri3_home(override: str | Path | None = None) -> PathResolution:
    candidates: list[tuple[Path, str]] = []
    for root in _iter_configured_search_roots():
        candidates.append((root, "configured-search-root"))

    return _resolve_from_candidates(
        "CIRI3 home",
        override=override if override is not None else os.environ.get("CIRCYTO_CIRI3_HOME"),
        override_source="override" if override is not None else "env:CIRCYTO_CIRI3_HOME",
        candidates=candidates,
    )


def find_ciri3_jar(override: str | Path | None = None) -> PathResolution:
    home = find_ciri3_home()
    candidates: list[tuple[Path, str]] = []

    if home.resolved_path is not None:
        root = home.resolved_path
        direct = sorted(root.glob("CIRI3*.jar"))
        for match in sorted(direct, key=lambda path: (_parse_version_tuple(path), path.name), reverse=True):
            candidates.append((match, "ciri3-home"))
    for root in _iter_configured_search_roots():
        direct = sorted(root.glob("CIRI3*.jar"))
        for match in sorted(direct, key=lambda path: (_parse_version_tuple(path), path.name), reverse=True):
            candidates.append((match, "configured-search-root"))

    return _resolve_from_candidates(
        "CIRI3 jar",
        override=override if override is not None else os.environ.get("CIRCYTO_CIRI3_JAR"),
        override_source="override" if override is not None else "env:CIRCYTO_CIRI3_JAR",
        candidates=candidates,
    )


def find_ciri3_bin(override: str | Path | None = None) -> PathResolution:
    home = find_ciri3_home()
    candidates: list[tuple[Path, str]] = []
    wrapper_names = ("CIRI3", "ciri3", "CIRI3.sh", "ciri3.sh", "CIRI3.pl", "ciri3.pl")

    if home.resolved_path is not None:
        for name in wrapper_names:
            candidates.append((home.resolved_path / "bin" / name, "ciri3-home"))
            candidates.append((home.resolved_path / "scripts" / name, "ciri3-home"))

    return _resolve_from_candidates(
        "CIRI3 wrapper",
        override=override if override is not None else os.environ.get("CIRCYTO_CIRI3_BIN"),
        override_source="override" if override is not None else "env:CIRCYTO_CIRI3_BIN",
        candidates=[
            *candidates,
            *[
                (_normalize_path(found), f"path:{name}")
                for name in wrapper_names
                for found in ([shutil.which(name)] if shutil.which(name) else [])
            ],
        ],
    )


def resolve_ciri3_installation() -> Ciri3Resolution:
    return Ciri3Resolution(
        home=find_ciri3_home(),
        jar=find_ciri3_jar(),
        bin=find_ciri3_bin(),
        java=find_java_executable(),
    )


def find_ciri_full_jar(override: str | Path | None = None) -> PathResolution:
    tools_dir = get_tools_dir()
    candidates: list[tuple[Path, str]] = []

    if tools_dir is not None:
        tools_source = _tools_source_label(tools_dir)
        candidates.extend(
            [
                (tools_dir / "CIRI-full_v2.0" / "CIRI-full.jar", tools_source),
                (tools_dir / "CIRI-full_v2.0" / "CIRI_Full.jar", tools_source),
                (tools_dir / "CIRI-full.jar", tools_source),
                (tools_dir / "CIRI_Full.jar", tools_source),
            ]
        )
        if tools_dir.exists():
            for match in sorted(tools_dir.rglob("CIRI-full*.jar")):
                candidates.append((match, f"{tools_source}-scan"))
            for match in sorted(tools_dir.rglob("CIRI_Full*.jar")):
                candidates.append((match, f"{tools_source}-scan"))

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
        candidates.append(
            (
                tools_dir / "CIRI-full_v2.0" / "bin" / "ciri_full_adapter.sh",
                _tools_source_label(tools_dir),
            )
        )

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
        candidates.append(
            (
                tools_dir / "CIRI-full_v2.0" / "bin" / "ciri2_adapter.sh",
                _tools_source_label(tools_dir),
            )
        )

    return _resolve_from_candidates(
        "CIRI2 adapter",
        override=override if override is not None else os.environ.get("CIRCYTO_CIRI2_ADAPTER"),
        override_source="override" if override is not None else "env:CIRCYTO_CIRI2_ADAPTER",
        candidates=candidates,
    )


def format_missing_resolution(label: str, resolution: PathResolution) -> str:
    checked = "\n".join(f"  - {path}" for path in resolution.checked_paths) or "  - <no candidates checked>"
    return f"{label} not found.\nChecked:\n{checked}"
