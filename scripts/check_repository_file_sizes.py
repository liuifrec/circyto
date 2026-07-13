#!/usr/bin/env python
from __future__ import annotations

import argparse
import fnmatch
import subprocess
from pathlib import Path


DEFAULT_ALLOW_GLOBS = [
    "circyto/resources/smoke_demo/**",
    "testdata/**",
]


def git_z(args: list[str]) -> list[str]:
    result = subprocess.run(["git", *args], check=True, capture_output=True)
    text = result.stdout.decode("utf-8", errors="replace")
    return [item for item in text.split("\0") if item]


def tracked_files() -> list[str]:
    return git_z(["ls-files", "-z"])


def added_files(base_ref: str) -> set[str]:
    try:
        return set(git_z(["diff", "--name-only", "--diff-filter=A", "-z", f"{base_ref}...HEAD"]))
    except subprocess.CalledProcessError:
        print(f"[WARN] Could not compare against {base_ref}; treating no files as newly added.")
        return set()


def allowed(path: str, patterns: list[str]) -> bool:
    return any(fnmatch.fnmatch(path, pattern) for pattern in patterns)


def human_bytes(value: int) -> str:
    units = ["B", "KiB", "MiB", "GiB"]
    number = float(value)
    for unit in units:
        if number < 1024 or unit == units[-1]:
            return f"{number:.2f} {unit}" if unit != "B" else f"{int(number)} B"
        number /= 1024
    return f"{value} B"


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Report tracked repository files above size thresholds.")
    parser.add_argument("--warn-bytes", type=int, default=5 * 1024 * 1024)
    parser.add_argument("--max-bytes", type=int, default=10 * 1024 * 1024)
    parser.add_argument("--base-ref", default="origin/main")
    parser.add_argument("--fail-scope", choices=["added", "all", "none"], default="added")
    parser.add_argument(
        "--allow-glob",
        action="append",
        default=[],
        help="Approved glob for large files. Repeatable. Defaults cover tiny fixtures only.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    allow_globs = DEFAULT_ALLOW_GLOBS + list(args.allow_glob)
    files = tracked_files()
    added = added_files(args.base_ref) if args.fail_scope == "added" else set()

    rows: list[tuple[str, int, bool, bool]] = []
    failures: list[tuple[str, int]] = []
    for file_name in files:
        path = Path(file_name)
        if not path.exists() or not path.is_file():
            continue
        size = path.stat().st_size
        is_allowed = allowed(file_name, allow_globs)
        is_added = file_name in added
        if size >= args.warn_bytes:
            rows.append((file_name, size, is_added, is_allowed))
        should_fail = size > args.max_bytes and not is_allowed
        if args.fail_scope == "added":
            should_fail = should_fail and is_added
        elif args.fail_scope == "none":
            should_fail = False
        if should_fail:
            failures.append((file_name, size))

    if rows:
        print("Tracked files above warning threshold:")
        for file_name, size, is_added, is_allowed in sorted(rows, key=lambda item: item[1], reverse=True):
            flags = []
            if is_added:
                flags.append("new")
            if is_allowed:
                flags.append("allowed")
            flag_text = f" ({', '.join(flags)})" if flags else ""
            print(f"  {human_bytes(size):>12}  {file_name}{flag_text}")
    else:
        print("No tracked files exceed the warning threshold.")

    if failures:
        print("\n[ERROR] Unapproved tracked files exceed the maximum size threshold:")
        for file_name, size in failures:
            print(f"  {human_bytes(size):>12}  {file_name}")
        print(
            "\nLarge manuscript-scale data should be archived in a versioned release/Zenodo record "
            "with SHA-256 checksums committed to Git."
        )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
