from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("circyto")
except PackageNotFoundError:
    __version__ = "0.9.0"

__all__ = ["__version__"]
