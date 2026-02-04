#!/usr/bin/env bash
set -euo pipefail
tree -a -L 4 -I ".git|.venv|__pycache__|*.pyc|*.pyo|*.egg-info|.pytest_cache|.mypy_cache|dist|build|work|logs|.ruff_cache|node_modules"
