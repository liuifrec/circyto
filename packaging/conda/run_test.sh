#!/usr/bin/env bash
set -euo pipefail
unset PYTHONPATH
export PYTHONNOUSERSITE=1
# MuData is a test dependency; mulled tests in meta.yaml use only run dependencies.
python -I examples/installed_release_check.py --fixtures examples/native_ciri3 --out acceptance --multimodal
