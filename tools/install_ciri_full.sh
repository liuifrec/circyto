#!/usr/bin/env bash
set -euo pipefail
echo "This script installs third-party CIRI-full v2.0 into tools/. It requires git, perl, java, an aligner (e.g., BWA), etc."
echo "It is NOT run automatically by pip. Use at your own discretion."
# Example scaffold; replace with your licensed source:
# git clone https://.../CIRI-full_v2.0 tools/CIRI-full_v2.0
# chmod -R a+rX tools/CIRI-full_v2.0
echo "Done. Wrapper will use tools/CIRI-full_v2.0 automatically."
