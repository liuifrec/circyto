mkdir -p tools docs
cat > tools/gen_repo_fingerprint.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

mkdir -p docs

{
  echo "# Circyto repo fingerprint (generated)"
  echo "# Generated at: $(date -Is)"
  echo

  echo "## Git"
  echo "TOPLEVEL: $(git rev-parse --show-toplevel 2>/dev/null || echo NA)"
  echo "BRANCH:   $(git branch --show-current 2>/dev/null || echo NA)"
  echo "COMMIT:   $(git rev-parse HEAD 2>/dev/null || echo NA)"
  echo "STATUS:"
  git status -sb 2>/dev/null || true
  echo
  echo "REMOTE:"
  git remote -v 2>/dev/null || true
  echo

  echo "## Python"
  echo "PYTHON:   $(python -V 2>&1)"
  echo "WHICH_PY: $(command -v python || echo NA)"
  echo "VENV:     ${VIRTUAL_ENV:-NA}"
  echo

  echo "## Installed CLI resolution"
  echo "WHICH_CIRCYTO: $(command -v circyto || echo NA)"
  echo -n "CIRCYTO_MODULE: "
  python -c "import circyto; print(circyto.__file__)" 2>/dev/null || echo NA
  echo

  echo "## Key tools (optional)"
  for t in bwa samtools; do
    printf "%-10s %s\n" "${t^^}:" "$(command -v $t 2>/dev/null || echo NA)"
  done
} > docs/REPO_FINGERPRINT.txt

echo "[ok] wrote docs/REPO_FINGERPRINT.txt"
EOF

chmod +x tools/gen_repo_fingerprint.sh
./tools/gen_repo_fingerprint.sh
