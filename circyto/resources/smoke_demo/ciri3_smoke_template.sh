#!/usr/bin/env bash
set -euo pipefail

alignment="${1:?alignment path required}"
raw_output="${2:?raw output path required}"
cell_id="${3:?cell id required}"

support="$(samtools view -c "$alignment" || echo 0)"
mkdir -p "$(dirname "$raw_output")"
cat > "$raw_output" <<EOF
circRNA_ID	chr	start	end	strand	bsj_reads
${cell_id}_smoke	chr21	100	200	+	${support}
EOF
