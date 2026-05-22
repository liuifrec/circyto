#!/usr/bin/env bash
set -euo pipefail

circyto summarize-benchmark \
  --workdirs /user/ifrec/liuyuchen/circyto_redo/emtab8735/work/diySpike_workflow_all192 \
  --workdirs /user/ifrec/liuyuchen/circyto_redo/scrr_imr90/work_hg38_fullstack \
  --workdirs /user/ifrec/liuyuchen/circyto_redo/scrr_hap1/work_hg38_fullstack \
  --output benchmark_summary.tsv \
  --json benchmark_summary.json
