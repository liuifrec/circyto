# Reproducing the A100 status replay

These are execution records, not permission to repeat raw analysis. The original
and repair budgets are separate; a future replay needs an unexpired approved
window and a fresh destination. The helper refuses existing workflow/import
outputs. Never run the one-cell or four-cell workflow into these directories.

Bindings used for the completed replay:

```sh
A02_ROOT=/home/liuyuchen/data/circyto_cmpb_a02_20261015
REPO=/tmp/circyto-cmpb-evaluation
REPAIR_ROOT="$A02_ROOT/status_revalidation_20261016"
CANDIDATE_ROOT="$REPAIR_ROOT/final_candidate"
REPAIR_COMMIT=e6843a38469d49f96d44758a4e37893ebbe2734b
```

The verified base is `b9726b2a9ab2cffcacc7bd22a2eae5f969b09234`. Its original
failed acceptance is unchanged. The production-only diff is:

```sh
git -C "$REPO" diff b9726b2a9ab2cffcacc7bd22a2eae5f969b09234 \
  "$REPAIR_COMMIT" -- circyto/pipeline/run_ciri3.py tests/test_alignment_status_handoff.py
```

All managed launch specs and logs reside under `$REPAIR_ROOT/specs` and
`$REPAIR_ROOT/logs`. The machine-readable
[commands.json](status_revalidation/commands.json) retains exact argv and hashes
with `${A02_ROOT}` and `${REPO}` placeholders replacing workstation paths. The
original literal specs remain on disk. Every spec declares `category=firstcell`
and `footprint_root=$A02_ROOT`; all stages share the same two-hour ledger and the
same total evaluation storage cap. `run_stage.py` applies the inherited 20 GiB
process-tree, CPU-affinity, wall-time and disk controls.

The final build command, executed inside the managed stage, was:

```sh
"$REPAIR_ROOT/test-env/bin/python" "$REPO/scripts/cmpb_a02/build_status_wheel.py" \
  --root "$CANDIDATE_ROOT" --repo "$REPO" --commit "$REPAIR_COMMIT"
```

It stages only versioned build inputs from that exact commit, builds wheel/sdist
in isolation, inspects every distribution, rebuilds a wheel from the sdist, and
creates `$CANDIDATE_ROOT/installed-env`. Dependencies are installed from
[constraints.txt](status_revalidation/constraints.txt), then the exact local wheel
is installed with `--no-deps`. No editable install or source checkout import is
used for installed acceptance or A100 replay. The staging tar is a local build
input; only the inspected wheel/sdist are distribution candidates. No uploads,
registry route assertions, tags or version changes are part of this record.

The actual postprocessing command was:

```sh
"$CANDIDATE_ROOT/installed-env/bin/python" \
  "$REPO/scripts/cmpb_a02/replay_alignment_status.py" \
  --original "$A02_ROOT" --root "$CANDIDATE_ROOT"
```

It is launched through `run_stage.py launch`, not directly without safeguards.
Before production functions run, the helper verifies installed member hashes,
blocks process creation (including BWA/STAR/Java/indexing), checks all original
bound evidence, and creates independent copies. It aggregates genuine child
records through the fixed production writer and uses the real QC/object/report
functions. The existing independent CSV/exact-integer/CSR checker is unmodified.
Both original and regenerated object contents are compared; only alignment
status/duration fields may differ in workflow AnnData.

The final default-suite command is `python -m pytest -q .` with
`CIRCYTO_SKIP_INTEGRATION=1` and a JUnit output path, followed by the existing
`.github/scripts/check_acceptance_junit.py`. Its exact source-path bootstrap and
environment are retained in the stage spec. Both successful and unsuccessful
attempts have distinct stage names and retained logs.

Actual report locations:

- Original failed workflow: `$A02_ROOT/workflow_firstcell/qc/report.html`.
- Revalidated workflow: `$CANDIDATE_ROOT/workflow_firstcell/qc/report.html`.
- Byte-preserved native import: `$CANDIDATE_ROOT/import_firstcell/qc/report.html`.
- Offline browser previews and DOM audit: `$CANDIDATE_ROOT/browser_attempt2/`.

Chrome's first temporary socket path was too long. For the successful offline
retry, `/tmp/circyto-a02-r` pointed to `$REPAIR_ROOT/tmp`, so temporary storage was
still persistent and included in the shared disk watchdog. The browser opened
byte-verified copies of both HTML files with DNS disabled. No workflow report
command was applied to the import bundle.

Native/raw/large scientific files remain outside Git. Per-candidate audit values
are in [native_value_checks.tsv](status_revalidation/native_value_checks.tsv) and
[report_value_checks.tsv](status_revalidation/report_value_checks.tsv). Source
hashes, complete resolved dependencies, output hashes and limits accompany them.
This validates preserved A100 evidence only; A101–A103 remain not started.
