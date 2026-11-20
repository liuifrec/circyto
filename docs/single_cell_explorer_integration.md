# Single-cell Explorer software integration

Validated locally on 20 November 2026. The software baseline is
`fix/external-detector-distribution@6c280426ebf10abcbcaeab72f87c9d34060d6ee8`.
Only `faf1996033b58abc15a0d139cc16f8cfdeb263d3` was cherry-picked from
`feature/single-cell-explorer`, whose parent is `8f729b3`. The resulting feature
commit is `19f6758`; no intervening analysis or manuscript commits were merged.
The package version remains **0.10.0**.

Both worktrees had no tracked edits. The source worktree was fully clean; the
software worktree contained 26 existing untracked files. Their paths did not
collide with the feature, and all 26 hashes were preserved. No files were
stashed, moved or deleted. The temporary integration branch was created in the
existing software worktree; no additional worktree was needed.

## Compatibility

The Explorer runtime and its four reporting modules are unchanged from the
accepted feature commit. Existing detector, matrix, alignment, multimodal and
scientific-analysis code stays at the software baseline. The manifest conflict
was resolved without importing the research recipe entries or the curation TSV.
The optional KYSE180 replay now takes explicit `--source` and
`--candidate-metadata` paths, so it can read accepted external inputs without a
research checkout being part of the software distribution. The CLI and public
reporting APIs retain their existing default behavior.

## Verification

| Check | Result |
| --- | --- |
| Full package suite, Python 3.10.20 | **517 passed, 10 skipped** |
| Critical acceptance-module execution gate | **13 modules passed** |
| Installed-wheel reporting, importer and multimodal suite | **174 passed, no skips** |
| Explorer tests within both suites | **37 passed** |
| Existing QC tests within both suites | **55 passed** |
| Default QC comparison | **8 cases, 16 files byte-identical**, same normalized CLI output |
| Offline synthetic and KYSE180 replay | HTML and JSON byte-identical to accepted feature outputs |
| Browser smoke | **5 cases, 92 checks passed**, no HTTP requests |
| Installed CLI, doctor, native import and package-content checks | Passed |
| Wheel and sdist | Built; unresolved Java assets and research materials excluded |

The ten package skips are six explicitly gated external-detector integrations,
two optional Scanpy tests, and two missing-MuData negative tests because MuData
is installed. No cohort, alignment, counting or detector analysis was rerun.
The package command is `CIRCYTO_SKIP_INTEGRATION=1 python -m pytest -q tests`;
this runs all package tests without collecting other local analysis worktrees.
The installed suite runs outside both checkouts and verifies that every loaded
`circyto` module comes from the new wheel environment's `site-packages`.

The wheel was installed offline into a new venv using existing environment
dependencies. User-site packages were disabled; an unrelated user-site Kaggle
installation otherwise caused `pip check` failures. With that isolation,
`pip check` passes. This is not a fresh dependency-resolver validation. Source
tests used SciPy 1.15.3; installed tests used 1.15.2. Both used AnnData 0.11.4,
MuData 0.3.10, NumPy 2.2.6, pandas 2.3.3 and Typer 0.27.0.

The KYSE180 replay checks all **181 irradiated cells**, six RNA genes and seven
native junctions against the accepted objects by exact cell ID. Coordinates,
stored expression, support, explicit cohort order, both completed-zero cells
and condition summaries are unchanged. The selected input objects, source QC
and curation table retain their SHA-256 hashes. Synthetic and unit tests
separately verify failed, missing, unknown and unprocessed states and their
exclusion from completion denominators.

Chrome checks cover synthetic and real desktop reports, a 500-pixel narrow
view, missing UMAP and the combined QC/Explorer document. The narrow and final
two cases use real completion polling over loopback because virtual-time DOM
capture returned before asynchronous exports completed. The accepted renderer
was not changed. HTTP name resolution was disabled, and no HTTP/CDN resource
requests occurred. Browser timings are not performance benchmarks.

## Preservation and resources

Jobs ran serially under the existing bounded runner: `memory.max=4294967296`,
`memory.swap.max=0`, and CPU affinity **0–1** for every observed descendant.
The CPU quota controller is not delegated on this workstation; the two-CPU
bound is enforced by verified affinity. No OOM events occurred. Initial
availability was about 26 GiB RAM and 233 GiB disk. Recorded service peaks for
the substantive suite, installation, replay and browser checks stayed below
0.4 GiB.

The nine tracked manuscript/review entries, all 26 original untracked files,
the accepted scientific inputs and every pre-existing branch other than the
requested software branch remain unchanged. This includes the SRRC manuscript.
No push, worktree deletion, history rewrite, publication or submission occurs
as part of this integration.

The machine-readable [verification record](examples/explorer/integration_20261120/verification.json)
contains the complete changed-file list, test/skip counts, module hashes,
preservation checks and compact resource receipts. The
[package JUnit](examples/explorer/integration_20261120/package-tests.xml) and
[installed-wheel JUnit](examples/explorer/integration_20261120/installed-tests.xml)
retain individual test results. Detailed local logs, replay HTML, browser
receipts and build artifacts are under
`/tmp/circyto-explorer-integration-20261120`.
