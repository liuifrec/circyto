# QC reporting validation

Validation date: 2026-10-09. Feature branch: `feature/qc-report-usability`.
This records the initial feature pass; the later
[reviewer-readiness audit](reviewer_readiness.md) supersedes its acceptance
assessment and documents a manifest-example validation correction.
Base: `44697355bcab1c525ca7ef9b130e2ad0094d9e1b` (v0.10.0). The package version
is unchanged. The manuscript branch remains at
`c29831c3f6e7a31b708cce493c0be8036176077a`; no merge or manuscript edits were made.

## Automated tests

The complete default suite passed: **396 passed, 10 skipped, 5 warnings**.
The command ran with the existing `circyto` Python environment activated:

```bash
CIRCYTO_SKIP_INTEGRATION=1 python -m pytest -q .
```

Test environment: Python 3.10.20, pytest 9.1.1, NumPy 2.2.6, pandas 2.3.3,
SciPy 1.15.3, AnnData 0.11.4, MuData 0.3.10, and Typer 0.27.0.

This is the repository's CI-style gate: external detector integrations remain
gated, and optional capability tests can skip depending on installed extras.
Warnings were from existing AnnData/MuData tests (index conversion/nonunique
variable names). No raw analysis or manuscript evidence regeneration was run.
The initial invocation without an activated environment failed shell-based
tests because `python` was absent from `PATH`; activating the environment
resolved those failures without changing those tests or scripts.

The reporting tests cover exact integer counts (including values above float
precision), source QC/matrix agreement, missing columns and values, empty
tables, malformed JSON/TSV, duplicate identities, failed/incomplete/dry-run
status, contradictory summaries, source checksums, deterministic regeneration,
privacy omissions, HTML escaping, offline assets, paths containing spaces,
unsafe symlink targets, source aliases, hardlinks, and write errors. The
documentation example is checked against its Matrix Market file and committed
metrics. Both workflow tests exercise automatic reporting, `--no-report`, and
simulated reporting errors after scientific outputs have been persisted.

## Executable examples and visual review

- Synthetic report generation through the actual CLI passed: 12 cells,
  6 candidates, 77 support, median 2 candidates/cell, 2/12 zero-count cells,
  and 4/6 recorded host-gene coverage. Every per-cell/per-candidate count agrees
  with both the source TSVs and the synthetic matrix.
- The later independent-user check found that the example's legacy manifest
  columns were rejected by `manifest validate --strict`, despite being accepted
  by the workflow. The reviewer pass corrects the example and verifies both.
- The real `workflow full-length-circrna --dry-run` command passed using the
  generated toy manifest and reference, with installed BWA/SAMtools on `PATH`.
  CIRI3 runtime was unavailable for detector-command preview; no detector ran.
  Standalone reporting on the plan correctly reported `dry_run` with unavailable
  detection metrics.
- The 1440 × 1400 Chrome screenshot was rendered from the local HTML and
  visually reviewed. The synthetic-data warning is visible before the status
  and metrics. CSS and SVG are embedded; there are no script or remote-asset
  dependencies. Screenshot reproduction is documented in [the guide](qc_reporting.md).

No completed public-data workflow directories were available in the current
workspace/local data roots. Validation therefore uses hermetic fixtures and
existing workflow tests, not a claim of revalidation on manuscript datasets.

## Distribution checks

Wheel and source distribution builds passed (build 1.5.0, setuptools 83.0.0,
wheel 0.47.0). A fresh virtual environment without system site packages installed
the wheel and all declared dependencies; `pip check` passed. Checks ran from
`/tmp`, outside the repository, and confirmed that Python imported the installed
package from that virtual environment. Installed `--version`, top-level help,
`report --help`, both workflow help commands, legacy `run --help`, and `doctor`
all passed. The installed CLI regenerated the demonstration HTML and metrics
byte-for-byte identically to the committed example. Clean-install Typer was
0.27.3, so CLI checks also cover a newer resolver-selected version.

Archive inspection confirmed the report modules and CLI entry point in the
wheel, and all demo/guide assets in the source distribution. `git diff --check`
passed. The build command is:

```bash
python -m build --no-isolation --outdir /tmp/circyto-qc-packages
```

The package embeds styling in Python modules, so it needs no HTML/CSS asset
loader or new dependency. The source distribution also includes the demo
generator, source fixture, report, metrics, screenshot, and reporting guide.
Build artifacts keep version `0.10.0` and must not be published as a replacement
for the existing v0.10.0 release.

## Preservation and submission decision

The initial tracked worktree was clean on the manuscript branch. That branch
and the main baseline were preserved before creating the feature branch from
`origin/main`. Four local duplicate manuscript figure files in
`manuscript/results_v2_recheck/` were ignored on the manuscript branch and became
visible as untracked files after switching to main's ignore rules. They remain
untouched and are excluded from the feature commit. BioRender Figure 1A and
all scientific figures/analyses are outside this change.

Recommendation: **post-review usability enhancement**. The feature is ready for
code review and laboratory-facing usability evaluation, but should not be
represented as part of the frozen Application Note software baseline. Including
it in the upcoming submission would require a separately reviewed feature
snapshot and explicit manuscript wording about its scope; it is not new
scientific evidence.
