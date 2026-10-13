# Application Note reviewer-readiness audit

Date: 2026-10-09. Branch: `feature/qc-report-usability`.
Starting commit: `a0d5f35602235c8b5b59bcd908190c23d2ea29ff`.
The audited changes are in the commit introducing this file; obtain its exact
SHA with `git log -1 --format=%H -- docs/reviewer_readiness.md`.
Package version remains 0.10.0. This is not a release or a change to the frozen
manuscript software baseline.

## Decision and remaining acceptance gate

**Keep reporting as a post-review usability enhancement for now.** It is a
candidate for the paper's eventual software release after original public-data
workflow report acceptance and the existing detector distribution gaps below
are resolved. This experiment does not gate RERF internal-review initiation.

No accessible completed Smart-seq3, IMR90 scRR, or HAP1 scRR workflow directory
with the original summary and both QC TSVs was found. The search covered local
home/data/work and temporary roots, available mounts, relevant development
archives/source distributions, and tracked filenames across local Git history.
The available report inputs were synthetic; other summaries were dry runs.
Historical server directories documented in
[the RamDA run record](ramda_full_stack_server_run.md) are not mounted here.
GitHub Actions returned zero retained artifacts; the listed releases supplied
no downloadable assets. This does not establish that the originals no longer
exist elsewhere.

Consequently **real-workflow HTML acceptance remains unperformed**. Historical
processed H5MU objects are separate manuscript evidence and were not converted
into pretend workflow summaries/QC. No detector, alignment, raw-data analysis,
or manuscript regeneration ran during this pass. Runtime `-H`/version checks
and a planning-only workflow are identified separately below.

## Changes supported by independent testing

- Fixed a scientific-denominator defect: failed, unprocessed, missing-status,
  conflicting, and dry-run cohorts cannot supply candidate prevalence or
  zero-count rates. All original QC values remain visible. Rates require the
  whole QC cohort to have completed alignment/detector evidence; aggregate
  candidate counts cannot be restricted to a successful-cell subset.
- Separated malformed optional RNA summaries (warning) from unreadable detector
  evidence (completion cannot be established). No calling, matrix, annotation,
  cell-identity or object semantics changed.
- Corrected the demo's legacy manifest columns: the workflow accepted them but
  the documented strict validator rejected them. The canonical example now
  passes both. Added the missing BWA indexing instruction and explained that
  a wheel does not install the source-archive demonstration script.
- Added a tested optional external-tool recipe, a five-question
  [reviewer guide](reviewer_guide.md), valid [software citation](../CITATION.cff)
  based on the existing James Liu authorship metadata (no invented DOI),
  [third-party notices](../THIRD_PARTY_NOTICES.md), and a
  [primary-source comparison](circrna_workflow_comparison.md).

## Independent installation and runtime checks

A fresh Python 3.10.20 virtual environment outside the repository installed
the built wheel without editable imports or system-site packages. All six
requested commands exited 0: `circyto --help`, `doctor`, `detectors`,
`manifest --help`, `workflow full-length-circrna --help`, and `report --help`.
`pip check` passed; imports came from the new environment's `site-packages`.
Each command took about 0.42–0.47 seconds. Wheel/dependency installation took
7.55 seconds using a warm pip cache; this is not an uncached download benchmark.

The lightweight install correctly reported core Python readiness alongside
missing BWA/STAR and an unsuitable system Java 8. Reports still worked.
Python installation does not install or validate references/external tools.

The separate Linux x86_64 micromamba recipe installed BWA 0.7.19-r1273,
SAMtools 1.22.1 (HTSlib 1.24), and OpenJDK 17.0.18. Initial external environment
creation took 19.46 seconds with 81 packages / 276 MB planned downloads; a
fresh cached recreation with user-site isolation took 4.33 seconds. The first
conda Python attempt exposed unrelated user-site packages; the delivered recipe
sets `PYTHONNOUSERSITE=1`. Its recreated environment passes `pip check` and
imports all dependencies locally. BWA's no-argument help returns its expected
exit 1; SAMtools/Java version commands and bundled CIRI3 `-H` return 0.
`doctor` marks CIRI3 ready. Strict toy-manifest validation and the documented
single-end workflow `--dry-run` pass, with a Java/BWA command preview.

This tests setup and planning, **not new biological end-to-end performance**.
BWA and Java >=12 are required for the single-end direct-SAM route; Java 17
was tested. SAMtools supports BAM handling; STAR is needed for the established
paired-end routes. Real execution also needs a matching FASTA/GTF, prepared BWA
indices, real per-cell FASTQs, and adequate storage/memory. Dry-run success does
not check BWA indices. The recipe pins key tools but is not a complete lockfile
and has not been tested on macOS/Windows.

## Verification results

| Check | Result |
| --- | --- |
| Full default suite: `CIRCYTO_SKIP_INTEGRATION=1 python -m pytest -q .` | **402 passed, 10 skipped, 5 warnings**, 17.52 s |
| Warnings / skips | Existing AnnData string-index conversion and MuData nonunique-variable warnings; external/optional integrations remain gated. No warnings suppressed. |
| Wheel + sdist: `python -m build --no-isolation` | Both build successfully; CLI, report modules, notices, and distributed demo/guide/recipe assets verified. No manuscript recheck files or H5MU inputs packaged. |
| Installed-only demo | Source-archive generator + installed wheel; 12 cells, 6 candidates, 77 support, median 2, zero-count fraction 2/12, host-gene coverage 4/6. Every cell/candidate count and prevalence numerator checked against TSV; regression tests also compare the synthetic matrix. |
| Determinism and preservation | HTML/JSON byte-identical across regeneration and equal to the distributed sample; source SHA-256 values, bytes, sizes, mtimes and modes unchanged. |
| Installed failed/missing-cell examples | Raw zero counts preserved; status failed/incomplete respectively; rates and denominators unavailable, with explanatory HTML. These are synthetic diagnostic variants. |
| Browser | Local `file://` report rendered in Chrome with hostname resolution disabled; no scripts/CDN assets. Refreshed and visually inspected 1440 × 1400 [screenshot](images/qc_report.png). |
| Citation | `cffconvert --validate --infile CITATION.cff`: valid CFF 1.2.0. Validator is test tooling only. |
| Git/data guards | `git diff --check` passes. All 71 entries under `manuscript/results_v2_recheck/` retain exact hashes/metadata; none staged. No scientific files or large inputs added. |

The Python test stack was pytest 9.1.1, NumPy 2.2.6, pandas 2.3.3, SciPy 1.15.3,
AnnData 0.11.4, MuData 0.3.10 and Typer 0.27.0; clean installation resolved
Typer 0.27.3. Local logs/artifacts are under `/tmp/circyto-review/` and are not
required by the installed package or committed documentation.

## Release and reproducibility concerns

The existing bundled CIRI-full v2.0 JAR lacks a verified applicable license.
CIRI3's supplied GPLv2 notice is retained, but corresponding-source distribution
and embedded HTSJDK notices need completion before a release. The audit
documents exact binary identity; it does not certify licensing completeness
or silently apply circyto's MIT license to third-party software.

The frozen manuscript reproduction route requires two exact historical H5MU
objects and a pinned scientific stack. The manuscript guide provides Git
archive recovery, but a dedicated processed-data deposit remains pending.
Reporting reads summary/QC files only; it cannot verify reference content from
basenames, validate matrix contents independently, or establish biological
truth. Those boundaries are explicit in the guide and report.

Changed files: `circyto/pipeline/qc_report.py`, `tests/test_qc_report.py`,
`examples/qc_report_demo.py`, `README.md`, `MANIFEST.in`, `pyproject.toml`,
`CITATION.cff`, `THIRD_PARTY_NOTICES.md`, `environment.ramda-se.yml`,
`docs/qc_reporting.md`, `docs/qc_reporting_validation.md`,
`docs/reviewer_guide.md`, `docs/reviewer_readiness.md`,
`docs/circrna_workflow_comparison.md`, and the synthetic example's
`qc/report.html`, `qc/metrics.json`, and `docs/images/qc_report.png`.
The main baseline remains `44697355bcab1c525ca7ef9b130e2ad0094d9e1b`; the
manuscript branch remains `c29831c3f6e7a31b708cce493c0be8036176077a`.
