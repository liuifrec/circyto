# Pre-release integration review

Date: 2026-10-13. Review baseline: main
`44697355bcab1c525ca7ef9b130e2ad0094d9e1b` through native feature
`9a80ec537aed2a81495c2e14b58f9f9c20e51375` (56 files, 5,976 insertions, 8 deletions).
Integration fixes are on `review/pre-release-integration`, leaving all original
branch tips unchanged. Obtain the integration commit with
`git log -1 --format=%H -- docs/pre_release_readiness.md`; PR checks identify their
own exact tested head/merge ref. Version remains 0.10.0, with no release or merge.

**GO for draft formal code review; NO-GO for public release.** RERF internal review
is independent of this feature branch and these release gates. The frozen
manuscript is not revised and no new novelty or biological claim is authorized.

## Review of the cumulative difference

Review covered the complete changed-file list: CLI registration/options, two
workflow hooks, reporter/style, strict import and native adapter, tests/examples,
fixture bytes/provenance, packaging metadata and all added guides/records.
No calling algorithm or biological threshold changed. Existing exporters and
historical objects are unchanged. The integration review found:

| Finding | Severity / disposition |
| --- | --- |
| Base CI installed no MuData, permitting importer/export tests to skip | Review blocker fixed: install dev extra, import dependency preflight, JUnit audit of nine critical modules; Python 3.10/3.12 matrix |
| Workflow `report` overwrote import-circ's HTML and machine-readable provenance because both use qc/metrics.json | Fixed: refuse recognized `circyto.import_circ.*` bundles with actionable instructions; regression checks every bundle file hash |
| Native pooled tables may have identical headers; genome, producer version and cell mapping cannot be inferred from syntax | Explicit user assertions remain necessary; genuine per-cell + matching RNA gate pending |
| Workflow status summaries lack independent cell/result binding and aggregate candidate counts cannot be subset by completion | Conservative whole-cohort denominator; report is source interpretation, not an audit of matrix values or detector success |
| Reports now run by default after scientific outputs are persisted; report failure makes the command fail despite retained results | Documented behavior; `--no-report` opt-out and rerun instructions; workflow regression coverage present |
| Optional MuData and broad dependency ranges; older outputs have no additive observation-status fields | Local 3.10 stack and remote matrix tested separately; older objects must remain unknown, not be reclassified from zero counts |
| Native parser reads raw/normalized tables in memory; exports keep sparse matrices | No dense cell-by-candidate construction found; sparse guard test present; large-cohort throughput/memory acceptance remains unmeasured |
| Source copies retain metadata/path strings; copied manifest paths are original and may not be relocatable | Documented evidence-versus-runnable-manifest boundary; review before sharing |
| Minimal/minimum dependency stack, other OSes, Python versions beyond CI matrix and actual external detector execution are not certified | Compatibility limits, not evidence supplied by green default CI |

The importer reuses CIRI3 normalization, the strict long-table CSR builder,
host-gene schema normalization, and MuData synchronization/writing. Direct
AnnData creation intentionally preserves int64 support; routing through legacy
int32/coercing exports would weaken that contract. Report HTML uses shared CSS,
but import and workflow report schemas/denominators differ intentionally:
imports can restrict a matrix to explicitly completed cells; workflow aggregate
QC cannot. Combining these scientific semantics merely to remove similar
rendering code would be unsafe. Generic TSV/ID/integer validation is shared by
the native adapter; its additional strict format checks are not a second import
pipeline. No additional detector parser or statistical model was introduced.

The recognized-schema overwrite guard does not authenticate a bundle whose
metrics were removed or corrupted. Source bundles should remain immutable;
report regeneration for imports requires rerunning from original inputs into a
new directory. Full native evidence fields stay in source copies; only selected
candidate annotations become `var` columns.

## Local validation in this integration pass

- Full default suite: **475 passed, 10 skipped, 9 warnings**, 19.19 seconds,
  Python 3.10.20, `CIRCYTO_SKIP_INTEGRATION=1 python -m pytest -q .`.
- Critical modules: general import 40, native import 33, QC report 55, MuData
  export 5, downstream 2, synchronization 8, multimodal export 1, full-length
  workflow 19 and Smart-seq3 workflow 12 executed. Two missing-MuData negative
  tests skip only because MuData is installed. No success-path skip is allowed.
- Remaining skips: six external-detector integration tests and two Scanpy tests.
  Warnings include AnnData index/duplicate-name and MuData future behavior notices.
  These are not new detector runs or proof of detector accuracy.
- The JUnit gate was exercised against deliberate missing-module, skipped-import
  and failed-test XML variants; each fails. It rejects absent/failed critical
  modules and unexpected skips, rather than trusting a green pytest exit alone.
- Wheel/sdist built locally; isolated wheel installation and examples run outside
  the repository. The PR validation record supplies installed acceptance results
  and exact build/CI identifiers; previous milestone measurements remain linked
  below and are not presented as a fresh biological benchmark.

Remote CI is a separate acceptance record: the draft PR runs the complete suite
on Python 3.10 and 3.12 with MuData, audits JUnit, builds wheel/sdist, installs the
wheel in a separate venv and runs the authentic native demo from outside the
checkout. CI uploads XML/dependency/acceptance records only, no detector binaries.
Consult the PR checks and retained artifacts for actual outcomes, not this
workflow description as a claim that CI passed. External detector integrations
remain intentionally gated; real paired-cell acceptance remains outside CI.

## Release-readiness matrix

| Category | Evidence / state | Remaining gate |
| --- | --- | --- |
| Validated implementation | Strict long/matrix/native import, CSR int64 preservation, statuses, optional MuData, HTML/JSON, overwrite protection; local suite above | Formal code review and green remote checks on exact PR head |
| Additional acceptance required | Authentic upstream native sample: 6 candidates / 50 support; technical RNA test; simulation annotation contract | Genuine per-cell native CIRI3 + same-cell RNA, independent reviewer and representative memory/runtime measurements; use the [checklist](external_ciri3_acceptance.md) |
| Additional acceptance required | Synthetic workflow report cases, failure/denominator regressions | Original completed public-data workflow summary + both QC TSVs for real HTML acceptance |
| Third-party distribution/legal review | CIRI3 source/notice route incomplete; CIRI-full v2.0 permission unresolved; exact binaries identified | Reviewed packaging migration or completed permission/source/notice obligations; [narrow recommendation](../THIRD_PARTY_NOTICES.md#pre-release-packaging-recommendation-2026-10-13) |
| Manuscript evidence and public archive | Checksum-pinned historical Smart-seq3/IMR90 object round trips; frozen results unchanged | Durable processed-data deposit/access instructions and manuscript release baseline decision; no inferred new biological evidence |
| Remaining release blockers | Review, distribution, real-data acceptance and archive gates are distinct | Do not equate CI success with public-release clearance; no version bump, merge, tag or publication in this task |

Historical processed-object round trips validate 192 × 2,503 / 2,659 nnz / 5,711
support (Smart-seq3) and 23 × 2,443 / 2,914 nnz / 7,422 support (IMR90). They are
not independent native CIRI3 inputs. Detailed [import evidence](import_circ_validation.md),
[native evidence](native_ciri3_validation.md), [report evidence](qc_reporting_validation.md)
and [artifact search/checklist](external_ciri3_acceptance.md) retain that boundary.
The additions strengthen reusable representation and auditable observation
semantics, not a new circRNA detection method. Any eventual manuscript novelty
revision requires a separate decision and evidence review.

All 72 protected manuscript recheck entries (including the directory itself)
are checked for bytes, size, mode and mtime. No scientific dataset, Figure 1A,
manuscript file, original local branch tip or main/frozen branch is changed.
