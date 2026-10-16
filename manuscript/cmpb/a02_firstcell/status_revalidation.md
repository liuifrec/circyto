# A100 alignment-status repair and derived-output revalidation

**A100 STATUS REPAIR REVALIDATED — scientific values unchanged; remaining cells
not started.** This is a separate record. The original failed acceptance and its four FAIL assertions
remain unchanged in `b9726b2a9ab2cffcacc7bd22a2eae5f969b09234` and on disk.

## Authorized amendment before replay

The verified evidence head is `b9726b2a9ab2cffcacc7bd22a2eae5f969b09234`, parent
`17b733c8f21bb128a035e7b100e704099e2bdb0e`. A new local branch,
`fix/workflow-alignment-status`, preserves all original branch tips. Evaluated
biological execution used `6c280426ebf10abcbcaeab72f87c9d34060d6ee8`.

The original launch deadline expired at 2026-10-15 22:14:54 UTC. The user explicitly
authorized a **separate two-hour repair/revalidation window**, conservatively
starting 2026-10-15 23:40:36 UTC and ending 2026-10-16 01:40:36 UTC. The original
budget is not reset. The new ledger is under
`$A02_ROOT/status_revalidation_20261016`; its disk watchdog charges the entire
`$A02_ROOT`, including retained originals, against the existing 80 GiB cap
(automatic stop at 79 GiB). Process-tree memory remains 20 GiB with swap disabled,
CPUs 0–7 and all-thread affinity monitoring, and the original user-service wall
and descendant termination controls. CPU quota is unavailable on this host.

Only A100 is eligible for status/QC/object/report replay. No alignment, detection,
indexing, RNA recounting, other cells or A03–A06 work is authorized here. Python
audit hooks reject external process creation and writes to original evidence;
tiny tests exercise those safeguards before the real replay. Scientific inputs
and outputs are hashed before/after. Mutable derived QC/HDF5 files are independent
copies, never hard links to originals. The registered six-file study design is
unchanged; this document records the narrow execution amendment.

## Root cause and repair boundary

The real per-aligner writer in `align_manifest.prepare_alignment_cache` records
A100 as `aligned`. `run_ciri3._write_alignment_summary` previously retained only
routing/planning fields. The full-length workflow's `_per_cell_record_map` and
`_build_execution_summary` therefore supplied no alignment record to
`workflow_reporting.build_cell_qc_table`; QC and AnnData received `missing`.
`qc_report.generate_report` correctly withheld completion and a denominator.

The production repair is confined to the combined writer. It preserves actual
child records and links, orders identities by the source manifest, derives status
counts from those records, and excludes the merged manifest from reaggregation.
Missing/aggregate-only legacy evidence does not establish per-cell completion;
record-only legacy evidence remains usable. Duplicate, foreign, conflicting
count/layout/reuse and dry-run-success evidence is rejected. Reused states remain
reused; dry-run cells remain planned. Report safeguards and all scientific
processing code are unchanged.

Existing full-length tests replaced `run_ciri3_workflow` with a fake that directly
wrote a correct combined summary. The new test uses the real per-aligner writer
with a tiny fake executable, then the real combined writer, QC/AnnData exporter
and report consumer. It failed before repair (`missing` versus `aligned`).

Local production commits:

- `bb1b8559d86c34ac955546e5381959f2f0593704`: preserve/reconcile child evidence.
- `e6843a38469d49f96d44758a4e37893ebbe2734b`: reject contradictory reuse metadata;
  this is the final artifact source.

The wheel retains local metadata `0.10.0` solely because this task forbids a
version bump. It is an unpublished repair artifact identified by commit and
SHA-256, not another frozen manuscript release. Neither original installation
nor original artifact is replaced.

## Original execution versus repaired derived outputs

The replay used the final installed wheel and genuine preserved child summary,
after checking the original exit records and tool configuration, rechecking
FASTQ/reference/JAR/native/SAM hashes and the alignment-to-detector binding. It invoked `_write_alignment_summary`,
`_build_execution_summary` (real QC and AnnData export), then `generate_report`.
It did not invoke the full raw-read workflow. Four attempted tool-launch probes
were blocked before execution; postprocessing made **zero** external execution
attempts. See [execution_block.json](status_revalidation/execution_block.json).

| Original failed assertion | Original observation | Repaired derived observation |
| --- | --- | --- |
| Recorded alignment versus report | `missing` | `aligned` — PASS |
| Recorded alignment versus workflow AnnData | `missing` | `aligned` — PASS |
| Bound completed workflow report state | `incomplete` | `completed` — PASS |
| Bound completed-cohort denominator | unavailable | 1 — PASS |

[before_after.tsv](status_revalidation/before_after.tsv) retains both outcomes.
Original files and all four FAIL records remain unchanged. Report denominators
use the unchanged contract; failed/unknown/missing controls remain non-completed.
The regenerated workflow report has no completion reasons or report warnings.

All **142 candidate identities and per-candidate integer values, 474 total support
units, coordinates, strands, feature/cell order, annotations, layers and other
scientific metadata are unchanged**. The workflow X remains int32; native-import
X remains int64. The independent native comparison permits only lossless integer
widening for workflow values. No migration or scientific transformation occurred.
The regenerated AnnData differs only in `obs.alignment_status` (`missing` to
`aligned`) and `obs.alignment_seconds` (empty to the recorded `495.261`). Its
SHA-256 is `e5b804d3733f222e6ff018ec872b994b8f37371e1f522abb7b63162d05e0d323`;
the original remains
`9f5c8176420230810cd752aacb7b7f709f3d658ff6a5463f70338ab50b21cc98`.
The circRNA QC TSV is byte-identical. See
[scientific_preservation.json](status_revalidation/scientific_preservation.json).

The genuine original A100 RNA profile (1 × 63,187) and native-import bundle were
copied unchanged and rechecked with the repaired installation, including ordinary
AnnData/MuData read/write/read. The separately generated same-run overlap RNA
files were also copied byte-for-byte; no RNA quantification was repeated. Its
historical reference identity remains unverified. The upstream six-candidate,
50-support example and synthetic RNA remain separately labeled installation
tests; neither substitutes for the genuine A100 checks.

The unchanged independent checker passed **149 native/RNA checks**, **479
numeric/HTML checks**, and the four original status assertions now pass separately
(483 report/status checks in the new table). The extra numeric/HTML assertion
checks the now-eligible candidate histogram; original assertions were not relaxed.
Both reports opened offline and all browser tables/cards/SVG values matched their
HTML. Original per-read ancillary native evidence remains outside Git.

## Build, tests, retained failures and resources

| Artifact, from `e6843a38469d49f96d44758a4e37893ebbe2734b` | SHA-256 |
| --- | --- |
| Wheel | `28c7cbde3c31cf90dd7c247487ba5f410c5b82f6f97d7a2281a0ce52912f7383` |
| Clean sdist | `e7b8ac5202bda683c9401fcefeb4d539f0442d8f7b12984c40f95835a0efcff7` |
| Wheel rebuilt from that sdist | `aeb90fcd4246051cf336ecf52c7be8d6b2b08b47f1300537ba002e8faf1149b1` |

All distributions passed `twine check`, archive integrity/exclusion inspection and
required-source/notice checks. The two wheels contain identical Python member
bytes. JARs, classes and the unresolved CIRI-full manual remain excluded. A fresh
environment outside the checkout passed `pip check`, installed-origin checks,
actual CLI help/import/report acceptance and native MuData loading. All 87
installed circyto Python members match the final wheel. The original evaluated
wheel and its 87 installed Python members are unchanged.

Final default suite: **523 passed, 10 skipped, 9 warnings**. Twenty new handoff
regressions passed, with real producer/exporter/report boundaries. The 10 skips
are six gated external-detector integrations, two optional Scanpy tests (Scanpy
absent), and two negative tests applicable only when MuData is absent. The critical
import/MuData/report/packaging skip audit passed. Warnings concern intentional
duplicate-name fixtures, MuData future defaults and empty-index string coercion.
Validation used Linux x86_64 / Python 3.10.22 and the original pinned scientific
stack; no new cross-platform claim is made.

Retained unsuccessful attempts, all resolved for this scoped acceptance:

- Expected original-writer regression: one failure, `missing` versus `aligned`.
- Initial new regression scaffold: two failures comparing entire path lists after
  deliberately changing repeated-manifest inputs; record/count assertions passed.
  The test was corrected to compare reconciled evidence and counts.
- First full-suite environment: six subprocess import failures because circyto
  was not installed in the isolated test environment. Installing the built wheel
  resolved them; dependencies and scientific assertions were unchanged.
- Expected reuse-conflict regression: one failure before adding the conflict
  guard, then passing in the final 20-test set.
- First offline browser attempt: Chrome's Unix socket path was too long. A short
  alias into persistent repair storage resolved it. The failed log and outputs
  remain; report HTML was not changed.

The original four A100 failures remain historical FAILs. There are **no unresolved
repair/revalidation checks**. All failed stage logs and hashes are referenced in
[tests_and_artifacts.json](status_revalidation/tests_and_artifacts.json).

Managed repair work took 176.13 s, including both builds, retained failed attempts,
tests and preservation checks; the replay itself occupied 33.98 s including source
hashing. The calendar interval from authorized start through the last managed
stage was 16 min 9.56 s, within the two-hour window. Peak sampled repair cgroup
memory was 7.06 GiB, summed process-tree RSS 1.21 GiB, and the entire retained
evaluation footprint was 20.49 GiB. The original 5,046.32 s of active A02 work is
retained separately; combined active stage time is 5,222.45 s. Sampling limits,
all stage outcomes and the original expired deadline are in
[resource_usage.json](status_revalidation/resource_usage.json).

## Preservation, review files and next decision

Before/after hashes and path sets match for all 219 files in the bounded original
input/reference/workflow/import/audit trees. All originally recorded branch tips,
the 21 pre-existing CMPB files, the actual protected recheck/blueprint file set and
both original H5MU hashes match their preserved records. Counts here describe
files actually compared, not assumed historical directory counts. See
[source_preservation.json](status_revalidation/source_preservation.json) and
[protected_work_preservation.json](status_revalidation/protected_work_preservation.json).

Changed production path: `circyto/pipeline/run_ciri3.py` only. Added regression:
`tests/test_alignment_status_handoff.py`. Analysis-only changes comprise a bounded
parent-directory disk-accounting extension to `run_stage.py`, the wheel-build and
replay helpers, and their six tiny safety tests. No reporting, detector, annotation,
RNA or count-matrix algorithms changed. The new audit packet adds no raw FASTQ,
SAM/BAM, HDF5, JAR, credentials, private logs or manuscript figure.

Exact argv, environment, hashes and helper snapshots are linked from
[commands.json](status_revalidation/commands.json),
[environment_manifest.json](status_revalidation/environment_manifest.json),
[output_artifacts.json](status_revalidation/output_artifacts.json), and
[reproduce_status_revalidation.md](reproduce_status_revalidation.md).
The original workflow timing is retained in the derived summary, with a separate
`postprocessing_replay` provenance block and source workflow UUID. New report
paths and regenerated metadata do not imply a second biological execution.

[pilot_cell_trace.tsv](status_revalidation/pilot_cell_trace.tsv) and
[pilot_progress.json](status_revalidation/pilot_progress.json) show A100 repaired
acceptance and A101/A102/A103 **not started**. This is technical acceptance, not
biological accuracy, independent-laboratory validation, four-cell completion or
historical provenance recovery. No original-run Methods gap is closed by replay.

Next decision: review the local repair and packet, then separately authorize a
bounded remaining-three batch and its resource allowance if desired. The current
two-hour authorization covers no additional raw-read execution. Nothing was
pushed, merged, tagged or published.
