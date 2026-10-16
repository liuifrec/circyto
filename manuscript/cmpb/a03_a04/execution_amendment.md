# A03/A04 execution amendment — 2026-10-16

The supplied owner instruction authorizes this bounded session. Original
registration `17b733c8f21bb128a035e7b100e704099e2bdb0e` is unchanged.
Base: `307758f57a39c009393d03fc0f69e4638b2498c6`; new branch
`analysis/cmpb-a03-a04`. Evaluated production is the repaired
`e6843a38469d49f96d44758a4e37893ebbe2734b`, replacing the old registration's
C2 evaluation target for these new analyses only. Reused installed wheel SHA-256:
`28c7cbde3c31cf90dd7c247487ba5f410c5b82f6f97d7a2281a0ce52912f7383`.
All 87 Python members, module origins, dependency versions and pip check passed
preflight. No package/environment upgrade or production edit is authorized.

## Current status, without changing historical records

| Milestone | Current evidence and boundary |
| --- | --- |
| A02 | [Completed packet](../a02_complete/cohort_acceptance.md), unchanged at the base commit |
| A100 | Original execution `6c280426…`; status replay `e6843a3…`; 142/474 |
| A101/A102/A103 | Normal repaired execution; candidates/support 95/317, 112/419, 123/489 |
| Joint control | 4 × 435 circ, 472 nonzeros, support 1699; authentic native inputs |
| Historical RNA | 4 × 63187; A103,A102,A101,A100 order; same-cell lightweight overlap profile, historical reference unverified |
| Reports | Remaining-three WORKFLOW denominator 3; four-cell IMPORT denominator 4 |
| Original historical detection status | Still unknown; the new pilot does not recover historical records |
| A03/A04 | Authorized here; results recorded separately after execution |
| Named-tool execution / A05 / A06 | Deferred; no Nextflow/plugin installation, CNV work or scaling |

## Budget and storage

Preflight finished and the one session clock started at
2026-10-16T02:39:09.646192Z (11:39:09 JST). A03 preparation/measurement must finish
by 03:39:09.646192Z; the overall deadline is 05:39:09.646192Z. A04 has at most
two hours from its first managed stage, also bounded by that overall deadline.
All setup/checks and human review pauses consume calendar allowance; only actual
process elapsed time is reported as execution time. Old budgets are unchanged.

`${EVAL_ROOT}` is new `a03_a04_20261016` below the existing `${A02_ROOT}`.
Reuse the A02 runner with only scope/cap configuration: CPUs 0–1, 4 GiB cgroup
memory, zero swap, process-tree wall termination, all-thread affinity watchdog,
10 GiB new storage (stop at 9.8), and whole-root 80 GiB (stop at 79). Re-exercise
the existing tiny stop probes. No privileged quota/cache changes. All runs are
sequential. No raw-tool executable is invoked.

No A02 tree copy is made. Necessary inputs are hash-checked against accepted
manifests; a current full path/size/mtime inventory supplements content hashes
for untouched old files. Protected manuscript/recheck content is hashed anew.
New perturbations retain full RNA fields. Each completed case's new input/output
files may be losslessly archived after verification, with member hashes checked
before removing only its redundant uncompressed generated copies. Original
inputs and failed evidence are never removed. Stop at the cap rather than
reduce the registered scientific grid after seeing outcomes.

## A03 matched endpoint

Both arms consume identical native mappings/status metadata and historical RNA.
Reference contract and source snapshots precede timing. Tiny valid/invalid tests
and the unchanged real inputs precede the frozen main evaluation. One unmeasured
warm-up per arm, then fresh processes in balanced pairs C/R, R/C, C/R (C=circyto,
R=reference). Startup, input parsing/validation, source hashes/evidence copies,
sparse AnnData/MuData, QC JSON/TSV and offline HTML are measured in both arms.
Independent post-output checking and archiving are outside both timings.
One BLAS/OpenMP thread, two allowed CPUs; ordinary warm OS cache, no flush.
Record wall/RSS, bytes, exact argv, exit status and all failed attempts.

## A04 fixed realization of registered scenarios

Use all 23 cells and all circ features from the checksum-matched original IMR90
object. Historical statuses remain unknown. Separate DERIVED/SIMULATED fixture
metadata starts with invented completed labels; these do not certify original
assays. RNA remains the source modality (all features/available fields).
Support is preserved detector junction-read support, not UMIs or molecules.

The old YAML key `stochastic_rounding` means **nearest_integer_count**:
`floor(fraction * n_cells + 0.5)`. At n=23, levels 0/5/10/25/50% yield
0/1/2/6/12 selected cells. Record nominal and realized fractions, selected IDs,
and actual moved IDs. Twenty seeds are exactly 2026101500–2026101519.

R01 uses the full seed/level grid with uniform selection: shuffle selected RNA
rows without changing their biological names or values. k=1 may move zero IDs.
R02, R05 and R06 each use the full grid for uniform sampling and weights
`1/(1 + stored RNA overlap-profile total)` without replacement. Validate totals
as finite/nonnegative; label this profile-associated sampling, not verified depth.
R02 alternates selected IDs between RNA removal (even selection positions) and
circ membership removal (odd positions), retaining RNA-only/circ-only members.
R05 changes only selected status, cycling failed/unknown/unprocessed/skipped;
all support remains. R06 marks selected rows failed and makes positive partial
support `max(1, floor(old_support/2))` for positive entries. This is a synthetic
partial-result perturbation, not read downsampling. The oracle reads the supplied
perturbed records for the same completed cohort; full-cohort composition changes
are recorded separately from errors.

Deterministic controls run once with a clean control: R03 duplicate metadata,
RNA IDs and native mapping; R04 literal `001` and `NA`; R07 mixed build,
coordinate and support declarations; R08 duplicated cell/candidate and native
records; R09 explicit sparse zeros; R10 all unknown/no completed denominator;
R11 completed cell without native mapping or explicit matrix source column.
Native cases use tiny explicitly synthetic native files; H5MU-derived values
use the supported matrix/long-table route and are never called recovered native
results. The unsafe positional/zero-fill arm is an illustration only.

Each arm's availability, valid rejection, invalid acceptance, identity/count/
annotation/state errors, same-cohort numerator/denominator and displayed report
agreement are recorded. Undefined prevalence stays null/NA. Output ordering is
checked against each arm's contract; scientific comparisons join explicit IDs.
Source tables include all planned and attempted configurations. A production
defect or source mutation stops execution for owner review; ordinary anticipated
invalid input does not stop the grid. No favorable outcome is required.

The version-pinned documentation comparison in the unchanged study design and
Nikaido audit is reused as documentation evidence. A future named-tool task could
compare preservation of actual library/sample IDs through each tool's supported
count-output stage plus explicit conversion to the same declared modality;
input-stage validation alone is not a full pipeline benchmark. That task needs
separate compatible inputs, endpoints and resources; it is not executed here.
