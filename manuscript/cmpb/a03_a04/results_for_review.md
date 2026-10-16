# A03/A04 execution review — 2026-10-16

**Verdict: COMPLETE WITH LIMITATIONS.** The A03 matched endpoint and all 715 A04
configurations (2145 primary arm attempts) completed. Both evaluated arms preserved
scientific values and correctly rejected the registered invalid inputs. The
additional export/HTML audit and source-preservation checks passed. These are
bounded technical results, not manuscript interpretation or biological accuracy.

## Identity and immutable prior evidence

Branch `analysis/cmpb-a03-a04` starts at accepted A02 evidence
`307758f57a39c009393d03fc0f69e4638b2498c6`. Production is unchanged:
`e6843a38469d49f96d44758a4e37893ebbe2734b`, installed from wheel SHA-256
`28c7cbde3c31cf90dd7c247487ba5f410c5b82f6f97d7a2281a0ce52912f7383`.
All 87 installed Python members matched both that wheel and production source;
module origins, the accepted dependency stack and `pip check` passed preflight.
The unchanged development version string `0.10.0` is not a public release ID.

[A02](../a02_complete/cohort_acceptance.md) remains a completed, separate milestone.
Its real control supplies four native files: A100 142/474, A101 95/317,
A102 112/419, A103 123/489 (positive candidates/junction-read support).
The union is 435 candidates, 472 nonzero cell/candidate entries and 1699 support.
The historical same-cell RNA slice has 4 × 63187 values/feature identities with
its own A103,A102,A101,A100 order. Circ metadata order is A100,A101,A102,A103.
A100 scientific execution was under `6c280426…`, with status replay under
`e6843a3…`; A101–A103 used the normal repaired workflow. The historical RNA is a
lightweight overlap profile with unverified historical reference identity.
It was not replaced by the separate new-run RNA. Native gene/type annotations
retain their native provenance; no GTF computation was added to either arm.

| Milestone | Current status | Unchanged evidence / new record |
| --- | --- | --- |
| A02 real four-cell acceptance | Completed; immutable | [A02 packet](../a02_complete/cohort_acceptance.md), [source manifest](../a02_complete/source_manifest.json) |
| A100 original failed report and repair | Both preserved | [Original failure](../a02_firstcell/first_cell_acceptance.md), [derived status revalidation](../a02_firstcell/status_revalidation.md) |
| A03 reference endpoint and measurements | Completed; endpoint and output checks passed | New `reference_comparison.tsv` and per-run records |
| A04 controlled robustness | Completed; all registered cases accounted for | New `robustness_outcomes.tsv` and complete candidate-level source table |

## A03: executed matched endpoint

Both arms passed the independent unchanged-real-input control, each separate
unmeasured warm-up, and all three fresh-process measurements. The required
endpoint includes native parsing and validation, exact sparse counts, circ
AnnData, matched RNA/circ MuData with its own axes, source hashes and copied
non-RNA evidence, QC TSV/JSON and offline HTML. Every candidate value,
coordinate, strand and available native annotation, RNA field and identity/map,
report denominator and displayed table value was checked after measured execution.
Cosmetic markup and nonessential provenance fields differ.

| Arm | Three wall times (s) | Median wall (s) | Three peak RSS values (KiB) | Median peak RSS (KiB) |
| --- | --- | ---: | --- | ---: |
| circyto | 1.243331, 1.199263, 1.190269 | 1.199263 | 188036, 187572, 187336 | 187572 |
| Independent reference | 0.653512, 0.632020, 0.669997 | 0.653512 | 156668, 156000, 155860 | 156000 |

The reference was faster and used less peak RSS in this tiny-input comparison.
This is a bespoke study converter, not an established competing package.
It does not establish scaling, human usability or superiority over named tools.
Runs were balanced/interleaved C1/R1, R2/C2, C3/R3 with startup included, one
BLAS/OpenMP thread, two allowed CPUs and warm OS cache. No privileged cache flush
was used. Identical input logical bytes were 10862525; output footprint was
27762688 for circyto and 15622144 for the reference. Independent acceptance and
retention work were outside both measured processes.

Exact per-run values, argv, dependencies, source/configuration hashes and output
hashes are retained in the source tables/manifests and local execution records.
The operation inventory describes actual supplied metadata and transformations;
it does not score code length or invent human effort. Repository manifests use
`${EVAL_ROOT}`, `${A02_ROOT}`, `${ANALYSIS_WORKTREE}` and
`${ORIGINAL_OBJECT_ROOT}` aliases to avoid publishing workstation paths. The
retained local records contain the exact emitted argv and original paths;
source hashes bind those records and the archived per-case evidence.

## A04: executed outcomes

| Arm | Valid inputs accepted | Valid cases passing every value/state check | Valid cases failing checks | Invalid inputs correctly rejected | False rejection of valid input |
| --- | ---: | ---: | ---: | ---: | ---: |
| circyto | 705 | 705 | 0 | 10 | 0 |
| Independent reference | 705 | 705 | 0 | 10 | 0 |
| Unsafe illustration | 705 | 172 | 533 | 10 | 0 |

No required field was unavailable in these outputs. Both evaluated arms had zero
count, annotation, identity, state, source, report and contract errors. The unsafe
arm's invalid rejections use the shared strict validator and are not independent
evidence of heuristic validation ability. Its failures were: native clean control
1, R01 52, R02 159, R05 160, R06 160 and R10 1. Its 172 valid passes are retained.
Error totals count assertion/field violations; they are not counts of distinct
misassigned cells or error probabilities.

| Scenario | Configurations per arm | Expected input class |
| --- | ---: | --- |
| Shared matrix/native controls | 2 | Valid; native zero case is synthetic |
| R01: RNA order | 100 | Valid |
| R02: modality membership | 200 | Valid |
| R03: duplicate cell identities | 3 | Invalid |
| R04: literal IDs | 1 | Valid |
| R05: processing states | 200 | Valid |
| R06: failed rows with partial positive support | 200 | Valid |
| R07: conflicting declarations | 3 | Invalid |
| R08: duplicate records | 2 | Invalid |
| R09: explicit sparse zero | 1 | Valid |
| R10: no completed denominator | 1 | Valid |
| R11: missing result/column for declared completion | 2 | Invalid |

The complete candidate-level table contains **5,159,619 rows** (1,719,873 per arm).
Both evaluated arms had zero exact same-cohort prevalence errors and retained all
2443 undefined candidate results in the no-denominator control. The unsafe arm had
781,440 nonzero same-cohort candidate/configuration errors and replaced all 2443
undefined results with defined values. These are table rows, not independent
biological observations. The repository extract contains the first candidate of
each valid case/arm (2115 rows); it is a deterministic review extract. The complete
lossless table remains outside Git at
`${EVAL_ROOT}/a04_combined/completed_cohort_prevalence.tsv.gz` (71,757,173 bytes;
SHA-256 `8163178301f702fb4217fc0e48d4390ed82aa39244a8d1efd2cc681ad9152b8c`).
All three source gzip streams passed integrity checks; 4886 earlier-attempt rows
remain in the original stream and are excluded from the primary combined table.

The additional independent output audit checked 2127 retained accepted-output
bundles: 711 circyto and 711 reference bundles all passed; the unsafe illustration
had 441 failures and 264 passes in this additional subset of checks. This includes
ten A03 bundles and two retained administrative-attempt bundles. It checks the
separate Matrix Market export, JSON declarations/states/overlap, HTML summary text
and absence of remote rendering resources. The primary oracle separately checks
all source values, RNA fields/maps and first-100 HTML table values. Browser-based
interaction/usability was not tested.

## A04: provenance and interpretation

The fixed plan has 715 configurations and 2145 arm attempts. The full 23-cell
IMR90 source object matched
`bb2e12f7c3b36f9fa72d98cd71e8bea905a67f50e22af1d6b713550ee92b60c8`.
Its 23 × 2443 circ and 23 × 63187 RNA fields are retained in separate derived
fixtures. Historical completion stays unknown. Invented completed labels apply
only to the explicitly DERIVED/SIMULATED experimental fixtures. Native edge
fixtures are separately synthetic. A02 supplies no observed real failed,
unknown or completed-zero example.

The [execution amendment](execution_amendment.md) records nearest-integer counts
0/1/2/6/12 at the five nominal levels, realized fractions, the twenty exact seeds,
uniform and overlap-profile-associated sampling, and each change to membership,
status or partial support. Selecting one identity in a permutation may move none.
The source overlap-profile total is not verified sequencing depth. Complete
same-cohort numerators/denominators are distinct from full-cohort composition.
Undefined prevalence is null/NA, not zero; fractions are not population prevalence.

The unsafe arm is a labelled positional/status-zero-fill illustration. It shares
the competent reference's strict validator and writer; its invalid-input
rejections therefore are inherited checks, not evidence of a heuristic's own
ability. Its deliberate identity/status mistakes are not attributed to named
tools or typical analysts. Error counts are assertion/field violations, not
cell-error probabilities. Output availability is separate from wrong values.

Every planned case/arm is present in `robustness_outcomes.tsv`; none is pending.
The two owner-authorized administrative checkpoints (exit 130) applied the storage
and time amendments. Complete archived cases were reused; partial cases were
repeated in separate directories. Four earlier arm records remain separately in
`administrative_prior_attempts.tsv`: two completed checks, one executed arm whose
check was interrupted, and one started arm without a completion record. None was
silently discarded or counted twice. Original source streams, partial directories,
logs and verified lossless archives are retained.

## Final validation, resources and retained limitations

The final focused helper run passed **28 tests plus 4 subtests**, with zero skips,
failures or errors and 43 warnings. Forty-two warnings describe the future MuData
obs/var update default; one is the intentional duplicate-ID fixture warning.
The final record is `records/a04-final-review-tests.xml` with its stage/log linked
from `resource_usage.json`. Syntax checks and `git diff --check` accompany the
reviewed commit. Empty final TSV fields in repository review copies are quoted;
`review_formatting.json` records before/after hashes and exact parsed-field equality.
The previous 523-test production suite remains historical
validation; it was not rerun for this analysis-only branch.

A [bookkeeping correction](bookkeeping_correction.md) cleared a stale `file absent`
initializer beside already successful gzip checks. The original generated manifest,
executed helper snapshot and before/after hashes remain retained. The focused
regression and final source suite passed. No evaluated-reference or production
defect was observed, and no scientific outcome changed in that correction.

Preservation passed: all 12 necessary input hashes, 112 protected-file hashes,
both original H5MU hashes and all 87 installed Python-member hashes matched.
All original branch tips matched. The prior A02 tree had 81,588 unchanged
path/size/mtime/symlink entries; this is an inventory comparison, not a claim of
freshly hashing every byte of that tree. A02 remains an immutable completed
milestone, including its original failures and mixed execution lineage.

The authorized amendments raised new storage to 16 GiB and reallocated 20 minutes
of unused A03 time, making A04's deadline 05:11:18 UTC. The overall deadline remained
05:39:09 UTC. All managed evaluation/test stages finished before the A04 deadline.
Active managed-stage time was 7300.362 seconds, separately recorded from human
review and A03 per-process timings. New persistent storage was 14.181 GiB; the
shared A02 tree was 53.345 GiB. Two CPUs, 4 GiB cgroup memory and disabled swap were
enforced. The sampled process-tree RSS peak was 665.03 MiB; the cgroup peak reached
4 GiB including cache. All six tiny limit probes passed and left no test descendants.
CPU enforcement used inherited affinity and the all-thread watchdog; the host did
not expose a CPU quota controller. Disk enforcement was a watchdog, not a quota.
Original A02/repair resource ledgers were not reset or included as new allocation.

`reference_comparison.tsv`, `operation_inventory.tsv`, `robustness_outcomes.tsv`,
`completed_cohort_prevalence.tsv`, `additional_output_checks.tsv` and
`figure_source_summary.tsv` provide factual source tables for the proposed
Figures 3/4. `input_output_manifest.json`, `environment.json` and
`resource_usage.json` bind inputs, source snapshots, output archives, exact local
records and resource accounting. No polished figure or manuscript interpretation
was generated. No production edit, raw-data task, named-tool installation,
A05/A06 analysis, release or manuscript rewrite occurred.

The version-pinned documentation comparison remains documentation-based.
The deferred nf-core input-stage harness is not a pipeline benchmark. A future
named-tool task could start each documented single-end route with the fixed four
IMR90 libraries, retain their biological IDs through each tool's own count
outputs, and audit an explicit conversion to the same historical RNA modality.
The existing pins are nf-core/circrna `23fdf1d22698a34e06efa99bf5af68d1bb963b78`
and CirComPara2 v0.1.2.1 `303c3621579c1a412493e40772c2ff2660148942` (inspected
2026-10-15 in the unchanged Nikaido documentation audit at `42459975`). Settings,
supported caller/quantification semantics and the conversion endpoint would
need separate review before execution. Count differences would not measure
accuracy, and complete-pipeline time would not be compared with import-only time.
No biological accuracy, historical provenance recovery, independent
laboratory adoption or methodological novelty is established here.
