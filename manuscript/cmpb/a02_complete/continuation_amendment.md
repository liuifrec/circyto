# A02 remaining-three continuation amendment

2026-10-16. Owner explicitly authorized the bounded continuation in this session.
Registration `17b733c8f21bb128a035e7b100e704099e2bdb0e` remains unchanged.
Base/review head: `05e5171419bffb96a7890c492b7e4ab7531ba18a`.
Production source: `e6843a38469d49f96d44758a4e37893ebbe2734b`.
Accepted wheel SHA-256:
`28c7cbde3c31cf90dd7c247487ba5f410c5b82f6f97d7a2281a0ce52912f7383`.
All 87 installed Python members were compared with this wheel and production
source before continuation. The repair branch was pushed without force; the new
`analysis/cmpb-a02-complete` worktree starts at the review head.

## Scope and immutable lineage

Run only IMR90_A_101/GSM8558853/SRR30918117,
IMR90_A_102/GSM8558854/SRR30918116, and
IMR90_A_103/GSM8558855/SRR30918115. Keep failures and genuine zeros; no replacement.
A100 retains its original scientific execution under `6c280426...`, original
failed acceptance under `b9726b2...`, and status-only replay under `e6843a3...`.
No A100 alignment, detection or RNA counting is repeated. No production edit,
new detector, reference change, A03–A06, manuscript/figure change or publication.

The persistent continuation root is `${A02_ROOT}/continuation_20261016`, initially
absent. It contains `inputs/manifest_remaining3.tsv`, `workflow_remaining3`,
`inputs/native_four.tsv`, `inputs/cells_four.tsv`,
`inputs/rna/imr90_original_overlap_profile.four.h5ad`, `import_four`, and optionally
`import_four_repeat`. All are new paths. Original workdirs/manifests stay read-only.
The three-cell workflow and four-cell mixed-lineage import are separate cohorts.

## One new bounded window

The six-hour calendar window starts immediately before the first managed
continuation safety stage. Its exact start/deadline are written once to
`records/budget.json` and appended below before launch. Setup, acquisition,
processing, checking and intervening review consume this same window. Previous
deadlines and the reported 5,222.45 s of active work remain separately
recorded and will be reconciled with original stage ledgers.

Reuse the tested runner with only a `continuation` category selecting 20 GiB
memory and this explicit aggregate deadline, instead of the expired first-cell
window. The original A02 deadline is expired; the prior repair session has
finished even if its recorded calendar deadline overlaps this separately
authorized continuation. No earlier window authorizes these new cells.
Eight inherited CPUs (0–7), all-descendant-thread affinity monitoring,
swap disabled, process-tree termination and the existing whole-A02 79 GiB disk
stop remain. The 80 GiB cap includes every old and new retained byte. CPU quota
is not delegated on this host. The existing tiny memory/wall/disk/affinity tests
are repeated, with disk accounting pointed at the whole A02 root. No host
configuration change is permitted. Java retains `-Xmx16g -XX:ActiveProcessorCount=8`.

## Fixed commands and acceptance

Reuse the accepted installed environment, original split BWA/Java prefixes,
external JAR, full GENCODE v45 primary FASTA/GTF and BWA index after hashes match.
Refresh ENA metadata only for the three fixed runs, compare every registered
field, then acquire complete files and run the reviewed gzip/FASTQ/mate validator.
No trimming, selection, subsampling or indexing. New audit wrappers use the
unchanged tool-audit helper with a new ledger root; old wrappers are untouched.

Invoke the installed `circyto workflow full-length-circrna` using the new manifest
and output, `--protocol ramda --skip-demux --detector ciri3 --threads 8 --parallel 1
--chunk-size 1 --export-h5ad --gene-expression-method simple-overlap --no-report`,
with the original registered FASTA/GTF paths. Follow with `circyto report
--workdir workflow_remaining3`. Actual emitted BWA/CIRI3 argv and hash-bound
outcomes are checked against the registered flags. No replay shortcut or manual
status assignment substitutes for the repaired normal workflow.

Joint import uses the unchanged A100 native file plus authentic new native files,
explicit bound cell metadata in A100,A101,A102,A103 order, and the historical RNA
slice in original A103,A102,A101,A100 order. Preserve all 63,187 features and fields.
The historical profile remains lightweight simple-overlap RNA with an unverified
historical reference identity. New same-run RNA stays separate.

Independent CSV/exact-integer/CSR and report checks extend only cohort/path scope.
Freeze helper hashes before checking; no production parser is an oracle. Compare
every value/coordinate/strand/annotation and exact-ID RNA matching, ordinary
object round trips, declared report populations and offline displayed values.
Range-checked lossless widening is allowed only for workflow int32 comparison;
native import stays int64. A new identity/value/schema/status defect stops work
for review. An annotation conflict is retained, not silently harmonized.

Only small audit records and helpers are committed/pushed. Original files,
historical failures and source objects remain unchanged. Exact commands and
private local logs stay in the durable ledger; repository records use logical
path placeholders. This is technical validation, not biological accuracy,
historical provenance recovery or external-laboratory adoption.

## Recorded launch window

Start: 2026-10-16T00:25:19.991361+00:00. Deadline: 2026-10-16T06:25:19.991361+00:00.
Budget record: `${CONTINUATION_ROOT}/records/budget.json`; SHA-256 `98132292dd9a89d078394ad7fc48591d1d600941d0a3ab0aa01cbe5772a95a1c`.

## Audit-format correction before a second independent check

The first independent audit passed 526 native/RNA/round-trip/repeat-output checks,
then stopped at displayed candidate row 20: the inherited one-cell checker
assumed six decimal places, expecting `1.000000`; Pandas rendered the mixed
four-cell fraction column with `1.00`. The source numerator/denominator are 4/4.
This is a checker presentation assumption, not an observed scientific-value,
identity, schema or status-propagation defect. No production output is changed.

Retain `records/acceptance/`, its FAIL row, stage log and launch-time helper
snapshot unchanged. The correction checks the IMPORT fraction using exact
`fractions.Fraction` equality, with no float coercion or tolerance. The four-cell
fractions terminate exactly in decimal; altered values such as 0.26 versus 0.25
must fail. Every count, identity, annotation, denominator and workflow assertion
remains unchanged. Tiny tests cover zero padding and unequal values before the
corrected helper hashes are frozen in a separate v2 record. A new
`records/acceptance_v2/` holds new round trips/checks/browser evidence. Existing
workflow/import directories and the two joint imports are reused read-only.

## RNA-source-location correction before a third independent check

The v2 audit passed 526 native/RNA checks and 4,539 report checks, then raised
`FileNotFoundError` for `import_four/evidence/rna.h5ad`. The inspected importer
deliberately hashes the supplied RNA input and preserves its contents in MuData
without making a redundant H5AD source copy (`import_circ.py`, reviewed source
`e6843a3`, source-copy loop). The checker incorrectly treated every source hash
as a file copied under `evidence/`; this is not a missing supplied RNA input.

Preserve `records/acceptance_v2/` and its exception/log unchanged. For RNA only,
verify the original extracted H5AD against the extraction hash, the actual
executed import `--rna` argument and the report's source hash. Native evidence
copies retain their byte checks. Tiny fixtures reject a changed source or a
different executed RNA path. All matrix/RNA-field/report-value checks remain.
Freeze the v3 helper hashes and run into new `records/acceptance_v3/`; no workflow,
import or prior audit output is modified. Both checker failures remain visible.
