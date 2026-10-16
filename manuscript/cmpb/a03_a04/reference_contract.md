# Independent reference conversion contract

Analysis-only bespoke Python reference; not an established package or a caller.
No circyto imports, copied production implementation or A02 acceptance-oracle
conversion code. CSV parsing, Python integer arithmetic, SciPy sparse matrices,
ordinary AnnData/MuData APIs and independent HTML generation implement the
following explicit contract. Verification is a separate per-record reader.

Inputs are native CIRI3 3.0.1/unknown single-sample 13-column files with an
explicit cell/sample/file map, the documented long table, or sparse integer
candidate × cell Matrix Market plus feature/cell indices. No coordinate change,
GTF lookup, liftover, threshold or normalization. Native IDs must agree with their
coordinates; long-table IDs may be opaque. Candidate identity is the supplied ID
bound to chromosome/start/end/strand. Conflicting coordinates or supplied
annotations under one ID are errors. Native gene/type are retained as native
annotations; no invented GTF/name attribution. Other ancillary evidence remains
in byte-exact copied source files.

IDs are literal nonempty trimmed strings (`001` and `NA` included), unique on
each axis. Cell metadata defines circ order; candidates follow first appearance
for native/long inputs and supplied matrix feature order otherwise. RNA retains
its own original axis and every field. MuData union order follows circ then
RNA-only cells; standard maps bind identities. No positional joins. No duplicate
cell/candidate pairs, duplicate source mappings or silent aggregation.

Counts are nonnegative exact integers through Python int into int64 CSR; both
entries and total must fit signed int64. Explicit zeros are removed from sparse
nonzero statistics without changing values. Coordinates obey declared closed
one-based or half-open zero-based bounds. Build, coordinate, support and detector
declarations must be internally uniform; native CIRI3 additionally requires its
reviewed closed-coordinate/junction-read declarations. RNA conflicting declared
builds are errors; absent historical declarations remain unverified.

Statuses: completed/failed/unknown/unprocessed/skipped. An absent status column
means unknown, not inferred completion. Completed positive/zero observations
are distinguished; failed positive support remains unchanged. Missing mapped
native files or matrix cell columns cannot certify completed-zero cells. A real
header-only native file or explicit matrix zero column can. Partial modality
membership stays explicit and is not treated as a completed assay.

Per-cell candidate/support totals and per-candidate positive-cell/support totals
include supplied counts. Completed-cohort numerators include only explicitly
completed cells; denominator is their count. Denominator zero gives JSON null,
TSV empty and object NaN, never zero. Fractions are not biological absence or
population prevalence. HTML exposes first 100 rows and declares display rounding
(six decimal places for reference rates); full exact numerators/denominators are
in JSON/TSV. Cosmetic markup and optional provenance field names may differ.

Required endpoint: sparse circ AnnData, optional ordinary MuData retaining RNA,
cell/candidate TSV, JSON summaries, offline HTML, copied non-RNA input evidence
and source SHA-256. RNA input is hashed and preserved in MuData without a
redundant source H5AD copy, matching the evaluated endpoint. Hash all inputs
before parsing and after output; check copied evidence hashes. Never replace an
existing output directory. Fail with a nonzero exit and retained error record.
Post-output oracle checks run outside measured execution for both arms.

Limitations: explicit supported schema only, no workflow orchestration, biology,
new detector compatibility or GTF annotation. Reference configuration and code
are bespoke study work. Output availability is scored separately from incorrect
values. Timing is allowed only after these endpoint requirements pass on the
unaltered four-cell control. Source/configuration hashes are frozen before the
main evaluation; any later correction is dated and failed outputs retained.
