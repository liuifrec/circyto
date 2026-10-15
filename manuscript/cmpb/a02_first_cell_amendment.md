# A02 amendment: first-cell gate

2026-10-15; parent registration `17b733c8f21bb128a035e7b100e704099e2bdb0e`.
The six registered files remain byte-for-byte unchanged. This amendment records
the subsequent first-cell instruction before any acquisition or detector job.
Evaluated implementation remains `6c280426ebf10abcbcaeab72f87c9d34060d6ee8`.

- Route A: the named object-bound historical IMR90 workdir is not mounted; no
  authorized remote endpoint or mounted replacement was supplied. No SSH, new
  archival sweep or historical-run reconstruction is performed. Original
  Methods/completion gaps remain open.
- Route B: NEW validation of **IMR90_A_100 / GSM8558852 / SRR30918126 only**.
  A101/A102/A103 remain in the registered cohort, **not started**. No substitution
  for low/zero signal. Full library and full registered reference are required.
- `EVALUATION_ROOT` is the persistent local directory
  `/home/liuyuchen/data/circyto_cmpb_a02_20261015`. Large inputs, external JAR,
  isolated environment, caches, temporary files and all outputs stay there.
  Only scripts and small audit records are candidates for Git.
- Replace the proposed four-cell execution manifest with
  `inputs/manifest_firstcell.tsv`, output `workflow_firstcell`; later native
  import uses `inputs/native_firstcell.tsv`, `inputs/cells_firstcell.tsv`,
  `inputs/rna/imr90_original_overlap_profile.firstcell.h5ad`, and
  `import_firstcell`. The four-cell selection remains separate. No merged
  workflow summary or changed-cohort resume is allowed.
- All other P01–P05 flags, BWA `-t 8 -k 15 -T 15`, CIRI3 `-T 8 -S 0 -Ma 0`,
  reference URLs, scientific constraints, simple-overlap and `--no-report` stay
  registered. Workflow HTML is generated separately from genuine workflow QC;
  IMPORT HTML belongs to its distinct native import bundle.
- The historical RNA slice retains its exact cell name and all available RNA
  fields. It remains a genuine same-cell **lightweight simple-overlap profile**,
  with incompletely identified historical reference. The new overlap profile
  remains separate. No RNA extraction or native acceptance is yet claimed here.

The original caps remain shared across the pilot: eight CPUs, one biological
cell at a time, 24 GiB process-tree memory, Java heap at most 16 GiB, 80 GiB total
evaluation storage and 16 hours including preparation. First-cell work has a
four-hour budget and a stricter 20 GiB cgroup memory ceiling, including page
cache as well as processes; this can stop earlier than a 20 GiB RSS trigger.
No new budget is granted to each stage. The current runner also imposes a
conservative absolute deadline from evaluation initialization; later review
must reconcile actual consumed time without silently resetting the allowance.

Resource preflight is a gate, not a scientific experiment. User-service cgroups
apply memory/swap limits and process-tree timeout/termination. The CPU quota
controller is **not delegated** on this host: an accepted `CPUQuota` property
alone is not evidence of enforcement. The candidate alternative is inherited
eight-CPU affinity plus a monitor of every descendant thread; it must pass a
deliberate affinity-widening stop test. A one-second disk monitor stops at 79 GiB,
leaving 1 GiB before the 80 GiB ceiling; this is an automatic stop, not a filesystem
quota. Tiny tests must demonstrate descendant termination for memory, wall,
disk and affinity violations before any long-running work. Any failed control
test blocks execution. No privileged host configuration is changed.

Environment amendment, recorded before installation succeeds: a single prefix
containing Python 3.10.22, BWA 0.7.17 and Java 17 failed dependency resolution.
The available BWA 0.7.17 builds require zlib below 1.3, whereas Python 3.10.22
requires libzlib at least 1.3.2. Retain the solver failure; use `env` for Python
and Java 17 and `env-bwa` for the unchanged BWA 0.7.17. Both are new prefixes
under the shared storage budget, using conda-forge then bioconda with strict
priority. The system Java 8 is not the selected runtime. No scientific software
version or detector flag is substituted to resolve this packaging conflict.
`JAVA_TOOL_OPTIONS` records `-Xmx16g -XX:ActiveProcessorCount=8`; the latter
keeps Java runtime thread-pool sizing within the allocated CPU mask and does
not change the registered CIRI3 scientific flags.

Input-audit clarification before alignment: the acquired ENA headers retain an
original flow-cell identifier ending in `/1` in their second field. The initial
validator counted unique SRR identifiers and verified all reads, bases, CRC and
hashes, but did not count this particular mate-suffix location. Retain its record
and add a complete second validation that recognizes this suffix. A bare numeric
SRA spot ordinal is not a mate flag; any explicit `/2` blocks the registered
single-end route. This is an audit-helper correction, not a read transformation.

Only the analysis-scoped A02 audit helpers are added. This amendment authorizes
no A03–A06 work, revised biological claim, production change, new release, merge
or publication. Later-design notes about nf-core, rounding and CNV remain deferred.
The outcome and performed-versus-pending operations will be recorded in
`a02_firstcell/first_cell_acceptance.md`; this amendment itself is not acceptance.

Execution-audit amendment before completion binding: the installed workflow sends
BWA stdout through a Python pipe, copies it to `alignment.sam.partial`, checks
the BWA exit, and renames it to `alignment.sam`. The wrapper's stdout-file probe
therefore recorded no output hash despite BWA exit 0. Preserve that original
record. The additional binding check may use the fresh one-cell alignment
manifest (`reused_alignment=false`), exact SAM path, unchanged input hashes,
BWA exit 0, and the CIRI3 input hash recorded before detector execution. It must
verify the SAM still matches that hash and preserve the alignment manifest/log
hashes. Missing/conflicting evidence leaves completion unbound. This correction
is tested with tiny synthetic manifests before use; it neither reruns BWA nor
changes the detector/production workflow or any scientific threshold.

Independent-check amendment after the first real check: all native-import/RNA
and round-trip checks passed, but the workflow comparison rejected its `int32`
CSR dtype against the native reader's `int64`. The reviewed legacy workflow
writer explicitly uses `int32`; both actual matrices have the same 1 × 142 axes.
Retain the failed stage/log and snapshot its audit outputs before a new check.
The workflow-only comparison now requires integer sparse data, checks range,
widens losslessly to int64 in memory, and requires every value/index to match.
Floats and changed counts still fail. Native-import and round-trip dtype checks
stay strict. No output or scientific setting is changed, and this one-cell check
does not establish safety for counts outside the legacy int32 range.
