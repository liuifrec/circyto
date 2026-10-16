# A02 four-cell technical acceptance

**A02 FOUR-CELL TECHNICAL ACCEPTANCE PASSED — documented mixed execution lineage;
scientific values preserved; no biological-accuracy or historical-recovery claim.**

Executed 2026-10-16 after explicit owner authorization. A101–A103 completed the
normal installed repaired workflow. A100 was not realigned, redetected or
recounted: its unchanged native result and previously accepted status replay
were reused. The four-cell import is joint postprocessing, not a homogeneous
four-cell raw-read workflow. No production code changed in this continuation.

## Identities and actual outcomes

Registration: `17b733c8f21bb128a035e7b100e704099e2bdb0e`.
Analysis base: `05e5171419bffb96a7890c492b7e4ab7531ba18a`.
The commit containing this packet is the analysis review head; its full hash is
reported in the handoff. Original failed A100 evidence remains at
`b9726b2a9ab2cffcacc7bd22a2eae5f969b09234`.

| Cell | Full-library reads | Native candidates | Detector support | Bound new execution / report |
| --- | ---: | ---: | ---: | --- |
| IMR90_A_100 | 3,824,354 | 142 | 474 | Original aligned/completed; repaired derived report completed |
| IMR90_A_101 | 3,935,595 | 95 | 317 | Aligned/completed; remaining-three report completed |
| IMR90_A_102 | 4,002,138 | 112 | 419 | Aligned/completed; remaining-three report completed |
| IMR90_A_103 | 3,954,925 | 123 | 489 | Aligned/completed; remaining-three report completed |

The joint coordinate/strand-compatible union contains **435 candidates**, with
**472 nonzero cell/candidate entries and 1,699 support units**. Per-cell feature
counts were not summed to obtain the union. The separate remaining-three
workflow has 305 unique candidates and 1,225 support units. Support is CIRI3
`#junction_reads`, not a molecule or UMI count. No cell was replaced or omitted.
Accessions, complete base counts, source hashes and both execution identities
are in [pilot_cell_trace.tsv](pilot_cell_trace.tsv).

A100 scientific execution used
`6c280426ebf10abcbcaeab72f87c9d34060d6ee8`; its status replay and the new
A101–A103 execution used `e6843a38469d49f96d44758a4e37893ebbe2734b`.
The reused wheel SHA-256 is
`28c7cbde3c31cf90dd7c247487ba5f410c5b82f6f97d7a2281a0ce52912f7383`.
Its 87 installed Python members match the wheel and repaired source. Metadata
`0.10.0` identifies neither a new public release nor the frozen historical build.

## Input, execution and integration evidence

ENA responses for the three fixed runs matched every registered identity,
layout, URL, byte-size and MD5 field. Complete gzip/FASTQ validation confirmed
read/base totals, unique read identifiers and `/1` mate suffixes without `/2`.
The deposited runs are SINGLE; the generic GEO paired-end prose is retained as
a metadata discrepancy, not used to invent mates. No read/chromosome subset,
trimming, demultiplexing or outcome-based selection was applied.

Existing complete GENCODE v45 primary FASTA/GTF, five BWA index files, external
CIRI3 JAR and pinned split environments were hash-verified and reused. There
was no new reference download, index build or environment solve. Actual emitted
commands retained BWA `mem -t 8 -k 15 -T 15` and CIRI3 `-T 8 -S 0 -Ma 0`,
with Java `-Xmx16g -XX:ActiveProcessorCount=8`. The high-level workflow used
`--parallel 1 --threads 8 --chunk-size 1 --gene-expression-method simple-overlap
--export-h5ad --no-report`, followed by the standalone workflow report command.

Actual BWA/CIRI3 exit records, reference/tool/input hashes, fresh manifests,
completion logs and CIRI3's pre-execution SAM hash bind each new result to its
cell. Biological command intervals do not overlap. The normal combined writer
automatically propagated `aligned` for all three cells into QC, AnnData and
the report; no manual status repair or replay shortcut was used.

The genuine historical RNA slice is **4 × 63,187**, retaining original order
A103,A102,A101,A100. Circ metadata order is A100,A101,A102,A103. Exact biological
IDs, rather than positions, align all four cells in MuData. Every available RNA
value, feature, annotation, layer, embedding and modality metadata field was
preserved; parent metadata remains in the checksum-matched original source.
This is a historical lightweight simple-overlap profile, with historical
reference identity still unverified (`rna_genome_build_verified=false`). The
new remaining-three same-run RNA profile, 3 × 63,241, stays separate.

## Acceptance and retained failures

| Check | Observed result |
| --- | --- |
| Independent native/RNA/object checks | 1,750 PASS; 0 FAIL |
| Independent QC/JSON/TSV/HTML checks | 6,133 PASS; 0 FAIL |
| Offline browser | Both reports opened; 13 tables, 6 cards and 2 SVGs preserved |
| Focused analysis-helper tests | 36 passed; 0 skipped; 0 warnings |
| Resource-control probes | Six passed, including expected nonzero stops; no surviving test descendants |
| Source/manuscript preservation | 1,643 original/repair files, 64 protected files and 40 earlier CMPB files verified unchanged; original branch tips unchanged |

The CSV/exact-integer/CSR oracle does not use circyto parsing or normalization.
It checks every native value, identity, coordinate, strand and available native
annotation, including retained original ancillary/read-ID bytes. Workflow int32
values are compared only by range-checked lossless widening; native import
remains int64. Ordinary AnnData/MuData read/write/read preserves scientific
contents, maps and missingness. Repeated joint import preserves text bytes and
HDF5 scientific contents; HDF5 byte equality is not required.

Report denominators are correctly **3** for the genuine remaining-three workflow
and **4** for the mixed-lineage import. Display checks cover all shown values
and the first-100-row limit, with full tables independently checked. Native
gene/type annotations and separate GTF-derived workflow host annotations retain
their distinct meanings. All four cells have positive support; completed-zero,
failed, unknown and unprocessed states were **not observed in this real pilot**.
These four cells do not establish population prevalence, biological absence or
independent experimental replication.

Original A100's four FAIL assertions remain unchanged: alignment versus report,
alignment versus workflow AnnData, report completion, and completed denominator.
Their repaired results remain in the
[separate status revalidation](../a02_firstcell/status_revalidation.md).
Two continuation checker failures are also retained with original helpers,
logs and outputs: a fixed-decimal-padding assumption (`1.000000` versus `1.00`),
and an incorrect expectation that the importer copies RNA into `evidence/`.
The [amendment](continuation_amendment.md) records both narrow corrections before
rechecking. Exact rational comparison and command-bound original RNA hashes
resolved them without changing scientific assertions or production outputs.
There is no unresolved acceptance failure. Ordinary MuData round trips emitted
eight future-default warnings in each of the three audit attempts; final
workflow-report warnings are empty. The previous **523-test** production-suite
result is reused evidence, not a new final-branch test run.

## Resources and review boundary

The single continuation window was 09:25:19.991–15:25:19.991 JST on October 16.
Managed work finished within 5,297.26 s of its start (about 88.29 minutes).
The final ledger includes the preservation/packet stage and all failed checks.

| Phase | Active elapsed seconds |
| --- | ---: |
| Original shared preparation | 4,344.31 |
| Original first-cell stages | 691.77 |
| Original other controls | 10.24 |
| Prior repair/replay | 176.13 |
| New acquisition and complete FASTQ validation | 1,748.12 |
| New remaining-three normal workflow | 1,987.85 |
| New standalone report / native evidence binding | 1.69 / 34.19 |
| Joint input/RNA preparation / import / repeat import | 14.13 / 2.18 / 1.89 |
| Other new controls, checks and preservation | 124.70 |

New active work totaled **3,914.74 s**, separately from the reconciled prior
**5,222.45 s**; all A02/repair active work totals **9,137.19 s**. These are managed
stage elapsed measurements, including hashing and orchestration, not detector
benchmarks. Human review pauses are not assigned to an algorithm.

Final whole-A02 retained storage was **39.16 GiB**. Maximum sampled summed RSS
was **12.39 GiB**. Cgroup accounting reached its **20 GiB** cap including file
cache; the contemporaneous memory event record shows reclamation and zero
OOM/kill events, and the workflow exited successfully. Swap remained disabled.
Eight-CPU affinity and descendant monitoring, memory/wall termination and the
79 GiB disk-stop watchdog were exercised. No delegated CPU quota or filesystem
quota is claimed. Original deadlines and ledgers remain intact.

[source_manifest.json](source_manifest.json) contains exact argv, environment and
dependency identities, source/output/helper hashes, evidence links, retained
failures and browser hashes. `${CONTINUATION_ROOT}` is the new
`continuation_20261016` directory under the existing authorized `${A02_ROOT}`;
`${ANALYSIS_WORKTREE}` denotes this branch's separate worktree. Detailed command
logs and large artifacts remain there, outside Git. Do not rerun recorded
commands into these occupied output directories. See also
[resource_usage.json](resource_usage.json), [pilot_progress.json](pilot_progress.json),
[native_value_checks.tsv](native_value_checks.tsv) and
[report_value_checks.tsv](report_value_checks.tsv).

No A03–A06, new figures, manuscript interpretation, publication or merge was
performed. Historical detection provenance, biological accuracy and independent
external-laboratory adoption remain separate gaps. The next decision is owner
review of this completed A02 packet before authorizing another work package.
