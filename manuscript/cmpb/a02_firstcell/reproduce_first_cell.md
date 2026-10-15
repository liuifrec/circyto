# A02 first-cell execution record and reproduction

This is a **NEW technical validation**, restricted to IMR90_A_100 /
GSM8558852 / SRR30918126. It does not recover the original IMR90 run.
The registered cohort still contains A101–A103, which require separate review.
See `first_cell_acceptance.md` for the actual outcome; commands alone do not
establish completion or acceptance.

The registration base is `17b733c8f21bb128a035e7b100e704099e2bdb0e`;
the evaluated circyto source is `6c280426ebf10abcbcaeab72f87c9d34060d6ee8`.
The local reviewed wheel SHA-256 is
`eca5f45f35b351ec44df4e5b0f73eb4e5e0564e2691631f8eb417b0e8a7ca03b`.
Its development metadata says 0.10.0; this is neither a new public release nor
the frozen historical main implementation. The installed Python members were
compared with the wheel. No editable install or checkout import is used.

The durable evaluation root is
`/home/liuyuchen/data/circyto_cmpb_a02_20261015` (`EVAL_ROOT` below).
It retains inputs, reference downloads and indices, both isolated environments,
external JAR, logs, source snapshots, partial/error evidence and outputs.
`commands.json` records the actual argv, selected environment, times, exits,
log paths and hashes for every managed stage and external tool invocation.
`input_manifest.json` records provider evidence, input hashes, cell mappings,
index hashes and the exact one-cell manifests. `environment_manifest.json`
records resolved Python distributions, module origins, external executables,
configuration and references to the full conda explicit lists and pip freeze.
`output_artifact_manifest.json` identifies retained outputs by bytes and hash.

No stage may be replayed into this root. The runner rejects a recorded stage
name and holds a single-stage lock. A separately approved reproduction needs a
fresh persistent root and corresponding paths, with the same immutable inputs
and source hashes; it is another new validation. A future remaining-three run
also needs its own manifest and workdir. Do not resume this one-cell workdir
with a four-cell manifest.

## Execution mechanism

The exact individual request specifications are retained in `EVAL_ROOT/specs`.
For a reviewed, fresh stage the launcher syntax is:

```bash
python3 scripts/cmpb_a02/run_stage.py launch \
  --root "$EVAL_ROOT" --spec /absolute/path/to/reviewed-stage.json
```

This requires a functioning user systemd manager and delegated cgroup memory
control. It does not require privileged host changes. Tiny local controls must
pass before acquisition/alignment: inherited eight-CPU affinity, descendant
termination on wall timeout, memory failure, disk growth and affinity widening,
plus successful and failing exit propagation. This host did not delegate the
CPU quota controller. The enforced CPU restriction is affinity 0–7 plus the
all-thread watchdog; an accepted `CPUQuota` property is not evidence of a quota.
Memory/swap and wall limits apply to the service process tree. Disk monitoring
stops at 79 GiB with headroom below 80 GiB; it is not a filesystem quota.

The allowance is shared: 16 hours including preparation, 80 GiB total retained
storage, eight CPUs, 24 GiB memory; first-cell stages use the stricter shared
four-hour deadline and 20 GiB cgroup cap. Java receives `-Xmx16g` and
`-XX:ActiveProcessorCount=8`. The first-cell budget starts with original RNA
extraction, so this envelope is conservative. Do not reset the allowance when
reviewing remaining cells. Timing definitions and sampling limits are in
`resource_usage.json`.

## Registered operations

The workflow invocation uses the installed executable under `EVAL_ROOT/env`:

```bash
"$EVAL_ROOT/env/bin/circyto" workflow full-length-circrna \
  --manifest "$EVAL_ROOT/inputs/manifest_firstcell.tsv" \
  --outdir "$EVAL_ROOT/workflow_firstcell" \
  --protocol ramda --skip-demux \
  --genome-fasta "$EVAL_ROOT/ref/GRCh38.primary_assembly.genome.fa" \
  --gtf "$EVAL_ROOT/ref/gencode.v45.primary_assembly.annotation.gtf" \
  --detector ciri3 --threads 8 --parallel 1 --chunk-size 1 \
  --export-h5ad --gene-expression-method simple-overlap --no-report
```

The environment recorded in each specification selects BWA 0.7.17 from
`env-bwa`, Java 17 from `env`, and the external hash-matched CIRI3 JAR via audit
wrappers. The wrappers forward argv and standard streams unchanged and record
input hashes, file-output hashes where available, and exits. Expected BWA flags
are `mem -t 8 -k 15 -T 15`;
CIRI3 receives `-T 8 -S 0 -Ma 0` with the full FASTA/GTF. Actual emitted argv,
rather than this abbreviated description, are authoritative in `commands.json`.
Wrapper elapsed times include post-execution hashing and exclude initial input
hashing; they are not isolated algorithm benchmarks.

After a successful workflow, a separate binding stage requires zero tool exits,
unchanged hashed inputs, the fresh workflow alignment manifest, the recorded
CIRI3 pre-execution SAM hash (still matching the file), and the native output
hash before declaring the **new** detector status completed. BWA stdout passes
through Python's byte-copy pipe; its wrapper record correctly retains a null
output hash. The added binding method and limitation are explicit in
`input_manifest.json`, rather than rewriting that original audit record.
It does not assign completion to any historical cell or to an unstarted cell.
The real workflow report is generated only from `workflow_firstcell`:

```bash
"$EVAL_ROOT/env/bin/circyto" report --workdir "$EVAL_ROOT/workflow_firstcell"
```

The original H5MU must match SHA-256
`bb2e12f7c3b36f9fa72d98cd71e8bea905a67f50e22af1d6b713550ee92b60c8`.
The independent helper extracts A100 by its exact biological identity, preserving
all available RNA fields. Original parent metadata is retained separately with
its evidence class. The historical lightweight simple-overlap profile retains
an unverified reference identity. The newly generated overlap profile remains
separate; neither profile is an independently generated transcript estimate.

The native import is a distinct output bundle, with its own IMPORT HTML:

```bash
"$EVAL_ROOT/env/bin/circyto" import-circ \
  --ciri3-manifest "$EVAL_ROOT/inputs/native_firstcell.tsv" \
  --cells "$EVAL_ROOT/inputs/cells_firstcell.tsv" \
  --rna "$EVAL_ROOT/inputs/rna/imr90_original_overlap_profile.firstcell.h5ad" \
  --detector CIRI3 --genome-build GRCh38 \
  --coordinate-convention 1-based-closed \
  --support-definition CIRI3_junction_reads \
  --outdir "$EVAL_ROOT/import_firstcell"
```

Never run the workflow-report command on this import bundle. The independent
CSV/exact-integer/CSR checker imports no circyto parsing or reporting functions.
It compares every native candidate/value/coordinate/strand/annotation, retained
source bytes and cell axes, all available RNA fields, ordinary AnnData/MuData
read/write/read behavior, and JSON/TSV/HTML values. A separate headless browser
opens each local report with DNS disabled and compares rendered DOM tables,
cards and SVG values. Screenshots remain outside Git.

The native importer preserves original native gene/type fields and source read
IDs. The workflow collector has its own GTF-derived host annotation and retains
matrix cell columns under its registered cell-total filter. Omitted zero matrix
columns and expanded selected-cell AnnData axes must be explained, not changed
to manufacture agreement. Equality with historical circRNA yields is not tested.

The executed helper snapshot and its hashes are recorded in
`executed_helper_manifest.json`; individual stage records also capture source
hashes. Durable source snapshots remain under `EVAL_ROOT/helpers`. The later
repository cleanup of unused imports in `tool_audit.py` and `probe_controls.py`
does not change the frozen executed files. The package implementation is unchanged.
The completion-binding helper has a separately tested pipe-output amendment;
its actual later source hash is in the `bind-firstcell-native` stage record.
The workflow-only value checker has a lossless integer-width amendment; its
initial failure and outputs are retained. After rendering checks, the separate
`workflow-status-acceptance` stage compares recorded alignment completion with
workflow QC/AnnData/report semantics and deliberately exits nonzero on the
observed mismatch. That failure prevents overall acceptance despite preserved
counts and successful browser rendering. Production outputs are unchanged.

## Review boundary

The public inputs are the complete deposited ENA single-end file and the
[GENCODE v45 primary assembly and matching annotation](https://www.gencodegenes.org/human/release_45.html).
Provider checksums, local compressed/uncompressed hashes and the layout conflict
with generic GEO prose are recorded rather than inferred from filenames.
No raw inputs, alignments, external JAR, indices or large scientific objects
belong in Git or distributions. No remaining cell, A03–A06 analysis, new release,
merge, tag or publication is authorized by these reproduction notes.
