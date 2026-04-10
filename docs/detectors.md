# Detectors

This page explains which detectors circyto can orchestrate, what each one needs, and when to choose which.

Use `circyto detectors` for the live readiness report and this page for behavior notes and workflow guidance.

---

## Execution tracks

circyto now supports two execution tracks:

- `per-cell-fastq`: the original manifest-driven FASTQ -> align -> detector flow
- `alignment-first`: a shared preprocessing flow where alignments are prepared once, cached, and reused by alignment-capable detectors

The alignment-first flow is designed for large demultiplexed datasets where rerunning alignments per detector or per rerun is not computationally practical.

Typical alignment-first flow:

```bash
circyto prepare-alignment-cache \
  --manifest manifest.tsv \
  --aligner reuse-existing \
  --detector ciri3 \
  --outdir work/alignment_cache

circyto run-detector-from-alignments \
  --detector ciri3 \
  --manifest work/alignment_cache/alignment_manifest.tsv \
  --outdir work/ciri3_from_alignments \
  --ref-fa ref/genome.fa
```

The normalized detector TSV schema remains:

```text
circ_id    chr    start    end    strand    support
```

---

## Detector summary

### `find-circ3`

**Best for**
- quick smoke testing
- junction/anchor-based circRNA detection workflows

**External requirements**
- `bowtie2`
- `samtools`
- `find-circ3` CLI installed and available on `PATH`

**Typical run**
```bash
circyto run-detector find-circ3 \
  --manifest manifest.tsv \
  --outdir work/find_circ3_run \
  --ref-fa ref/genome.fa \
  --threads 8 \
  --parallel 4
```

Notes:
- A GTF is not required for the detector call itself.
- You may still want a GTF later for host-gene annotation (`annotate-host-genes`).

---

### `ciri-full`

**Best for**
- full-length short-read protocols (e.g., Smart-seq2-like data)
- workflows where you want CIRI-full-style calling through the unified detector interface

**External requirements**
- `bwa`
- `java` (JRE)
- the CIRI-full JAR placed in `tools/` per this repo’s `tools/` layout

**Typical run**
```bash
circyto run-detector ciri-full \
  --manifest manifest.tsv \
  --outdir work/ciri_full_run \
  --ref-fa ref/genome.fa \
  --gtf ref/genes.gtf \
  --threads 8 \
  --parallel 1
```

Notes:
- `circyto` reports `ciri-full` as READY when the required external tools and bundled assets are available. That status does **not** imply one identical internal execution path for every read layout.
- **Paired-end input (`r1` + `r2`)** uses the upstream bundled **CIRI-full Java Pipeline**.
- **Single-end input (`r1` only)** does **not** use the upstream CIRI-full Pipeline directly, because upstream Pipeline mode is paired-end only. `circyto` falls back to the bundled **CIRI2-based detection path** and normalizes the result to the same TSV schema.
- The normalized TSV schema is the same for both layouts:

```text
circ_id    chr    start    end    strand    support
```

- This preserves compatibility with downstream collection and multi-detector workflows, but **single-end `ciri-full` runs are not equivalent to upstream CIRI-full full-length reconstruction**.
- If you do not have a GTF, the detector may still run with limited annotation context, but `--gtf` is the recommended public workflow.

### `ciri2`

**Best for**
- direct CIRI2-based circRNA detection
- single-end short-read datasets where you want the explicit CIRI2 path rather than the `ciri-full` compatibility wrapper

**External requirements**
- `bwa`
- `perl`
- bundled `tools/CIRI-full_v2.0/bin/CIRI_v2.0.6/CIRI2.pl`

**Typical run**
```bash
circyto run-detector ciri2 \
  --manifest manifest.tsv \
  --outdir work/ciri2_run \
  --ref-fa ref/genome.fa \
  --gtf ref/genes.gtf \
  --threads 8 \
  --parallel 4
```

Notes:
- `ciri2` is the direct bundled CIRI2 integration.
- For short single-end reads, the wrapper applies more permissive BWA/CIRI2 settings than the default paired-end path.

---

### `ciri3`

**Best for**
- alignment-first workflows on large demultiplexed datasets
- detector runs where BAM/SAM reuse is more important than repeated FASTQ alignment
- multisample-alignment backends with reusable cached alignments

**External requirements**
- `java` (mandatory)
- a local CIRI3 jar under `tools/CIRI3/` or configured via `CIRCYTO_CIRI3_HOME` / `CIRCYTO_CIRI3_JAR`
- `samtools`
- `bwa` for BWA mode
- `STAR` only for STAR mode

**Typical run**
```bash
circyto prepare-alignment-cache \
  --manifest manifest.tsv \
  --aligner bwa-mem \
  --detector ciri3 \
  --ref-fa ref/genome.fa \
  --outdir work/alignment_cache

circyto run-detector-from-alignments \
  --detector ciri3 \
  --manifest work/alignment_cache/alignment_manifest.tsv \
  --outdir work/ciri3_run \
  --ref-fa ref/genome.fa
```

Notes:
- `ciri3` is an alignment-first detector with a real backend.
- Direct execution is `java -jar` based and is reported as `READY` when circyto finds both a usable CIRI3 jar and Java.
- `NOT READY` means circyto could not find a usable CIRI3 jar or execution path.
- `PARTIAL` means local CIRI3 assets were found but the runtime contract is incomplete.
- Supported runtime environment variables include `CIRCYTO_CIRI3_HOME`, `CIRCYTO_CIRI3_JAR`, `CIRCYTO_CIRI3_JAVA`, and `CIRCYTO_CIRI3_EXTRA_ARGS`.
- `CIRCYTO_CIRI3_CMD_TEMPLATE` and `--command-template` remain available for explicit template execution.
- BWA mode requires unsorted SAM input, is alignment-sensitive, and is the validated local production path.
- Validated local BWA settings are `-S 0`, `-Ma 0`, and `bwa mem -k 15 -T 15`.
- STAR mode is supported in code for alignment-first workflows, also alignment-sensitive, and uses `-Ma 1`.
- STAR mode is not yet fully validated end-to-end in this release.
- The backend normalizes tabular CIRI3 output to the stable circyto TSV schema.

Verification example:

```bash
export CIRCYTO_CIRI3_HOME=/path/to/CIRI3
circyto doctor
circyto detectors
```

---

## Other detectors (future)

Additional STAR-based or alignment-native detectors are still planned as optional integrations.

Examples:
- DCC
- circhunter

---

## `circyto detectors`

Example output shape:

```text
$ circyto detectors

NAME        TYPE        NEEDS                     NOTES
find-circ3   CLI tool    bowtie2, samtools         good first-run smoke test
ciri-full    JAR tool    bwa, java                 requires tools/CIRI-full*.jar
```

Optional flags under consideration:
- `--json` for tooling
- `--verbose` to show expected output layout and required inputs

---

## `circyto doctor`

Current checks:
- detect executables on PATH: `bowtie2`, `samtools`, `bwa`, `java`, `STAR` (optional)
- detect the presence of required detector assets under `tools/`
- print one actionable line per missing dependency

Example:

```text
$ circyto doctor

[OK] python: 3.11.7
[OK] bowtie2: /usr/bin/bowtie2
[OK] samtools: /usr/bin/samtools
[MISSING] bwa: not found on PATH (needed for ciri-full)
[OK] java: 17.0.10
[WARN] STAR: not found (only needed for STAR-based detectors)
[WARN] tools/: CIRI-full jar not found (needed for ciri-full)
```
