# Detectors

This page explains which detectors circyto can orchestrate, what each one needs, and when to choose which.

> Planned: `circyto detectors` will print a live, authoritative list. Until then, this doc is the reference.

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
circyto run-batch \
  --detector find-circ3 \
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
- workflows where you want CIRI-full’s reconstruction/quantification behavior

**External requirements**
- `bwa`
- `java` (JRE)
- the CIRI-full JAR placed in `tools/` per this repo’s `tools/` layout

**Typical run**
```bash
circyto run-batch \
  --detector ciri-full \
  --manifest manifest.tsv \
  --outdir work/ciri_full_run \
  --ref-fa ref/genome.fa \
  --gtf ref/genes.gtf \
  --threads 8 \
  --parallel 4
```

Notes:
- CIRI-full’s own documentation indicates it is Java-based and requires `bwa` for generating SAM inputs. It bundles CIRI2 and CIRI-AS as part of the CIRI-full software package.
- If you do not have a GTF, you can still run on a reference FASTA, but downstream annotation will be limited.

---

## STAR-based detectors (future)

STAR-based detectors are planned/under consideration and will be treated as **optional integrations** (because STAR is a large external dependency).

Examples discussed for future integration:
- DCC
- circhunter
- CIRI3

---

## Planned UX: `circyto detectors`

Proposed output (example):

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

## Planned UX: `circyto doctor`

Proposed checks:
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
