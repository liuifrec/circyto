Multi-Detector Integration (v0.6.0+)

circyto v0.6.0 introduces the detector engine API and the first working implementation of multi-detector orchestration.

This document explains how detectors work, how to run them individually, and how to run multiple detectors in parallel on the same set of single-cell FASTQ files.

1. Detector Engine Overview

A detector is defined by a Python class that implements:

class DetectorBase:
    name: str
    input_type: str                # e.g., "fastq"
    supports_paired_end: bool

    def is_available(self) -> bool: ...
    def version(self) -> Optional[str]: ...
    def run(self, inputs: DetectorRunInputs) -> DetectorResult: ...


DetectorRunInputs provides:

cell_id

r1 / r2

ref_fa

gtf

threads

outdir

Each detector must write a normalized TSV with the columns:

circ_id   chr   start   end   strand   support


This is the format consumed by collect, convert, and downstream multimodal export.

Implemented detectors (v0.6.0)
Detector	Status	Notes
ciri-full	✓ stable	full CIRI-full Pipeline via ciri_full_adapter.sh
ciri2	✓ experimental	wrapper for CIRI2.pl; strict filtering often gives sparse output
2. Running a Single Detector

Use:

circyto run-detector <detector> \
  --manifest manifest.tsv \
  --outdir work/<detector> \
  --ref-fa ref/genome.fa \
  --gtf ref/genes.gtf \
  --threads 8 \
  --parallel 4


Example (CIRI-full on chr21 subset):

circyto run-detector ciri-full \
  --manifest manifest_2.tsv \
  --outdir work/ciri-full \
  --ref-fa ref/chr21.fa \
  --gtf ref/chr21.gtf \
  --threads 8 --parallel 1


⚠ Note:
CIRI-full is not thread-safe when running multiple jobs concurrently on the same reference directory.
--parallel 1 is recommended for Codespaces or local testing.

3. Running Multiple Detectors (NEW)

You can now run any number of detectors in one unified command.

circyto run-multidetector <detector1> <detector2> ... \
  --manifest manifest.tsv \
  --outdir work/multi \
  --ref-fa ref/genome.fa \
  --gtf ref/genes.gtf \
  --threads 8 --parallel 1


Example using two detectors on a 2-cell manifest:

circyto run-multidetector ciri-full ciri2 \
  --manifest manifest_2.tsv \
  --outdir work/multi \
  --ref-fa ref/chr21.fa \
  --gtf ref/chr21.gtf \
  --threads 8 --parallel 1

Output Structure
work/multi/
├── ciri-full/
│   ├── ERR2139486.tsv
│   ├── ERR2139559.tsv
│   └── <run dirs + logs>
├── ciri2/
│   ├── ERR2139486.tsv
│   ├── ERR2139559.tsv
│   └── <run dirs + logs>
└── summary.json


summary.json contains a simple machine-readable record:

{
  "ciri-full": {
    "n_cells": 2,
    "detector": "ciri-full",
    "outdir": "work/multi/ciri-full"
  },
  "ciri2": {
    "n_cells": 2,
    "detector": "ciri2",
    "outdir": "work/multi/ciri2"
  }
}

4. Notes on CIRI2 Strictness

CIRI2.pl is designed for bulk RNA-seq and applies very strict filtering, often removing circRNAs even when candidate signals are present.

Common reasons for empty outputs:

small number of junction reads (<2)

insufficient PCC support

multi-mapping reads

single-cell sparse coverage

partial support only from R1/R2

Even with low-stringency flags (-w), many single-cell circRNAs do not pass.

This is expected.

circyto reports these cleanly as header-only TSVs.