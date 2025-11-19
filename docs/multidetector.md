# Multidetector Workflow (v0.7.0)

The v0.7 multidetector pipeline enables:
- Running multiple circRNA detectors in one command
- Merging their outputs into a unified circRNA union matrix
- Computing detector similarity and sensitivity metrics

---

## 1. Run multiple detectors

```bash
circyto run-multidetector ciri-full ciri2 \
  --manifest manifest.tsv \
  --outdir work/multi \
  --ref-fa ref/genome.fa \
  --gtf ref/genes.gtf \
  --threads 8 --parallel 1
```

Produces:

```
work/multi/
  ├── ciri-full/
  │     ├── <cell>.tsv
  │     └── <run dirs + logs>
  ├── ciri2/
  │     ├── <cell>.tsv
  │     └── <run dirs + logs>
  └── summary.json
```

---

## 2. Merge detector outputs

```bash
circyto merge-detectors work/multi work/multi_merged
```

Produces:

```
work/multi_merged/
  ├── circ_union.tsv
  ├── circ_by_detector.tsv
  └── metadata.json
```

---

## 3. Compare detectors

```bash
circyto compare-detectors work/multi_merged work/multi_compare
```

Produces:

```
work/multi_compare/
  ├── jaccard.tsv
  ├── detector_summary.tsv
  └── compare_metadata.json
```

This supports benchmarking, agreement scoring, and detector sensitivity analysis.

---

## Notes
- `ciri2` defaults to relaxed stringency (`-0`) for sparse scRNA-seq.
- `ciri-full` uses the adapter wrapper and requires proper reference indices.
