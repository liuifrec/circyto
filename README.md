# circCyto

[![CI](https://github.com/liuifrec/circyto/actions/workflows/test.yml/badge.svg?branch=main)](https://github.com/liuifrec/circyto/actions/workflows/test.yml)
[![release](https://img.shields.io/github/v/release/liuifrec/circyto?display_name=tag)](https://github.com/liuifrec/circyto/releases)

<p align="center">
  <img src="assets/circCyto_logo.png" alt="circCyto logo" width="450">
</p>

---

A **Python CLI toolkit** for **single-cell circRNA detection, merging, and conversion**.  
Supports modular detectors (`CIRI-full`, `CIRI-long`, `find_circ`, `CIRCexplorer2`) and  
outputs sparse matrices ready for downstream analysis in **Scanpy**, **scVelo**, etc.

---

## 🧬 Installation

```
git clone https://github.com/liuifrec/circyto.git
cd circyto
pip install -e .
```

Optional extras:

```
pip install -e .[detectors]
```

---

## ⚙️ Commands

| Command | Description |
|----------|--------------|
| `circyto prepare` | 10x BAM → FASTQs (`--chemistry tenx-3p|tenx-5p`) |
| `circyto run` | Run CIRI-full on batches produced by `prepare` |
| `circyto run-manifest` | Run CIRI-full using a manifest TSV (per-cell FASTQs) |
| `circyto collect` | Merge all CIRI outputs into a circ × cell matrix |
| `circyto convert` | Convert `.mtx` + index files to `.loom` or `.h5ad` |
| `circyto make` | All-in-one pipeline (detects manifest or BAM input) |

---

## 📄 Manifest format

```
cell_id    r1                     r2
A01        /data/A01_R1.fastq.gz  /data/A01_R2.fastq.gz
A02        /data/A02_R1.fastq.gz  /data/A02_R2.fastq.gz
```

For single-end, omit `r2`.

---

## 🚀 Example: plate protocols

```
circyto make   --manifest cells.tsv   --ref-fa ref.fa   --gtf genes.gtf   --outdir work   --threads 32   --cmd-template 'CIRI-full -r {ref_fa} -a {gtf} -1 {r1} -2 {r2} -o {out_tsv}'
```

---

## 🧬 Example: 10x Genomics (5′ or 3′)

```
circyto make   --bam possorted_genome_bam.bam   --chemistry tenx-5p   --ref-fa ref.fa   --gtf genes.gtf   --outdir work   --threads 16   --cmd-template 'CIRI-full -r {ref_fa} -a {gtf} -1 {r1} -2 {r2} -o {out_tsv}'
```

> For 10x BAM, R1/R2 reconstruction is heuristic — manifest mode is most accurate.

---

## 🧪 Example workflow (using dummy FASTQs)

```
# 1. Prepare fake FASTQ batches
circyto prepare --bam tests/mini_R1.fastq --outdir tests/fastq_batches

# 2. Run placeholder detector
circyto run   --fastq-dir tests/fastq_batches   --outdir tests/out_ciri-full   --ref-fa tests/mini_R1.fastq   --gtf tests/mini_R2.fastq   --cmd-template 'bash -lc '''printf "circ_id	chr	start	end	strand	support\n" > {out_tsv}''''

# 3. Merge outputs
circyto collect   --cirifull-dir tests/out_ciri-full   --matrix tests/circ.mtx   --circ-index tests/circ_ids.txt   --cell-index tests/cell_ids.txt

# 4. Convert to .loom
circyto convert   --matrix tests/circ.mtx   --circ-index tests/circ_ids.txt   --cell-index tests/cell_ids.txt   --loom tests/circ_full.loom
```

---

## ⚡ Quick-start sanity check

```
bash -lc 'circyto prepare --bam tests/mini_R1.fastq --outdir tests/fastq_batches && circyto run --fastq-dir tests/fastq_batches --outdir tests/out_ciri-full   --ref-fa tests/mini_R1.fastq --gtf tests/mini_R2.fastq   --cmd-template "bash -lc '''printf "circ_id\tchr\tstart\tend\tstrand\tsupport\n" > {out_tsv}'''" && circyto collect --cirifull-dir tests/out_ciri-full --matrix tests/circ.mtx   --circ-index tests/circ_ids.txt --cell-index tests/cell_ids.txt && circyto convert --matrix tests/circ.mtx --circ-index tests/circ_ids.txt   --cell-index tests/cell_ids.txt --loom tests/circ_full.loom && echo ✅ "circCyto CLI OK"'
```

---

## 📂 Project structure

```
circyto/
├── cli/
├── detectors/
├── parsers/
├── pipeline/
├── writers/
├── tests/
└── assets/
```

---

## 🔄 Continuous Integration

Runs on every push to `main`:
- installs package with extras  
- runs CLI smoke tests (`bash tests/test_cli.sh`)  
- verifies core commands (`prepare`, `run`, `collect`, `convert`) work

View → [CI Dashboard](https://github.com/liuifrec/circyto/actions/workflows/test.yml)

---

## 📖 Citation

Liu, Y.-C. *et al.* **circCyto** — a modular toolkit for single-cell circRNA profiling and integration with scRNA-seq workflows.  
*(Manuscript in preparation, 2025).*

---

## 📜 License

<<<<<<< HEAD
MIT License © 2025 **Yusei (James) Liu**
=======
MIT License © 2025 **Yu-Chen (James) Liu**
>>>>>>> origin/main
