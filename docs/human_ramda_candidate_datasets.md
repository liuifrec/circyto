# Human scRR RNA Candidate Datasets

This note resolves public human RNA-side scRepli-RamDA/scRR runs for the next `circyto` validation step.

Scope:

- human only
- RNA-side runs only
- `GSE278952` and `GSE278958`
- first server pilot constrained to 2 cells

Non-goals:

- no DNA-side runs
- no mouse biological validation
- no full production download

## Series-level resolution

### `GSE278952`

- title: `Genome-wide DNA replication profiling and full-length total RNA sequencing from the same single cell [HAP1 mid-S]`
- BioProject: `PRJNA1169834`
- organism: `Homo sapiens`
- assay family: `scRR-seq` / `scRepli-RamDA-seq`
- RNA reference clues from GEO: HAP1 RNA supplementary files are labeled `hg38`
- practical status for `circyto`: public, biologically important, and now aligned with the validated paired-end `STAR+CIRI3` full-length workflow route

### `GSE278958`

- title: `Genome-wide DNA replication profiling and full-length total RNA sequencing from the same single cell [IMR-90 G1]`
- BioProject: `PRJNA1169833`
- organism: `Homo sapiens`
- assay family: `scRR-seq` / `scRepli-RamDA-seq`
- RNA reference clues from GEO: sample processing says `Assembly: hg38`
- practical status for `circyto`: public and immediately usable for a first 2-cell server pilot because confirmed RNA runs are single-end

## Confirmed RNA-side candidates

| dataset_id | gsm | srx | srr | assay_side | organism | cell_type | protocol | read_layout | expected_reference | recommended_route | status | notes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `GSE278958` | `GSM8558852` | `SRX26321174` | `SRR30918126` | RNA | `Homo sapiens` | IMR-90, aphidicolin-treated G1 | scRR-seq / scRamDA-seq RNA side | single-end, 151 bp | `hg38` + human GENCODE v38 | execute with `circyto workflow full-length-circrna` using the validated single-end RamDA route | `ready_for_execution` | Best first pilot cell. SRA XML shows `LIBRARY_LAYOUT SINGLE`, `is_public=true`, one original FASTQ file. |
| `GSE278958` | `GSM8558853` | `SRX26321183` | `SRR30918117` | RNA | `Homo sapiens` | IMR-90, aphidicolin-treated G1 | scRR-seq / scRamDA-seq RNA side | single-end, 151 bp | `hg38` + human GENCODE v38 | execute with `circyto workflow full-length-circrna` using the validated single-end RamDA route | `ready_for_execution` | Good second pilot cell. Same single-end public layout as `GSM8558852`. |
| `GSE278952` | `GSM8558630` | `SRX26315002` | `SRR30911454` | RNA | `Homo sapiens` | HAP1 mid-S | scRR-seq / scRamDA-seq RNA side | paired-end, 2 x 151 bp | `hg38` + human GENCODE v38 | execute with `circyto workflow full-length-circrna --star-index ... --allow-paired-ramda` using the validated STAR+CIRI3 paired-end path | `validated_chr21_subset_route` | SRA XML shows `LIBRARY_LAYOUT PAIRED`, two original FASTQs, `is_public=true`. Local chr21 subset execution on `SRR30911454` completed end-to-end with STAR, BWA rescue, CIRI3, matrix, and h5ad. |
| `GSE278952` | `GSM8558631` | `SRX26315003` | `SRR30911453` | RNA | `Homo sapiens` | HAP1 mid-S | scRR-seq / scRamDA-seq RNA side | paired-end, 2 x 151 bp | `hg38` + human GENCODE v38 | execute with `circyto workflow full-length-circrna --star-index ... --allow-paired-ramda` using the validated STAR+CIRI3 paired-end path | `ready_for_paired_execution` | Slightly smaller than `GSM8558630`; still paired-end. |
| `GSE278952` | `GSM8558632` | `SRX26315028` | `SRR30911559` | RNA | `Homo sapiens` | HAP1 mid-S | scRR-seq / scRamDA-seq RNA side | paired-end, 2 x 151 bp | `hg38` + human GENCODE v38 | execute with `circyto workflow full-length-circrna --star-index ... --allow-paired-ramda` using the validated STAR+CIRI3 paired-end path | `ready_for_paired_execution` | Another public HAP1 RNA run suitable for later paired-end validation at larger scale. |

## Explicit exclusions

Exclude these from the first circRNA pilot:

- `DNA_*` samples from the same series, for example `GSM8558695 / SRX26315019 / SRR30911568`
- replication-only DNA libraries
- any non-RNA assay rows with `Sample_molecule = genomic DNA`, `LIBRARY_SOURCE = GENOMIC`, or `LIBRARY_STRATEGY = OTHER`

## Immediate recommendation

Use `GSE278958` first for the real server pilot.

Reason:

- RNA-side runs are confirmed public
- actual SRA layout is single-end
- single-end RamDA aligns with the already validated `BWA-MEM -> direct SAM -> CIRI3` route
- 2-cell download size is modest relative to the current disk limit

Do not start biological validation with `GSE278952` before the smaller IMR-90 pilot.

Reason:

- the RNA-side HAP1 runs are clearly paired-end in SRA
- the paired-end route now has a real `GSE278952 / SRR30911454` chr21 subset execution proof
- HAP1 remains the intended next human validation target once the smaller IMR-90 server pilot has confirmed the environment

## Storage and execution notes

Observed public object sizes from SRA XML:

- `SRR30918126`: one original FASTQ object about `242 MB`; normalized SRA about `183 MB`
- `SRR30918117`: one original FASTQ object about `248 MB`; normalized SRA about `190 MB`
- `SRR30911454`: original FASTQs about `424 MB` and `473 MB`; normalized SRA about `688 MB`
- `SRR30911453`: original FASTQs about `382 MB` and `425 MB`; normalized SRA about `613 MB`
- `SRR30911559`: original FASTQs about `433 MB` and `478 MB`; normalized SRA about `696 MB`

Practical implication:

- the 2-cell IMR-90 pilot is low-risk for `~290 GB` free disk
- the 2-cell HAP1 paired-end pull is still feasible, and is now a realistic validated-route follow-up after the IMR-90 pilot

## Source notes

Primary records inspected:

- GEO series text for `GSE278952` and `GSE278958`
- GEO sample text for `GSM8558630`, `GSM8558631`, `GSM8558632`, `GSM8558852`, `GSM8558853`
- SRA XML for `SRX26315002`, `SRX26315003`, `SRX26315028`, `SRX26321174`, `SRX26321183`

Key inference boundary:

- use SRA XML `LIBRARY_LAYOUT` as the deciding layout field when generic protocol prose and per-run layout disagree
