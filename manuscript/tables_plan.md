# Application Note tables plan

No main-paper table is planned. Keep the supplement to the smallest set needed
to audit datasets, claims, and reproducibility.

## Supplementary Table S1: dataset/workflow inventory

Columns:

- dataset and accession;
- protocol and read layout;
- route;
- reference;
- output object/modalities;
- validation scale;
- manuscript role;
- limitation.

Rows: Smart-seq3 / E-MTAB-8735, IMR90 scRR / GSE278958, and HAP1 scRR /
GSE278952. Long-read examples may be added as clearly separated optional rows.

## Supplementary Table S2: software and reproducibility validation

Columns:

- workflow or capability;
- dataset/test fixture;
- cells;
- circRNAs;
- annotation recovery;
- validation status;
- evidence source;
- limitation.

Use “not regenerated” or “not applicable” rather than importing stale values.
The HAP1 manuscript-scale counts remain blank until reconciled; the validated
10-cell route may be reported as protocol execution evidence.

## Supplementary Table S3: versions, commands, and provenance

Columns:

- component or workflow;
- circyto version and baseline commit;
- Python/dependency or external-tool version;
- command/script;
- input accession/object SHA-256;
- output/provenance record;
- status.

Include the package build, wheel-only install, H5MU round trip, resource lookup,
MuData synchronization audit, and exact manuscript regeneration commands.

## Optional long-read section

If included, append a short, visibly secondary block to S1-S3 rather than
creating additional tables. Report generic SRR4048177 interoperability and the
official CIRI-long demo separately.

## Deferred biological follow-up

Former tables for HAP1 RT regression, IMR90 CNV-high programs, known/novel
biology, cross-dataset host-gene overlap, and local CNV-at-circRNA loci are not
submission requirements.
