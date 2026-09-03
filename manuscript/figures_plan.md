# Application Note figures plan

## Main Figure 1: workflow and Smart-seq3 output

Figure 1 is the only main-paper figure. It should read from architecture to
real-data output in one visual sequence.

### Panel A: circyto architecture

Primary path:

```text
full-length scRNA-seq
    -> protocol-aware alignment / circRNA detection
    -> cell x circRNA matrix + QC
    -> host-gene / known-circRNA annotation
    -> AnnData / MuData
```

Show that detector calls enter circyto; do not depict circyto as inventing a
new back-splice detector. RNA and circRNA are the primary output modalities.
Optional branches for processed CNV, processed RT, and exploratory candidate
variant signals should be smaller, lighter, and labelled optional.

### Panel B: RNA-derived Smart-seq3 UMAP with circRNA burden

- Dataset: E-MTAB-8735, 192 cells.
- Construct the embedding from RNA features only.
- Colour cells by detected circRNAs per cell; include median = 12 only after
  regeneration from the checksum-matched object.
- State in the legend that the circRNA metric is an overlay and did not define
  the embedding.

### Panel C: representative MAN1A2-associated circRNA detection

- Use candidate `chr1:117402186|117420649`, host gene `MAN1A2`.
- Reuse the Panel B RNA-derived coordinates.
- Prefer a binary detected/not-detected overlay for clarity; support may be a
  continuous overlay only if the legend says “support,” not expression.
- Historical export: detected in 6/192 cells, total support 15. Regenerate
  before putting either number in the legend.
- Treat the feature as an illustrative output, not evidence for MAN1A2 biology.

## Supplementary Figure S1: protocol/workflow validation

One compact panel set covering:

- Smart-seq3 pooled/demultiplexed paired-end route;
- IMR90 single-end scRR route;
- HAP1 paired-end scRR route.

For each, show protocol, alignment/detector route, completed output class, and
validation status. Use the validated HAP1 10-cell batch; do not substitute the
unresolved historical 63-cell/RT object.

## Supplementary Figure S2: multimodal object/schema demonstration

Show modality structure and cell-axis contracts rather than biological
associations:

- RNA+circ AnnData/MuData;
- validated IMR90 RNA+circ+processed-CNV MuData (23-cell overlap);
- optional RT import/merge contract, labelled implemented with synthetic tests
  and pending reconciled real processed-file validation.

If space permits, include the scoped `pull_on_update=True` compatibility policy
as a small provenance note, not a separate panel.

## Optional long-read inset

Only add a compact supplementary inset or table row if it does not compete
with the short-read single-cell story. Separate generic ONT cDNA alignment/QC
from the chemistry-gated bulk CIRI-long adapter and retain all non-claims in
`application_note_evidence.md`.

## Deferred biological follow-up

The former HAP1 RT association, IMR90 CNV biology, and cross-dataset host-gene
program figures are not part of this Application Note.
