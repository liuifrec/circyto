# CIRCYTO CMPB A05 review packet

Start with `a05_results_for_review.md`, `cnv_eligibility_decision.md` and `figure2c_rna_circ_only_alternative.png`. All source tables and draft SVG/PNG figures are included. This packet stops after A05. The CNV gate is NOT VERIFIED; no heatmap or fake empty display table is included. Publication figure choice and manuscript prose await review.

Scientific analysis helper: `6ec0ad2f6977a593fabdc5978a3043511214bb61`. Figure renderer helper: `640980f9bf2035a01070a43522f8e8f6b989eae6`. The ZIP's top-level HANDOFF_README.txt records the final committed packet SHA and branch. `SHA256SUMS` in the ZIP covers every other archive member; the outer checksum file covers the ZIP. `inputs_outputs_sha256.json` binds original source hashes, outputs and scripts, explicitly excluding its own hash to avoid self-reference.

To recreate only the vector drafts from the included committed tables, unpack the ZIP, run from its root, and use a new empty destination (Python standard library only):

```bash
python3 scripts/cmpb_a05/plot_drafts.py \
  --tables manuscript/cmpb/a05_descriptive \
  --out /tmp/a05-recreated-figures \
  --renderer-commit 640980f9bf2035a01070a43522f8e8f6b989eae6
```

The SVG/figure-observation output should match the packet byte-for-byte with the same script and pinned renderer commit. Optional PNG conversion uses existing ImageMagick 6.9.11-60 via `python3 scripts/cmpb_a05/render_png.py /tmp/a05-recreated-figures`; raster pixels may vary with fonts/renderer version on other machines. No dependency installation is required to recreate SVGs.

Full scientific recomputation requires the two original checksum-matched H5MU files (not in this ZIP) and the accepted Python/AnnData environment. `analyze.py --original ORIGINAL_LOAD_WORK --out NEW_OUTPUT`, `cnv_gate.py A05_ROOT ORIGINAL_LOAD_WORK NEW_OUTPUT`, and `verify.py --original ORIGINAL_LOAD_WORK --tables NEW_OUTPUT` document the calculation/check interfaces. The exact executed argv, control properties, input URLs/hashes and environment versions are included in execution_details, execution_record.json and source_object_inventory.json. Running resource controls on a new machine requires verifying equivalent limits first; no uncontrolled rerun is recommended by this packet.

Only small A05 outputs/scripts are packaged. The original H5MU, downloaded full CNV source tables, raw reads, native A02 outputs, full A04 prevalence archive, environments, executables, indices and secrets are excluded. Public CNV table URLs and hashes are provided for provenance. The unchanged A02/A03/A04 evidence is linked by accepted commits and recorded preservation checks.
