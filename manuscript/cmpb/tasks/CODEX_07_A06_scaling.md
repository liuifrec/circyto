# CODEX A06 — Controlled synthetic scaling (CMPB; run only after reading)

**Owner review gate:** This is the next proposed work package, not an instruction to start automatically. The owner starts it by telling Codex to execute this file. Do not execute A07, rerun A02–A05, merge/release/publish anything, or rewrite the circulated manuscript.

## 0. Provenance and boundaries

- Start from **accepted A05 analysis head** \`90381340cc0365bc66c815d06ff01b2978dc3d34\`; the A03/A04 parent is \`02ee3d4acfed63c5b2f02c4baef6d22c822d4dac\`. The task file lives on a separate planning branch. Verify the actual local and origin commit graph before using it.
- Evaluated production implementation: \`e6843a38469d49f96d44758a4e37893ebbe2734b\`, repaired installed wheel SHA-256 \`28c7cbde3c31cf90dd7c247487ba5f410c5b82f6f97d7a2281a0ce52912f7383\`. Its unpublished \`0.10.0\` metadata is **not** the frozen historical public release. Verify installed member digests and \`pip check\`; if a different artifact is present, stop and report before changing the software.
- First read \`manuscript/cmpb/benchmark_plan.yaml\` **A06/scaling**, \`manuscript/cmpb/figure_plan.md\`, \`manuscript/cmpb/a03_a04/results_for_review.md\`, \`manuscript/cmpb/a05_descriptive/a05_results_for_review.md\` and the A05 CNV decision.
- Use a **new separate worktree and analysis branch** \`analysis/cmpb-a06-scaling\` based on A05, leaving existing branch tips, uncommitted work, original H5MU objects, A02–A05 artifacts and protected manuscript/figure/recheck material unchanged. Do not reuse occupied analysis output paths. If that branch exists, inspect and do not overwrite.
- This is a **synthetic import / AnnData–MuData / reporting overhead evaluation**, not raw-read alignment, circRNA detection, clinical accuracy, biological prevalence, human usability, or a named-tool pipeline benchmark. Do not alter production functions, thresholds, scientific normalization or UI to obtain better measurements.

## 1. Fixed workloads: no unapproved expansion

Create deterministic synthetic **native CIRI3-compatible 13-column per-cell files** using the documented native schema, **one file per synthetic cell**, seed \`20261015\`, one synthetic read ID of length 32 bytes (unless explicitly varied), support 1 per nonzero record. Generated identifiers/coordinates/strands/ancillary fields must remain valid for the accepted importer. Every declared feature must occur at least once; all counts, unique IDs, and input UTF-8 bytes must be validated independently. Do not misrepresent synthetic inputs as genuine CIRI3 calls.

| ID | Cells | Circ features | Circ nnz | Additional controlled variable |
| --- | ---: | ---: | ---: | --- |
| S01 | 100 | 1,000 | 10,000 | baseline |
| S02 | 1,000 | 10,000 | 100,000 | baseline |
| S03 | 10,000 | 50,000 | 1,000,000 | **conditional only**, see size gate |
| S04 | 1,000 | 10,000 | 100,000 | 128-byte read IDs, hold other dimensions fixed |
| S05 | 1,000 | 10,000 | 100,000 | add synthetic RNA, 2,000 features, RNA nnz 100,000 |

The controls are **not** a full Cartesian experiment. Each scenario needs actual generated input bytes, source-file count, string/read-ID sizes, circ and RNA axes/nnz, metadata/status declarations and repeat-specific outputs. No raw FASTQs, SAM/BAM, reference download, index build or detector JAR is required or allowed.

## 2. Resource gate — enforce before work

The registration's **upper bounds for the A06 job**, including setup, checks, output and handoff, are **2 CPUs, 8 GiB cgroup RAM including cache, swap disabled, 10 GiB new persistent storage and 2 hours total wall-clock**. These are ceilings, not performance expectations. Start a new explicit time/disk ledger; do not reset previous studies' budgets. Reuse previously tested descendant monitoring, affinity enforcement and real cgroup limits only after small stop-probe verification. An accepted scheduler property does not by itself prove a CPU quota. If caps cannot be enforced, stop before generating the large data.

Execute S01 then S02. **Before generating or running S03**, estimate native-file/parsed string overhead, peak RAM and **all** output/repeat disk bytes from observed S01/S02. If projected memory exceeds 8 GiB, disk cannot stay within the added 10 GiB, or the shared deadline is infeasible, record S03 as \`SKIPPED_RESOURCE_GATE\` with measured evidence. Do not truncate the 10,000-cell workload and call it S03, silently increase limits, or remove previous study assets. Continue smaller eligible S04/S05 only within the same budget. Failure/timeout/OOM are reportable outcomes, not cases to replace.

Keep A06 large inputs/outputs in a **dedicated durable local root outside Git**. Deduplicate generated inputs when safe; retain the original measurements, failures, summary outputs and their hashes. Treat temporary-output cleanup as requiring a narrowly documented reproducibility-preserving policy. Never alter original research/A02–A05 directories.

## 3. Execution and verification

- One separate unmeasured warm-up per configuration, then **three fresh-process measured repeats** in separate output directories, startup included, warm filesystem cache, no privileged cache flush. Interleave ordering if both equivalent operations are timed; do not compare unlike endpoints.
- Measure exact wall seconds, user/system time when feasible, peak process-tree RSS and cgroup peak, input/output bytes, CPU/thread settings, dependency and installed-wheel identifiers, number of files, string widths, sparse shape and nnz, failure/exit code and wall/resource termination.
- Run the same documented \`circyto import-circ\` endpoint for S01–S04, including its output/report contract; S05 adds the explicitly synthetic RNA modality and MuData. Do not silently compare S05 with circ-only scenarios as if only cell count changed. If the importer requires fixed schema/support declarations, register and preserve them, never patch production code to accommodate generated records.
- Verify exact sparse integer values and feature/cell identities using an **independent streaming/CSR oracle** not calling circyto parsing/normalization. Independently check cell reorder, observation status, feature annotations, report denominators, TSV/JSON and HTML-displayed statistics, and ordinary scientific-object round trips. For the largest workload use a justified streaming exact or chunked digest strategy rather than densifying the matrix. A handful of sampled rows is **not** a substitute for a global exact count/hash check.
- Add a test guard against full sparse \`.toarray()\` / \`.todense()\` and inspect potentially dense allocations; document what the guard does **not** prove.
- Report any warning, unsupported shape, incomplete/report withholding or unknown status as actually observed. No invented scalability ceiling, perfect linear claim, or asymptotic estimate from three points.

## 4. Compact manuscript-ready handoff

Produce at least:
\`scaling_runs.tsv\` (all warmup/measured/failed/skipped runs),
\`scaling_inputs.tsv\`,
\`scaling_value_checks.tsv\`,
\`scaling_figure_sources.tsv\`,
\`environment.json\`,
\`source_output_manifest.json\`,
\`resource_usage.json\`,
\`results_for_review.md\`,
and draft Figure 4 scaling SVG/PNG **only if data permit**.

Show **measured** wall time and RSS against both workload and actual input bytes. Separate effects of (i) cells/features/nnz, (ii) read-ID byte length, and (iii) adding RNA. Label every configuration **synthetic**. This benchmark is **not** direct extension of A03's four-real-cell arm-to-arm comparison and cannot be used to claim a speed advantage over the bespoke reference or nf-core/CirComPara2. If S03 is gated, show a visibly missing/gated point, not an interpolated value. Archive all meaningful negative outcomes and avoid selecting favorable repetitions.

Keep analysis scripts, aggregate tables, figures, manifests and small review files on the analysis branch and push only **after** checking the diff; no automated merge or release. Exclude large native files and 10,000-file directories from Git. Make one **verified ZIP under \`~/Downloads/CIRCYTO_CMPB_A06_<timestamp>/\`** containing the results/report/code/plots plus checksums and a README distinguishing complete data from review extracts. Do not package FASTQs, environments, JARs, sensitive local paths or bulk synthetic files. Verify ZIP CRC and member hashes; return the exact ZIP path, SHA-256, final commit, table of **all five scenarios** with status and measured resource statistics, as well as any limitations and next manuscript decision.

**Stop after A06.** A07 publishing, named-tool execution, figure submission and final manuscript rewriting remain owner-review decisions.
