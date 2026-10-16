# Additional output verification — 2026-10-16

Recorded during review of the running frozen grid, before the additional check
is executed. The frozen per-record oracle already compares all sparse AnnData
counts, source annotations, RNA fields/maps, QC TSV/JSON rates and displayed
first-100 HTML table values. It does not independently parse the separate output
Matrix Market file or every prose summary paragraph. Those are a coverage limit,
not an observed production or reference defect.

`output_audit.py` adds a separate standard CSV, exact-integer Matrix Market and
HTML parser over retained outputs. It checks output matrix axes/every integer,
JSON counts/states/rates/declarations/overlap, HTML summary counts/denominators
and applicable declarations/overlap, and absence of remote rendering resources.
It reads only the RNA observation axis needed for overlap; the frozen oracle
retains responsibility for every RNA value/field. Reference and unsafe HTML have
a simpler summary format; unavailable prose sections are not called wrong values.

The added checker is tested with tiny valid and deliberately corrupted generated
outputs and its source hash is recorded before use. It does not change either
converter, their timing, the perturbations, frozen source hashes, expected grid
outcomes or earlier results. Every retained valid-arm output is eligible; there
is no outcome-based selection. A new evaluated-arm inconsistency stops the
additional check and remains a failure for review. Unsafe illustrative mistakes
remain recorded. This is outside measured A03 execution for both arms.

The same resource/time ceilings apply. If insufficient allowance remains, this
additional coverage stays explicitly incomplete; it will not be asserted as
passed. Any interrupted archive or output remains retained without repair.
