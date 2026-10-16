# A04 implementation details fixed before the grid

The full original 23 × 63,187 RNA modality is retained in the derived fixtures;
there is no feature downselection. Original circ annotations are supplied through
the documented matrix route. Long/native invalid cases are separate schema
controls. The clean matrix control and clean synthetic native/zero control are
shared controls for the deterministic scenarios, avoiding redundant clean runs.
The resulting plan is 715 cases × three arms = 2,145 arm attempts.

The unsafe arm deliberately overwrites parsed statuses with completed and assigns
RNA IDs from circ row positions (extra RNA rows receive labelled unsafe IDs).
It shares the reference's strict input validation and output writer. Its invalid
input rejections therefore are inherited validation, not independent evidence
of a heuristic's detection ability. Its incorrect identity/status manipulations
are an illustration and are never attributed to a named package or scientists.

A04 uses a separate extension of the A03 oracle to count missing global membership
rows in deliberately unsafe outputs as errors rather than raising a lookup
exception. A03's already frozen oracle is unchanged. Error counts are assertion/
field violations, not a cell-error probability; affected configurations and exact
candidate numerators/denominators remain available for review. All RNA fields
are compared; RNA value/ID mismatch and missing mapping are reported separately
from unavailable fields. Undefined prevalence is explicitly flagged.

Prevalence tables retain exact error numerator/denominator and separately store
the full-original-cohort numerator and denominator (23). The latter is a
composition reference, not truth for a changed completed cohort. Any display-rate
comparison applies the declared six-place decimal rounding and ignores only
equivalent trailing zeros. Reference, generator and oracle remain separate.

Every case records its selected IDs, exact order/membership changes and all arm
commands. Completed case directories are losslessly tar/gzip archived with a
verified SHA-256 for every member. Only the redundant newly generated directory
is then removed; failed/unexpected cases remain unpacked. This retention policy
does not touch any A02 or historical file. Source/configuration hashes are frozen
in `records/frozen_a04_sources.json` before the first grid configuration.
