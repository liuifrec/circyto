# Owner-authorized time reallocation — 2026-10-16

The owner authorized reallocating 20 minutes of unused A03 allowance to A04.
The A04 deadline changes from 2026-10-16T04:51:18.414622Z to
2026-10-16T05:11:18.414622Z. Its original start remains 02:51:18.414610Z.
The overall session deadline remains 05:39:09.646192Z. This is not a fresh clock
or an added overall session window. CPUs, memory, no-swap policy, 16 GiB new
storage and overall 80 GiB ceiling remain unchanged.

Preserve the original A04 budget bytes and record old/new hashes in
`${EVAL_ROOT}/records/time_amendment.json`. The active continued stage retains
its original launch timer. If it reaches that timer with cases pending, retain
its complete and interrupted evidence, then run only pending registered cases
in a separate `a04_final` output tree under the amended deadline. A brief recorded
checkpoint just before that timer may allow open output streams to close.
The scientific converters, generator, oracles, inputs and grid do not change.

Completion still requires final output and preservation checks within the
amended allowance. Budget exhaustion remains an incomplete result, not a pass.

The continued grid was administratively checkpointed at 04:08:09 UTC (exit 130),
before its old timer expired. The final segment uses a distinct `a04_final`
output tree and skips only fully archived, three-arm completed cases from both
prior segments. Partial attempts are retained and repeated as a whole case.
Only orchestration paths/checkpoint selection changed; converter, perturbation
generator, oracle and scientific case definitions remain frozen.
