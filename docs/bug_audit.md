# Bug audit: resolved findings

Updated on 11 November 2026. All eight confirmed bug groups from the
post-backend-integration audit have been repaired. The original logs, scripts,
fingerprints, benchmark data, and frozen BURY reference are retained as
historical evidence. The observation scripts sometimes assert the old erroneous
behavior; passing those assertions is no longer the desired result.

The [usage guide](usage.md) describes the current public API and validation
policy. The [audit README](../experiments/bug_audit/README.md) distinguishes the
historical probes from current regressions.

## Repairs

| Finding | Current behavior | Regression coverage |
| --- | --- | --- |
| One-part matching summed an empty `Any[]` | Integer accumulator returns zero, an empty witness, and black edge colors; both matching plotting entry points work | `test/known_bugs.jl`, `test/cost_regressions.jl` |
| Recursive matching traversal overflowed on a 50,000-vertex alternating path | Bipartite sides come directly from part labels; no recursive traversal | Independent cardinality/witness checks in `test/cost_regressions.jl`; bounded `test/run_matching_stress.py` |
| Range search discarded scores above `nv(g)` | Initialize from the first candidate, minimize without a graph-size bound, and retain the last equal-score candidate | `test/bury_validation.jl`, `test/bury_backends.jl` |
| Saran–Vazirani skipped original components | Attempt packing before any cut and record that check in the trace; report solver/final assignment failure explicitly | `test/baseline_legacy_regressions.jl` |
| Legacy constructor aliased its dictionary | Copy the registry and validate counts, exact vertex coverage, and integer labels; preserve empty parts | `test/baseline_legacy_regressions.jl` |
| Legacy displays ignored IO | All three `show` methods and their nested displays use the supplied stream | `test/baseline_legacy_regressions.jl` |
| Invalid targets, partition counts, and empty seed searches behaved inconsistently | Shared validation across backends/balanced entry points, before sampling or solver execution | `test/bury_validation.jl`, `test/baseline_legacy_regressions.jl`, `test/api.jl` |
| Research comparisons mixed tuples/scalars and objectives | Grid and sandbox partitioner comparisons share scalar pairwise GF(2) cut rank for every competitor and label it explicitly | `resources/comparison_scores.jl`, `test/research_scores.jl` |

The Saran–Vazirani regression includes the ten-vertex graph with edges
`(1,2), (2,3), (4,5), (5,6), (7,8), (9,10)`. Its components of sizes
3, 3, 2, 2 now yield a zero-cut 5+5 assignment. Edgeless graphs and `k=1`
are also covered. This repair does not promise optimal heuristic results.

## Intentional boundary changes

- Every BURY backend requires integer targets in `0:nv(g)`. Legacy no longer
  accepts integral floats such as `2.0`, and oversized targets now raise
  `ArgumentError` rather than entering a nonterminating loop.
- Seeds are integers in `1:nv(g)`. Seed searches and seeded BURY reject empty
  graphs even for target zero, before consuming any RNG draws. Ordinary BURY
  on an empty graph with target zero remains supported.
- Range bounds must be integers with `0 ≤ min_v ≤ max_v ≤ nv(g)`, validated
  before callbacks. Built-in callbacks receive the backend; custom callbacks
  retain their two-argument contract.
- Balanced generation requires a nonempty graph and integer `k` in `1:nv(g)`
  dividing `nv(g)`. Invalid calls raise `ArgumentError` before RNG/solver use.
- Legacy `quantumLand` counts require positive integer hamlet count and
  nonnegative integer vertex count. The dictionary must contain exactly
  integer keys `1:n`, with integer labels in `1:k`. Empty membership with
  positive `k` and declared empty parts are supported.
- Saran–Vazirani traces now begin with the original component-packing attempt.
  Matching witnesses may select different optimal edges; cardinality and
  pairwise matching validity define the contract.

Supported integer BURY choices, deterministic ties, seed policies, return
formats, verbose traces, and valid-call RNG behavior are checked against the
unchanged `test/reference/bury_original.jl`. The multipartition driver and
legacy algorithm bodies have not been rewritten.

## Validation and limits

The five former `@test_broken` cases are now ordinary passing regressions.
The package suite adds independent exhaustive matching checks through four
vertices, including one part, disconnected graphs, isolated vertices, and
empty declared parts. It checks cardinality and witness validity instead of
requiring identical selected edges.

Run the package suite with:

```sh
julia --startup-file=no --project=. -e 'using Pkg; Pkg.test(; allow_reresolve=false)'
```

The optional large matching regression runs in a disposable Julia process
with external startup and per-stage watchdogs:

```sh
python3 test/run_matching_stress.py
```

Validation completed on Julia 1.12.6:

- The standard `Pkg.test(; allow_reresolve=false)` command passed **457,209
  assertions**, with no broken tests. This includes **379,348** comparisons
  against the frozen original, **67,734** independent matching regressions,
  **9,191** cost-semantics checks, and **166** packed-rank regressions. The final
  restoration of the unchanged legacy edge-cut loop after benchmarking was
  separately checked with all **9,191** cost-semantics assertions.
- The existing independent cost audit passed again: **127,098 assertions on
  7,122 fixtures**, using exhaustive matching and independent GF(2) XOR
  elimination. New output was saved outside the repository; historical logs
  were not overwritten.
- The durable bounded stress runner passed: the alternating 50,000-vertex
  path returned cost **25,000** and a witness of 25,000 valid edges with
  50,000 distinct endpoints. Scalar and witness calls took approximately
  **18.7 s** and **19.1 s** after the bug fixes alone. With the subsequent
  bulk-construction/unweighted-matching optimization, the same bounded check
  passed in approximately **15 ms** and **19 ms**. These are single warmed
  checks, not benchmark medians. The initial pre-optimization five-second
  large-case probe timed out before verification under a longer bound.
  No larger case was attempted.
- Tiny research-score equivalents passed scalar, `mean`/`std`, objective
  consistency, and plotting checks. `Statistics` was declared as a test-only
  dependency; runtime dependencies, versions, and the manifest are unchanged.
- `git diff --check` passed. The frozen reference and saved historical data,
  logs, and fingerprints remain unchanged.

The original saved 407 API observation assertions describe the pre-fix
snapshot and are not a substitute for the new desired-behavior regressions.

Remaining limitations: BURY and Saran–Vazirani remain heuristics; greedy
vertex-cover nonoptimality is intentional. The old vertex-cover plotting
callback remains incompatible because its second return value is a numeric
estimate, not edge colors. Legacy nested membership objects remain mutable
and should be treated as read-only. Large matching solver workloads can still
be slow despite eliminating recursion. No complete research sweep, large
automatic-seed search, or `resources/graphml*` workload was launched.

The subsequent [cost performance experiment](../experiments/cost_performance/README.md)
records the optimized evaluator measurements against a separate, unchanged
post-bug-fix snapshot. This preserves the distinction between correctness
repairs and later performance changes.
