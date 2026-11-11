# Cost evaluator performance comparison

This experiment compares the optimized evaluators with the **post-bug-fix,
pre-optimization** `src/costs.jl` snapshot in `reference_costs.jl` (SHA-256
`0932a6ebe59db1a46a8fe7586b772d6600b89663f67b663d25373ea6a7f52211`). It does not
modify the older BURY benchmark records or audit logs.

Run from the package root:

```sh
python3 experiments/cost_performance/run.py
```

The runner uses the existing `experiments/bury_performance` environment; no new
dependency installation is needed. Each invocation creates a fresh timestamped
result directory and refuses to overwrite an existing directory. An explicit
new destination can be supplied with `--output-dir PATH`. To remeasure just one
changed evaluator, use e.g. `--operation edge_cut_legacy`; repeat `--operation`
to select several. The default still measures all eight API variants.

`results.tsv` contains numerical results, and `output.log` also records the
Julia version, CPU, thread count, source hashes, and progress markers. The
external watchdog allows 180 seconds for package loading, 60 seconds per
operation (including its first compilation), and 480 seconds overall. A warmed
call over two seconds aborts instead of repeatedly launching a slow workload.
No GraphML resources, automatic seed searches, or large research sweeps run.

## Method

The fixed fixtures are a six-vertex path, an alternating 1,000-vertex path,
two 400-vertex grid memberships, a seeded sparse 300-vertex random graph, a
seeded dense 120-vertex random graph, a twelve-vertex graph with fourteen
empty declared parts, and an empty graph with three parts. They cover sparse,
dense, balanced, uneven, multipartite, disconnected-boundary, and empty inputs.

Every fixture measures all four original tuple-returning evaluators and all
four scalar/named-tuple `GraphPartition` cost APIs. The baseline scalar adapters
include their original `quantumLand` conversion and unused color construction;
legacy evaluator timings receive a preconstructed land on both sides.

Each pair is warmed and calibrated separately to target approximately two
milliseconds per timing sample, with between 1 and 1,000 evaluations per
sample. Benchmark objects are constructed once per pair. One explicit full garbage
collection precedes measurement; per-trial and per-sample forced collections
are disabled, while normal automatic GC remains enabled. Nine rounds alternate
before/after order. Reported times, allocation
counts, and bytes are medians across those nine samples. Graph generation,
partition/land setup, verification, and compilation are outside measured calls.
Allocation data counts Julia allocations; native igraph/Nemo allocations are
not fully represented in those counters.

Every measured operation checks exact numerical agreement with the snapshot
and graph/partition preservation. Edge-cut colors are compared exactly;
matching witnesses are checked for cardinality, valid cross-part edges, and
pairwise disjoint endpoints, allowing different optimal selected edges.
`speedup` is `before_ns / after_ns`; values below one mean a slowdown. These
small controlled fixtures do not establish a universal performance guarantee.

## Files

- `reference_costs.jl`: unchanged evaluator snapshot.
- `reference.jl`: isolated reference module and original scalar adapters.
- `benchmark.jl`: fixtures, correctness checks, and interleaved measurements.
- `run.py`: external watchdog and non-overwriting output capture.

## Results (11 November 2026)

The completed eight-operation run is [run_20261111T195403Z/results.tsv](run_20261111T195403Z/results.tsv),
with metadata and progress in [output.log](run_20261111T195403Z/output.log).
All **64 before/after pairs** passed the numerical and preservation checks.
Measurements used Julia 1.12.6, one thread, and the Apple M1 CPU target, after
the package tests and independent oracle worker had finished.

Selected medians below are **microseconds**. “Legacy” includes the original
return format; “scalar” includes the public `GraphPartition` adapter.

| Evaluator | Fixture | Before μs | After μs | Speedup |
| --- | --- | ---: | ---: | ---: |
| Matching, legacy | Alternating path, 1,000 vertices | 8,737 | 177 | 49.2× |
| Matching, scalar | Four-stripe grid, 400 vertices | 3,609 | 21.2 | 170× |
| Matching, legacy | Dense graph, 120 vertices | 18,181 | 230 | 79.1× |
| Rank, legacy | Four-stripe grid, 400 vertices | 22,939 | 6.32 | 3,631× |
| Rank, legacy | Dense graph, 120 vertices | 1,842 | 51.2 | 36.0× |
| Edge cut, scalar | Alternating path, 1,000 vertices | 44.6 | 4.48 | 9.94× |
| Cover estimates, scalar | Four-stripe grid, 400 vertices | 53.2 | 22.3 | 2.39× |
| Cover estimates, legacy | Dense graph, 120 vertices | 161.5 | 95.3 | 1.69× |

Across the five moderate fixtures (excluding the tiny and empty cases), legacy
matching improved 19.5–136× and legacy rank improved 36–4,549×. Scalar edge-cut
improved 2.79–9.94× with **zero measured Julia allocations** on every fixture.
Scalar cover estimates improved 1.24–2.39×.

The changes are not uniformly faster. In the initial eight-operation run,
legacy edge-cut was approximately unchanged overall and was 8.2% slower on
the interleaved grid (6.02→6.52 μs). Its helper-based implementation was then
reverted to the exact frozen function body; those legacy edge-cut rows are
superseded by the focused follow-up below. Scalar edge-cut remains optimized,
and all other production evaluator code is unchanged from that full run.

Legacy cover was 5.5% slower on the alternating path (136→143 μs), 3.9% slower
on the sparse random graph (61.8→64.2 μs), and 3.1% slower on the empty-parts
fixture (1.276→1.315 μs). Cover allocations also increased in bytes on the
sparse/dense fixtures even when allocation counts or runtime decreased. The
full TSV retains these results rather than reporting only improvements.

For the 1,000-vertex path, rank's legacy Julia allocations dropped from
approximately 82.3 MB and 2.01 million allocations to 155 kB and 34 allocations.
Matching's legacy allocations dropped from 451 kB to 191 kB there. Native
allocation sizes are not included; dense-case matching does not reduce Julia
allocated bytes in this run despite its large runtime improvement.

`run_20261111T195153Z` is an **incomplete preliminary run**, deliberately
interrupted after five comparisons to remove repeated benchmark construction
and forced-GC overhead from the harness. It is retained as evidence but is
not used in this summary. The completed run uses the method described above.


### Focused follow-up after restoring the legacy edge-cut body

[run_20261111T195712Z/results.tsv](run_20261111T195712Z/results.tsv) remeasures
**only** `edge_cut_legacy` across the eight existing fixtures after restoring
that function's exact pre-optimization body. All eight numerical/preservation
checks passed, with identical allocation counts and bytes for each pair.
No rank, matching, cover, or scalar workload was repeated for this follow-up.

No legacy edge-cut speedup is claimed. Most observed median ratios were
0.98–1.02×; the sparse fixture was 1.08× and the dense fixture was 0.81×
(26.1→32.3 μs), despite the identical function bodies. The experiment does not
estimate statistical confidence or independently repeat runs, so these
residual timing differences are retained without assigning a statistical
cause. The scalar edge-cut improvement comes from its separate direct-label,
allocation-free path and is unaffected by the restoration.

The final source hashes are recorded in the focused follow-up's
[output.log](run_20261111T195712Z/output.log). The original reference snapshot
hash is unchanged; the other seven evaluator/API variants are identical to
the completed eight-operation run.

## Correctness and bounded stress check

The full package suite passed **457,209 assertions**, with no broken tests.
After restoring the original legacy edge-cut body, the focused **9,191**
cost-semantics assertions passed again. The independent cost audit passed
**127,098 assertions on 7,122 fixtures** with exhaustive matching and separate
GF(2) elimination. The new rank suite covers multiple UInt64 words, rectangular
and deficient matrices, empty parts, and large declared partition counts.

The [bounded stress log](matching_stress.log) records the 50,000-vertex
alternating path: matching cost 25,000 and a valid 25,000-edge witness, in
approximately 15 ms for the scalar call and 19 ms for details. These are single
warmed validation calls, not benchmark medians. The corresponding pre-optimization
checks took approximately 18.7 s and 19.1 s. No larger case was attempted.
