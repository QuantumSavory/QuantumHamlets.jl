# Bipartition target-range performance

This bounded experiment compares the complete ordinary-BURY range workflow:

1. The legacy tuple wrapper at commit `a05adb5`, before the range changes.
2. The current `bury_not_balanced_bipartition` tuple wrapper.
3. The new `bipartition_graph` convenience API returning `GraphPartition`.

`reference_range.jl` contains the exact old function body, isolated in a module.
Its SHA-256 is `45025a957267b2ad048dd055fa66cb82be240ac71516ecf3af473c0af013a0fd`.
The reference uses the same current BURY algorithms and previously optimized
cost evaluators as production. This isolates range selection and representation
changes rather than repeating the earlier cost-kernel comparison. It does not
alter historical benchmark records or dependency versions.

Run from the package root, after other tests and benchmarks have stopped:

```sh
python3 experiments/bipartition_range_performance/run.py
```

The runner uses the existing `experiments/bury_performance` environment and
creates a new timestamped result directory, refusing to overwrite a directory.
An explicit fresh destination can be supplied with `--output-dir PATH`.
The external watchdog permits 180 seconds for loading, 60 seconds per comparison
(including compilation), and 480 seconds overall. A warmed call over two seconds
aborts rather than repeating the workload. No GraphML resources, automatic seed
searches, approximation prototypes, or research sweeps are run.

## Method

The fixed fixtures are an 11-vertex path, a 15×15 grid, a 201-vertex sparse random
graph (`p=0.025`), and a 61-vertex dense random graph (`p=0.35`). Random graph
generation uses a local seeded `Xoshiro` and is outside measurement. Each graph
uses a singleton target, the two nearest-balanced odd targets, and six consecutive
targets from `fld(n,2)-2` through `fld(n,2)+3`.

All three backends and matching, edge-cut, and rank objectives are measured:
108 built-in comparisons. Another 12 comparisons wrap the heap partitioner and
matching tuple objective in closures to exercise the legacy generic path, and
wrap the helper's scalar objective in a closure. These closures are called on
every candidate; the helper still uses its ordinary-BURY implementation.

All calls include every target's partition generation, scoring, winner
selection, and returned representation. The two wrappers return dictionaries
and Boolean vectors; the helper returns `GraphPartition`, so its comparison
also reflects the different required output construction. All methods must
produce identical memberships, scores, and winning targets, and leave the
graph unchanged. Legacy output types and registry/color agreement are checked.

Each method is warmed and calibrated independently to approximately two
milliseconds per timing sample, between 1 and 1,000 evaluations. One explicit
full garbage collection precedes each comparison. Nine rounds rotate the three
methods through all timing positions; per-sample forced collection is disabled,
while normal automatic GC remains enabled. The reported times, allocated bytes,
and allocation counts are medians of those nine samples. Compilation, graph
generation, verification, and benchmark construction are excluded.

`results.tsv` contains the medians and `samples.tsv` the individual samples.
`output.log` records the Julia version, CPU, thread count, source hashes, progress,
and both tables. Allocation counters measure Julia allocation traffic rather
than peak memory, and do not fully count native igraph allocations. Speedup is
the old wrapper's median divided by the new method's median; below one means
a slowdown. These controlled fixtures and nine samples do not establish a
universal speedup or statistical confidence interval.

## Results (11 November 2026)

The complete run is [run_20261111T212033Z/results.tsv](run_20261111T212033Z/results.tsv),
with [individual samples](run_20261111T212033Z/samples.tsv) and
[metadata/progress](run_20261111T212033Z/output.log). All **120 three-way
comparisons** passed. Measurements used Julia 1.12.6, one thread, and the Apple
M1 CPU target, with no other package tests or benchmarks running concurrently.

Selected median elapsed times are in **microseconds**; “before” and “after”
are the legacy tuple wrapper, and “helper” is `bipartition_graph`.

| Fixture and targets | Backend | Objective | Before μs | After μs | Helper μs |
| --- | --- | --- | ---: | ---: | ---: |
| Path11, `5:6` | scan | Matching | 8.69 | 5.82 | 5.05 |
| Path11, `5:6` | scan | Edge cut | 3.04 | 0.855 | 0.415 |
| Grid225, `112:112` | heap | Matching | 49.4 | 27.8 | 16.2 |
| Grid225, `112:113` | heap | Matching | 96.4 | 42.7 | 31.2 |
| Grid225, `110:115` | heap | Matching | 306 | 107 | 92.8 |
| Grid225, `110:115` | heap | Rank | 242 | 75.8 | 64.8 |
| Sparse201, `98:103` | heap | Edge cut | 237 | 59.6 | 47.3 |
| Dense61, `30:31` | heap | Matching | 91.3 | 67.1 | 63.7 |
| Dense61, `30:31` | heap | Rank | 38.7 | 26.5 | 24.3 |

Across the 108 recognized built-in comparisons, the wrapper improved
**1.30–4.24× with heap**, **1.33–5.61× with scan**, and **1.01–1.62× with
legacy**. The helper improved **1.42–5.32×**, **1.44–7.81×**, and **1.30–2.15×**,
respectively, relative to the old wrapper. Every built-in median was faster in
this run; near-unity differences are not evidence of a reliable general gain.

Allocation reductions include the returned representation. For grid225's
six-target matching search with heap, Julia allocated bytes fell from **464,816
to 103,808** for the wrapper and to **74,352** for the helper; allocation counts
fell from **6,889 to 618 and 208**, respectively. For the singleton matching
case, the same comparison was **77,488 → 45,280 → 15,824 bytes** and
**1,149 → 448 → 38 allocations**. Sparse201's six-target edge-cut search used
**412,256 → 65,440 → 35,184 bytes** and **6,283 → 464 → 54 allocations**.

There are two measured tradeoffs. First, the generic legacy callback path is
**2.0–12.0% slower** in all 12 callback comparisons, because retained results
now require independent copies. On path11's two targets, it changes from
**8.66 to 9.70 μs**, **12,288 to 15,280 bytes**, and **201 to 222 allocations**.
On grid225's six targets it changes from **278 to 305 μs**, **464,816 to
573,120 bytes**, and **6,889 to 6,962 allocations**. This is the cost of fixing
shared-buffer ownership; arbitrary callbacks still receive the legacy objects
and are invoked for every target.

Second, keeping the legacy algorithm intact means it still constructs a
registry for each candidate, while the optimized selection path builds the
returned registry from the winning membership. Six legacy-backend singleton
rows (all three objectives on grid225 and sparse201) therefore allocate more
bytes despite fewer allocation events and lower measured elapsed times. The
largest increase is sparse201 rank: **105,168 to 112,112 bytes** (**6.6%**),
with **1,365 to 1,155 allocations**. The helper uses **81,856 bytes** there.
All rows, including callback slowdowns and allocation increases, remain in
the TSV. Seeded search policies and their performance were outside this pass.


## Package validation

After measurement, the required package command
`julia --startup-file=no --project=. -e 'using Pkg; Pkg.test(; allow_reresolve=false)'`
passed **459,086 assertions** with no broken tests, including **1,877** new
target-range assertions. The [full output](validation_pkg_test.log) is retained
separately from timing evidence. The suite includes frozen-original BURY
compatibility, cost semantics, and independent matching cardinality checks.
