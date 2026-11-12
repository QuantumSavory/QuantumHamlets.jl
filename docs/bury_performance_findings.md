# BURY performance findings

Measured on 11 November 2026 against the working tree at that time. This is a
historical prototype study: references to “current” code, source line numbers,
and proposed changes describe that snapshot. See the [usage guide](usage.md)
for the maintained implementation and the [archive index](../archive/README.md#bury-performance)
for the complete study and restoration instructions.

Yes: the BURY implementations can be substantially faster while preserving the tested partitions. The strongest result is an indexed-heap implementation for larger sparse graphs. On a 4,096-vertex grid it reduced ordinary BURY from **9.44 ms to 0.350 ms** and explicitly seeded BURY from **29.4 ms to 0.444 ms**. A simpler scan implementation is also useful, particularly for seeded BURY and dense graphs. The heap is not uniformly better: it was more than twice as slow as the current code on the tested dense random graph.

I recommend retaining both approaches for further evaluation: start with the simpler seeded scan, then consider a heap path for larger sparse workloads. Do not replace every selector with a heap or choose a density threshold from these few fixtures alone.

This study adds only new experimental files and this report. It does not install the prototypes into the library or fix the previously identified bugs. The existing source, tests, documentation, `Project.toml`, and `Manifest.toml` are preserved. The experiment includes hashes of all 187 pre-existing files and a verification script.

## What was measured

The reference is the current working tree, including the earlier source reorganization, rather than an older Git checkout. The experiment uses Julia 1.12.6, Graphs 1.13.4, BenchmarkTools 1.6.3, one Julia thread, and an Apple M1 Pro with 16 GiB of memory. Its separate environment preserves the dependency versions from the existing package manifest.

Each competitor is warmed before timing. Three rounds rotate measurement order, with one call per sample and the RNG reset in untimed setup. Graph generation, GraphML loading, and plotting are outside the timings. Ordinary and explicitly seeded timings include final dictionary and color-vector construction. Multipartition timings include graph extraction and registry assembly. Seed/range-search timings include the existing matching-cost computation.

Tables report pooled sample medians across the three rotated rounds. The broad run covers 53 workloads and 159 implementation/workload combinations, allowing 101 samples or 0.4 seconds per trial; raw results include actual sample counts and quartiles. The largest grid receives a separate longer run because the short budget admits too few slow baseline calls. Allocation measurements describe bytes allocated per call, not peak resident memory. These are local workstation measurements, not performance guarantees for other machines or graph families. The warmup, sampling, and allocation distinctions follow the [BenchmarkTools guidance](https://juliaci.github.io/BenchmarkTools.jl/stable/manual/).

## Single bipartitions

Every row selects half the vertices, rounded down. “Seeded” uses explicit seed vertex 1, so it excludes automatic seed search. Times are milliseconds; speedup is current time divided by heap time.

| Graph | Variant | Current | Scan prototype | Heap prototype | Heap speedup |
| --- | --- | ---: | ---: | ---: | ---: |
| 6×6 grid, 36 vertices | Ordinary | 0.00200 | 0.00104 | 0.00117 | 1.7× |
| 6×6 grid, 36 vertices | Seeded | 0.00317 | 0.000791 | 0.00158 | 2.0× |
| 64×64 grid, 4,096 vertices | Ordinary | 9.44 | 10.1 | 0.350 | 27× |
| 64×64 grid, 4,096 vertices | Seeded | 29.4 | 4.33 | 0.444 | 66× |
| 128×128 grid, 16,384 vertices | Ordinary | 150 | 165 | 1.80 | 83× |
| 128×128 grid, 16,384 vertices | Seeded | 452 | 65.7 | 2.23 | 203× |
| Random 3-regular graph, 4,096 vertices | Ordinary | 7.06 | 7.56 | 0.396 | 18× |
| Random 3-regular graph, 4,096 vertices | Seeded | 20.6 | 3.91 | 0.535 | 38× |
| QAOA fixture `qaoa_60`, 1,890 vertices | Ordinary | 1.95 | 2.03 | 0.168 | 12× |
| QAOA fixture `qaoa_60`, 1,890 vertices | Seeded | 5.97 | 1.14 | 0.287 | 21× |
| Dense random graph, 512 vertices, p=0.25 | Ordinary | 0.0901 | 0.0581 | 0.190 | 0.47× |
| Dense random graph, 512 vertices, p=0.25 | Seeded | 0.149 | 0.0424 | 0.364 | 0.41× |
| Complete graph, 512 vertices | Ordinary | 0.382 | 0.230 | 0.393 | 0.97× |
| Complete graph, 512 vertices | Seeded | 0.657 | 0.136 | 0.683 | 0.96× |

The simpler prototype is about 6.8× faster than the current seeded implementation on the 4,096-vertex grid, but approximately 7% slower for ordinary BURY there. It still performs a full minimum scan each round. Removing allocations and writing a manual minimum loop do not by themselves remove that dominant quadratic work. Conversely, on the dense random graph, the scan is about 1.6× faster for ordinary BURY and 3.5× faster for seeded BURY; maintaining heap order after many neighbor updates is more expensive.

The 128×128 rows use the additional four-second trials. The seeded baseline has 26 timed samples, with an interquartile range of 430–476 ms; the heap has 303 samples, with an interquartile range of 2.217–2.248 ms. Its cumulative allocation traffic falls from 2.14 GB per call to 1.23 MB. This is a strong result for that workload, not a claim that every BURY call becomes 203× faster.

The bundled QAOA filename refers to circuit size, not the number of graph vertices. The experiment reads the declared node order and uses that same graph for every competitor. All three implementations produce identical membership on each timed input.

## Complete multipartition and search calls

These times include the full workflow, in milliseconds. Multipartition prototype rows combine the revised driver with the selected kernel. Searches retain the original scoring and candidate policy; typed registry construction is also included in their improvement.

| Workload | Current | Scan prototype | Heap prototype | Heap speedup |
| --- | ---: | ---: | ---: | ---: |
| 64×64 grid, ordinary `bury_k`, k=2 | 15.0 | 10.5 | 0.670 | 22× |
| 64×64 grid, ordinary `bury_k`, k=8 | 13.7 | 13.0 | 3.09 | 4.4× |
| 64×64 grid, ordinary `bury_k`, k=16 | 16.3 | 14.7 | 5.88 | 2.8× |
| 64×64 grid, seeded multipartition, k=8 | 37.1 | 6.15 | 3.41 | 11× |
| QAOA fixture, 1,890 vertices, ordinary k=7 | 4.10 | 3.57 | 1.63 | 2.5× |
| QAOA fixture, 1,890 vertices, seeded k=7 | 9.71 | 3.07 | 0.558 | 17× |
| 16×16 grid, default automatic seed selection and final BURY | 5.71 | 3.86 | 3.94 | 1.45× |
| 8×8 grid, `best_seed` over all 64 seeds | 4.46 | 3.88 | 4.08 | 1.09× |
| 8×8 grid, target search over 16:32 | 1.26 | 1.17 | 1.17 | 1.07× |

Large single-kernel speedups do not automatically carry over to many-part or search workflows. More parts require more graph extraction and mapping work, while seed searches repeatedly evaluate matchings. In the automatic-seed profile, matching evaluation appears in more sampled seed-search frames than the candidate BURY call. Improving BURY alone therefore gives a much smaller end-to-end gain there. The scan is slightly faster than the heap on all three small search fixtures.

## Why the implementation can be faster

### Avoid repeated counting and temporary vectors

In the measured `src/bury.jl` snapshot, ordinary BURY repeatedly calls `sum(colored)` and `findmin(weights)` at lines 62–63. Seeded BURY does this at lines 252–255 and additionally constructs `.!region * Inf + weights` before selecting a vertex. Those full-size temporary arrays are especially expensive when many rounds are needed.

The [archived scan prototype](../archive/README.md#bury-performance) (`experiments/bury_performance/scan.jl`) maintains the colored count incrementally, scans eligible vertices directly, visits the chosen vertex and neighbors without `vcat`, and uses explicit neighbor-update loops. It also uses typed, pre-sized registries. The remaining budget stays separate from the colored count, as in the reference implementation.

The sampling profiles support this diagnosis: the ordinary 128×128-grid run spends nearly all its sampled BURY frames under the minimum-search line, while the seeded run concentrates at the masked-selection expression, including array arithmetic and allocation. The profiles are supporting evidence; the timing tables come from separate unprofiled runs.

On the 4,096-vertex seeded grid, allocation traffic falls from **132.5 MB per call to 0.181 MB** for the scan and **0.308 MB** for the heap. Allocation counts fall from **36,591 to 16 and 25**, respectively. These figures do not imply that the heap has the smallest live workspace: it keeps additional integer arrays for heap positions and regional candidates.

### Replace repeated global searches with indexed heaps

The [archived heap prototype](../archive/README.md#bury-performance) (`experiments/bury_performance/heap.jl`) keeps the next candidate in an indexed binary min-heap ordered by **weight, then vertex ID**. Weight decreases update the heap incrementally. Seeded BURY uses a second heap for regional candidates and falls back to the global heap only when no unprocessed regional candidate remains.

For a single valid bipartition, repeated full scans give the current implementation O(n²+m) worst-case work. A vertex can be processed at most once, and each newly colored vertex updates its own weight and its neighbors. An indexed heap changes this to O((n+m) log n), with O(n) additional storage. The asymptotic advantage is compelling for large sparse graphs, but dense graphs pay for many logarithmic updates. The measured dense-graph slowdown is consistent with that tradeoff.

Two details are essential to correctness. First, **colored does not mean processed**: a colored neighbor must remain eligible until it is selected as the choice. Removing every colored vertex would change the algorithm. Second, ties must retain the smallest current vertex ID. Both prototypes preserve these rules and the seeded region's global fallback.

### Avoid unnecessary multipartition work

The [archived multipartition driver](../archive/README.md#bury-performance) (`experiments/bury_performance/multi.jl`) uses a color-only kernel, so intermediate rounds do not construct registries that `bury_k` immediately discards. It carries the active-to-original vertex map directly, extracts the next graph from the current remainder, and assigns the last remainder without another BURY call or another graph rebuild.

Vertices remain in ascending original order, preserving the reference tie behavior after relabeling. The scan and heap rows for multipartitioning include both this driver and their respective kernel; they are not isolated measurements of the heap alone.

## What the tests establish

The standalone suite passed **253,988 assertions**. It covers all 1,099 labeled simple graphs with 1–5 vertices and 14 larger fixtures: paths, cycles, stars, edgeless and disconnected graphs, grids, random regular graphs, sparse/dense random graphs, and a complete graph.

- **78,012 single-bipartition comparisons** check every registry entry and Boolean color flag, including every target from 0 to n and every explicit seed on the exhaustive small graphs.
- **9,388 multipartition comparisons** check exact registries and balance for both kernels and seeded/unseeded drivers.
- **36 seed/range-search comparisons** cover three additional small graph fixtures. Sampled and automatic searches also check the next RNG draw after identical reseeding, verifying preserved RNG consumption.
- Graph preservation is checked, and every timed competitor is compared against the reference before timing.

Identical membership on an identical graph preserves the partition objectives; the experiment does not trade partition quality for speed on these cases. Finite tests do not prove equivalence on every possible input or prove BURY optimality. The existing package suite also remains at **115 passing tests and five expected broken tests**.

The prototypes deliberately differ in some observable API details: registries are `Dict{Int,Int}` rather than `Dict{Any,Any}`, dictionary iteration order is not promised, diagnostic output is not implemented by the single-bipartition kernels, and their input guards differ from legacy failures. They are experimental alternatives, not drop-in replacements with complete behavioral equivalence. Automatic seed search is handled by a separate adapter; omitting `seed` from a prototype kernel selects ordinary BURY.

## Recommended next steps

1. **Separate the reusable color-only kernel from registry construction.** Keep compatibility wrappers if the original dictionary type or logging must be retained, and benchmark those wrappers before adopting the reported speedups as package-level expectations.
2. **Consider the mask-free seeded scan as a small first change.** It removes substantial allocation traffic without a custom data structure. Retain the current ordinary selector unless a measured improvement justifies changing it.
3. **Add the heap implementation as a separately reviewable option for larger sparse graphs.** Keep deterministic tie order and colored-versus-processed invariants explicit in its tests. Avoid a universal heap default based only on this sample of workloads.
4. **Evaluate the multipartition driver independently.** Its map handling, color-only calls, and final-remainder assignment do not require changing the partition objective or seed policy.
5. **Optimize seed-search scoring separately if it dominates the application.** These experiments retain the current matching evaluator, candidate order, repeated sampled seeds, and tie rules. Caching repeated seed evaluations may help, but sampling without replacement would change the search policy and needs a separate quality comparison.

Further ideas remain unmeasured here: reusable workspace across seed trials, active-vertex masks that avoid induced-subgraph rebuilding, and replacing Floyd–Warshall in `antipode_seeds` with unweighted-graph traversal. The last change needs an explicit disconnected-graph policy. Do not assume that all target sizes in `bury_not_balanced_bipartition` can reuse a single greedy prefix: the remaining budget affects neighborhood-versus-singleton choices.

## Files and reproduction

The [complete study archive](../archive/README.md#bury-performance) preserves
the experiment README, commands, environment setup, prototypes, scripts,
validation and package-test logs, individual timing samples, sampling profiles,
and original-file verification script. Restore it before running the historical
commands. Two summary tables remain directly available:

- [Broad benchmark results](results/bury_performance/results.tsv).
- [Longer-trial results](results/bury_performance/slow_results.tsv) for the largest grid.

Run measurements serially after correctness checks. The isolated environment adds BenchmarkTools only within the experiment; it does not change the library's dependencies. The allocation and type-specialization improvements are consistent with [Julia's performance guidance](https://docs.julialang.org/en/v1/manual/performance-tips/). The measured results determine which prototype is useful for a workload.
