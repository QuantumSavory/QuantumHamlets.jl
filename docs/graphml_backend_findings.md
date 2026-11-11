# BURY backend performance on the GraphML resources

Measured on 11 November 2026. These results support **heap as the default,
scan as a useful explicit alternative, and legacy as a reference backend**.
Heap had the lowest median for ordinary BURY on 39 of 45 graphs. Explicitly
seeded BURY was more evenly divided: heap won 24 graphs and scan won 21.
Legacy had no lowest median in the completed comparisons.

The distinction is broader than sparse versus dense. The reduced graphs often
favor scan when seeded, while their ordinary BURY calls usually favor heap.
Changing the requested target can reverse the winner on the same graph.
The proposed `backend=:heap/:scan/:legacy` interface therefore fits the evidence
better than a `dense=true` switch. No production code has been changed.

## Graph properties and initial expectations

All 45 files in the three `resources/graphml*` directories were inspected.
Every graph is connected, simple, and undirected, with no self-loops, duplicate
edges, or isolated vertices. Density is `2m / (n(n-1))`.

| Directory | Files | Vertices | Mean degree | Density | Expectation from properties |
| --- | ---: | ---: | ---: | ---: | --- |
| `graphml_initial` | 20 | 67–652 | 2.39–2.67 | 0.394–3.62% | Low degree should favor heap as size grows. |
| `graphml_reduced` | 20 | 25–294 | 3.49–8.07 | 2.75–16.67% | Small size and more neighbor updates make scan competitive. |
| `graphml_evolution` | 5 | 25–1,772 | 3.36–13.66 | 0.771–14.0% | Small cases may favor scan; larger cases should favor heap, with hubs adding update work. |

Reduction removes vertices but does not always remove edges. For example,
`gf2^6_mult` changes from 652 vertices / 855 edges to 294 vertices / 1,186 edges.
The reduced version has maximum degree 88. Evolution `k5` has 1,772 vertices,
12,102 edges, and maximum degree 222. These graphs can have high-degree hubs
despite low overall density.

The [full inventory](../experiments/graphml_bury/inventory/inventory.tsv) gives
each file's size, edge count, degree distribution, density, components, and hash.
The [per-graph guide](../experiments/graphml_bury/per_graph_guide.md) combines
those properties with the observed backend choices.

## Half-size partition results

The target is `floor(n/2)`. Explicitly seeded calls use the first declared
GraphML node, mapped to Julia vertex 1. All three backends return identical
partition contents on every completed workload. Automatic searches use one
fixed RNG stream, reset before each evaluation. These measurements do not
characterize variation across seed choices.

| Directory | Ordinary: heap / scan wins | Explicit seeded: heap / scan wins |
| --- | ---: | ---: |
| Initial | 20 / 0 | 17 / 3 |
| Reduced | 15 / 5 | 4 / 16 |
| Evolution | 4 / 1 | 3 / 2 |
| **Total** | **39 / 6** | **24 / 21** |

These count the lowest sample median, including close results. Of the 90
ordinary/seeded comparisons, ten have heap and scan within 10% of each other.
The strongest practical scan case is small, explicitly seeded input, rather
than a general rule that every reduced graph should use scan.

Selected medians below are **microseconds**. “Seeded” excludes automatic seed
selection; allocation and seed-search results follow separately.

| Graph | Vertices | Variant | Legacy µs | Scan µs | Heap µs |
| --- | ---: | --- | ---: | ---: | ---: |
| Initial `barenco_tof_3` | 83 | Ordinary | 8.33 | 5.04 | **2.96** |
| Initial `barenco_tof_3` | 83 | Seeded | 13.8 | **3.17** | 4.25 |
| Reduced `barenco_tof_3` | 25 | Ordinary | 1.25 | **0.583** | 0.709 |
| Reduced `barenco_tof_3` | 25 | Seeded | 2.04 | **0.542** | 1.08 |
| Initial `barenco_tof_10` | 571 | Ordinary | 208 | 202 | **30.6** |
| Initial `barenco_tof_10` | 571 | Seeded | 572 | 101 | **42.7** |
| Reduced `barenco_tof_10` | 123 | Seeded | 28.1 | **7.29** | 7.63 |
| Initial `gf2^6_mult` | 652 | Ordinary | 295 | 271 | **46.3** |
| Reduced `gf2^6_mult` | 294 | Ordinary | 65.6 | 53.4 | **20.9** |
| Reduced `gf2^6_mult` | 294 | Seeded | 189 | 37.1 | **28.6** |
| Evolution `k2` | 159 | Seeded | 36.3 | **9.58** | 11.4 |
| Evolution `k5` | 1,772 | Ordinary | 1,471 | 1,512 | **194** |
| Evolution `k5` | 1,772 | Seeded | 5,943 | 1,018 | **322** |

Thus heap is about **6.8× faster than legacy for ordinary initial
`barenco_tof_10`**, and **13.4× faster for its explicitly seeded call**. On
evolution `k5`, the corresponding gains are **7.6× and 18.5×**. For the tiny
reduced `barenco_tof_3`, scan is approximately twice as fast as heap when seeded,
although the absolute difference is only about half a microsecond.

![Scan versus heap timing across the 45 GraphML files](../experiments/graphml_bury/backend_comparison.png)

Values above 1 favor heap; values below 1 favor scan. The gray band marks
differences within 10%, not a statistical confidence interval. The ordinary
and seeded panels have different crossover behavior.

The seeded `k5` legacy call allocates **26.2 MB per call**, versus **0.0879 MB**
for scan and **0.144 MB** for heap. Scan can allocate less yet run slower,
because it still scans all vertices each round. These are cumulative allocated
bytes, not peak memory usage.

## Target size can reverse the choice

For reduced `gf2^6_mult`, the graph properties and explicit seed are unchanged:

| Target | Scan µs | Heap µs | Lower median |
| ---: | ---: | ---: | --- |
| 73, approximately one quarter | 11.96 | 18.13 | Scan |
| 147, one half | 37.13 | 28.63 | Heap |
| 220, approximately three quarters | 51.71 | 40.33 | Heap |

This is a direct counterexample to selecting solely by density. Target size
changes the sequence and number of greedy selection rounds. The experiment
completed 40 additional ordinary/seeded quarter- and three-quarter-target
comparisons on representative graphs; all retained exact output agreement.

## Automatic seed search and the runtime stop

Automatic selection runs multiple candidate partitions and evaluates matching
costs. Its timing should not be inferred from the explicitly seeded kernel.
Here the heap and scan medians are within **5% on 41 of 44 completed graphs**
and within **10% on 42 of 44**. Neither has a decisive general advantage for
these complete searches.

Selected automatic-search medians are **milliseconds**:

| Graph | Legacy ms | Scan ms | Heap ms |
| --- | ---: | ---: | ---: |
| Initial `barenco_tof_3` | 1.20 | 1.08 | 1.10 |
| Reduced `barenco_tof_3` | 0.230 | 0.215 | 0.215 |
| Initial `barenco_tof_10` | 59.8 | 44.0 | 43.2 |
| Initial `gf2^6_mult` | 81.2 | 59.4 | 57.3 |
| Reduced `gf2^6_mult` | 36.9 | 33.7 | 33.8 |
| Evolution `k4` | ≈921 | ≈844 | ≈820 |

The `k4` automatic row has only **three samples per backend** because one call
exceeds the short trial budget. Its small heap/scan difference is only an
approximation, not a confident ranking. Search and scoring overhead limit the
benefit of accelerating the individual partition calls.

The Barenco files were tested in parameter order, and the smaller calls
finished comfortably before larger cases were attempted. No long-running
Barenco call was reproduced for these targets and seed policies.

The stop rule did trigger for **automatic seed selection on evolution `k5`**:
a single warmed **legacy preflight** took **4.15 seconds** and allocated
**1.09 GB**. That is not a benchmark median. Scan and heap automatic searches
were not attempted afterward, and the four remaining quarter/three-quarter
core workloads for that graph were skipped. Its ordinary and explicitly
seeded half-target calls had already passed. The skipped core workloads have
not themselves been shown to be slow.

## Recommendation for the backend API

Use `:heap` by default, keep `:scan` directly selectable, and retain `:legacy`
as a reference. For these resources:

- **Ordinary BURY:** heap is a strong default, including all 20 initial graphs
  and most reduced graphs.
- **Explicitly seeded BURY:** prefer scan for many small reduced graphs;
  use the per-graph guide when the workload is known. Heap remains the stronger
  option on the larger initial and evolution graphs.
- **Automatic seed search:** expect much smaller gains from switching heap
  and scan. Treat matching/scoring optimization as a separate task.

Do not advertise scan only as a dense-graph option, or infer an automatic
density threshold from these fixtures. Explicit backend selection lets users
account for graph size, degree distribution, target, and seeded behavior.

## Measurement and validation

Measurements used Julia 1.12.6, Graphs 1.13.4, BenchmarkTools 1.6.3, and one
Julia thread on an Apple M1 Pro. Each backend was warmed, and three timing
rounds rotated their order. Medians pool raw samples across the rounds. Timing
includes registry/color construction and, for automatic calls, unchanged seed
search and scoring. Loading, compilation, and graph-property computation are
excluded. A full GC runs once per workload group outside timing; automatic GC
remains enabled throughout.

There are **174 completed workload groups**, **522 backend result rows**, and
**121,851 timing samples**. The 348 nontrivial optimized-versus-legacy output
comparisons passed. Every checked call preserves the graph and reaches the
target; the 44 complete automatic comparisons also preserve the next RNG draw
for each optimized backend. A separate audit recomputed medians and sample
counts from raw measurements. All 210 pre-existing repository files remain
byte-for-byte unchanged.

The files contain node type, phase, layout, and edge-type attributes. BURY
here operates on their simple unweighted topology. Node declaration order is
preserved when mapping noncontiguous XML IDs to `1:n`; circuit attributes do
not affect these timings. This study covers single bipartitions and automatic
seed search, not complete multipartition workflows or exhaustive seed search.

The prototypes return typed dictionaries while legacy returns `Dict{Any,Any}`,
so their construction costs are included in the measured differences.
Partitions agree on tested inputs; logging, invalid-input behavior, dictionary
types, and complete API compatibility still need reconciliation before the
prototypes become production backends.

The [experiment README](../experiments/graphml_bury/README.md) contains commands,
watchdog limits, and file descriptions. See [all measured comparisons](../experiments/graphml_bury/comparison.tsv),
[raw results](../experiments/graphml_bury/results.tsv),
[timing samples](../experiments/graphml_bury/samples.tsv), and
[completion/skip records](../experiments/graphml_bury/status.tsv).
