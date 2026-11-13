# Using QuantumHamlets.jl

The numerical workflow is to construct a graph, assign its vertices to partitions, and evaluate a cost. Visualization is optional. BURY uses the heap backend by default; the original qualified BURY function names still work and accept an explicit backend.

## Graphs and partition membership

```julia
using Graphs, QuantumHamlets

g = grid([6, 6])
p = partition_graph(g, 3)
partition_sizes(p)  # [12, 12, 12]
```

The supported graph input is a simple undirected `Graphs.Graph` with vertices numbered `1:nv(g)`. `partition_graph(g, k)` requires a nonempty graph, `1 ≤ k ≤ nv(g)`, and exact divisibility of `nv(g)` by `k`. It returns a `GraphPartition` rather than a dictionary or a tuple with algorithm metadata.

The algorithm is selected with a keyword:

| Call | Implementation |
| --- | --- |
| `partition_graph(g, k)` | BURY, equivalent to `algorithm=:bury, backend=:heap` |
| `partition_graph(g, k; algorithm=:bury_seeded)` | Seeded multipartition BURY using the existing minimum-weight seed rule |
| `partition_graph(g, k; algorithm=:random)` | Balanced random assignment using Julia's default RNG |
| `partition_graph(g, k; algorithm=:saran_vazirani)` | The existing Saran–Vazirani baseline, including its optimization solver |

`verbose=true` forwards BURY's progress output. The random method retains the existing default-RNG behavior and has no `rng` keyword. Use `Random.seed!(seed)` before a call when reproducing a run in the same Julia environment. The Saran–Vazirani baseline is experimental. It attempts to pack the original connected components before removing edges, so feasible zero-cut component assignments are retained. Its trace includes that initial attempt. Solver failure or exhaustion without a complete balanced assignment raises an error.

### Near-balanced bipartitions, including odd graphs

Use `bipartition_graph` to search ordinary BURY target sizes and return a
`GraphPartition` with two nonempty parts:

```julia
odd_graph = path_graph(11)
p = bipartition_graph(odd_graph; backend=:heap)
partition_sizes(p)  # [5, 6]
p = bipartition_graph(odd_graph; targets=4:7, objective=pairwise_cut_rank)
```

The default targets are `fld(nv(g),2):cld(nv(g),2)`: one candidate for an even
graph, both nearest-balanced candidates for an odd graph. A target is the size
of label 2. Every target gets an independent BURY run, in ascending order;
equal scores retain the last target. Complementary sizes can produce different
quality partitions, and their colored sets need not be nested.

The graph must have at least two vertices. Explicit `targets` must be a nonempty
contiguous ascending integer range within `1:nv(g)-1`, such as `4:7` or `5:5`.
Empty, reversed, floating-point, and non-unit-step ranges are rejected before
the objective runs. To require at least `m` vertices on each side, use
`m:nv(g)-m`. Wider intervals allow more imbalance and can favor smaller cuts;
this remains a heuristic search among BURY's candidates.

The default objective is `matching_cost`; `edge_cut_cost` and
`pairwise_cut_rank` are also supported directly. A custom scalar objective
uses the same graph-first contract and is invoked once for every target:

```julia
objective = (graph, part) -> edge_cut_cost(graph, part)
p = bipartition_graph(odd_graph; targets=4:7, objective, backend=:scan)
```

All three backends are supported. This helper uses ordinary BURY. Seeded and
custom partitioners remain available through the qualified range function
below. `partition_graph(g, k)` continues to require exact divisibility and
equal-sized parts, including when `k=2`.

### Simulated annealing and METIS baselines

The baseline functions return the same `GraphPartition` values as BURY:

```julia
using Graphs, QuantumHamlets
g = grid([6, 6])

p = min_cutrank(g; seed=3)        # random balanced start, then cut-rank SA
p = bury_sa(g; seed=3)           # ordinary BURY, then SA
p = bury_seeded_sa(g; seed=3)    # sqrt-sampled BURY, then SA
p = metis_partition(g)           # native METIS KWAY, 3% imbalance tolerance
p = metis_sa(g; seed=3)          # explicitly balance METIS, then SA

pairwise_cut_rank(g, p)
partition_sizes(p)
registry(p)                      # dictionary for legacy quantumLand calls
```

SA defaults to ten temperatures from 1.0 to 0.1 and retains the best rank seen,
including its initializer. Use `temperatures=range(1.0, 0.1; length=100)` to
change the schedule, or `keep_best=false` for the upstream final-state policy.
There is no time cap or automatic restart. Supply your own partition with
`min_cutrank(g; initial=p, seed=3)`: it preserves that partition's two nonempty
sizes. An empty schedule (`temperatures=[]`) returns the initializer alone.
All new methods use local seeded RNGs without resetting the caller's streams.
The seeded BURY initializer samples `ceil(sqrt(n))` seeds with replacement,
selects inner seeds by matching cost, and selects floor/ceiling targets by rank.

`metis_partition(g, k; balance=true)` explicitly repairs native METIS output to
floor/ceiling capacities using minimum edge-cut-increase vertex moves. This is
postprocessing, not a guarantee from the METIS tolerance. The extra capacities
go to the largest native parts, ties by label. Native output is the default for
`metis_partition`; repaired output is the default initializer for `metis_sa`.
Set `metis_sa(...; balance=false)` to preserve native sizes instead. Empty native
parts are rejected by SA. Both choices leave METIS's caller options unchanged.

The [bundled upstream source and MIT license](../src/min_cutrank/README.md) live
under `src/min_cutrank/`; PythonCall loads it directly with no separate Python
package installation. Metis and PythonCall are package dependencies; CondaPkg
manages a Python 3.12 interpreter. Call SA from the main Julia task unless you
manage PythonCall's GIL explicitly.

For the Toffoli grouped histogram, use the existing parent research environment
(`julia --project=..` from the package root), which supplies GraphIO and EzXML:

```julia
include("resources/toffoli_baselines.jl")
fig, data = plot_toffoli_baselines()  # reduced Toffoli + Barenco, all four sizes
display(fig)
save("toffoli_baselines.svg", fig)
# plot_toffoli_baselines(; collection="graphml_initial", repeats=20)
```

The single plotting function returns raw ranks, seeds, and actual sizes as well
as the figure. Bars show medians, whiskers show observed min–max, and METIS is
explicitly repaired so all plotted bipartitions have the same capacities.
Including the file does not run the comparison.

### Choose a BURY backend

The algorithm determines the partitioning procedure; the backend determines how BURY selects and updates its candidates:

| Backend | Implementation and intended use |
| --- | --- |
| `:heap` (default) | Indexed minimum heaps; generally strongest on larger sparse graphs |
| `:scan` | Minimum scans with incremental colored counts and fewer temporary arrays; useful for small or dense inputs and many small seeded workloads |
| `:legacy` | Original ordinary and seeded implementations, retained for compatibility and reference comparisons |

```julia
p = partition_graph(g, 3; backend=:heap)
p = partition_graph(g, 3; backend=:scan)
p = partition_graph(g, 3; backend=:legacy)
p = partition_graph(g, 3; algorithm=:bury_seeded, backend=:scan)
```

Backend selection preserves the greedy rule, smallest-vertex-ID ties, and seeded-region behavior. There is no automatic density threshold: graph size, target size, degree distribution, and seed can all change which backend is fastest. The [GraphML measurements](graphml_backend_findings.md) and [per-graph guide](results/graphml_bury/per_graph_guide.md) give examples.

`backend` applies to BURY. The `:random` and `:saran_vazirani` algorithms accept the default `backend=:heap` without using it; explicitly requesting `:scan` or `:legacy` with either raises `ArgumentError`. An unknown backend name also raises `ArgumentError`.

The original qualified BURY calls accept the same keyword:

```julia
const Q = QuantumHamlets
reg, colored = Q.bury(g, nv(g) ÷ 2; backend=:scan)
reg, colored = Q.bury_seeded(g, nv(g) ÷ 2, 1; backend=:heap)
reg = Q.bury_k(g, 3; backend=:legacy)
reg = Q.bury_k_seeded(g, 3; backend=:scan)
```

These functions retain their existing registry and tuple returns: `bury` and `bury_seeded` return `(Dict{Any,Any}, Vector{Bool})`, and the multipartition functions return `Dict{Any,Any}`. A colored vertex receives label 2 in a single bipartition. `verbose=true` retains the original trace for supported integer inputs on every backend.

Omitting the positional seed in `Q.bury_seeded(g, target; backend=...)` runs the existing sampled automatic seed search using that backend for every candidate and the final partition. Candidate order, repeated sampled seeds, matching scoring, and RNG consumption are preserved. `Q.best_seed` and `Q.best_sqrt_seed` also accept `backend`; switching it does not change their seed-selection policy. Multipartition calls carry the backend through every part. Seeded multipartition calls retain their deterministic minimum-weight seed rule.

`Q.bury_not_balanced_bipartition(g, min_v, max_v; backend=:scan)` uses scan for every target in the inclusive range. If `bury_v` is the built-in `Q.bury` or `Q.bury_seeded`, the keyword is forwarded. A custom `bury_v(graph, target)` callback keeps its two-argument calling convention and controls its own backend:

```julia
custom = (graph, target) -> Q.bury(graph, target; backend=:scan)
reg, colored = Q.bury_not_balanced_bipartition(g, 10, 18; bury_v=custom)
```

The legacy range function minimizes the first element of the tuple returned by
`cost_f(land, graph)`, retaining the last target on ties. It copies retained
registries and color vectors, so callbacks may reuse their output buffers or
mutate them after return without changing the winner. Custom callbacks retain
the existing `quantumLand` validation and are never inferred to be built-ins.
Ordinary BURY with the exact built-in matching, edge-cut, or rank function uses
the scalar cost kernels; heap/scan construct only the winning registry. The
legacy backend keeps its original algorithm. The [range performance
experiment archive](../archive/README.md#bipartition-range-performance) includes
scoring and returned representations in its before/after measurements.

For example, `bury_v=Q.bury_seeded, cost_f=Q.cut_rank_cost` invokes a sampled
seed search at each target that **always selects by matching cost**. The outer
search then selects the lowest rank among those matching-selected candidates.
An edge-cut outer objective behaves analogously. This preserves sampled order,
RNG consumption, first-best seed ties, and last-best target ties. Objective-aware
seed selection, common or unique seed pools, and reuse of the chosen inner
partition/score are deferred compatibility decisions. The seeded multipartition
driver separately uses a minimum-degree seed per residual graph.

Every backend follows the same input policy, with `ArgumentError` checks before sampling, callbacks, or solver execution:

- Single-partition targets must have an integer type and lie in `0:nv(g)`. Integral floating-point targets such as `2.0` are rejected, including by legacy.
- Explicit seeds must have an integer type and lie in `1:nv(g)`. All seeded calls and seed searches require a nonempty graph, even for target zero.
- Ordinary `bury(SimpleGraph(0), 0)` remains supported and returns an empty registry and color vector.
- Target ranges require integer bounds `0 ≤ min_v ≤ max_v ≤ nv(g)`. Ordinary range search supports the empty graph with `0:0`. Scores are initialized from the first candidate, so arbitrary score scales work; equal scores retain the last candidate.
- The low-level range still permits empty sides. Default matching search over `0:nv(g)` returns target `nv(g)` under last-target ties. `bipartition_graph` excludes these endpoints and requires two nonempty parts.
- Balanced calls (`partition_graph`, `bury_k`, `bury_k_seeded`, `k_partition_random`, and `k_partition_saran_vazirani`) require a nonempty graph, integer `k` in `1:nv(g)`, and exact divisibility of `nv(g)` by `k`.

Validation failures consume no RNG draws. Supported integer representations, valid-call RNG consumption, and the legacy algorithm's partition choices are preserved.

### Supply your own partition

A membership vector gives the partition label for each vertex. A dictionary maps vertex IDs to labels:

```julia
p = GraphPartition([1, 1, 2, 2])
p_from_dict = GraphPartition(Dict(1 => 1, 2 => 1, 3 => 2, 4 => 2))

membership(p)       # [1, 1, 2, 2]
registry(p)         # Dict(1 => 1, 2 => 1, 3 => 2, 4 => 2), in arbitrary order
partition_sizes(p)  # [2, 2]
length(p)           # 4
```

Labels must be integers in `1:nparts`, and dictionaries must contain exactly the vertex keys `1:n`. The partition count defaults to the largest label, or to one for empty membership. Specify `nparts` to represent empty parts:

```julia
p = GraphPartition([1, 1, 3, 3]; nparts=3)
partition_sizes(p)  # [2, 0, 2]
```

Unlike balanced partition generation, a supplied partition may be uneven or contain empty parts. Cost evaluation checks that its membership length matches the graph's vertex count.

Constructors copy their input. `membership(p)` and `registry(p)` also return copies, so changing either result does not change `p`. Treat the object's fields as internal and read-only; construct a new `GraphPartition` to change an assignment. This avoids synchronizing the separate villager and hamlet structures used by the original API.

## Numerical costs

All scalar cost calls take the graph first and partition second:

| Function | Result and meaning |
| --- | --- |
| `edge_cut_cost(g, p)` | Integer number of edges whose endpoints have different partition labels |
| `matching_cost(g, p)` | Integer sum of maximum matching sizes over unordered pairs of partitions |
| `pairwise_cut_rank(g, p)` | Integer sum of GF(2) biadjacency ranks over unordered pairs of partitions |
| `vertex_cover_costs(g, p)` | Named tuple `(ghz_free=..., ghz_n_minus_one=...)` calculated from the existing greedy vertex cover |

One declared part has matching cost zero, an empty matching witness, and no highlighted edges. Empty graphs and declared empty parts are supported by the cost APIs. Matching derives bipartite sides directly from part labels, without recursive graph traversal.

Scalar costs operate directly on `GraphPartition` membership, without constructing a legacy `quantumLand` or unused edge-color arrays. Matching groups cross edges by part pair, builds each native graph in bulk, and uses unweighted maximum-cardinality matching. Rank uses exact packed GF(2) elimination on biadjacency matrices, omitting isolated endpoints. The cover method keeps the existing greedy selection and counts neighboring part labels with reusable scratch storage. Legacy tuple-returning functions retain their formats.

The [cost performance archive](../archive/README.md#cost-performance) preserves bounded before/after timings and allocations against the implementation saved after the bug fixes. No graph partitioning or seed-selection policy changed.

The matching and rank values are exact for a fixed partition and their stated objectives. BURY is a heuristic for selecting the partition; it does not guarantee an optimal assignment. The cover routine uses a greedy cover, so its two estimates need not be minimum costs.

Matching and rank are different objectives. A complete bipartite graph with two vertices on each side illustrates the distinction:

```julia
g = complete_bipartite_graph(2, 2)
p = GraphPartition([1, 1, 2, 2])
edge_cut_cost(g, p)     # 4
matching_cost(g, p)     # 2
pairwise_cut_rank(g, p) # 1
```

For more than two partitions, `pairwise_cut_rank` is a sum over pairs, not one bipartite cut rank. The physical interpretation of each objective depends on the preparation model; see the [paper](https://arxiv.org/abs/2603.06387) for those assumptions.

### Inspect the matching edges

```julia
details = matching_details(g, p)
details.cost   # 2
details.edges  # matching edges with the original graph's vertex IDs
```

The result has named fields `cost` and `edges`. These edges form a valid maximum matching within each pair of partitions. A vertex may occur in matchings for several different pairs, so their union need not be one global matching. Different maximum matchings can have the same cost.

## Visualization

```julia
g = grid([6, 6])
p = partition_graph(g, 3)

figure, cost = plot_partition(g, p; highlight=:matching)
figure
```

`plot_partition` returns `(figure, cost)`. Nodes are colored by partition; the highlight controls the edge overlay:

| Highlight | Overlay | Returned cost |
| --- | --- | --- |
| `:matching` (default) | Pairwise matching edges | `matching_cost(g, p)` |
| `:cut_edges` | All cross-partition edges | `edge_cut_cost(g, p)` |
| `:none` | No edge overlay | `nothing` |

Use `highlight=:none` to display membership while reporting a rank or cover estimate separately. Save a figure with the plotting backend:

```julia
using CairoMakie
save("partition.png", figure)
```

The plotting dependencies still load with the package in this revision. Separating them into an optional extension is future work.

For the 6×6 grid example, all three BURY backends produce three groups of 12 vertices. The cut-edge score is 15 and the matching score is 9:

![The original unpartitioned 6×6 grid](../resources/figs/unpartioned6_6_grid.png)

![BURY partition with 15 cut edges highlighted](../resources/figs/bury_naive_6_6_grid.png)

![BURY partition with pairwise matching edges highlighted](../resources/figs/bury_vcg_6_6_grid.png)

These figures are the original example assets. They illustrate one graph and partition, not a general performance guarantee.

## Existing API and migration

The original API remains available through `QuantumHamlets`. Existing registries, tuple returns, and `quantumLand` constructors have not been replaced. Its BURY calls now default to heap; add `backend=:legacy` to request the original implementation. A calculation with the original API can still be written as:

```julia
using Graphs, QuantumHamlets
const Q = QuantumHamlets

g = grid([6, 6])
reg = Q.bury_k(g, 3)
land = Q.quantumLand(reg, 3, nv(g))
cost, edge_colors = Q.matching_cost_of_partition(land, g)
figure, cost = Q.visualize_graph_on_land(
    land, g; method=Q.matching_cost_of_partition
)
```

Use the conversion methods when integrating new code with existing scripts:

```julia
p = GraphPartition(land)
new_land = Q.quantumLand(p)
reg = registry(p)
```

| Original expression | Simpler interface |
| --- | --- |
| `Q.bury_k(g, k)` | `partition_graph(g, k)` |
| Ordinary BURY range with two nonempty parts | `bipartition_graph(g; targets=lo:hi, objective=matching_cost)` |
| `Q.bury_k_seeded(g, k)` | `partition_graph(g, k; algorithm=:bury_seeded)` |
| `Q.k_partition_random(g, k)` | `partition_graph(g, k; algorithm=:random)` |
| `Q.quantumLand(reg, k, nv(g))` | `GraphPartition(reg; nparts=k)` |
| `first(Q.naive_cost_of_partition(land, g))` | `edge_cut_cost(g, p)` |
| `first(Q.matching_cost_of_partition(land, g))` | `matching_cost(g, p)` |
| `first(Q.cut_rank_cost(land, g))` | `pairwise_cut_rank(g, p)` |
| `Q.vertex_cover_cost_of_partition(land, g)` | `vertex_cover_costs(g, p)` with named fields |
| `Q.visualize_graph_on_land(land, g; method=Q.matching_cost_of_partition)` | `plot_partition(g, p; highlight=:matching)` |

The original cost functions keep their return formats: edge-cut and matching return `(cost, edge_colors)`, rank returns `(cost, nothing)`, and cover returns two numerical costs. In particular, the cover tuple is not a score plus plotting colors. The new plotting interface uses explicit overlays instead of accepting an arbitrary cost function.

The legacy `quantumLand(reg, k, n)` constructor copies its dictionary, requires positive integer `k`, nonnegative integer `n`, exactly the integer keys `1:n`, and integer labels in `1:k`. Declared empty parts and an empty registry with positive `k` are supported. Conversions from `GraphPartition` preserve them. Treat legacy fields and nested membership objects as read-only: directly mutating one representation can still make it disagree with the others. All legacy display methods honor the supplied IO stream.

### Compare an external partitioner

Metis is an example dependency, not a QuantumHamlets dependency. Install it in an environment that also develops this checkout, then convert its labels directly:

```julia
using Graphs, QuantumHamlets, Metis

g = grid([6, 6])
bury = partition_graph(g, 3)
metis = GraphPartition(Metis.partition(g, 3); nparts=3)

matching_cost(g, bury)
matching_cost(g, metis)
```

Use the same objective for both algorithms, and inspect `partition_sizes` when comparing their balance. External partitioners can vary with their version or random choices. The original METIS example figure is retained here:

![Original METIS grid example with matching edges highlighted](../resources/figs/metis_vcg_6_6_grid.png)

## Tests and source organization

Run the suite from the package directory:

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
```

The tests cover deterministic fixtures, exhaustive small graphs, valid partition sizes, numerical costs, input graph preservation, conversions, and equivalence against a frozen original BURY implementation. Backend checks also cover automatic seed selection, RNG consumption, custom callbacks, return types, and verbose traces. The five originally deferred failures are now ordinary passing regressions. Additional tests cover independent matching cardinalities and witnesses, one-part plots, validation before RNG/callbacks, disconnected zero-cut packing, constructor ownership, captured display output, and scalar research scores. Large matching checks run separately in disposable subprocesses with external time limits (`python3 test/run_matching_stress.py`). The full suite after the target-range improvements passed 459,086 assertions with no broken tests, including 1,877 new range assertions; its validation log is preserved in the [range-performance archive](../archive/README.md#bipartition-range-performance).

The source is organized within one module, with seven top-level Julia files
and three BURY backend files:

| File | Responsibility |
| --- | --- |
| `src/QuantumHamlets.jl` | Imports, exported API, and includes |
| `src/quantum_land.jl` | Original villagers, hamlets, and `quantumLand`, their display methods, and legacy edge removal |
| `src/partition.jl` | `GraphPartition`, conversions, accessors, and graph/partition size validation |
| `src/costs.jl` | Shared cost kernels, graph-first scalar costs and matching details, and legacy tuple adapters |
| `src/bury.jl` | BURY kernels' dispatch, seed searches, vertex-map helpers, `partition_graph`, and `bipartition_graph` |
| `src/bury/heap.jl` | Indexed-heap color kernel |
| `src/bury/scan.jl` | Optimized scan color kernel |
| `src/bury/legacy.jl` | Original ordinary and seeded implementations |
| `src/baselines.jl` | Random, Saran–Vazirani, METIS, and min-cutrank/SA baselines and initializers |
| `src/visualization.jl` | Graph plotting, legacy plotting adapters, and `plot_partition` |
| `src/min_cutrank/` | Vendored upstream Python implementation, adapter, provenance, and license |

The November 2026 source consolidation placed each public function beside its
implementation. Existing function names, arguments, return types, exports,
and algorithm policies are preserved; normal `using QuantumHamlets` calls need
no migration. The [streamlining record](streamlining_plan.md#3-source-consolidation)
maps the former files to their current locations. `GraphPartition` is loaded
before cost and plotting methods that use it in their signatures.

After consolidation, the full package suite passed 459,283 assertions with
`Pkg.test(; allow_reresolve=false)`. An independent parsed-source comparison
also confirmed that the moved definitions and docstrings were unchanged.

Maintained research examples and original graph datasets remain under `resources/`. Completed studies, earlier plotting workflows, and historical design reviews are preserved through the [archive index](../archive/README.md). Reports and selected results remain directly readable under `docs/`.

## Known limitations

The repairs preserve valid-input BURY partition choices and scoring/seed policies across all backends. Matching witnesses may choose different optimal edges; their cardinality and validity define the contract.

- BURY and Saran–Vazirani remain heuristics; neither guarantees globally optimal balanced partitions. Saran–Vazirani depends on an optimization solver and reports failures explicitly.
- The greedy vertex-cover estimates need not be minimum costs. The legacy cover callback returns two numeric estimates and remains incompatible with the old plot callback's `(score, edge_colors)` convention; use the explicit `plot_partition` overlays.
- Copying constructor inputs prevents input aliasing, but legacy nested objects remain mutable. Construct a new partition rather than mutating their fields independently.
- Matching no longer has a recursive traversal limit, but large solver workloads can still be expensive. Plotting dependencies still load with the package.
- The grid and sandbox partitioner comparisons now use scalar pairwise GF(2) cut rank for every competitor and label that objective. Their defaults are full research sweeps, not quickstarts; only tiny equivalents were tested. Historical results and audit logs were not regenerated.

The resolved audit and validation record is preserved in the [historical documentation archive](../archive/README.md#historical-docs); the limitations above describe the maintained API.
