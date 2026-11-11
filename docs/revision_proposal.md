# QuantumHamlets.jl revision proposal

Draft reviewed on 11 November 2026 against the current working tree, including local edits, based on commit `9b0ad4f`. Source locations below refer to that snapshot.

The most useful first step is a small regression suite around partition validity and cost evaluation, followed by fixes to the confirmed defects below. Then separate the public numerical API from plotting, split the source into focused files, and optimize the measured bottlenecks. This order protects the existing scientific behavior while making the library easier to learn and extend.

The current package has one 682-line source module, eleven direct dependencies including `Random`, and no tracked test suite or CI workflow. Its size does not require a large framework: one module with several included files should be sufficient. The core BURY and cost calculations passed the small valid-input checks described below; the clearest problems concern boundary conditions, inconsistent interfaces, and research-script maintenance.

## Suggested priorities

| Order | Change | Benefit | Relative effort |
| --- | --- | --- | --- |
| 1 | Add regressions and fix invalid targets, missing results, and the one-hamlet matching case | Prevent hangs and invalid partitions | Small to medium |
| 2 | Repair experiment cost handling and define numerical return contracts | Make comparisons and downstream analysis reliable | Small to medium |
| 3 | Validate membership, fix display, and document supported inputs | Give users predictable behavior and useful errors | Medium |
| 4 | Split files while preserving existing calls | Make algorithms and tests easier to navigate | Small |
| 5 | Introduce a simpler partition API and optional plotting | Reduce setup requirements and conceptual overhead | Medium |
| 6 | Benchmark and optimize matching, pair extraction, rank, and BURY | Improve throughput without changing objectives accidentally | Medium to large |

“Small” means a localized change; “medium” means coordinated API, documentation, and test work. These are scope estimates, not timing estimates.

## Confirmed defects and regression targets

The examples use `using QuantumHamlets, Graphs` and `const Q = QuantumHamlets`.

### Invalid partition targets can hang or produce misleading results

In [the BURY implementation](../src/QuantumHamlets.jl), lines 296–345, `bury(path_graph(4), 5)` does not terminate. A subprocess probe was stopped after 15 seconds. Once all four vertices are colored, every weight is `Inf`, but the loop still requires five colored vertices. The seeded implementation at lines 473–532 has the same missing upper-bound guard.

Other observed input failures are less severe but should use the same validation policy: `bury_k(path_graph(4), -2)` returns an empty dictionary, `k=0` raises `DivideError`, a negative BURY target returns an all-uncolored result, and a fractional target is accepted. `k_partition_random(path_graph(5), 2)` raises `InexactError` through `Int(n/k)` at line 245.

Require integer targets with `0 ≤ r0 ≤ nv(g)`, valid seeds, and positive integer `k`. For balanced nonempty partitions, require `1 ≤ k ≤ nv(g)` and `nv(g) % k == 0`. Keep exact divisibility initially and explain it in an `ArgumentError`; support unequal capacities later as an explicit feature. Decide and test the empty-graph policy separately. Validate the public entry points before indexing or looping.

### Searching valid candidates can return no partition

`bury_not_balanced_bipartition` initializes `lowest_cost = n` at line 358, but its configurable cost can exceed `n`:

```julia
Q.bury_not_balanced_bipartition(
    complete_graph(8), 4, 4; cost_f=Q.naive_cost_of_partition
)
# Current result: (nothing, nothing)
# The valid 4+4 candidate has cut-edge cost 16.
```

Initialize the best result from the first valid candidate or use an appropriate unbounded sentinel. Reject inverted or out-of-range search intervals. Test a valid cost greater than `nv(g)` so this does not regress when new objectives are added.

### One hamlet has no matching cost but currently throws

At lines 181 and 216, `cost_all = []` followed by `sum(cost_all)` fails when `k=1`: the empty vector has element type `Any`, so Julia cannot obtain its additive identity. A three-vertex path entirely within one hamlet raises `MethodError: zero(::Type{Any})`; it should return zero with all edges unhighlighted. An integer accumulator also removes the need to retain per-pair costs when they are not requested.

### Saran–Vazirani can return an empty registry for a feasible problem

`k_partition_saran_vazirani(SimpleGraph(4), 2)` returns an empty registry and trace. All component assignment occurs inside the edge-cut loop at lines 588–632, which is never entered for an edgeless graph.

Attempt to pack the original connected components before removing edges. If the routine cannot produce a partition, return an explicit failure or raise a documented exception instead of returning an incomplete registry. Add cases for edgeless graphs, already disconnected graphs, `k=1`, indivisible sizes, and solver failure. Keep solver status distinct from a valid partition result.

### Membership can become internally inconsistent

The constructor at lines 55–69 stores an input `Dict{Int,Int}` by reference while independently creating hamlet/villager objects:

```julia
registry = Dict(1 => 1, 2 => 2)
land = Q.quantumLand(registry, 2, 2)
registry[1] = 2
# land.registry[1] is now 2, but villager 1 remains stored in hamlet 1.
```

Copy constructor inputs as an immediate repair. Longer term, use one authoritative membership representation and derive groups and counts from it. Merely copying the input does not prevent callers from mutating the public nested objects inconsistently. The current constructor also accepts negative dimensions and extra registry entries. Validate exact vertex coverage, label bounds, counts, and agreement with the graph passed to a cost function.

### Display ignores the requested output stream

All three `Base.show` methods at lines 668–680 call `print` or `println` without `io`. Consequently, `sprint(show, Q.villager(1, 1))` returns `""` while writing to stdout. Pass `io` through every nested display call. Test captured output, and consider a compact ordinary `show` plus a more detailed `MIME"text/plain"` display for interactive use.

### Comparison scripts no longer consume compatible scores

In [grid_plots.jl](../resources/grid_plots.jl), lines 46–47 replace the BURY matching score with the complete `(rank, nothing)` tuple from `cut_rank_cost`. METIS still receives a scalar matching score at lines 59–60. This creates incompatible plot data and selects different objectives for the two methods.

In [sandbox.jl](../resources/sandbox.jl), lines 165, 172, 184, and 191 also collect cut-rank tuples where subsequent `mean` and `std` calls need numerical scores. These are source-confirmed integration defects; the large experiment runs were not executed during this review.

Give comparison functions one explicit scalar objective argument, apply it to every algorithm, and put the objective name in the plot label and saved data. For an immediate repair, destructure the cut-rank return consistently. Test a tiny comparison with two algorithms and one graph before launching a sweep.

## Make the public interface smaller and more explicit

The intended user workflow can be expressed as **graph → partition membership → numerical cost → optional visualization**. Ordinary graph users should not need to construct villagers and hamlets to evaluate a partition.

An illustrative future API, not currently implemented, is:

```julia
using Graphs, QuantumHamlets
g = grid([6, 6])
p = partition_graph(g, 3; algorithm=:bury)
edge_cut_cost(g, p)
matching_cost(g, p)
pairwise_cut_rank(g, p)
```

Use a small validated `Partition` representation with membership stored as a `Vector{Int}` for vertices `1:n`, together with the declared partition count if empty parts are supported. Derive `n` from the vector length. Accept dictionaries at the boundary and provide conversion back to the existing registry format. Specify ownership and mutation rules; an immutable Julia struct containing a vector does not make the vector immutable.

Keep `quantumLand` and its quantum terminology available for existing users. Introduce new calls through wrappers rather than simultaneously changing every name, return type, and algorithm. A first structural change can preserve all existing behavior; the new API can follow with migration examples and compatibility tests.

The module currently has no `export` or `public` declarations. Identify a supported surface explicitly, using `public` if qualified names are preferred, and label seed-search helpers and experimental baselines accordingly. Generalize `Graphs.Graph` arguments only where the implementation supports the broader Graphs interface. State the assumptions about simple undirected graphs, vertex numbering, and weights before promising additional graph types.

### Separate scores from plotting data

Current return values have three different meanings:

| Function | Current result |
| --- | --- |
| `naive_cost_of_partition` | `(cost, edge_colors)` |
| `matching_cost_of_partition` | `(cost, edge_colors)` |
| `cut_rank_cost` | `(cost, nothing)`, despite its integer-return docstring |
| `vertex_cover_cost_of_partition` | `(ghz_free_cost, ghz_n_minus_one_cost)` |

`visualize_graph_on_land` at lines 659–664 treats every second value as an edge color. Passing the vertex-cover function actually constructs a figure in the smoke probe, but feeds its second numerical cost into the color parameter. That is a semantic mismatch, even when it does not raise an error.

Make ordinary cost calls return integers. Provide a details function or structured result when witnesses are needed, such as `(cost=..., matching_edges=...)`; plotting should convert those edges to colors. Give the two GHZ estimates named fields. Keep the existing tuple-returning methods as compatibility wrappers. A shared internal evaluator can compute the score and witness together without rerunning matching.

For rank, a plot may reasonably show membership without a matching overlay. Make this an explicit supported behavior. Reject incompatible overlays with a useful message instead of interpreting arbitrary values as colors.

### State precisely what each metric measures

`matching_cost_of_partition` sums maximum matching sizes over unordered pairs of hamlets. A vertex may participate in matchings for different pairs; the combined highlighted edges need not form one global matching.

`cut_rank_cost` sums the GF(2) ranks of the pairwise biadjacency matrices. Its present division by two is mathematically correct: each bipartite adjacency matrix has block form `[0 B; Bᵀ 0]`, whose rank is twice `rank(B)`. For more than two hamlets, call this a sum of pairwise cut ranks. Avoid conflating it with one bipartite cut rank or an unrestricted optimum preparation cost.

The vertex-cover routine uses `DegreeVertexCover`, which Graphs documents as a greedy heuristic. Its result is an estimate based on a particular cover, not necessarily a minimum. On the concrete 13-vertex fixture in the validation notes below, it returns 5 although a cover of size 4 exists. Document this distinction; an exact option can be added where scientifically useful. [Graphs vertex-cover documentation](https://juliagraphs.org/Graphs.jl/stable/algorithms/vertexcover/).

## Proposed organization

Keep a single Julia module and use `include` for these responsibilities:

```text
src/
  QuantumHamlets.jl          imports, public declarations, includes
  partition.jl              membership, validation, conversions
  bury.jl                   BURY variants and seed strategies
  baselines.jl              random and Saran–Vazirani partitioning
  costs.jl                  numerical objectives
  bipartite.jl              pair extraction and matching adapter
  display.jl                text display
  compatibility.jl          existing names and return formats
ext/
  QuantumHamletsMakieExt.jl  optional visualization
test/
  runtests.jl
  partition.jl
  bury.jl
  costs.jl
  regressions.jl
  plotting.jl
examples/
  Project.toml
  grid_partition.jl
  compare_metis.jl
benchmark/
  Project.toml
  benchmarks.jl
experiments/
  Project.toml
  ...research scripts...
docs/
  ...user guide and API documentation...
```

Start by moving existing definitions without changing behavior. Share cross-partition extraction between matching and rank only after tests establish the current contracts. Move research scripts and their paths together in a later change; preserve publication assets and intentional datasets.

Plotting is the clearest optional dependency boundary. Move graph visualization into a Makie/GraphMakie extension and let users load their preferred backend, with CairoMakie used in examples. Review the unused direct `Colors` import. Consider a separate solver extension for Saran–Vazirani later; keep IGraphs and Nemo required while the supported numerical API depends on them. Julia package extensions support this separation without requiring every plotting dependency during ordinary package loading. [Pkg extension guide](https://pkgdocs.julialang.org/v1.11/creating-packages/#Conditional-loading-of-code-in-packages-(Extensions)).

## Testing plan

Add `Test` through `[extras]` and `[targets]`, with `test/runtests.jl` as the entry point. This works with the declared Julia 1.11 minimum. Use separate declared environments for examples and benchmarks; do not require a newer workspace mechanism implicitly. [Pkg test dependency guide](https://pkgdocs.julialang.org/v1.11/creating-packages/#Test-specific-dependencies).

| Area | Tests that provide useful protection |
| --- | --- |
| Membership | Missing/extra vertices, invalid labels and dimensions, graph-size mismatch, copied inputs, and explicit mutation behavior |
| Partition algorithms | Exact coverage, label range, promised sizes, disconnected/edgeless graphs, `k=1`, `k=n`, valid target endpoints, invalid targets, and unchanged input graph |
| BURY | Exhaustive small-graph cardinality checks for ordinary and seeded variants; fixed fixtures for intentional tie behavior; search with a score greater than `n` |
| Matching | Independent exhaustive matching oracle on tiny graphs; unmatched vertices; one hamlet; isolated vertices; highlight count equals cost; matching validity within each hamlet pair |
| Rank | Independent Boolean XOR elimination; known GF(2) examples; local-edge invariance; pairwise aggregation |
| Cover estimates | Every cross edge is covered; the documented greedy counterexample; named GHZ cost meanings |
| Public behavior | Capture `show`; quiet defaults; explicit RNG reproducibility; old and new API equivalence |
| Integration | Numerical README example in a clean environment; a small saved headless plot; a tiny experiment using one common objective |

For a fixed partition, test `0 ≤ pairwise_rank ≤ matching_sum ≤ cut_edges`. These three exact objectives should be invariant under consistent vertex relabeling, hamlet relabeling, and adding intra-hamlet edges. A greedy cover estimate can change under vertex relabeling because tie choices can change; test cover validity rather than requiring identical estimates. Similarly, do not require a partition heuristic to choose the same partition after vertex relabeling. Alternative maximum matchings can also have different highlighted edges while sharing the same cost.

Start CI with the declared minimum Julia version and current stable, including Linux and macOS coverage. Run `Pkg.test()` and the small example smoke checks. Add a fresh dependency-resolution job to complement the existing pinned manifest, without modifying the research environment. Keep large data sweeps and strict wall-clock performance thresholds out of routine CI.

## Performance opportunities to benchmark

These are source-grounded opportunities, not measured speedup claims. Measure time and allocations separately for partitioning, graph extraction, numerical scoring, and plotting.

| Priority | Current work and location in `src/QuantumHamlets.jl` | Candidate change and validation |
| --- | --- | --- |
| High | Matching creates all-one weights and calls the weighted backend, lines 194–198 | Use the unweighted backend through the installed wrapper's `IGNull()` support. Verify cardinality and witness validity against the current path. Different optimal witnesses are acceptable. |
| High | Matching/rank scan vertices for every hamlet pair; matching rescans all original edges for every pair, lines 120–140 and 154–214 | Group vertices once, scan cross edges once into pair buckets, retain original edge IDs, skip empty pairs, and process pairs without retaining all temporary graphs. This targets the current O(k²n) vertex scans and O(k²m) matching color scans, excluding solver work. |
| High | Matching reconstructs bipartite sides using recursive DFS and repeated `sum(visited)`, lines 170–192 | Assign sides directly from the two hamlet labels. This removes unnecessary traversal, repeated scans for disconnected components, and deep-recursion risk. |
| High | Cost evaluation always builds visualization data, including inside seed searches, lines 78, 180–214, 394–396, 417–419 | Provide a score-only path. Consume the backend's returned cardinality and request original-edge mapping only when details are needed. |
| Medium | Rank builds full square adjacency matrices, lines 136–139 | Rank the rectangular cross-part matrix B directly over GF(2). With equal side size s, this represents s² entries instead of 4s² in the full block matrix. Compare results against the existing implementation before considering bit-packed elimination. |
| Medium | Registries and several temporary collections use `Dict()` and `[]`, throughout | Use `Dict{Int,Int}`, typed edge/group collections, or dense membership vectors. Remove unused `maps` in rank and `matchings` in matching. Check inference and allocations at function boundaries. |
| Medium | BURY repeatedly scans colors and weights and allocates neighborhood temporaries, lines 315–338 and 497–524 | Track colored count incrementally; update weights with explicit neighbor loops; avoid constructing a closed-neighborhood vector. Profile before adopting a priority queue, and preserve tie order and the distinct meanings of colored versus processed vertices. |
| Medium | `bury_k` and its seeded variant rebuild remainder graphs after every part, lines 260–284 and 445–470 | Assign the final remainder directly; then investigate reusable maps or active masks. Keep original-vertex mapping tests. |
| Later | Seed searches rerun BURY and matching; antipodes use Floyd–Warshall, lines 387–442 | Reuse workspace; sample without replacement if that is the intended policy; compare sparse-graph BFS alternatives for antipodes. Reducing the search budget changes heuristic quality and needs separate evaluation. |
| Later | Saran–Vazirani allocates a dense capacity matrix and solves many cuts/packing problems, lines 581–614 | Profile solver and preprocessing costs separately. Investigate supported sparse capacities and reusable solver state; make the retained trace optional. |

The unweighted matching candidate has unusually direct support: the installed IGraphs wrapper accepts `IGNull()` and returned `(2, 2.0)` on `K₂,₂`, matching the all-one weighted result. The underlying igraph 0.10.15 source dispatches these inputs to different implementations and documents O(√V E) for unweighted versus O(VE) for weighted matching. Benchmark end-to-end package workloads before claiming that complexity difference translates to a particular speedup. [igraph matching implementation](https://raw.githubusercontent.com/igraph/igraph/0.10.15/src/misc/matching.c).

Typed containers, reduced allocations, and separating preparation from specialized numerical kernels follow Julia's performance guidance. Profile first; adding concrete annotations to every argument or replacing a proven numerical backend is not a substitute for measurement. [Julia performance tips](https://docs.julialang.org/en/v1/manual/performance-tips/).

Use grids of increasing size, paths/stars, disconnected and edgeless graphs, seeded random regular graphs, a dense bipartite case, and a small representative QAOA fixture. Vary both `n` and `k` because many small hamlets stress different work than two large ones. Include balanced and skewed user-provided partitions for cost benchmarks.

Warm each function before steady-state timing, use identical graphs and memberships for before/after scoring comparisons, and report median time, allocations, and objective values. Measure package-load and first-call latency in separate processes. Record Julia/dependency versions, machine, thread count, and RNG seed. Partition algorithm changes should report balance and objective quality as well as runtime; compare both matching and rank if both matter to the intended use.

## Documentation and experiment usability

Put a copyable numerical quickstart before visualization in the README. Explain that a villager is a graph vertex and a hamlet is a partition, then introduce the quantum interpretation. Include installation/development commands, input restrictions, result formats, and the exact versus heuristic distinction for every metric.

Repair specific drift: remove the expired August 2026 documentation promise and fill or remove the empty quantum-user section; correct the claim that `cut_rank_cost` lives in the sandbox; update the visualization docstring that says only naive cost is supported. Keep METIS comparisons in an example environment declaring Metis explicitly. Document illustrated results as examples, including graph/seed details, rather than as general performance guarantees.

Make randomness an argument: accept `rng` in `k_partition_random`, `best_sqrt_seed`, and callers. Define seed-sampling and tie policies. Gate the unconditional improvement messages at lines 402 and 425 behind `verbose` or logging.

Separate small runnable examples from expensive experiments. The scripts import packages absent from the library project, including Metis, GraphIO, EzXML, and QuantumClifford. Declare only the dependencies each environment needs. Put expensive runs behind explicit entry points so including helper definitions does not launch a sweep.

Use `joinpath(@__DIR__, ...)` for bundled resources, accept input/output paths as arguments, filter GraphML inputs, and avoid string concatenation that requires a trailing slash. In particular, [qaoa.py](../resources/qaoa.py) contains an author-specific absolute output path. Store experiment parameters, algorithm/objective names, seeds, and versions with results. The two tracked QAOA directories occupy roughly 145 MB on this checkout; retain a small test fixture and consider separate distribution for the full research data if clone/install size becomes burdensome.

Add a focused `.gitignore`, stop tracking `.DS_Store` files in a separate cleanup change, choose and document the intended license, and provide short contribution/test/benchmark instructions. Preserve figures and data that are intentional research outputs.

## Suggested sequence of changes

1. **Regression and correctness change.** Establish the test entry point; fix target validation, candidate initialization, one-hamlet matching, edgeless Saran–Vazirani handling, and display. Repair comparison scripts. Complete when the new regression tests and numerical README example pass.
2. **Source organization change.** Split the existing module by responsibility with unchanged supported calls. Complete when the same suite and examples pass without altered numerical results.
3. **API and membership change.** Introduce validated membership, scalar cost methods, explicit details, RNG handling, and compatibility wrappers. Complete when invalid inputs fail predictably and both interfaces agree on supported inputs.
4. **Usability and dependency change.** Extract optional plotting; add declared examples and a concise guide. Complete when a numerical-only session works and a plotting-enabled session can save an example figure.
5. **Measured optimization changes.** First benchmark unweighted matching and pair extraction, then rank and BURY. Make each optimization independently reviewable, with before/after allocations, runtime, and correctness/quality checks.

The first implementation change can remain narrow; deciding on every future API name or optional backend need not delay the immediate regressions and fixes.

## Validation notes and reproducible examples

The package imported successfully with Julia 1.12.6 using the existing project and manifest. No production source, dependency files, or existing research edits were changed for this review. The review did not execute the full research datasets, audit the paper's proofs, benchmark speedups, or establish support across the entire declared Julia/dependency range.

Observed passing results:

| Fixture | Cut edges | Matching sum | Pairwise GF(2) rank |
| --- | --- | --- | --- |
| README 6×6 grid, `bury_k(g, 3)`, sizes 12/12/12 | 15 | 9 | 9 |
| `K₂,₂`, its natural two parts | 4 | 2 | 1 |
| Two disjoint cross edges | 2 | 2 | 2 |
| Empty graph with two declared hamlets | 0 | 0 | 0 |
| One edge with three declared hamlets | 1 | 1 | 1 |

All 64 simple graphs on four vertices and all 16 assignments to two hamlets were checked: 1,024 cases passed comparison against an independent exhaustive matching enumerator and independent Boolean XOR rank calculation. This included the metric inequality and highlighted-edge count checks. A separate exhaustive BURY probe over simple graphs with 1–5 vertices, all valid targets, explicit seeds, and valid divisors for `k` passed 43,136 cardinality/coverage checks. These finite checks support the tested properties; they do not prove heuristic optimality.

The following small examples use the current API. They are diagnostic observations, not assertions to preserve after fixing the defects:

```julia
using QuantumHamlets, Graphs
const Q = QuantumHamlets

# One hamlet: currently throws; desired matching cost is zero.
land = Q.quantumLand(Dict(1 => 1, 2 => 1, 3 => 1), 1, 3)
Q.matching_cost_of_partition(land, path_graph(3))

# Run separately after the preceding exception.
# Currently returns no candidate despite a feasible balanced partition.
Q.bury_not_balanced_bipartition(
    complete_graph(8), 4, 4; cost_f=Q.naive_cost_of_partition
)

# Currently returns an empty registry for four assignable vertices.
Q.k_partition_saran_vazirani(SimpleGraph(4), 2)

# Currently prints to stdout and captures an empty string.
sprint(show, Q.villager(1, 1))
```

Do not run `bury(path_graph(4), 5)` in an interactive session before adding the range check; it was tested with an external timeout.

For a metric distinction worth keeping in the guide:

```julia
g = complete_bipartite_graph(2, 2)
land = Q.quantumLand(Dict(1 => 1, 2 => 1, 3 => 2, 4 => 2), 2, 4)
first(Q.naive_cost_of_partition(land, g))     # 4
first(Q.matching_cost_of_partition(land, g))  # 2
first(Q.cut_rank_cost(land, g))               # 1
```

For the greedy-cover counterexample, construct 13 vertices with hamlet 1 containing `1:4`, hamlet 2 containing `5:13`, and edges `(i,5)`, `(i,4+2i)`, `(i,5+2i)` for each `i in 1:4`. The current cover result is `(5,5)`, while `1:4` itself is a valid cover of size 4 and the matching cost is 4.
