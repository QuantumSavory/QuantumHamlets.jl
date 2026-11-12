# BURY target range review and implementation

Originally reviewed on 11 November 2026 against commit `a05adb5`. The scoped
first implementation pass is now complete in the working source. The historical
findings below distinguish the reviewed behavior from the implemented API.

`bury_not_balanced_bipartition` is a useful exhaustive search over **target sizes**, with one candidate partition generated for each size. It already supports graphs with an odd number of vertices. The first pass adds a graph-first API and removes per-trial representation work for recognized ordinary built-ins. It also fixes the confirmed ownership issue: reused callback buffers no longer overwrite the retained winner.

## Implemented first pass

- Exported `bipartition_graph(g; targets=fld(nv(g),2):cld(nv(g),2),
  objective=matching_cost, backend=:heap)` returns two nonempty declared parts.
  It validates an ascending contiguous integer range before callbacks and
  keeps ascending evaluation and last-target ties.
- Ordinary BURY with the exact matching, edge-cut, or rank functions uses the
  existing private scalar kernels on Boolean membership. Heap/scan avoid
  trial registries, legacy lands, dense integer labels, and unused edge colors.
  The legacy algorithm body is unchanged and still constructs trial registries.
- Generic legacy callbacks keep their two-positional-argument partitioner and
  land-first tuple cost contracts, including `quantumLand` validation. Retained
  registries and colors are copied. New scalar callbacks receive a separate
  `GraphPartition` per target and are invoked on every target.
- Seed selection, RNG consumption, explicit seeds, traces, and the exactly
  balanced multipartition driver are unchanged. Automatic seeded range search
  still selects by matching internally, whatever the outer objective.

See [usage](usage.md#near-balanced-bipartitions-including-odd-graphs),
[regressions](../test/bipartition_range.jl), and the new
[bounded performance evidence](../experiments/bipartition_range_performance/README.md).
The full `Pkg.test(; allow_reresolve=false)` suite passed **459,086 assertions**,
including **1,877** new range assertions. All **120** benchmark comparisons
passed membership, scoring, and preservation checks. On the 225-vertex grid,
six-target heap matching improved from **306 to 107 μs** for the wrapper and
**92.8 μs** for the helper; allocated bytes fell from **464,816 to 103,808 and
74,352**, respectively. Generic callback copies cost **2–12%** in the measured
cases, and six legacy singleton rows allocated more bytes (up to **6.6%**).
The linked experiment retains every row and measurement limitation.

Objective-aware seed selection, shared/unique seed pools, inner-result reuse,
workspace reuse, pruning, and multipartition balance policies remain follow-up
work. The approximation study and prototypes are outside this pass.

## Odd graphs already work

The target is the number of colored vertices, which receive label 2. Label 1 contains the other vertices. There is no parity restriction in either `bury(g, target)` or the range function.

```julia
using Graphs, QuantumHamlets
const Q = QuantumHamlets

g = path_graph(5)
reg, colored = Q.bury_not_balanced_bipartition(
    g, fld(nv(g), 2), cld(nv(g), 2); backend=:heap,
)
p = GraphPartition(reg; nparts=2)
partition_sizes(p)  # [2, 3] for this fixture
matching_cost(g, p) # 1
```

For odd `n`, the inclusive range `fld(n,2):cld(n,2)` contains two targets. For even `n`, it contains one. A single `bury(g, fld(n,2))` also gives the nearest-balanced sizes, but considering both targets can produce a better partition.

That distinction is substantive. On five vertices with edges `(1,4)` and `(2,3)`, leaving vertex 5 isolated:

| Target | Colored vertices from ordinary BURY | Matching cost |
| --- | --- | ---: |
| 2 | `{1,5}` | 1 |
| 3 | `{1,4,5}` | 0 |

Searching `2:3` finds the zero-cut partition. Although matching cost is invariant under complementing a fixed bipartition, independently running BURY with complementary targets need not produce complementary sets.

The public `partition_graph(g, 2)` rejects odd `n`, because it delegates to the exactly balanced multipartition driver, whose contract requires `nv(g) % k == 0`. This is an API policy difference, not a failure of the bipartition algorithm. The new `bipartition_graph` preserves that established exact-balance contract in `partition_graph`.

## The role of the range function

The relevant layers are:

| Layer | Responsibility |
| --- | --- |
| Heap, scan, and legacy kernels | Generate one BURY candidate at a specified target and optional seed |
| `bury` and explicit-seed `bury_seeded` | Return the legacy `(registry, colored)` representation |
| Automatic `bury_seeded` | Sample seeds, select one by matching cost, and generate its partition |
| `bury_not_balanced_bipartition` | Generate and compare candidates across an inclusive target interval |
| `bury_k` and `bury_k_seeded` | Repeatedly extract exactly equal-sized parts from residual graphs |
| `partition_graph` | Return a `GraphPartition` for the balanced multipartition algorithms |
| Scalar cost API | Evaluate a supplied `GraphPartition`, including uneven and empty parts |

Range search is therefore a natural building block for near-balanced bipartitioning. It is not a replacement for the multipartition driver: choosing a locally favorable first cut does not establish that the final multipartition objective will be favorable, and arbitrary balanced `k`-partitioning needs an explicit capacity policy.

The function minimizes cost among its generated candidates, not over all partitions of those sizes. The [approximation study](bury_approximation.md) also cautions against giving BURY a general optimality or approximation guarantee.

## Balance constraints must remain explicit

The interval constrains the size of **label 2**, rather than directly constraining the smaller part. To require both sides to contain at least `m` vertices, use targets `m:n-m`. To allow an absolute size difference at most `d`, the feasible targets satisfy `abs(n - 2r) <= d`.

The low-level API intentionally permits targets zero and `n`. With default nonnegative matching scores, either yields a zero-cost empty-sided partition. In particular, searching `0:n` returns the target-`n` result because equal scores retain the last candidate. This is expected under that contract, but it is a poor default for a user requesting two nonempty parts.

Even after excluding the endpoints, a wide range trades balance against the objective. Matching cost is at most `min(r, n-r)`, so highly unequal sizes can have inherently small scores. Prefer a narrow, explicit balance interval. A balance penalty or a different tie rule would be a distinct policy change.

## Efficiency findings at the reviewed revision

The search loop is straightforward and retains only the current winner; it does not store all candidates. For `q = max_v - min_v + 1`, its ordinary work is approximately

`sum(T_bury(r) + T_score(r) + O(n) representation work for each target r)`.

For the two nearest-balanced odd targets, this is only two trials. For a wide interval, repeated runs and scoring scale with `q`.

At `a05adb5`, the loop in [src/bury.jl](../src/bury.jl) gave every candidate a `Dict{Any,Any}` registry, then a validating/copying `quantumLand` with hamlets and villagers. The cost wrapper converts that registry back to dense labels. Matching also constructs an edge-color vector that the range function discards. The newer scalar APIs and private kernels in [src/costs.jl](../src/costs.jl) already avoid much of this work.

The implemented first improvement is a **specialized path for recognized ordinary built-in partitioners and objectives**:

1. For ordinary heap/scan BURY, obtain the colored vector directly from the kernel.
2. Score dense membership with `_matching_cost`, `_edge_cut_cost`, or `_pairwise_cut_rank` as appropriate, without plotting colors or legacy membership objects.
3. Keep an independent winning membership and construct the legacy registry or `GraphPartition` only for the returned result.
4. Keep the generic callback path and its validation, calling convention, order, and side effects intact. Keep legacy available through its existing algorithm.

The new [range experiment](../experiments/bipartition_range_performance/README.md) measures these changes; the [cost benchmarks](../experiments/cost_performance/README.md) measure cost evaluators, not the full range workflow. They establish that a scalar path exists, not a measured speedup for this wrapper.

Reusing degree arrays or workspace between independent trials is a possible later improvement. All mutable per-run state must be reset, and retained winners must be copied before the workspace is reused.

## A confirmed ownership defect in custom callbacks

The reviewed winner update assigned `best_reg = bury_reg` and `best_colors = colors`. It kept references to the callback's outputs. If that callback fills the same buffers on the next trial, the saved winner changes even when the next score is worse.

Before the ownership fix, this reproducer returned target 3 although target 1 uniquely minimizes its score. It now retains target 1:

```julia
using Graphs, QuantumHamlets
const Q = QuantumHamlets

g = path_graph(4)
shared_reg = Dict(v => 1 for v in vertices(g))
shared_colors = fill(false, nv(g))

candidate = function (graph, target)
    for v in vertices(graph)
        shared_colors[v] = v <= target
        shared_reg[v] = shared_colors[v] ? 2 : 1
    end
    return shared_reg, shared_colors
end
score = (land, graph) -> (count(==(2), values(land.registry)), nothing)

reg, colors = Q.bury_not_balanced_bipartition(
    g, 1, 3; bury_v=candidate, cost_f=score,
)
count(colors) # 1; formerly 3 because later trials overwrote the winner
```

Each callback result is a valid registry/color pair of the requested size when returned. The defect is retaining it across later calls without taking ownership. The generic path now copies the winning registry and colors; the same ownership rule applies to any future reusable workspace. Built-in BURY currently creates fresh outputs and does not trigger this problem.

## Seeded range searches have two objectives

With `bury_v=bury_seeded`, every target invokes an automatic search of `ceil(sqrt(n))` sampled seeds. That inner search always scores with matching, then reruns BURY using the winning seed. The outer range search subsequently evaluates `cost_f`.

Consequently, `cost_f=cut_rank_cost` selects the best rank among **matching-selected candidates**, not the lowest-rank candidate across all sampled target-seed combinations. This is a limitation of the composition, not an incorrect rank calculation. The usage guide and range docstring now document it.

An objective-aware seed search would be useful as an explicit extension. Retaining the best trial partition and score could also avoid rerunning and rescoring it when the inner and outer objectives agree. Preserve the existing sampled-seed order, RNG draws, first-best seed ties, improvement messages, and last-best target ties unless a new policy is explicitly selected. Shared seed pools, distinct-seed sampling, and early stopping change the search policy and belong in a separate experiment.

## Optimizations that would change results

**One run cannot provide all target partitions as prefixes.** On four vertices with edges `(1,3), (1,4), (2,3), (2,4), (3,4)`, ordinary BURY colors `{1,2}` at target 2 and `{1,3,4}` at target 3. The remaining budget changes neighborhood-versus-singleton decisions. The sets are not nested.

**Do not skip complementary target sizes.** The five-vertex example above proves they can have different costs. Complement symmetry also cannot be assumed for arbitrary custom scores.

**Do not assume the target-score curve is monotone or unimodal.** No such property is established by this implementation. Binary/ternary search or a coarse target grid would need a separate approximate-search contract.

**Do not stop at cost zero under the legacy contract.** A later zero-cost candidate wins ties, and skipped calls can change callback side effects and seeded RNG consumption.

**Pruning must respect last-candidate ties.** For minimizing matching cost, a partial or greedy matching is a lower bound. An outer candidate may be discarded when that bound is strictly greater than the incumbent. Equality cannot be discarded: an exact tie must replace the incumbent, whereas a higher exact cost must not. The inner seed search has the opposite tie policy and can discard at equality. Do not reuse its pruning condition blindly.

## Implemented additive API

The exported graph-first helper implements the following contract:

```julia
g = path_graph(11)
p = bipartition_graph(g; backend=:heap)
p = bipartition_graph(g; targets=4:7, objective=pairwise_cut_rank)
```

- Returns `GraphPartition`, with two declared parts.
- Defaults to `fld(n,2):cld(n,2)` and evaluates both candidates on odd graphs.
- Requires `n >= 2` and targets within `1:n-1`, guaranteeing two nonempty parts. Zero/empty-side behavior is preserved in the existing low-level calls.
- Accepts an inclusive contiguous ascending integer target range, validating it before computation.
- Uses a graph-first scalar `objective(g, p)` callback. The existing `cost_f(land, g)` tuple contract is preserved in the legacy wrapper.
- Retains ascending target order and last-candidate ties, and supports all three BURY backends.
- Exposes ordinary BURY. Seeded/custom workflows remain available through the existing range function; extending the new helper to seeded search should address its objective policy explicitly.

This is now an exported API. A separate opt-in near-balance policy on `partition_graph(g, 2)` could later delegate to it. Do not silently relax the existing exact-divisibility rule for all `k`.

## Priorities and validation

The implemented first pass fixes callback result ownership, adds the ordinary built-in scalar path and graph-first helper, and documents the seeded objective limitation. Changes to seed sampling, multipartition capacities, and approximation/certificate prototypes remain deferred.

The original review's odd-path probes passed 28 assertions over sizes 3, 5, and 7 and all three backends, including cardinalities and agreement with explicit candidate enumeration. Tiny bounded probes confirmed the two target counterexamples and the reused-buffer defect. That original review did not measure performance; the subsequent implementation measurements are linked above.

The implementation regressions cover both counterexamples, odd/even/empty/singleton boundaries, zero/high/negative scores, last ties, buffer reuse, callback invocation order, invalid-input RNG preservation, backend forwarding, and comparison with an explicit exhaustive loop over target sizes. Test custom callbacks as callbacks; do not infer purity from their outputs. Run the standard package suite and benchmark only small/moderate fixed-target versus range workloads, with external time limits for seeded cases.

The [handoff prompt](bipartition_range_handoff.md) records the scope and preservation requirements for this implementation pass.
