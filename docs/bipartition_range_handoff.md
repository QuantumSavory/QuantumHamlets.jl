# Handoff prompt for BURY target range improvements

Copy the prompt below into a new agent conversation.

---

Work in `/Users/micciche/Research/QuantumInfo23/JuliaProjects/QuantumHamlets.jl`.

Improve `bury_not_balanced_bipartition` and make near-balanced bipartitions, including odd-vertex graphs, easy to request through the graph-first API. Start by reading:

- `docs/bipartition_range_review.md` — findings, confirmed counterexamples, and proposed API.
- `docs/usage.md`, `docs/bug_audit.md`.
- `src/bury.jl`, `src/bury/heap.jl`, `src/bury/scan.jl`, `src/bury/legacy.jl`.
- `src/costs.jl`, `src/api.jl`, `src/partition.jl`, `src/quantum_land.jl`.
- `test/bury_backends.jl`, `test/bury_validation.jl`, `test/api.jl`, and applicable `AGENTS.md` files.
- `experiments/cost_performance/README.md` for current scalar cost performance.

The review was made at commit `a05adb5`. Inspect the current working tree before editing; preserve all unrelated work. Do not reset files, change dependency versions, rewrite `test/reference/bury_original.jl`, or overwrite historical benchmark/audit logs. The bug fixes and cost optimizations are already implemented. Do not reintroduce the old bounded-score initialization, recursive matching traversal, or weighted unit-weight matching path.

## Current behavior to preserve

`bury_not_balanced_bipartition(g, min_v, max_v; backend, bury_v, cost_f)` validates integer ordered bounds in `0:nv(g)`, evaluates every target in ascending order, and returns `(registry, colored)`. Colors identify label 2. It minimizes the first component of `cost_f(land,g)` and keeps the LAST candidate on equal scores.

The exact built-in `bury` and `bury_seeded` callbacks receive the backend keyword. A custom `bury_v(graph,target)` receives exactly two positional arguments and controls its own backend. Preserve callback count/order and exceptions; do not catch `MethodError` and retry a different signature. Custom cost callbacks use the old land-first tuple contract.

All three backends are supported. Preserve valid-input greedy choices, explicit-seed behavior, return formats, traces, and default-RNG consumption. Ordinary empty BURY and the low-level empty range `0:0` are valid; seeded empty calls reject before sampling.

`partition_graph(g,k)` currently promises exactly equal parts and requires divisibility. `bury_k_seeded` uses a minimum-degree seed per residual graph, not sampled seed search. Do not rewrite the multipartition driver.

## Implement this first pass

1. **Fix retained callback output ownership.** The range function currently stores mutable callback outputs by reference. The review has a reproducer whose custom partitioner reuses a dictionary and color vector: scores uniquely favor target 1, but the returned result describes target 3 after later mutation. Snapshot retained winners on the generic path. Any reused internal workspace must also have independent winning state. Add a desired-behavior regression, including registry/color agreement and mutation after the search returns.

2. **Use the scalar cost kernels for recognized built-ins.** The existing loop constructs a registry, a copied/validated `quantumLand`, hamlets/villagers, dense labels, and unused edge colors for every trial. Ordinary heap/scan kernels can produce colored membership directly. Score known matching/edge-cut/rank objectives through existing private kernels, retain the winning membership, and materialize only the returned representation. Keep the generic callback path and its existing validation. Keep legacy's algorithm intact. Do not silently treat arbitrary closures as pure built-ins or skip their invocations.

3. **Add an exported graph-first convenience helper.** Recommended spelling:

   ```julia
   g = path_graph(11)
   p = bipartition_graph(g; backend=:heap)
   p = bipartition_graph(g; targets=4:7, objective=pairwise_cut_rank)
   ```

   Return `GraphPartition(...; nparts=2)`. Default targets are `fld(n,2):cld(n,2)`, evaluating one candidate for even `n` and both candidates for odd `n`. Require `n >= 2` and explicit targets to be a nonempty contiguous integer range in `1:n-1`, giving two nonempty parts. Validate before callbacks. Accept graph-first scalar `objective(g,p)`, defaulting to matching. Preserve ascending target order and last-target ties. Support all three backends.

   Scope this new helper to ordinary BURY initially. Seeded/custom workflows continue through the existing qualified range function. Share a small internal selection mechanism where useful, without turning this into a broad representation or partition-driver rewrite. Preserve the legacy public name, tuple return, and zero/empty-side boundary behavior. Do not silently change `partition_graph`'s exact-divisibility contract.

4. **Document the seeded objective composition.** `bury_v=bury_seeded` invokes a sampled seed search that always selects by matching cost, even when the outer `cost_f` is rank or edge cut. Explain that it selects among matching-selected candidates per target. Do not change that established policy in this pass. Objective-aware seed selection, common seed pools, unique sampling, and reuse of the selected inner partition/score are follow-up work requiring explicit compatibility choices.

5. **Measure and document the result.** Compare the legacy wrapper before/after and the new convenience helper on small/moderate graphs, with odd near-balanced ranges, a singleton target, and a modest wider interval. Report elapsed time and allocations including scoring/representation work. Use the existing experiment environment and external time limits. Record slowdowns as well as gains. Keep new evidence in a new directory.

## Counterexamples that must remain correct

- **Odd targets can differ in quality:** `SimpleGraph(5)` with edges `(1,4)` and `(2,3)` has isolated vertex 5. Ordinary BURY at target 2 colors `{1,5}` with matching cost 1; at target 3 it colors `{1,4,5}` with cost 0. The new default helper must choose the latter, with sizes `[2,3]`. Test all backends.
- **Target results are not nested:** `SimpleGraph(4)` with edges `(1,3),(1,4),(2,3),(2,4),(3,4)` gives colored sets `{1,2}` for target 2 and `{1,3,4}` for target 3. Do not derive all targets from one greedy prefix.
- **Endpoints are degenerate:** default matching search over `0:n` chooses target `n`, a zero-cost empty-sided partition, under last-candidate ties. Preserve this low-level behavior while making the new convenience API reject empty-sided ranges.
- **Winner aliasing:** use the shared-buffer reproducer from the review, then mutate the callback's buffers after return. The winning registry and colors must remain an independent, consistent target-1 result.

Do not skip complementary target sizes: complement symmetry of an objective does not imply complementary BURY outputs. Do not assume monotonicity/unimodality of target costs. Do not early-return at zero: it changes last-target ties, callbacks, and seeded RNG consumption. If investigating matching lower-bound pruning later, the OUTER last-tie search can reject only when the lower bound is strictly GREATER than the incumbent; equality must still be resolved. Inner seed search uses first-best ties, so its condition differs.

## Tests and completion criteria

- Add tests for odd/even graphs, `n=0/1/2`, explicit nonempty ranges, invalid/reversed/floating bounds, and graph nonmutation.
- Compare results to an explicit loop over the exact same targets, applying the documented last-tie policy. Test built-in objectives, custom graph-first scalar objectives for the new helper, and legacy tuple callbacks.
- Preserve and extend callback call-order/count tests, backend forwarding, negative/high/zero scores, seeded RNG and output comparisons, and original public types.
- Test the shared-buffer ownership defect and the two small graph counterexamples above.
- Validate failures before sampling or invoking callbacks. Preserve the current legacy callback validation rather than bypassing it in a fast path.
- Run `julia --startup-file=no --project=. -e 'using Pkg; Pkg.test(; allow_reresolve=false)'`. The existing experiment environment can also run `julia --startup-file=no --project=experiments/bury_performance test/runtests.jl`.
- The latest full suite passed 457,209 assertions; the independent cost audit passed 127,098 assertions on 7,122 fixtures. Assertion counts may increase. A focused 9,191 cost-semantics rerun passed after the final cost-loop restoration.
- Avoid full research sweeps, large automatic-seed searches, or `resources/graphml*` workloads. Preserve `docs/bury_approximation.md` and `experiments/bury_approximation/`; their prototypes are not part of this implementation pass.
- Update usage documentation and the review's status so proposed APIs are distinguished from implemented ones. Summarize code changes, measured gains/tradeoffs, tests, and any deferred items.

Proceed autonomously with routine implementation choices. The current task is scoped implementation, tests, and documentation; no dependency upgrades, publishing, or unrelated refactoring.
