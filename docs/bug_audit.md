# Current bug review

Reviewed on 11 November 2026 after production backend integration. The earlier
findings still matter, but their scope has changed: heap and scan now reject
invalid single-partition targets, and the simpler API protects several legacy
boundaries. The range-search score bound, one-part matching failure, legacy
state/display problems, experimental baseline behavior, and research-script
score handling remain unresolved.

This review also confirms a **stack overflow in matching on a valid large
path**, a nonexistent seed returned for an empty graph, and a broader
consequence of the Saran–Vazirani initialization problem. These are shared or
pre-existing issues, not a newly observed heap/scan partition disagreement.
Only new review scripts, logs, and this report were added; no production code
or existing tests were changed.

Examples below use `using Graphs, QuantumHamlets` and
`const Q = QuantumHamlets`.

## Earlier findings and their current status

### Range search can discard every valid candidate

**Impact:** a valid search returns `(nothing, nothing)` instead of a partition.
This affects all three backends.

```julia
Q.bury_not_balanced_bipartition(
    complete_graph(8), 4, 4;
    cost_f=Q.naive_cost_of_partition,
)
# (nothing, nothing)
```

The only candidate has a valid 4+4 partition and cut-edge cost 16. The search
initializes `lowest_cost` to the number of vertices, 8, so it never accepts
that candidate. See [bury.jl](../src/bury.jl), lines 141–152.

Initialize from the first candidate or use an unbounded initial score. Also
validate the requested interval: an inverted range currently returns no result
silently. The problem is especially relevant when using a custom objective
whose scale differs from matching cardinality.

### Matching cost fails for one declared part

**Impact:** the valid all-in-one partition raises an exception although its
matching cost is zero. The simpler `matching_cost` and `matching_details`
interfaces inherit this evaluator failure.

```julia
matching_cost(path_graph(3), GraphPartition([1, 1, 1]))
# MethodError involving zero(::Type{Any}); expected 0
```

There are no pairs of parts to evaluate. The evaluator sums an empty `Any[]`
instead of an integer accumulator, which has no defined additive identity.
See [costs.jl](../src/costs.jl), lines 110 and 145. An integer accumulator or
explicit zero result should preserve all edges as unhighlighted in this case.

### Saran–Vazirani skips the initial component assignment

**Impact:** an edgeless graph gets no partition, and some already-disconnected
graphs receive unnecessary cuts.

```julia
Q.k_partition_saran_vazirani(SimpleGraph(4), 2)
# Empty registry and trace, despite a feasible 2+2 partition.
```

Assignments are attempted only after entering the edge-cut loop. With no
edges, that loop never runs. The simpler `partition_graph(...;
algorithm=:saran_vazirani)` catches the incomplete registry and raises an error;
it does not repair the missing partition. See [baselines.jl](../src/baselines.jl),
lines 78–88.

**Newly confirmed consequence:** take ten vertices with edges
`(1,2), (2,3), (4,5), (5,6), (7,8), (9,10)`. Their components have sizes
3, 3, 2, 2. The two sets `{1,2,3,7,8}` and `{4,5,6,9,10}` give a balanced
zero-cut partition. The observed result instead cuts edge `(1,2)`, costing 1.
The function removes an edge before checking whether the original components
already fit. This extends the same initialization defect; it is not a claim
that every suboptimal heuristic result is a bug.

Attempt component packing before removing edges, and represent infeasibility
or solver failure explicitly. The counterexample and five other disconnected
fixtures are saved in the [integration audit](../experiments/bug_audit/integration/findings.txt).

### Legacy membership can become inconsistent through aliasing

**Impact:** mutating an input dictionary can change the registry while leaving
the stored hamlets and villagers unchanged.

```julia
r = Dict(1 => 1, 2 => 2)
land = Q.quantumLand(r, 2, 2)
r[1] = 2
land.registry[1]                       # 2
land.hamlets[1].villagers[1].address    # 1
```

The legacy constructor retains a `Dict{Int,Int}` input by reference while
building separate membership objects. See [quantum_land.jl](../src/quantum_land.jl),
line 57. Copying the input addresses this example; choosing one authoritative
membership representation addresses the wider synchronization problem.

`GraphPartition` copies inputs and its public accessors return copies, so it
avoids this particular aliasing issue. Its internal fields remain documented
as read-only. Legacy constructor validation remains weaker than that of
`GraphPartition`.

### Display writes to the wrong stream

**Impact:** captured display output is empty while text unexpectedly goes to
stdout. This affects logging, notebooks, and embedding in other displays.

```julia
sprint(show, Q.villager(1, 1))
# Returns "", while printing "Qubit 1, who lives in hamlet 1" to stdout.
```

All three methods in [display.jl](../src/display.jl) ignore their `io` argument.
Every nested `print`, `println`, or `show` must use that stream.

### Validation remains inconsistent at legacy entry points

The original oversized-target hang is now restricted to `backend=:legacy`:
`Q.bury(path_graph(4), 5; backend=:legacy)` can loop after every vertex has been
colored. The seeded legacy implementation has the same problem. The original
review already demonstrated the hang; this audit did not launch it again.
The retained bodies in [bury/legacy.jl](../src/bury/legacy.jl) still lack the
upper-bound guard. Heap and scan now reject these calls with `ArgumentError`.

Other low-level validation issues remain on all applicable backends:

- `Q.bury_k(path_graph(4), -2)` returns an empty dictionary.
- `Q.bury_k(path_graph(4), 0)` raises `DivideError` instead of an input error.
- `Q.k_partition_random(path_graph(5), 2)` raises `InexactError` while converting
  a fractional capacity to `Int`.

These are unsupported-input handling defects, distinct from wrong partitions
on valid input. `partition_graph` already validates positive, divisible `k`
before dispatch. A later cleanup should give the qualified entry points a
consistent validation policy as well.

### Research scripts mishandle scores and compare different objectives

Two problems remain in [grid_plots.jl](../resources/grid_plots.jl), lines 46–59:
the BURY scalar matching score is overwritten by the complete `(rank, nothing)`
tuple, and METIS continues to receive a matching score. Plotting the tuple data
fails in a tiny equivalent probe. Extracting the rank scalar alone would still
leave the two algorithms measured with different objectives.

In [sandbox.jl](../resources/sandbox.jl), lines 165, 172, 184, and 191, cut-rank
tuples are similarly collected where subsequent `mean`/`std` calls need scalar
numbers. The tiny equivalent `mean` probe raises a `MethodError` involving
division of `Tuple{Int64,Nothing}` by an integer.

Use one scalar objective consistently for every competitor, and label saved
results with that objective. The scalar APIs now make this easier. No complete
research sweep was run or existing saved research result reinterpreted.

## Newly confirmed issues

### Recursive matching traversal overflows on a valid graph

**Impact:** matching evaluation can fail on large sparse inputs before the
matching solver runs. It is shared by all BURY backends and can also affect
seed-search scoring, matching details, and matching-highlight plots.

```julia
g = path_graph(50_000)
p = GraphPartition([isodd(v) ? 1 : 2 for v in vertices(g)])
matching_cost(g, p)
# StackOverflowError in the tested environment; correct cost is 25_000.
```

The recursive `flip_nodes` traversal in [costs.jl](../src/costs.jl), lines
99–120, follows the entire path and exhausts the stack. The public call was
confirmed in a fresh Julia process, with no preceding stack overflow in that
worker. A 1,000-vertex call returned the correct cost 500 first. The 50,000-node
example is an observed failure, not a portable threshold.

The confirmation used a five-second external stage watchdog, which did not
fire; the worker exited after the overflow. See the
[fresh-worker log](../experiments/bug_audit/costs/recursion_confirm.log).
Avoid repeating this probe in an interactive research session: Julia warns
that a stack overflow can make further execution unreliable.

The two sides of each pairwise graph are already known from the part labels,
so the traversal can be eliminated. An iterative traversal is another repair.

### Empty seed searches return a nonexistent vertex

**Impact:** seed helpers report success with an invalid seed.

```julia
Q.best_seed(SimpleGraph(0), 0)       # 1
Q.best_sqrt_seed(SimpleGraph(0), 0)  # 1
```

Both initialize the result to 1 and perform zero iterations. This occurs on
all three backends. The subsequent `bury_seeded` call rejects that seed with
`ArgumentError` for heap/scan or raises `BoundsError` for legacy. See
[bury.jl](../src/bury.jl), lines 177–226. Reject empty seed-search input
explicitly or define a consistent empty-graph result.

A related lower-priority issue is that calling `best_sqrt_seed` directly with
an invalid target on heap/scan consumes a random draw before throwing. The omitted-seed
`bury_seeded` wrapper validates first, so the two entry points differ on failed
calls. This is an invalid-input/RNG policy inconsistency, not a wrong valid
partition. Both observations are reproduced in the
[API probe log](../experiments/bug_audit/api/probes.log).

## Checks that passed and limits of the review

The independent numerical audit passed **127,098 assertions on 7,122
fixtures**, using exhaustive matching and independent GF(2) XOR elimination
instead of treating the original evaluator as its own reference. It checked
matching witnesses, edge cuts, rank/matching/cut relationships, cover bounds,
consistent relabeling, and nonmutation. Empty graphs and declared empty parts
worked when at least two parts were declared. The known one-part matching
failure was checked separately and excluded from that sweep.

Another **407 assertions** covered supported integer representations,
zero-target/empty ordinary BURY, explicit seeds, and `GraphPartition` copy and
conversion behavior. The earlier-findings script passed 26 observation
assertions. Passing those assertions means the recorded findings are
reproduced; it does not mean the known bugs are fixed.

Twelve tiny Saran–Vazirani cases returned complete balanced partitions; all
six disconnected fixtures also remained balanced, with the zero-cut issue
described above in one case. Normal no-highlight plots, legacy rank plots,
and empty no-highlight plots constructed successfully.

I did not classify greedy vertex-cover nonoptimality as a bug: it is now
documented as an estimate. Passing the vertex-cover evaluator to the old
plotting callback still misinterprets its second numeric cost as an edge
color, but that combination is already documented as incompatible. The newer
plotting interface avoids it. A suspected `Figure(yticks=...)` failure was
ruled out in the installed plotting stack.

The checks used Julia 1.12.6 and the existing isolated experiment environment.
No new valid-input heap/scan partition disagreement or additional small-graph
numerical error was found. Finite tests do not establish correctness for every
graph, and no large automatic-seed or research-script sweep was attempted.

## Suggested repair order

1. Fix the one-part matching result and replace the recursive matching traversal.
2. Fix range-search initialization and input validation, including empty seed
   searches and low-level `k` guards.
3. Attempt Saran–Vazirani component packing before cutting edges.
4. Repair legacy constructor ownership and display streams.
5. Repair scalar objective handling before using the affected comparison scripts.

Reproducers, logs, and preservation checks are listed in the
[audit README](../experiments/bug_audit/README.md).
