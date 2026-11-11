# Independent numerical cost audit

The audit passed **127,098 assertions on 7,122 fixtures** with no additional
small-graph cost defect found. Production source and existing tests were not
modified. Run from the package root with:

```sh
python3 experiments/bug_audit/costs/run.py
```

`run.py` starts Julia in the existing isolated experiment environment and
enforces a 120-second external process-group timeout. Results are recorded in
`oracles.log`; the final successful run exited normally. The test portion took
1.5 seconds on the current environment, excluding package loading/compilation.

Coverage:

- Every labeled simple graph on 1–4 vertices, with every membership assignment
  to two or three declared parts, including empty parts.
- Empty graphs with two, three, or four declared parts.
- Three hundred deterministic random graphs on five/six vertices, with two
  through five declared parts, plus consistently relabeled copies.
- Exact edge-cut counts against direct edge filtering.
- Maximum matching cardinalities against an exhaustive bipartite matching
  search for each pair of parts; witness edges are checked for membership,
  disjoint endpoints within each pair, and correct cardinality.
- Pairwise GF(2) ranks against independent packed-row XOR elimination.
- Greedy cover sizes bounded by exhaustive minimum vertex cover and vertex
  count; the second cover estimate satisfies the matching lower bound and
  `(nparts-1)*cover_size` upper bound.
- Input-graph and partition preservation; exact cost invariance under vertex
  and part relabeling. Greedy cover values are not required to be invariant.

The existing one-part matching defect is reproduced separately:

```julia
matching_cost(path_graph(3), GraphPartition([1, 1, 1]))
# MethodError involving zero(::Type{Any}); expected cost is 0.
```

This comes from summing the empty `Any[]` at the end of the original matching
evaluator. The oracle sweep excludes one-part matching afterward so the known
failure does not obscure independent coverage. The audit does not test large
graphs, guarantee heuristic optimality, or establish the physical optimality
of the GHZ estimates.

## Confirmed recursion failure on a valid larger graph

A separate bounded check confirmed the earlier recursion risk as an actual
failure in the public matching API:

```julia
g = path_graph(50_000)
p = GraphPartition([isodd(v) ? 1 : 2 for v in vertices(g)])
matching_cost(g, p)
# StackOverflowError; expected cost is 25_000.
```

All path edges cross the two parts, and pairing `(1,2), (3,4), ...` attains
the upper bound of 25,000 matching edges. The recursive `flip_nodes` function
in `src/costs.jl` traverses a path to depth proportional to its vertex count.
The copied exact DFS completed for 20,000 vertices and overflowed at 50,000;
this is an observed case, not an exact or portable overflow threshold.

The full public call was then confirmed in a **fresh Julia worker**, without
any preceding stack overflow in that worker. A 1,000-vertex full call first
returned the expected 500 in 0.0086 seconds. The 50,000-vertex full call raised
`StackOverflowError` before the external five-second stage deadline. The
worker was closed afterward because Julia warns that state after a stack
overflow may be unreliable. Logs are `recursion.log` and
`recursion_confirm.log`.

The shared cost evaluator affects matching details, matching plots, and
seed-search scoring as well as `matching_cost`; choosing another BURY backend
does not remove this failure. An iterative traversal would avoid recursive
stack growth. In this evaluator, the two hamlet labels already define the
bipartite sides, so traversal can instead be removed entirely.

Bounded reproduction:

```sh
python3 experiments/bug_audit/costs/recursion_probe.py
python3 experiments/bug_audit/costs/recursion_probe.py 50000
```

The first command locates a copied-DFS failure with a maximum of 100,000
vertices, then tries the full evaluator. The second directly confirms the
full evaluator in a fresh worker after a small successful check. Both use an
external process-group watchdog: 120 seconds for loading/warmup and five
seconds per test stage. A small full call exceeding one second prevents
escalation. No GraphML seed search or large weighted matching sweep is run.
