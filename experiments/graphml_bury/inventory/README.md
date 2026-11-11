# GraphML inventory and property-based expectations

This inventory covers all 45 files directly inside `resources/graphml*`: 20
initial graphs, 20 reduced graphs, and five evolution graphs. The supplied XML
files total 4,696,087 bytes; the largest is 1,329,846 bytes. Inventory computation
is inexpensive. This does not imply that every partition or seed-search
workflow is inexpensive.

All 45 files describe connected, simple, undirected graphs: no self loops,
duplicate undirected edges, directed overrides, or isolated vertices were
found. Node IDs are numeric and ascending in declaration order, but none is
exactly `1:n`. Declaration order is retained when mapping to Julia vertex IDs.

| Directory | Vertices | Mean degree | Edge density | Structural expectation before timing |
| --- | ---: | ---: | ---: | --- |
| `graphml_initial` | 67–652 | 2.39–2.67 | 0.394%–3.62% | Low degree favors the heap as size grows; the smallest graphs may favor scan because setup and constant costs matter. |
| `graphml_reduced` | 25–294 | 3.49–8.07 | 2.75%–16.67% | Small vertex counts and higher degree make scan a plausible winner; heap is still worth testing at the larger end. |
| `graphml_evolution` | 25–1,772 | 3.36–13.66 | 0.771%–14.0% | Scan is plausible for the first small cases; heap becomes a stronger candidate as vertex count grows, but hubs increase update work. |

These are hypotheses about the selector, not measured winners. They do not
predict the cost of automatic seed search or the matching objective. Ordinary
and explicitly seeded BURY must be evaluated separately. The optimized scan
removes allocations from seeded selection, while the legacy ordinary minimum
scan can have favorable constants on small inputs.

Reduction always makes these paired graphs smaller in vertex count, but it
does not necessarily reduce the edge count. For example, the `gf2` multiplier
and `hwb6` reduced graphs have more edges than their initial counterparts.
Their degree distributions are heterogeneous: the reduced `gf2^6_mult` has
294 vertices, mean degree 8.07, and maximum degree 88. Evolution `k4` has 1,159
vertices and maximum degree 228; `k5` has 1,772 vertices and maximum degree 222.
A scalar density threshold loses this degree-distribution information and
ignores the requested target size, which affects how many selection rounds
BURY performs.

The XML includes node type, phase, and layout attributes plus an edge-type
attribute. No numeric edge-weight field is declared. Explicit edge-type data
is `2`; the declared default is `1`. The topology exports intentionally omit
these attributes because the tested BURY implementations accept simple
unweighted graphs. They remain recorded in `metadata.json`.

## Reproduction and files

From the package root, run:

```sh
python3 experiments/graphml_bury/inventory/inventory.py
```

- `inventory.tsv`: counts, density, degree statistics, components, input hashes,
  family/parameter ordering, XML paths, and topology/map paths.
- `initial_reduced.tsv`: paired count/density/mean-degree comparisons.
- `metadata.json`: key definitions, observed attribute-value counts, and ID
  samples/ranges.
- `topology/*.tsv`: the first line is the numeric pair `n<TAB>m`; subsequent
  lines are unique undirected edges with 1-based mapped endpoints.
- `node_maps/*.tsv`: each mapped vertex and its original GraphML ID.
- `inventory.log`: compact property output from the inventory run.

Degree percentiles use the nearest-rank definition. Density is
`2m / (n*(n-1))`. Connected components are computed with union-find. XML is
parsed with the Python standard library's streaming element iterator; source
files are read only. The script validates unique node IDs and declared edge
endpoints. Family parameters include `barenco_tof` 3/4/5/10, `tof` 3/4/5/10,
`gf2_mult` 4/5/6, and `evolution` 1/2/3/4/5 for progressive benchmark limits.
