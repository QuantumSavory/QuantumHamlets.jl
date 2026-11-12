# GraphML METIS versus BURY comparisons

For the redesigned family figures comparing **BURY, seeded BURY, and METIS**,
use the [family comparison guide](graphml_family_plots_README.md). It includes
fixed-seed ranges, complete-family findings, and counterexamples. The original
all-graph figures and their evidence below are retained as historical output.

`graphml_plots.jl` compares METIS and ordinary BURY on every `.graphml` file
directly inside `resources/graphml*` directories. It produces scientific
figures in PNG and SVG, together with the underlying scores and actual part
sizes. The presentation follows the METIS/BURY comparisons in `grid_plots.jl`
and `qaoa_plots.jl`, with named graph detail panels and explicit balance checks.

## Completed run (12 November 2026)

[run_20261112T142638Z](figs/graphml_comparison/run_20261112T142638Z/) completed
with the default heap backend, both objectives, and part counts 2 through 8.
All **45 input graphs** were processed. The **315 status rows** comprise
**78 completed comparisons**, **236 unsupported capacity combinations**, and
**one METIS empty-part result**. The **316 raw score rows** include the two
flagged scores for that empty-part METIS output; those points are excluded
from the figures while BURY's valid scores are retained.

| Figure | PNG | SVG |
| --- | --- | --- |
| Matching overview | [View](figs/graphml_comparison/run_20261112T142638Z/overview_matching.png) | [Vector](figs/graphml_comparison/run_20261112T142638Z/overview_matching.svg) |
| Pairwise rank overview | [View](figs/graphml_comparison/run_20261112T142638Z/overview_rank.png) | [Vector](figs/graphml_comparison/run_20261112T142638Z/overview_rank.svg) |
| Initial graph bipartitions | [View](figs/graphml_comparison/run_20261112T142638Z/graphml_initial_bipartitions.png) | [Vector](figs/graphml_comparison/run_20261112T142638Z/graphml_initial_bipartitions.svg) |
| Reduced graph bipartitions | [View](figs/graphml_comparison/run_20261112T142638Z/graphml_reduced_bipartitions.png) | [Vector](figs/graphml_comparison/run_20261112T142638Z/graphml_reduced_bipartitions.svg) |
| Evolution graph bipartitions | [View](figs/graphml_comparison/run_20261112T142638Z/graphml_evolution_bipartitions.png) | [Vector](figs/graphml_comparison/run_20261112T142638Z/graphml_evolution_bipartitions.svg) |

All ten PNG/SVG exports are present. The five PNGs were visually inspected;
labels, axes, and legends are visible without clipping. Numerical evidence is
available in [results.tsv](figs/graphml_comparison/run_20261112T142638Z/results.tsv),
[status.tsv](figs/graphml_comparison/run_20261112T142638Z/status.tsv), and
[inputs.tsv](figs/graphml_comparison/run_20261112T142638Z/inputs.tsv).

Balance differences affect interpretation: **48 of the 78 nonempty METIS
graph/part-count results** fall outside BURY's exact or nearest-balanced
constraint and are marked with red rings. The empty-part case is
`graphml_evolution/k1.graphml` with five requested parts: METIS returns sizes
**`0,6,6,6,7`**. Its raw matching and rank costs are preserved for diagnosis,
not counted as valid five-part comparisons.

The [metadata](figs/graphml_comparison/run_20261112T142638Z/metadata.txt) records
Julia 1.12.6, one thread, Metis 1.5.0, package versions, options, and source
hashes. The [log](figs/graphml_comparison/run_20261112T142638Z/output.log)
records successful completion in **29.691 seconds**, including loading and
rendering; this elapsed time is not a partitioner benchmark.

The final focused validation passed **1,379 assertions in 7.7 seconds**; see
[validation.log](figs/graphml_comparison/run_20261112T142638Z/validation.log).

## Run

From the package root, using the already installed parent `JuliaProjects`
environment:

```sh
python3 resources/run_graphml_plots.py
```

Defaults are `--objectives matching,rank --backend heap --parts 2,3,4,5,6,7,8`.
For a narrower run:

```sh
python3 resources/run_graphml_plots.py --objectives rank --parts 2 --backend scan
```

The runner creates a fresh timestamped directory under
`resources/figs/graphml_comparison/`. Supply `--output-dir PATH` to select a
different fresh directory; existing directories are rejected. It launches
Julia with `--startup-file=no --threads=1 --project=..` from the package root.
It does not install, resolve, or upgrade dependencies. Existing figures and
historical measurements are preserved.

An external watchdog allows 180 seconds for package loading, 90 seconds per
load/comparison/plot stage, and 900 seconds overall. A failure or timeout kills
the worker process group and preserves the partial artifacts and `output.log`.
A complete run ends with `COMPLETE` in that log. Including `graphml_plots.jl`
defines `GraphMLComparisonPlots` and its helpers without launching the sweep.

## Inputs and partition policies

The current inventory contains **45 graphs**: 20 initial, 20 reduced, and five
evolution graphs. GraphIO loads their simple, undirected, unweighted topology.
XML node declaration order maps vertices to `1:n`; node phase/type/layout and
edge-type attributes do not enter these graph objectives. Each input's path,
vertex/edge counts, and SHA-256 are saved in `inputs.tsv`.

For **two parts**, every graph with at least two vertices is included. BURY
uses `bipartition_graph(g; objective, backend)`: on odd graphs it independently
generates both targets `floor(n/2)` and `ceil(n/2)`, scores them using the
selected objective, and retains the last target on ties. Thus matching and
rank can select different BURY bipartitions. Even graphs have one target.

For **three or more parts**, the comparison includes only graph/part-count
pairs where `n % k == 0` and `k <= n`. BURY uses `partition_graph(g, k)` and
returns exactly equal parts. The same BURY membership is scored by both
objectives. Unsupported pairs receive an `unsupported_balance` row in
`status.tsv`; they are not silently padded, truncated, or repartitioned under
a different capacity policy. With the current inputs and default part counts,
there are **79 supported pairs** (45 bipartitions and 34 multipartitions) and
**236 unsupported pairs**. These inventory counts specify requested feasible
coverage, not evidence that a particular run finished or that METIS returned
nonempty parts for every supported pair.

METIS uses native **KWAY**, equal target part weights, native **seed 0**, and
`UFACTOR=30`, the standard **3% imbalance tolerance**. The wrapper restores the
caller's METIS options after each invocation. One METIS partition per graph
and part count is reused for all requested scores. METIS optimizes **edge
cut** internally; it does not optimize matching or GF(2) rank. Both algorithms'
returned partitions are evaluated with the same selected scalar cost function:

- `matching`: `matching_cost(g, p)`, the sum of maximum matching sizes across
  unordered pairs of parts.
- `rank`: `pairwise_cut_rank(g, p)`, the sum of GF(2) biadjacency ranks across
  unordered pairs of parts. For more than two parts, this is a sum of pairwise
  ranks, not a single bipartite cut rank.

Lower scores are better for these objectives. A cost comparison must also
account for balance: METIS's tolerance can produce different capacities from
BURY's exact or nearest-balanced constraint. `results.tsv` records every part
size and `max_size_over_mean = maximum(part_sizes)/(n/k)`. Red rings in the
figures identify METIS partitions outside BURY's balance requirement: a size
difference greater than one for `k=2`, or unequal sizes for divisible `k>=3`.
The comparison checks that graphs are unchanged.

METIS can return an **empty part** even for a supported graph/part-count pair.
Such a return is not a valid comparison with BURY's nonempty-part requirement.
The runner preserves its raw costs and sizes in `results.tsv` with
`nonempty=false`, records `metis_empty_part` in `status.tsv`, and excludes that
METIS point from the figures. BURY's valid results remain in both the table and
plots. Raw empty-part scores must not be counted as successful METIS results
when summarizing or comparing algorithms. Nonempty METIS results that miss
BURY's stricter balance constraint remain plotted with red rings.

## Outputs and figure reading

The default two-objective run produces **five figures**, each in PNG and SVG:

| File stem | Content |
| --- | --- |
| `overview_matching` | Matching cost versus vertex count; one panel per collection |
| `overview_rank` | Pairwise cut rank versus vertex count; one panel per collection |
| `graphml_initial_bipartitions` | Named initial graphs, paired METIS/BURY scores for two parts |
| `graphml_reduced_bipartitions` | Named reduced graphs, paired METIS/BURY scores for two parts |
| `graphml_evolution_bipartitions` | Named evolution graphs, paired METIS/BURY scores for two parts |

Overview colors indicate the number of parts; circles indicate BURY and
crosses indicate METIS. Named bipartition panels use blue BURY circles and
orange METIS crosses, connected within each graph for comparison. Graphs are
ordered by vertex count and then name. Both plot types use red rings for
METIS balance differences. Restricting `--objectives` creates only the selected
overview and objective panels; omitting `2` from `--parts` omits bipartition
detail figures.

| Evidence file | Content |
| --- | --- |
| `results.tsv` | Collection, graph, vertices, edges, part count, objective, algorithm, raw cost, actual sizes, balance flag, imbalance ratio, and `nonempty` flag |
| `status.tsv` | One row per requested graph/part-count pair: `completed`, `unsupported_balance`, or `metis_empty_part`, with reason |
| `inputs.tsv` | Source paths, topology sizes, and source-file hashes |
| `metadata.txt` | Julia/package versions, thread count, algorithm settings, per-objective bipartition selection policy, topology interpretation, and implementation/backend hashes |
| `output.log` | Exact command, watchdog limits, stage progress, errors if any, and completion marker |

A complete default run is expected to contain **316 raw score rows**: 79
supported pairs × two algorithms × two objectives, including any flagged
empty-part METIS outputs. The count of `completed` pairs can therefore be less
than 79; inspect `status.tsv` and `nonempty` before calculating comparison
statistics. The TSV files are the authoritative numerical results; plots can
overlap points with identical vertex counts or scores. Elapsed wall time in
`output.log` includes loading, scoring, and figure rendering and is not a
partitioner performance benchmark.

The initial run, [run_20261112T142507Z](figs/graphml_comparison/run_20261112T142507Z/),
is preserved as incomplete evidence. It stopped when METIS returned an empty
part for `graphml_evolution/k1.graphml` with five requested parts. The current
policy records this condition and continues while keeping the invalid METIS
comparison out of the figures. This initial run is not a completed result set.

## Focused validation

```sh
julia --startup-file=no --project=.. resources/graphml_plots_validation.jl
```

The validation checks all 45 source files against the existing inventory's
hashes, independently stored topologies, and graph statistics. Small fixtures
check supported capacities, odd-target selection, objective scoring, METIS
repeatability and option restoration, graph nonmutation, and preservation of
existing output tables. It loads the plotting source without running the full
comparison. Successful validation and visual inspection of generated figures
should be recorded with the completed run rather than inferred from the
expected counts above.
