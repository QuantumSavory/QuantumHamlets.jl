# QuantumHamlets.jl

QuantumHamlets partitions graphs for distributed graph-state preparation. It provides the BURY heuristic, comparison algorithms, and costs based on cut edges, pairwise maximum matchings, and GF(2) cut ranks. The same partitioning tools can be used for ordinary graph applications.

A graph vertex is called a **villager** and a partition is called a **hamlet** in the original API. The simpler API works directly with a graph and a `GraphPartition`; the original functions and `quantumLand` remain available.

## Set up a local checkout

Julia 1.11 or later is required. From the package directory:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=.
```

To use this checkout from another Julia project, run `Pkg.develop(path="/path/to/QuantumHamlets.jl")` in that project's environment.

## Partition a graph and evaluate its cost

```julia
using Graphs, QuantumHamlets

g = grid([6, 6])
p = partition_graph(g, 3)

partition_sizes(p)     # [12, 12, 12]
edge_cut_cost(g, p)    # 15
matching_cost(g, p)    # 9
pairwise_cut_rank(g, p) # 9
```

`partition_graph` defaults to BURY with the heap backend. Its balanced partitioners require a nonempty simple undirected graph and a positive partition count that divides the number of vertices. You can also construct `GraphPartition` from your own membership vector or vertex-to-partition dictionary.

Choose a BURY backend explicitly when needed:

```julia
p = partition_graph(g, 3)                  # backend=:heap
p = partition_graph(g, 3; backend=:scan)    # Optimized minimum scans
p = partition_graph(g, 3; backend=:legacy)  # Original implementation
p = partition_graph(g, 3; algorithm=:bury_seeded, backend=:scan)
```

All three backends preserve BURY's partition choices on supported inputs. Heap is a useful default for larger sparse graphs; scan can be faster for small or dense inputs, particularly with explicit seeds. Legacy remains available for comparison and compatibility. The choice also applies to the original `bury` functions and their seed searches; see [backend selection](docs/usage.md#choose-a-bury-backend).

Plotting is a separate step:

```julia
figure, cost = plot_partition(g, p; highlight=:matching)
figure
```

See the [usage guide](docs/usage.md) for algorithms, metric definitions, custom partitions, plotting, and migration from the original API. The [revision proposal](docs/revision_proposal.md) records the broader review and planned follow-up work.

## Tests and current scope

From the package directory:

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
```

The tests compare all BURY backends with a frozen original implementation, including exhaustive small graphs, seeded behavior, RNG consumption, and verbose output. They also check partition and cost results and the simpler API. The confirmed audit defects now have passing regressions, including validation, matching witnesses, component packing, and constructor/display behavior. See [known limitations](docs/usage.md#known-limitations).

The [BURY performance study](docs/bury_performance_findings.md) and [GraphML study](docs/graphml_backend_findings.md) record the prototype measurements behind the backend choice. Their figures describe the measured workloads, rather than a universal speedup.

The [production backend checks](experiments/bury_backends/README.md) verify the integrated API and record timings with the original public return types.

## Citing

If you use this work, please cite [Quantum Hamlets: Distributed Compilation of Large Algorithmic Graph States](https://arxiv.org/abs/2603.06387).

```bibtex
@misc{micciche2026quantumhamletsdistributedcompilation,
      title={Quantum Hamlets: Distributed Compilation of Large Algorithmic Graph States},
      author={Anthony Micciche and Naphan Benchasattabuse and Andrew McGregor and Michal Hajdušek and Rodney Van Meter and Stefan Krastanov},
      year={2026},
      eprint={2603.06387},
      archivePrefix={arXiv},
      primaryClass={quant-ph},
      url={https://arxiv.org/abs/2603.06387},
}
```
