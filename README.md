# QuantumHamlets (Currently Repo is a bit unstable as it is under active improvement - we should have most things tidied up by the end of the month (November 2026))
Please feel free to email us questions if you would like to use this library and the current state of documentation is stopping you. More examples and API details are in the [usage guide](docs/usage.md).

To use a local checkout, with Julia 1.11 or later, run these commands from the package directory:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=.
```

## Citing
If you would like to cite this work, please cite our paper: https://arxiv.org/abs/2603.06387

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
## For a non-quantum user, interested in regular graph partitioning for other applications
Our paper, "Quantum Hamlets: Distributed Compilation of Large Algorithmic Graph States", looks at partitioning certain classes of quantum states that can be represented as graphs for use in distributed quantum computing.

In that context, rather than minimizing edges cut, it is better to minimize the sum of the maximum matching sizes between pairs of partitions (for the purpose of reducing the amount of long range quantum entanglement needed).

In this pursuit, we developed a new heuristic algorithm for partitioning graphs with this goal in mind. To our knowledge, we were the first to consider this problem of graph partitioning with matchings as the objective, and the NP hardness proof is also present in our paper.

We hope that our abstractions were clear enough that Example 1 below shows you how to use our heuristic, BURY, for your non-quantum application.

## For a user interested in quantum graph state partitioning

In the language of our paper, graph vertices are "villagers" and each group of vertices is a "hamlet" hosted on a QPU. The example below shows how to assign villagers to hamlets and compare the inter-hamlet entanglement costs of preparing the desired graph state.

## Example 1: Partitioning a Grid into $k$ pieces
First, we need to define how many hamlets ($k$) and villagers per hamlet ($n/k$) we are considering:
```julia
numHamlets = 3
villagersPerHamlet = 12
n = numHamlets * villagersPerHamlet
```
We also need a graph to consider:
```julia
using Graphs
g = grid([6, 6])
```
QuantumHamlets.jl has a simple built in plotting function for graphs. Most of its functionality is for visualizing partitions, but it works just as well for looking at unpartitioned graphs. It is internally called with more parameters by `plot_partition`, which shows up later in this example.
```julia
using QuantumHamlets
QuantumHamlets.print_graph(g)
```
![Our example graph before we partition it](./resources/figs/unpartioned6_6_grid.png)
Now we need to create a "registry of villagers". That is, an assignment that takes integers from 1 to $n$, and maps each to one of the $k$ hamlets. QuantumHamlets.jl has a few different partition methods to choose from, and the current interface returns this assignment as a `GraphPartition`. Here is how one would use the BURY algorithm from our paper:
```julia
bury_part = partition_graph(g, numHamlets)
```

This object tracks which villager node is in which hamlet. We can inspect the assignment with `membership`: entry `v` is the hamlet containing vertex `v`. We can also check how many villagers ended up in each hamlet:
```julia
membership(bury_part)
partition_sizes(bury_part)  # [12, 12, 12]
```
`partition_graph` requires a nonempty simple undirected graph and a positive number of hamlets that divides the number of vertices, so each hamlet has the same number of villagers.

Finally, we can visualize the inter-hamlet connections needed to generate a graph state on the configuration specified by our `GraphPartition`. A few ways of evaluating the cost of the partitioning are implemented in QuantumHamlets.jl:

- `edge_cut_cost(g, p)` counts the total number of cut edges.
- `matching_cost(g, p)` sums the maximum matching sizes between each unordered pair of hamlets.
- `pairwise_cut_rank(g, p)` sums the GF(2) ranks of the biadjacency matrices between each unordered pair of hamlets.
- `vertex_cover_costs(g, p)` gives two greedy vertex-cover estimates; see the [usage guide](docs/usage.md#numerical-costs) for their meanings.

`plot_partition` takes a graph, its partition, and a choice of edges to highlight, and returns a figure and the corresponding cost. We can start by highlighting all cut edges:
```julia
f_bury, bury_cost = plot_partition(g, bury_part; highlight=:cut_edges)
f_bury
```
![BURY partitioned 6×6 grid with cut edges highlighted](./resources/figs/bury_naive_6_6_grid.png)
The `bury_cost` variable here is equal to the number of highlighted edges, which is 15 in this case. `edge_cut_cost(g, bury_part)` gives the same count of cut edges, a metric that our paper argues is quite suboptimal for *graph state* partitioning. Furthermore, BURY is designed to perform well on minimizing the matching sizes that cross the cuts. Here is how we visualize that:
```julia
f_bury, bury_cost = plot_partition(g, bury_part; highlight=:matching)
f_bury
```
![BURY partitioned 6×6 grid with matching edges highlighted](./resources/figs/bury_vcg_6_6_grid.png)
Now `bury_cost` is equal to 9, and as we describe in the paper, this corresponds to being able to couple villager qubits to their mayors and perform $Y$ measurements in order to generate the same desired graph state. The cost here is the number of entanglements needed between hamlets (QPUs) under that preparation scheme. This cost is also equal to the sum of the sizes of the maximum matchings between each pair of hamlets.

We can calculate the costs directly, without plotting:
```julia
edge_cut_cost(g, bury_part)     # 15
matching_cost(g, bury_part)     # 9
pairwise_cut_rank(g, bury_part) # 9
```
Matching cost and GF(2) rank happen to agree here, but they are different objectives and need not agree in general. With more than two hamlets, `pairwise_cut_rank` is a sum of pairwise ranks, not a single bipartite cut rank. Rank has no corresponding set of highlighted edges; use `plot_partition(g, bury_part; highlight=:none)` to show membership and calculate rank separately. The highlighted matching edges may also differ between runs while giving the same matching cost.

For comparison, out of the existing methods we considered, we found METIS to be most competitive. Here is how this graph compares against METIS:
```julia
metis_part = metis_partition(g, numHamlets; seed=0)
partition_sizes(metis_part)
edge_cut_cost(g, metis_part)
matching_cost(g, metis_part)
pairwise_cut_rank(g, metis_part)
f_metis, metis_cost = plot_partition(g, metis_part; highlight=:matching)
f_metis
```
![METIS partitioned 6×6 grid with matching edges highlighted](./resources/figs/metis_vcg_6_6_grid.png)
The original figure above has matching cost 10. In the current environment, the call with `seed=0` also gives matching cost 10, edge-cut cost 10, and pairwise rank 10, with three groups of 12 vertices. Notice that METIS does better on cut edges (10 versus 15), while BURY does better at reducing matching sizes across the cuts (9 versus 10) in this example. METIS results can vary with the version and seed, so calculate the costs and check the partition sizes for your own run. Native METIS allows imbalance; use `metis_partition(g, numHamlets; balance=true)` if you need explicit balancing.

## A few more comparisons

For partitioning a graph into two pieces, we also provide a few baselines. The simulated annealing (SA) methods minimize GF(2) cut rank:
```julia
p = min_cutrank(g; seed=3)       # random balanced start, then cut-rank SA
p = bury_sa(g; seed=3)           # ordinary BURY, then SA
p = bury_seeded_sa(g; seed=3)    # sqrt-seeded BURY, then SA
p = metis_partition(g)          # native METIS bipartition, minimizing edge cut
p = metis_sa(g; seed=3)         # explicitly balance METIS, then SA
```
All of these calls return a `GraphPartition`. Sqrt-seeded BURY samples `ceil(sqrt(n))` seeds with replacement for each target size, selects the inner seeds by matching cost, and compares the floor/ceiling target sizes by rank. SA retains its best observed partition, including its initializer, by default. Native METIS allows imbalance, while `metis_sa` explicitly balances its initializer by default. See the [usage guide](docs/usage.md#simulated-annealing-and-metis-baselines) for options and further examples.
