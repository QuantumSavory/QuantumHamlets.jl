module QuantumHamlets
using CairoMakie
using Graphs
using GraphsFlows
using GraphMakie
using Colors
using MultiwayNumberPartitioning, HiGHS
using Random
using IGraphs
using Makie
using Nemo
import Metis
import PythonCall

export GraphPartition, partition_graph, bipartition_graph, membership, registry, partition_sizes,
       edge_cut_cost, matching_cost, pairwise_cut_rank, vertex_cover_costs,
       matching_details, plot_partition,
       min_cutrank, bury_sa, bury_seeded_sa, metis_partition, metis_sa

# Existing qualified entry points remain available with their original contracts.
public quantumLand, villager, hamlet, defaultVillagerRegistry,
       bury, bury_k, bury_seeded, bury_k_seeded, bury_not_balanced_bipartition,
       best_seed, best_sqrt_seed,
       k_partition_random, k_partition_saran_vazirani,
       naive_cost_of_partition, matching_cost_of_partition, cut_rank_cost,
       vertex_cover_cost_of_partition, print_graph, visualize_graph_on_land

include("quantum_land.jl")
include("graph_utils.jl")
include("costs.jl")
include("bury.jl")
include("visualization.jl")
include("display.jl")
include("partition.jl")
include("api.jl")
include("baselines.jl")

end # module QuantumHamlets
