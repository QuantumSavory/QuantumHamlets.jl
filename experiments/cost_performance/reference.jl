# Frozen post-bug-fix / pre-optimization evaluator snapshot.
module ReferenceCosts
import QuantumHamlets
using QuantumHamlets: quantumLand, inverse_map
using QuantumHamlets.Graphs
using QuantumHamlets.Nemo
using QuantumHamlets.IGraphs
include("reference_costs.jl")

# Preserve the original scalar API's validation/conversion and edge-color work.
function _land_for(g, p)
    nv(g) == length(p) || throw(ArgumentError("partition length does not match graph"))
    return QuantumHamlets.quantumLand(p)
end
edge_cut_cost(g, p) = first(naive_cost_of_partition(_land_for(g, p), g))
matching_cost(g, p) = first(matching_cost_of_partition(_land_for(g, p), g))
pairwise_cut_rank(g, p) = first(cut_rank_cost(_land_for(g, p), g))
function vertex_cover_costs(g, p)
    free, paid = vertex_cover_cost_of_partition(_land_for(g, p), g)
    return (ghz_free=free, ghz_n_minus_one=paid)
end
end
