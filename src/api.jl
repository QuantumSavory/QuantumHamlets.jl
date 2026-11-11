# These adapters expose partition values and forward BURY backend selection.
# Cost routines retain the legacy tuple formats behind scalar adapters.

function _labels_for(g::Graphs.Graph, p::GraphPartition)
    nv(g) == length(p) ||
        throw(ArgumentError("partition has $(length(p)) vertices, but graph has $(nv(g))"))
    return p.labels
end

function _land_for(g::Graphs.Graph, p::GraphPartition)
    _labels_for(g, p)
    return quantumLand(p)
end

"""
    partition_graph(g, k; algorithm=:bury, backend=:heap, verbose=false) -> GraphPartition

Partition a nonempty simple undirected graph into `k` equally sized parts.
Require `1 <= k <= nv(g)` and `nv(g) % k == 0`.

Algorithms delegate to the existing implementations:
- `:bury` calls [`bury_k`](@ref).
- `:bury_seeded` calls [`bury_k_seeded`](@ref), with its deterministic
  minimum-weight seeder for each part.
- `:random` calls [`k_partition_random`](@ref), using the default Random RNG.
- `:saran_vazirani` calls [`k_partition_saran_vazirani`](@ref); its trace is
  available through that legacy function instead.

BURY accepts `backend=:heap` (default), `:scan`, or `:legacy`, forwarded through
every partition round. Seed and partition policies are unchanged. `verbose`
is forwarded to BURY. Other algorithms only accept the default backend value,
which does not apply to them, and retain their existing output behavior.
An incomplete registry raises an error. The Saran–Vazirani baseline first
tries to pack the original components and reports failure explicitly.
"""
function partition_graph(g::Graphs.Graph, k;
                         algorithm::Symbol=:bury, backend::Symbol=:heap,
                         verbose::Bool=false)
    k = _check_partition_count(g, k)
    _bury_check_backend(backend)
    if algorithm in (:random, :saran_vazirani) && backend !== :heap
        throw(ArgumentError("backend selection applies only to :bury and :bury_seeded"))
    end

    reg = if algorithm === :bury
        bury_k(g, k; backend=backend, verbose=verbose)
    elseif algorithm === :bury_seeded
        bury_k_seeded(g, k; backend=backend, verbose=verbose)
    elseif algorithm === :random
        k_partition_random(g, k)
    elseif algorithm === :saran_vazirani
        first(k_partition_saran_vazirani(g, k))
    else
        throw(ArgumentError("unknown algorithm $algorithm; use :bury, :bury_seeded, :random, or :saran_vazirani"))
    end

    length(reg) == nv(g) ||
        error("$algorithm returned an incomplete registry")
    return GraphPartition(reg; nparts=k)
end

"""
    edge_cut_cost(g, p::GraphPartition) -> Int

Count edges with endpoints in different parts. This is the scalar form of
[`naive_cost_of_partition`](@ref).
"""
edge_cut_cost(g::Graphs.Graph, p::GraphPartition) =
    _edge_cut_cost(_labels_for(g, p), g)

"""
    matching_cost(g, p::GraphPartition) -> Int

Sum maximum matching sizes over unordered pairs of parts. One declared part
has cost zero. Uses the same unweighted matching kernel as
[`matching_cost_of_partition`](@ref), without constructing edge colors. Use
[`matching_details`](@ref) when the matching edges are also needed.
"""
matching_cost(g::Graphs.Graph, p::GraphPartition) =
    _matching_cost(_labels_for(g, p), g)

"""
    pairwise_cut_rank(g, p::GraphPartition) -> Int

Sum the GF(2) ranks of the biadjacency matrices between unordered pairs of
parts. This is the scalar form of [`cut_rank_cost`](@ref). For more than two
parts this is a sum of pairwise ranks, not a single bipartite cut rank.
"""
pairwise_cut_rank(g::Graphs.Graph, p::GraphPartition) =
    _pairwise_cut_rank(_labels_for(g, p), p.nparts, g)

"""
    vertex_cover_costs(g, p::GraphPartition)

Return `(ghz_free=..., ghz_n_minus_one=...)`, the two estimates from
[`vertex_cover_cost_of_partition`](@ref). The underlying cover is greedy,
not necessarily minimum. The second model charges `m-1` for an `m`-partite
GHZ state. Both estimates retain the original implementation's semantics.
"""
function vertex_cover_costs(g::Graphs.Graph, p::GraphPartition)
    free, paid = _vertex_cover_costs(_labels_for(g, p), p.nparts, g)
    return (ghz_free=free, ghz_n_minus_one=paid)
end

"""
    matching_details(g, p::GraphPartition)

Return `(cost=..., edges=...)` from one evaluation of the matching
routine. Edges use the original graph's vertex IDs and follow `edges(g)`
order. They form a maximum matching within each pair of parts; their union
need not be a single matching across the whole graph.
"""
function matching_details(g::Graphs.Graph, p::GraphPartition)
    labels = _labels_for(g, p)
    colors = fill(:black, ne(g))
    cost = _matching_cost(labels, g, colors)
    selected = [e for (e, color) in zip(edges(g), colors) if color === :red]
    return (cost=cost, edges=selected)
end

"""
    plot_partition(g, p::GraphPartition; highlight=:matching)

Return `(figure, cost)` using the existing plotting implementation.
`highlight=:matching` highlights pairwise matching edges; `:cut_edges`
highlights every cross-part edge. With `:none`, show membership without
highlighted edges and return `nothing` for the cost.

Plotting dependencies still load with the package. To display rank, evaluate
[`pairwise_cut_rank`](@ref) separately and use `highlight=:none`.
"""
function plot_partition(g::Graphs.Graph, p::GraphPartition; highlight::Symbol=:matching)
    land = _land_for(g, p)
    if highlight === :matching
        return visualize_graph_on_land(land, g; method=matching_cost_of_partition)
    elseif highlight === :cut_edges
        return visualize_graph_on_land(land, g; method=naive_cost_of_partition)
    elseif highlight === :none
        return visualize_graph_on_land(land, g; method=(_, _) -> (nothing, nothing))
    else
        throw(ArgumentError("unknown highlight $highlight; use :matching, :cut_edges, or :none"))
    end
end
