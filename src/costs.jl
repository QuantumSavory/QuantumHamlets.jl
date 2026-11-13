# Cost kernels work on dense labels. Legacy wrappers copy dictionary values
# once; the GraphPartition scalar API can use its already-dense membership.
_cost_labels(land::quantumLand, g) = Int[land.registry[v] for v in vertices(g)]

function _edge_cut_cost(labels, g, edgecolors=nothing)
    cost = 0
    for (index, e) in enumerate(edges(g))
        if labels[src(e)] != labels[dst(e)]
            cost += 1
            isnothing(edgecolors) || (edgecolors[index] = :red)
        end
    end
    return cost
end

"""Return `(cut_edge_count, edge_colors)` for a legacy partition.
Cross-partition edges are red and all other edges are black.
"""
function naive_cost_of_partition(land::quantumLand, graph::Graphs.Graph) 
    edgecolors = [:black for _ in 1:ne(graph)]

    cost = 0
    for (index, e) in enumerate(edges(graph))
        if land.registry[src(e)] != land.registry[dst(e)]
            cost += 1
            edgecolors[index] = :red
        end
    end

    return cost, edgecolors
end

function _vertex_cover_costs(labels, k, g)
    # Construct the same sorted adjacency as copy-and-remove, without copying
    # internal edges. Keep Graphs' greedy algorithm and its tie order intact.
    cut_graph = SimpleGraph{eltype(g)}(nv(g))
    for e in edges(g)
        labels[src(e)] == labels[dst(e)] || add_edge!(cut_graph, e)
    end
    cover = vertex_cover(cut_graph, DegreeVertexCover())
    ghz_n_minus_one_cost = 0
    seen = zeros(Int, k)
    for (stamp, v) in enumerate(cover)
        seen[labels[v]] = stamp
        for u in neighbors(cut_graph, v)
            label = labels[u]
            if seen[label] != stamp
                seen[label] = stamp
                ghz_n_minus_one_cost += 1
            end
        end
    end
    return length(cover), ghz_n_minus_one_cost
end

"""Return the existing greedy vertex-cover estimates `(ghz_free, ghz_n_minus_one)`.
The second model charges `m-1` for an `m`-partite GHZ state. These are greedy
estimates, not minimum vertex-cover costs.
"""
vertex_cover_cost_of_partition(land::quantumLand, g::Graphs.Graph) =
    _vertex_cover_costs(_cost_labels(land, g), land.numHamlets, g)

# Store only nonempty pairs. Orient endpoints by part label and retain each
# original edge index for optional plotting colors. Each graph edge is read once.
function _cut_edge_groups(labels, g)
    groups = Dict{Tuple{Int,Int},Vector{NTuple{3,Int}}}()
    for (index, e) in enumerate(edges(g))
        u, v = Int(src(e)), Int(dst(e))
        a, b = labels[u], labels[v]
        a == b && continue
        if a > b
            a, b = b, a
            u, v = v, u
        end
        push!(get!(Vector{NTuple{3,Int}}, groups, (a, b)), (u, v, index))
    end
    return groups
end

# Each Julia matrix column stores one packed binary row, keeping wordwise XOR
# contiguous in memory. Only this private work matrix is modified.
function _gf2_packed_rank!(rows::Matrix{UInt64}, ncols::Int)
    nwords, nrows = size(rows)
    pivots = zeros(Int, ncols)
    result = 0
    @inbounds for row in 1:nrows
        word = 1
        while word <= nwords
            if iszero(rows[word, row])
                word += 1
                continue
            end
            col = 64 * (word - 1) + trailing_zeros(rows[word, row]) + 1
            pivot = pivots[col]
            if iszero(pivot)
                pivots[col] = row
                result += 1
                result == ncols && return result
                break
            end
            for w in word:nwords
                rows[w, row] ⊻= rows[w, pivot]
            end
        end
    end
    return result
end

function _pairwise_cut_rank(labels::AbstractVector{<:Integer}, nparts::Int,
                            g::Graphs.Graph)
    nparts <= 1 && return 0
    pairs = _cut_edge_groups(labels, g)
    isempty(pairs) && return 0
    positions = zeros(Int, nv(g))
    touched = Int[]
    result = 0
    for pair_edges in values(pairs)
        nleft = nright = 0
        empty!(touched)
        # Omit isolated vertices of both parts: they contribute only zero rows
        # or columns, which cannot affect rank.
        for (u, v, _) in pair_edges
            if iszero(positions[u])
                nleft += 1
                positions[u] = nleft
                push!(touched, u)
            end
            if iszero(positions[v])
                nright += 1
                positions[v] = nright
                push!(touched, v)
            end
        end
        # Rank is invariant under transpose. Pack the smaller side into bits.
        ncols, nrows = minmax(nleft, nright)
        packed = zeros(UInt64, cld(ncols, 64), nrows)
        for (u, v, _) in pair_edges
            row, col = nleft <= nright ? (positions[v], positions[u] - 1) :
                                        (positions[u], positions[v] - 1)
            packed[(col >>> 6) + 1, row] |= UInt64(1) << (col & 63)
        end
        result += _gf2_packed_rank!(packed, ncols)
        for v in touched
            positions[v] = 0
        end
    end
    return result
end

"""Return `(pairwise_cut_rank, nothing)` for a legacy partition.
Compute exact GF(2) ranks of the cross-part biadjacency matrices, using packed
binary elimination. Empty parts contribute zero.
"""
cut_rank_cost(land::quantumLand, g::Graphs.Graph) =
    (_pairwise_cut_rank(_cost_labels(land, g), land.numHamlets, g), nothing)

function _matching_cost(labels, g, edgecolors=nothing)
    groups = _cut_edge_groups(labels, g)
    isempty(groups) && return 0
    local_ids = zeros(Int, nv(g))
    total_cost = 0
    for pair_edges in values(groups)
        # Isolated vertices cannot affect a matching. Number only incident
        # vertices, avoiding a graph copy and a full vertex scan for every pair.
        types = Bool[]
        endpoints = Vector{Int}(undef, 2length(pair_edges))
        for (i, (u, v, _)) in enumerate(pair_edges)
            if local_ids[u] == 0
                push!(types, false)
                local_ids[u] = length(types)
            end
            if local_ids[v] == 0
                push!(types, true)
                local_ids[v] = length(types)
            end
            endpoints[2i - 1] = local_ids[u] - 1
            endpoints[2i] = local_ids[v] - 1
        end
        # IGraph(g) inserts edges one at a time. Bulk construction avoids
        # repeated rebuilding of igraph's edge indices.
        ig = IGraphs.IGraph(_uninitialized=Val(true))
        IGraphs.LibIGraph.create(ig, IGraphs.IGVectorInt(endpoints), length(types), false)
        matching = isnothing(edgecolors) ? IGraphs.IGNull() : IGraphs.IGVectorInt()
        # NULL weights request maximum cardinality directly. Unit weights
        # unnecessarily select the more expensive weighted algorithm.
        cardinality, _ = IGraphs.LibIGraph.maximum_bipartite_matching(
            ig, IGraphs.IGVectorBool(types), matching, IGraphs.IGNull(), 0.0)
        total_cost += cardinality
        if !isnothing(edgecolors)
            for (u, v, index) in pair_edges
                matching[local_ids[u]] == local_ids[v] - 1 && (edgecolors[index] = :red)
            end
        end
        for (u, v, _) in pair_edges
            local_ids[u] = 0
            local_ids[v] = 0
        end
    end
    return total_cost
end

"""Return `(matching_cost, edge_colors)` for a legacy partition.
The cost sums maximum matching cardinalities over unordered pairs of parts.
Selected matching edges are red, all others black. One declared part costs zero.
Alternative maximum matching witnesses may be selected by the unweighted solver.
"""
function matching_cost_of_partition(land::quantumLand, g::Graphs.Graph)
    colors = fill(:black, ne(g))
    return _matching_cost(_cost_labels(land, g), g, colors), colors
end

# Graph-first scalar costs and matching details.
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
