"""Validate the common nonempty, exactly balanced partitioning boundary."""
function _check_partition_count(g, k)
    k isa Integer && 1 <= k <= nv(g) ||
        throw(ArgumentError("k must be an integer satisfying 1 <= k <= nv(g)"))
    nv(g) % k == 0 ||
        throw(ArgumentError("k must divide the number of vertices"))
    return Int(k)
end

"""Helper function to some of the other methods"""
function inverse_map(map, size=maximum(map))
    inv_map = [0 for _ in 1:size]
    for (index, e) in enumerate(map)
        inv_map[e] = index
    end
    return inv_map
end

"""Takes a land and a graph, and returns a graph with the intra-hamlet edges removed"""
function remove_local_edges(land::quantumLand, graph::Graphs.Graph)
    graph_copy = copy(graph)
    for e in edges(graph)
        if land.registry[src(e)] == land.registry[dst(e)]
            rem_edge!(graph_copy, e)
        end
    end
    return graph_copy
end
