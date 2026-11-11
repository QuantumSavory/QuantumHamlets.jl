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
