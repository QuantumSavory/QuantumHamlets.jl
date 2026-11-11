"""Calculates the edges cut to establish a graph state on the given quantumLand\n
Input: a quantum land object which specifies how nodes are distributed, and a graph from Graphs.jl\n
Output: First the number edges cut, then a vector of colors (cut edges are set to red and others black). 
The vector of colors is used by "visualize_graph_on_land" to illustrate the partition and cost.
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

"""Returns first the size of the vertex cover, corresponding to GHZ states being 100% free.
Second returns the cost if an m partite GHZ states has cost m-1."""
function vertex_cover_cost_of_partition(land::quantumLand, g::Graphs.Graph)
    k = land.numHamlets
    local_removed_g = QuantumHamlets.remove_local_edges(land, g)
    cover = vertex_cover(local_removed_g, DegreeVertexCover())
    ghz_free_cost = length(cover)

    ghz_n_minus_one_cost = 0
    partite_counter = [false for _ in 1:k]
    for v in cover
        partite_counter[land.registry[v]] = true
        for n in neighbors(local_removed_g, v)
            partite_counter[land.registry[n]] = true
        end
        ghz_n_minus_one_cost += (sum(partite_counter) - 1)
        partite_counter = [false for _ in 1:k]
    end
    return ghz_free_cost, ghz_n_minus_one_cost
end

"""Calculates the cost of a k-partition as a sum of pairwise GF(2) cut ranks.\n
Input: a quantum land object, a graph from Graphs.jl\n
Output: `(cost, nothing)`, where cost is an integer. For a scalar result use `pairwise_cut_rank(g, p)`.
"""
function cut_rank_cost(land::quantumLand, g::Graphs.Graph)
    k = land.numHamlets
    local_removed_g = QuantumHamlets.remove_local_edges(land, g)

    bipartite_decomposition = []
    maps = []
    for i in 1:k
        for j in (i+1):k
            subgraph_vertices = Array{Int64}([])
            for v in vertices(local_removed_g)
                if land.registry[v] == i || land.registry[v] == j
                    push!(subgraph_vertices,v)
                end
            end
            g_bi, map = induced_subgraph(local_removed_g, subgraph_vertices)
            push!(bipartite_decomposition, g_bi)
            push!(maps, map)
        end
    end

    R, = residue_ring(ZZ, 2)
    sum_ranks = 0
    for graph in bipartite_decomposition
        sum_ranks += rank(matrix(R, adjacency_matrix(graph)))
    end    

    return sum_ranks÷2, nothing
end

"""Calculates the cost of generating a given graph state by only using bell pairs and Y measurements.\n
Input: a quantum land object, a graph from Graphs.jl\n
Output: First the sum of the matching sizes between each pair of partitions,
then a vector of colors (selected matching edges are red and others black). 
The vector of colors is used by "visualize_graph_on_land" to illustrate the partition and cost."""
function matching_cost_of_partition(land::quantumLand, g::Graphs.Graph)
    k = land.numHamlets
    local_removed_g = QuantumHamlets.remove_local_edges(land, g)

    bipartite_decomposition = []
    maps = []
    for i in 1:k
        for j in (i+1):k
            subgraph_vertices = Array{Int64}([])
            for v in vertices(local_removed_g)
                if land.registry[v] == i || land.registry[v] == j
                    push!(subgraph_vertices,v)
                end
            end
            g_bi, map = induced_subgraph(local_removed_g, subgraph_vertices)
            push!(bipartite_decomposition, g_bi)
            push!(maps, map)
        end
    end

    function flip_nodes(v::Int, g::Graphs.Graph, visted::Vector{Bool}, vtypes::Vector{Bool})
        for neighbor in neighbors(g,v)
            if !visted[neighbor]
                visted[neighbor] = true
                vtypes[neighbor] = !vtypes[v]
                flip_nodes(neighbor, g, visted, vtypes)
            end
        end
    end

    edgecolors = [:black for _ in 1:ne(g)]
    cost_all = []
    matchings = []
    for (index, bi_g) in enumerate(bipartite_decomposition)
        cost = 0

        vtypes = [false for _ in 1:nv(bi_g)]
        visited = [false for _ in 1:nv(bi_g)]
        while sum(visited) < length(visited)
            node = findfirst(==(false),visited)  
            visited[node] = true      
            flip_nodes(node, bi_g, visited, vtypes)
        end

        vtypes = IGraphs.IGVectorBool(vtypes)
        matching = IGraphs.IGVectorInt()
        weights = IGraphs.IGVectorFloat(ones(Graphs.ne(bi_g)))
        ig = IGraphs.IGraph(bi_g)
        IGraphs.LibIGraph.maximum_bipartite_matching(ig, vtypes, matching, weights, 1e-8)

        matching .+= 1
       
        inv_map = inverse_map(maps[index],nv(g))

        for (index, e) in enumerate(edges(g))
            if inv_map[src(e)]!= 0 && inv_map[dst(e)] != 0
                if matching[inv_map[src(e)]] == inv_map[dst(e)] 
                    cost += 1
                    edgecolors[index] = :red
                end
            end
        end
        push!(cost_all, cost)
        push!(matchings, matching)
    end

    return sum(cost_all), edgecolors
end
