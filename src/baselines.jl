"""Randomly partitions the graph into k pieces\n
Input: a graph from Graphs.jl, and integer k denoting the number of pieces\n
Output: a registry (dictionary mapping villager ids to hamlet ids)"""
function k_partition_random(g::Graphs.Graph, k)
    nodes = shuffle(collect(1:nv(g)))
    new_reg = Dict()
    cap = Int(nv(g)/k) # Assuming the graph can be equally partitioned
    for i in 1:k
        for j in 1:cap
            new_reg[nodes[(i-1)*cap+j]] = i
        end
    end
    return new_reg
end

"""Implementation of the n/2 balanced partitioning algorithm from https://ics.uci.edu/~vazirani/k-cut.pdf \n
Input: a graph and number of pices to partition into.\n
Output: `(registry, trace)`, with component snapshots and balance checks in the trace.
This experimental implementation can return an incomplete registry; see the usage guide.
"""
function k_partition_saran_vazirani(g::Graphs.Graph, k)
    # Takes an input to multipartion S, and indices returned by the optimizer, inds. 
    # This function returns a bool stating whether the partiion is perfectly balanced or not.
    function partitionable(S, inds, k)
        ans = true
        val = sum(S[inds .== 1])
        for i in 2:k
            val_prime = sum(S[inds .== i])
            if val_prime != val
                ans = false
                break
            end
        end
        return ans
    end
    # Given a graph g, and two vectors corresponding the vertices within that graph, returns the set of edges the crosses the partition
    function edges_spanning_partition(g::Graphs.Graph, v1, v2)
        if length(v1) > length(v2)
            v3 = v2
            v2 = v1
            v1 = v3
        end

        edges = []
        for v in v1
            neighbors = Graphs.neighbors(g,v)
            for n in neighbors
                if n ∉ v1
                    push!(edges,Edge(v,n))
                end
            end
        end

        return edges
    end

    d = Graphs.DiGraph(g)
    capacity_matrix = zeros(Int, nv(g), nv(g))
    for e in edges(d)
        u, v = src(e), dst(e)
        capacity_matrix[u, v] = 1
    end

    min_cuts = []
    for e in edges(g)
        cut = GraphsFlows.mincut(d,src(e),dst(e),capacity_matrix,GraphsFlows.BoykovKolmogorovAlgorithm())
        push!(min_cuts, cut)
    end

    unique!(min_cuts)
    sort!(min_cuts, by = x-> x[3], rev=false)

    edge_set = []
    for c in min_cuts
        push!(edge_set, edges_spanning_partition(g,c[1],c[2]))
    end

    dump = []
    g_prime = copy(g)
    new_reg = Dict()
    for cut in edge_set
        for e in cut
            rem_edge!(g_prime, e)
        end
    
        components = Graphs.connected_components(g_prime)
        S = [length(x) for x in components]
        inds = partition(S, k; optimizer = HiGHS.Optimizer);
        ans = partitionable(S, inds, k)
        push!(dump, (copy(components), ans))

        
        if ans == true
            # Turn the found solution into a dictionary/ vector of vectors (i.e. [[Nodes 1 4 5], [Nodes 2 3 6]])
            for i in 1:k
                for j in S[inds .== i]
                    index = findfirst(length.(components).==j)
                    chosen_nodes = components[index]
                    deleteat!(components, index)
                    for node in chosen_nodes
                        new_reg[node] = i
                    end
                end
            end 

            break
        end
    end

    return new_reg, dump
end
