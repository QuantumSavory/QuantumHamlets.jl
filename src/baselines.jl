"""Randomly partitions the graph into k pieces\n
Input: a nonempty graph from Graphs.jl, and integer k in 1:nv(g) dividing nv(g).\n
Output: a registry (dictionary mapping villager ids to hamlet ids)"""
function k_partition_random(g::Graphs.Graph, k)
    k = _check_partition_count(g, k)
    nodes = shuffle(collect(1:nv(g)))
    new_reg = Dict()
    cap = nv(g) ÷ k
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
The first snapshot checks the original connected components before removing edges.
Require a nonempty graph and integer `k` in `1:nv(g)` dividing `nv(g)`.
Solver failure or failure to find a complete balanced assignment raises an error.
"""
function k_partition_saran_vazirani(g::Graphs.Graph, k)
    k = _check_partition_count(g, k)
    dump = []
    function pack_components(graph)
        components = Graphs.connected_components(graph)
        sizes = length.(components)
        inds = try
            partition(sizes, k; optimizer=HiGHS.Optimizer, strict=true)
        catch err
            err isa InterruptException && rethrow()
            error("Saran–Vazirani component packing failed: $(sprint(showerror, err))")
        end
        length(inds) == length(components) && all(i -> i isa Integer && 1 <= i <= k, inds) ||
            error("Saran–Vazirani component packing returned invalid labels")
        balanced = all(i -> sum(sizes[inds .== i]) == nv(g) ÷ k, 1:k)
        push!(dump, (copy(components), balanced))
        balanced || return nothing

        new_reg = Dict()
        for (component, label) in zip(components, inds), node in component
            new_reg[node] = label
        end
        return new_reg
    end

    # A disconnected input may already admit a balanced zero-cut partition.
    initial = pack_components(g)
    initial === nothing || return initial, dump

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

    g_prime = copy(g)
    for cut in edge_set
        for e in cut
            rem_edge!(g_prime, e)
        end
    
        new_reg = pack_components(g_prime)
        new_reg === nothing || return new_reg, dump
    end

    error("Saran–Vazirani did not find a complete balanced partition")
end
