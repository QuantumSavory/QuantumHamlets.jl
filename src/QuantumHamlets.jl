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

"""Struct to hold villager qubits. Their address should be the hamlet their in, and their id is unique identifier amongst all villagers"""
mutable struct villager
    address::Int 
    id::Int
end

mutable struct hamlet
    numVillagers::Int
    numMayors::Int
    numFarms::Int
    villagers::Vector{villager}
    id::Int
end
function hamlet(hamlet_ID::Int)
    return hamlet(0,1,0,villager[],hamlet_ID)
end

function add_villager!(h::hamlet, p::villager)
    push!(h.villagers,p)
    h.numVillagers+=1
    return h
end

# Default map from villagers to their hamlets
function defaultVillagerRegistry(numHamlets::Int,hamletSize::Int)
    villagerLookup = Dict()
    currentCapacity = 0
    currentHamlet = 1
    for p in 1:numHamlets*hamletSize
        villagerLookup[p] = currentHamlet
        currentCapacity += 1
        if currentCapacity == hamletSize
            currentCapacity = 0
            currentHamlet += 1
        end
    end
    return villagerLookup
end

struct quantumLand
    hamlets::Vector{hamlet}
    registry::Dict{Int,Int}
    numHamlets::Int
    numVillagers::Int
end
function quantumLand(citizenRegistry::Dict, numHamlets::Int, numVillagers::Int)
    land = quantumLand([hamlet(i) for i in 1:numHamlets],citizenRegistry, numHamlets, numVillagers)

    for i in 1:numVillagers
        hamlet_id = citizenRegistry[i]
        add_villager!(land.hamlets[hamlet_id],villager(hamlet_id,i))
    end

    return land
end

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

"""Calculates the the cost of a k-partition in terms of cut rank.\n
Input: a quantum land object, a graph from Graphs.jl\n
Output: cut rank returned as an integer
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

    return sum_ranks÷2
end

"""Calculates the cost of generating a given graph state by only using bell pairs and Y measurements.\n
Input: a quantum land object, a graph from Graphs.jl\n
Output: First the sum of the matching sizes between each pair of partitions,
then a vector of colors (cut edges are set to red and others black). 
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

"""
Implementation of k-partition BURY algorithm from our work, QuantumHamlets: https://arxiv.org/abs/2603.06387
This function assumes that k divides the number of vertices in g. \n
Input: a graph from Graphs.jl, and integer k denoting the number of pieces\n
Output: a registry (dictionary mapping villager ids to hamlet ids)
"""
function bury_k(g::Graphs.Graph, k; verbose=false)
    @assert nv(g)%k ==0

    current_graph = copy(g)
    new_reg = Dict()
    inv_map = collect(1:nv(g))
    total_colored = [false for _ in 1:nv(g)]
    for i in 1:k
        _, colored = bury(current_graph, nv(g)÷k, verbose=verbose)

        for v in 1:nv(g)
            v_in_subgraph = inv_map[v]
            if v_in_subgraph != 0
                if colored[v_in_subgraph]
                    new_reg[v] = i
                    total_colored[v] = true
                end
            end
        end

        current_graph, map = induced_subgraph(g, findall(==(false),total_colored))
        inv_map = inverse_map(map, nv(g))
    end

    return new_reg
end

"""
Implementation of BURY algorithm from our work, QuantumHamlets: https://arxiv.org/abs/2603.06387
This version colors a number of nodes in the graph equal to the second parameter in the input, and has k set to 2.
To use BURY for k partition, use the "bury_k" function in this library.\n
Input: A graph to be bipartitioned, and a number of nodes that will be in the first partition \n
Output: A registry (see README example for explanation), and a Boolean vector marking which nodes where colored
(k=2 case has only colored/not colored in the bipartition). For most users, the registry is all you need. 
The vector is primarily metadata for "bury_k".
"""
function bury(g::Graphs.Graph, r0=nv(g)÷2; verbose=false)
    r = r0

    # Node weight initialization
    weights = [0.0 for _ in 1:nv(g)]
    colored = [false for _ in 1:nv(g)]
    for (index,v) in enumerate(vertices(g))
        weights[index] = Graphs.degree(g,v) + 1
    end

    if verbose
        println("Initial weights and colors:")
        println("\nWeights table")
        println(weights)
        println("Colored table")
        println(colored)
    end

    # Global greedy selection
    while sum(colored)<r0
        cost, choice = findmin(weights)

        if verbose
            println("\n\n###########################\nStarting round. r=", r)
            println("Chose vertex ", choice, " for the price of ", cost)
        end

        if r>= cost
            nbd = vcat(choice,neighbors(g,choice))
            r -= cost
        else
            nbd = choice
            r -= 1
        end

        for v in nbd
            if !colored[v]
                weights[v] -= 1
                weights[neighbors(g,v)] .-= 1
                colored[v] = true
            end
        end
        weights[choice] = Inf
        if verbose
            println("\nUpdated weights table")
            println(weights)
            println("Updated Colored table")
            println(colored)
        end
    end

    new_reg = Dict()
    for i in 1:nv(g)
        new_reg[i] = colored[i]+1
    end
    return new_reg, colored
end

# Should maybe add a check that the seed can be chosen - not a problem with this minweight seeder
function seeder_1(g::Graphs.Graph)
    # Node weight initialization
    weights = [0.0 for _ in 1:nv(g)]
    for (index,v) in enumerate(vertices(g))
        weights[index] = Graphs.degree(g,v) + 1
    end
    _, index = findmin(weights)
    #println("seed is", index)
    return index
end

function best_seed(g::Graphs.Graph, r0)
    best_score = Inf
    best_seed = 1
    for seed in 1:nv(g)
        reg, _ = bury_seeded(g, r0, seed)
        land = quantumLand(reg, 2, nv(g))
        cost, _ = matching_cost_of_partition(land, g)

        if cost< best_score
            best_seed = seed
            best_score = cost
        end
        #println(cost, " seed ", seed)
    end
    #println("using seed ", best_seed, "\n")
    return best_seed
end

# this version should work for any k
function bury_k_seeded(g::Graphs.Graph, k; verbose=false)
    @assert nv(g)%k ==0

    current_graph = copy(g)
    new_reg = Dict()
    inv_map = collect(1:nv(g))
    total_colored = [false for _ in 1:nv(g)]
    for i in 1:k
        returned_reg, colored = bury_seeded(current_graph, nv(g)÷k, seeder_1(current_graph), verbose=verbose)
        #returned_reg, colored = bury_seeded(current_graph, nv(g)÷k, verbose=verbose)

        for v in 1:nv(g)
            v_in_subgraph = inv_map[v]
            if v_in_subgraph != 0
                if colored[v_in_subgraph]
                    new_reg[v] = i
                    total_colored[v] = true
                end
            end
        end

        current_graph, map = induced_subgraph(g, findall(==(false),total_colored))
        inv_map = inverse_map(map, nv(g))
    end

    return new_reg
end
function bury_seeded(g::Graphs.Graph, r0=nv(g)÷2, seed = best_seed(g, r0); verbose=false)
    r = r0

    # Node weight initialization
    weights = [0.0 for _ in 1:nv(g)]
    colored = [false for _ in 1:nv(g)]

    region = [false for _ in 1:nv(g)]
    region[neighbors(g,seed)] .= true
    region[seed] = true

    for (index,v) in enumerate(vertices(g))
        weights[index] = Graphs.degree(g,v) + 1
    end

    if verbose
        println("Initial weights and colors:")
        println("\nWeights table")
        println(weights)
        println("Colored table")
        println(colored)
    end

    # Global greedy selections
    while sum(colored)<r0
        cost, choice = findmin(.!region*Inf + weights)
        if cost == Inf
            cost, choice = findmin(weights)
        end

        if verbose
            println("\n\n###########################\nStarting round. r=", r)
            println("Chose vertex ", choice, " for the price of ", cost)
        end

        if r>= cost
            nbd = vcat(choice,neighbors(g,choice))
            r -= cost
        else
            nbd = choice
            r -= 1
        end

        for v in nbd
            if !colored[v]
                weights[v] -= 1
                weights[neighbors(g,v)] .-= 1
                colored[v] = true
            end
            region[v] = true
        end
        weights[choice] = Inf

        if verbose
            println("\nUpdated weights table")
            println(weights)
            println("Updated Colored table")
            println(colored)
        end
    end

    new_reg = Dict()
    for i in 1:nv(g)
        new_reg[i] = colored[i]+1
    end
    return new_reg, colored
end

"""Implementation of the n/2 balanced partitioning algorithm from https://ics.uci.edu/~vazirani/k-cut.pdf \n
Input: a graph and number of pices to partition into.\n
Output:
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

######################## Visualization functions ######################## 
function print_graph(g::Graph; edgecolors=nothing, nodecolors=[:white for _ in 1:nv(g)], nodelabels = true)
    labels = collect(1:nv(g))
    if nodelabels == false
        labels = nothing
    end
    thicken_lines = edgecolors.==:red
    f, ax, p = graphplot(g,
                        node_color=nodecolors,
                        edge_color = edgecolors,
                        edge_width= [2 for _ in 1:ne(g)] .+ 2*thicken_lines,
                        nlabels_align=(:center,:center);
                        ilabels= labels, 
                        size = (900,700))
    hidedecorations!(ax); hidespines!(ax); ax.aspect = DataAspect()

    return f
end

"""Returns a graph with nodes colored to match the hamlet they belong to. Edges are colored red
to indicate cost. Cost is default any cross hamlet edge. Using a different method can color edges differently.
For now, only the naive_cost_of_partition() method is supported."""
function visualize_graph_on_land(land::quantumLand, g::Graphs.Graph; method=QuantumHamlets.naive_cost_of_partition)
    cost, edgecolors = method(land,g)
    hamlet_colors = Makie.wong_colors()

    nodecolors = [hamlet_colors[(land.registry[i]%7)+1] for i in 1:land.numVillagers]
    f = print_graph(g, edgecolors=edgecolors, nodecolors=nodecolors)
    return f, cost
end

function Base.show(io::IO, land::quantumLand)
    for i in land.hamlets
        println(i)
    end
end
function Base.show(io::IO, p::villager)
    print("Qubit $(p.id), who lives in hamlet $(p.address)")
end
function Base.show(io::IO, hamlet::hamlet)
    print("Village $(hamlet.id) has $(hamlet.numVillagers) villagers, $(hamlet.numFarms) farms, and $(hamlet.numMayors) mayors.\n")
    for p in hamlet.villagers
        println(p)
    end
end
end # module QuantumHamlet
