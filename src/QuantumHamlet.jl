module QuantumHamlet
using QuantumClifford
using CairoMakie
using Graphs
using GraphsFlows
using GraphMakie
using Colors
using MultiwayNumberPartitioning, HiGHS
using Random
using IGraphs

greet() = print("(|To be> + |not to be>)/√2")

mutable struct villager
    address::Int 
    id::Int
end

# Maybe add a layout, or a state to represent what is and isn't busy
# Maybe also think about other codes? this could instead hold data on the physical qubits instead
# THis could hold information on the tiles?
mutable struct quantumVillage
    numVillagers::Int
    numMayors::Int
    numFarms::Int
    villagers::Vector{villager}
    id::Int
end
function quantumVillage(村番号::Int)
    return quantumVillage(0,0,0,villager[],村番号)
end

function add_villager!(v::quantumVillage, p::villager)
    push!(v.villagers,p)
    v.numVillagers+=1
    return v
end

# Default map from villagers to their villages
function defaultVillagerRegistry(numVillages::Int,logicalBitsPerVillage::Int)
    villagerLookup = Dict()
    currentCapacity = 0
    currentVillage = 1
    for p in 1:numVillages*logicalBitsPerVillage
        villagerLookup[p] = currentVillage
        currentCapacity += 1
        if currentCapacity == logicalBitsPerVillage
            currentCapacity = 0
            currentVillage += 1
        end
    end
    return villagerLookup
end

struct quantumLand
    villages::Vector{quantumVillage}
    registry::Dict{Int,Int}
    numVillages::Int
    numVillagers::Int
end
function quantumLand(citizenRegistry::Dict, numVillages::Int, numVillagers::Int)
    land = quantumLand([quantumVillage(i) for i in 1:numVillages],citizenRegistry, numVillages, numVillagers)

    for i in 1:numVillagers
        village_id = citizenRegistry[i]
        add_villager!(land.villages[village_id],villager(village_id,i))
    end

    return land
end

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

# TODO this should calculate the cost when in the GHZ model
function vertex_cover_cost_of_partition()
end

"""Calculates the cost of generating a given graph state by only using bell pairs and YY measurement.
Does not look for bicliques"""
function matching_cost_of_partition(land::quantumLand, g::Graphs.Graph)
    k = land.numVillages
    local_removed_g = QuantumHamlet.remove_local_edges(land, g)

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
    # debugging return statement
    #return cost_all, matchings, maps, edgecolors
end

# TODO this function should exist inside matching_cost_of_partition, but its useful to leave out for debugging for now
function inverse_map(map, size=maximum(map))
    inv_map = [0 for _ in 1:size]
    for (index, e) in enumerate(map)
        inv_map[e] = index
    end
    return inv_map
end

"""Takes a land and a graph, and returns a graph with the intra-village edges removed"""
function remove_local_edges(land::quantumLand, graph::Graphs.Graph)
    graph_copy = copy(graph)
    for e in edges(graph)
        if land.registry[src(e)] == land.registry[dst(e)]
            rem_edge!(graph_copy, e)
        end
    end
    return graph_copy
end

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

# TODO currently only supported for k=2
function bury_heuristic_global(g::Graphs.Graph, k=2; verbose=false)
    r = nv(g)÷k

    # Node weight initialization
    weights = [0.0 for _ in 1:nv(g)]
    colored = [false for _ in 1:nv(g)]
    for (index,v) in enumerate(vertices(g))
        weights[index] = degree(g,v) + 1
    end

    if verbose
        println("Initial weights and colors:")
        println("\nWeights table")
        println(weights)
        println("Colored table")
        println(colored)
    end

    # Global greedy selection
    while sum(colored)<nv(g)÷k
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
    return new_reg
end

# TODO currently only supported for k=2
# This version should run much faster than the gloabal version
function bury_heuristic_seeded(g::Graphs.Graph, k, v)
    return
end

# TODO performance on grid graphs is so bad that I think this might not actually be implemented correctly
function k_partition_saran_vazirani(g::Graphs.Graph, k)
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
        push!(edge_set, QuantumHamlet.edges_spanning_partition(g,c[1],c[2]))
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

"""Helper function for Saran-Vazirani. 
Takes an input to multipartion S, and indices returned by the optimizer, inds. 
This function returns a bool stating whether the partiion is perfectly balanced or not."""
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

"""Helper function for Saran-Vazirani. 
Given a graph g, and two vectors corresponding the vertices within that graph, returns the set of edges the crosses the partition"""
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

######################## Visualization functions ######################## 
function print_graph(g::Graph; edgecolors=nothing, nodecolors=[:white for _ in 1:nv(g)])
    f, ax, p = graphplot(g,
                        node_color=nodecolors,
                        edge_color = edgecolors,
                        nlabels_align=(:center,:center);
                        ilabels= collect(1:nv(g)))
    hidedecorations!(ax); hidespines!(ax); ax.aspect = DataAspect()

    return f
end

"""Returns a graph with nodes colored to match the village they belong to. Edges are colored red
to indicate cost. Cost is default any cross village edge. Using a different method can color edges differently.
For now, only the naive_cost_of_partition() method is supported."""
function visualize_graph_on_land(land::quantumLand, g::Graphs.Graph; method=QuantumHamlet.naive_cost_of_partition)
    # This subfunction was generated with ChatGPT
    function distinct_colors(n::Int)
        # Generate `n` distinct colors using HSV space
        return [HSV(i * 360 / n, 0.8, 0.9) for i in 0:n-1] .|> RGB
    end

    cost, edgecolors = method(land,g)

    village_colors = distinct_colors(land.numVillages)
    nodecolors = [village_colors[land.registry[i]] for i in 1:land.numVillagers]
    f = print_graph(g, edgecolors=edgecolors, nodecolors=nodecolors)
    return f, cost
end

function Base.show(io::IO, land::quantumLand)
    for i in land.villages
        println(i)
    end
end
function Base.show(io::IO, p::villager)
    print("Qubit $(p.id), who lives in hamlet $(p.address)")
end
function Base.show(io::IO, village::quantumVillage)
    print("Village $(village.id) has $(village.numVillagers) villagers, $(village.numFarms) farms, and $(village.numMayors) mayors.\n")
    for p in village.villagers
        println(p)
    end
end
end # module QuantumHamlet
