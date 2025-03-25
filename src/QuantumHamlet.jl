module QuantumHamlet
using QuantumClifford
using CairoMakie
using Graphs
using GraphsFlows
using GraphMakie
using Colors
using MultiwayNumberPartitioning, HiGHS

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
end
function quantumLand(citizenRegistry::Dict, numVillages::Int, numVillagers::Int)
    land = quantumLand([quantumVillage(i) for i in 1:numVillages],citizenRegistry)

    for i in 1:numVillagers
        village_id = citizenRegistry[i]
        add_villager!(land.villages[village_id],villager(village_id,i))
    end

    return land
end

function naive_cost_of_partition(land::quantumLand, graph::Graphs.Graph) 
    println(land)
    println(graph)
    
    edgecolors = [:black for _ in 1:ne(graph)]

    cost = 0
    for (index, e) in enumerate(edges(graph))
        if land.registry[src(e)] != land.registry[dst(e)]
            cost += 1
            edgecolors[index] = :red
        end
    end

    println("Naive cost is ", cost)
    return cost, edgecolors
end

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
