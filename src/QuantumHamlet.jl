module QuantumHamlet
using QuantumClifford
using CairoMakie
using Graphs
using GraphsFlows
using GraphMakie
using Colors

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

function k_partition_saran_vazirani(g::Graphs.Graph)
    d = Graphs.DiGraph(g)
    capacity_matrix = zeros(Int, 9, 9)
    for e in edges(d)
        u, v = src(e), dst(e)
        capacity_matrix[u, v] = 1
    end

    min_cuts = []
    for e in edges(g)
        push!(min_cuts, GraphsFlows.mincut(d,src(e),dst(e),capacity_matrix,GraphsFlows.BoykovKolmogorovAlgorithm()))
    end

    sort!(min_cuts, by = x-> x[3])

    return min_cuts
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
