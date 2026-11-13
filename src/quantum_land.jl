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

"""
    quantumLand(citizenRegistry, numHamlets, numVillagers)

Legacy representation of villagers grouped into hamlets, constructed from a
vertex-to-hamlet dictionary and explicit counts. Existing cost and plotting
functions accept this object. For the simpler graph-first API, use
[`GraphPartition`](@ref), with `quantumLand(p)` available for conversion.

The dictionary is copied. Counts must be integers, with positive `numHamlets`
and nonnegative `numVillagers`; keys must cover exactly `1:numVillagers` and
labels must be integers in `1:numHamlets`. Empty hamlets are supported.
Treat the registry and stored membership objects as read-only after construction.
"""
struct quantumLand
    hamlets::Vector{hamlet}
    registry::Dict{Int,Int}
    numHamlets::Int
    numVillagers::Int
end
function quantumLand(citizenRegistry::AbstractDict, numHamlets, numVillagers)
    numHamlets isa Integer && 1 <= numHamlets <= typemax(Int) ||
        throw(ArgumentError("numHamlets must be a positive integer representable as Int"))
    numVillagers isa Integer && 0 <= numVillagers <= typemax(Int) ||
        throw(ArgumentError("numVillagers must be a nonnegative integer representable as Int"))
    n, k = Int(numVillagers), Int(numHamlets)
    length(citizenRegistry) == n &&
        all(v -> v isa Integer && 1 <= v <= n, keys(citizenRegistry)) &&
        all(v -> haskey(citizenRegistry, v), 1:n) ||
        throw(ArgumentError("registry must contain exactly the vertex keys 1:numVillagers"))
    all(label -> label isa Integer && 1 <= label <= k, values(citizenRegistry)) ||
        throw(ArgumentError("registry labels must be integers in 1:numHamlets"))
    owned_registry = Dict{Int,Int}(citizenRegistry)
    land = quantumLand([hamlet(i) for i in 1:k], owned_registry, k, n)

    for i in 1:n
        hamlet_id = owned_registry[i]
        add_villager!(land.hamlets[hamlet_id],villager(hamlet_id,i))
    end

    return land
end

# Legacy graph helpers.
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

# Text display for the legacy partition objects.
function Base.show(io::IO, land::quantumLand)
    for i in land.hamlets
        println(io, i)
    end
end
function Base.show(io::IO, p::villager)
    print(io, "Qubit $(p.id), who lives in hamlet $(p.address)")
end
function Base.show(io::IO, hamlet::hamlet)
    print(io, "Village $(hamlet.id) has $(hamlet.numVillagers) villagers, $(hamlet.numFarms) farms, and $(hamlet.numMayors) mayors.\n")
    for p in hamlet.villagers
        println(io, p)
    end
end
