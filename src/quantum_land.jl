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
