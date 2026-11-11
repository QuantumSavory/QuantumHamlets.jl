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

This constructor retains its original ownership and validation behavior;
see the usage guide's known limitations.
"""
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
