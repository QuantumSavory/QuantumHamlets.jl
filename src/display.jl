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
