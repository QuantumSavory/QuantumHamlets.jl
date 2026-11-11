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
