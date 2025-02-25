using QuantumHamlet
using QuantumHamlet: quantumLand
using QuantumClifford
using Graphs
using GraphMakie
using CairoMakie
using Colors

function random_graphstate(n::Int)
    rand_stab = QuantumClifford.random_stabilizer(n)
    return QuantumClifford.graphstate(rand_stab)
end

numVillages = 3
logicalBitsPerVillage = 4
reg = QuantumHamlet.defaultVillagerRegistry(numVillages,logicalBitsPerVillage)

land = quantumLand(reg, numVillages, numVillages*logicalBitsPerVillage)

rand_graphstate = random_graphstate(logicalBitsPerVillage*numVillages)
g = rand_graphstate[1]

cost, edgecolors = QuantumHamlet.naive_cost_of_partition(land,g)

function distinct_colors(n::Int)
    # Generate `n` distinct colors using HSV space
    return [HSV(i * 360 / n, 0.8, 0.9) for i in 0:n-1] .|> RGB
end

function print_graph(g::Graph; edgecolors=nothing, nodecolors=nothing)
    f, ax, p = graphplot(g,
                        node_color=nodecolors,
                        edge_color = edgecolors,
                        nlabels_align=(:center,:center);
                        ilabels= collect(1:nv(g)))
    hidedecorations!(ax); hidespines!(ax); ax.aspect = DataAspect()

    return f
end

village_colors = distinct_colors(numVillages)
nodecolors = [village_colors[reg[i]] for i in 1:numVillages*logicalBitsPerVillage]
f = print_graph(g, edgecolors=edgecolors, nodecolors=nodecolors)
