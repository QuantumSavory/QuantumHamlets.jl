using QuantumHamlet
using QuantumHamlet: quantumLand
using QuantumClifford
using Graphs
using GraphMakie
using CairoMakie
using Colors
using Statistics

function random_graphstate(n::Int)
    rand_stab = QuantumClifford.random_stabilizer(n)
    return QuantumClifford.graphstate(rand_stab)
end

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

#### Init
numVillages = 6
logicalBitsPerVillage = 4
reg = QuantumHamlet.defaultVillagerRegistry(numVillages,logicalBitsPerVillage)
rand_graphstate = random_graphstate(logicalBitsPerVillage*numVillages)
g = rand_graphstate[1]
village_colors = distinct_colors(numVillages)

#### Random initial labeling 
land = quantumLand(reg, numVillages, numVillages*logicalBitsPerVillage)

cost, edgecolors = QuantumHamlet.naive_cost_of_partition(land,g)

nodecolors = [village_colors[reg[i]] for i in 1:numVillages*logicalBitsPerVillage]
f = print_graph(g, edgecolors=edgecolors, nodecolors=nodecolors)

#### n/2 approximation algorithm for balanced k partitioning
new_reg, dump = QuantumHamlet.k_partition_saran_vazirani(g, numVillages)

# TODO instead of creating a new land, maybe add a function to relabel qubits? maybe?
improvedLand = quantumLand(new_reg, numVillages, numVillages*logicalBitsPerVillage)

cost2, edgecolors2 = QuantumHamlet.naive_cost_of_partition(improvedLand,g)

nodecolors2 = [village_colors[new_reg[i]] for i in 1:numVillages*logicalBitsPerVillage]
f2 = print_graph(g, edgecolors=edgecolors2, nodecolors=nodecolors2)

#### Random balanced partitionings 
samples = 50
worst_cost = 100000000
best_cost =  0
f_random_best = nothing
random_costs = []

for _ in 1:samples
    random_reg = QuantumHamlet.k_partition_random(g, numVillages)

    random_land = quantumLand(random_reg, numVillages, numVillages*logicalBitsPerVillage)

    random_cost, edgecolors = QuantumHamlet.naive_cost_of_partition(random_land,g)

    nodecolors = [village_colors[random_reg[i]] for i in 1:numVillages*logicalBitsPerVillage]
    f_random = print_graph(g, edgecolors=edgecolors, nodecolors=nodecolors)

    # if random_cost < best_cost
    #     best_cost = random_cost
    # end

    # if random_cost > worst_cost
    #     worst_cost = random_cost
    # end
    push!(random_costs, random_cost)
end

mean(random_costs)