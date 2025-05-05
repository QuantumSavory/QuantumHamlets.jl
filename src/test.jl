using QuantumHamlet
using QuantumHamlet: quantumLand
using QuantumClifford
using Graphs
using GraphMakie
using CairoMakie
using Colors
using Statistics
using IGraphs

function random_graphstate(n::Int)
    rand_stab = QuantumClifford.random_stabilizer(n)
    return QuantumClifford.graphstate(rand_stab)
end

#### Init
numVillages = 2
logicalBitsPerVillage = 10
reg = QuantumHamlet.defaultVillagerRegistry(numVillages,logicalBitsPerVillage)
# rand_graphstate = random_graphstate(logicalBitsPerVillage*numVillages)
# g = rand_graphstate[1]
#g = grid([8,8])
g = random_regular_graph(20, 7)

#### n/2 approximation algorithm for balanced k partitioning
saran_reg, _ = QuantumHamlet.k_partition_saran_vazirani(g, numVillages)

# TODO instead of creating a new land, maybe add a function to relabel qubits? maybe?
saranLand = quantumLand(saran_reg, numVillages, numVillages*logicalBitsPerVillage)

saran_cost, _ = QuantumHamlet.naive_cost_of_partition(saranLand,g)
f_saran, saran_cost = QuantumHamlet.visualize_graph_on_land(saranLand, g, method=QuantumHamlet.matching_cost_of_partition)

#### Random balanced partitionings 
function random_sample(g::Graphs.Graph, numVillages, numVillagers; method=QuantumHamlet.naive_cost_of_partition, samples=100)
    best_cost =  [10000000]
    f_random_best = [Figure()]
    random_costs = []
    best_land = [quantumLand(Dict(), 0, 0)]
    problem_land = [quantumLand(Dict(), 0, 0)]

    for _ in 1:samples
        random_reg = QuantumHamlet.k_partition_random(g, numVillages)
        random_land = quantumLand(random_reg, numVillages, numVillages*numVillagers)
        #random_cost, _ = method(random_land,g)
        try
            f_random, random_cost = QuantumHamlet.visualize_graph_on_land(random_land, g, method=method)

            if random_cost < best_cost[1]
                println("FOUND BETTER")
                best_cost[1] = random_cost
                f_random_best[1] = f_random
                best_land[1] = random_land
            end
    
            push!(random_costs, random_cost)
        catch KeyError
            println("problem_land founrd")
            problem_land = random_land
        end
    end

    return best_land[1], f_random_best[1], random_costs, problem_land
end

best_randomLand_naive, f_random_best_naive, random_costs_naive, _ = random_sample(g, numVillages, logicalBitsPerVillage)
best_randomLand_matching, f_random_best_matching, random_costs_matching, problem_land = random_sample(g, numVillages, logicalBitsPerVillage, method=QuantumHamlet.matching_cost_of_partition)

# TODO Only defined currently when numVillages = 2
bury_reg = QuantumHamlet.bury_heuristic_global(g, numVillages)
buryLand = quantumLand(bury_reg, numVillages, numVillages*logicalBitsPerVillage)
f_bury, bury_cost = QuantumHamlet.visualize_graph_on_land(buryLand, g, method=QuantumHamlet.matching_cost_of_partition)

#QuantumHamlet.visualize_graph_on_land(best_randomLand_naive, g, method=QuantumHamlet.matching_cost_of_partition)