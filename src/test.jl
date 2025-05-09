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
numVillages = 3
logicalBitsPerVillage = 12
#rand_graphstate = random_graphstate(logicalBitsPerVillage*numVillages)
#g = rand_graphstate[1]
g = grid([6,6])
#g = random_regular_graph(numVillages*logicalBitsPerVillage, 4)

#### n/2 approximation algorithm for balanced k partitioning
saran_reg, _ = QuantumHamlet.k_partition_saran_vazirani(g, numVillages)

# TODO instead of creating a new land, maybe add a function to relabel qubits? maybe?
saranLand = quantumLand(saran_reg, numVillages, numVillages*logicalBitsPerVillage)

saran_cost, _ = QuantumHamlet.naive_cost_of_partition(saranLand,g)
f_saran, saran_cost = QuantumHamlet.visualize_graph_on_land(saranLand, g, method=QuantumHamlet.matching_cost_of_partition)

#### Random balanced partitionings 
function random_sample(g::Graphs.Graph, numVillages, numVillagers; method=QuantumHamlet.matching_cost_of_partition, samples=50, vis=true)
    best_cost =  [10000000]
    f_random_best = [Figure()]
    random_costs = []
    best_land = [quantumLand(Dict(), 0, 0)]

    for _ in 1:samples
        random_reg = QuantumHamlet.k_partition_random(g, numVillages)
        random_land = quantumLand(random_reg, numVillages, numVillages*numVillagers)

        if vis
            f_random, random_cost = QuantumHamlet.visualize_graph_on_land(random_land, g, method=method)
        else
            random_cost, _ = method(random_land, g)
        end

        if random_cost < best_cost[1]
            best_cost[1] = random_cost
            best_land[1] = random_land
            if vis
                println("FOUND BETTER")
                f_random_best[1] = f_random
            end
        end

        push!(random_costs, random_cost)
    end

    return best_land[1], f_random_best[1], random_costs
end

#best_randomLand_naive, f_random_best_naive, random_costs_naive = random_sample(g, numVillages, logicalBitsPerVillage, method=QuantumHamlet.naive_cost_of_partition)
best_randomLand_matching, f_random_best_matching, random_costs_matching = random_sample(g, numVillages, logicalBitsPerVillage, method=QuantumHamlet.matching_cost_of_partition)

# TODO Only defined currently when numVillages = 2
bury_reg = QuantumHamlet.bury_heuristic_global_v2(g, numVillages)
buryLand = quantumLand(bury_reg, numVillages, numVillages*logicalBitsPerVillage)
f_bury, bury_cost = QuantumHamlet.visualize_graph_on_land(buryLand, g, method=QuantumHamlet.matching_cost_of_partition)

#QuantumHamlet.visualize_graph_on_land(best_randomLand_naive, g, method=QuantumHamlet.matching_cost_of_partition)
function random_regular_graph_deg3(numVertices)
    random_regular_graph(numVertices, 3)
end

function random_regular_graph_deg4(numVertices)
    random_regular_graph(numVertices, 4)
end

function random_regular_graph_deg5(numVertices)
    random_regular_graph(numVertices, 5)
end

function grid_graph(numVertices)
    a = sqrt(numVertices)
    b = floor(a) |> Int
    while numVertices%b != 0
        b -= 1
    end

    a = numVertices÷b
    return Graphs.grid([a,b])
end

function compare_partition_methods(graph_generator, numVillages; numSamples=20)
    pt = 4/3
    f = Figure(size=(900, 700),px_per_unit = 5.0, fontsize = 15pt)

    ax = f[1,1] = Axis(f[1,1],  xlabel="Logical Qubits per Village",ylabel="Required Bell pairs",title=string(numVillages)*" Villages on "*string(graph_generator))
    
    sizes = 4:2:30
    bury_costs = []
    bury_stds = []
    bury_ghz_costs_optimistic = []
    bury_ghz_costs_pessimistic = []

    saran_costs = []
    saran_stds = []

    randommin_costs = []
    randommin_stds = []

    for logicalBitsPerVillage in sizes
        bury_sample_arr = []; bsa_cheap = []; bsa_expensive = [];
        saran_sample_arr = []
        randommin_sample_arr = []
        for _ in 1:numSamples
            g = graph_generator(numVillages*logicalBitsPerVillage)

            bury_reg = QuantumHamlet.bury_heuristic_global_v2(g, numVillages)
            buryLand = quantumLand(bury_reg, numVillages, numVillages*logicalBitsPerVillage)
            bury_cost, _ = QuantumHamlet.matching_cost_of_partition(buryLand, g)
            cheap, expen = QuantumHamlet.vertex_cover_cost_of_partition(buryLand,g)
            push!(bury_sample_arr, bury_cost)
            push!(bsa_cheap, cheap)
            push!(bsa_expensive,expen)

            saran_reg, _ = QuantumHamlet.k_partition_saran_vazirani(g, numVillages)
            saranLand = quantumLand(saran_reg, numVillages, numVillages*logicalBitsPerVillage)
            saran_cost, _ = QuantumHamlet.matching_cost_of_partition(saranLand, g)
            push!(saran_sample_arr, saran_cost)

            _, _, random_costs = random_sample(g, numVillages, logicalBitsPerVillage, method=QuantumHamlet.matching_cost_of_partition, vis=false)
            push!(randommin_sample_arr, minimum(random_costs))
        end
        push!(bury_costs, mean(bury_sample_arr))
        push!(bury_stds, std(bury_sample_arr))
        push!(bury_ghz_costs_optimistic, mean(bsa_cheap))
        push!(bury_ghz_costs_pessimistic, mean(bsa_expensive))

        push!(saran_costs, mean(saran_sample_arr))
        push!(saran_stds, std(saran_sample_arr))

        push!(randommin_costs, mean(randommin_sample_arr))
        push!(randommin_stds, std(randommin_sample_arr))
    end
    scatter!(ax, sizes, bury_costs, color=:blue, marker=:circle)
    errorbars!(ax, sizes, bury_costs, bury_stds, color = :blue, whiskerwidth=10)

    scatter!(ax, sizes, bury_ghz_costs_optimistic, color=:blue, marker=:utriangle)
    scatter!(ax, sizes, bury_ghz_costs_pessimistic, color=:blue, marker=:cross)

    scatter!(ax, sizes, saran_costs, color=:green, marker=:circle)
    errorbars!(ax, sizes, saran_costs, saran_stds, color = :green, whiskerwidth=10)

    scatter!(ax, sizes, randommin_costs, color=:orange, marker=:circle)
    errorbars!(ax, sizes, randommin_costs, randommin_stds, color = :orange, whiskerwidth=10)

    lines!(ax, [0,0], [0,0], label="n(k-1)/2 Bound", color=:black)
    lines!(ax, [0,0], [0,0], label="Best partition found\nover 50 random samples", color=:orange)
    lines!(ax, [0,0], [0,0], label="Saran Vazirani n/2\napprox for edge min", color=:green)
    lines!(ax, [0,0], [0,0], label="Bury heuristic - Global", color=:blue)

    scatter!(ax, [0,0], [0,0], label="m party GHZ states\nhave cost m-1", color=:gray, marker=:cross)
    scatter!(ax, [0,0], [0,0], label="Cheap GHZ states", color=:gray, marker=:utriangle)

    lines!(ax, [0, maximum(sizes)], [0, (maximum(sizes)*numVillages)*(numVillages-1)/2], color=:gray)
 
    f[1,2] = Legend(f, ax, "Partition Method")
    return f
end

function compare_generation_methods(graph_generator, numVillages)
    pt = 4/3
    f = Figure(size=(900, 700),px_per_unit = 5.0, fontsize = 15pt)

    ax = f[1,1] = Axis(f[1,1],  xlabel="Logical Qubits per Village",ylabel="Required Bell pairs",title=string(numVillages)*" Villages on "*string(graph_generator))
    
    sizes = 5:5:100

    bury_costs = []; bury_ghz_costs_optimistic = []; bury_ghz_costs_pessimistic = [];
    for logicalBitsPerVillage in sizes
        g = graph_generator(numVillages*logicalBitsPerVillage)

        bury_reg = QuantumHamlet.bury_heuristic_global_v2(g, numVillages)
        buryLand = quantumLand(bury_reg, numVillages, numVillages*logicalBitsPerVillage)
        bury_cost, _ = QuantumHamlet.matching_cost_of_partition(buryLand, g)
        cheap, expen = QuantumHamlet.vertex_cover_cost_of_partition(buryLand,g)
        push!(bury_costs, bury_cost)
        push!(bury_ghz_costs_optimistic, cheap)
        push!(bury_ghz_costs_pessimistic,expen)
    end
    scatter!(ax, sizes, bury_costs, color=:blue, marker=:circle)
    scatter!(ax, sizes, bury_ghz_costs_optimistic, color=:blue, marker=:utriangle)
    scatter!(ax, sizes, bury_ghz_costs_pessimistic, color=:blue, marker=:cross)

    #lines!(ax, [0,0], [0,0], label="n(k-1)/2 Bound", color=:black)
    lines!(ax, [0,0], [0,0], label="Bury heuristic - Global", color=:blue)

    scatter!(ax, [0,0], [0,0], label="m party GHZ states\nhave cost m-1", color=:gray, marker=:cross)
    scatter!(ax, [0,0], [0,0], label="Bell pairs only", color=:gray, marker=:circle)
    scatter!(ax, [0,0], [0,0], label="Cheap GHZ states", color=:gray, marker=:utriangle)

    #lines!(ax, [0, maximum(sizes)], [0, (maximum(sizes)*numVillages)*(numVillages-1)/2], color=:gray)
 
    f[1,2] = Legend(f, ax, "Graph Preparation method")
    return f
end