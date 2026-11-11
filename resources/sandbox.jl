using QuantumHamlets
using QuantumHamlets: quantumLand
using Graphs
using GraphMakie
using CairoMakie
using Colors
using Statistics
using IGraphs
using Metis

using GraphIO
using GraphIO.GraphML
using EzXML
#### Init
# numHamlets = 2
# hamletSize = 98

# n = hamletSize * numHamlets
#g = random_regular_graph_deg3(n)
#g = graphstate(ghz(9))[1]
#g = grid([2, 2])
# g = grid_graph(n)
#g = graphstate(random_stabilizer(9))[1]

g = loadgraph("./QuantumHamlets.jl/resources/graphml_evolution/k1.graphml", GraphIO.GraphML.GraphMLFormat())
n = nv(g)

bury_reg, _  = QuantumHamlets.bury_not_balanced_bipartition(g, ceil(Int, n/3), floor(Int,n/2))
buryLand = QuantumHamlets.quantumLand(bury_reg, 2, n)
f_bury, bury_cost = QuantumHamlets.visualize_graph_on_land(buryLand, g, method=QuantumHamlets.matching_cost_of_partition)
bury_cost_cr, _ = QuantumHamlets.cut_rank_cost(buryLand, g)

bury_reg_seeded, _  = QuantumHamlets.bury_not_balanced_bipartition(g, ceil(Int, n/3), floor(Int,n/2), bury_v=QuantumHamlets.bury_seeded)
buryLand_seeded = QuantumHamlets.quantumLand(bury_reg_seeded, 2, n)
f_bury_seeded, bury_cost_seeded = QuantumHamlets.visualize_graph_on_land(buryLand_seeded, g, method=QuantumHamlets.matching_cost_of_partition)
bury_cost_cr_seeded, _ = QuantumHamlets.cut_rank_cost(buryLand_seeded, g)

#### Random balanced partitionings 
function random_sample(g::Graphs.Graph, numHamlets, numVillagers; method=QuantumHamlets.matching_cost_of_partition, samples=50, vis=true)
    best_cost =  [10000000]
    f_random_best = [Figure()]
    random_costs = []
    best_land = [quantumLand(Dict(), 0, 0)]

    for _ in 1:samples
        random_reg = QuantumHamlets.k_partition_random(g, numHamlets)
        random_land = quantumLand(random_reg, numHamlets, numHamlets*numVillagers)

        if vis
            f_random, random_cost = QuantumHamlets.visualize_graph_on_land(random_land, g, method=method)
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

#best_randomLand_naive, f_random_best_naive, random_costs_naive = random_sample(g, numHamlets, hamletSize, method=QuantumHamlets.naive_cost_of_partition)

# this one is the random one to use?
#best_randomLand_matching, f_random_best_matching, random_costs_matching = random_sample(g, numHamlets, hamletSize, method=QuantumHamlets.matching_cost_of_partition)

# bury_reg = QuantumHamlets.bury_k(g, numHamlets)
# buryLand = quantumLand(bury_reg, numHamlets, numHamlets*hamletSize)
# f_bury, bury_cost = QuantumHamlets.visualize_graph_on_land(buryLand, g, method=QuantumHamlets.matching_cost_of_partition)

# bury_seeded_reg = QuantumHamlets.bury_k_seeded(g, numHamlets)
# burySeededlLand = quantumLand(bury_seeded_reg, numHamlets, numHamlets*hamletSize)
# f_bury_seeded, bury_cost_seeded = QuantumHamlets.visualize_graph_on_land(burySeededlLand, g, method=QuantumHamlets.matching_cost_of_partition)

# metis_reg = Dict()
# metis_part = Metis.partition(g,numHamlets)
# for a in eachindex(metis_part)
#     metis_reg[a] = metis_part[a]
# end

# metisLand = quantumLand(metis_reg, numHamlets, numHamlets*hamletSize)
# f_metis, metis_cost = QuantumHamlets.visualize_graph_on_land(metisLand, g, method=QuantumHamlets.matching_cost_of_partition)

function random_regular_graph_deg3(numVertices)
    random_regular_graph(numVertices, 3)
end

function random_regular_graph_deg4(numVertices)
    random_regular_graph(numVertices, 4)
end

function random_regular_graph_deg5(numVertices)
    random_regular_graph(numVertices, 5)
end

function random_regular_graph_deg6(numVertices)
    random_regular_graph(numVertices, 6)
end

function my_erdos(n, p=4/n)
    erdos_renyi(n, p)
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

function compare_partition_methods(graph_generator, numHamlets; name=string(graph_generator), numSamples=50)
    pt = 4/3
    f = Figure(size=(900, 700),px_per_unit = 5.0, fontsize = 21pt)

    ax = f[1,1] = Axis(f[1,1],  xlabel="Logical Qubits per Hamlet",ylabel="Cut Rank",title=string(numHamlets)*" Hamlets partitioning "*name)
    
    # Should maybe make this less hardcoded
    colors = Makie.wong_colors()

    sizes = 20:4:100
    bury_costs = []
    bury_stds = []

    bury_costs_local = []
    bury_stds_local = []

    metis_costs = []
    metis_stds = []

    saran_costs = []
    saran_stds = []
    
    kl_costs = []
    kl_stds = []

    randommin_costs = []
    randommin_stds = []

    for hamletSize in sizes
        bury_sample_arr = []
        bury_sample_arr_local = []

        metis_sample_arr = []

        saran_sample_arr = []
        randommin_sample_arr = []
        kl_sample_arr = []
        for _ in 1:numSamples
            g = graph_generator(numHamlets*hamletSize)

            bury_reg = QuantumHamlets.bury_k(g, numHamlets)
            buryLand = quantumLand(bury_reg, numHamlets, numHamlets*hamletSize)
            #bury_cost, _ = QuantumHamlets.matching_cost_of_partition(buryLand, g)
            bury_cost = QuantumHamlets.cut_rank_cost(buryLand, g)

            push!(bury_sample_arr, bury_cost)

            bury_reg_local = QuantumHamlets.bury_k_seeded(g, numHamlets)
            buryLand_local = quantumLand(bury_reg_local, numHamlets, numHamlets*hamletSize)
            #bury_cost_local, _ = QuantumHamlets.matching_cost_of_partition(buryLand_local, g)
            bury_cost_local = QuantumHamlets.cut_rank_cost(buryLand_local, g)

            push!(bury_sample_arr_local, bury_cost_local)

            metis_reg = Dict()
            metis_part = Metis.partition(g,numHamlets)
            for a in eachindex(metis_part)
                metis_reg[a] = metis_part[a]
            end

            metisLand = quantumLand(metis_reg, numHamlets, numHamlets*hamletSize)
            #metis_cost, _ = QuantumHamlets.matching_cost_of_partition(metisLand, g)
            metis_cost = QuantumHamlets.cut_rank_cost(metisLand, g)

            push!(metis_sample_arr, metis_cost)

            saran_reg, _ = QuantumHamlets.k_partition_saran_vazirani(g, numHamlets)
            saranLand = quantumLand(saran_reg, numHamlets, numHamlets*hamletSize)
            #saran_cost, _ = QuantumHamlets.matching_cost_of_partition(saranLand, g)
            saran_cost = QuantumHamlets.cut_rank_cost(saranLand, g)

            push!(saran_sample_arr, saran_cost)

            # if numHamlets == 2
            #     py_g = nx.Graph()
            #     py_g.add_edges_from([(src(e), dst(e)) for e in edges(g)])
            #     kl_part = nxcomm.kernighan_lin.kernighan_lin_bisection(py_g, partition=nothing, max_iter=10)
            #     a, b = pyconvert(Vector{Int}, kl_part[0]), pyconvert(Vector{Int}, kl_part[1])
            #     kl_reg = Dict()
            #     for i in a
            #         kl_reg[i] = 1
            #     end
            #     for i in b
            #         kl_reg[i] = 2
            #     end
            #     klLand = quantumLand(kl_reg, numHamlets, numHamlets*hamletSize)
            #     #kl_cost, _ = QuantumHamlets.matching_cost_of_partition(klLand, g)
            #     kl_cost = cut_rank_cost(klLand, g)

            #     push!(kl_sample_arr, kl_cost)
            # end

            # _, _, random_costs = random_sample(g, numHamlets, hamletSize, method=QuantumHamlets.matching_cost_of_partition, samples=100, vis=false)
            # push!(randommin_sample_arr, minimum(random_costs))
        end
        push!(bury_costs, mean(bury_sample_arr))
        push!(bury_stds, std(bury_sample_arr))

        push!(bury_costs_local, mean(bury_sample_arr_local))
        push!(bury_stds_local, std(bury_sample_arr_local))


        push!(metis_costs, mean(metis_sample_arr))
        push!(metis_stds, std(metis_sample_arr))
        push!(saran_costs, mean(saran_sample_arr))
        push!(saran_stds, std(saran_sample_arr))

        # if numHamlets == 2
        #     push!(kl_costs, mean(kl_sample_arr))
        #     push!(kl_stds, std(kl_sample_arr))
        # end

        # push!(randommin_costs, mean(randommin_sample_arr))
        # push!(randommin_stds, std(randommin_sample_arr))
    end
    scatter!(ax, sizes, bury_costs, color=colors[1], marker=:circle)
    errorbars!(ax, sizes, bury_costs, bury_stds, color=colors[1], whiskerwidth=10)

    scatter!(ax, sizes, bury_costs_local, color=colors[2], marker=:circle)
    errorbars!(ax, sizes, bury_costs_local, bury_stds_local, color = colors[2], whiskerwidth=10)

    scatter!(ax, sizes, saran_costs, color=colors[3], marker=:circle)
    errorbars!(ax, sizes, saran_costs, saran_stds, color = colors[3], whiskerwidth=10)

    # if numHamlets == 2
    #     scatter!(ax, sizes, kl_costs, color=colors[4], marker=:circle)
    #     errorbars!(ax, sizes, kl_costs, kl_stds, color = colors[4], whiskerwidth=10)
    # end

    # scatter!(ax, sizes, randommin_costs, color=colors[5], marker=:circle)
    # errorbars!(ax, sizes, randommin_costs, randommin_stds, color = colors[5], whiskerwidth=10)

    scatter!(ax, sizes, metis_costs, color=colors[6], marker=:circle)
    errorbars!(ax, sizes, metis_costs, metis_stds, color = colors[6], whiskerwidth=10)

    #lines!(ax, [0,0], [0,0], label="n(k-1)/2 Bound", color=:black)
    
    lines!(ax, [0,0], [0,0], label="BURY", color=colors[1], linewidth=12)
    lines!(ax, [0,0], [0,0], label="BURY - \nSEEDED", color=colors[2], linewidth=12)
    #lines!(ax, [0,0], [0,0], label="KL", color=colors[4], linewidth=12)
    lines!(ax, [0,0], [0,0], label="Saran-Vazirani", color=colors[3], linewidth=12)
    #lines!(ax, [0,0], [0,0], label="Random\nSampling 100", color=colors[5], linewidth=12)
    
    lines!(ax, [0,0], [0,0], label="METIS", color=colors[6], linewidth=12)
    #lines!(ax, [0, maximum(sizes)], [0, (maximum(sizes)*numHamlets)*(numHamlets-1)/2], color=:gray)
 
    xlims!(ax, low=18)
    ylims!(ax, low=5)
    f[1,2] = Legend(f, ax, "Partition Method")
    return f
end

function compare_generation_methods(graph_generator, numHamlets)
    pt = 4/3
    f = Figure(size=(900, 700),px_per_unit = 5.0, fontsize = 15pt)

    ax = f[1,1] = Axis(f[1,1],  xlabel="Logical Qubits per Hamlet",ylabel="Required Bell pairs",title=string(numHamlets)*" Hamlets on "*string(graph_generator))
    
    sizes = 5:5:100

    bury_costs = []; bury_ghz_costs_optimistic = []; bury_ghz_costs_pessimistic = [];
    for hamletSize in sizes
        g = graph_generator(numHamlets*hamletSize)

        bury_reg = QuantumHamlets.bury_k(g, numHamlets)
        buryLand = quantumLand(bury_reg, numHamlets, numHamlets*hamletSize)
        bury_cost, _ = QuantumHamlets.matching_cost_of_partition(buryLand, g)
        cheap, expen = QuantumHamlets.vertex_cover_cost_of_partition(buryLand,g)
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

    #lines!(ax, [0, maximum(sizes)], [0, (maximum(sizes)*numHamlets)*(numHamlets-1)/2], color=:gray)
 
    f[1,2] = Legend(f, ax, "Graph Preparation method")
    return f
end

# ### KL portion
# using CondaPkg
# CondaPkg.add("networkx")
# using PythonCall
# nx = pyimport("networkx")
# nxcomm = pyimport("networkx.algorithms.community")

# py_g = nx.Graph()
# py_g.add_edges_from([(src(e), dst(e)) for e in edges(g)])
# kl_part = nxcomm.kernighan_lin.kernighan_lin_bisection(py_g, partition=nothing, max_iter=10)
# a, b = pyconvert(Vector{Int}, kl_part[0]), pyconvert(Vector{Int}, kl_part[1])

# kl_reg = Dict()
# for i in a
#     kl_reg[i] = 1
# end
# for i in b
#     kl_reg[i] = 2
# end
# klLand = quantumLand(kl_reg, numHamlets, numHamlets*hamletSize)
# f_kl, kl_cost = QuantumHamlets.visualize_graph_on_land(klLand, g, method=QuantumHamlets.matching_cost_of_partition)

#f_3reg = compare_partition_methods(random_regular_graph_deg3, 2, name="Random 3-Regular Graphs", numSamples=5)
#f_6reg = compare_partition_methods(random_regular_graph_deg6, 2, name="Random 6-Regular Graphs", numSamples=50)