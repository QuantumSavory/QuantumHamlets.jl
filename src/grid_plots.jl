using QuantumHamlet
using QuantumHamlet: quantumLand
using QuantumClifford
using Graphs
using GraphMakie
using CairoMakie
using Colors
using Statistics
using IGraphs
using Metis
### KL portion
using CondaPkg
CondaPkg.add("networkx")
using PythonCall
nx = pyimport("networkx")
nxcomm = pyimport("networkx.algorithms.community")

using GraphIO.GraphML
using EzXML

function grid_graph(numVertices)
    a = sqrt(numVertices)
    b = floor(a) |> Int
    while numVertices%b != 0
        b -= 1
    end

    a = numVertices÷b
    return Graphs.grid([a,b])
end

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

# This function was generated with ChatGPT
function distinct_colors(n::Int)
    # Generate `n` distinct colors using HSV space
    return [HSV(i * 360 / n, 0.8, 0.9) for i in 0:n-1] .|> RGB
end

# Plots performance of METIS versus BURY on compiled graph states provided by a path to a directory with graphml files
function deterministic_comparison(graph_generator; name="", maxHamlets=10)
    pt = 4/3
    f = Figure(size=(900, 700),px_per_unit = 5.0, fontsize = 20pt)
    ax = f[1,1] = Axis(f[1,1],  xlabel="Vertices in Graph", ylabel="VCG Required Bell Pairs",title="METIS vs BURY on "*name)
    
    bury_costs = [[] for _ in 1:maxHamlets]
    metis_costs = [[] for _ in 1:maxHamlets]

    for i in 5:40
        g = graph_generator(i*i)
        n = nv(g)
        for k in 2:maxHamlets
            if n%k != 0
                continue
            end

            bury_reg = QuantumHamlet.bury_heuristic_global_v2(g, k)
            buryLand = quantumLand(bury_reg, k, n)
            bury_cost, _ = QuantumHamlet.matching_cost_of_partition(buryLand, g)

            push!(bury_costs[k], [n,bury_cost])

            metis_reg = Dict()
            metis_part = Metis.partition(g,k)
            for a in eachindex(metis_part)
                metis_reg[a] = metis_part[a]
            end

            metisLand = quantumLand(metis_reg, k, n)
            metis_cost, _ = QuantumHamlet.matching_cost_of_partition(metisLand, g)
            push!(metis_costs[k], [n,metis_cost])
        end
    end

    k_colors = distinct_colors(maxHamlets-1)
    for k in 2:maxHamlets
        sort!(bury_costs[k], by = x -> x[1])
        sort!(metis_costs[k], by = x -> x[1])
        x_bury = [data[1] for data in bury_costs[k]]
        y_bury = [data[2] for data in bury_costs[k]]

        x_metis = [data[1] for data in metis_costs[k]]
        y_metis = [data[2] for data in metis_costs[k]]

        scatterlines!(ax, x_bury, y_bury, color=k_colors[k-1], marker=:circle, markersize=12)
        scatterlines!(ax, x_metis, y_metis, color=k_colors[k-1], marker=:xcross, markersize=12)

        lines!(ax, [0,0], [0,0], label=string(k), color=k_colors[k-1], linewidth = 12)
    end

    scatter!(ax, [0,0], [0,0], label="METIS", color=:gray, marker=:xcross, markersize= 16)
    scatter!(ax, [0,0], [0,0], label="BURY", color=:gray, marker=:circle, markersize = 16)
 
    f[1,2] = Legend(f, ax, "Number of\nHamlets")
    return f, bury_costs, metis_costs
end

f_grid, bury_costs, metis_costs = deterministic_comparison(grid_graph, name="Grid Graphs")
#f_3reg, bury_costs, metis_costs = grid_plots(random_regular_graph_deg3, name="3-regular Graphs")