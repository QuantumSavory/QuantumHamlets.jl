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

using GraphIO
using GraphIO.GraphML
using EzXML

#### 
# g = loadgraph("./QuantumHamlet/QAOAgraphs/qaoa_100.graphml", GraphIO.GraphML.GraphMLFormat())
# n = nv(g)

# numVillages = 2
# @assert n%numVillages == 0
# logicalBitsPerVillage = n ÷ numVillages

# bury_reg = QuantumHamlet.bury_heuristic_global_v2(g, numVillages)
# buryLand = quantumLand(bury_reg, numVillages, numVillages*logicalBitsPerVillage)
# bury_cost, _ = QuantumHamlet.matching_cost_of_partition(buryLand, g)

# metis_reg = Dict()
# metis_part = Metis.partition(g,numVillages)
# for a in eachindex(metis_part)
#     metis_reg[a] = metis_part[a]
# end

# metisLand = quantumLand(metis_reg, numVillages, numVillages*logicalBitsPerVillage)
# metis_cost, _ = QuantumHamlet.matching_cost_of_partition(metisLand, g)

# Plots performance of METIS versus BURY on compiled graph states provided by a path to a directory with graphml files
function metis_bury_comparison_from_dir(graphml_dir; dir_name = "", hamletSizes = 2:8)
    pt = 4/3
    f = Figure(size=(900, 700),px_per_unit = 6.0, fontsize = 21pt)
    ax = f[1,1] = Axis(f[1,1],  xlabel="Pre-compilation QAOA Circuit Size ",ylabel="VCG Required Bell Pairs",title="METIS vs BURY on "*string(dir_name))
    
    bury_costs = Dict()
    metis_costs = Dict()
    for k in hamletSizes
        bury_costs[k] = []
        metis_costs[k] = []
    end

    for file in readdir(graphml_dir)
        println(file)
        g = loadgraph(graphml_dir*file, GraphIO.GraphML.GraphMLFormat())
        n = nv(g)
        circuit_size = parse(Int, split(file, ['_', '.'])[2])
        for k in hamletSizes
            if n%k != 0
                continue
            end

            bury_reg = QuantumHamlet.bury_heuristic_global_v2(g, k)
            buryLand = quantumLand(bury_reg, k, n)
            bury_cost, _ = QuantumHamlet.matching_cost_of_partition(buryLand, g)

            push!(bury_costs[k], [circuit_size,bury_cost])

            metis_reg = Dict()
            metis_part = Metis.partition(g,k)
            for a in eachindex(metis_part)
                metis_reg[a] = metis_part[a]
            end

            metisLand = quantumLand(metis_reg, k, n)
            metis_cost, _ = QuantumHamlet.matching_cost_of_partition(metisLand, g)
            push!(metis_costs[k], [circuit_size,metis_cost])
        end
    end

    wong_colors = Makie.wong_colors()
    colorCount = 0
    for k in hamletSizes
        if length(bury_costs[k]) == 0
            continue
        end
        sort!(bury_costs[k], by = x -> x[1])
        sort!(metis_costs[k], by = x -> x[1])
        x_bury = [data[1] for data in bury_costs[k]]
        y_bury = [data[2] for data in bury_costs[k]]

        x_metis = [data[1] for data in metis_costs[k]]
        y_metis = [data[2] for data in metis_costs[k]]

        scatterlines!(ax, x_bury, y_bury, color=wong_colors[(colorCount%7)+1], marker=:circle, markersize=15)
        scatterlines!(ax, x_metis, y_metis, color=wong_colors[(colorCount%7)+1], marker=:xcross, markersize=15)

        lines!(ax, [0,0], [0,0], label=string(k), color=wong_colors[(colorCount%7)+1], linewidth=12)
        colorCount+= 1
    end

    scatter!(ax, [0,0], [0,0], label="METIS", color=:gray, marker=:xcross, markersize=16)
    scatter!(ax, [0,0], [0,0], label="BURY", color=:gray, marker=:circle, markersize=16)
 
    xlims!(low = 25)
    f[1,2] = Legend(f, ax, "Number of\nHamlets")
    return f
end

f_qaoa = metis_bury_comparison_from_dir("./QuantumHamlet/QAOAgraphs/", dir_name="QAOA", hamletSizes= [4, 8, 10, 13, 18, 20])
