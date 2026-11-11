using CairoMakie, Colors, Graphs, QuantumHamlets
include("counterexamples/path_hubs.jl")

g, hubs, connectors, cliques, optimum = path_hub_obstruction(3)
reg, _ = QuantumHamlets.bury(g, nv(g) ÷ 2)
candidate = GraphPartition(reg; nparts=2)
positions = Dict{Int,Point2f}()
for j in eachindex(hubs)
    x = 2.0f0 * (j - 1)
    positions[hubs[j]] = Point2f(x, 0)
    for (v, offset) in zip(cliques[j], ((0,-0.8),(-0.65,-2.05),(0.65,-2.05),(0,-1.6)))
        positions[v] = Point2f(x + offset[1], offset[2])
    end
end

fig = Figure(size=(1100, 650), fontsize=17, backgroundcolor=:white)
palette = [colorant"#167080", colorant"#EAA548"]
for (row, p, title) in ((1, candidate, "BURY: matching cost 5"),
                        (2, optimum, "Optimal: matching cost 1"))
    labels = membership(p)
    ax = Axis(fig[row, 1]; title, aspect=DataAspect(), titlesize=22)
    hidedecorations!(ax)
    hidespines!(ax)
    for e in edges(g)
        crossing = labels[src(e)] != labels[dst(e)]
        lines!(ax, [positions[src(e)], positions[dst(e)]];
               color=crossing ? colorant"#B7353F" : colorant"#BBC1C9",
               linewidth=crossing ? 2.5 : 1.6)
    end
    scatter!(ax, [positions[v] for v in vertices(g)];
             color=palette[labels], markersize=17,
             strokecolor=:white, strokewidth=1.2)
    limits!(ax, -1, 11, -2.35, 0.35)
end
Label(fig[3, 1], "Same 30-vertex graph; both partitions have 15 vertices per side.\nRed edges cross the partition. Each four-vertex group is a complete graph.";
      fontsize=16, padding=(0, 0, 4, 4))
out = joinpath(@__DIR__, "figures")
mkpath(out)
save(joinpath(out, "bury_obstruction.png"), fig; px_per_unit=2)
save(joinpath(out, "bury_obstruction.pdf"), fig)
println("Saved counterexample figure; costs ", matching_cost(g, candidate), " and ", matching_cost(g, optimum))
