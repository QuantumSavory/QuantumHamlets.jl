using QuantumHamlets, Graphs, CairoMakie, Statistics
using GraphIO.GraphML, EzXML

"""
    plot_toffoli_baselines(; collection="graphml_reduced", repeats=10,
        parameters=(3, 4, 5, 10), families=("tof", "barenco_tof"),
        temperatures=range(1.0, 0.1; length=10), output=nothing)

Compare BURY, sqrt-seeded BURY, METIS, and their simulated-annealing variants
on the saved Toffoli graphs. Return `(figure, data)`, where `data` contains one
row per actual partition, including its seed, rank, and part sizes.

The grouped bars show median GF(2) cut rank; whiskers show the observed minimum
and maximum over seeds `0:repeats-1`, not confidence intervals. Ordinary BURY
is evaluated once. Every partition is exactly balanced (off by one when odd).
METIS uses the explicit capacity repair in `metis_partition(...; balance=true)`;
SA retains its best encountered partition, including its initializer.

Run from the package root with `julia --project=..`, then:
```julia
include("resources/toffoli_baselines.jl")
fig, data = plot_toffoli_baselines()
display(fig)
save("toffoli_baselines.svg", fig)
```
Change `collection` to `"graphml_initial"` for the original circuits, or pass
`families=("tof",)` for Toffoli alone. Including this file runs no experiment.
"""
function plot_toffoli_baselines(; collection="graphml_reduced", repeats=10,
        parameters=(3, 4, 5, 10), families=("tof", "barenco_tof"),
        temperatures=range(1.0, 0.1; length=10), output=nothing)
    repeats isa Integer && repeats > 0 || throw(ArgumentError("repeats must be a positive integer"))
    isempty(parameters) && throw(ArgumentError("parameters must not be empty"))
    !isempty(families) && all(f -> f in ("tof", "barenco_tof"), families) ||
        throw(ArgumentError("families must contain tof and/or barenco_tof"))
    graphdir = joinpath(@__DIR__, collection)
    collectionlabel = replace(basename(normpath(graphdir)), "graphml_" => "")
    for family in families, parameter in parameters
        path = joinpath(graphdir, "$(family)_$(parameter).graphml")
        isfile(path) || throw(ArgumentError("graph not found: $path"))
    end
    temperatures = collect(temperatures)
    names = ["BURY", "sqrt BURY", "Random SA", "BURY + SA", "sqrt BURY + SA",
             "METIS (balanced)", "METIS + SA"]
    colors = Makie.wong_colors()
    data = NamedTuple[]
    fig = Figure(size=(1400, 360 * length(families) + 160), fontsize=17)
    legendplots = Any[]

    for (panel, family) in enumerate(families)
        vertexcounts = Int[]
        for parameter in parameters
            g = loadgraph(joinpath(graphdir, "$(family)_$(parameter).graphml"), GraphMLFormat())
            push!(vertexcounts, nv(g))
            p = bipartition_graph(g; objective=pairwise_cut_rank)
            push!(data, (family=family, parameter=parameter, n=nv(g), method=names[1],
                         seed=missing, rank=pairwise_cut_rank(g, p), sizes=partition_sizes(p)))
            for seed in 0:repeats-1
                # Zero SA sweeps returns the sqrt-selected initializer itself.
                sqrt_partition = bury_seeded_sa(g; seed, temperatures=Float64[])
                partitions = (sqrt_partition,
                    min_cutrank(g; seed, temperatures),
                    bury_sa(g; seed, temperatures),
                    bury_seeded_sa(g; seed, temperatures),
                    metis_partition(g; seed, balance=true),
                    metis_sa(g; seed, temperatures, balance=true))
                for (name, partition) in zip(names[2:end], partitions)
                    sizes = partition_sizes(partition)
                    maximum(sizes) - minimum(sizes) <= 1 || error("$name returned an unbalanced partition")
                    push!(data, (family=family, parameter=parameter, n=nv(g), method=name,
                                 seed=seed, rank=pairwise_cut_rank(g, partition), sizes=sizes))
                end
            end
        end

        title = family == "tof" ? "Toffoli" : "Barenco Toffoli"
        ticks = ["$parameter\n(n=$n)" for (parameter, n) in zip(parameters, vertexcounts)]
        ax = Axis(fig[panel+1, 1], title="$title — $collectionlabel graphs", ylabel="GF(2) cut rank ↓",
                  xlabel="Circuit parameter", xticks=(1:length(parameters), ticks))
        for (method, name) in enumerate(names)
            samples = [[row.rank for row in data if row.family == family &&
                        row.parameter == parameter && row.method == name] for parameter in parameters]
            centers = median.(samples)
            x = collect(1:length(parameters)) .+ (method - 4) * 0.11
            bars = barplot!(ax, x, centers; width=0.105, color=colors[method])
            errorbars!(ax, x, centers, centers .- minimum.(samples), maximum.(samples) .- centers;
                       color=:black, whiskerwidth=5, linewidth=1)
            panel == 1 && push!(legendplots, bars)
        end
        ylims!(ax, 0, nothing)
    end
    Legend(fig[1, 1], legendplots, names; orientation=:horizontal, nbanks=2, framevisible=false)
    Label(fig[length(families)+2, 1],
        "Medians and observed min–max; $repeats stochastic seeds. SA: $(length(temperatures)) sweeps, best partition retained.\n" *
        "Exact floor/ceiling capacities; METIS is repaired by minimum edge-cut-increase moves when needed.";
        fontsize=13, tellwidth=false)
    isnothing(output) || save(output, fig)
    return fig, data
end
