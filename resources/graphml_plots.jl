"""Reproducible METIS/BURY cost comparisons for resources/graphml*.

Including this file defines helpers without starting a sweep. Run through
run_graphml_plots.py for fresh output directories and external time limits.
"""
module GraphMLComparisonPlots

using QuantumHamlets, Graphs, CairoMakie, GraphIO.GraphML, EzXML, Metis, SHA
import GraphIO, Makie

const OBJECTIVE_LABELS = Dict(:matching => "Matching cost", :rank => "Pairwise cut rank (GF(2))")
const COLLECTION_LABELS = Dict("graphml_initial" => "Initial graphs",
                               "graphml_reduced" => "Reduced graphs",
                               "graphml_evolution" => "Evolution graphs")
const BURY_COLOR = "#0072B2"
const METIS_COLOR = "#D55E00"
const IMBALANCE_COLOR = "#C23B50"
collection_label(name) = get(COLLECTION_LABELS, name, replace(name, "graphml_" => ""))
collection_order(name) = (get(Dict("graphml_initial" => 1, "graphml_reduced" => 2,
                                 "graphml_evolution" => 3), name, 4), name)

objective_function(objective::Symbol) = objective === :matching ? matching_cost :
    objective === :rank ? pairwise_cut_rank : throw(ArgumentError("objective must be matching or rank"))

function graphml_inputs(resources_dir=@__DIR__)
    files = String[]
    for directory in sort(readdir(resources_dir; join=true))
        isdir(directory) && startswith(basename(directory), "graphml") || continue
        append!(files, sort(filter(p -> isfile(p) && endswith(lowercase(p), ".graphml"),
                                  readdir(directory; join=true))))
    end
    isempty(files) && throw(ArgumentError("no graphml*/*.graphml inputs in $resources_dir"))
    return files
end

supported_parts(n, k) = n >= 2 && 2 <= k <= n && (k == 2 || n % k == 0)

# Metis.jl 1.5 exposes native options as a global vector. Pin the native seed
# and standard KWAY imbalance tolerance, then restore every caller option.
function metis_partition(g, k; seed=0)
    saved = copy(Metis.options)
    try
        fill!(Metis.options, -1)
        Metis.options[Int(Metis.LibMetis.METIS_OPTION_NUMBERING) + 1] = 1
        Metis.options[Int(Metis.LibMetis.METIS_OPTION_SEED) + 1] = seed
        Metis.options[Int(Metis.LibMetis.METIS_OPTION_UFACTOR) + 1] = 30
        return GraphPartition(Metis.partition(g, k; alg=:KWAY); nparts=k)
    finally
        copyto!(Metis.options, saved)
    end
end

function balance_matches(p, n, k)
    sizes = partition_sizes(p)
    return all(>(0), sizes) && (k == 2 ? maximum(sizes) - minimum(sizes) <= 1 :
                               all(==(n ÷ k), sizes))
end

function write_tsv_row(io, values)
    println(io, join(values, '\t'))
    flush(io)
end

function build_comparisons(output; resources_dir=@__DIR__, parts=2:8,
                           objectives=[:matching, :rank], backend=:heap)
    parts = collect(parts)
    !isempty(parts) && all(k -> k isa Integer && k >= 2, parts) ||
        throw(ArgumentError("parts must contain integers at least 2"))
    parts = sort(unique(parts))
    objectives = unique(collect(objectives))
    !isempty(objectives) || throw(ArgumentError("choose at least one objective"))
    foreach(objective_function, objectives)
    backend in (:heap, :scan, :legacy) || throw(ArgumentError("invalid BURY backend"))
    files = graphml_inputs(resources_dir)
    mkpath(output)
    for name in ("results.tsv", "status.tsv", "inputs.tsv")
        ispath(joinpath(output, name)) && error("refusing to overwrite $name")
    end
    rows = NamedTuple[]
    open(joinpath(output, "results.tsv"), "w") do results
        open(joinpath(output, "status.tsv"), "w") do status
            open(joinpath(output, "inputs.tsv"), "w") do inputs
                write_tsv_row(results, ("collection", "graph", "vertices", "edges", "parts", "objective",
                    "algorithm", "cost", "part_sizes", "matches_bury_balance", "max_size_over_mean", "nonempty"))
                write_tsv_row(status, ("collection", "graph", "vertices", "parts", "status", "reason"))
                write_tsv_row(inputs, ("collection", "graph", "vertices", "edges", "sha256", "path"))
                for path in files
                    collection, name = basename(dirname(path)), splitext(basename(path))[1]
                    println("STAGE\tload\t", collection, '\t', name); flush(stdout)
                    g = loadgraph(path, GraphMLFormat())
                    g isa Graphs.Graph || error("expected simple undirected graph: $path")
                    original = copy(g)
                    n = nv(g)
                    write_tsv_row(inputs, (collection, name, n, ne(g), bytes2hex(sha256(read(path))),
                                           relpath(path, resources_dir)))
                    for k in parts
                        if !supported_parts(n, k)
                            write_tsv_row(status, (collection, name, n, k, "unsupported_balance",
                                k > n ? "more parts than vertices" : "exact multipartition requires n divisible by k"))
                            continue
                        end
                        println("STAGE\tcompare\t", collection, '\t', name, '\t', k); flush(stdout)
                        metis = metis_partition(g, k)
                        length(metis) == n || error("METIS returned incomplete membership")
                        metis_nonempty = all(>(0), partition_sizes(metis))
                        ordinary = k == 2 ? nothing : partition_graph(g, k; backend)
                        for objective in objectives
                            score = objective_function(objective)
                            bury = k == 2 ? bipartition_graph(g; objective=score, backend) : ordinary
                            balance_matches(bury, n, k) || error("BURY violated its balance contract")
                            for (algorithm, p) in (("BURY", bury), ("METIS", metis))
                                sizes = partition_sizes(p)
                                row = (collection=collection, graph=name, n=n, m=ne(g), k=k,
                                       objective=objective, algorithm=algorithm, cost=score(g, p),
                                       sizes=sizes, balanced=balance_matches(p, n, k),
                                       imbalance=maximum(sizes) / (n / k), nonempty=all(>(0), sizes))
                                push!(rows, row)
                                write_tsv_row(results, (collection, name, n, ne(g), k, objective, algorithm,
                                    row.cost, join(sizes, ','), row.balanced, row.imbalance, row.nonempty))
                            end
                        end
                        g == original || error("comparison modified graph $path")
                        write_tsv_row(status, (collection, name, n, k,
                            metis_nonempty ? "completed" : "metis_empty_part",
                            metis_nonempty ? "" : "METIS sizes=" * join(partition_sizes(metis), ',') * "; excluded from plots"))
                    end
                end
            end
        end
    end
    return rows
end

function save_plot(output, stem, figure)
    for extension in ("png", "svg")
        path = joinpath(output, stem * "." * extension)
        ispath(path) && error("refusing to overwrite $path")
        extension == "png" ? save(path, figure; px_per_unit=1.6) : save(path, figure)
    end
end

function plot_overview(rows, objective, output)
    selected = filter(r -> r.objective == objective && r.nonempty, rows)
    collections = sort(unique(r.collection for r in selected); by=collection_order)
    parts = sort(unique(r.k for r in selected))
    colors = Makie.wong_colors()
    fig = Figure(size=(1500, 620), fontsize=17, backgroundcolor=:white)
    Label(fig[0, 1:3], "METIS vs BURY · " * OBJECTIVE_LABELS[objective], fontsize=25, font=:bold)
    for (column, collection) in enumerate(collections)
        rs = filter(r -> r.collection == collection, selected)
        ax = Axis(fig[1, column], title=collection_label(collection),
                  xlabel="Vertices in graph", ylabel=column == 1 ? OBJECTIVE_LABELS[objective] : "",
                  xgridvisible=false, ygridcolor=(:gray, 0.15), titlesize=20)
        for (index, k) in enumerate(parts)
            color = colors[mod1(index, length(colors))]
            for (algorithm, marker) in (("BURY", :circle), ("METIS", :xcross))
                values = sort(filter(r -> r.k == k && r.algorithm == algorithm, rs); by=r -> (r.n, r.graph))
                isempty(values) && continue
                scatter!(ax, [r.n for r in values], [r.cost for r in values]; color, marker, markersize=10)
                imbalanced = filter(r -> !r.balanced, values)
                isempty(imbalanced) || scatter!(ax, [r.n for r in imbalanced], [r.cost for r in imbalanced];
                    marker=:circle, markersize=16, color=:transparent,
                    strokecolor=IMBALANCE_COLOR, strokewidth=1.3)
            end
        end
        ylims!(ax, low=0)
        hidespines!(ax, :t, :r)
    end
    part_elements = [MarkerElement(color=colors[mod1(i, length(colors))], marker=:circle, markersize=12)
                     for i in eachindex(parts)]
    method_elements = [MarkerElement(color=:gray35, marker=:circle, markersize=12),
                       MarkerElement(color=:gray35, marker=:xcross, markersize=12)]
    Legend(fig[2, 1:3], [part_elements, method_elements], [string.(parts), ["BURY", "METIS"]],
           ["Number of parts", "Algorithm"]; orientation=:horizontal, framevisible=false, tellheight=true)
    coverage = 2 in parts ? "k = 2: includes odd n  •  k ≥ 3: exactly divisible n only" : "Exactly divisible n only"
    Label(fig[3, 1:3], coverage * "  •  Lower cost is better  •  Empty-part METIS results excluded", fontsize=14)
    Label(fig[4, 1:3], "Red rings: METIS outside BURY's balance constraint. METIS optimizes edge cut with 3% balance tolerance.", fontsize=14, color=:gray35)
    save_plot(output, "overview_" * string(objective), fig)
end

function plot_bipartitions(rows, collection, objectives, output)
    selected = filter(r -> r.collection == collection && r.k == 2, rows)
    isempty(selected) && return
    names = sort(unique(r.graph for r in selected); by=name ->
        (first(r.n for r in selected if r.graph == name), name))
    labels = [replace(name, '_' => ' ') * "  (n=" * string(first(r.n for r in selected if r.graph == name)) * ")"
              for name in names]
    fig = Figure(size=(max(850, 650length(objectives)), 230 + 30length(names)), fontsize=16, backgroundcolor=:white)
    Label(fig[0, 1:length(objectives)], collection_label(collection) * " · near-balanced bipartitions",
          fontsize=24, font=:bold)
    for (column, objective) in enumerate(objectives)
        ax = Axis(fig[1, column], xlabel=OBJECTIVE_LABELS[objective],
                  yticks=(collect(eachindex(names)), column == 1 ? labels : fill("", length(names))),
                  yreversed=true, ygridvisible=false, xgridcolor=(:gray, 0.15),
                  title=OBJECTIVE_LABELS[objective], titlesize=19)
        for (index, name) in enumerate(names)
            pair = filter(r -> r.graph == name && r.objective == objective, selected)
            bury = only(filter(r -> r.algorithm == "BURY", pair))
            metis = only(filter(r -> r.algorithm == "METIS", pair))
            scatter!(ax, [bury.cost], [index]; color=BURY_COLOR, marker=:circle, markersize=11)
            metis.nonempty || continue
            lines!(ax, [bury.cost, metis.cost], [index, index]; color=:gray75, linewidth=1.7)
            scatter!(ax, [metis.cost], [index]; color=METIS_COLOR, marker=:xcross, markersize=12)
            if !metis.balanced
                scatter!(ax, [metis.cost], [index]; color=:transparent, marker=:circle,
                         strokecolor=IMBALANCE_COLOR, strokewidth=1.4, markersize=18)
            end
        end
        xlims!(ax, low=0)
        ylims!(ax, 0.4, length(names) + 0.6)
        hidespines!(ax, :t, :r)
    end
    elements = [MarkerElement(color=BURY_COLOR, marker=:circle, markersize=12),
                MarkerElement(color=METIS_COLOR, marker=:xcross, markersize=12),
                MarkerElement(color=:transparent, strokecolor=IMBALANCE_COLOR, strokewidth=1.4,
                              marker=:circle, markersize=16)]
    Legend(fig[2, 1:length(objectives)], elements,
           ["BURY", "METIS", "METIS outside BURY balance"]; orientation=:horizontal, framevisible=false)
    Label(fig[3, 1:length(objectives)], "BURY selects floor(n/2) or ceil(n/2) separately for each objective; METIS reuses one partition. Lower is better.", fontsize=13)
    Label(fig[4, 1:length(objectives)], "METIS uses edge-cut optimization and 3% balance tolerance; actual part sizes are in results.tsv.", fontsize=13, color=:gray35)
    save_plot(output, collection * "_bipartitions", fig)
end

function write_metadata(output, parts, objectives, backend)
    path = joinpath(output, "metadata.txt")
    ispath(path) && error("refusing to overwrite $path")
    open(path, "w") do io
        println(io, "Julia: ", VERSION, "\nThreads: ", Threads.nthreads())
        for package in (QuantumHamlets, Graphs, CairoMakie, GraphIO, EzXML, Metis)
            println(io, nameof(package), ": ", Base.pkgversion(package))
        end
        println(io, "BURY backend: ", backend, "\nParts: ", join(parts, ','), "\nObjectives: ", join(objectives, ','))
        println(io, "METIS: KWAY, edge-cut objective, equal target weights, native seed=0, UFACTOR=30 (3%)")
        println(io, "At k=2 BURY selects floor/ceil targets separately for each objective; METIS membership is reused across objectives.")
        println(io, "METIS empty-part outputs retain raw costs/sizes in tables but are excluded from plots (status=metis_empty_part).")
        println(io, "GraphML: simple undirected unweighted topology; XML node declaration order maps to 1:n.")
        println(io, "Node phase/type/layout and edge-type attributes do not enter these graph objectives.")
        for source in ("graphml_plots.jl", "run_graphml_plots.py", "../src/bury.jl", "../src/api.jl", "../src/costs.jl",
                       "../src/bury/heap.jl", "../src/bury/scan.jl", "../src/bury/legacy.jl")
            println(io, "sha256\t", source, '\t', bytes2hex(sha256(read(joinpath(@__DIR__, source)))))
        end
    end
end

function main(args=ARGS)
    options = Dict("--objectives" => "matching,rank", "--backend" => "heap", "--parts" => "2,3,4,5,6,7,8")
    iseven(length(args)) || error("arguments must be --name value pairs")
    for i in 1:2:length(args)
        args[i] in ("--output-dir", "--objectives", "--backend", "--parts") || error("unknown option $(args[i])")
        options[args[i]] = args[i + 1]
    end
    haskey(options, "--output-dir") || error("supply --output-dir, or use run_graphml_plots.py")
    output = abspath(options["--output-dir"])
    objectives = unique(Symbol.(split(options["--objectives"], ',')))
    parts = sort(unique(parse.(Int, split(options["--parts"], ','))))
    backend = Symbol(options["--backend"])
    println("READY"); flush(stdout)
    rows = build_comparisons(output; parts, objectives, backend)
    isempty(rows) && error("no feasible graph/part-count pairs; see status.tsv")
    write_metadata(output, parts, objectives, backend)
    for objective in objectives
        println("STAGE\tplot overview\t", objective); flush(stdout)
        plot_overview(rows, objective, output)
    end
    for collection in sort(unique(r.collection for r in rows); by=collection_order)
        println("STAGE\tplot bipartitions\t", collection); flush(stdout)
        plot_bipartitions(rows, collection, objectives, output)
    end
    println("RESULT\t", length(rows), " score rows; output: ", output)
    println("COMPLETE"); flush(stdout)
    return rows
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    GraphMLComparisonPlots.main()
end
