"""Family-focused plots of BURY, seeded BURY, and METIS.

Use run_graphml_family_plots.py. --input-dir rerenders saved data without
rerunning partitioners. Existing figures and measurements are never overwritten.
"""
module GraphMLFamilyPlots

using CairoMakie, Statistics, SHA
include("graphml_family_analysis.jl")
using .GraphMLFamilyAnalysis

const METHODS = ("BURY", "BURY_seeded", "METIS")
const METHOD_NAMES = ("BURY", "Seeded BURY", "METIS")
const COLORS = ("#4169B1", "#359887", "#D98C3B")
const INK = "#24313C"
const MUTED = "#64727B"
const METRICS = Dict("matching" => "Matching cost", "rank" => "Pairwise cut rank (GF(2))")

function read_results(input)
    lines = readlines(joinpath(input, "results.tsv"))
    header = split(first(lines), '\t')
    return [Dict(header .=> split(line, '\t'; keepempty=true)) for line in lines[2:end]]
end

function case_stats(rows, collection, graph, k, objective, method)
    selected = filter(r -> r["collection"] == collection && r["graph"] == graph &&
        parse(Int, r["parts"]) == k && r["objective"] == objective && r["algorithm"] == method, rows)
    valid = filter(r -> r["nonempty"] == "true", selected)
    values = parse.(Float64, [r["cost"] for r in valid])
    balanced = count(r -> r["matches_bury_balance"] == "true", valid)
    return (center=isempty(values) ? NaN : median(values),
            low=isempty(values) ? NaN : minimum(values),
            high=isempty(values) ? NaN : maximum(values),
            valid=length(valid), total=length(selected), balanced=balanced)
end

function cases_for(rows; collection, family, k=nothing, multipart=false)
    selected = filter(r -> r["collection"] == collection && r["family"] == family &&
        r["algorithm"] == "BURY" && r["objective"] == "matching" &&
        (isnothing(k) || parse(Int, r["parts"]) == k) &&
        (!multipart || parse(Int, r["parts"]) >= 3), rows)
    cases = [(graph=r["graph"], n=parse(Int, r["n"]), k=parse(Int, r["parts"]),
              parameter=parse(Int, r["parameter"])) for r in selected]
    return sort(unique(cases); by=r -> (r.parameter, r.k))
end

function new_figure(title, subtitle; width=1120, height=770, columns=2)
    fig = Figure(size=(width, height), fontsize=16, backgroundcolor=:white,
                 figure_padding=(25, 25, 20, 18))
    Label(fig[0, 1:columns], title; fontsize=27, font=:bold, color=INK, halign=:left)
    Label(fig[1, 1:columns], subtitle; fontsize=15, color=MUTED, halign=:left)
    elements = [PolyElement(color=color) for color in COLORS]
    Legend(fig[2, 1:columns], elements, collect(METHOD_NAMES); orientation=:horizontal,
           framevisible=false, labelsize=15, halign=:left, patchsize=(19, 12),
           colgap=26, padding=(0, 0, 2, 8))
    rowgap!(fig.layout, 13)
    colgap!(fig.layout, 34)
    return fig
end

format_score(value) = isinteger(value) ? string(Int(value)) : string(round(value; digits=1))

function bar_panel!(position, rows, collection, cases, objective;
                    title, ticklabels, xlabel="", annotate=true, ymax=nothing)
    ax = Axis(position; title, titlealign=:left, titlesize=17, titlecolor=INK,
              xlabel, ylabel=METRICS[objective], xticks=(collect(eachindex(cases)), ticklabels),
              xticklabelsize=13, yticklabelsize=13, xgridvisible=false,
              ygridcolor=("#95A3AA", 0.16), xgridwidth=0.7, ygridwidth=0.7,
              xticksvisible=false, yticksvisible=false, xlabelsize=14, ylabelsize=14,
              spinewidth=0.7, bottomspinecolor="#C7D0D5")
    hidespines!(ax, :l, :r, :t)
    top = 0.0
    for case in cases, method in METHODS
        stat = case_stats(rows, collection, case.graph, case.k, objective, method)
        isnan(stat.high) || (top = max(top, stat.high))
    end
    top = max(top, 1.0)
    limit = isnothing(ymax) ? top * 1.23 : ymax
    for (i, case) in enumerate(cases), (j, method) in enumerate(METHODS)
        s = case_stats(rows, collection, case.graph, case.k, objective, method)
        x = i + (j - 2) * 0.235
        if s.valid == 0
            text!(ax, x, limit * 0.04; text="No valid\npartition", align=(:center, :bottom),
                  fontsize=10, color=MUTED, rotation=pi/2)
            continue
        end
        barplot!(ax, [x], [s.center]; color=COLORS[j], width=0.205)
        if s.high > s.low
            errorbars!(ax, [x], [s.center], [s.center - s.low], [s.high - s.center];
                       color=INK, linewidth=1.15, whiskerwidth=8)
        end
        if annotate
            suffix = method == "METIS" && s.balanced < s.valid ? "†" : ""
            s.valid < s.total && (suffix *= "‡")
            text!(ax, x, s.high + limit * 0.027; text=format_score(s.center) * suffix,
                  fontsize=11, align=(:center, :bottom), color=INK)
        end
    end
    xlims!(ax, 0.48, length(cases) + 0.52)
    ylims!(ax, 0, limit)
    return ax
end

function footer!(fig, row, columns, rows; sampled)
    runs = length(unique(r["seed"] for r in rows if r["algorithm"] == "METIS"))
    description = sampled ? "Bars: median cost. Whiskers: observed min–max across up to $runs valid runs; ordinary BURY is deterministic." :
        "BURY methods are deterministic here. METIS: median and observed min–max across up to $runs valid runs."
    Label(fig[row, 1:columns], description; fontsize=12, color=MUTED, halign=:left)
    Label(fig[row + 1, 1:columns], "† Some valid METIS runs exceed BURY’s balance bound.  ‡ Some METIS runs returned an empty part and were excluded.";
          fontsize=11, color=MUTED, halign=:left)
end

function save_figure(output, stem, fig)
    for extension in ("png", "svg")
        path = joinpath(output, stem * "." * extension)
        ispath(path) && error("refusing to overwrite $path")
        extension == "png" ? save(path, fig; px_per_unit=1.7) : save(path, fig)
    end
end

function plot_multipliers_bipartitions(rows, output)
    fig = new_figure("GF(2ᵐ) multipliers", "Two parts · lower cost is better · seeded BURY uses matching-selected samples")
    for (i, (collection, title)) in enumerate((("graphml_initial", "Initial"), ("graphml_reduced", "Reduced")))
        cases = cases_for(rows; collection, family="gf_multiplier", k=2)
        for (j, objective) in enumerate(("matching", "rank"))
            bar_panel!(fig[i + 2, j], rows, collection, cases, objective;
                title=title * " · " * METRICS[objective],
                ticklabels=["m = $(c.parameter)\n$(c.n) vertices" for c in cases])
        end
    end
    footer!(fig, 5, 2, rows; sampled=true)
    save_figure(output, "multipliers_bipartitions", fig)
end

function plot_multipliers_multipartitions(rows, objective, output)
    fig = new_figure("GF(2ᵐ) multipliers · " * METRICS[objective],
        "Feasible counts within 3:8 parts · seeded BURY uses minimum-degree seeding · lower cost is better";
        width=1330, height=800, columns=3)
    for (i, (collection, title)) in enumerate((("graphml_initial", "Initial"), ("graphml_reduced", "Reduced")))
        allcases = cases_for(rows; collection, family="gf_multiplier", multipart=true)
        for (j, exponent) in enumerate(4:6)
            cases = filter(c -> c.parameter == exponent, allcases)
            isempty(cases) && continue
            bar_panel!(fig[i + 2, j], rows, collection, cases, objective;
                title="$title · m = $exponent · $(first(cases).n) vertices",
                ticklabels=["$(c.k) parts" for c in cases])
        end
    end
    footer!(fig, 5, 3, rows; sampled=false)
    save_figure(output, "multipliers_multipartitions_" * objective, fig)
end

function plot_toffoli(rows, objective, output)
    fig = new_figure("Toffoli families · " * METRICS[objective],
        "Two parts · all available sizes · seeded BURY uses matching-selected samples"; width=1220, height=790)
    for (i, (collection, title)) in enumerate((("graphml_initial", "Initial"), ("graphml_reduced", "Reduced")))
        for (j, (family, name)) in enumerate((("toffoli", "Toffoli"), ("barenco_toffoli", "Barenco Toffoli")))
            cases = cases_for(rows; collection, family, k=2)
            bar_panel!(fig[i + 2, j], rows, collection, cases, objective;
                title="$title · $name", ticklabels=["$(c.parameter)\n$(c.n) vertices" for c in cases],
                xlabel="Circuit parameter")
        end
    end
    footer!(fig, 5, 2, rows; sampled=true)
    save_figure(output, "toffoli_bipartitions_" * objective, fig)
end

function plot_evolution(rows, output)
    cases = cases_for(rows; collection="graphml_evolution", family="evolution", k=2)
    fig = new_figure("Evolution family", "Two parts · all five stages · seeded BURY uses matching-selected samples"; width=1250, height=510)
    for (j, objective) in enumerate(("matching", "rank"))
        bar_panel!(fig[3, j], rows, "graphml_evolution", cases, objective; title=METRICS[objective],
                   ticklabels=["Stage $(c.parameter)\n$(c.n) vertices" for c in cases])
    end
    footer!(fig, 4, 2, rows; sampled=true)
    save_figure(output, "evolution_bipartitions", fig)
end

function render_all(input, output)
    rows = read_results(input)
    isempty(rows) && error("empty results table")
    println("STAGE\tplot\tmultiplier bipartitions"); flush(stdout)
    plot_multipliers_bipartitions(rows, output)
    for objective in ("matching", "rank")
        println("STAGE\tplot\tmultiplier multipartitions\t", objective); flush(stdout)
        plot_multipliers_multipartitions(rows, objective, output)
        println("STAGE\tplot\tToffoli families\t", objective); flush(stdout)
        plot_toffoli(rows, objective, output)
    end
    println("STAGE\tplot\tevolution"); flush(stdout)
    plot_evolution(rows, output)
    open(joinpath(output, "figure_provenance.txt"), "w") do io
        println(io, "Input results: ", abspath(joinpath(input, "results.tsv")))
        println(io, "Input SHA-256: ", bytes2hex(sha256(read(joinpath(input, "results.tsv")))))
        println(io, "Plot source SHA-256: ", bytes2hex(sha256(read(@__FILE__))))
        println(io, "Figures use all available sizes within their named families; full tables cover every input family.")
        println(io, "Bars are medians of valid runs; whiskers are observed min/max, not confidence intervals.")
        println(io, "Bipartition seeded BURY samples seeds and selects internally by matching, even for rank panels.")
        println(io, "Multipartition seeded BURY uses the deterministic minimum-degree residual seeder.")
    end
    return rows
end

function main(args=ARGS)
    options = Dict("--backend" => "heap", "--parts" => "2,3,4,5,6,7,8", "--seeds" => "0,1,2,3,4,5,6,7,8,9")
    iseven(length(args)) || error("expected --name value pairs")
    for i in 1:2:length(args)
        args[i] in ("--output-dir", "--input-dir", "--backend", "--parts", "--seeds") || error("unknown option $(args[i])")
        options[args[i]] = args[i + 1]
    end
    output = abspath(options["--output-dir"])
    println("READY"); flush(stdout)
    input = if haskey(options, "--input-dir")
        abspath(options["--input-dir"])
    else
        GraphMLFamilyAnalysis.run_analysis(output;
            parts=parse.(Int, split(options["--parts"], ',')),
            seeds=parse.(Int, split(options["--seeds"], ',')), backend=Symbol(options["--backend"]))
        output
    end
    render_all(input, output)
    println("COMPLETE"); flush(stdout)
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    GraphMLFamilyPlots.main()
end
