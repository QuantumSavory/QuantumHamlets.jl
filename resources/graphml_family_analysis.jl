"""Fixed-replicate, three-algorithm GraphML data collection; no plotting or sweep on include."""
module GraphMLFamilyAnalysis

using QuantumHamlets, Graphs, GraphIO.GraphML, EzXML, Random, Statistics, SHA, Dates
import Metis
include("graphml_plots.jl")
const GCP = GraphMLComparisonPlots
const Q = QuantumHamlets
const OBJECTIVES = (:matching, :rank)
const RESULT_FIELDS = (:collection, :graph, :family, :parameter, :n, :m, :parts,
    :objective, :algorithm, :replicate, :seed, :variant, :cost, :part_sizes,
    :matches_bury_balance, :nonempty)

"""Classify every filename before examining scores; return `(family, parameter)`."""
function family_parameter(collection, name)
    if collection == "graphml_evolution"
        found = match(r"^k(\d+)$", name)
        isnothing(found) || return ("evolution", parse(Int, only(found.captures)))
    end
    for (pattern, family) in ((r"^gf2\^(\d+)_mult$", "gf_multiplier"),
                              (r"^barenco_tof_(\d+)$", "barenco_toffoli"),
                              (r"^tof_(\d+)$", "toffoli"))
        found = match(pattern, name)
        isnothing(found) || return (family, parse(Int, only(found.captures)))
    end
    found = match(r"^(.*?)[_]?([0-9]+)$", name)
    isnothing(found) && return (String(name), missing)
    return (rstrip(found.captures[1], '_'), parse(Int, found.captures[2]))
end

tsv_value(value::Missing) = ""
tsv_value(value::AbstractVector) = join(value, ',')
tsv_value(value) = string(value)
function write_row(io, row)
    println(io, join(tsv_value.(collect(row)), '\t'))
    flush(io)
end

function write_table(path, rows)
    isempty(rows) && error("cannot infer empty table schema for $path")
    open(path, "w") do io
        write_row(io, propertynames(first(rows)))
        foreach(row -> write_row(io, values(row)), rows)
    end
end

function stage(collection, graph, parts, objective, algorithm, replicate)
    println("STAGE\t", join((collection, graph, parts, objective, algorithm, replicate), '\t'))
    flush(stdout)
end

function result_row(collection, name, family, parameter, g, parts, objective,
                    algorithm, replicate, seed, variant, p)
    sizes = partition_sizes(p)
    return (collection=collection, graph=name, family=family, parameter=parameter,
            n=nv(g), m=ne(g), parts=parts, objective=objective, algorithm=algorithm,
            replicate=replicate, seed=seed, variant=variant,
            cost=GCP.objective_function(objective)(g, p), part_sizes=sizes,
            matches_bury_balance=GCP.balance_matches(p, nv(g), parts),
            nonempty=all(>(0), sizes))
end

"""Run established sampled inner matching selection and outer target scoring.

Reset the default RNG to `seed` before each objective. This preserves the
existing sampled-seed policy while making objective comparisons reproducible.
Only established improvement-print output is suppressed; exceptions propagate.
"""
function sampled_bipartition(g, objective, seed; backend=:heap)
    legacy_cost = objective === :matching ? Q.matching_cost_of_partition :
                  objective === :rank ? Q.cut_rank_cost :
                  throw(ArgumentError("objective must be matching or rank"))
    Random.seed!(seed)
    reg, _ = redirect_stdout(devnull) do
        Q.bury_not_balanced_bipartition(g, fld(nv(g), 2), cld(nv(g), 2);
            bury_v=Q.bury_seeded, cost_f=legacy_cost, backend)
    end
    return GraphPartition(reg; nparts=2)
end

function score_summary(rows; balanced_only=false)
    selected = filter(row -> row.nonempty && (!balanced_only || row.matches_bury_balance), rows)
    costs = [row.cost for row in selected]
    return (median=isempty(costs) ? missing : median(costs),
            minimum=isempty(costs) ? missing : minimum(costs),
            maximum=isempty(costs) ? missing : maximum(costs),
            valid=length(costs), total=length(rows),
            balanced=count(row -> row.nonempty && row.matches_bury_balance, rows))
end

comparison(a, b) = ismissing(a) || ismissing(b) ? "unavailable" :
                    a < b ? "lower" : a > b ? "higher" : "tied"

"""Summarize every graph/part-count/objective case without selecting successes.

Raw empty-part outcomes are excluded from cost statistics, and total/valid/
balanced replicate counts remain explicit. Both all-nonempty and balance-matched
METIS statistics are retained. Deterministic variants have one observation.
"""
function summarize_rows(rows)
    groups = Dict{Tuple{String,String,Int,Symbol},Vector{NamedTuple}}()
    for row in rows
        key = (row.collection, row.graph, row.parts, row.objective)
        push!(get!(Vector{NamedTuple}, groups, key), row)
    end
    summary = NamedTuple[]
    for key in sort(collect(keys(groups)); by=key -> (key[1], key[2], key[3], string(key[4])))
        group = groups[key]
        base = first(group)
        ordinary = only(filter(row -> row.algorithm == "BURY", group))
        seeded = score_summary(filter(row -> row.algorithm == "BURY_seeded", group))
        metis_rows = filter(row -> row.algorithm == "METIS", group)
        metis = score_summary(metis_rows)
        balanced_metis = score_summary(metis_rows; balanced_only=true)
        bury_cost = ordinary.nonempty ? ordinary.cost : missing
        push!(summary, (collection=base.collection, graph=base.graph, family=base.family,
            parameter=base.parameter, n=base.n, m=base.m, parts=base.parts,
            objective=base.objective, bury_cost=bury_cost,
            seeded_median=seeded.median, seeded_min=seeded.minimum, seeded_max=seeded.maximum,
            seeded_total_replicates=seeded.total, seeded_valid_replicates=seeded.valid,
            seeded_balanced_replicates=seeded.balanced,
            metis_median=metis.median, metis_min=metis.minimum, metis_max=metis.maximum,
            metis_total_replicates=metis.total, metis_valid_replicates=metis.valid,
            metis_balanced_replicates=metis.balanced,
            metis_balanced_median=balanced_metis.median,
            metis_balanced_min=balanced_metis.minimum, metis_balanced_max=balanced_metis.maximum,
            bury_vs_metis_median=comparison(bury_cost, metis.median),
            seeded_vs_metis_median=comparison(seeded.median, metis.median),
            bury_vs_balanced_metis_median=comparison(bury_cost, balanced_metis.median),
            seeded_vs_balanced_metis_median=comparison(seeded.median, balanced_metis.median)))
    end
    return summary
end

"""Count lower/tied/higher/unavailable cases for every family and collection.

Cases are graph/part-count combinations, not independent statistical trials.
This table includes all classified families and both objectives.
"""
function summarize_families(summary)
    keys = unique((row.collection, row.family, row.objective) for row in summary)
    families = NamedTuple[]
    for key in sort(collect(keys); by=key -> (key[1], key[2], string(key[3])))
        selected = filter(row -> (row.collection, row.family, row.objective) == key, summary)
        for (comparator, field) in (("all_nonempty", :bury_vs_metis_median),
                                    ("bury_balance", :bury_vs_balanced_metis_median)),
            (algorithm, selected_field) in (("BURY", field),
                ("BURY_seeded", comparator == "all_nonempty" ?
                 :seeded_vs_metis_median : :seeded_vs_balanced_metis_median))
            outcomes = [getproperty(row, selected_field) for row in selected]
            push!(families, (collection=key[1], family=key[2], objective=key[3],
                algorithm=algorithm, metis_comparator=comparator,
                graphs=length(unique(row.graph for row in selected)), cases=length(selected),
                lower=count(==("lower"), outcomes), tied=count(==("tied"), outcomes),
                higher=count(==("higher"), outcomes), unavailable=count(==("unavailable"), outcomes)))
        end
    end
    return families
end

function write_metadata(output, parts, seeds, backend)
    open(joinpath(output, "metadata.txt"), "w") do io
        println(io, "UTC start: ", Dates.now(Dates.UTC), "\nJulia: ", VERSION,
                "\nThreads: ", Threads.nthreads())
        for package in (QuantumHamlets, Graphs, Metis)
            println(io, nameof(package), ": ", Base.pkgversion(package))
        end
        println(io, "Parts: ", join(parts, ','), "\nFixed replicate seeds: ", join(seeds, ','),
                "\nBURY backend: ", backend, "\nObjectives: matching,rank")
        println(io, "Ordinary k=2: floor/ceil target search selected separately by each objective; last-target ties.")
        println(io, "Seeded k=2: existing sampled matching-selected inner seeds, outer matching/rank target selection; RNG reset per objective and replicate.")
        println(io, "Seeded k>=3: deterministic minimum-degree residual seeds, one observation; not sampled seed search.")
        println(io, "METIS: KWAY edge-cut, native fixed seeds, equal target weights, UFACTOR=30 (3%); each membership reused for both scores.")
        println(io, "Deterministic rows: replicate=0, empty seed; sampled rows: replicate=1:length(seeds), seed=specified seed.")
        println(io, "Family classification uses filenames before observing costs; all supported cases retained.")
        println(io, "Empty-part raw scores retained but excluded from summaries; balance-matched METIS summaries also reported.")
        println(io, "Ranges are observed replicate min/max, not confidence intervals; family cases share graphs and are not independent trials.")
        println(io, "GraphML: simple undirected unweighted topology; node phase/type/layout and edge-type attributes ignored.")
        println(io, "Established seeded improvement messages suppressed. Timing includes scoring and data writing; not a partitioner benchmark.")
        for source in ("graphml_family_analysis.jl", "graphml_plots.jl", "../src/bury.jl", "../src/api.jl",
                       "../src/costs.jl", "../src/bury/heap.jl", "../src/bury/scan.jl", "../src/bury/legacy.jl")
            println(io, "sha256\t", source, '\t', bytes2hex(sha256(read(joinpath(@__DIR__, source)))))
        end
    end
end

"""
    run_analysis(output; parts=2:8, seeds=0:9, resources_dir=@__DIR__, backend=:heap)

Collect all feasible GraphML cases using fixed seeds chosen before scoring.
Return `Vector{NamedTuple}` with fields `RESULT_FIELDS`; `objective` is a Symbol,
`algorithm` is "BURY", "BURY_seeded", or "METIS", `part_sizes` is `Vector{Int}`,
and `parameter`/`seed` can be `missing`. The TSV uses empty cells for missing.
Deterministic rows have replicate zero; sampled rows use 1-based replicate IDs.

Write results, inputs, status, per-case summary, full-family summary, and
metadata without overwriting existing artifacts. Emits STAGE markers for an
external watchdog. The caller emits READY/COMPLETE and owns plotting/output.log.
"""
function run_analysis(output; parts=2:8, seeds=0:9, resources_dir=@__DIR__, backend=:heap)
    parts, seeds = collect(parts), collect(seeds)
    !isempty(parts) && all(k -> k isa Integer && 2 <= k <= typemax(Int), parts) ||
        throw(ArgumentError("parts must be nonempty integers at least two"))
    !isempty(seeds) && all(seed -> seed isa Integer && 0 <= seed <= typemax(Int32), seeds) ||
        throw(ArgumentError("seeds must be nonempty nonnegative Int32-compatible integers"))
    length(unique(seeds)) == length(seeds) || throw(ArgumentError("replicate seeds must be distinct"))
    parts, seeds = sort(unique(Int.(parts))), Int.(seeds)
    backend in (:heap, :scan, :legacy) || throw(ArgumentError("invalid BURY backend"))
    files = GCP.graphml_inputs(resources_dir)
    artifacts = ("results.tsv", "status.tsv", "inputs.tsv", "summary.tsv", "family_summary.tsv", "metadata.txt")
    for name in artifacts
        ispath(joinpath(output, name)) && error("refusing to overwrite $name")
    end
    mkpath(output)
    started = time_ns()
    write_metadata(output, parts, seeds, backend)
    rows = NamedTuple[]
    open(joinpath(output, "results.tsv"), "w") do results
        open(joinpath(output, "status.tsv"), "w") do status
            open(joinpath(output, "inputs.tsv"), "w") do inputs
                write_row(results, RESULT_FIELDS)
                write_row(status, (:collection, :graph, :n, :parts, :status,
                                   :metis_valid_replicates, :metis_total_replicates, :reason))
                write_row(inputs, (:collection, :graph, :family, :parameter, :n, :m, :sha256, :path))
                for path in files
                    collection, name = basename(dirname(path)), splitext(basename(path))[1]
                    family, parameter = family_parameter(collection, name)
                    stage(collection, name, "", "load", "", "")
                    g = loadgraph(path, GraphMLFormat())
                    g isa Graphs.Graph || error("expected simple undirected graph: $path")
                    original = copy(g)
                    n = nv(g)
                    write_row(inputs, (collection, name, family, parameter, n, ne(g),
                                       bytes2hex(sha256(read(path))), relpath(path, resources_dir)))
                    for k in parts
                        if !GCP.supported_parts(n, k)
                            write_row(status, (collection, name, n, k, "unsupported_balance", 0, 0,
                                k > n ? "more parts than vertices" : "exact multipartition requires n divisible by k"))
                            continue
                        end
                        ordinary = seeded_deterministic = nothing
                        if k >= 3
                            stage(collection, name, k, "both", "BURY", 0)
                            ordinary = partition_graph(g, k; backend)
                            stage(collection, name, k, "both", "BURY_seeded", 0)
                            seeded_deterministic = partition_graph(g, k; algorithm=:bury_seeded, backend)
                        end
                        record = function (objective, algorithm, replicate, seed, variant, p)
                            length(p) == n || error("$algorithm returned incomplete membership")
                            algorithm == "METIS" || GCP.balance_matches(p, n, k) ||
                                error("$algorithm violated BURY balance")
                            row = result_row(collection, name, family, parameter, g, k,
                                             objective, algorithm, replicate, seed, variant, p)
                            push!(rows, row)
                            write_row(results, values(row))
                        end
                        for objective in OBJECTIVES
                            stage(collection, name, k, objective, "BURY", 0)
                            p = k == 2 ? bipartition_graph(g; objective=GCP.objective_function(objective), backend) : ordinary
                            record(objective, "BURY", 0, missing,
                                   k == 2 ? "ordinary_target_range" : "ordinary_equal_multipartition", p)
                            if k == 2
                                for (replicate, seed) in enumerate(seeds)
                                    stage(collection, name, k, objective, "BURY_seeded", replicate)
                                    p = sampled_bipartition(g, objective, seed; backend)
                                    record(objective, "BURY_seeded", replicate, seed,
                                           "sampled_matching_inner_target_range", p)
                                end
                            else
                                record(objective, "BURY_seeded", 0, missing,
                                       "minimum_degree_residual_multipartition", seeded_deterministic)
                            end
                        end
                        invalid_seeds = Int[]
                        for (replicate, seed) in enumerate(seeds)
                            stage(collection, name, k, "both", "METIS", replicate)
                            p = GCP.metis_partition(g, k; seed)
                            all(>(0), partition_sizes(p)) || push!(invalid_seeds, seed)
                            for objective in OBJECTIVES
                                record(objective, "METIS", replicate, seed, "kway_edgecut_ufactor30", p)
                            end
                        end
                        g == original || error("comparison modified graph $path")
                        write_row(status, (collection, name, n, k,
                            isempty(invalid_seeds) ? "completed" : "metis_empty_part",
                            length(seeds) - length(invalid_seeds), length(seeds),
                            isempty(invalid_seeds) ? "" : "empty parts at native seeds=" * join(invalid_seeds, ',')))
                    end
                end
            end
        end
    end
    isempty(rows) && error("no feasible graph/part-count pairs; see status.tsv")
    stage("all", "all", "", "summary", "", "")
    summary = summarize_rows(rows)
    write_table(joinpath(output, "summary.tsv"), summary)
    write_table(joinpath(output, "family_summary.tsv"), summarize_families(summary))
    open(joinpath(output, "metadata.txt"), "a") do io
        println(io, "Data collection complete: true\nResult rows: ", length(rows),
                "\nPer-case summary rows: ", length(summary),
                "\nElapsed data collection seconds (not benchmark): ", (time_ns() - started) / 1e9)
    end
    return rows
end

end # module GraphMLFamilyAnalysis
