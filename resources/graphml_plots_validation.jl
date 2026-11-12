# Focused, bounded validation; no research sweeps or source graph mutations.
# From the package root:
#   julia --startup-file=no --project=.. resources/graphml_plots_validation.jl
using Test, SHA, Statistics, Graphs, GraphIO.GraphML, EzXML, Metis, QuantumHamlets

include("graphml_plots.jl")
const GCP = GraphMLComparisonPlots
const GRAPHML_PACKAGE_ROOT = dirname(@__DIR__)

function read_inventory_tsv(path)
    lines = readlines(path)
    header = split(first(lines), '\t')
    return [Dict(header .=> split(line, '\t'; keepempty=true)) for line in Iterators.drop(lines, 1)]
end

function read_reference_topology(path)
    open(path) do io
        n, m = parse.(Int, split(readline(io), '\t'))
        g = SimpleGraph(n)
        for line in eachline(io)
            u, v = parse.(Int, split(line, '\t'))
            add_edge!(g, u, v) || error("duplicate reference edge in $path")
        end
        ne(g) == m || error("incomplete reference topology in $path")
        return g
    end
end

@testset "GraphML comparison plotting validation" begin
    inventory = read_inventory_tsv(joinpath(GRAPHML_PACKAGE_ROOT,
        "experiments", "graphml_bury", "inventory", "inventory.tsv"))
    inputs = GCP.graphml_inputs()
    historical_inputs = [joinpath(GRAPHML_PACKAGE_ROOT, row["path"]) for row in inventory]
    @test length(inputs) == 45
    @test length(unique(inputs)) == 45
    @test Set(abspath.(inputs)) == Set(historical_inputs)
    @test length(filter(path -> basename(dirname(path)) == "graphml_initial", inputs)) == 20
    @test length(filter(path -> basename(dirname(path)) == "graphml_reduced", inputs)) == 20
    @test length(filter(path -> basename(dirname(path)) == "graphml_evolution", inputs)) == 5

    @testset "All source hashes, topologies, and graph statistics" begin
        for row in inventory
            @testset "$(row["directory"])/$(row["graph"])" begin
                path = joinpath(GRAPHML_PACKAGE_ROOT, row["path"])
                @test filesize(path) == parse(Int, row["bytes"])
                @test bytes2hex(open(sha256, path)) == row["sha256"]
                g = loadgraph(path, GraphMLFormat())
                expected = read_reference_topology(joinpath(GRAPHML_PACKAGE_ROOT, row["topology_path"]))
                @test g == expected
                @test !is_directed(g)
                @test nv(g) == parse(Int, row["n"])
                @test ne(g) == parse(Int, row["m"])
                @test ne(g) == parse(Int, row["declared_edges"])
                @test density(g) ≈ parse(Float64, row["density"])
                degrees = degree(g)
                @test minimum(degrees) == parse(Int, row["degree_min"])
                @test median(degrees) == parse(Float64, row["degree_median"])
                @test mean(degrees) ≈ parse(Float64, row["degree_mean"])
                @test sort(degrees)[ceil(Int, 0.95nv(g))] == parse(Int, row["degree_p95"])
                @test maximum(degrees) == parse(Int, row["degree_max"])
                @test count(iszero, degrees) == parse(Int, row["isolated"])
                components = connected_components(g)
                @test length(components) == parse(Int, row["components"])
                @test maximum(length, components) == parse(Int, row["largest_component"])
                for k in 2:8
                    @test GCP.supported_parts(nv(g), k) == (k <= nv(g) && (k == 2 || rem(nv(g), k) == 0))
                end
            end
        end
    end

    @testset "Objective and balance policies" begin
        @test GCP.objective_function(:matching) === matching_cost
        @test GCP.objective_function(:rank) === pairwise_cut_rank
        @test_throws ArgumentError GCP.objective_function(:unknown)
        for n in 0:9, k in -1:11
            expected = n >= 2 && 2 <= k <= n && (k == 2 || rem(n, k) == 0)
            @test GCP.supported_parts(n, k) == expected
        end
        odd = SimpleGraph(5)
        add_edge!(odd, 1, 4)
        add_edge!(odd, 2, 3)
        for g in (path_graph(2), odd, path_graph(6)), objective in (:matching, :rank)
            original = copy(g)
            score = GCP.objective_function(objective)
            candidates = [GraphPartition(first(QuantumHamlets.bury(g, target)); nparts=2)
                          for target in fld(nv(g), 2):cld(nv(g), 2)]
            scores = [score(g, p) for p in candidates]
            expected = candidates[findlast(==(minimum(scores)), scores)]
            actual = bipartition_graph(g; objective=score)
            @test membership(actual) == membership(expected)
            @test sort(partition_sizes(actual)) == [fld(nv(g), 2), cld(nv(g), 2)]
            @test g == original
        end
        @test matching_cost(odd, bipartition_graph(odd)) == 0
    end

    @testset "METIS options and output ownership" begin
        original_options = copy(Metis.options)
        try
            seed_index = Int(Metis.METIS_OPTION_SEED) + 1
            imbalance_index = Int(Metis.METIS_OPTION_UFACTOR) + 1
            for g in (path_graph(5), path_graph(6), complete_graph(6)), k in (2, 3)
                original = copy(g)
                Metis.options[seed_index] = 17
                Metis.options[imbalance_index] = 7
                caller_options = copy(Metis.options)
                p = GCP.metis_partition(g, k; seed=0)
                @test p isa GraphPartition
                @test p.nparts == k
                @test length(p) == nv(g)
                @test all(label -> 1 <= label <= k, membership(p))
                @test sum(partition_sizes(p)) == nv(g)
                @test Metis.options == caller_options
                @test g == original
                again = GCP.metis_partition(g, k; seed=0)
                @test membership(again) == membership(p)
                @test Metis.options == caller_options
                labels = membership(p)
                fill!(labels, 1)
                @test membership(p) == membership(again)
            end
            before_error = copy(Metis.options)
            @test_throws InexactError GCP.metis_partition(path_graph(5), 2;
                seed=big(typemax(Int)) + 1)
            @test Metis.options == before_error
        finally
            copyto!(Metis.options, original_options)
        end
        @test Metis.options == original_options
    end

    @testset "Small comparison builder and output preservation" begin
        mktempdir() do scratch
            resources = joinpath(scratch, "resources")
            mkpath(resources)
            @test_throws ArgumentError GCP.graphml_inputs(resources)
            fixtures = Dict("odd" => path_graph(5), "even" => path_graph(6),
                            "singleton" => SimpleGraph(1))
            paths = String[]
            for (collection, name) in (("graphml_initial", "odd"),
                                       ("graphml_reduced", "even"),
                                       ("graphml_evolution", "singleton"))
                directory = joinpath(resources, collection)
                mkpath(directory)
                path = joinpath(directory, name * ".graphml")
                savegraph(path, fixtures[name], GraphMLFormat())
                push!(paths, path)
            end
            # Directories outside graphml* and non-GraphML files are excluded.
            mkpath(joinpath(resources, "other"))
            savegraph(joinpath(resources, "other", "ignored.graphml"), path_graph(3), GraphMLFormat())
            write(joinpath(resources, "graphml_initial", "ignored.txt"), "not a graph")
            @test Set(GCP.graphml_inputs(resources)) == Set(paths)
            hashes = Dict(path => bytes2hex(open(sha256, path)) for path in paths)
            output = joinpath(scratch, "output")
            for arguments in ((parts=Int[],), (parts=[1],), (parts=[2.0],),
                              (parts=[2, "3"],), (parts=[2, nothing],),
                              (objectives=Symbol[],), (objectives=[:unknown],),
                              (backend=:unknown,))
                @test_throws ArgumentError GCP.build_comparisons(output;
                    resources_dir=resources, arguments...)
            end
            @test !ispath(output)
            # Repeated objectives/part counts must not duplicate comparisons;
            # mixed-type invalid parts above must fail before sorting.
            rows = GCP.build_comparisons(output; resources_dir=resources,
                parts=(3, 2, 2), objectives=[:matching, :rank, :matching], backend=:scan)
            @test length(rows) == 12
            @test Set(row.graph for row in rows) == Set(("odd", "even"))
            @test all(row -> row.n == nv(fixtures[row.graph]) && row.m == ne(fixtures[row.graph]), rows)
            @test count(row -> row.graph == "odd", rows) == 4
            @test all(row -> row.graph != "odd" || row.k == 2, rows)
            for row in rows
                g = fixtures[row.graph]
                objective = GCP.objective_function(row.objective)
                expected = row.algorithm == "METIS" ? GCP.metis_partition(g, row.k) :
                    row.k == 2 ? bipartition_graph(g; backend=:scan, objective) :
                    partition_graph(g, row.k; backend=:scan)
                sizes = partition_sizes(expected)
                @test row.cost == objective(g, expected)
                @test row.sizes == sizes
                @test row.nonempty == all(>(0), sizes)
                @test row.imbalance == maximum(sizes) / (nv(g) / row.k)
                @test row.balanced == (all(>(0), sizes) &&
                    (row.k == 2 ? maximum(sizes) - minimum(sizes) <= 1 : all(==(nv(g) ÷ row.k), sizes)))
            end
            status = read_inventory_tsv(joinpath(output, "status.tsv"))
            @test length(status) == 6
            @test count(row -> row["status"] == "completed", status) == 3
            @test count(row -> row["status"] == "unsupported_balance", status) == 3
            result_table = read_inventory_tsv(joinpath(output, "results.tsv"))
            @test length(result_table) == length(rows)
            @test Set(row["objective"] for row in result_table) == Set(("matching", "rank"))
            @test all(row -> !isempty(row["part_sizes"]), result_table)
            @test all(row -> parse(Bool, row["nonempty"]) ==
                all(>(0), parse.(Int, split(row["part_sizes"], ','))), result_table)
            input_table = read_inventory_tsv(joinpath(output, "inputs.tsv"))
            @test length(input_table) == 3
            for row in input_table
                @test row["sha256"] == hashes[joinpath(resources, row["path"])]
            end
            for path in paths
                @test bytes2hex(open(sha256, path)) == hashes[path]
            end
            before = Dict(name => read(joinpath(output, name)) for name in readdir(output))
            @test_throws ErrorException GCP.build_comparisons(output; resources_dir=resources,
                parts=(2, 3), objectives=[:matching, :rank])
            @test before == Dict(name => read(joinpath(output, name)) for name in readdir(output))
        end
    end

    @testset "Real METIS empty-part result remains explicit" begin
        # Native KWAY with the documented seed/tolerance returns an empty part
        # on this supplied 25-vertex graph when five parts are requested. Keep
        # the raw result available for audit without losing valid BURY rows.
        source = joinpath(@__DIR__, "graphml_evolution", "k1.graphml")
        source_hash = bytes2hex(open(sha256, source))
        g = loadgraph(source, GraphMLFormat())
        original = copy(g)
        caller_options = copy(Metis.options)
        metis = GCP.metis_partition(g, 5; seed=0)
        @test nv(g) == 25
        @test any(iszero, partition_sizes(metis))
        @test !GCP.balance_matches(metis, nv(g), 5)
        @test Metis.options == caller_options
        @test g == original
        mktempdir() do scratch
            resources = joinpath(scratch, "resources")
            collection = joinpath(resources, "graphml_evolution")
            mkpath(collection)
            copy_path = joinpath(collection, "k1.graphml")
            cp(source, copy_path)
            output = joinpath(scratch, "output")
            rows = GCP.build_comparisons(output; resources_dir=resources,
                parts=(5,), objectives=[:matching, :rank])
            @test length(rows) == 4
            @test all(row -> row.graph == "k1" && row.k == 5, rows)
            for objective in (:matching, :rank)
                score = GCP.objective_function(objective)
                metis_row = only(filter(row -> row.algorithm == "METIS" && row.objective == objective, rows))
                @test !metis_row.nonempty
                @test !metis_row.balanced
                @test metis_row.sizes == partition_sizes(metis)
                @test metis_row.cost == score(g, metis)
                bury_row = only(filter(row -> row.algorithm == "BURY" && row.objective == objective, rows))
                @test bury_row.nonempty
                @test bury_row.balanced
                @test bury_row.sizes == fill(5, 5)
                @test bury_row.cost == score(g, partition_graph(g, 5))
            end
            status = only(read_inventory_tsv(joinpath(output, "status.tsv")))
            @test status["status"] == "metis_empty_part"
            @test occursin("excluded from plots", status["reason"])
            @test occursin(join(partition_sizes(metis), ','), status["reason"])
            result_table = read_inventory_tsv(joinpath(output, "results.tsv"))
            @test count(row -> row["nonempty"] == "false", result_table) == 2
            @test count(row -> row["algorithm"] == "BURY" && row["nonempty"] == "true", result_table) == 2
            @test bytes2hex(open(sha256, copy_path)) == source_hash
            @test Metis.options == caller_options
        end
        @test bytes2hex(open(sha256, source)) == source_hash
    end
end
