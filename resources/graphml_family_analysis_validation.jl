# Bounded data validation using tiny fixtures, not a full GraphML sweep.
# julia --startup-file=no --project=.. resources/graphml_family_analysis_validation.jl
using Test, Graphs, GraphIO.GraphML, EzXML, QuantumHamlets, Random, SHA, Metis
include("graphml_family_analysis.jl")
const FA = GraphMLFamilyAnalysis
const FQ = QuantumHamlets

function family_table(path)
    lines = readlines(path)
    header = split(first(lines), '\t')
    return [Dict(header .=> split(line, '\t'; keepempty=true)) for line in Iterators.drop(lines, 1)]
end

function explicit_sampled_partition(g, objective, seed; backend=:heap)
    legacy = objective === :matching ? FQ.matching_cost_of_partition : FQ.cut_rank_cost
    Random.seed!(seed)
    best = nothing
    score = nothing
    redirect_stdout(devnull) do
        for target in fld(nv(g), 2):cld(nv(g), 2)
            reg, _ = FQ.bury_seeded(g, target; backend)
            value = first(legacy(FQ.quantumLand(reg, 2, nv(g)), g))
            if best === nothing || value <= score
                best = copy(reg)
                score = value
            end
        end
    end
    return GraphPartition(best; nparts=2)
end

function synthetic_family_row(name, algorithm, cost; family="toffoli", parts=2,
                              objective=:matching, nonempty=true, balanced=true,
                              replicate=0)
    return (collection="graphml_reduced", graph=name, family=family, parameter=3,
        n=12, m=14, parts=parts, objective=objective, algorithm=algorithm,
        replicate=replicate, seed=missing, variant="synthetic", cost=cost,
        part_sizes=fill(12 ÷ parts, parts), matches_bury_balance=balanced, nonempty=nonempty)
end

@testset "GraphML family analysis validation" begin
    @testset "Classification precedes scores" begin
        for (collection, name, expected) in (
                ("graphml_initial", "barenco_tof_10", ("barenco_toffoli", 10)),
                ("graphml_reduced", "tof_3", ("toffoli", 3)),
                ("graphml_reduced", "gf2^5_mult", ("gf_multiplier", 5)),
                ("graphml_evolution", "k4", ("evolution", 4)),
                ("graphml_initial", "hwb6", ("hwb", 6)),
                ("graphml_initial", "mod_red_21", ("mod_red", 21)))
            @test FA.family_parameter(collection, name) == expected
        end
        @test isequal(FA.family_parameter("graphml_initial", "named_fixture"), ("named_fixture", missing))
    end

    @testset "Seeded objectives reset the same sampled matching search" begin
        g = path_graph(7)
        original = copy(g)
        for backend in (:heap, :scan), seed in (0, 9), objective in (:rank, :matching)
            expected = explicit_sampled_partition(g, objective, seed; backend)
            expected_next = rand(UInt64)
            # Deliberately perturb the incoming RNG. Each helper invocation
            # must start with the requested fixed seed independently.
            Random.seed!(123456)
            rand(19)
            actual = FA.sampled_bipartition(g, objective, seed; backend)
            @test membership(actual) == membership(expected)
            @test rand(UInt64) == expected_next
            @test g == original
        end
        @test_throws ArgumentError FA.sampled_bipartition(g, :unknown, 0)
    end

    @testset "Tiny collection, exact rows, and preservation" begin
        mktempdir() do scratch
            resources = joinpath(scratch, "resources")
            mkpath(resources)
            fixtures = Dict("tof_3" => path_graph(5), "gf2^4_mult" => path_graph(6),
                            "k1" => SimpleGraph(1))
            originals = Dict(name => copy(g) for (name, g) in fixtures)
            paths = String[]
            for (collection, name) in (("graphml_initial", "tof_3"),
                                      ("graphml_reduced", "gf2^4_mult"),
                                      ("graphml_evolution", "k1"))
                directory = joinpath(resources, collection)
                mkpath(directory)
                path = joinpath(directory, name * ".graphml")
                savegraph(path, fixtures[name], GraphMLFormat())
                push!(paths, path)
            end
            hashes = Dict(path => bytes2hex(open(sha256, path)) for path in paths)
            output = joinpath(scratch, "output")
            for arguments in ((parts=Int[],), (parts=[1],), (parts=[2.0],),
                              (parts=[2, "3"],), (seeds=Int[],), (seeds=[-1],),
                              (seeds=[0.0],), (seeds=[0, 0],),
                              (seeds=[big(typemax(Int32)) + 1],), (backend=:unknown,))
                Random.seed!(11)
                expected_next = rand(UInt64)
                Random.seed!(11)
                @test_throws ArgumentError FA.run_analysis(output; resources_dir=resources, arguments...)
                @test rand(UInt64) == expected_next
                @test !ispath(output)
            end
            options = copy(Metis.options)
            rows = FA.run_analysis(output; resources_dir=resources,
                parts=(3, 2, 2), seeds=(9, 0), backend=:scan)
            @test length(rows) == 28
            @test all(row -> propertynames(row) == FA.RESULT_FIELDS, rows)
            @test Metis.options == options
            for row in rows
                g = fixtures[row.graph]
                objective = FA.GCP.objective_function(row.objective)
                expected = if row.algorithm == "BURY"
                    @test row.replicate == 0 && ismissing(row.seed)
                    @test row.variant == (row.parts == 2 ? "ordinary_target_range" : "ordinary_equal_multipartition")
                    row.parts == 2 ? bipartition_graph(g; objective, backend=:scan) :
                                     partition_graph(g, row.parts; backend=:scan)
                elseif row.algorithm == "BURY_seeded" && row.parts == 2
                    @test row.seed == (9, 0)[row.replicate]
                    @test row.variant == "sampled_matching_inner_target_range"
                    explicit_sampled_partition(g, row.objective, row.seed; backend=:scan)
                elseif row.algorithm == "BURY_seeded"
                    @test row.replicate == 0 && ismissing(row.seed)
                    @test row.variant == "minimum_degree_residual_multipartition"
                    GraphPartition(FQ.bury_k_seeded(g, row.parts; backend=:scan); nparts=row.parts)
                else
                    @test row.algorithm == "METIS"
                    @test row.seed == (9, 0)[row.replicate]
                    @test row.variant == "kway_edgecut_ufactor30"
                    FA.GCP.metis_partition(g, row.parts; seed=row.seed)
                end
                @test row.cost == objective(g, expected)
                @test row.part_sizes == partition_sizes(expected)
                @test row.nonempty == all(>(0), row.part_sizes)
                @test row.matches_bury_balance == FA.GCP.balance_matches(expected, nv(g), row.parts)
                @test row.n == nv(g) && row.m == ne(g)
                @test fixtures[row.graph] == originals[row.graph]
            end
            @test count(row -> row.algorithm == "BURY_seeded" && row.parts == 3, rows) == 2
            @test count(row -> row.algorithm == "BURY_seeded" && row.parts == 2, rows) == 8
            @test length(family_table(joinpath(output, "results.tsv"))) == 28
            @test length(family_table(joinpath(output, "summary.tsv"))) == 6
            @test length(family_table(joinpath(output, "family_summary.tsv"))) == 16
            status = family_table(joinpath(output, "status.tsv"))
            @test length(status) == 6
            @test count(row -> row["status"] == "unsupported_balance", status) == 3
            @test all(row -> row["status"] == "unsupported_balance" || row["metis_total_replicates"] == "2", status)
            @test occursin("Fixed replicate seeds: 9,0", read(joinpath(output, "metadata.txt"), String))
            for path in paths
                @test bytes2hex(open(sha256, path)) == hashes[path]
            end
            before = Dict(name => read(joinpath(output, name)) for name in readdir(output))
            @test_throws ErrorException FA.run_analysis(output; resources_dir=resources, parts=(2,), seeds=(0,))
            @test before == Dict(name => read(joinpath(output, name)) for name in readdir(output))
        end
    end

    @testset "Summaries retain all cases and honest replicate denominators" begin
        rows = NamedTuple[]
        add(name, algorithm, cost; kwargs...) = push!(rows, synthetic_family_row(name, algorithm, cost; kwargs...))
        add("tof_3", "BURY", 8)
        add("tof_3", "BURY_seeded", 4; replicate=1)
        add("tof_3", "BURY_seeded", 8; replicate=2)
        add("tof_3", "METIS", 0; replicate=1, nonempty=false, balanced=false)
        add("tof_3", "METIS", 4; replicate=2, balanced=false)
        add("tof_3", "METIS", 10; replicate=3)
        for (algorithm, cost) in (("BURY", 10), ("BURY_seeded", 5), ("METIS", 10))
            add("tof_3", algorithm, cost; parts=3)
        end
        for (algorithm, cost) in (("BURY", 5), ("BURY_seeded", 7))
            add("tof_4", algorithm, cost)
        end
        add("tof_4", "METIS", 0; nonempty=false, balanced=false)
        for (algorithm, cost) in (("BURY", 4), ("BURY_seeded", 2), ("METIS", 4))
            add("tof_5", algorithm, cost; objective=:rank)
        end
        for (algorithm, cost) in (("BURY", 3), ("BURY_seeded", 5), ("METIS", 4))
            add("gf2^4_mult", algorithm, cost; family="gf_multiplier")
        end
        summary = FA.summarize_rows(rows)
        @test length(summary) == 5
        a = only(filter(r -> r.graph == "tof_3" && r.parts == 2, summary))
        @test (a.metis_total_replicates, a.metis_valid_replicates, a.metis_balanced_replicates) == (3, 2, 1)
        @test (a.metis_min, a.metis_median, a.metis_max) == (4, 7, 10)
        @test (a.metis_balanced_min, a.metis_balanced_median, a.metis_balanced_max) == (10, 10, 10)
        @test (a.seeded_min, a.seeded_median, a.seeded_max) == (4, 6, 8)
        @test a.bury_vs_metis_median == "higher"
        @test a.bury_vs_balanced_metis_median == "lower"
        @test a.seeded_vs_metis_median == a.seeded_vs_balanced_metis_median == "lower"
        b = only(filter(r -> r.graph == "tof_4", summary))
        @test (b.metis_total_replicates, b.metis_valid_replicates, b.metis_balanced_replicates) == (1, 0, 0)
        @test ismissing(b.metis_median) && ismissing(b.metis_balanced_median)
        @test b.bury_vs_metis_median == b.seeded_vs_metis_median == "unavailable"
        families = FA.summarize_families(summary)
        @test length(families) == 12
        @test Set(row.family for row in families) == Set(("toffoli", "gf_multiplier"))
        @test Set(row.objective for row in families) == Set((:matching, :rank))
        @test all(r -> r.lower + r.tied + r.higher + r.unavailable == r.cases, families)
        ordinary = only(filter(r -> r.family == "toffoli" && r.objective == :matching &&
            r.algorithm == "BURY" && r.metis_comparator == "all_nonempty", families))
        @test (ordinary.graphs, ordinary.cases) == (2, 3)
        @test (ordinary.lower, ordinary.tied, ordinary.higher, ordinary.unavailable) == (0, 1, 1, 1)
        balanced = only(filter(r -> r.family == "toffoli" && r.objective == :matching &&
            r.algorithm == "BURY" && r.metis_comparator == "bury_balance", families))
        @test (balanced.lower, balanced.tied, balanced.higher, balanced.unavailable) == (1, 1, 0, 1)
    end

    @testset "Actual empty METIS replicate stays in raw data only" begin
        source = joinpath(@__DIR__, "graphml_evolution", "k1.graphml")
        original_hash = bytes2hex(open(sha256, source))
        mktempdir() do scratch
            resources = joinpath(scratch, "resources")
            directory = joinpath(resources, "graphml_evolution")
            mkpath(directory)
            cp(source, joinpath(directory, "k1.graphml"))
            output = joinpath(scratch, "output")
            rows = FA.run_analysis(output; resources_dir=resources, parts=(5,), seeds=(0, 1))
            @test length(rows) == 8
            @test all(r -> !r.nonempty, filter(r -> r.algorithm == "METIS" && r.seed == 0, rows))
            @test all(r -> r.nonempty && r.matches_bury_balance, filter(r -> r.algorithm != "METIS", rows))
            for summary in FA.summarize_rows(rows)
                metis = filter(r -> r.algorithm == "METIS" && r.objective == summary.objective, rows)
                @test summary.metis_total_replicates == 2
                @test summary.metis_valid_replicates == count(r -> r.nonempty, metis)
                @test summary.metis_balanced_replicates == count(r -> r.nonempty && r.matches_bury_balance, metis)
            end
            @test only(family_table(joinpath(output, "status.tsv")))["status"] == "metis_empty_part"
        end
        @test bytes2hex(open(sha256, source)) == original_hash
    end
end
