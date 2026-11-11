# This module contains the pre-backend algorithms, not calls through the new
# dispatch wrappers. It remains an independent reference for every backend.
include("reference/bury_original.jl")

const BURY_BACKENDS = (:heap, :scan, :legacy)

function capture_bury_output(f)
    mktemp() do _, io
        result = redirect_stdout(f, io)
        flush(io)
        seekstart(io)
        return result, read(io, String)
    end
end

function capture_bury_rng(f, rngseed)
    Random.seed!(rngseed)
    result, output = capture_bury_output(f)
    return result, output, rand(UInt64)
end

function check_bury_result(actual, expected, target)
    @test actual == expected
    # Preserve the original public registry and color-vector types as well as
    # values; the standalone prototypes originally used a different Dict type.
    @test typeof(actual) === Tuple{Dict{Any,Any},Vector{Bool}}
    @test count(actual[2]) == target
end

function backend_disconnected_graph(n=12)
    g = SimpleGraph(n)
    for (u, v) in ((1, 2), (2, 3), (5, 6), (6, 7), (7, 8), (8, 5), (10, 11))
        v <= n && add_edge!(g, u, v)
    end
    return g
end

@testset "BURY backends against the frozen implementation" begin
    @testset "Every labeled simple graph on 1 to 5 vertices" begin
        for n in 1:5
            @testset "$n vertices" begin
                possible_edges = [(u, v) for u in 1:n for v in (u + 1):n]
                for mask in 0:((1 << length(possible_edges)) - 1)
                    g = SimpleGraph(n)
                    for (i, (u, v)) in enumerate(possible_edges)
                        iszero(mask & (1 << (i - 1))) || add_edge!(g, u, v)
                    end
                    original = copy(g)
                    for target in 0:n
                        expected = OriginalBury.bury(g, target)
                        for backend in BURY_BACKENDS
                            check_bury_result(Q.bury(g, target; backend), expected, target)
                        end
                        for seed in 1:n
                            expected = OriginalBury.bury_seeded(g, target, seed)
                            for backend in BURY_BACKENDS
                                check_bury_result(Q.bury_seeded(g, target, seed; backend), expected, target)
                            end
                        end
                    end
                    for k in 1:n
                        n % k == 0 || continue
                        for (partitioner, reference) in ((Q.bury_k, OriginalBury.bury_k),
                                                         (Q.bury_k_seeded, OriginalBury.bury_k_seeded))
                            expected = reference(g, k)
                            for backend in BURY_BACKENDS
                                actual = partitioner(g, k; backend)
                                @test actual == expected
                                @test actual isa Dict{Any,Any}
                            end
                        end
                    end
                    @test g == original
                end
            end
        end
    end

    @testset "Representative graphs, defaults, and boundary targets" begin
        rng = Xoshiro(0x42555259)
        fixtures = (path_graph(128), cycle_graph(120), star_graph(128),
                    SimpleGraph(64), grid([10, 12]), backend_disconnected_graph(),
                    random_regular_graph(128, 3; rng),
                    erdos_renyi(96, 0.4; rng), complete_graph(48), grid([32, 32]))
        for g in fixtures
            original = copy(g)
            n = nv(g)
            for target in unique([0, 1, n ÷ 4, n ÷ 2, 3n ÷ 4, n - 1, n])
                expected = OriginalBury.bury(g, target)
                check_bury_result(Q.bury(g, target), expected, target)
                for backend in BURY_BACKENDS
                    check_bury_result(Q.bury(g, target; backend), expected, target)
                    for seed in unique([1, n ÷ 2, n])
                        seeded_expected = OriginalBury.bury_seeded(g, target, seed)
                        check_bury_result(Q.bury_seeded(g, target, seed; backend), seeded_expected, target)
                    end
                end
            end
            @test Q.bury(g) == OriginalBury.bury(g)
            @test Q.bury_seeded(g, n ÷ 2, 1) == OriginalBury.bury_seeded(g, n ÷ 2, 1)
            @test g == original
        end
        empty = SimpleGraph(0)
        for backend in BURY_BACKENDS
            check_bury_result(Q.bury(empty, 0; backend), OriginalBury.bury(empty, 0), 0)
        end

        # These fixtures specifically exercise colored-but-unprocessed choices
        # and the seeded global fallback after exhausting a disconnected region.
        @test last(Q.bury(path_graph(6), 3)) == [true, true, true, false, false, false]
        for backend in BURY_BACKENDS
            @test Q.bury_seeded(SimpleGraph(6), 4, 6; backend) ==
                  OriginalBury.bury_seeded(SimpleGraph(6), 4, 6)
        end
    end

    @testset "Multipartition and simplified API routing" begin
        for g in (path_graph(12), grid([4, 4]), backend_disconnected_graph(),
                  complete_graph(12)), k in (1, 2, 3, 4, nv(g))
            nv(g) % k == 0 || continue
            original = copy(g)
            for (algorithm, partitioner, reference) in
                    ((:bury, Q.bury_k, OriginalBury.bury_k),
                     (:bury_seeded, Q.bury_k_seeded, OriginalBury.bury_k_seeded))
                expected = reference(g, k)
                @test partitioner(g, k) == expected
                @test Q.registry(Q.partition_graph(g, k; algorithm)) == expected
                for backend in BURY_BACKENDS
                    @test partitioner(g, k; backend) == expected
                    p = Q.partition_graph(g, k; algorithm, backend)
                    @test Q.registry(p) == expected
                    @test Q.partition_sizes(p) == fill(nv(g) ÷ k, k)
                end
            end
            @test g == original
        end
        @test Q._bury_color_kernel(:heap) === Q._bury_heap_colors
        @test Q._bury_color_kernel(:scan) === Q._bury_scan_colors
    end

    @testset "Verbose traces preserve public behavior" begin
        for g in (path_graph(6), complete_graph(6), star_graph(8),
                  backend_disconnected_graph(), SimpleGraph(5))
            for target in unique([0, 1, nv(g) ÷ 2, nv(g)])
                expected = capture_bury_output(() -> OriginalBury.bury(g, target; verbose=true))
                for backend in BURY_BACKENDS
                    @test capture_bury_output(() -> Q.bury(g, target; backend, verbose=true)) == expected
                end
                for seed in unique([1, nv(g)])
                    expected = capture_bury_output(() -> OriginalBury.bury_seeded(g, target, seed; verbose=true))
                    for backend in BURY_BACKENDS
                        @test capture_bury_output(() -> Q.bury_seeded(g, target, seed; backend, verbose=true)) == expected
                    end
                end
            end
        end
        g = path_graph(6)
        for (partitioner, reference) in ((Q.bury_k, OriginalBury.bury_k),
                                        (Q.bury_k_seeded, OriginalBury.bury_k_seeded))
            expected = capture_bury_output(() -> reference(g, 3; verbose=true))
            for backend in BURY_BACKENDS
                @test capture_bury_output(() -> partitioner(g, 3; backend, verbose=true)) == expected
            end
        end
    end

    @testset "Seed search, automatic defaults, and RNG consumption" begin
        for g in (path_graph(12), grid([4, 4]), backend_disconnected_graph())
            original = copy(g)
            target = nv(g) ÷ 2
            expected = capture_bury_output(() -> OriginalBury.best_seed(g, target; verbose=true))
            for backend in BURY_BACKENDS
                @test capture_bury_output(() -> Q.best_seed(g, target; backend, verbose=true)) == expected
            end
            for rngseed in (71, 991)
                sampled = capture_bury_rng(() -> OriginalBury.best_sqrt_seed(g, target; verbose=true), rngseed)
                automatic = capture_bury_rng(() -> OriginalBury.bury_seeded(g, target; verbose=true), rngseed)
                omitted_target = capture_bury_rng(() -> OriginalBury.bury_seeded(g), rngseed)
                @test capture_bury_rng(() -> Q.bury_seeded(g), rngseed) == omitted_target
                for backend in BURY_BACKENDS
                    @test capture_bury_rng(() -> Q.best_sqrt_seed(g, target; backend, verbose=true), rngseed) == sampled
                    @test capture_bury_rng(() -> Q.bury_seeded(g, target; backend, verbose=true), rngseed) == automatic
                    @test capture_bury_rng(() -> Q.bury_seeded(g; backend), rngseed) == omitted_target
                end
            end
            @test g == original
        end
    end

    @testset "Range search and two-argument callbacks" begin
        g = path_graph(8)
        original = copy(g)
        expected = capture_bury_output(() -> OriginalBury.bury_not_balanced_bipartition(g, 2, 4))
        @test capture_bury_output(() -> Q.bury_not_balanced_bipartition(g, 2, 4)) == expected
        for backend in BURY_BACKENDS
            @test capture_bury_output(() -> Q.bury_not_balanced_bipartition(g, 2, 4; backend)) == expected
            @test capture_bury_output(() -> Q.bury_not_balanced_bipartition(g, 2, 4; backend, bury_v=Q.bury)) == expected

            # Equal scores retain the last target, as in the original search.
            seen = Int[]
            custom = function (graph, target)
                @test graph === g
                push!(seen, target)
                return OriginalBury.bury(graph, target)
            end
            zero_cost = (_, _) -> (0, nothing)
            @test Q.bury_not_balanced_bipartition(g, 2, 4; backend, bury_v=custom, cost_f=zero_cost) ==
                  OriginalBury.bury(g, 4)
            @test seen == [2, 3, 4]

            # An explicitly supplied built-in seeded partitioner still honors
            # backend, while retaining its automatic-seed RNG stream and logs.
            expected_seeded = capture_bury_rng(
                () -> OriginalBury.bury_not_balanced_bipartition(g, 2, 4;
                      bury_v=OriginalBury.bury_seeded, cost_f=zero_cost), 101)
            @test capture_bury_rng(
                () -> Q.bury_not_balanced_bipartition(g, 2, 4;
                      backend, bury_v=Q.bury_seeded, cost_f=zero_cost), 101) == expected_seeded
        end
        @test g == original
    end

    @testset "Backend and input validation" begin
        g = path_graph(6)
        unknown_calls = (
            () -> Q.bury(g, 3; backend=:unknown),
            () -> Q.bury_seeded(g, 3, 1; backend=:unknown),
            () -> Q.bury_seeded(g, 3; backend=:unknown),
            () -> Q.bury_k(g, 2; backend=:unknown),
            () -> Q.bury_k_seeded(g, 2; backend=:unknown),
            () -> Q.best_seed(g, 3; backend=:unknown),
            () -> Q.best_sqrt_seed(g, 3; backend=:unknown),
            () -> Q.bury_not_balanced_bipartition(g, 2, 3; backend=:unknown),
            () -> Q.partition_graph(g, 2; backend=:unknown),
        )
        for f in unknown_calls
            Random.seed!(451)
            expected_next = rand(UInt64)
            Random.seed!(451)
            @test_throws ArgumentError f()
            @test rand(UInt64) == expected_next
        end
        for backend in BURY_BACKENDS
            for target in (-1, nv(g) + 1, 2.0, 2.5)
                @test_throws ArgumentError Q.bury(g, target; backend)
                @test_throws ArgumentError Q.bury_seeded(g, target, 1; backend)
                # Reject before automatic sampling consumes RNG state.
                Random.seed!(451)
                expected_next = rand(UInt64)
                Random.seed!(451)
                @test_throws ArgumentError Q.bury_seeded(g, target; backend)
                @test rand(UInt64) == expected_next
            end
            for seed in (0, nv(g) + 1, 1.0)
                @test_throws ArgumentError Q.bury_seeded(g, 3, seed; backend)
            end
        end
        # All backends now reject Float64 targets, even integral ones. The
        # frozen original remains unchanged; its permissive invalid-input
        # behavior is intentionally no longer a compatibility requirement.
        for backend in BURY_BACKENDS
            for f in (
                    () -> Q.bury(g, 2.0; backend),
                    () -> Q.bury_seeded(g, 2.0; backend),
                    () -> Q.best_seed(g, 2.0; backend),
                    () -> Q.best_sqrt_seed(g, 2.0; backend),
                    () -> Q.bury_not_balanced_bipartition(g, 2.0, 3.0; backend))
                @test_throws ArgumentError f()
            end
        end
        for backend in (:scan, :legacy), algorithm in (:random, :saran_vazirani)
            @test_throws ArgumentError Q.partition_graph(g, 2; algorithm, backend)
        end
        Random.seed!(120)
        expected = Q.k_partition_random(g, 2)
        Random.seed!(120)
        @test Q.registry(Q.partition_graph(g, 2; algorithm=:random, backend=:heap)) == expected
    end
end
