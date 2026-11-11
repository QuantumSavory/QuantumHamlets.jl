# Failed calls must validate before either sampling or invoking user callbacks.
function check_bury_argument_error_without_rng(f)
    Random.seed!(0x76616c)
    next_draw = rand(UInt64)
    Random.seed!(0x76616c)
    @test_throws ArgumentError f()
    @test rand(UInt64) == next_draw
end

@testset "BURY validation and range-search regressions" begin
    g = path_graph(4)
    empty = SimpleGraph(0)
    for backend in (:heap, :scan, :legacy)
        @testset "$backend" begin
            # Includes legacy's former infinite-loop targets and direct seed
            # helpers, which formerly sampled before rejecting bad targets.
            for target in (-1, 5, 2.0, 2.5, NaN, Inf, "2", big(typemax(Int)) + 1)
                for call in (
                        () -> Q.bury(g, target; backend),
                        () -> Q.bury_seeded(g, target, 1; backend),
                        () -> Q.bury_seeded(g, target; backend),
                        () -> Q.best_seed(g, target; backend),
                        () -> Q.best_sqrt_seed(g, target; backend))
                    check_bury_argument_error_without_rng(call)
                end
            end
            for seed in (0, 5, 1.0, "1")
                check_bury_argument_error_without_rng(
                    () -> Q.bury_seeded(g, 0, seed; backend))
            end
            for call in (
                    () -> Q.bury_seeded(empty, 0; backend),
                    () -> Q.bury_seeded(empty, 0, 1; backend),
                    () -> Q.best_seed(empty, 0; backend),
                    () -> Q.best_sqrt_seed(empty, 0; backend))
                check_bury_argument_error_without_rng(call)
            end
            @test Q.bury(empty, 0; backend) == (Dict{Any,Any}(), Bool[])
            @test Q.bury(empty; backend) == (Dict{Any,Any}(), Bool[])

            for graph in (g, empty), k in (-2, 0, 3, 5, 2.0, "2")
                for partitioner in (Q.bury_k, Q.bury_k_seeded)
                    check_bury_argument_error_without_rng(() -> partitioner(graph, k; backend))
                end
            end
            # Supported integer representations and target endpoints retain
            # the same assignments; zero target still needs a valid seed.
            for T in (Int8, UInt8, Int128, BigInt)
                @test Q.bury(g, T(2); backend) == Q.bury(g, 2; backend)
                @test Q.bury_seeded(g, T(0), T(1); backend) == Q.bury_seeded(g, 0, 1; backend)
                @test Q.bury_k(g, T(2); backend) == Q.bury_k(g, 2; backend)
                @test Q.bury_k_seeded(g, T(2); backend) == Q.bury_k_seeded(g, 2; backend)
            end

            complete = complete_graph(8)
            result = Q.bury_not_balanced_bipartition(complete, 4, 4;
                backend, cost_f=Q.naive_cost_of_partition)
            @test result == Q.bury(complete, 4; backend)
            @test first(Q.naive_cost_of_partition(Q.quantumLand(result[1], 2, 8), complete)) == 16

            calls = Int[]
            custom = function (graph, target)
                push!(calls, target)
                return Q.bury(graph, target; backend)
            end
            scores = function (land, graph)
                count_at_two = count(==(2), values(land.registry))
                return (count_at_two == 1 ? 101 : 100, nothing)
            end
            @test Q.bury_not_balanced_bipartition(g, 1, 3;
                backend, bury_v=custom, cost_f=scores) == Q.bury(g, 3; backend)
            @test calls == [1, 2, 3] # high scores still minimize; final tie wins

            empty!(calls)
            for (lo, hi) in ((2, 1), (-1, 2), (0, 5), (0.0, 2), (0, 2.0))
                check_bury_argument_error_without_rng(() ->
                    Q.bury_not_balanced_bipartition(g, lo, hi;
                        backend, bury_v=custom, cost_f=scores))
            end
            @test isempty(calls)
            @test Q.bury_not_balanced_bipartition(empty, 0, 0; backend) ==
                  (Dict{Any,Any}(), Bool[])
            check_bury_argument_error_without_rng(() ->
                Q.bury_not_balanced_bipartition(empty, 0, 0; backend, bury_v=Q.bury_seeded))
        end
    end
end
