# Run from the package root:
# julia --startup-file=no --project=. experiments/bury_performance/runtests.jl
#
# This standalone suite compares experimental code with the library's legacy backend.
# It never edits or includes changes in the package's production source.
using Test
using Graphs
using Random
using QuantumHamlets

include("scan.jl")
include("heap.jl")
include("multi.jl")
include("search.jl")

const Q = QuantumHamlets
const SINGLE_IMPLEMENTATIONS = ((:scan, BuryScan.bury_scan),
                                (:heap, BuryHeap.bury_heap))
const COLOR_KERNELS = ((:scan, BuryScan._colors), (:heap, BuryHeap._colors))

function check_single(g, target, seed, context, counters)
    expected = isnothing(seed) ? Q.bury(g, target; backend=:legacy) : Q.bury_seeded(g, target, seed; backend=:legacy)
    for (name, implementation) in SINGLE_IMPLEMENTATIONS
        actual = implementation(g, target; seed=seed)
        # Equality includes every registry entry and every colored flag, not
        # just balance or objective value. Dictionary iteration order is ignored.
        @test (context, name, actual) == (context, name, expected)
        @test count(actual[2]) == target
        @test length(actual[1]) == nv(g)
        counters[:single] += 1
    end
end

function check_multi(g, k, context, counters)
    for seeded in (false, true)
        expected = seeded ? Q.bury_k_seeded(g, k; backend=:legacy) : Q.bury_k(g, k; backend=:legacy)
        for (name, kernel) in COLOR_KERNELS
            actual = BuryMulti.bury_k_fast(g, k; seeded=seeded, kernel=kernel)
            @test (context, seeded, name, actual) == (context, seeded, name, expected)
            @test [count(==(i), values(actual)) for i in 1:k] == fill(nv(g) ÷ k, k)
            counters[:multi] += 1
        end
    end
end

function exhaustive_checks(counters)
    @testset "Every labeled simple graph on 1 to 5 vertices" begin
        for n in 1:5
            @testset "$n vertices" begin
                possible_edges = [(u, v) for u in 1:n for v in (u + 1):n]
                for mask in 0:((1 << length(possible_edges)) - 1)
                    g = SimpleGraph(n)
                    for (i, (u, v)) in enumerate(possible_edges)
                        if !iszero(mask & (1 << (i - 1)))
                            add_edge!(g, u, v)
                        end
                    end
                    original = copy(g)
                    for target in 0:n
                        context = (n=n, graph_mask=mask, target=target, seed=nothing)
                        check_single(g, target, nothing, context, counters)
                        for seed in 1:n
                            context = (n=n, graph_mask=mask, target=target, seed=seed)
                            check_single(g, target, seed, context, counters)
                        end
                    end
                    for k in 1:n
                        if n % k == 0
                            check_multi(g, k, (n=n, graph_mask=mask, k=k), counters)
                        end
                    end
                    @test g == original
                    counters[:graphs] += 1
                end
            end
        end
    end
end

function disconnected_fixture()
    g = SimpleGraph(96)
    for v in 1:23
        add_edge!(g, v, v + 1)
    end
    for v in 34:48
        add_edge!(g, 33, v)
    end
    for v in 65:95
        add_edge!(g, v, v + 1)
    end
    # Several isolated vertices remain between these components.
    return g
end

function representative_checks(counters)
    rng = Xoshiro(0x42555259)
    fixtures = [
        ("path128", path_graph(128)),
        ("cycle120", cycle_graph(120)),
        ("star128", star_graph(128)),
        ("edgeless64", SimpleGraph(64)),
        ("grid6x6", grid([6, 6])),
        ("grid10x12", grid([10, 12])),
        ("disconnected96", disconnected_fixture()),
        ("regular3_128", random_regular_graph(128, 3; rng=rng)),
        ("regular6_120", random_regular_graph(120, 6; rng=rng)),
        ("sparse_random128", erdos_renyi(128, 0.03; rng=rng)),
        ("dense_random96", erdos_renyi(96, 0.5; rng=rng)),
        ("complete48", complete_graph(48)),
        ("grid32x32", grid([32, 32])),
        ("regular3_1024", random_regular_graph(1024, 3; rng=rng)),
    ]
    @testset "Representative graphs and boundary targets" begin
        for (name, g) in fixtures
            @testset "$name" begin
                original = copy(g)
                n = nv(g)
                targets = unique([0, 1, n ÷ 4, n ÷ 2, 3n ÷ 4, n - 1, n])
                seeds = unique([1, max(1, n ÷ 2), n])
                for target in targets
                    check_single(g, target, nothing, (graph=name, target=target, seed=nothing), counters)
                    for seed in seeds
                        check_single(g, target, seed, (graph=name, target=target, seed=seed), counters)
                    end
                end
                for k in unique([1, 2, 3, 4, 8, 16, n])
                    if k <= n && n % k == 0
                        check_multi(g, k, (graph=name, k=k), counters)
                    end
                end
                @test g == original
                counters[:graphs] += 1
            end
        end
    end
end

function search_checks(counters)
    disconnected = SimpleGraph(12)
    for (u, v) in ((1, 2), (2, 3), (5, 6), (6, 7), (7, 8), (8, 5), (10, 11))
        add_edge!(disconnected, u, v)
    end
    fixtures = (("path12", path_graph(12)),
                ("grid4x4", grid([4, 4])),
                ("disconnected12", disconnected))
    # Both original and experimental search implementations intentionally retain
    # the same improvement messages. Discard both outside the comparisons.
    quiet(f) = redirect_stdout(f, devnull)
    @testset "Search wrappers and RNG consumption" begin
        for (name, g) in fixtures
            original = copy(g)
            target = nv(g) ÷ 2
            @testset "$name" begin
                expected_seed = quiet(() -> Q.best_seed(g, target; backend=:legacy))
                expected_range = quiet(() -> Q.bury_not_balanced_bipartition(g, cld(nv(g), 3), target; backend=:legacy))
                for (_, kernel) in COLOR_KERNELS
                    @test quiet(() -> BurySearch.best_seed_fast(g, target; kernel=kernel)) == expected_seed
                    @test quiet(() -> BurySearch.bury_range_fast(g, cld(nv(g), 3), target; kernel=kernel)) == expected_range
                    counters[:search] += 2
                end
                for rngseed in (71, 991)
                    Random.seed!(rngseed)
                    sampled_seed = quiet(() -> Q.best_sqrt_seed(g, target; backend=:legacy))
                    after_sample = rand(UInt64)
                    Random.seed!(rngseed)
                    automatic = quiet(() -> Q.bury_seeded(g, target; backend=:legacy))
                    after_auto = rand(UInt64)
                    for (_, kernel) in COLOR_KERNELS
                        Random.seed!(rngseed)
                        @test quiet(() -> BurySearch.best_seed_fast(g, target; sampled=true, kernel=kernel)) == sampled_seed
                        @test rand(UInt64) == after_sample
                        Random.seed!(rngseed)
                        @test quiet(() -> BurySearch.bury_auto_fast(g, target; kernel=kernel)) == automatic
                        @test rand(UInt64) == after_auto
                        counters[:search] += 2
                    end
                end
                @test g == original
            end
        end
    end
end

function main()
    counters = Dict(:single => 0, :multi => 0, :graphs => 0, :search => 0)
    @testset "BURY experimental compatibility" begin
        exhaustive_checks(counters)
        representative_checks(counters)
        search_checks(counters)
    end
    println("VALIDATION_COUNTS graphs=$(counters[:graphs]) single_comparisons=$(counters[:single]) multipartition_comparisons=$(counters[:multi]) search_graphs=3 search_comparisons=$(counters[:search])")
    println("Exact registry and color equivalence against backend=:legacy passed.")
end

# Invalid oversized targets are deliberately excluded: the existing baseline
# can hang, so calling it would not be a useful equivalence test.
main()
