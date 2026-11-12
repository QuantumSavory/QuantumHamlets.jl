import Metis
import PythonCall

@testset "Julia annealing baselines" begin
    same_cut(a, b) = membership(a) == membership(b) ||
                     membership(a) == 3 .- membership(b)
    schedule = [1.0, 0.1]

    @testset "Random starts, ownership, and supplied capacities" begin
        for g in (path_graph(2), path_graph(9), grid([3, 3]), SimpleGraph(7)), seed in 0:2
            original = copy(g)
            p = min_cutrank(g; temperatures=schedule, seed)
            @test p isa GraphPartition
            @test sort(partition_sizes(p)) == [fld(nv(g), 2), cld(nv(g), 2)]
            @test 0 <= pairwise_cut_rank(g, p) <= fld(nv(g), 2)
            @test membership(min_cutrank(g; temperatures=schedule, seed)) == membership(p)
            @test g == original
        end

        g = path_graph(9)
        labels = [1, 1, 2, 2, 2, 2, 2, 2, 2]
        initial = GraphPartition(labels; nparts=2)
        for supplied in (initial, labels, labels .== 1, registry(initial))
            before = deepcopy(supplied)
            p = min_cutrank(g; initial=supplied, temperatures=schedule, seed=37)
            @test sort(partition_sizes(p)) == [2, 7]
            @test pairwise_cut_rank(g, p) <= pairwise_cut_rank(g, initial)
            @test same_cut(min_cutrank(g; initial=supplied, temperatures=Float64[]), initial)
            @test supplied isa GraphPartition ? membership(supplied) == membership(before) : supplied == before
        end

        # Complementing labels preserves the cut and canonical row orientation.
        p = min_cutrank(g; initial, temperatures=schedule, seed=37)
        complemented = GraphPartition(3 .- labels; nparts=2)
        @test membership(min_cutrank(g; initial=complemented, temperatures=schedule, seed=37)) == membership(p)
    end

    @testset "Incumbents and RNG isolation" begin
        g = path_graph(9)
        initial = GraphPartition(vcat(ones(Int, 4), fill(2, 5)))
        for seed in 0:5
            best = min_cutrank(g; initial, temperatures=[1.0], seed)
            final = min_cutrank(g; initial, temperatures=[1.0], seed, keep_best=false)
            @test pairwise_cut_rank(g, best) <= pairwise_cut_rank(g, initial)
            @test pairwise_cut_rank(g, best) <= pairwise_cut_rank(g, final)
            @test partition_sizes(final) == partition_sizes(best)
            if seed == 0
                @test pairwise_cut_rank(g, best) == 1
                @test pairwise_cut_rank(g, final) == 2
            end
        end

        pyrandom = PythonCall.pyimport("random")
        saved = pyrandom.getstate()
        try
            pyrandom.seed(1829)
            expected = PythonCall.pyconvert(Float64, pyrandom.random())
            pyrandom.seed(1829)
            min_cutrank(g; initial, temperatures=schedule, seed=81)
            @test PythonCall.pyconvert(Float64, pyrandom.random()) == expected
        finally
            pyrandom.setstate(saved)
        end
        for method in (min_cutrank, bury_sa, bury_seeded_sa)
            Random.seed!(923)
            expected = rand(UInt)
            Random.seed!(923)
            method(g; temperatures=schedule, seed=81)
            @test rand(UInt) == expected
        end
    end

    @testset "BURY initializers compose with the same annealer" begin
        g = grid([3, 3])
        for backend in (:heap, :scan, :legacy)
            ordinary = bipartition_graph(g; objective=pairwise_cut_rank, backend)
            sqrt_seeded = Q._sqrt_bury_bipartition(g; seed=23, backend)
            for (method, initial) in ((bury_sa, ordinary), (bury_seeded_sa, sqrt_seeded))
                @test same_cut(method(g; temperatures=Float64[], seed=23, backend), initial)
                expected = min_cutrank(g; initial, temperatures=schedule, seed=23)
                actual = method(g; temperatures=schedule, seed=23, backend)
                @test membership(actual) == membership(expected)
                @test pairwise_cut_rank(g, actual) <= pairwise_cut_rank(g, initial)
            end
        end

        # Reconstruct the existing sqrt search using an explicit independent RNG.
        rng = Random.Xoshiro(23)
        reg, _ = Q.bury_not_balanced_bipartition(g, fld(nv(g), 2), cld(nv(g), 2);
            bury_v=(graph, target) -> Q.bury_seeded(graph, target; rng, log_improvements=false),
            cost_f=Q.cut_rank_cost)
        @test same_cut(Q._sqrt_bury_bipartition(g; seed=23), GraphPartition(reg; nparts=2))
    end

    @testset "METIS and explicit capacity repair" begin
        original_options = copy(Metis.options)
        g = path_graph(11)
        original = copy(g)
        for k in (1, 2, 3, 7, 11)
            p = metis_partition(g, k; balance=true, seed=5)
            sizes = partition_sizes(p)
            @test length(p) == nv(g)
            @test length(sizes) == k
            @test maximum(sizes) - minimum(sizes) <= 1
            @test minimum(sizes) >= 1
            @test membership(metis_partition(g, k; balance=true, seed=5)) == membership(p)
            @test Metis.options == original_options
        end
        @test g == original
        @test metis_partition(g; seed=5) isa GraphPartition
        @test Metis.options == original_options
        try
            Metis.options[Int(Metis.LibMetis.METIS_OPTION_SEED) + 1] = 913
            custom_options = copy(Metis.options)
            metis_partition(g; seed=5, balance=true)
            @test Metis.options == custom_options
        finally
            copyto!(Metis.options, original_options)
        end

        # A deliberately empty native part is repaired using cut cost and ID ties.
        @test Q._balance_metis(path_graph(6), ones(Int, 6), 2) == [2, 2, 2, 1, 1, 1]

        initial = metis_partition(g; balance=true, seed=5)
        @test same_cut(metis_sa(g; seed=5, temperatures=Float64[]), initial)
        expected = min_cutrank(g; initial, seed=5, temperatures=schedule)
        actual = metis_sa(g; seed=5, temperatures=schedule)
        @test membership(actual) == membership(expected)
        @test sort(partition_sizes(actual)) == [5, 6]
        @test pairwise_cut_rank(g, actual) <= pairwise_cut_rank(g, initial)
        @test Metis.options == original_options

        for k in (0, -1, 12, 2.0), balance in (false, true)
            @test_throws ArgumentError metis_partition(g, k; balance)
            @test Metis.options == original_options
        end
        @test_throws ArgumentError metis_partition(SimpleGraph(0))
        @test_throws ArgumentError metis_partition(g; seed=-1)
        @test_throws ArgumentError metis_partition(g; seed=big(typemax(Int32)) + 1)
        @test Metis.options == original_options
    end

    @testset "Invalid annealing inputs fail before solving" begin
        g = path_graph(5)
        for n in (0, 1)
            @test_throws ArgumentError min_cutrank(SimpleGraph(n))
        end
        for initial in (ones(Int, 5), [1, 1, 2], [1, 1, 1, 2, 3],
                        GraphPartition([1, 1, 1, 2, 2]; nparts=3),
                        Dict(1 => 1, 3 => 2), [1.0, 1.0, 1.0, 2.0, 2.0])
            @test_throws ArgumentError min_cutrank(g; initial)
        end
        for temperatures in ([0.0], [-0.1], [NaN], [Inf])
            @test_throws ArgumentError min_cutrank(g; temperatures)
        end
        for seed in (-1, 0.5)
            @test_throws ArgumentError min_cutrank(g; seed)
        end
        @test_throws ArgumentError bury_sa(g; backend=:unknown)
        @test_throws ArgumentError bury_seeded_sa(g; backend=:unknown)
        looped = copy(g)
        add_edge!(looped, 1, 1)
        @test_throws ArgumentError min_cutrank(looped)
    end
end
