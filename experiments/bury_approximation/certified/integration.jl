using Graphs, QuantumHamlets, Test
include("certified_bury.jl")
using .CertifiedBury

module PathHubFixtures
include("../counterexamples/path_hubs.jl")
end

@testset "Component guard repairs zero-cost packing" begin
    # The original greedy order consumes both length-two paths, then must
    # split a length-three path. Packing one of each gives a zero-cost cut.
    g = SimpleGraph(10)
    for (u, v) in ((1,2), (3,4), (5,6), (6,7), (8,9), (9,10))
        add_edge!(g, u, v)
    end
    before = collect(edges(g))
    reg, _ = QuantumHamlets.bury(g, 5)
    @test matching_cost(g, GraphPartition(reg; nparts=2)) == 1
    for backend in (:heap, :scan, :legacy), seed in (nothing, 1, 10)
        p = component_guarded_bury(g, 5; seed=seed, backend=backend)
        @test partition_sizes(p) == [5,5]
        @test matching_cost(g, p) == 0
        @test edge_cut_cost(g, p) == 0
    end
    @test collect(edges(g)) == before
end

@testset "Component guard boundaries and validation before successful prepass" begin
    empty = SimpleGraph(0)
    @test isempty(membership(component_guarded_bury(empty, 0)))
    for backend in (:heap, :scan, :legacy), seed in (nothing, 1)
        g = path_graph(4)
        @test partition_sizes(component_guarded_bury(g, 0; seed, backend)) == [4,0]
        @test partition_sizes(component_guarded_bury(g, 4; seed, backend)) == [0,4]
    end
    @test_throws ArgumentError component_guarded_bury(empty, 0; seed=1)
    @test_throws ArgumentError component_guarded_bury(empty, 0; backend=:missing)
    @test_throws ArgumentError component_guarded_bury(path_graph(4), 0; seed=0)
    @test_throws ArgumentError component_guarded_bury(path_graph(4), 4; seed=5)
    @test_throws ArgumentError component_guarded_bury(path_graph(4), 0; seed=1.0)
    @test_throws ArgumentError component_guarded_bury(path_graph(4), 0.0)
    @test_throws ArgumentError component_guarded_bury(path_graph(4), -1)
    @test_throws ArgumentError component_guarded_bury(path_graph(4), 5)
    @test_throws ArgumentError component_guarded_bury(path_graph(4), 4; backend=:missing)
end

@testset "Fallback is exactly the requested BURY variant" begin
    for g in (path_graph(8), cycle_graph(8), complete_graph(8))
        before = collect(edges(g))
        for backend in (:heap, :scan, :legacy), seed in (nothing, 1, 4, 8)
            p = component_guarded_bury(g, 4; seed, backend)
            _, expected = seed === nothing ?
                QuantumHamlets.bury(g, 4; backend) :
                QuantumHamlets.bury_seeded(g, 4, seed; backend)
            @test membership(p) == Int.(expected) .+ 1
            @test partition_sizes(p) == [4,4]
        end
        @test collect(edges(g)) == before
    end
end

@testset "One-vertex certification repairs bounded-degree counterexample" begin
    println("r\tn\toriginal_matching\trepaired_matching\tlower_bound\texamined")
    for r in (3, 4, 8)
        g, _, _, _, witness = PathHubFixtures.path_hub_obstruction(r; b=4)
        target = nv(g) ÷ 2
        before = collect(edges(g))
        result = certified_bury(g, target; max_cover=1)
        @test result.original_cost >= r+1
        @test matching_cost(g, witness) == 1
        @test result.cost == result.lower_bound == 1
        @test result.optimal
        @test result.factor == 1
        @test partition_sizes(result.partition) == [target,target]
        @test length(result.separator) == 1
        @test matching_cost(g, result.partition) == 1
        labels = membership(result.partition)
        @test all(e -> labels[src(e)] == labels[dst(e)] ||
                       src(e) in result.separator || dst(e) in result.separator,
                  edges(g))
        @test collect(edges(g)) == before
        println(join((r,nv(g),result.original_cost,result.cost,
                      result.lower_bound,result.examined), '\t'))
    end
end
