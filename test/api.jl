@testset "GraphPartition membership" begin
    labels = [1, 1, 2, 2]
    p = Q.GraphPartition(labels)
    labels[1] = 2
    @test Q.membership(p) == [1, 1, 2, 2]
    @test length(p) == 4
    @test Q.partition_sizes(p) == [2, 2]

    copied_labels = Q.membership(p)
    copied_registry = Q.registry(p)
    copied_labels[1] = 2
    copied_registry[1] = 2
    @test Q.membership(p) == [1, 1, 2, 2]

    input_registry = Dict(1 => 1, 2 => 2)
    from_registry = Q.GraphPartition(input_registry; nparts=3)
    input_registry[1] = 2
    @test Q.membership(from_registry) == [1, 2]
    @test Q.partition_sizes(from_registry) == [1, 1, 0]
    @test Q.registry(p) == Dict(1 => 1, 2 => 1, 3 => 2, 4 => 2)
    @test Q.membership(Q.GraphPartition(Q.quantumLand(p))) == Q.membership(p)
    @test Q.partition_sizes(Q.GraphPartition(Int[]; nparts=2)) == [0, 0]

    @test_throws ArgumentError Q.GraphPartition([0, 1])
    @test_throws ArgumentError Q.GraphPartition([1, 3]; nparts=2)
    @test_throws ArgumentError Q.GraphPartition([1]; nparts=0)
    @test_throws ArgumentError Q.GraphPartition(Dict(1 => 1, 3 => 2))
    @test_throws ArgumentError Q.GraphPartition(Dict(1 => 1.5))
end

@testset "Simplified partition and cost API" begin
    g = grid([2, 3])
    original = copy(g)
    for (algorithm, legacy) in ((:bury, Q.bury_k), (:bury_seeded, Q.bury_k_seeded))
        p = Q.partition_graph(g, 2; algorithm)
        @test Q.registry(p) == legacy(g, 2)
        @test Q.partition_sizes(p) == [3, 3]
        @test g == original
    end
    Random.seed!(123)
    random_partition = Q.partition_graph(g, 2; algorithm=:random)
    Random.seed!(123)
    @test Q.registry(random_partition) == Q.k_partition_random(g, 2)
    @test Q.partition_sizes(random_partition) == [3, 3]
    @test Q.partition_sizes(Q.partition_graph(g, 1)) == [6]
    @test Q.partition_sizes(Q.partition_graph(g, 6)) == ones(Int, 6)
    for k in (0, -1, 4, 7, 2.0, 2.5)
        @test_throws ArgumentError Q.partition_graph(g, k)
    end
    @test_throws ArgumentError Q.partition_graph(g, 2; algorithm=:unknown)
    @test Q.partition_sizes(Q.partition_graph(SimpleGraph(4), 2; algorithm=:saran_vazirani)) == [2, 2]

    g = complete_bipartite_graph(2, 2)
    p = Q.GraphPartition([1, 1, 2, 2])
    land = Q.quantumLand(p)
    @test_throws ArgumentError Q.edge_cut_cost(path_graph(3), p)
    @test Q.edge_cut_cost(g, p) == first(Q.naive_cost_of_partition(land, g)) == 4
    @test Q.matching_cost(g, p) == first(Q.matching_cost_of_partition(land, g)) == 2
    @test Q.pairwise_cut_rank(g, p) == first(Q.cut_rank_cost(land, g)) == 1
    free, paid = Q.vertex_cover_cost_of_partition(land, g)
    @test Q.vertex_cover_costs(g, p) == (ghz_free=free, ghz_n_minus_one=paid)

    details = Q.matching_details(g, p)
    @test details.cost == 2
    @test length(details.edges) == 2
    @test all(e -> has_edge(g, e), details.edges)
    @test length(Set(v for e in details.edges for v in (src(e), dst(e)))) == 4
    @test g == complete_bipartite_graph(2, 2)
end
