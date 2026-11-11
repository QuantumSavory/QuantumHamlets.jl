@testset "Existing partition API" begin
    # Recorded against the original single-file implementation before moving it.
    fixtures = (
        (path_graph(4), [1, 1, 2, 2]),
        (cycle_graph(6), [1, 1, 2, 2, 2, 1]),
        (grid([2, 3]), [1, 1, 1, 2, 2, 2]),
    )
    for (g, expected) in fixtures
        original = copy(g)
        for partitioner in (Q.bury_k, Q.bury_k_seeded)
            registry = partitioner(g, 2)
            @test [registry[v] for v in vertices(g)] == expected
            @test g == original
        end
        _, colored = Q.bury_seeded(g, nv(g) ÷ 2, 1)
        @test colored == (expected .== 1)
    end

    g = path_graph(4)
    for target in (0, 2, 4)
        for partitioner in (Q.bury, (g, r) -> Q.bury_seeded(g, r, 1))
            registry, colored = partitioner(g, target)
            @test count(colored) == target
            @test Set(keys(registry)) == Set(vertices(g))
        end
    end
    for k in (1, 2, 4), partitioner in (Q.bury_k, Q.bury_k_seeded, Q.k_partition_random)
        registry = partitioner(g, k)
        @test Set(keys(registry)) == Set(vertices(g))
        @test [count(==(i), values(registry)) for i in 1:k] == fill(4 ÷ k, k)
        @test g == path_graph(4)
    end
    @test Q.inverse_map([4, 2], 5) == [0, 2, 0, 1, 0]
end

@testset "Existing cost API" begin
    for (g, expected) in ((path_graph(4), (1, 1, 1)),
                          (complete_bipartite_graph(2, 2), (4, 2, 1)))
        original = copy(g)
        land = Q.quantumLand(Dict(1 => 1, 2 => 1, 3 => 2, 4 => 2), 2, 4)
        edge_cost, edge_colors = Q.naive_cost_of_partition(land, g)
        match_cost, match_colors = Q.matching_cost_of_partition(land, g)
        rank_cost, metadata = Q.cut_rank_cost(land, g)
        @test (edge_cost, match_cost, rank_cost) == expected
        @test count(==(:red), edge_colors) == edge_cost
        @test count(==(:red), match_colors) == match_cost
        @test length(match_colors) == ne(g)
        @test metadata === nothing
        @test g == original
        cross_graph = Q.remove_local_edges(land, g)
        @test ne(cross_graph) == edge_cost
        @test g == original
    end

    g = grid([6, 6])
    registry = Q.bury_k(g, 3)
    land = Q.quantumLand(registry, 3, nv(g))
    @test [count(==(i), values(registry)) for i in 1:3] == [12, 12, 12]
    @test first(Q.naive_cost_of_partition(land, g)) == 15
    @test first(Q.matching_cost_of_partition(land, g)) == 9
    @test Q.cut_rank_cost(land, g) == (9, nothing)
end
