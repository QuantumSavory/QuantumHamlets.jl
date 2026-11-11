# Fixed post-bug-fix implementations, saved before cost optimizations. The
# shared remove_local_edges helper is intentionally retained unchanged.
module CostSemanticsReference
using Graphs
import QuantumHamlets
using QuantumHamlets: quantumLand, inverse_map, IGraphs, residue_ring, ZZ, rank, matrix
include(joinpath(@__DIR__, "..", "experiments", "cost_performance", "reference_costs.jl"))
end

function cost_with_next_draw(f, land, graph)
    Random.seed!(9143)
    result = f(land, graph)
    return result, rand(UInt)
end

function check_cost_semantics(graph, labels, parts)
    p = GraphPartition(labels; nparts=parts)
    land = Q.quantumLand(p)
    original_graph = copy(graph)
    original_registry = copy(land.registry)

    expected, next_draw = cost_with_next_draw(CostSemanticsReference.naive_cost_of_partition, land, graph)
    actual, actual_next = cost_with_next_draw(Q.naive_cost_of_partition, land, graph)
    @test actual == expected
    @test actual_next == next_draw
    @test edge_cut_cost(graph, p) == first(expected)

    expected, next_draw = cost_with_next_draw(CostSemanticsReference.cut_rank_cost, land, graph)
    actual, actual_next = cost_with_next_draw(Q.cut_rank_cost, land, graph)
    @test actual == expected
    @test actual_next == next_draw
    @test pairwise_cut_rank(graph, p) == first(expected)

    expected, next_draw = cost_with_next_draw(CostSemanticsReference.vertex_cover_cost_of_partition, land, graph)
    actual, actual_next = cost_with_next_draw(Q.vertex_cover_cost_of_partition, land, graph)
    @test actual == expected
    @test actual_next == next_draw
    @test vertex_cover_costs(graph, p) == (ghz_free=expected[1], ghz_n_minus_one=expected[2])

    expected, next_draw = cost_with_next_draw(CostSemanticsReference.matching_cost_of_partition, land, graph)
    actual, actual_next = cost_with_next_draw(Q.matching_cost_of_partition, land, graph)
    @test first(actual) == first(expected)
    @test actual_next == next_draw
    @test matching_cost(graph, p) == first(expected)
    cost, colors = actual
    @test length(colors) == ne(graph)
    @test all(color -> color in (:red, :black), colors)
    selected = [e for (e, color) in zip(edges(graph), colors) if color === :red]
    @test length(selected) == cost
    @test all(e -> labels[src(e)] != labels[dst(e)], selected)
    for i in 1:parts, j in (i + 1):parts
        pair_edges = [e for e in selected if minmax(labels[src(e)], labels[dst(e)]) == (i, j)]
        @test length(Set(v for e in pair_edges for v in (src(e), dst(e)))) == 2length(pair_edges)
    end
    @test graph == original_graph
    @test land.registry == original_registry
    @test membership(p) == labels
end

@testset "Cost semantics after allocation reductions" begin
    @testset "Tiny exhaustive graphs and memberships" begin
        for n in 0:3
            possible = [(u, v) for u in 1:n for v in (u + 1):n]
            for mask in 0:((1 << length(possible)) - 1)
                graph = SimpleGraph(n)
                for (i, (u, v)) in enumerate(possible)
                    !iszero(mask & (1 << (i - 1))) && add_edge!(graph, u, v)
                end
                for parts in 1:3, code in 0:(parts^n - 1)
                    labels = [((code ÷ parts^(v - 1)) % parts) + 1 for v in 1:n]
                    check_cost_semantics(graph, labels, parts)
                end
            end
        end
    end

    @testset "Varied densities, repeated neighbor labels, and empty parts" begin
        fixture_rng = MersenneTwister(7301)
        for trial in 1:80
            n = rand(fixture_rng, 4:14)
            parts = rand(fixture_rng, 1:7)
            probability = rand(fixture_rng, (0.0, 0.1, 0.35, 0.7, 1.0))
            graph = SimpleGraph(n)
            for u in 1:n, v in (u + 1):n
                rand(fixture_rng) < probability && add_edge!(graph, u, v)
            end
            labels = rand(fixture_rng, 1:parts, n)
            check_cost_semantics(graph, labels, parts)
        end

        # Several neighbors in the same part must contribute only once to a
        # selected cover vertex's GHZ m-minus-one term.
        check_cost_semantics(star_graph(9), [1, 2, 2, 2, 3, 3, 4, 4, 4], 6)
        check_cost_semantics(complete_graph(6), [1, 1, 2, 2, 3, 3], 5)
        check_cost_semantics(SimpleGraph{Int32}(path_graph(7)), [1, 2, 3, 1, 2, 3, 1], 4)
    end
end
