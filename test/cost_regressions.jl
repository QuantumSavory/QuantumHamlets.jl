# Independent cardinality oracle: enumerate assignments of right-side vertices
# rather than comparing selected edges with another call to the same solver.
function exhaustive_pair_matching(g, left, right)
    function search(i, used)
        i > length(left) && return 0
        best = search(i + 1, used)
        for (j, v) in enumerate(right)
            bit = UInt(1) << (j - 1)
            if iszero(used & bit) && has_edge(g, left[i], v)
                best = max(best, 1 + search(i + 1, used | bit))
            end
        end
        return best
    end
    return search(1, UInt(0))
end

function check_pairwise_matching(g, labels, k)
    p = GraphPartition(labels; nparts=k)
    original = copy(g)
    details = matching_details(g, p)
    expected_total = 0
    for i in 1:k, j in (i + 1):k
        expected = exhaustive_pair_matching(g, findall(==(i), labels), findall(==(j), labels))
        selected = [e for e in details.edges if minmax(labels[src(e)], labels[dst(e)]) == (i, j)]
        expected_total += expected
        @test length(selected) == expected
        @test length(Set(v for e in selected for v in (src(e), dst(e)))) == 2length(selected)
    end
    @test details.cost == matching_cost(g, p) == expected_total
    @test length(details.edges) == expected_total
    @test all(e -> has_edge(g, e) && labels[src(e)] != labels[dst(e)], details.edges)
    @test g == original
    @test membership(p) == labels
end

@testset "Matching cost regressions" begin
    @testset "One part has no matching edges" begin
        for g in (SimpleGraph(0), SimpleGraph(1), path_graph(3))
            p = GraphPartition(ones(Int, nv(g)))
            @test matching_cost(g, p) === 0
            details = matching_details(g, p)
            @test details.cost === 0
            @test isempty(details.edges)
            cost, colors = Q.matching_cost_of_partition(Q.quantumLand(p), g)
            @test cost === 0
            @test colors == fill(:black, ne(g))
        end
        g = path_graph(3)
        p = GraphPartition([1, 1, 1])
        figure, cost = plot_partition(g, p; highlight=:matching)
        @test figure isa Q.Makie.Figure
        @test cost === 0
        figure, cost = Q.visualize_graph_on_land(Q.quantumLand(p), g;
                                                method=Q.matching_cost_of_partition)
        @test figure isa Q.Makie.Figure
        @test cost === 0
    end

    @testset "Independent cardinality and witnesses" begin
        # All labeled simple graphs and one/two/three-part assignments through
        # four vertices exercise disconnected sides, isolated vertices, empty
        # parts, repeated component sizes, and non-unique maximum matchings.
        for n in 0:4
            possible = [(u, v) for u in 1:n for v in (u + 1):n]
            for mask in 0:((1 << length(possible)) - 1)
                g = SimpleGraph(n)
                for (i, (u, v)) in enumerate(possible)
                    !iszero(mask & (1 << (i - 1))) && add_edge!(g, u, v)
                end
                for k in 1:3, code in 0:(k^n - 1)
                    labels = n == 0 ? Int[] : [((code ÷ k^(v - 1)) % k) + 1 for v in 1:n]
                    check_pairwise_matching(g, labels, k)
                end
            end
        end
    end
end
