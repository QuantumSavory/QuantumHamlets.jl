module ResearchScoreRegression
using Test, Graphs, QuantumHamlets, Statistics
include("../resources/comparison_scores.jl")

@testset "Scalar research comparison objective" begin
    # Rank differs from matching here, so a mixed-objective comparison fails.
    g = complete_bipartite_graph(2, 2)
    assignments = ([1, 1, 2, 2], [1, 2, 1, 2])
    scores = [comparison_score(g, labels, 2) for labels in assignments]
    @test scores == [1, 2]
    @test all(score -> score isa Int, scores)
    @test matching_cost(g, GraphPartition(first(assignments))) == 2
    @test mean(scores) == 1.5
    @test std(scores) ≈ sqrt(0.5)
    # The same registry/label forms used by BURY and METIS give the same score.
    for labels in assignments
        @test comparison_score(g, Dict(v => labels[v] for v in vertices(g)), 2) ==
              comparison_score(g, labels, 2)
    end
    f = QuantumHamlets.CairoMakie.Figure()
    ax = QuantumHamlets.CairoMakie.Axis(f[1, 1]; ylabel=COMPARISON_OBJECTIVE_LABEL)
    plot = QuantumHamlets.CairoMakie.scatterlines!(ax, [4, 4], scores)
    @test plot isa QuantumHamlets.Makie.ScatterLines
    @test ax.ylabel[] == "Pairwise cut rank (GF(2))"
end
end
