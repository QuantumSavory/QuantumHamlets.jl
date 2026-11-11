# Desired behavior for pre-existing bugs, deliberately deferred during reorganization.
# @test_broken reports these explicitly and fails if a fix makes one unexpectedly pass.
# See docs/revision_proposal.md for explanations and the later repair plan.
@testset "Known bugs deferred" begin
    @testset "One-hamlet matching should cost zero" begin
        land = Q.quantumLand(Dict(1 => 1, 2 => 1, 3 => 1), 1, 3)
        @test_broken first(Q.matching_cost_of_partition(land, path_graph(3))) == 0
    end

    @testset "A valid candidate may cost more than the number of vertices" begin
        registry, _ = Q.bury_not_balanced_bipartition(
            complete_graph(8), 4, 4; cost_f=Q.naive_cost_of_partition)
        @test_broken registry !== nothing && length(registry) == 8
    end

    @testset "Edgeless graphs still have balanced partitions" begin
        registry, _ = Q.k_partition_saran_vazirani(SimpleGraph(4), 2)
        @test_broken length(registry) == 4
    end

    @testset "Display should write to its IO argument" begin
        captured = redirect_stdout(devnull) do
            sprint(show, Q.villager(1, 1))
        end
        @test_broken occursin("Qubit 1", captured)
    end

    @testset "Constructor input should not alias land membership" begin
        registry = Dict(1 => 1, 2 => 2)
        land = Q.quantumLand(registry, 2, 2)
        registry[1] = 2
        @test_broken land.registry[1] == 1
    end
end

# Do not execute bury(path_graph(4), 5; backend=:legacy): the original
# oversized-target bug hangs. Optimized-backend guards are covered separately.
