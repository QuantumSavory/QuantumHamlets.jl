using Test, Graphs, QuantumHamlets
const Q = QuantumHamlets

function observe(name, f)
    try
        value = f()
        println(name, "\tRETURN\t", repr(value))
        return value
    catch ex
        println(name, "\tERROR\t", nameof(typeof(ex)), "\t", sprint(showerror, ex))
        return ex
    end
end

@testset "Current status of earlier findings" begin
    g = complete_graph(8)
    for backend in (:heap, :scan, :legacy)
        result = observe("range_score_bound_$backend", () ->
            Q.bury_not_balanced_bipartition(g, 4, 4; backend, cost_f=Q.naive_cost_of_partition))
        @test result == (nothing, nothing)
        candidate, _ = Q.bury(g, 4; backend)
        @test first(Q.naive_cost_of_partition(Q.quantumLand(candidate, 2, 8), g)) == 16
    end

    result = observe("one_part_matching", () -> matching_cost(path_graph(3), GraphPartition([1,1,1])))
    @test result isa MethodError

    reg = Dict(1=>1, 2=>2)
    land = Q.quantumLand(reg, 2, 2)
    reg[1] = 2
    @test land.registry[1] == 2 && land.hamlets[1].villagers[1].address == 1
    println("constructor_alias\tregistry=", land.registry[1], "\tstored_address=", land.hamlets[1].villagers[1].address)
    p = GraphPartition(Dict(1=>1,2=>2))
    copied = registry(p); copied[1] = 2
    @test membership(p) == [1,2]

    printed = mktemp() do _, io
        captured = redirect_stdout(io) do
            sprint(show, Q.villager(1, 1))
        end
        seekstart(io)
        (captured, read(io, String))
    end
    @test printed == ("", "Qubit 1, who lives in hamlet 1")
    println("show_wrong_stream\tcaptured=", repr(printed[1]), "\tstdout=", repr(printed[2]))

    g = path_graph(4)
    for backend in (:heap, :scan, :legacy)
        @test isempty(observe("negative_k_$backend", () -> Q.bury_k(g, -2; backend)))
        @test observe("zero_k_$backend", () -> Q.bury_k(g, 0; backend)) isa DivideError
        @test Q.bury_not_balanced_bipartition(g, 3, 2; backend) == (nothing, nothing)
    end
    for backend in (:heap, :scan)
        @test observe("oversized_target_guard_$backend", () -> Q.bury(g, 5; backend)) isa ArgumentError
        @test observe("oversized_seeded_target_guard_$backend", () -> Q.bury_seeded(g, 5, 1; backend)) isa ArgumentError
    end
    @test_throws ArgumentError partition_graph(g, -2)
    @test_throws ArgumentError partition_graph(g, 0)
    @test observe("random_indivisible", () -> Q.k_partition_random(path_graph(5), 2)) isa InexactError
end

# Deliberately do not execute oversized legacy targets: the infinite loop is
# already demonstrated in the original review and unchanged in src/bury/legacy.jl.
