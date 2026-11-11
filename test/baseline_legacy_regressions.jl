@testset "Baseline and legacy object regressions" begin
    @testset "Balanced baseline validation and random draw policy" begin
        g = path_graph(6)
        for partitioner in (Q.k_partition_random, Q.k_partition_saran_vazirani)
            for k in (-2, 0, 4, 7, 2.0, 2.5)
                Random.seed!(572)
                expected_next = rand(UInt)
                Random.seed!(572)
                @test_throws ArgumentError partitioner(g, k)
                @test rand(UInt) == expected_next
            end
            @test_throws ArgumentError partitioner(SimpleGraph(0), 1)
            @test_throws ArgumentError partitioner(path_graph(5), 2)
        end

        # Supported random calls still shuffle once and fill parts in that order.
        for k in (1, 2, 3, 6)
            Random.seed!(831)
            nodes = shuffle(collect(vertices(g)))
            expected = Dict(v => (i - 1) ÷ (nv(g) ÷ k) + 1 for (i, v) in enumerate(nodes))
            expected_next = rand(UInt)
            Random.seed!(831)
            @test Q.k_partition_random(g, k) == expected
            @test rand(UInt) == expected_next
        end
    end

    @testset "Pack original connected components" begin
        disconnected = SimpleGraph(10)
        for (u, v) in ((1, 2), (2, 3), (4, 5), (5, 6), (7, 8), (9, 10))
            add_edge!(disconnected, u, v)
        end
        for (g, k) in ((SimpleGraph(4), 2), (disconnected, 2), (path_graph(4), 1))
            original = copy(g)
            reg, trace = Q.k_partition_saran_vazirani(g, k)
            p = Q.GraphPartition(reg; nparts=k)
            @test Q.partition_sizes(p) == fill(nv(g) ÷ k, k)
            @test Q.edge_cut_cost(g, p) == 0
            @test length(trace) == 1
            @test trace[1] == (connected_components(g), true)
            @test g == original
            @test Q.partition_sizes(Q.partition_graph(g, k; algorithm=:saran_vazirani)) == fill(nv(g) ÷ k, k)
        end

        g = path_graph(4)
        reg, trace = Q.k_partition_saran_vazirani(g, 2)
        @test trace[1] == (connected_components(g), false)
        @test last(trace)[2]
        @test Q.partition_sizes(Q.GraphPartition(reg; nparts=2)) == [2, 2]
        @test g == path_graph(4)
    end

    @testset "Legacy constructor owns and validates membership" begin
        supplied = Dict(1 => 1, 2 => 3)
        land = Q.quantumLand(supplied, 3, 2)
        supplied[1] = 3
        delete!(supplied, 2)
        @test land.registry == Dict(1 => 1, 2 => 3)
        @test [h.numVillagers for h in land.hamlets] == [1, 0, 1]
        @test [[(v.id, v.address) for v in h.villagers] for h in land.hamlets] == [[(1, 1)], [], [(2, 3)]]
        @test land.numHamlets == 3
        @test land.numVillagers == 2

        for p in (Q.GraphPartition([1, 3]; nparts=4), Q.GraphPartition(Int[]; nparts=3))
            converted = Q.quantumLand(p)
            @test Q.membership(Q.GraphPartition(converted)) == Q.membership(p)
            @test Q.partition_sizes(Q.GraphPartition(converted)) == Q.partition_sizes(p)
            @test [h.numVillagers for h in converted.hamlets] == Q.partition_sizes(p)
        end
        p = Q.GraphPartition([1, 2])
        converted = Q.quantumLand(p)
        converted.registry[1] = 2
        @test Q.membership(p) == [1, 2]

        for integer in (Int8, Int32, UInt, BigInt)
            land = Q.quantumLand(Dict(integer(1) => integer(2)), integer(3), integer(1))
            @test land.registry == Dict(1 => 2)
            @test [h.numVillagers for h in land.hamlets] == [0, 1, 0]
        end
        for (reg, parts, n) in (
            (Dict(), 0, 0), (Dict(), -1, 0), (Dict(), 1.0, 0),
            (Dict(), 1, -1), (Dict(), 1, 0.0),
            (Dict(), big(typemax(Int)) + 1, 0), (Dict(), 1, big(typemax(Int)) + 1),
            (Dict(1 => 1), 1, 2), (Dict(1 => 1, 2 => 1), 1, 1),
            (Dict(0 => 1, 2 => 1), 1, 2), (Dict(1 => 1, 3 => 1), 1, 2),
            (Dict(1.0 => 1), 1, 1), (Dict("1" => 1), 1, 1),
            (Dict(1 => 0), 2, 1), (Dict(1 => 3), 2, 1),
            (Dict(1 => 1.0), 1, 1), (Dict(1 => "1"), 1, 1),
        )
            @test_throws ArgumentError Q.quantumLand(reg, parts, n)
        end
    end

    @testset "Legacy display honors IO throughout nested objects" begin
        land = Q.quantumLand(Dict(1 => 1, 2 => 2), 2, 2)
        villager_text = "Qubit 1, who lives in hamlet 1"
        hamlet_text = "Village 1 has 1 villagers, 0 farms, and 1 mayors.\n" * villager_text * "\n"
        land_text = hamlet_text * "\nVillage 2 has 1 villagers, 0 farms, and 1 mayors.\nQubit 2, who lives in hamlet 2\n\n"
        mktemp() do _, stdout_file
            for (object, expected) in ((land.hamlets[1].villagers[1], villager_text),
                                       (land.hamlets[1], hamlet_text), (land, land_text))
                captured = redirect_stdout(stdout_file) do
                    sprint(show, object)
                end
                @test captured == expected
            end
            flush(stdout_file)
            seekstart(stdout_file)
            @test isempty(read(stdout_file, String))
        end
    end
end
