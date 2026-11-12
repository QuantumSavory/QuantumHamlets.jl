# Deliberately enumerate independent BURY runs: target-dependent greedy choices
# are not prefixes of one run, and equal scores select the final target.
function range_reference(g, targets, partitioner, objective)
    best = nothing
    best_score = nothing
    for target in targets
        candidate = partitioner(g, target)
        score = objective(g, candidate)
        if best === nothing || score <= best_score
            best = deepcopy(candidate)
            best_score = score
        end
    end
    return best
end

function range_capture_rng(f, seed)
    Random.seed!(seed)
    mktemp() do _, io
        result = redirect_stdout(f, io)
        flush(io)
        seekstart(io)
        return result, read(io, String), rand(UInt64)
    end
end

function range_invalid_without_rng(f)
    Random.seed!(0x72616e6765)
    next_draw = rand(UInt64)
    Random.seed!(0x72616e6765)
    @test_throws ArgumentError f()
    @test rand(UInt64) == next_draw
end

@testset "Bipartition target ranges" begin
    backends = (:heap, :scan, :legacy)
    objectives = ((matching_cost, Q.matching_cost_of_partition),
                  (edge_cut_cost, Q.naive_cost_of_partition),
                  (pairwise_cut_rank, Q.cut_rank_cost))

    odd = SimpleGraph(5)
    add_edge!(odd, 1, 4)
    add_edge!(odd, 2, 3)
    nonnested = SimpleGraph(4)
    for (u, v) in ((1, 3), (1, 4), (2, 3), (2, 4), (3, 4))
        add_edge!(nonnested, u, v)
    end

    @testset "Independent enumeration, built-in costs, and return types" begin
        for g in (SimpleGraph(0), SimpleGraph(1), path_graph(2), path_graph(5),
                  path_graph(6), grid([3, 3]), complete_graph(8), odd, nonnested)
            original = copy(g)
            n = nv(g)
            ranges = unique([0:n, fld(n, 2):fld(n, 2), fld(n, 2):cld(n, 2),
                             max(0, fld(n, 2) - 2):min(n, cld(n, 2) + 2)])
            for backend in backends, (scalar, legacy) in objectives, targets in ranges
                partitioner = (graph, target) -> Q.bury(graph, target; backend)
                score = (graph, result) -> scalar(graph, GraphPartition(result[1]; nparts=2))
                expected = range_reference(g, targets, partitioner, score)
                actual = Q.bury_not_balanced_bipartition(g, first(targets), last(targets);
                    backend, cost_f=legacy)
                @test actual == expected
                @test typeof(actual) === Tuple{Dict{Any,Any},Vector{Bool}}
                @test [actual[1][v] == 2 for v in vertices(g)] == actual[2]

                if n >= 2 && first(targets) >= 1 && last(targets) < n
                    p = bipartition_graph(g; targets, backend, objective=scalar)
                    @test p isa GraphPartition
                    @test p.nparts == 2
                    @test registry(p) == expected[1]
                    @test all(>(0), partition_sizes(p))
                end
            end
            if n >= 2
                for backend in backends, (scalar, _) in objectives
                    targets = fld(n, 2):cld(n, 2)
                    expected = range_reference(g, targets,
                        (graph, target) -> Q.bury(graph, target; backend),
                        (graph, result) -> scalar(graph, GraphPartition(result[1]; nparts=2)))
                    @test registry(bipartition_graph(g; backend, objective=scalar)) == expected[1]
                end
                @test registry(bipartition_graph(g)) ==
                      registry(bipartition_graph(g; backend=:heap, objective=matching_cost))
            end
            @test g == original
        end
    end

    @testset "Complementary targets and target-dependent choices" begin
        for backend in backends
            target2 = Q.bury(odd, 2; backend)
            target3 = Q.bury(odd, 3; backend)
            @test findall(target2[2]) == [1, 5]
            @test findall(target3[2]) == [1, 4, 5]
            @test matching_cost(odd, GraphPartition(target2[1]; nparts=2)) == 1
            @test matching_cost(odd, GraphPartition(target3[1]; nparts=2)) == 0
            p = bipartition_graph(odd; backend)
            @test registry(p) == target3[1]
            @test partition_sizes(p) == [2, 3]
            @test matching_cost(odd, p) == 0

            @test findall(last(Q.bury(nonnested, 2; backend))) == [1, 2]
            @test findall(last(Q.bury(nonnested, 3; backend))) == [1, 3, 4]
            seen = Vector{Vector{Int}}()
            objective = function (_, candidate)
                push!(seen, findall(==(2), membership(candidate)))
                return 0
            end
            p = bipartition_graph(nonnested; targets=2:3, backend, objective)
            @test seen == [[1, 2], [1, 3, 4]]
            @test registry(p) == first(Q.bury(nonnested, 3; backend))

            # Empty sides remain legal only on the qualified low-level API.
            @test Q.bury_not_balanced_bipartition(odd, 0, nv(odd); backend) ==
                  Q.bury(odd, nv(odd); backend)
            @test Q.bury_not_balanced_bipartition(SimpleGraph(0), 0, 0; backend) ==
                  (Dict{Any,Any}(), Bool[])
            @test partition_sizes(bipartition_graph(path_graph(2); backend)) == [1, 1]
            @test_throws ArgumentError partition_graph(odd, 2; backend)
        end
    end

    @testset "Graph-first scalar callbacks and target validation" begin
        for backend in backends
            for n in (2, 5, 6)
                g = path_graph(n)
                original = copy(g)
                seen = Int[]
                objective = function (graph, p)
                    @test graph === g
                    @test p isa GraphPartition
                    @test p.nparts == 2
                    push!(seen, partition_sizes(p)[2])
                    return 0
                end
                p = bipartition_graph(g; backend, objective)
                @test seen == collect(fld(n, 2):cld(n, 2))
                @test registry(p) == first(Q.bury(g, cld(n, 2); backend))
                @test g == original
            end

            g = path_graph(7)
            for scores in ([101, 100, 100], [-3, -2, -1], [0, 2, 0], [2, 0, 1])
                seen = Int[]
                objective = function (graph, p)
                    @test graph === g
                    target = partition_sizes(p)[2]
                    push!(seen, target)
                    return scores[target - 1]
                end
                p = bipartition_graph(g; targets=2:4, backend, objective)
                expected_target = findlast(==(minimum(scores)), scores) + 1
                @test registry(p) == first(Q.bury(g, expected_target; backend))
                @test seen == [2, 3, 4]
            end

            # Integer representations and explicit unit-step ranges are valid.
            for targets in (Int8(2):Int8(4), big(2):big(4), 2:1:4, 3:3)
                @test registry(bipartition_graph(g; targets, backend)) ==
                      first(Q.bury_not_balanced_bipartition(g, first(targets), last(targets); backend))
            end
            calls = Ref(0)
            objective = (_, _) -> (calls[] += 1; 0)
            for targets in (0:2, 2:7, -1:2, 4:3, 4:-1:2, 2:2:4,
                            2.0:4.0, 2:0.5:4, [2, 3, 4], (2, 4), 2,
                            big(2):(big(typemax(Int)) + 1))
                range_invalid_without_rng(() -> bipartition_graph(g; targets, backend, objective))
            end
            for graph in (SimpleGraph(0), SimpleGraph(1))
                range_invalid_without_rng(() -> bipartition_graph(graph; backend, objective))
                range_invalid_without_rng(() -> bipartition_graph(graph; targets=1:1, backend, objective))
            end
            range_invalid_without_rng(() -> bipartition_graph(g; backend=:unknown, objective))
            @test calls[] == 0
        end
    end

    @testset "Legacy callbacks retain validation, signatures, order, and errors" begin
        g = path_graph(5)
        for backend in backends
            for scores in ([101, 100, 100], [-3, -2, -1], [0, 2, 0], [2, 0, 1])
                events = Tuple{Symbol,Int}[]
                # This callback accepts exactly two positional arguments and
                # deliberately controls its own backend.
                partitioner = function (graph, target)
                    @test graph === g
                    push!(events, (:partition, target))
                    return Q.bury(graph, target; backend=:scan)
                end
                objective = function (land, graph)
                    @test land isa Q.quantumLand
                    @test graph === g
                    @test land.numHamlets == 2
                    @test land.numVillagers == nv(g)
                    target = count(==(2), values(land.registry))
                    @test land.hamlets[2].numVillagers == target
                    push!(events, (:score, target))
                    return scores[target], :unused
                end
                result = Q.bury_not_balanced_bipartition(g, 1, 3;
                    backend, bury_v=partitioner, cost_f=objective)
                target = findlast(==(minimum(scores)), scores)
                @test result == Q.bury(g, target; backend=:scan)
                @test events == [(:partition, 1), (:score, 1), (:partition, 2),
                                 (:score, 2), (:partition, 3), (:score, 3)]
            end

            for (scalar, legacy) in objectives
                seen = Int[]
                # Wrapping a known objective is still a callback invocation.
                wrapped = function (land, graph)
                    push!(seen, count(==(2), values(land.registry)))
                    return legacy(land, graph)
                end
                result = Q.bury_not_balanced_bipartition(g, 1, 4; backend, cost_f=wrapped)
                @test result == Q.bury_not_balanced_bipartition(g, 1, 4; backend, cost_f=legacy)
                @test seen == [1, 2, 3, 4]
                empty!(seen)
                wrapped_scalar = function (graph, p)
                    push!(seen, partition_sizes(p)[2])
                    return scalar(graph, p)
                end
                @test registry(bipartition_graph(g; targets=1:4, backend, objective=wrapped_scalar)) == result[1]
                @test seen == [1, 2, 3, 4]
            end

            for invalid in (Dict(1 => 1), Dict(v => 3 for v in vertices(g)),
                            Dict(v => 1.0 for v in vertices(g)))
                calls = Ref(0)
                partitioner = (_, _) -> (calls[] += 1; (invalid, falses(nv(g))))
                for (_, legacy) in objectives
                    calls[] = 0
                    @test_throws ArgumentError Q.bury_not_balanced_bipartition(g, 1, 3;
                        backend, bury_v=partitioner, cost_f=legacy)
                    @test calls[] == 1
                end
                scores = Ref(0)
                @test_throws ArgumentError Q.bury_not_balanced_bipartition(g, 1, 3;
                    backend, bury_v=partitioner, cost_f=(_, _) -> (scores[] += 1; (0, nothing)))
                @test scores[] == 0
            end

            # A MethodError thrown inside a callback must propagate unchanged;
            # it must never trigger a retry with another calling convention.
            for location in (:partition, :legacy_cost, :scalar_cost)
                calls = Ref(0)
                failure = MethodError(identity, (:callback_failure, :sentinel))
                callback = (_, _) -> (calls[] += 1; throw(failure))
                actual = try
                    if location === :partition
                        Q.bury_not_balanced_bipartition(g, 1, 3; backend, bury_v=callback)
                    elseif location === :legacy_cost
                        Q.bury_not_balanced_bipartition(g, 1, 3; backend, cost_f=callback)
                    else
                        bipartition_graph(g; targets=1:3, backend, objective=callback)
                    end
                catch err
                    err
                end
                @test actual === failure
                @test calls[] == 1
            end
        end
    end

    @testset "Retained callback results own their buffers" begin
        g = path_graph(4)
        for backend in backends
            shared_reg = Dict(v => 1 for v in vertices(g))
            shared_colors = fill(false, nv(g))
            calls = Int[]
            candidate = function (graph, target)
                push!(calls, target)
                for v in vertices(graph)
                    shared_colors[v] = v <= target
                    shared_reg[v] = shared_colors[v] ? 2 : 1
                end
                return shared_reg, shared_colors
            end
            score = (land, _) -> (count(==(2), values(land.registry)), nothing)
            reg, colors = Q.bury_not_balanced_bipartition(g, 1, 3;
                backend, bury_v=candidate, cost_f=score)
            expected_reg = Dict(1 => 2, 2 => 1, 3 => 1, 4 => 1)
            expected_colors = [true, false, false, false]
            @test calls == [1, 2, 3]
            @test reg == expected_reg
            @test colors == expected_colors
            @test [reg[v] == 2 for v in vertices(g)] == colors
            @test reg !== shared_reg
            @test colors !== shared_colors
            @test typeof(reg) === typeof(shared_reg)
            @test typeof(colors) === typeof(shared_colors)
            # Mutations after search completion cannot alter either winner.
            empty!(shared_reg)
            fill!(shared_colors, false)
            @test reg == expected_reg
            @test colors == expected_colors
            @test [reg[v] == 2 for v in vertices(g)] == colors

            # A custom scalar callback may retain the objects it sees. They
            # must not expose storage shared with the returned partition.
            observed = GraphPartition[]
            objective = function (_, p)
                push!(observed, p)
                return partition_sizes(p)[2]
            end
            p = bipartition_graph(g; targets=1:3, backend, objective)
            expected = membership(p)
            for candidate in observed
                fill!(candidate.labels, 1)
            end
            @test membership(p) == expected
            @test partition_sizes(p) == [3, 1]
        end
    end

    @testset "Seeded target composition preserves matching seed search and RNG" begin
        g = path_graph(7)
        original = copy(g)
        for backend in backends, (_, legacy) in objectives
            expected = range_capture_rng(0x73656564) do
                range_reference(g, 2:4,
                    (graph, target) -> Q.bury_seeded(graph, target; backend),
                    (graph, result) -> first(legacy(Q.quantumLand(result[1], 2, nv(graph)), graph)))
            end
            actual = range_capture_rng(0x73656564) do
                Q.bury_not_balanced_bipartition(g, 2, 4;
                    backend, bury_v=Q.bury_seeded, cost_f=legacy)
            end
            @test actual == expected # includes logs and the next RNG draw
            @test typeof(actual[1]) === Tuple{Dict{Any,Any},Vector{Bool}}
        end
        @test g == original
    end
end
