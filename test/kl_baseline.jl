_kl_test_cut(g, labels) = count(e -> labels[src(e)] != labels[dst(e)], edges(g))

# Deliberately recompute the full cut for every hypothetical swap. This oracle
# does not share the implementation's gain formula or incremental updates.
function _kl_test_pass(g, initial)
    labels = copy(initial)
    left, right = findall(==(1), initial), findall(==(2), initial)
    pairs = Tuple{Int,Int}[]
    costs = Int[_kl_test_cut(g, labels)]
    while !isempty(left) && !isempty(right)
        best = (typemax(Int), typemax(Int), typemax(Int))
        for a in left, b in right
            trial = copy(labels)
            trial[a], trial[b] = trial[b], trial[a]
            candidate = (_kl_test_cut(g, trial), a, b)
            candidate < best && (best = candidate)
        end
        cost, a, b = best
        labels[a], labels[b] = labels[b], labels[a]
        filter!(!=(a), left)
        filter!(!=(b), right)
        push!(pairs, (a, b))
        push!(costs, cost)
    end
    return pairs, costs
end

function _kl_test_reference(g, initial; max_iter=10)
    labels = copy(initial)
    for _ in 1:max_iter
        pairs, costs = _kl_test_pass(g, labels)
        prefix = argmin(costs) - 1 # First minimum; zero means no improvement.
        prefix == 0 && break
        for (a, b) in pairs[1:prefix]
            labels[a], labels[b] = labels[b], labels[a]
        end
    end
    return labels
end

function _kl_test_graph(n, pairs; vertex_type=Int)
    g = SimpleGraph{vertex_type}(n)
    for (u, v) in pairs
        add_edge!(g, u, v)
    end
    return g
end

@testset "Classic Kernighan–Lin baseline" begin
    @test kl_partition === Q.kl_partition
    @test :kl_partition in names(Q)

    @testset "Exhaustive small-graph oracle" begin
        pairs = [(u, v) for u in 1:4 for v in (u + 1):4]
        starts = [[isodd(mask >> (v - 1)) ? 1 : 2 for v in 1:4]
                  for mask in 0:15 if count_ones(mask) == 2]
        append!(starts, [[1, 2, 2, 2], [1, 1, 1, 2]])
        for mask in 0:63
            g = _kl_test_graph(4, [pairs[i] for i in eachindex(pairs) if isodd(mask >> (i - 1))])
            for initial in starts, limit in (1, 5)
                expected = _kl_test_reference(g, initial; max_iter=limit)
                p = Q.kl_partition(g; initial=initial, max_iter=limit)
                @test Q.membership(p) == expected
                @test Q.partition_sizes(p) == [count(==(1), initial), count(==(2), initial)]
                @test _kl_test_cut(g, Q.membership(p)) <= _kl_test_cut(g, initial)
            end
        end
    end

    @testset "Negative moves, rollback, and first best prefix" begin
        initial = [1, 1, 1, 1, 2, 2, 2, 2]
        g = _kl_test_graph(8, [(1, 2), (1, 7), (2, 6), (3, 4), (3, 8), (4, 5),
                               (4, 8), (5, 6), (5, 7), (5, 8), (6, 7), (6, 8), (7, 8)])
        pairs, costs = _kl_test_pass(g, initial)
        @test pairs == [(1, 8), (2, 5), (3, 6), (4, 7)]
        @test costs == [5, 6, 4, 6, 5]
        @test Q.membership(Q.kl_partition(g; initial=initial, max_iter=1)) == [2, 2, 1, 1, 1, 2, 2, 1]

        zero_first = _kl_test_graph(8, [(2, 8), (5, 7), (6, 8)])
        @test last(_kl_test_pass(zero_first, initial)) == [1, 1, 0, 1, 1]
        @test Q.edge_cut_cost(zero_first, Q.kl_partition(zero_first; initial=initial, max_iter=1)) == 0

        tied = _kl_test_graph(6, [(1, 4)])
        start = [1, 1, 1, 2, 2, 2]
        @test last(_kl_test_pass(tied, start)) == [1, 0, 0, 1]
        @test Q.membership(Q.kl_partition(tied; initial=start, max_iter=1)) == [2, 1, 1, 2, 1, 2]
    end

    @testset "Initializers, unequal sizes, and input ownership" begin
        g = path_graph(6)
        labels = [1, 2, 2, 2, 2, 2]
        partition = Q.GraphPartition(labels; nparts=2)
        initializers = (labels, labels .== 1, Dict(i => labels[i] for i in 1:6), partition)
        for initial in initializers
            saved = initial isa Q.GraphPartition ? Q.membership(initial) : copy(initial)
            original = copy(g)
            p = Q.kl_partition(g; initial=initial)
            @test Q.membership(p) == _kl_test_reference(g, labels)
            @test Q.partition_sizes(p) == [1, 5]
            @test g == original
            @test (initial isa Q.GraphPartition ? Q.membership(initial) : initial) == saved
            @test Q.membership(Q.kl_partition(g; initial=initial, max_iter=0)) == labels
        end
        copied = Q.kl_partition(g; initial=labels, max_iter=0)
        labels[1] = 2
        @test Q.membership(copied) == [1, 2, 2, 2, 2, 2]

        plain = _kl_test_graph(7, [(v, v + 1) for v in 1:6]; vertex_type=Int32)
        looped = copy(plain)
        foreach(v -> add_edge!(looped, v, v), vertices(looped))
        start = [1, 1, 1, 2, 2, 2, 2]
        @test Q.membership(Q.kl_partition(looped; initial=start)) == _kl_test_reference(plain, start)
        @test Q.partition_sizes(Q.kl_partition(plain)) == [3, 4]
        @test Q.partition_sizes(Q.kl_partition(path_graph(2))) == [1, 1]
    end

    @testset "Local randomness and validation" begin
        g = path_graph(7)
        expected = fill(2, 7)
        expected[randperm(Xoshiro(91), 7)[1:3]] .= 1
        @test Q.membership(Q.kl_partition(g; seed=91, max_iter=0)) == expected
        Random.seed!(654)
        next_draw = rand(UInt)
        Random.seed!(654)
        first_result = Q.membership(Q.kl_partition(g; seed=91))
        @test rand(UInt) == next_draw
        @test Q.membership(Q.kl_partition(g; seed=91)) == first_result
        @test Q.membership(Q.kl_partition(g; initial=expected, seed=0)) ==
              Q.membership(Q.kl_partition(g; initial=expected, seed=99))
        for n in (0, 1)
            @test_throws ArgumentError Q.kl_partition(SimpleGraph(n))
        end
        for seed in (-1, 1.0, true, false, "1")
            @test_throws ArgumentError Q.kl_partition(g; seed=seed)
        end
        for limit in (-1, 1.0, 1.5, true, false, "1")
            @test_throws ArgumentError Q.kl_partition(g; max_iter=limit)
        end
        for initial in (ones(Int, 7), fill(2, 7), trues(7), [1, 2], [0, 1, 1, 2, 2, 2, 2],
                        [1, 1, 1, 2, 2, 2, 3], Float64.(expected), Dict(1 => 1, 2 => 2),
                        Q.GraphPartition(expected; nparts=3), false)
            @test_throws ArgumentError Q.kl_partition(g; initial=initial)
        end
    end
end
