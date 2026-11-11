using Graphs, QuantumHamlets, Random, Test
include("certified_bury.jl")
using .CertifiedBury

# Independent exhaustive matching oracle. It enumerates whether each left
# vertex is unmatched or paired to an unused right neighbor, without relying
# on the package's IGraphs matching routine or a vertex-cover theorem.
function matching_oracle(g, colored)
    left = findall(colored)
    used = falses(nv(g))
    function search(i)
        i > length(left) && return 0
        best = search(i + 1)
        for v in neighbors(g, left[i])
            if !colored[v] && !used[v]
                used[v] = true
                best = max(best, 1 + search(i + 1))
                used[v] = false
            end
        end
        return best
    end
    return search(1)
end

function optimum_oracle(g, target)
    best = nv(g)
    for mask in 0:(2^nv(g) - 1)
        count_ones(mask) == target || continue
        colored = [!iszero(mask & (1 << (v - 1))) for v in vertices(g)]
        best = min(best, matching_oracle(g, colored))
    end
    return best
end

function check_case(g, target, opt, cap)
    before = collect(edges(g))
    result = certified_bury(g, target; max_cover=cap)
    labels = membership(result.partition)
    @test count(==(2), labels) == target
    @test result.cost == matching_oracle(g, labels .== 2)
    @test result.lower_bound <= opt <= result.cost <= result.original_cost
    @test result.factor == (result.cost == 0 ? 1.0 : result.cost / result.lower_bound)
    @test result.optimal == (result.lower_bound == result.cost)
    @test !result.optimal || result.cost == opt
    @test opt > cap || result.optimal
    @test collect(edges(g)) == before
    return nothing
end

@testset "Certified matching BURY: every labeled graph through five vertices" begin
    for n in 0:5
        pairs = [(u, v) for u in 1:n for v in (u + 1):n]
        for mask in 0:(2^length(pairs) - 1)
            g = SimpleGraph(n)
            for (i, (u, v)) in enumerate(pairs)
                iszero(mask & (1 << (i - 1))) || add_edge!(g, u, v)
            end
            for target in 0:n
                opt = optimum_oracle(g, target)
                for cap in 0:min(n, 2)
                    check_case(g, target, opt, cap)
                end
            end
        end
    end
end

@testset "Random balanced graphs and larger cover cap" begin
    rng = MersenneTwister(20261111)
    for n in (6, 8), trial in 1:40
        g = SimpleGraph(n)
        density = rand(rng)
        for u in 1:n, v in (u + 1):n
            rand(rng) < density && add_edge!(g, u, v)
        end
        target = n ÷ 2
        opt = optimum_oracle(g, target)
        for cap in 0:3
            check_case(g, target, opt, cap)
        end
    end
end

@testset "Boundary validation" begin
    @test_throws ArgumentError certified_bury(path_graph(4), -1)
    @test_throws ArgumentError certified_bury(path_graph(4), 5)
    @test_throws ArgumentError certified_bury(path_graph(4), 2.0)
    @test_throws ArgumentError certified_bury(path_graph(4); max_cover=-1)
    @test_throws ArgumentError certified_bury(path_graph(4); max_cover=5)
    @test_throws ArgumentError certified_bury(path_graph(4); max_cover=1.0)
end
