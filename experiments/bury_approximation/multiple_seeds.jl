using Graphs, QuantumHamlets, Test, Random

module MultiSeedFixtures
include("counterexamples/path_hubs.jl")
end

"""
Experimental multiple-seed adaptation: initially make every vertex in the
union of the seeds' closed neighborhoods eligible, then use the unchanged
scan BURY selection/update rule. Seeds are not precolored. This is separate
from independently grown regions or seeds committed to opposite partitions.

On the path-hub K4 construction (r>=3), every nonempty seed set still gives
matching cost >=ceil(n/6), against optimum one. More generally, any nonempty
initial eligibility region fails: the first choice is a hub, a nonconnector
(fills one clique), or a connector (fills one clique and its hub). An eligible
hub frontier of residual price <=2 then dominates untouched clique candidates
of price >=3. All such neighborhoods fit in the union of the hubs, connectors,
and at most three initial clique interiors: 4r+3<=5r. Thus all hubs/connectors
are colored, and the original split-clique counting proof applies unchanged.
This statement does not apply to seeds preassigned to opposing partitions,
front-specific quotas, or a different growth priority.
"""
function multiple_seed_colors(g, target, seeds)
    target isa Integer && 0 <= target <= nv(g) || throw(ArgumentError("invalid target"))
    isempty(seeds) && throw(ArgumentError("at least one seed is required"))
    all(s -> s isa Integer && 1 <= s <= nv(g), seeds) || throw(ArgumentError("invalid seed"))
    weights = Float64[Graphs.degree(g, v) + 1 for v in vertices(g)]
    colored = fill(false, nv(g))
    region = fill(false, nv(g))
    for seed in seeds
        region[seed] = true
        region[neighbors(g, seed)] .= true
    end
    remaining = target
    ncolored = 0
    while ncolored < target
        cost, choice = QuantumHamlets._bury_scan_minimum(weights, region)
        whole = remaining >= cost
        remaining -= whole ? cost : 1
        ncolored += QuantumHamlets._bury_scan_color!(g, choice, weights, colored)
        region[choice] = true
        if whole
            for v in neighbors(g, choice)
                ncolored += QuantumHamlets._bury_scan_color!(g, v, weights, colored)
                region[v] = true
            end
        end
        weights[choice] = Inf
    end
    return colored
end

function run_multiple_seed_checks()
    @testset "Singleton agrees with production; every vertex agrees with ordinary" begin
        for n in 1:4
            pairs = [(u,v) for u in 1:n for v in u+1:n]
            for mask in 0:((1 << length(pairs)) - 1)
                g = SimpleGraph(n)
                for (i,(u,v)) in enumerate(pairs)
                    iszero(mask & (1 << (i-1))) || add_edge!(g,u,v)
                end
                for target in 0:n
                    @test multiple_seed_colors(g,target,collect(vertices(g))) == last(QuantumHamlets.bury(g,target))
                    for seed in vertices(g)
                        @test multiple_seed_colors(g,target,[seed]) == last(QuantumHamlets.bury_seeded(g,target,seed))
                    end
                end
            end
        end
    end
    rng = MersenneTwister(20261112)
    @testset "Multiple simultaneous seeds retain the counterexample" begin
        println("n\tseed_sets_checked\tminimum_matching\tmaximum_matching\toptimal_matching")
        for r in (3,4,8)
            g,hubs,connectors,_,optimum = MultiSeedFixtures.path_hub_obstruction(r)
            n = nv(g)
            before = collect(edges(g))
            sets = Vector{Int}[[v] for v in vertices(g)]
            append!(sets, [[u,v] for u in 1:n for v in u+1:n])
            if r == 3
                append!(sets, [[u,v,w] for u in 1:n for v in u+1:n for w in v+1:n])
            end
            append!(sets, [hubs,connectors,collect(vertices(g))])
            for _ in 1:100
                push!(sets, randperm(rng,n)[1:rand(rng,1:n)])
            end
            scores = Int[]
            for seeds in sets
                colors = multiple_seed_colors(g,n÷2,seeds)
                @test count(colors) == n÷2
                @test all(colors[hubs])
                @test all(colors[connectors])
                cost = matching_cost(g,GraphPartition(Int.(colors).+1; nparts=2))
                @test cost >= cld(n,6)
                push!(scores,cost)
            end
            @test matching_cost(g,optimum) == 1
            @test collect(edges(g)) == before
            println(join((n,length(sets),minimum(scores),maximum(scores),1),'\t'))
            flush(stdout)
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_multiple_seed_checks()
end
