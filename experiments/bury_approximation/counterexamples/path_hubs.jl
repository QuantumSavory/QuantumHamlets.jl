using Graphs, QuantumHamlets, Test, Random
const Q = QuantumHamlets

"""
A path of 2r hub vertices, each attached to a separate K_b by one connector.
The two halves of the hub path (with their cliques) form equally sized modules
joined by one edge. Hubs have IDs1:2r; each following group of b IDs is a clique.
"""
function path_hub_obstruction(r::Integer; b::Integer=4)
    r >= 2 || throw(ArgumentError("r must be at least two"))
    b >= 4 || throw(ArgumentError("b must be at least four"))
    b == 4 && r < 3 && throw(ArgumentError("b=4 requires r at least three"))
    g = SimpleGraph(2r*(b+1))
    hubs = collect(1:2r)
    cliques = [collect(2r+1+(j-1)*b:2r+j*b) for j in 1:2r]
    connectors = first.(cliques)
    labels = zeros(Int, nv(g))
    for j in 1:2r
        j < 2r && add_edge!(g, j, j+1)
        add_edge!(g, j, connectors[j])
        labels[j] = j <= r ? 1 : 2
        labels[cliques[j]] .= labels[j]
        for u in cliques[j], v in cliques[j]
            u < v && add_edge!(g, u, v)
        end
    end
    return g, hubs, connectors, cliques, GraphPartition(labels)
end

function path_hub_checks()
    rng = MersenneTwister(271828)
    println("r\tb\tn\tordinary_matching\tordinary_rank\tordinary_cut\tseed_matching_min\tseed_matching_max\tseed_rank_min\tseed_rank_max")
    @testset "Bounded-degree all-seed BURY obstruction" begin
        for (r, b) in [(r,b) for b in (4,5,6) for r in (2,3,4,8) if b != 4 || r >= 3]
            g, hubs, connectors, cliques, optimum = path_hub_obstruction(r; b)
            target = nv(g) ÷ 2
            @test is_connected(g)
            @test maximum(Graphs.degree(g)) == b
            @test partition_sizes(optimum) == [target, target]
            @test matching_cost(g, optimum) == 1
            @test pairwise_cut_rank(g, optimum) == 1
            @test edge_cut_cost(g, optimum) == 1
            function check_output(colored)
                @test count(colored) == target
                @test all(colored[hubs])
                @test all(colored[connectors])
                partial = count(c -> 0 < count(colored[c]) < b, cliques)
                @test partial >= r+1
                p = GraphPartition(Int.(colored) .+ 1; nparts=2)
                mc, rc, ec = matching_cost(g, p), pairwise_cut_rank(g, p), edge_cut_cost(g, p)
                @test mc >= r+1
                @test rc == partial
                @test ec >= (r+1)*(b-1)
                return mc, rc, ec
            end
            _, ordinary = Q.bury(g, target)
            ordinary_costs = check_output(ordinary)
            for backend in (:scan, :legacy)
                _, colored = Q.bury(g, target; backend)
                @test colored == ordinary
                check_output(colored)
            end
            mc_list, rc_list = Int[], Int[]
            for seed in vertices(g)
                _, colored = Q.bury_seeded(g, target, seed)
                mc, rc, _ = check_output(colored)
                push!(mc_list, mc)
                push!(rc_list, rc)
                if seed in (1, 2, connectors[1], connectors[1]+1, nv(g))
                    for backend in (:scan, :legacy)
                        _, alternative = Q.bury_seeded(g, target, seed; backend)
                        @test alternative == colored
                    end
                end
            end
            for trial in 1:2
                permutation = randperm(rng, nv(g))
                permuted, _ = induced_subgraph(g, permutation)
                for seed in vertices(permuted)
                    _, colored = Q.bury_seeded(permuted, target, seed)
                    original_colors = falses(nv(g))
                    original_colors[permutation] = colored
                    check_output(original_colors)
                end
            end
            println(join((r, b, nv(g), ordinary_costs..., minimum(mc_list), maximum(mc_list),
                          minimum(rc_list), maximum(rc_list)), '\t'))
            flush(stdout)
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    path_hub_checks()
end
