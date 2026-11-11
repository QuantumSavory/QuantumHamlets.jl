using Graphs, QuantumHamlets, Test, Random
const Q = QuantumHamlets

"""
Two equal modules: one hub and r cliques of size b, each attached to its hub
by a designated connector. The hubs are joined by one bridge. Each module
has r*b+1 vertices. No labeling assumption is needed for the obstruction.
"""
function all_seed_obstruction(r::Integer; b::Integer=4r)
    r >= 2 || throw(ArgumentError("r must be at least two"))
    b >= r+3 || throw(ArgumentError("b must be at least r+3"))
    module_size = r*b+1
    g = SimpleGraph(2module_size)
    hubs = [1, module_size+1]
    connectors = Int[]
    cliques = Vector{Int}[]
    labels = zeros(Int, nv(g))
    for mod in 1:2
        offset = (mod-1)*module_size
        labels[offset+1:offset+module_size] .= mod
        for j in 1:r
            clique = collect(offset+2+(j-1)*b:offset+1+j*b)
            push!(cliques, clique)
            push!(connectors, first(clique))
            add_edge!(g, hubs[mod], first(clique))
            for u in clique, v in clique
                u < v && add_edge!(g, u, v)
            end
        end
    end
    add_edge!(g, hubs...)
    return g, hubs, connectors, cliques, GraphPartition(labels)
end

function all_seed_checks()
    rng = MersenneTwister(314159)
    println("r\tb\tn\tordinary_matching\tordinary_rank\tordinary_cut\tseed_matching_min\tseed_matching_max\tseed_rank_min\tseed_rank_max")
    @testset "All-seed BURY obstruction" begin
        for r in (2, 3, 4, 8)
            b = 4r
            g, hubs, connectors, cliques, optimum = all_seed_obstruction(r; b)
            target = nv(g) ÷ 2
            @test is_connected(g)
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
                mc = matching_cost(g, p)
                rc = pairwise_cut_rank(g, p)
                ec = edge_cut_cost(g, p)
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
            seeds = r <= 4 ? collect(vertices(g)) : [hubs; connectors; 3; nv(g)]
            for seed in seeds
                _, colored = Q.bury_seeded(g, target, seed)
                mc, rc, _ = check_output(colored)
                push!(mc_list, mc)
                push!(rc_list, rc)
                if seed in (1, 2, 3, target+1, target+2, target+3)
                    for backend in (:scan, :legacy)
                        _, alternative = Q.bury_seeded(g, target, seed; backend)
                        @test alternative == colored
                    end
                end
            end
            # Verify that ID changes do not rescue the rule on bounded cases.
            if r <= 3
                for trial in 1:3
                    permutation = randperm(rng, nv(g))
                    permuted, _ = induced_subgraph(g, permutation)
                    for seed in vertices(permuted)
                        _, colored = Q.bury_seeded(permuted, target, seed)
                        original_colors = falses(nv(g))
                        original_colors[permutation] = colored
                        check_output(original_colors)
                    end
                end
            end
            println(join((r, b, nv(g), ordinary_costs..., minimum(mc_list), maximum(mc_list),
                          minimum(rc_list), maximum(rc_list)), '\t'))
            flush(stdout)
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    all_seed_checks()
end
