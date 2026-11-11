using Graphs, QuantumHamlets, Test
const Q = QuantumHamlets

"""
Two equal connected modules. Each has t copies of K4,4 minus its diagonal,
joined along A vertices. Every core vertex has a private leaf. The modules
are joined by one B--B bridge. IDs put all A leaves, then all B leaves,
then all A cores, then all B cores.
"""
function obstruction_graph(t::Integer)
    t >= 1 || throw(ArgumentError("t must be positive"))
    q = 8t                   # number of A vertices across both modules
    g = SimpleGraph(4q)
    a_leaf(mod, block, i) = (mod - 1) * 4t + 4(block - 1) + i
    b_leaf(mod, block, i) = q + a_leaf(mod, block, i)
    a_core(mod, block, i) = 2q + a_leaf(mod, block, i)
    b_core(mod, block, i) = 3q + a_leaf(mod, block, i)
    module_labels = zeros(Int, nv(g))
    for mod in 1:2, block in 1:t, i in 1:4
        a, b = a_core(mod, block, i), b_core(mod, block, i)
        add_edge!(g, a, a_leaf(mod, block, i))
        add_edge!(g, b, b_leaf(mod, block, i))
        for j in 1:4
            i == j || add_edge!(g, a, b_core(mod, block, j))
        end
        for v in (a, b, a_leaf(mod, block, i), b_leaf(mod, block, i))
            module_labels[v] = mod
        end
    end
    for mod in 1:2, block in 1:t-1
        add_edge!(g, a_core(mod, block, 1), a_core(mod, block+1, 1))
    end
    add_edge!(g, b_core(1, 1, 1), b_core(2, 1, 1))
    expected = vcat(fill(true, q), fill(false, q), fill(true, q), fill(false, q))
    return g, expected, GraphPartition(module_labels)
end

function check_family()
    println("t\tn\tm\tmaxdegree\tbury_matching\tbury_rank\tbury_cut\topt_witness\tseed_min\tseed_max")
    @testset "Connected bounded-degree BURY obstruction" begin
        for t in (1, 2, 3, 4, 8)
            g, expected, optimum = obstruction_graph(t)
            n = nv(g)
            @test is_connected(g)
            @test maximum(Graphs.degree(g)) <= 6
            @test partition_sizes(optimum) == [n ÷ 2, n ÷ 2]
            @test matching_cost(g, optimum) == 1
            @test pairwise_cut_rank(g, optimum) == 1
            @test edge_cut_cost(g, optimum) == 1
            observed = nothing
            for backend in (:heap, :scan, :legacy)
                reg, colored = Q.bury(g, n ÷ 2; backend)
                @test colored == expected
                p = GraphPartition(reg; nparts=2)
                @test partition_sizes(p) == [n ÷ 2, n ÷ 2]
                @test matching_cost(g, p) == 8t
                @test pairwise_cut_rank(g, p) == 8t
                @test edge_cut_cost(g, p) == 24t
                p2 = partition_graph(g, 2; backend)
                # Multipartition labels are reversed relative to bury.
                @test membership(p2) == 3 .- membership(p)
                @test matching_cost(g, p2) == 8t
                observed = p
            end
            seed_costs = Int[]
            # Bound all-seed study to the first three small family members.
            if t <= 3
                for seed in vertices(g)
                    reg, _ = Q.bury_seeded(g, n ÷ 2, seed; backend=:heap)
                    push!(seed_costs, matching_cost(g, GraphPartition(reg; nparts=2)))
                end
            end
            seed_min = isempty(seed_costs) ? "not_run" : string(minimum(seed_costs))
            seed_max = isempty(seed_costs) ? "not_run" : string(maximum(seed_costs))
            println(join((t, n, ne(g), maximum(Graphs.degree(g)), matching_cost(g, observed),
                          pairwise_cut_rank(g, observed), edge_cut_cost(g, observed), 1,
                          seed_min, seed_max), '\t'))
            flush(stdout)
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    check_family()
end
