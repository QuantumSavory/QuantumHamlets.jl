# Frozen before backend integration. Keep independent of production BURY dispatch.
module OriginalBury
using Graphs, Random
using QuantumHamlets: inverse_map, quantumLand, matching_cost_of_partition

"""
Implementation of k-partition BURY algorithm from our work, QuantumHamlets: https://arxiv.org/abs/2603.06387
This function assumes that k divides the number of vertices in g. \n
Input: a graph from Graphs.jl, and integer k denoting the number of pieces\n
Output: a registry (dictionary mapping villager ids to hamlet ids)
"""
function bury_k(g::Graphs.Graph, k; verbose=false)
    @assert nv(g)%k ==0

    current_graph = copy(g)
    new_reg = Dict()
    inv_map = collect(1:nv(g))
    total_colored = [false for _ in 1:nv(g)]
    for i in 1:k
        _, colored = bury(current_graph, nv(g)÷k, verbose=verbose)

        for v in 1:nv(g)
            v_in_subgraph = inv_map[v]
            if v_in_subgraph != 0
                if colored[v_in_subgraph]
                    new_reg[v] = i
                    total_colored[v] = true
                end
            end
        end

        current_graph, map = induced_subgraph(g, findall(==(false),total_colored))
        inv_map = inverse_map(map, nv(g))
    end

    return new_reg
end

"""
Implementation of BURY algorithm from our work, QuantumHamlets: https://arxiv.org/abs/2603.06387
This version colors a number of nodes in the graph equal to the second parameter in the input, and has k set to 2.
To use BURY for k partition, use the "bury_k" function in this library.\n
Input: A graph to be bipartitioned, and a number of nodes that will be in the first partition \n
Output: A registry (see README example for explanation), and a Boolean vector marking which nodes where colored
(k=2 case has only colored/not colored in the bipartition). For most users, the registry is all you need. 
The vector is primarily metadata for "bury_k".
"""
function bury(g::Graphs.Graph, r0=nv(g)÷2; verbose=false)
    r = r0

    # Node weight initialization
    weights = [0.0 for _ in 1:nv(g)]
    colored = [false for _ in 1:nv(g)]
    for (index,v) in enumerate(vertices(g))
        weights[index] = Graphs.degree(g,v) + 1
    end

    if verbose
        println("Initial weights and colors:")
        println("\nWeights table")
        println(weights)
        println("Colored table")
        println(colored)
    end

    # Global greedy selection
    while sum(colored)<r0
        cost, choice = findmin(weights)

        if verbose
            println("\n\n###########################\nStarting round. r=", r)
            println("Chose vertex ", choice, " for the price of ", cost)
        end

        if r>= cost
            nbd = vcat(choice,neighbors(g,choice))
            r -= cost
        else
            nbd = choice
            r -= 1
        end

        for v in nbd
            if !colored[v]
                weights[v] -= 1
                weights[neighbors(g,v)] .-= 1
                colored[v] = true
            end
        end
        weights[choice] = Inf
        if verbose
            println("\nUpdated weights table")
            println(weights)
            println("Updated Colored table")
            println(colored)
        end
    end

    new_reg = Dict()
    for i in 1:nv(g)
        new_reg[i] = colored[i]+1
    end
    return new_reg, colored
end

"""
    bury_not_balanced_bipartition(g, min_v, max_v;
                                 bury_v=bury, cost_f=matching_cost_of_partition)

Search the inclusive range of bipartition targets `min_v:max_v`, returning
the best `(registry, colored)` produced by `bury_v`. `cost_f(land, g)` must
return a tuple whose first element is the score to minimize. Ties retain
the last candidate. This legacy search retains the known score-bound and
input-validation limitations described in the usage guide.
"""
function bury_not_balanced_bipartition(g, min_v, max_v; bury_v = bury, cost_f = matching_cost_of_partition)
    n = nv(g)
    lowest_cost = n
    best_reg = nothing
    best_colors = nothing
    for v in min_v:max_v
        bury_reg, colors  = bury_v(g, v)
        buryLand = quantumLand(bury_reg, 2, n)
        cost, _ = cost_f(buryLand, g)
        if cost <= lowest_cost
            best_reg = bury_reg
            best_colors = colors
            lowest_cost = cost
        end
    end

    return best_reg, best_colors
end

# Should maybe add a check that the seed can be chosen - not a problem with this minweight seeder
function seeder_1(g::Graphs.Graph)
    # Node weight initialization
    weights = [0.0 for _ in 1:nv(g)]
    for (index,v) in enumerate(vertices(g))
        weights[index] = Graphs.degree(g,v) + 1
    end
    _, index = findmin(weights)
    return index
end

"""Finds the seed vertex that producing the smallest matching"""
function best_seed(g::Graphs.Graph, r0; verbose=false)
    best_score = Inf
    best_seed = 1
    for seed in 1:nv(g)
        if verbose
            println("Trying seed ", seed, " out of ", nv(g))
        end
        reg, _ = bury_seeded(g, r0, seed)
        land = quantumLand(reg, 2, nv(g))
        cost, _ = matching_cost_of_partition(land, g)
        #cost, _ = cut_rank_cost(land, g)

        if cost< best_score
            best_seed = seed
            best_score = cost
            println("Found better. Seed: ", best_seed, ". Cost: ", best_score)
        end
    end
    return best_seed
end

"""Finds the seed vertex, out of sqrt(n) samples, that producing the smallest matching"""
function best_sqrt_seed(g::Graphs.Graph, r0; verbose=false)
    best_score = Inf
    best_seed = 1
    for i in 1:ceil(Int,sqrt(nv(g)))
        seed = rand(vertices(g))
        if verbose
            println("Trying seed ", i, " out of ", ceil(Int,sqrt(nv(g))))
        end
        reg, _ = bury_seeded(g, r0, seed)
        land = quantumLand(reg, 2, nv(g))
        cost, _ = matching_cost_of_partition(land, g)
        #cost, _ = cut_rank_cost(land, g)

        if cost< best_score
            best_seed = seed
            best_score = cost
            println("Found better. Seed: ", best_seed, ". Cost: ", best_score)
        end
    end
    return best_seed
end

"""Finds seeds based on maximum anitpodal vertices of the graphs"""
function antipode_seeds(g::Graphs.Graph, r0)
    dists = floyd_warshall_shortest_paths(g).dists
    antipodes = findall(==(maximum(dists)), dists)
    seed_list = []
    for pair in antipodes
        push!(seed_list, pair[1])
    end
    sort!(unique!(seed_list))

    return seed_list
end

# this version should work for any k
function bury_k_seeded(g::Graphs.Graph, k; verbose=false)
    @assert nv(g)%k ==0

    current_graph = copy(g)
    new_reg = Dict()
    inv_map = collect(1:nv(g))
    total_colored = [false for _ in 1:nv(g)]
    for i in 1:k
        returned_reg, colored = bury_seeded(current_graph, nv(g)÷k, seeder_1(current_graph), verbose=verbose)
        #returned_reg, colored = bury_seeded(current_graph, nv(g)÷k, verbose=verbose)

        for v in 1:nv(g)
            v_in_subgraph = inv_map[v]
            if v_in_subgraph != 0
                if colored[v_in_subgraph]
                    new_reg[v] = i
                    total_colored[v] = true
                end
            end
        end

        current_graph, map = induced_subgraph(g, findall(==(false),total_colored))
        inv_map = inverse_map(map, nv(g))
    end

    return new_reg
end
#  
function bury_seeded(g::Graphs.Graph, r0=nv(g)÷2, seed =  best_sqrt_seed(g, r0); verbose=false)
    r = r0

    # Node weight initialization
    weights = [0.0 for _ in 1:nv(g)]
    colored = [false for _ in 1:nv(g)]

    region = [false for _ in 1:nv(g)]
    region[neighbors(g,seed)] .= true
    region[seed] = true

    for (index,v) in enumerate(vertices(g))
        weights[index] = Graphs.degree(g,v) + 1
    end

    if verbose
        println("Initial weights and colors:")
        println("\nWeights table")
        println(weights)
        println("Colored table")
        println(colored)
    end

    # Global greedy selections
    while sum(colored)<r0
        cost, choice = findmin(.!region*Inf + weights)
        if cost == Inf
            cost, choice = findmin(weights)
        end

        if verbose
            println("\n\n###########################\nStarting round. r=", r)
            println("Chose vertex ", choice, " for the price of ", cost)
        end

        if r>= cost
            nbd = vcat(choice,neighbors(g,choice))
            r -= cost
        else
            nbd = choice
            r -= 1
        end

        for v in nbd
            if !colored[v]
                weights[v] -= 1
                weights[neighbors(g,v)] .-= 1
                colored[v] = true
            end
            region[v] = true
        end
        weights[choice] = Inf

        if verbose
            println("\nUpdated weights table")
            println(weights)
            println("Updated Colored table")
            println(colored)
        end
    end

    new_reg = Dict()
    for i in 1:nv(g)
        new_reg[i] = colored[i]+1
    end
    return new_reg, colored
end

end # module OriginalBury
