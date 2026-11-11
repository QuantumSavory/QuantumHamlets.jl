module BurySearch

using Graphs, Random, QuantumHamlets
using ..BuryHeap

const Q = QuantumHamlets

# Preserve candidate order, sampled seeds (including duplicates), score routine,
# strict tie rule, and legacy improvement messages. The BURY kernel and typed
# registry construction change; matching evaluation itself stays unchanged.
function best_seed_fast(g, target; sampled=false, kernel=BuryHeap._colors)
    best_score = Inf
    chosen = 1
    attempts = sampled ? ceil(Int, sqrt(nv(g))) : nv(g)
    for i in 1:attempts
        seed = sampled ? rand(vertices(g)) : i
        colors = kernel(g, target; seed=seed)
        reg = Dict(v => Int(colors[v]) + 1 for v in vertices(g))
        land = Q.quantumLand(reg, 2, nv(g))
        cost, _ = Q.matching_cost_of_partition(land, g)
        if cost < best_score
            chosen = seed
            best_score = cost
            println("Found better. Seed: ", chosen, ". Cost: ", best_score)
        end
    end
    return chosen
end

function bury_auto_fast(g, target=nv(g)÷2; kernel=BuryHeap._colors)
    seed = best_seed_fast(g, target; sampled=true, kernel=kernel)
    colors = kernel(g, target; seed=seed)
    return Dict(v => Int(colors[v]) + 1 for v in vertices(g)), colors
end

function bury_range_fast(g, min_v, max_v; kernel=BuryHeap._colors)
    lowest_cost = nv(g) # Deliberately retain the original bound and tie rule.
    best_reg = nothing
    best_colors = nothing
    for target in min_v:max_v
        colors = kernel(g, target)
        reg = Dict(v => Int(colors[v]) + 1 for v in vertices(g))
        land = Q.quantumLand(reg, 2, nv(g))
        cost, _ = Q.matching_cost_of_partition(land, g)
        if cost <= lowest_cost
            best_reg, best_colors, lowest_cost = reg, colors, cost
        end
    end
    return best_reg, best_colors
end

end # module BurySearch
