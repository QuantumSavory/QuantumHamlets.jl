module BuryMulti

using Graphs
using ..BuryHeap

"""
    bury_k_fast(g, k; seeded=false, kernel=BuryHeap._colors)

Experimental balanced multipartition driver. Keep vertices in ascending original
order, ask the kernel only for colors, and assign the final remainder directly.
No production methods are extended or changed. Valid inputs only; diagnostic
logging, return dictionary element types, and invalid-input behavior are not
intended to duplicate the legacy API.
"""
function bury_k_fast(g::Graphs.Graph, k::Integer; seeded::Bool=false,
                     kernel=BuryHeap._colors)
    n = nv(g)
    1 <= k <= n && n % k == 0 || throw(ArgumentError("require 1 <= k <= n and k dividing n"))
    target = n ÷ k
    active = collect(1:n)
    current = g
    result = Dict{Int,Int}()
    sizehint!(result, n)
    for part in 1:k-1
        seed = seeded ? argmin(degree(current)) : nothing
        colored = kernel(current, target; seed=seed)
        retained = Int[]
        sizehint!(retained, length(active) - target)
        for v in eachindex(active)
            if colored[v]
                result[active[v]] = part
            else
                push!(retained, v)
            end
        end
        current, map = induced_subgraph(current, retained)
        active = active[map]
    end
    for v in active
        result[v] = k
    end
    return result
end

end # module BuryMulti
