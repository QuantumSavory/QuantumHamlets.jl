module BuryScan

using Graphs

export bury_scan

# Isolated experiment: no QuantumHamlets definitions are extended or replaced.
# Supports valid simple undirected graphs with vertices 1:n, integer targets
# 0:n, and either no seed or one explicitly chosen vertex in 1:n.

@inline function _color!(g, v, weights, colored)
    if !colored[v]
        weights[v] -= 1
        for neighbor in neighbors(g, v)
            weights[neighbor] -= 1
        end
        colored[v] = true
        return 1
    end
    return 0
end

function _minimum(weights, region::Nothing)
    cost = Inf
    choice = 1
    for v in eachindex(weights)
        if weights[v] < cost
            cost = weights[v]
            choice = v
        end
    end
    return cost, choice
end

function _minimum(weights, region::Vector{Bool})
    cost = Inf
    choice = 1
    # Strict comparison retains the first (smallest-ID) minimizer.
    for v in eachindex(weights)
        if region[v] && weights[v] < cost
            cost = weights[v]
            choice = v
        end
    end
    # Match the legacy seeded fallback, including disconnected components.
    return cost == Inf ? _minimum(weights, nothing) : (cost, choice)
end

"""
    _colors(g, r0; seed=nothing)

Color-only BURY kernel for the isolated performance experiment. The returned
`Vector{Bool}` follows the original algorithm's choices and tie order.
"""
function _colors(g::Graphs.Graph, r0::Integer; seed::Union{Nothing,Int}=nothing)
    n = nv(g)
    0 <= r0 <= n || throw(ArgumentError("target must lie in 0:nv(g)"))
    isnothing(seed) || 1 <= seed <= n || throw(ArgumentError("invalid seed"))

    weights = Float64[degree(g, v) + 1 for v in vertices(g)]
    colored = fill(false, n)
    region = if isnothing(seed)
        nothing
    else
        eligible = fill(false, n)
        eligible[seed] = true
        for v in neighbors(g, seed)
            eligible[v] = true
        end
        eligible
    end

    remaining = Float64(r0)
    ncolored = 0
    while ncolored < r0
        cost, choice = _minimum(weights, region)
        whole_neighborhood = remaining >= cost
        remaining -= whole_neighborhood ? cost : 1

        # Preserve the original update order: choice, then its neighbors.
        ncolored += _color!(g, choice, weights, colored)
        isnothing(region) || (region[choice] = true)
        if whole_neighborhood
            for v in neighbors(g, choice)
                ncolored += _color!(g, v, weights, colored)
                isnothing(region) || (region[v] = true)
            end
        end

        # Colored and processed are different states. Only the chosen vertex
        # becomes unavailable; newly colored neighbors may still be chosen.
        weights[choice] = Inf
    end
    return colored
end

"""
    bury_scan(g, r0=nv(g)÷2; seed=nothing)

Return `(registry, colored)` with the same contents as ordinary `bury`, or
`bury_seeded(g, r0, seed)` when a seed is supplied, on supported valid inputs.
The registry is concretely typed as `Dict{Int,Int}`. Logging is not reproduced.
"""
function bury_scan(g::Graphs.Graph, r0::Integer=nv(g) ÷ 2;
                   seed::Union{Nothing,Int}=nothing)
    colored = _colors(g, r0; seed=seed)
    registry = Dict{Int,Int}()
    sizehint!(registry, length(colored))
    for v in eachindex(colored)
        registry[v] = Int(colored[v]) + 1
    end
    return registry, colored
end

end # module BuryScan
