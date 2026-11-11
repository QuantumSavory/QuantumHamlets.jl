# Allocation-light minimum scans. Only the chosen vertex is processed;
# colored neighbors remain candidates until selected in a later round.

@inline function _bury_scan_color!(g, v, weights, colored)
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

function _bury_scan_minimum(weights, ::Nothing)
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

function _bury_scan_minimum(weights, region::Vector{Bool})
    cost = Inf
    choice = 1
    # Strict comparison keeps the smallest vertex ID when weights tie.
    for v in eachindex(weights)
        if region[v] && weights[v] < cost
            cost = weights[v]
            choice = v
        end
    end
    # An exhausted seed region falls back to the global minimum, including
    # when the remaining vertices belong to disconnected components.
    return cost == Inf ? _bury_scan_minimum(weights, nothing) : (cost, choice)
end

function _bury_scan_colors(g::Graphs.Graph, r0::Integer;
                           seed::Union{Nothing,Integer}=nothing, verbose=false)
    n = nv(g)
    0 <= r0 <= n || throw(ArgumentError("target must lie in 0:nv(g)"))
    isnothing(seed) || 1 <= seed <= n || throw(ArgumentError("invalid seed"))

    weights = Float64[Graphs.degree(g, v) + 1 for v in vertices(g)]
    colored = fill(false, n)
    region = if isnothing(seed)
        nothing
    else
        eligible = fill(false, n)
        eligible[Int(seed)] = true
        for v in neighbors(g, Int(seed))
            eligible[v] = true
        end
        eligible
    end

    verbose && _bury_trace_initial(weights, colored)
    # Keep the initial integer budget for the legacy trace. It becomes a
    # Float64 only after subtracting a neighborhood's floating-point cost.
    remaining = r0
    ncolored = 0
    while ncolored < r0
        cost, choice = _bury_scan_minimum(weights, region)
        verbose && _bury_trace_choice(remaining, choice, cost)
        whole_neighborhood = remaining >= cost
        remaining -= whole_neighborhood ? cost : 1

        # Preserve the reference update order: choice, then its neighbors.
        ncolored += _bury_scan_color!(g, choice, weights, colored)
        isnothing(region) || (region[choice] = true)
        if whole_neighborhood
            for v in neighbors(g, choice)
                ncolored += _bury_scan_color!(g, v, weights, colored)
                isnothing(region) || (region[v] = true)
            end
        end

        weights[choice] = Inf
        verbose && _bury_trace_updated(weights, colored)
    end
    return colored
end
