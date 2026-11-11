module CertifiedBury

using Graphs
using QuantumHamlets

export certified_bury, component_guarded_bury, partition_after_deleting

# All connected components outside S must lie wholly on one side of the cut.
# The vertices in S may be split arbitrarily because they cover every cut edge.
function partition_after_deleting(g, target::Int, removed::Vector{Int})
    n = nv(g)
    excluded = falses(n)
    excluded[removed] .= true
    seen = copy(excluded)
    components = Vector{Int}[]
    for root in 1:n
        seen[root] && continue
        component = Int[]
        stack = [root]
        seen[root] = true
        while !isempty(stack)
            v = pop!(stack)
            push!(component, v)
            for u in neighbors(g, v)
                if !seen[u]
                    seen[u] = true
                    push!(stack, u)
                end
            end
        end
        push!(components, component)
    end

    # A subset of component sizes must sum to t in [target-|S|, target].
    # Entries are assigned once, so predecessor paths remain valid.
    reachable = falses(target + 1)
    reachable[1] = true
    previous = fill(-1, target + 1)
    used_component = zeros(Int, target + 1)
    for (i, component) in enumerate(components)
        size = length(component)
        for t in target:-1:size
            if !reachable[t + 1] && reachable[t - size + 1]
                reachable[t + 1] = true
                previous[t + 1] = t - size
                used_component[t + 1] = i
            end
        end
    end

    for t in max(0, target - length(removed)):target
        reachable[t + 1] || continue
        colored = falses(n)
        cursor = t
        while cursor > 0
            colored[components[used_component[cursor + 1]]] .= true
            cursor = previous[cursor + 1]
        end
        colored[removed[1:(target - t)]] .= true
        return colored
    end
    return nothing
end

function _first_feasible_separator(g, target, size)
    chosen = Int[]
    examined = Ref(0)
    function visit(first)
        if length(chosen) == size
            examined[] += 1
            candidate = partition_after_deleting(g, target, chosen)
            return candidate === nothing ? nothing : (candidate, copy(chosen))
        end
        missing = size - length(chosen)
        for v in first:(nv(g) - missing + 1)
            push!(chosen, v)
            result = visit(v + 1)
            pop!(chosen)
            result === nothing || return result
        end
        return nothing
    end
    return visit(1), examined[]
end

"""
    component_guarded_bury(g, target=nv(g)÷2; seed=nothing, backend=:heap)

Return a zero-cost partition if whole connected components can fill `target`;
otherwise call ordinary BURY, or explicit-seeded BURY when `seed` is supplied.
Return a `GraphPartition`. Validate every argument before the component check.

For simple graphs with maximum degree at most two, this is a matching-cost
2-approximation: the prepass solves every zero-cost instance, and otherwise
BURY cuts at most two edges while the optimum matching cost is at least one.
No approximation guarantee is claimed here for graphs of larger degree.
The subset-sum component check adds O(n*target+n+m) worst-case time.
"""
function component_guarded_bury(g::Graphs.Graph, target=nv(g)÷2;
                                seed=nothing, backend::Symbol=:heap)
    backend in (:heap, :scan, :legacy) ||
        throw(ArgumentError("backend must be :heap, :scan, or :legacy"))
    target isa Integer && 0 <= target <= nv(g) ||
        throw(ArgumentError("target must be an integer in 0:nv(g)"))
    if seed !== nothing
        seed isa Integer && 1 <= seed <= nv(g) ||
            throw(ArgumentError("seed must be an integer in 1:nv(g)"))
    end
    target = Int(target)
    colored = partition_after_deleting(g, target, Int[])
    if colored === nothing
        _, colored = seed === nothing ?
            QuantumHamlets.bury(g, target; backend=backend) :
            QuantumHamlets.bury_seeded(g, target, seed; backend=backend)
    end
    return GraphPartition(Int.(colored) .+ 1; nparts=2)
end

"""
    certified_bury(g, target=nv(g)÷2; max_cover=1, backend=:heap)

Run existing ordinary BURY and search exhaustively for a better cut whose
crossing edges have a vertex cover of size at most `max_cover`. Return a
partition, its exact matching cost, a proved lower bound, and a certified
approximation factor. A zero-cost optimum is reported with factor 1.
Require an integer cap in 0:nv(g); use max_cover=0 for an empty graph.

This is a research prototype, not a general polynomial-time constant-factor
algorithm. For fixed `max_cover=q`, the extra work is bounded by
O(sum(binomial(n,s), s=0:q) * (n+m+n*target)). Exhaustive search is deliberate.
It may be expensive even for q=2 on large inputs. It stops early if it proves
the incumbent optimal. The production BURY implementation is not changed.
"""
function certified_bury(g::Graphs.Graph, target=nv(g)÷2;
                        max_cover=1, backend::Symbol=:heap)
    target isa Integer && 0 <= target <= nv(g) ||
        throw(ArgumentError("target must be an integer in 0:nv(g)"))
    max_cover isa Integer && 0 <= max_cover <= nv(g) ||
        throw(ArgumentError("max_cover must be an integer in 0:nv(g)"))
    target = Int(target)
    _, original = QuantumHamlets.bury(g, target; backend=backend)
    partition = GraphPartition(Int.(original) .+ 1; nparts=2)
    upper = matching_cost(g, partition)
    original_cost = upper
    upper == 0 && return (partition=partition, cost=0, lower_bound=0,
                         factor=1.0, optimal=true, original_cost=0,
                         separator=nothing, examined=0)

    examined = 0
    last_size = min(Int(max_cover), upper - 1)
    for size in 0:last_size
        result, count = _first_feasible_separator(g, target, size)
        examined += count
        if result !== nothing
            colored, separator = result
            partition = GraphPartition(Int.(colored) .+ 1; nparts=2)
            # By Konig's theorem and the exhausted smaller sizes, cost=size.
            cost = matching_cost(g, partition)
            @assert cost == size
            return (partition=partition, cost=cost, lower_bound=cost,
                    factor=1.0, optimal=true, original_cost=original_cost,
                    separator=separator, examined=examined)
        end
    end
    lower = last_size + 1
    return (partition=partition, cost=upper, lower_bound=lower,
            factor=upper / lower, optimal=upper == lower,
            original_cost=original_cost, separator=nothing, examined=examined)
end

end
