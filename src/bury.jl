# Public BURY dispatch and shared diagnostics. Algorithm bodies live in bury/.
function _bury_check_backend(backend::Symbol)
    backend in (:heap, :scan, :legacy) ||
        throw(ArgumentError("unknown BURY backend $backend; use :heap, :scan, or :legacy"))
    return backend
end

function _bury_check_target(g, target)
    target isa Integer && 0 <= target <= nv(g) ||
        throw(ArgumentError("BURY target must be an integer in 0:nv(g)"))
    return Int(target)
end

function _bury_check_seed(g, seed)
    seed isa Integer && 1 <= seed <= nv(g) ||
        throw(ArgumentError("BURY seed must be an integer in 1:nv(g)"))
    return Int(seed)
end

function _bury_trace_initial(weights, colored)
    println("Initial weights and colors:")
    println("\nWeights table")
    println(weights)
    println("Colored table")
    println(colored)
end

function _bury_trace_choice(remaining, choice, cost)
    println("\n\n###########################\nStarting round. r=", remaining)
    println("Chose vertex ", choice, " for the price of ", cost)
end

function _bury_trace_updated(weights, colored)
    println("\nUpdated weights table")
    println(weights)
    println("Updated Colored table")
    println(colored)
end

include("bury/legacy.jl")
include("bury/scan.jl")
include("bury/heap.jl")

function _bury_color_kernel(backend::Symbol)
    backend === :heap && return _bury_heap_colors
    backend === :scan && return _bury_scan_colors
    throw(ArgumentError("color-only BURY kernels are :heap and :scan"))
end

function _bury_registry(colored)
    # Preserve the original public dictionary type and insertion order.
    reg = Dict{Any,Any}()
    for v in eachindex(colored)
        reg[v] = Int(colored[v]) + 1
    end
    return reg
end

# These kernels depend on equality/order of the two labels, so Bool membership
# is equivalent to labels 1/2 and needs no intermediate dense Int vector.
_bury_cut_rank(colored, g) = _pairwise_cut_rank(colored, 2, g)

function _bury_range_cost(cost_f)
    cost_f === matching_cost_of_partition && return _matching_cost
    cost_f === naive_cost_of_partition && return _edge_cut_cost
    cost_f === cut_rank_cost && return _bury_cut_rank
    return nothing
end

# Ordinary BURY only. Each target needs its own run: target results need not
# be nested. The legacy algorithm retains its original registry construction.
function _bury_select_colors(g, targets, backend, score)
    lowest_cost = nothing
    best_colors = nothing
    for target in targets
        colors = if backend === :legacy
            last(_bury_legacy(g, target))
        else
            _bury_color_kernel(backend)(g, target)
        end
        cost = score(colors, g)
        if best_colors === nothing || cost <= lowest_cost
            best_colors = copy(colors)
            lowest_cost = cost
        end
    end
    return best_colors
end

"""
    bury_k(g, k; backend=:heap, verbose=false)

Partition nonempty `g` into `k` equal parts with repeated ordinary BURY calls.
Require integer `k` in `1:nv(g)` dividing `nv(g)`. Forward `backend` and `verbose`
to every round; return the original dictionary mapping vertices to parts.
Graph extraction and vertex mapping retain their existing behavior.
"""
function bury_k(g::Graphs.Graph, k; backend::Symbol=:heap, verbose=false)
    _bury_check_backend(backend)
    k = _check_partition_count(g, k)

    current_graph = copy(g)
    new_reg = Dict()
    inv_map = collect(1:nv(g))
    total_colored = [false for _ in 1:nv(g)]
    for i in 1:k
        _, colored = bury(current_graph, nv(g)÷k; backend=backend, verbose=verbose)

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
    bury(g, target=nv(g)÷2; backend=:heap, verbose=false)

Color `target` vertices of a simple undirected graph using BURY. Return
`(registry, colored)`, with a `Dict{Any,Any}` registry and a `Vector{Bool}`.
Uncolored/colored vertices receive part labels 1/2. The graph is not modified.

Choose `:heap` (default indexed heap), `:scan` (optimized linear scan), or
`:legacy` (original implementation). Supported inputs retain the same greedy
choices and smallest-vertex-ID ties. Heap favors larger sparse workloads;
scan can help small or dense workloads, particularly seeded calls.

Every backend requires an integer target in `0:nv(g)`; an empty graph supports
target zero. Invalid inputs raise `ArgumentError` before starting the algorithm.
`verbose=true` prints the original round-by-round trace with every backend.
"""
function bury(g::Graphs.Graph, r0=nv(g)÷2; backend::Symbol=:heap, verbose=false)
    _bury_check_backend(backend)
    target = _bury_check_target(g, r0)
    backend === :legacy && return _bury_legacy(g, target; verbose=verbose)
    colored = _bury_color_kernel(backend)(g, target; verbose=verbose)
    return _bury_registry(colored), colored
end

"""
    bury_not_balanced_bipartition(g, min_v, max_v;
        backend=:heap, bury_v=bury, cost_f=matching_cost_of_partition)

Search the inclusive range of bipartition targets `min_v:max_v`, returning
the best `(registry, colored)` produced by `bury_v`. `cost_f(land, g)` must
return a tuple whose first element is the score to minimize. Ties retain
the last candidate. Built-in `bury_v=bury` or `bury_seeded` receives the
selected backend for every candidate. Custom callbacks keep the two-argument
call `bury_v(g, target)` and control their own backend choice. Require integer
bounds satisfying `0 <= min_v <= max_v <= nv(g)`, checked before any callback.
Retained callback outputs are copied, so later buffer reuse does not change
the winner. Ordinary BURY with the built-in matching, edge-cut, or rank cost
uses scalar kernels without constructing a legacy land for each target.

With `bury_v=bury_seeded`, each target's sampled seed search always selects
by matching cost. The outer `cost_f` therefore selects among matching-selected
candidates, even when it scores rank or edge cut. Use [`bipartition_graph`](@ref)
for ordinary BURY with a graph-first scalar objective and nonempty parts.
"""
function bury_not_balanced_bipartition(g, min_v, max_v; backend::Symbol=:heap,
                                      bury_v=bury, cost_f=matching_cost_of_partition)
    _bury_check_backend(backend)
    min_v = _bury_check_target(g, min_v)
    max_v = _bury_check_target(g, max_v)
    min_v <= max_v || throw(ArgumentError("minimum target must not exceed maximum target"))
    score = _bury_range_cost(cost_f)
    if bury_v === bury && score !== nothing
        colors = _bury_select_colors(g, min_v:max_v, backend, score)
        return _bury_registry(colors), colors
    end
    partitioner = if bury_v === bury
        (g, target) -> bury(g, target; backend=backend)
    elseif bury_v === bury_seeded
        (g, target) -> bury_seeded(g, target; backend=backend)
    else
        bury_v
    end
    n = nv(g)
    lowest_cost = nothing
    best_reg = nothing
    best_colors = nothing
    for v in min_v:max_v
        bury_reg, colors = partitioner(g, v)
        buryLand = quantumLand(bury_reg, 2, n)
        cost, _ = cost_f(buryLand, g)
        if best_reg === nothing || cost <= lowest_cost
            best_reg = copy(bury_reg)
            best_colors = copy(colors)
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

"""
    best_seed(g, target; backend=:heap, verbose=false)

Evaluate every vertex as a seed using the selected BURY backend and original
matching score. Return the first seed achieving the lowest score. Improvement
messages are printed regardless of `verbose`, as in the original function.
Require a nonempty graph and an integer target in `0:nv(g)`.
"""
function best_seed(g::Graphs.Graph, r0; backend::Symbol=:heap, verbose=false)
    _bury_check_backend(backend)
    target = _bury_check_target(g, r0)
    nv(g) > 0 || throw(ArgumentError("seed search requires a nonempty graph"))
    best_score = Inf
    best_seed = 1
    for seed in 1:nv(g)
        if verbose
            println("Trying seed ", seed, " out of ", nv(g))
        end
        reg, _ = bury_seeded(g, target, seed; backend=backend)
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

"""
    best_sqrt_seed(g, target; backend=:heap, verbose=false,
                   rng=Random.default_rng(), log_improvements=true)

Evaluate `ceil(sqrt(nv(g)))` sampled seeds using the selected BURY backend.
Sampling uses `rng` with replacement (the default RNG unless supplied).
Matching scores and first-best ties retain their original behavior.
Set `log_improvements=false` to suppress improvement messages.
Require a nonempty graph and an integer target in `0:nv(g)`; validate before
sampling so invalid calls leave the RNG unchanged.
"""
function best_sqrt_seed(g::Graphs.Graph, r0; backend::Symbol=:heap, verbose=false,
                        rng::AbstractRNG=Random.default_rng(), log_improvements::Bool=true)
    _bury_check_backend(backend)
    target = _bury_check_target(g, r0)
    nv(g) > 0 || throw(ArgumentError("seed search requires a nonempty graph"))
    best_score = Inf
    best_seed = 1
    for i in 1:ceil(Int,sqrt(nv(g)))
        seed = rand(rng, vertices(g))
        if verbose
            println("Trying seed ", i, " out of ", ceil(Int,sqrt(nv(g))))
        end
        reg, _ = bury_seeded(g, target, seed; backend=backend)
        land = quantumLand(reg, 2, nv(g))
        cost, _ = matching_cost_of_partition(land, g)
        #cost, _ = cut_rank_cost(land, g)

        if cost< best_score
            best_seed = seed
            best_score = cost
            log_improvements && println("Found better. Seed: ", best_seed, ". Cost: ", best_score)
        end
    end
    return best_seed
end

"""
    bury_k_seeded(g, k; backend=:heap, verbose=false)

Partition nonempty `g` into `k` equal parts using seeded BURY. Require integer
`k` in `1:nv(g)` dividing `nv(g)`. Every round uses the original minimum-degree
seed rule and the selected backend. Return the original registry dictionary.
"""
function bury_k_seeded(g::Graphs.Graph, k; backend::Symbol=:heap, verbose=false)
    _bury_check_backend(backend)
    k = _check_partition_count(g, k)

    current_graph = copy(g)
    new_reg = Dict()
    inv_map = collect(1:nv(g))
    total_colored = [false for _ in 1:nv(g)]
    for i in 1:k
        _, colored = bury_seeded(current_graph, nv(g)÷k, seeder_1(current_graph); backend=backend, verbose=verbose)

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
    bury_seeded(g, target=nv(g)÷2[, seed]; backend=:heap, verbose=false)

Seeded BURY with the backend and return conventions of [`bury`](@ref). An
explicit positional seed starts the original regional rule; it need not be
selected first. Omit the seed to run [`best_sqrt_seed`](@ref), using the chosen
backend for every candidate and the final partition. Sampling order, duplicate
samples, scoring, and tie policies are unchanged.

Every backend requires a nonempty graph, an integer target in `0:nv(g)`, and
an explicit integer seed in `1:nv(g)`. Invalid inputs raise `ArgumentError`
before sampling or running the algorithm. Automatic search's improvement
messages are preserved, including when `verbose=false`. Automatic searches
also accept `rng=Random.default_rng()` and `log_improvements=false` for a
caller-owned random stream and quiet operation.
"""
function bury_seeded(g::Graphs.Graph, r0=nv(g)÷2;
                     backend::Symbol=:heap, verbose=false,
                     rng::AbstractRNG=Random.default_rng(), log_improvements::Bool=true)
    _bury_check_backend(backend)
    target = _bury_check_target(g, r0)
    # Resolve omitted seeds here so trial calls receive the selected backend.
    # Outer verbose still controls only the final trace, as in the original.
    seed = best_sqrt_seed(g, target; backend=backend, rng=rng, log_improvements=log_improvements)
    return bury_seeded(g, target, seed; backend=backend, verbose=verbose)
end

function bury_seeded(g::Graphs.Graph, r0, seed;
                     backend::Symbol=:heap, verbose=false)
    _bury_check_backend(backend)
    target = _bury_check_target(g, r0)
    selected_seed = _bury_check_seed(g, seed)
    backend === :legacy && return _bury_seeded_legacy(g, target, selected_seed; verbose=verbose)
    colored = _bury_color_kernel(backend)(g, target; seed=selected_seed, verbose=verbose)
    return _bury_registry(colored), colored
end

# Shared partition validation and vertex mapping.
"""Validate the common nonempty, exactly balanced partitioning boundary."""
function _check_partition_count(g, k)
    k isa Integer && 1 <= k <= nv(g) ||
        throw(ArgumentError("k must be an integer satisfying 1 <= k <= nv(g)"))
    nv(g) % k == 0 ||
        throw(ArgumentError("k must divide the number of vertices"))
    return Int(k)
end

"""Helper function to some of the other methods"""
function inverse_map(map, size=maximum(map))
    inv_map = [0 for _ in 1:size]
    for (index, e) in enumerate(map)
        inv_map[e] = index
    end
    return inv_map
end

# Graph-first partitioning entry points.
"""
    partition_graph(g, k; algorithm=:bury, backend=:heap, verbose=false) -> GraphPartition

Partition a nonempty simple undirected graph into `k` equally sized parts.
Require `1 <= k <= nv(g)` and `nv(g) % k == 0`.

Algorithms delegate to the existing implementations:
- `:bury` calls [`bury_k`](@ref).
- `:bury_seeded` calls [`bury_k_seeded`](@ref), with its deterministic
  minimum-weight seeder for each part.
- `:random` calls [`k_partition_random`](@ref), using the default Random RNG.
- `:saran_vazirani` calls [`k_partition_saran_vazirani`](@ref); its trace is
  available through that legacy function instead.

BURY accepts `backend=:heap` (default), `:scan`, or `:legacy`, forwarded through
every partition round. Seed and partition policies are unchanged. `verbose`
is forwarded to BURY. Other algorithms only accept the default backend value,
which does not apply to them, and retain their existing output behavior.
An incomplete registry raises an error. The Saran–Vazirani baseline first
tries to pack the original components and reports failure explicitly.
"""
function partition_graph(g::Graphs.Graph, k;
                         algorithm::Symbol=:bury, backend::Symbol=:heap,
                         verbose::Bool=false)
    k = _check_partition_count(g, k)
    _bury_check_backend(backend)
    if algorithm in (:random, :saran_vazirani) && backend !== :heap
        throw(ArgumentError("backend selection applies only to :bury and :bury_seeded"))
    end

    reg = if algorithm === :bury
        bury_k(g, k; backend=backend, verbose=verbose)
    elseif algorithm === :bury_seeded
        bury_k_seeded(g, k; backend=backend, verbose=verbose)
    elseif algorithm === :random
        k_partition_random(g, k)
    elseif algorithm === :saran_vazirani
        first(k_partition_saran_vazirani(g, k))
    else
        throw(ArgumentError("unknown algorithm $algorithm; use :bury, :bury_seeded, :random, or :saran_vazirani"))
    end

    length(reg) == nv(g) ||
        error("$algorithm returned an incomplete registry")
    return GraphPartition(reg; nparts=k)
end

function _bipartition_targets(g, targets)
    n = nv(g)
    n >= 2 || throw(ArgumentError("bipartition_graph requires at least two vertices"))
    targets isa AbstractRange{<:Integer} && !isempty(targets) && step(targets) == 1 ||
        throw(ArgumentError("targets must be a nonempty contiguous ascending integer range"))
    1 <= first(targets) <= last(targets) < n ||
        throw(ArgumentError("targets must lie in 1:nv(g)-1 for two nonempty parts"))
    return Int(first(targets)):Int(last(targets))
end

function _bipartition_cost(objective)
    objective === matching_cost && return _matching_cost
    objective === edge_cut_cost && return _edge_cut_cost
    objective === pairwise_cut_rank && return _bury_cut_rank
    # Arbitrary callbacks must run once per target, with their public contract.
    return (colors, g) -> objective(g, GraphPartition(Int.(colors) .+ 1; nparts=2))
end

"""
    bipartition_graph(g; targets=fld(nv(g),2):cld(nv(g),2),
                     objective=matching_cost, backend=:heap) -> GraphPartition

Run ordinary BURY independently for each target and minimize the scalar
`objective(g, p)`. Targets specify the size of label 2, are evaluated in
ascending order, and ties retain the last target. The default evaluates one
target for even graphs and both nearest-balanced targets for odd graphs.

Require at least two vertices and a nonempty contiguous ascending integer
range within `1:nv(g)-1`, validated before invoking the objective. Return two
nonempty declared parts without modifying the graph. Built-in matching,
edge-cut, and pairwise-rank objectives use the scalar membership kernels.
All BURY backends (`:heap`, `:scan`, `:legacy`) are supported.

This helper uses ordinary BURY; seeded/custom partitioners remain available
through [`bury_not_balanced_bipartition`](@ref). [`partition_graph`](@ref)
continues to require exactly equal-sized parts and divisibility.
"""
function bipartition_graph(g::Graphs.Graph;
                           targets=fld(nv(g), 2):cld(nv(g), 2),
                           objective=matching_cost, backend::Symbol=:heap)
    _bury_check_backend(backend)
    targets = _bipartition_targets(g, targets)
    colors = _bury_select_colors(g, targets, backend, _bipartition_cost(objective))
    return GraphPartition(Int.(colors) .+ 1; nparts=2)
end
