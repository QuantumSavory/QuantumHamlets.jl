"""Randomly partitions the graph into k pieces\n
Input: a nonempty graph from Graphs.jl, and integer k in 1:nv(g) dividing nv(g).\n
Output: a registry (dictionary mapping villager ids to hamlet ids)"""
function k_partition_random(g::Graphs.Graph, k)
    k = _check_partition_count(g, k)
    nodes = shuffle(collect(1:nv(g)))
    new_reg = Dict()
    cap = nv(g) ÷ k
    for i in 1:k
        for j in 1:cap
            new_reg[nodes[(i-1)*cap+j]] = i
        end
    end
    return new_reg
end

# Load the bundled Python package on first use, without modifying sys.path.
const _MIN_CUTRANK = PythonCall.pynew()
const _SA_LOCK = ReentrantLock()
const _METIS_LOCK = ReentrantLock()

function _min_cutrank_module()
    if PythonCall.pyisnull(_MIN_CUTRANK)
        util = PythonCall.pyimport("importlib.util")
        path = joinpath(@__DIR__, "min_cutrank", "__init__.py")
        spec = util.spec_from_file_location("_quantumhamlets_min_cutrank", path)
        mod = util.module_from_spec(spec)
        PythonCall.pyimport("sys").modules[spec.name] = mod
        spec.loader.exec_module(mod)
        PythonCall.pycopy!(_MIN_CUTRANK, mod)
    end
    return _MIN_CUTRANK
end

function _baseline_seed(seed)
    seed isa Integer && !(seed isa Bool) && seed >= 0 ||
        throw(ArgumentError("seed must be a nonnegative integer"))
    return seed
end

function _sa_graph(g)
    nv(g) >= 2 || throw(ArgumentError("annealing requires at least two vertices"))
    any(v -> has_edge(g, v, v), vertices(g)) &&
        throw(ArgumentError("annealing requires a graph without self loops"))
    return nothing
end

"""
    min_cutrank(g; initial=nothing, temperatures=range(1.0, 0.1; length=10),
                seed=0, keep_best=true) -> GraphPartition

Run the bundled min-cutrank simulated annealer to minimize GF(2) cut rank.
By default start from a random floor/ceiling bipartition. Supply `initial` as a
`GraphPartition`, a registry, a vector of labels 1/2, or Boolean row flags.
Both parts must be nonempty; their sizes are preserved, even if unbalanced.
The smaller side is always used as the annealer's row side.

`temperatures` supplies one positive finite temperature per sweep; an empty
sequence returns the initializer without annealing. The default is upstream's
ten-sweep example schedule. A local RNG controlled by `seed` makes calls
repeatable without resetting Julia's or Python's global RNG.
Call from the main Julia task unless managing PythonCall's GIL explicitly.

By default return the first best rank seen, including the initializer.
`keep_best=false` returns upstream's final state, which can have a higher rank.
There is no time limit or restart policy. Inputs are not modified. Source and
MIT attribution are in `src/min_cutrank/`; PythonCall manages Python 3.12 and
the annealer needs no additional Python packages.
"""
function min_cutrank(g::Graphs.Graph; initial=nothing,
                     temperatures=range(1.0, 0.1; length=10), seed=0,
                     keep_best::Bool=true)
    _sa_graph(g)
    _baseline_seed(seed)
    temps = collect(temperatures)
    all(t -> t isa Real && isfinite(t) && t > 0, temps) ||
        throw(ArgumentError("temperatures must contain positive finite numbers"))
    temps = Float64.(temps)
    all(t -> isfinite(t) && t > 0, temps) ||
        throw(ArgumentError("temperatures must be positive finite Float64 values"))
    n = nv(g)
    p = if initial === nothing
        labels = fill(2, n)
        labels[randperm(Xoshiro(seed), n)[1:fld(n, 2)]] .= 1
        GraphPartition(labels; nparts=2)
    elseif initial isa GraphPartition
        initial
    elseif initial isa AbstractVector{Bool}
        GraphPartition(ifelse.(initial, 1, 2); nparts=2)
    elseif initial isa Union{AbstractVector{<:Integer},AbstractDict}
        GraphPartition(initial; nparts=2)
    else
        throw(ArgumentError("initial must be a GraphPartition, registry, labels, or Boolean flags"))
    end
    length(p) == n && p.nparts == 2 && all(>(0), partition_sizes(p)) ||
        throw(ArgumentError("initial must assign every vertex to two nonempty parts"))
    flags = membership(p) .== 1
    count(flags) > fld(n, 2) && (flags .= .!flags)
    isempty(temps) && return GraphPartition(ifelse.(flags, 1, 2); nparts=2)
    PythonCall.C.PyGILState_Check() == 1 ||
        throw(ArgumentError("call min_cutrank from Julia's main task, or acquire PythonCall's GIL"))
    adjacency = [Int[has_edge(g, u, v) for v in vertices(g)] for u in vertices(g)]
    result = lock(_SA_LOCK) do
        PythonCall.pyconvert(Vector{Bool},
            _min_cutrank_module().anneal(adjacency, flags, temps, seed, keep_best))
    end
    length(result) == n && count(result) == count(flags) ||
        error("min-cutrank returned invalid partition capacities")
    return GraphPartition(ifelse.(result, 1, 2); nparts=2)
end

"""
    bury_sa(g; backend=:heap, kwargs...) -> GraphPartition

Initialize [`min_cutrank`](@ref) with ordinary near-balanced BURY selected by
rank. Forward `temperatures`, `seed`, and `keep_best` to the annealer.
"""
function bury_sa(g::Graphs.Graph; backend::Symbol=:heap, kwargs...)
    p = bipartition_graph(g; objective=pairwise_cut_rank, backend=backend)
    return min_cutrank(g; initial=p, kwargs...)
end

function _sqrt_bury_bipartition(g; seed=0, backend::Symbol=:heap)
    _sa_graph(g)
    _baseline_seed(seed)
    rng = Xoshiro(seed)
    seeded = (graph, target) -> bury_seeded(graph, target; backend=backend,
                                           rng=rng, log_improvements=false)
    reg, _ = bury_not_balanced_bipartition(g, fld(nv(g), 2), cld(nv(g), 2);
        bury_v=seeded, cost_f=cut_rank_cost, backend=backend)
    return GraphPartition(reg; nparts=2)
end

"""
    bury_seeded_sa(g; seed=0, backend=:heap, kwargs...) -> GraphPartition

Initialize [`min_cutrank`](@ref) with sqrt-seeded BURY: ceil(sqrt(n)) sampled
seeds with replacement, inner matching selection, outer floor/ceiling rank
selection. Uses local RNGs. Set `temperatures=[]` to obtain just this initializer.
"""
function bury_seeded_sa(g::Graphs.Graph; seed=0, backend::Symbol=:heap, kwargs...)
    p = _sqrt_bury_bipartition(g; seed=seed, backend=backend)
    return min_cutrank(g; initial=p, seed=seed, kwargs...)
end

# Capacity repair for explicitly balanced METIS comparisons. Give any extra
# capacities to the largest existing parts to avoid unnecessary vertex moves.
function _balance_metis(g, labels, k)
    sizes = [count(==(i), labels) for i in 1:k]
    targets = fill(nv(g) ÷ k, k)
    order = sortperm(1:k; by=i -> (-sizes[i], i))
    targets[order[1:rem(nv(g), k)]] .+= 1
    while sizes != targets
        best = (typemax(Int), typemax(Int), typemax(Int))
        for v in vertices(g)
            source = labels[v]
            sizes[source] > targets[source] || continue
            counts = [count(u -> labels[u] == i, neighbors(g, v)) for i in 1:k]
            for dest in 1:k
                sizes[dest] < targets[dest] || continue
                candidate = (counts[source] - counts[dest], Int(v), dest)
                candidate < best && (best = candidate)
            end
        end
        _, v, dest = best
        sizes[labels[v]] -= 1
        labels[v] = dest
        sizes[dest] += 1
    end
    return labels
end

"""
    metis_partition(g, k=2; seed=0, balance=false) -> GraphPartition

Run METIS KWAY edge-cut partitioning with equal target weights, native seed,
and its standard 3% imbalance tolerance. Require `1 <= k <= nv(g)`. Native
outputs can be unbalanced or contain empty parts; inspect `partition_sizes`.

With `balance=true`, explicitly repair the result to floor/ceiling capacities:
move vertices from oversized to undersized parts, choosing the smallest edge-cut
increase, then smallest vertex ID and destination on ties. This postprocessing
can change the native cost. The extra capacities go to the largest native parts
(smallest label on ties). METIS's global options are restored after each call.
"""
function metis_partition(g::Graphs.Graph, k=2; seed=0, balance::Bool=false)
    n = nv(g)
    k isa Integer && !(k isa Bool) && 1 <= k <= n ||
        throw(ArgumentError("k must be an integer in 1:nv(g)"))
    _baseline_seed(seed)
    seed <= typemax(Cint) || throw(ArgumentError("METIS seed must fit in Int32"))
    labels = lock(_METIS_LOCK) do
        saved = copy(Metis.options)
        try
            fill!(Metis.options, -1)
            Metis.options[Int(Metis.LibMetis.METIS_OPTION_NUMBERING) + 1] = 1
            Metis.options[Int(Metis.LibMetis.METIS_OPTION_SEED) + 1] = seed
            Metis.options[Int(Metis.LibMetis.METIS_OPTION_UFACTOR) + 1] = 30
            Int.(Metis.partition(g, k; alg=:KWAY))
        finally
            copyto!(Metis.options, saved)
        end
    end
    p = GraphPartition(labels; nparts=k)
    length(p) == n || error("METIS returned an incomplete partition")
    balance && _balance_metis(g, labels, k)
    return GraphPartition(labels; nparts=k)
end

"""
    metis_sa(g; seed=0, balance=true, kwargs...) -> GraphPartition

Initialize [`min_cutrank`](@ref) with METIS. By default explicitly repair METIS
to near-balanced capacities before annealing, as in [`metis_partition`](@ref).
With `balance=false`, preserve native sizes; an empty native part is rejected.
"""
function metis_sa(g::Graphs.Graph; seed=0, balance::Bool=true, kwargs...)
    p = metis_partition(g, 2; seed=seed, balance=balance)
    return min_cutrank(g; initial=p, seed=seed, kwargs...)
end

"""Implementation of the n/2 balanced partitioning algorithm from https://ics.uci.edu/~vazirani/k-cut.pdf \n
Input: a graph and number of pices to partition into.\n
Output: `(registry, trace)`, with component snapshots and balance checks in the trace.
The first snapshot checks the original connected components before removing edges.
Require a nonempty graph and integer `k` in `1:nv(g)` dividing `nv(g)`.
Solver failure or failure to find a complete balanced assignment raises an error.
"""
function k_partition_saran_vazirani(g::Graphs.Graph, k)
    k = _check_partition_count(g, k)
    dump = []
    function pack_components(graph)
        components = Graphs.connected_components(graph)
        sizes = length.(components)
        inds = try
            partition(sizes, k; optimizer=HiGHS.Optimizer, strict=true)
        catch err
            err isa InterruptException && rethrow()
            error("Saran–Vazirani component packing failed: $(sprint(showerror, err))")
        end
        length(inds) == length(components) && all(i -> i isa Integer && 1 <= i <= k, inds) ||
            error("Saran–Vazirani component packing returned invalid labels")
        balanced = all(i -> sum(sizes[inds .== i]) == nv(g) ÷ k, 1:k)
        push!(dump, (copy(components), balanced))
        balanced || return nothing

        new_reg = Dict()
        for (component, label) in zip(components, inds), node in component
            new_reg[node] = label
        end
        return new_reg
    end

    # A disconnected input may already admit a balanced zero-cut partition.
    initial = pack_components(g)
    initial === nothing || return initial, dump

    # Given a graph g, and two vectors corresponding the vertices within that graph, returns the set of edges the crosses the partition
    function edges_spanning_partition(g::Graphs.Graph, v1, v2)
        if length(v1) > length(v2)
            v3 = v2
            v2 = v1
            v1 = v3
        end

        edges = []
        for v in v1
            neighbors = Graphs.neighbors(g,v)
            for n in neighbors
                if n ∉ v1
                    push!(edges,Edge(v,n))
                end
            end
        end

        return edges
    end

    d = Graphs.DiGraph(g)
    capacity_matrix = zeros(Int, nv(g), nv(g))
    for e in edges(d)
        u, v = src(e), dst(e)
        capacity_matrix[u, v] = 1
    end

    min_cuts = []
    for e in edges(g)
        cut = GraphsFlows.mincut(d,src(e),dst(e),capacity_matrix,GraphsFlows.BoykovKolmogorovAlgorithm())
        push!(min_cuts, cut)
    end

    unique!(min_cuts)
    sort!(min_cuts, by = x-> x[3], rev=false)

    edge_set = []
    for c in min_cuts
        push!(edge_set, edges_spanning_partition(g,c[1],c[2]))
    end

    g_prime = copy(g)
    for cut in edge_set
        for e in cut
            rem_edge!(g_prime, e)
        end
    
        new_reg = pack_components(g_prime)
        new_reg === nothing || return new_reg, dump
    end

    error("Saran–Vazirani did not find a complete balanced partition")
end
