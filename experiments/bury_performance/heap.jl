module BuryHeap

using Graphs

export bury_heap

# Standalone performance experiment. Nothing in this module replaces or adds
# methods to QuantumHamlets. Supported inputs are simple undirected graphs,
# integer targets 0:nv(g), and (for the seeded path) an explicit valid seed.

struct IndexedMinHeap
    vertices::Vector{Int}
    positions::Vector{Int}
    weights::Vector{Int}
end

@inline function _less(h::IndexedMinHeap, u::Int, v::Int)
    wu, wv = h.weights[u], h.weights[v]
    return wu < wv || (wu == wv && u < v)
end

@inline function _swap!(h::IndexedMinHeap, i::Int, j::Int)
    u, v = h.vertices[i], h.vertices[j]
    h.vertices[i], h.vertices[j] = v, u
    h.positions[u], h.positions[v] = j, i
    return nothing
end

function _up!(h::IndexedMinHeap, i::Int)
    while i > 1
        parent = i ÷ 2
        _less(h, h.vertices[i], h.vertices[parent]) || break
        _swap!(h, i, parent)
        i = parent
    end
    return nothing
end

function _down!(h::IndexedMinHeap, i::Int)
    n = length(h.vertices)
    while 2i <= n
        child = 2i
        if child < n && _less(h, h.vertices[child + 1], h.vertices[child])
            child += 1
        end
        _less(h, h.vertices[child], h.vertices[i]) || break
        _swap!(h, i, child)
        i = child
    end
    return nothing
end

function _heap(weights::Vector{Int})
    n = length(weights)
    h = IndexedMinHeap(collect(1:n), collect(1:n), weights)
    for i in n÷2:-1:1
        _down!(h, i)
    end
    return h
end

function _empty_heap(weights::Vector{Int})
    vertices = Int[]
    sizehint!(vertices, length(weights))
    return IndexedMinHeap(vertices, zeros(Int, length(weights)), weights)
end

function _insert!(h::IndexedMinHeap, v::Int)
    h.positions[v] == 0 || return nothing
    push!(h.vertices, v)
    h.positions[v] = length(h.vertices)
    _up!(h, length(h.vertices))
    return nothing
end

function _remove!(h::IndexedMinHeap, v::Int)
    i = h.positions[v]
    i == 0 && return nothing
    last = pop!(h.vertices)
    h.positions[v] = 0
    if i <= length(h.vertices)
        h.vertices[i] = last
        h.positions[last] = i
        if i > 1 && _less(h, last, h.vertices[i ÷ 2])
            _up!(h, i)
        else
            _down!(h, i)
        end
    end
    return nothing
end

@inline function _decrease!(h::IndexedMinHeap, v::Int)
    i = h.positions[v]
    i == 0 || _up!(h, i)
    return nothing
end

@inline function _decrement!(global_heap::IndexedMinHeap, ::Nothing, v::Int)
    # Processed vertices have infinite weight in the original implementation.
    # Removing them from the heap is equivalent; do not decrement them.
    global_heap.positions[v] == 0 && return nothing
    global_heap.weights[v] -= 1
    _decrease!(global_heap, v)
    return nothing
end

@inline function _decrement!(global_heap::IndexedMinHeap,
                            region_heap::IndexedMinHeap, v::Int)
    global_heap.positions[v] == 0 && return nothing
    global_heap.weights[v] -= 1
    _decrease!(global_heap, v)
    _decrease!(region_heap, v)
    return nothing
end

@inline _include_region!(::IndexedMinHeap, ::Nothing, ::Int) = nothing

@inline function _include_region!(global_heap::IndexedMinHeap,
                                 region_heap::IndexedMinHeap, v::Int)
    global_heap.positions[v] == 0 || _insert!(region_heap, v)
    return nothing
end

function _visit!(g, colored, global_heap, region_heap, v)
    newly_colored = !colored[v]
    if newly_colored
        _decrement!(global_heap, region_heap, v)
        for u in neighbors(g, v)
            _decrement!(global_heap, region_heap, Int(u))
        end
        colored[v] = true
    end
    # Region grows for every member of the selected neighborhood, even when
    # that vertex was already colored. Colored candidates are NOT removed.
    _include_region!(global_heap, region_heap, v)
    return Int(newly_colored)
end

@inline _selection_heap(global_heap, ::Nothing) = global_heap

@inline function _selection_heap(global_heap, region_heap::IndexedMinHeap)
    return isempty(region_heap.vertices) ? global_heap : region_heap
end

@inline _remove_region!(::Nothing, ::Int) = nothing
@inline _remove_region!(region_heap::IndexedMinHeap, v::Int) = _remove!(region_heap, v)

function _run(g, r0, global_heap, region_heap)
    colored = fill(false, nv(g))
    colored_count = 0
    remaining = r0

    while colored_count < r0
        selection = _selection_heap(global_heap, region_heap)
        choice = first(selection.vertices)
        cost = global_heap.weights[choice]
        take_neighborhood = remaining >= cost
        remaining -= take_neighborhood ? cost : 1

        # The choice cannot be selected again. Removing it before the updates
        # changes no remaining candidate's weight and avoids repairing its key.
        _remove!(global_heap, choice)
        _remove_region!(region_heap, choice)

        colored_count += _visit!(g, colored, global_heap, region_heap, choice)
        if take_neighborhood
            for v in neighbors(g, choice)
                colored_count += _visit!(g, colored, global_heap, region_heap, Int(v))
            end
        end
    end
    return colored
end

"""
    _colors(g, r0=nv(g)÷2; seed=nothing)

Experimental color-only BURY kernel with indexed minimum heaps. `seed=nothing`
uses ordinary BURY; an explicit seed uses the original seeded region rule.
Weight ties select the smallest vertex ID. Colored-but-unprocessed vertices
remain candidates, as in the original implementation.
"""
function _colors(g::Graphs.Graph, r0::Integer=nv(g)÷2;
                 seed::Union{Nothing,Int}=nothing)
    n = nv(g)
    0 <= r0 <= n || throw(ArgumentError("target must lie in 0:nv(g)"))
    isnothing(seed) || 1 <= seed <= n || throw(ArgumentError("invalid seed"))
    weights = [degree(g, v) + 1 for v in vertices(g)]
    global_heap = _heap(weights)
    if isnothing(seed)
        return _run(g, r0, global_heap, nothing)
    end

    region_heap = _empty_heap(weights)
    _insert!(region_heap, seed)
    for v in neighbors(g, seed)
        _insert!(region_heap, Int(v))
    end
    return _run(g, r0, global_heap, region_heap)
end

"""
    bury_heap(g, r0=nv(g)÷2; seed=nothing)

Return `(registry, colored)` using the isolated heap prototype. Registry values
and colors follow the original BURY convention (`false => 1`, `true => 2`).
This experiment does not implement verbose logging or automatic seed search.
"""
function bury_heap(g::Graphs.Graph, r0::Integer=nv(g)÷2;
                   seed::Union{Nothing,Int}=nothing)
    colored = _colors(g, r0; seed=seed)
    registry = Dict{Int,Int}()
    sizehint!(registry, nv(g))
    for v in vertices(g)
        registry[v] = Int(colored[v]) + 1
    end
    return registry, colored
end

end # module BuryHeap
