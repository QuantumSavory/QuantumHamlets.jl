# An indexed binary min-heap orders candidates by weight, then vertex ID.
# The two heaps in seeded BURY share weights but maintain separate positions.
struct _BuryMinHeap
    vertices::Vector{Int}
    positions::Vector{Int}
    weights::Vector{Int}
end

@inline function _bury_heap_less(h::_BuryMinHeap, u::Int, v::Int)
    wu, wv = h.weights[u], h.weights[v]
    return wu < wv || (wu == wv && u < v)
end

@inline function _bury_heap_swap!(h::_BuryMinHeap, i::Int, j::Int)
    u, v = h.vertices[i], h.vertices[j]
    h.vertices[i], h.vertices[j] = v, u
    h.positions[u], h.positions[v] = j, i
    return nothing
end

function _bury_heap_up!(h::_BuryMinHeap, i::Int)
    while i > 1
        parent = i ÷ 2
        _bury_heap_less(h, h.vertices[i], h.vertices[parent]) || break
        _bury_heap_swap!(h, i, parent)
        i = parent
    end
    return nothing
end

function _bury_heap_down!(h::_BuryMinHeap, i::Int)
    n = length(h.vertices)
    while 2i <= n
        child = 2i
        if child < n && _bury_heap_less(h, h.vertices[child + 1], h.vertices[child])
            child += 1
        end
        _bury_heap_less(h, h.vertices[child], h.vertices[i]) || break
        _bury_heap_swap!(h, i, child)
        i = child
    end
    return nothing
end

function _bury_heap(weights::Vector{Int})
    n = length(weights)
    h = _BuryMinHeap(collect(1:n), collect(1:n), weights)
    for i in n÷2:-1:1
        _bury_heap_down!(h, i)
    end
    return h
end

function _bury_heap_empty(weights::Vector{Int})
    vertices = Int[]
    sizehint!(vertices, length(weights))
    return _BuryMinHeap(vertices, zeros(Int, length(weights)), weights)
end

function _bury_heap_insert!(h::_BuryMinHeap, v::Int)
    h.positions[v] == 0 || return nothing
    push!(h.vertices, v)
    h.positions[v] = length(h.vertices)
    _bury_heap_up!(h, length(h.vertices))
    return nothing
end

function _bury_heap_remove!(h::_BuryMinHeap, v::Int)
    i = h.positions[v]
    i == 0 && return nothing
    last = pop!(h.vertices)
    h.positions[v] = 0
    if i <= length(h.vertices)
        h.vertices[i] = last
        h.positions[last] = i
        if i > 1 && _bury_heap_less(h, last, h.vertices[i ÷ 2])
            _bury_heap_up!(h, i)
        else
            _bury_heap_down!(h, i)
        end
    end
    return nothing
end

@inline function _bury_heap_decrease!(h::_BuryMinHeap, v::Int)
    i = h.positions[v]
    i == 0 || _bury_heap_up!(h, i)
    return nothing
end

@inline function _bury_heap_decrement!(global_heap::_BuryMinHeap, ::Nothing, v::Int)
    # Processed vertices have Inf weight in the original implementation.
    # Removing them is equivalent; their weight no longer needs updating.
    global_heap.positions[v] == 0 && return nothing
    global_heap.weights[v] -= 1
    _bury_heap_decrease!(global_heap, v)
    return nothing
end

@inline function _bury_heap_decrement!(global_heap::_BuryMinHeap,
                                     region_heap::_BuryMinHeap, v::Int)
    global_heap.positions[v] == 0 && return nothing
    global_heap.weights[v] -= 1
    _bury_heap_decrease!(global_heap, v)
    _bury_heap_decrease!(region_heap, v)
    return nothing
end

@inline _bury_heap_include_region!(::_BuryMinHeap, ::Nothing, ::Int) = nothing

@inline function _bury_heap_include_region!(global_heap::_BuryMinHeap,
                                          region_heap::_BuryMinHeap, v::Int)
    global_heap.positions[v] == 0 || _bury_heap_insert!(region_heap, v)
    return nothing
end

function _bury_heap_visit!(g, colored, global_heap, region_heap, v)
    newly_colored = !colored[v]
    if newly_colored
        _bury_heap_decrement!(global_heap, region_heap, v)
        for u in neighbors(g, v)
            _bury_heap_decrement!(global_heap, region_heap, Int(u))
        end
        colored[v] = true
    end
    # Visit every member of the selected neighborhood, even if already colored.
    # Colored neighbors remain candidates until they are chosen themselves.
    _bury_heap_include_region!(global_heap, region_heap, v)
    return Int(newly_colored)
end

@inline _bury_heap_selection(global_heap, ::Nothing) = global_heap

@inline function _bury_heap_selection(global_heap, region_heap::_BuryMinHeap)
    return isempty(region_heap.vertices) ? global_heap : region_heap
end

@inline _bury_heap_remove_region!(::Nothing, ::Int) = nothing
@inline _bury_heap_remove_region!(region_heap::_BuryMinHeap, v::Int) =
    _bury_heap_remove!(region_heap, v)

# Construct the legacy display representation only when tracing is requested.
function _bury_heap_trace_weights(h::_BuryMinHeap)
    return Float64[h.positions[v] == 0 ? Inf : h.weights[v] for v in eachindex(h.weights)]
end

function _bury_heap_run(g, r0, global_heap, region_heap, verbose)
    colored = fill(false, nv(g))
    colored_count = 0
    remaining = r0
    display_remaining = r0
    verbose && _bury_trace_initial(_bury_heap_trace_weights(global_heap), colored)

    while colored_count < r0
        selection = _bury_heap_selection(global_heap, region_heap)
        choice = first(selection.vertices)
        cost = global_heap.weights[choice]
        verbose && _bury_trace_choice(display_remaining, choice, Float64(cost))
        take_neighborhood = remaining >= cost
        remaining -= take_neighborhood ? cost : 1
        if verbose
            display_remaining -= take_neighborhood ? Float64(cost) : 1
        end

        # Remove only the chosen vertex. Doing so before updates changes no
        # remaining candidate's weight and avoids repairing an irrelevant key.
        _bury_heap_remove!(global_heap, choice)
        _bury_heap_remove_region!(region_heap, choice)

        colored_count += _bury_heap_visit!(g, colored, global_heap, region_heap, choice)
        if take_neighborhood
            for v in neighbors(g, choice)
                colored_count += _bury_heap_visit!(g, colored, global_heap, region_heap, Int(v))
            end
        end
        verbose && _bury_trace_updated(_bury_heap_trace_weights(global_heap), colored)
    end
    return colored
end

function _bury_heap_colors(g::Graphs.Graph, r0::Integer;
                           seed::Union{Nothing,Integer}=nothing, verbose=false)
    n = nv(g)
    0 <= r0 <= n || throw(ArgumentError("target must lie in 0:nv(g)"))
    isnothing(seed) || 1 <= seed <= n || throw(ArgumentError("invalid seed"))
    weights = [Graphs.degree(g, v) + 1 for v in vertices(g)]
    global_heap = _bury_heap(weights)
    if isnothing(seed)
        return _bury_heap_run(g, r0, global_heap, nothing, verbose)
    end

    region_heap = _bury_heap_empty(weights)
    _bury_heap_insert!(region_heap, Int(seed))
    for v in neighbors(g, Int(seed))
        _bury_heap_insert!(region_heap, Int(v))
    end
    return _bury_heap_run(g, r0, global_heap, region_heap, verbose)
end
