"""
    GraphPartition(labels; nparts=isempty(labels) ? 1 : maximum(labels))
    GraphPartition(registry::AbstractDict; nparts=...)
    GraphPartition(land::quantumLand)

Membership of vertices `1:n` in parts `1:nparts`. `labels[v]` is the part
containing vertex `v`; a dictionary must contain exactly the keys `1:n`.
Labels and `nparts` must be integers, and `nparts` must be positive. Empty
parts are allowed, including an empty graph with a declared partition count.

Construction copies the input. Use [`membership`](@ref), [`registry`](@ref),
and [`partition_sizes`](@ref) to inspect a partition; treat its fields as
internal and read-only. To change membership, construct a new partition.
The accessors return independent collections.

# Examples
```julia
p = GraphPartition([1, 1, 2, 2])
membership(p)       # [1, 1, 2, 2]
partition_sizes(p)  # [2, 2]
```
"""
struct GraphPartition
    labels::Vector{Int}
    nparts::Int

    function GraphPartition(labels::AbstractVector{<:Integer}, nparts::Integer)
        nparts >= 1 || throw(ArgumentError("nparts must be positive"))
        all(label -> 1 <= label <= nparts, labels) ||
            throw(ArgumentError("partition labels must lie in 1:nparts"))
        return new(collect(Int, labels), Int(nparts))
    end
end

GraphPartition(labels::AbstractVector{<:Integer};
               nparts::Integer=isempty(labels) ? 1 : maximum(labels)) =
    GraphPartition(labels, nparts)

function GraphPartition(reg::AbstractDict; nparts::Union{Nothing,Integer}=nothing)
    n = length(reg)
    all(key -> key isa Integer && 1 <= key <= n, keys(reg)) &&
        all(v -> haskey(reg, v), 1:n) ||
        throw(ArgumentError("registry must contain exactly the vertex keys 1:n"))
    all(label -> label isa Integer, values(reg)) ||
        throw(ArgumentError("partition labels must be integers"))
    parts = isnothing(nparts) ? (isempty(reg) ? 1 : maximum(values(reg))) : nparts
    return GraphPartition(Int[reg[v] for v in 1:n], parts)
end

function GraphPartition(land::quantumLand)
    length(land.registry) == land.numVillagers ||
        throw(ArgumentError("land registry does not match its vertex count"))
    return GraphPartition(land.registry; nparts=land.numHamlets)
end

"""
    membership(p::GraphPartition)

Return a copy of the part labels in vertex order. Modifying it does not change `p`.
"""
membership(p::GraphPartition) = copy(p.labels)

"""
    registry(p::GraphPartition)

Return an independent `Dict{Int,Int}` mapping each vertex to its part, for use
with the original QuantumHamlets API.
"""
registry(p::GraphPartition) = Dict{Int,Int}(v => label for (v, label) in enumerate(p.labels))

"""
    partition_sizes(p::GraphPartition)

Return the number of vertices in each of the `p.nparts` parts, including zeros
for empty parts.
"""
function partition_sizes(p::GraphPartition)
    sizes = zeros(Int, p.nparts)
    for label in p.labels
        sizes[label] += 1
    end
    return sizes
end

Base.length(p::GraphPartition) = length(p.labels)

"""
    quantumLand(p::GraphPartition)

Create an independent legacy `quantumLand` from `p`. Existing constructors
and their behavior are unchanged.
"""
quantumLand(p::GraphPartition) = quantumLand(registry(p), p.nparts, length(p))
