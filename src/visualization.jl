######################## Visualization functions ######################## 
function print_graph(g::Graph; edgecolors=nothing, nodecolors=[:white for _ in 1:nv(g)], nodelabels = true)
    labels = collect(1:nv(g))
    if nodelabels == false
        labels = nothing
    end
    thicken_lines = edgecolors.==:red
    f, ax, p = graphplot(g,
                        node_color=nodecolors,
                        edge_color = edgecolors,
                        edge_width= [2 for _ in 1:ne(g)] .+ 2*thicken_lines,
                        nlabels_align=(:center,:center);
                        ilabels= labels, 
                        size = (900,700))
    hidedecorations!(ax); hidespines!(ax); ax.aspect = DataAspect()

    return f
end

"""Returns a graph with nodes colored to match the hamlet they belong to. Edges are colored red
to indicate cost. Cost is default any cross hamlet edge. Using a different method can color edges differently.
Use `naive_cost_of_partition` or `matching_cost_of_partition` for edge highlights.
`cut_rank_cost` supplies no edge highlights. The vertex-cover method returns two
costs rather than edge colors and is not compatible with this interface.
The simpler `plot_partition(g, p; highlight=...)` API selects these modes explicitly."""
function visualize_graph_on_land(land::quantumLand, g::Graphs.Graph; method=QuantumHamlets.naive_cost_of_partition)
    cost, edgecolors = method(land,g)
    hamlet_colors = Makie.wong_colors()

    nodecolors = [hamlet_colors[(land.registry[i]%7)+1] for i in 1:land.numVillagers]
    f = print_graph(g, edgecolors=edgecolors, nodecolors=nodecolors)
    return f, cost
end

# GraphPartition plotting adapter.
function _land_for(g::Graphs.Graph, p::GraphPartition)
    _labels_for(g, p)
    return quantumLand(p)
end

"""
    plot_partition(g, p::GraphPartition; highlight=:matching)

Return `(figure, cost)` using the existing plotting implementation.
`highlight=:matching` highlights pairwise matching edges; `:cut_edges`
highlights every cross-part edge. With `:none`, show membership without
highlighted edges and return `nothing` for the cost.

Plotting dependencies still load with the package. To display rank, evaluate
[`pairwise_cut_rank`](@ref) separately and use `highlight=:none`.
"""
function plot_partition(g::Graphs.Graph, p::GraphPartition; highlight::Symbol=:matching)
    land = _land_for(g, p)
    if highlight === :matching
        return visualize_graph_on_land(land, g; method=matching_cost_of_partition)
    elseif highlight === :cut_edges
        return visualize_graph_on_land(land, g; method=naive_cost_of_partition)
    elseif highlight === :none
        return visualize_graph_on_land(land, g; method=(_, _) -> (nothing, nothing))
    else
        throw(ArgumentError("unknown highlight $highlight; use :matching, :cut_edges, or :none"))
    end
end
