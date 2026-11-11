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
