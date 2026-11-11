using QuantumHamlets, Graphs, Statistics
const Q = QuantumHamlets
using QuantumHamlets.CairoMakie

function probe(f, name)
    println("START\t", name); flush(stdout)
    try
        value = f()
        println("RESULT\t", name, "\t", repr(value))
    catch ex
        println("ERROR\t", name, "\t", typeof(ex), "\t", sprint(showerror, ex))
    end
    println("END\t", name); flush(stdout)
end

println("READY"); flush(stdout)

probe("saran_edgeless4_k2") do
    reg, trace = Q.k_partition_saran_vazirani(SimpleGraph(4), 2)
    (registry=reg, trace_length=length(trace))
end

disconnected = SimpleGraph(4)
add_edge!(disconnected,1,2);add_edge!(disconnected,3,4)
cases = (("path4",path_graph(4)),("cycle4",cycle_graph(4)),
         ("complete4",complete_graph(4)),("two_edges4",disconnected))
for (name,g) in cases, k in (1,2,4)
    probe("saran_$(name)_k$k") do
        original=copy(g)
        reg,trace=Q.k_partition_saran_vazirani(g,k)
        sizes=[count(==(i),values(reg)) for i in 1:k]
        (complete=Set(keys(reg))==Set(vertices(g)), sizes=sizes,
         balanced=sizes==fill(nv(g)÷k,k), unchanged=g==original,
         trace_length=length(trace), first_components=isempty(trace) ? nothing : first(trace)[1])
    end
end

probe("show_villager_captured") do
    mktemp() do _, out
        text=redirect_stdout(out) do
            sprint(show,Q.villager(1,1))
        end
        flush(out);seekstart(out)
        (requested_io=text,stdout=read(out,String))
    end
end

g=path_graph(4);p=Q.GraphPartition([1,1,2,2]);land=Q.quantumLand(p)
probe("plot_none") do
    figure,cost=Q.plot_partition(g,p;highlight=:none)
    (figure_type=string(typeof(figure)),cost=cost)
end
probe("legacy_plot_rank") do
    figure,cost=Q.visualize_graph_on_land(land,g;method=Q.cut_rank_cost)
    (figure_type=string(typeof(figure)),cost=cost)
end
probe("legacy_plot_cover_semantic_mismatch") do
    costs=Q.vertex_cover_cost_of_partition(land,g)
    figure,cost=Q.visualize_graph_on_land(land,g;method=Q.vertex_cover_cost_of_partition)
    (cover_returns=costs,second_value_type=string(typeof(costs[2])),
     figure_type=string(typeof(figure)),returned_cost=cost)
end
probe("figure_with_yticks_from_grid_plots") do
    f=Figure(size=(900,700),px_per_unit=6.0,fontsize=28.0,yticks=[0,20,100,200,300,400])
    string(typeof(f))
end

probe("sandbox_tuple_mean") do
    scores=[Q.cut_rank_cost(land,g),Q.cut_rank_cost(land,g)]
    mean(scores)
end
probe("grid_plots_tuple_scatter") do
    f=Figure();ax=Axis(f[1,1])
    scores=[Q.cut_rank_cost(land,g),Q.cut_rank_cost(land,g)]
    scatterlines!(ax,[4,4],scores)
    "constructed"
end

probe("empty_plot_none") do
    figure,cost=Q.plot_partition(SimpleGraph(0),Q.GraphPartition(Int[];nparts=2);highlight=:none)
    (figure_type=string(typeof(figure)),cost=cost)
end

println("INTEGRATION_PROBES_COMPLETE");flush(stdout)
