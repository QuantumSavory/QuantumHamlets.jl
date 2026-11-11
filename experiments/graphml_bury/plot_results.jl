# Uses the plotting dependency already present in QuantumHamlets.
using QuantumHamlets
using QuantumHamlets.CairoMakie

const OUT = @__DIR__
lines = readlines(joinpath(OUT, "comparison.tsv"))
header = split(first(lines), '\t')
rows = [Dict(zip(header, split(line, '\t'))) for line in lines[2:end]]
styles = (("graphml_initial", "Initial", "#2166ac", :circle),
          ("graphml_reduced", "Reduced", "#c76e13", :rect),
          ("graphml_evolution", "Evolution", "#3b8b5a", :utriangle))
fig = Figure(size=(1100, 500), fontsize=15)
axes = Axis[]
for (column, kind, title) in ((1, "ordinary", "Ordinary BURY"),
                            (2, "seeded", "Explicitly seeded BURY"))
    ax = Axis(fig[1, column], title=title, xscale=log10, yscale=log10,
              xlabel="Vertices (log scale)",
              ylabel=column==1 ? "Scan median time / heap median time" : "",
              xticks=([25,100,500,2000],["25","100","500","2,000"]),
              yticks=([.5,1,2,4,8],["0.5×","1×","2×","4×","8×"]),
              limits=((20,2200),(.28,12)))
    hspan!(ax, 1/1.1, 1.1, color=(:gray,.15))
    hlines!(ax, [1], color=:gray, linewidth=1)
    for (directory, label, color, marker) in styles
        rs = filter(r -> r["directory"]==directory && r["kind"]==kind && r["half"]=="True", rows)
        scatter!(ax, parse.(Int,[r["n"] for r in rs]),
                 parse.(Float64,[r["scan_over_heap"] for r in rs]);
                 label, color, marker, markersize=11, strokecolor=:white, strokewidth=.6)
    end
    text!(ax, 24, 10; text="Heap faster", fontsize=13, color=:gray40)
    text!(ax, 24, .34; text="Scan faster", fontsize=13, color=:gray40)
    push!(axes, ax)
end
Legend(fig[0, 1:2], axes[1], orientation=:horizontal, framevisible=false)
Label(fig[2, 1:2], "45 GraphML files • target = floor(n/2) • explicit seed = first declared node • gray band = within 10%", fontsize=13)
save(joinpath(OUT, "backend_comparison.png"), fig, px_per_unit=1.6)
save(joinpath(OUT, "backend_comparison.pdf"), fig)
println("Saved backend_comparison.png and .pdf")
