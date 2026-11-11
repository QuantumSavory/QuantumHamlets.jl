# Run through run_matching_stress.py: the external parent limits every stage.
# Kept outside runtests.jl so the 50,000-vertex regression always runs in a
# disposable process with an external watchdog, including on slower machines.
using Graphs, QuantumHamlets, Test

function check_path_matching(stage, n)
    g = path_graph(n)
    p = GraphPartition([isodd(v) ? 1 : 2 for v in vertices(g)])
    elapsed = @elapsed if stage == "scalar"
        @test matching_cost(g, p) == n ÷ 2
    elseif stage == "witness"
        details = matching_details(g, p)
        @test details.cost == n ÷ 2
        @test length(details.edges) == n ÷ 2
        @test all(e -> has_edge(g, e), details.edges)
        @test length(Set(v for e in details.edges for v in (src(e), dst(e)))) == n
    else
        error("unknown stage $stage")
    end
    println("RESULT\t", stage, '\t', n, '\t', elapsed)
    flush(stdout)
end

# Compile before the small-case speed guard or large-case watchdog starts.
check_path_matching("scalar", 8)
check_path_matching("witness", 8)
println("READY")
flush(stdout)
for line in eachline(stdin)
    line == "QUIT" && break
    stage, value = split(line)
    n = parse(Int, value)
    n in (1_000, 50_000) || error("unsupported stress size")
    check_path_matching(stage, n)
end
