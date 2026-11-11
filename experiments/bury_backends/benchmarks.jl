# Production API checks, separate from the historical prototype measurements.
# Run with --project=experiments/bury_performance from the package root.
include("../bury_performance/benchmarks.jl")

const INTEGRATION_OUT = @__DIR__
const PACKAGE_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const BACKENDS = (:legacy, :scan, :heap)

function load_topology(path)
    open(path) do io
        n, m = parse.(Int, split(readline(io), '\t'))
        g = SimpleGraph(n)
        for line in eachline(io)
            u, v = parse.(Int, split(line, '\t'))
            @assert u != v && add_edge!(g, u, v)
        end
        @assert ne(g) == m
        return g
    end
end

function fixture(directory, graph)
    return load_topology(joinpath(PACKAGE_ROOT, "experiments", "graphml_bury",
                                 "inventory", "topology", directory * "__" * graph * ".tsv"))
end

function production_calls(g, target, kind)
    return Tuple((String(backend),
                  kind == "ordinary" ? () -> Q.bury(g, target; backend=backend) :
                  kind == "seeded" ? () -> Q.bury_seeded(g, target, 1; backend=backend) :
                  kind == "multi" ? () -> Q.bury_k(g, target; backend=backend) :
                  () -> Q.bury_seeded(g, target; backend=backend)) for backend in BACKENDS)
end

function validate_graphml()
    path = joinpath(PACKAGE_ROOT, "experiments", "graphml_bury", "inventory", "inventory.tsv")
    lines = readlines(path)
    header = split(first(lines), '\t')
    column = findfirst(==("topology_path"), header)
    count = 0
    for line in lines[2:end]
        row = split(line, '\t')
        g = load_topology(joinpath(PACKAGE_ROOT, row[column]))
        original = copy(g)
        for kind in ("ordinary", "seeded")
            fs = production_calls(g, nv(g)÷2, kind)
            expected = fs[1][2]()
            for (_, f) in fs
                actual = f()
                @assert actual == expected
                @assert typeof(actual) === Tuple{Dict{Any,Any},Vector{Bool}}
                @assert g == original
            end
            count += 1
        end
    end
    println("GRAPHML_VALIDATED workloads=", count, " backends=3 exact_outputs_types_and_graph_preservation=true")
end

function integration_main()
    validate_graphml()
    fixtures = (("grid64x64", grid([64,64])),
                ("initial_barenco3", fixture("graphml_initial", "barenco_tof_3")),
                ("initial_barenco10", fixture("graphml_initial", "barenco_tof_10")),
                ("reduced_barenco3", fixture("graphml_reduced", "barenco_tof_3")),
                ("evolution_k5", fixture("graphml_evolution", "k5")))
    open(joinpath(INTEGRATION_OUT, "results.tsv"), "w") do io
        open(joinpath(INTEGRATION_OUT, "samples.tsv"), "w") do raw
            println(io, "graph\tkind\tn\tm\tparameter\timplementation\tmedian_ns\tp25_ns\tp75_ns\tmin_ns\tbytes\tallocations\tsamples")
            println(raw, "graph\tkind\tparameter\timplementation\tround\ttime_ns")
            for (name, g) in fixtures, kind in ("ordinary", "seeded")
                measure_group(io, raw, name, kind, g, nv(g)÷2,
                              production_calls(g, nv(g)÷2, kind); samples=51, seconds=0.3)
            end
            g = grid([32,32])
            measure_group(io, raw, "grid32x32", "multi", g, 4,
                          production_calls(g, 4, "multi"); samples=51, seconds=0.3)
            redirect_stdout(devnull) do
                g = grid([8,8])
                measure_group(io, raw, "grid8x8", "automatic", g, nv(g)÷2,
                              production_calls(g, nv(g)÷2, "automatic"); samples=21, seconds=0.3)
            end
        end
    end
    open(joinpath(INTEGRATION_OUT, "environment.txt"), "w") do io
        println(io, "Julia: ", VERSION, "\nCPU: ", Sys.cpu_info()[1].model,
                "\nThreads: ", Threads.nthreads(), "\nGraphs: ", pkgversion(Graphs),
                "\nBenchmarkTools: ", pkgversion(BenchmarkTools),
                "\nProduction public APIs with Dict{Any,Any} outputs; original multipartition driver.",
                "\nThree rotated rounds, evals=1, samples<=51 (automatic<=21), seconds/trial=0.3.",
                "\nFull untimed BenchmarkTools GC per trial; automatic GC remains enabled.",
                "\nCompilation/loading excluded, result construction included; automatic includes matching.",
                "\nAutomatic RNG seed reset before each sample: ", BENCH_SEED,
                "\nNo automatic seed search on large GraphML fixtures.")
    end
    println("PRODUCTION_BACKEND_BENCHMARKS_COMPLETE")
end

integration_main()
