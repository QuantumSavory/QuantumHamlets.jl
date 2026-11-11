# From the package root:
# julia --startup-file=no --project=experiments/bury_performance experiments/bury_performance/benchmarks.jl
using Graphs, QuantumHamlets, BenchmarkTools, Random, Statistics, Profile

include("scan.jl")
include("heap.jl")
include("multi.jl")
include("search.jl")

const Q = QuantumHamlets
const OUT = @__DIR__
const BENCH_SEED = 0x42555259
const ROUNDS = 3
const SAMPLES = 101
const SECONDS = 0.4

# This fixture reader is intentionally limited to the plain, unweighted GraphML
# files bundled in resources/QAOAgraphs_mini. Preserve their node-declaration order.
function qaoa_fixture(size)
    path = joinpath(OUT, "..", "..", "resources", "QAOAgraphs_mini", "qaoa_$(size).graphml")
    xml = read(path, String)
    ids = [m.captures[1] for m in eachmatch(r"<node\s+id=\"([^\"]+)\"", xml)]
    @assert length(unique(ids)) == length(ids)
    map = Dict(id => i for (i, id) in enumerate(ids))
    g = SimpleGraph(length(ids))
    for m in eachmatch(r"<edge\s+source=\"([^\"]+)\"\s+target=\"([^\"]+)\"", xml)
        @assert add_edge!(g, map[m.captures[1]], map[m.captures[2]])
    end
    @assert ne(g) == length(collect(eachmatch(r"<edge\s", xml)))
    return g
end

function measure_group(io, raw, name, kind, g, parameter, functions;
                       samples=SAMPLES, seconds=SECONDS)
    original = copy(g)
    Random.seed!(BENCH_SEED)
    expected = functions[1][2]()
    for (_, f) in functions
        Random.seed!(BENCH_SEED)
        @assert f() == expected "output mismatch in $name / $kind"
        @assert g == original
        # Warm all specializations before timing any competitor.
        f()
    end
    trials = Dict(label => Any[] for (label, _) in functions)
    for round in 1:ROUNDS
        # Rotate order to avoid giving one implementation a fixed position.
        for offset in 0:length(functions)-1
            i = mod1(round + offset, length(functions))
            label, f = functions[i]
            benchmark = @benchmarkable $f() setup=(Random.seed!($BENCH_SEED)) evals=1 samples=samples seconds=seconds
            trial = run(benchmark)
            push!(trials[label], trial)
            for time in trial.times
                println(raw, join((name, kind, parameter, label, round, time), '\t'))
            end
        end
    end
    for (label, _) in functions
        collected = trials[label]
        times = reduce(vcat, [t.times for t in collected])
        row = (name, kind, nv(g), ne(g), parameter, label, median(times),
               quantile(times, 0.25), quantile(times, 0.75), minimum(times),
               minimum(t.memory for t in collected), minimum(t.allocs for t in collected), length(times))
        println(io, join(row, '\t'))
    end
    flush(io); flush(raw)
    println(stderr, "Finished ", name, " / ", kind, " / ", parameter)
    flush(stderr)
    return nothing
end

function single_functions(g, target, seeded)
    return (
        ("current", seeded ? () -> Q.bury_seeded(g, target, 1; backend=:legacy) : () -> Q.bury(g, target; backend=:legacy)),
        ("scan", () -> BuryScan.bury_scan(g, target; seed=seeded ? 1 : nothing)),
        ("heap", () -> BuryHeap.bury_heap(g, target; seed=seeded ? 1 : nothing)),
    )
end

function multi_functions(g, k, seeded)
    return (
        ("current", seeded ? () -> Q.bury_k_seeded(g, k; backend=:legacy) : () -> Q.bury_k(g, k; backend=:legacy)),
        ("scan", () -> BuryMulti.bury_k_fast(g, k; seeded=seeded, kernel=BuryScan._colors)),
        ("heap", () -> BuryMulti.bury_k_fast(g, k; seeded=seeded, kernel=BuryHeap._colors)),
    )
end

function save_profile(f, name, io)
    f() # Warm the code before collecting samples.
    Profile.clear()
    deadline = time() + 1.0
    @profile while time() < deadline
        f()
    end
    println(io, "\nPROFILE ", name)
    Profile.print(io; format=:flat, sortedby=:count, mincount=5, C=false)
end

function main()
    rng = Xoshiro(BENCH_SEED)
    fixtures = [
        ("grid6x6", grid([6, 6])),
        ("grid32x32", grid([32, 32])),
        ("grid64x64", grid([64, 64])),
        ("grid128x128", grid([128, 128])),
        ("path4096", path_graph(4096)),
        ("star4096", star_graph(4096)),
        ("edgeless4096", SimpleGraph(4096)),
        ("regular3_4096", random_regular_graph(4096, 3; rng=rng)),
        ("regular12_4096", random_regular_graph(4096, 12; rng=rng)),
        ("dense512", erdos_renyi(512, 0.25; rng=rng)),
        ("complete512", complete_graph(512)),
        ("qaoa30", qaoa_fixture(30)),
        ("qaoa60", qaoa_fixture(60)),
    ]
    open(joinpath(OUT, "environment.txt"), "w") do io
        println(io, "Julia: ", VERSION)
        println(io, "CPU model: ", Sys.cpu_info()[1].model)
        println(io, "CPU target: ", Sys.CPU_NAME)
        println(io, "Architecture: ", Sys.ARCH)
        println(io, "Threads: ", Threads.nthreads())
        println(io, "Word size: ", Sys.WORD_SIZE)
        println(io, "Memory bytes: ", Sys.total_memory())
        println(io, "Graphs: ", pkgversion(Graphs))
        println(io, "BenchmarkTools: ", pkgversion(BenchmarkTools))
        println(io, "QuantumHamlets: ", pkgversion(Q))
        println(io, "RNG seed: ", BENCH_SEED)
        println(io, "Rounds: ", ROUNDS, "; samples/trial maximum: ", SAMPLES, "; seconds/trial: ", SECONDS)
        println(io, "evals=1; RNG reset in untimed setup; compiled code; all methods warmed; rotated order")
        println(io, "Graph construction/loading, plotting, and objective evaluation excluded from core partition timing.")
        println(io, "Search measurements include unchanged matching evaluation and identical improvement messages redirected to devnull.")
        println(io, "Memory/allocations are BenchmarkTools per-call estimates, not peak resident memory.")
    end
    open(joinpath(OUT, "results.tsv"), "w") do io
        open(joinpath(OUT, "samples.tsv"), "w") do raw
            println(io, "graph\tkind\tn\tm\tparameter\timplementation\tmedian_ns\tp25_ns\tp75_ns\tmin_ns\tbytes\tallocations\tsamples")
            println(raw, "graph\tkind\tparameter\timplementation\tround\ttime_ns")
            for (name, g) in fixtures, seeded in (false, true)
                target = nv(g) ÷ 2
                measure_group(io, raw, name, seeded ? "seeded" : "bury", g, target,
                              single_functions(g, target, seeded))
            end
            for target in (1024, 3072), seeded in (false, true)
                g = fixtures[3][2]
                measure_group(io, raw, "grid64x64", seeded ? "seeded" : "bury", g, target,
                              single_functions(g, target, seeded))
            end
            for (name, g, ks) in (("grid6x6", fixtures[1][2], [3]),
                                  ("grid64x64", fixtures[3][2], [2, 4, 8, 16]),
                                  ("regular3_4096", fixtures[8][2], [4, 16]),
                                  ("qaoa60", fixtures[13][2], [2, 7, 14]))
                for k in ks, seeded in (false, true)
                    measure_group(io, raw, name, seeded ? "multi_seeded" : "multi", g, k,
                                  multi_functions(g, k, seeded))
                end
            end
            # Preserve seed order and default RNG consumption. Silence both old
            # and experimental improvement messages outside the timed calls.
            redirect_stdout(devnull) do
                let g = grid([16, 16]), target = 128
                functions = (("current", () -> Q.bury_seeded(g, target; backend=:legacy)),
                             ("scan", () -> BurySearch.bury_auto_fast(g, target; kernel=BuryScan._colors)),
                             ("heap", () -> BurySearch.bury_auto_fast(g, target; kernel=BuryHeap._colors)))
                measure_group(io, raw, "grid16x16", "automatic_seed", g, target, functions)
                end
                let g = grid([8, 8]), target = 32
                functions = (("current", () -> Q.best_seed(g, target; backend=:legacy)),
                             ("scan", () -> BurySearch.best_seed_fast(g, target; kernel=BuryScan._colors)),
                             ("heap", () -> BurySearch.best_seed_fast(g, target; kernel=BuryHeap._colors)))
                measure_group(io, raw, "grid8x8", "all_seeds", g, target, functions)
                functions = (("current", () -> Q.bury_not_balanced_bipartition(g, 16, 32; backend=:legacy)),
                             ("scan", () -> BurySearch.bury_range_fast(g, 16, 32; kernel=BuryScan._colors)),
                             ("heap", () -> BurySearch.bury_range_fast(g, 16, 32; kernel=BuryHeap._colors)))
                measure_group(io, raw, "grid8x8", "range_search", g, 17, functions)
                end
            end
        end
    end
    open(joinpath(OUT, "profiles.txt"), "w") do io
        g = fixtures[4][2]
        save_profile(() -> Q.bury(g, nv(g)÷2; backend=:legacy), "current ordinary BURY grid128x128", io)
        save_profile(() -> Q.bury_seeded(g, nv(g)÷2, 1; backend=:legacy), "current seeded BURY grid128x128", io)
        redirect_stdout(devnull) do
            let g = grid([16, 16])
            save_profile(() -> Q.bury_seeded(g, nv(g)÷2; backend=:legacy), "current automatic seed grid16x16", io)
            end
        end
    end
    println("BENCHMARKS_COMPLETE")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
