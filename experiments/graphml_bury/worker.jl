# Run through run_benchmarks.py: its external watchdog can stop native code too.
using Graphs, QuantumHamlets, BenchmarkTools, Random, Statistics
include("../bury_performance/scan.jl")
include("../bury_performance/heap.jl")
include("../bury_performance/search.jl")

const Q = QuantumHamlets
const OUT = @__DIR__
const BENCH_SEED = 0x47524150484d4c
const LABELS = ("legacy", "scan", "heap")

function event(parts...)
    println(stdout, join(parts, '\t'))
    flush(stdout)
end

function topology(path)
    open(path) do io
        n, m = parse.(Int, split(readline(io), '\t'))
        g = SimpleGraph(n)
        count = 0
        for line in eachline(io)
            u, v = parse.(Int, split(line, '\t'))
            @assert 1 <= u <= n && 1 <= v <= n && u != v
            @assert add_edge!(g, u, v)
            count += 1
        end
        @assert count == m == ne(g)
        return g
    end
end

function competitors(g, target, kind, seed)
    if kind == "ordinary"
        return (() -> Q.bury(g, target; backend=:legacy),
                () -> BuryScan.bury_scan(g, target),
                () -> BuryHeap.bury_heap(g, target))
    elseif kind == "seeded"
        return (() -> Q.bury_seeded(g, target, seed; backend=:legacy),
                () -> BuryScan.bury_scan(g, target; seed=seed),
                () -> BuryHeap.bury_heap(g, target; seed=seed))
    elseif kind == "automatic"
        return (() -> Q.bury_seeded(g, target; backend=:legacy),
                () -> BurySearch.bury_auto_fast(g, target; kernel=BuryScan._colors),
                () -> BurySearch.bury_auto_fast(g, target; kernel=BuryHeap._colors))
    end
    error("unknown workload $kind")
end

function run_job(id, path, kind, target, seed)
    event("START", id, "load", "all")
    g = topology(path)
    original = copy(g)
    fs = competitors(g, target, kind, seed)
    expected = nothing
    expected_next = nothing
    slow = false
    jobdir = joinpath(OUT, "jobs", id)
    mkpath(jobdir)
    open(joinpath(jobdir, "validation.tsv"), "w") do io
        println(io, "implementation\telapsed_seconds\tbytes\tpartition_equal\tgraph_unchanged\tbalanced\trng_equal")
        redirect_stdout(devnull) do
            for (label, f) in zip(LABELS, fs)
                # Protocol messages go to stderr while legacy search output is silenced.
                println(stderr, "START\t", id, "\tpreflight\t", label); flush(stderr)
                Random.seed!(BENCH_SEED)
                timed = @timed f()
                result = timed.value
                next_draw = rand(UInt64)
                if expected === nothing
                    expected = result
                    expected_next = next_draw
                end
                equal = result == expected
                unchanged = g == original
                balanced = length(result[1]) == nv(g) && sum(result[2]) == target &&
                           all(result[1][v] == Int(result[2][v]) + 1 for v in vertices(g))
                rng_equal = kind != "automatic" || next_draw == expected_next
                @assert equal && unchanged && balanced && rng_equal
                println(io, join((label, timed.time, timed.bytes, equal, unchanged, balanced, rng_equal), '\t'))
                flush(io)
                println(stderr, "CHECK\t", id, '\t', label, '\t', timed.time); flush(stderr)
                if timed.time > 1.0
                    slow = true
                    break
                end
            end
        end
    end
    if slow
        event("SLOW", id, "preflight_over_one_second")
        return
    end

    rounds = 3
    maxsamples = kind == "automatic" ? 11 : 101
    seconds = 0.2
    trials = [Any[] for _ in LABELS]
    event("START", id, "gc", "all")
    GC.gc() # One untimed full collection per group, not nine per tiny workload.
    open(joinpath(jobdir, "samples.tsv"), "w") do raw
        println(raw, "job\timplementation\tround\ttime_ns")
        redirect_stdout(devnull) do
            for round in 1:rounds, offset in 0:2
                i = mod1(round + offset, 3)
                f, label = fs[i], LABELS[i]
                println(stderr, "START\t", id, "\ttrial", round, '\t', label); flush(stderr)
                benchmark = @benchmarkable $f() setup=(Random.seed!($BENCH_SEED)) evals=1 samples=maxsamples seconds=seconds gctrial=false
                trial = run(benchmark)
                push!(trials[i], trial)
                for time in trial.times
                    println(raw, join((id, label, round, time), '\t'))
                end
                flush(raw)
            end
        end
    end
    open(joinpath(jobdir, "results.tsv"), "w") do io
        println(io, "job\tkind\tn\tm\ttarget\tseed\timplementation\tmedian_ns\tp25_ns\tp75_ns\tmin_ns\tbytes\tallocations\tsamples")
        for (i, label) in enumerate(LABELS)
            ts = trials[i]
            times = reduce(vcat, (t.times for t in ts))
            row = (id, kind, nv(g), ne(g), target, seed, label, median(times),
                   quantile(times, 0.25), quantile(times, 0.75), minimum(times),
                   minimum(t.memory for t in ts), minimum(t.allocs for t in ts), length(times))
            println(io, join(row, '\t'))
        end
    end
    @assert g == original
    event("DONE", id)
end

function warmup()
    g = path_graph(12)
    redirect_stdout(devnull) do
        for kind in ("ordinary", "seeded", "automatic")
            for f in competitors(g, 6, kind, 1)
                Random.seed!(BENCH_SEED)
                f()
                @timed f()
                benchmark = @benchmarkable $f() setup=(Random.seed!($BENCH_SEED)) evals=1 samples=2 seconds=0.01
                run(benchmark)
            end
        end
    end
end

open(joinpath(OUT, "environment.txt"), "w") do io
    println(io, "Julia: ", VERSION, "\nCPU: ", Sys.cpu_info()[1].model,
            "\nArchitecture: ", Sys.ARCH, "\nThreads: ", Threads.nthreads(),
            "\nMemory bytes: ", Sys.total_memory(), "\nGraphs: ", pkgversion(Graphs),
            "\nBenchmarkTools: ", pkgversion(BenchmarkTools), "\nQuantumHamlets: ", pkgversion(Q),
            "\nRNG seed: ", BENCH_SEED,
            "\nThree rotated rounds; evals=1; samples/trial<=101 (automatic<=11); seconds/trial=0.2.",
            "\nOne untimed full GC per workload group; gctrial=false; automatic GC remains enabled.",
            "\nCompilation and loading excluded; result registry and colors included.",
            "\nAutomatic search includes matching scoring; all improvement output redirected outside timed calls.",
            "\nAllocation figures are total allocated bytes/call, not peak memory.")
end
event("WARMING")
warmup()
event("READY")
for line in eachline(stdin)
    line == "QUIT" && break
    id, path, kind, target, seed = split(line, '\t')
    try
        run_job(id, path, kind, parse(Int, target), parse(Int, seed))
    catch ex
        showerror(stderr, ex, catch_backtrace()); println(stderr); flush(stderr)
        event("ERROR", id, string(typeof(ex)))
    end
end
