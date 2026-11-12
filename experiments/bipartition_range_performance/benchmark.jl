# Use run.py for external deadlines and non-overwriting evidence capture.
using QuantumHamlets, Graphs, Random, Statistics, BenchmarkTools, SHA
const Q = QuantumHamlets
include("reference_range.jl")
const R = BeforeRange

function fixtures()
    rng = Xoshiro(0x52414e4745504552)
    return (
        ("path11", path_graph(11)),
        ("grid225", grid([15, 15])),
        ("sparse201", erdos_renyi(201, 0.025; rng)),
        ("dense61", erdos_renyi(61, 0.35; rng)),
    )
end

function workloads(g)
    center = fld(nv(g), 2)
    ranges = (("singleton", center:center),
              ("near_balanced", center:cld(nv(g), 2)),
              ("wider", center-2:center+3))
    costs = (("matching", Q.matching_cost_of_partition, matching_cost),
             ("edge_cut", Q.naive_cost_of_partition, edge_cut_cost),
             ("rank", Q.cut_rank_cost, pairwise_cut_rank))
    return ranges, costs
end

function methods(g, targets, backend, legacy_cost, scalar_cost, mode)
    # Exact callbacks exercise the fast path. Closures deliberately exercise
    # the generic callback path and its retained-output ownership work.
    partitioner = mode == "builtin" ? Q.bury : (g, r) -> Q.bury(g, r; backend)
    tuple_cost = mode == "builtin" ? legacy_cost : (land, g) -> legacy_cost(land, g)
    objective = mode == "builtin" ? scalar_cost : (g, p) -> scalar_cost(g, p)
    return (
        () -> R.bury_not_balanced_bipartition(g, first(targets), last(targets);
                                             backend, bury_v=partitioner, cost_f=tuple_cost),
        () -> Q.bury_not_balanced_bipartition(g, first(targets), last(targets);
                                             backend, bury_v=partitioner, cost_f=tuple_cost),
        () -> bipartition_graph(g; targets, backend, objective),
    )
end

function verify(calls, g, scalar_cost)
    before, after, helper = map(f -> f(), calls)
    @assert before == after
    @assert before isa Tuple{Dict{Any,Any},Vector{Bool}}
    @assert after isa Tuple{Dict{Any,Any},Vector{Bool}}
    @assert helper isa GraphPartition
    @assert helper.nparts == 2
    @assert membership(helper) == [before[1][v] for v in vertices(g)]
    @assert before[2] == (membership(helper) .== 2)
    @assert scalar_cost(g, helper) == scalar_cost(g, GraphPartition(before[1]; nparts=2))
    return scalar_cost(g, helper), count(before[2])
end

function compare(calls)
    foreach(f -> f(), calls)
    elapsed = map(f -> @elapsed(f()), calls)
    maximum(elapsed) <= 2.0 || error("warm call exceeded 2s; refusing repeated workload")
    evals = map(t -> clamp(floor(Int, 0.002 / max(t, 1e-9)), 1, 1000), elapsed)
    benchmarks = map(calls, evals) do f, count
        @benchmarkable $f() evals=count samples=1 seconds=1.0 gctrial=false gcsample=false
    end
    samples = [NamedTuple[] for _ in calls]
    # Nine rounds rotate all three methods through every timing position.
    GC.gc()
    for round in 1:9
        for offset in 0:2
            index = mod1(round + offset, 3)
            estimate = minimum(run(benchmarks[index]))
            push!(samples[index], (round=round, ns=estimate.time,
                                   bytes=estimate.memory, allocations=estimate.allocs))
        end
    end
    summaries = map(samples) do sample
        (ns=median([s.ns for s in sample]), bytes=median([s.bytes for s in sample]),
         allocations=median([s.allocations for s in sample]))
    end
    return summaries, samples, evals
end

println("META\tjulia\t", VERSION)
println("META\tthreads\t", Threads.nthreads())
println("META\tcpu\t", Sys.CPU_NAME)
println("META\tbaseline_commit\ta05adb5")
for filename in ("reference_range.jl", "benchmark.jl", "../../src/bury.jl",
                 "../../src/api.jl", "../../src/partition.jl", "../../src/costs.jl",
                 "../../src/bury/heap.jl", "../../src/bury/scan.jl", "../../src/bury/legacy.jl")
    println("META\tsha256\t", filename, '\t', bytes2hex(sha256(read(joinpath(@__DIR__, filename)))))
end
println("HEADER\tfixture\tvertices\tedges\trange\tfirst_target\tlast_target\tbackend\tobjective\tmode\twinner_target\tscore\tbefore_ns\tafter_ns\thelper_ns\tafter_speedup\thelper_speedup\tbefore_bytes\tafter_bytes\thelper_bytes\tbefore_allocs\tafter_allocs\thelper_allocs\tbefore_evals\tafter_evals\thelper_evals")
println("SAMPLE_HEADER\tfixture\trange\tbackend\tobjective\tmode\tmethod\tround\tns\tbytes\tallocations")
println("READY"); flush(stdout)
for (fixture, g) in fixtures()
    original = copy(g)
    ranges, costs = workloads(g)
    for (range_name, targets) in ranges, backend in (:heap, :scan, :legacy),
        (objective_name, legacy_cost, scalar_cost) in costs, mode in ("builtin", "callbacks")
        mode == "callbacks" && (backend != :heap || objective_name != "matching") && continue
        println("STAGE\t", join((fixture, range_name, backend, objective_name, mode), '\t')); flush(stdout)
        calls = methods(g, targets, backend, legacy_cost, scalar_cost, mode)
        score, winner = verify(calls, g, scalar_cost)
        summaries, samples, evals = compare(calls)
        @assert g == original
        old, new, helper = summaries
        fields = (fixture, nv(g), ne(g), range_name, first(targets), last(targets), backend,
                  objective_name, mode, winner, score, old.ns, new.ns, helper.ns,
                  old.ns/new.ns, old.ns/helper.ns, old.bytes, new.bytes, helper.bytes,
                  old.allocations, new.allocations, helper.allocations, evals...)
        println("RESULT\t", join(fields, '\t'))
        for (method, sample) in zip(("before", "after", "helper"), samples), s in sample
            println("SAMPLE\t", join((fixture, range_name, backend, objective_name, mode,
                                      method, s.round, s.ns, s.bytes, s.allocations), '\t'))
        end
        flush(stdout)
    end
end
println("COMPLETE"); flush(stdout)
