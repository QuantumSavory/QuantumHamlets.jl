# Run via run.py for external deadlines and non-overwriting result capture.
using QuantumHamlets, Graphs, Random, Statistics, BenchmarkTools, SHA
const Q = QuantumHamlets
include("reference.jl")
const R = ReferenceCosts
const OPERATION_NAMES = ("edge_cut_legacy", "cover_legacy", "rank_legacy", "matching_legacy",
                         "edge_cut_scalar", "cover_scalar", "rank_scalar", "matching_scalar")
const requested_operations = Set(ARGS)
issubset(requested_operations, Set(OPERATION_NAMES)) || error("unknown benchmark operation")

function fixtures()
    rng = Xoshiro(0x434f535450455246)
    entries = Tuple{String,SimpleGraph{Int},Vector{Int},Int}[]
    push!(entries, ("tiny_path6", path_graph(6), [1, 1, 1, 2, 2, 2], 2))
    push!(entries, ("path1000_alternating", path_graph(1000), [isodd(v) ? 1 : 2 for v in 1:1000], 2))
    push!(entries, ("grid400_stripes4", grid([20, 20]), [1 + (v - 1) ÷ 100 for v in 1:400], 4))
    push!(entries, ("grid400_interleaved4", grid([20, 20]), [1 + (v - 1) % 4 for v in 1:400], 4))
    push!(entries, ("sparse300_parts4", erdos_renyi(300, 0.02; rng), rand(rng, 1:4, 300), 4))
    push!(entries, ("dense120_parts2", erdos_renyi(120, 0.5; rng), [1 + (v - 1) ÷ 60 for v in 1:120], 2))
    push!(entries, ("empty_parts12_of16", path_graph(12), [isodd(v) ? 1 : 16 for v in 1:12], 16))
    push!(entries, ("empty_graph_parts3", SimpleGraph(0), Int[], 3))
    return entries
end

function operations(g, p, land)
    return (
        ("edge_cut_legacy", () -> R.naive_cost_of_partition(land, g), () -> Q.naive_cost_of_partition(land, g)),
        ("cover_legacy", () -> R.vertex_cover_cost_of_partition(land, g), () -> Q.vertex_cover_cost_of_partition(land, g)),
        ("rank_legacy", () -> R.cut_rank_cost(land, g), () -> Q.cut_rank_cost(land, g)),
        ("matching_legacy", () -> R.matching_cost_of_partition(land, g), () -> Q.matching_cost_of_partition(land, g)),
        ("edge_cut_scalar", () -> R.edge_cut_cost(g, p), () -> Q.edge_cut_cost(g, p)),
        ("cover_scalar", () -> R.vertex_cover_costs(g, p), () -> Q.vertex_cover_costs(g, p)),
        ("rank_scalar", () -> R.pairwise_cut_rank(g, p), () -> Q.pairwise_cut_rank(g, p)),
        ("matching_scalar", () -> R.matching_cost(g, p), () -> Q.matching_cost(g, p)),
    )
end

function verify_results(name, before, after, g, p)
    if name == "matching_legacy"
        @assert first(before) == first(after)
        for result in (before, after)
            @assert length(last(result)) == ne(g)
            selected = [e for (e, color) in zip(edges(g), last(result)) if color == :red]
            @assert length(selected) == first(result)
            labels = membership(p)
            for i in 1:p.nparts, j in i+1:p.nparts
                es = [e for e in selected if minmax(labels[src(e)], labels[dst(e)]) == (i, j)]
                @assert length(Set(v for e in es for v in (src(e), dst(e)))) == 2length(es)
            end
            @assert all(e -> labels[src(e)] != labels[dst(e)], selected)
        end
    else
        @assert before == after
    end
end

function bench_sample(benchmark)
    trial = run(benchmark)
    estimate = minimum(trial)
    return (ns=estimate.time, bytes=estimate.memory, allocations=estimate.allocs)
end

function compare(before, after)
    # Compilation is excluded. Calibrate each separately to target ~2ms/sample.
    before(); after()
    t_before = @elapsed before()
    t_after = @elapsed after()
    max(t_before, t_after) <= 2.0 || error("warm call exceeded 2s; refusing repeated workload")
    eval_before = clamp(floor(Int, 0.002 / max(t_before, 1e-9)), 1, 1000)
    eval_after = clamp(floor(Int, 0.002 / max(t_after, 1e-9)), 1, 1000)
    before_benchmark = @benchmarkable $before() evals=eval_before samples=1 seconds=1.0 gctrial=false gcsample=false
    after_benchmark = @benchmarkable $after() evals=eval_after samples=1 seconds=1.0 gctrial=false gcsample=false
    # Collect once before interleaving, rather than forcing 18 full collections.
    GC.gc()
    pre, post = NamedTuple[], NamedTuple[]
    for round in 1:9
        if isodd(round)
            push!(pre, bench_sample(before_benchmark))
            push!(post, bench_sample(after_benchmark))
        else
            push!(post, bench_sample(after_benchmark))
            push!(pre, bench_sample(before_benchmark))
        end
    end
    summarize(samples) = (median([s.ns for s in samples]),
                          median([s.bytes for s in samples]),
                          median([s.allocations for s in samples]))
    return summarize(pre), summarize(post), eval_before, eval_after
end

println("META\tjulia\t", VERSION)
println("META\tthreads\t", Threads.nthreads())
println("META\toperations\t", isempty(requested_operations) ? "all" : join(sort!(collect(requested_operations)), ","))
println("META\tcpu\t", Sys.CPU_NAME)
println("META\treference_sha256\t", bytes2hex(sha256(read(joinpath(@__DIR__, "reference_costs.jl")))))
println("META\tcurrent_costs_sha256\t", bytes2hex(sha256(read(joinpath(@__DIR__, "../../src/costs.jl")))))
println("META\tcurrent_api_sha256\t", bytes2hex(sha256(read(joinpath(@__DIR__, "../../src/api.jl")))))
println("HEADER\tfixture\tvertices\tedges\tparts\toperation\tbefore_ns\tafter_ns\tspeedup\tbefore_bytes\tafter_bytes\tbefore_allocs\tafter_allocs\tbefore_evals\tafter_evals")
println("READY"); flush(stdout)
for (fixture, g, labels, k) in fixtures()
    p = GraphPartition(labels; nparts=k)
    land = Q.quantumLand(p)
    original = copy(g)
    for (operation, before, after) in operations(g, p, land)
        isempty(requested_operations) || operation in requested_operations || continue
        println("STAGE\t", fixture, '\t', operation); flush(stdout)
        verify_results(operation, before(), after(), g, p)
        old, new, pre_evals, post_evals = compare(before, after)
        @assert g == original
        @assert membership(p) == labels
        fields = (fixture, nv(g), ne(g), k, operation, old[1], new[1], old[1]/new[1], old[2], new[2], old[3], new[3], pre_evals, post_evals)
        println("RESULT\t", join(fields, '\t')); flush(stdout)
    end
end
println("COMPLETE"); flush(stdout)
