# Additional samples for the two largest cases, after the main benchmark process
# has finished. The initial 0.4-second budget admits too few slow baseline calls.
include("benchmarks.jl")

open(joinpath(OUT, "slow_results.tsv"), "w") do io
    open(joinpath(OUT, "slow_samples.tsv"), "w") do raw
        println(io, "graph\tkind\tn\tm\tparameter\timplementation\tmedian_ns\tp25_ns\tp75_ns\tmin_ns\tbytes\tallocations\tsamples")
        println(raw, "graph\tkind\tparameter\timplementation\tround\ttime_ns")
        g = grid([128, 128])
        for seeded in (false, true)
            measure_group(io, raw, "grid128x128", seeded ? "seeded" : "bury", g, nv(g)÷2,
                          single_functions(g, nv(g)÷2, seeded); samples=101, seconds=4.0)
        end
    end
end
println("SLOW_CASES_COMPLETE")
