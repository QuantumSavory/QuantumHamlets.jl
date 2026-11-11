using Graphs, QuantumHamlets, Random, Test
const Q = QuantumHamlets

function probe(name, f)
    print(name, "\t")
    try
        result = redirect_stdout(f, devnull)
        println("RETURN\t", repr(result))
    catch e
        println("THROW\t", typeof(e), "\t", sprint(showerror,e))
    end
    flush(stdout)
end

g0 = SimpleGraph(0)
for backend in (:heap, :scan, :legacy)
    probe("best_seed empty $backend", () -> Q.best_seed(g0,0;backend))
    probe("best_sqrt_seed empty $backend", () -> Q.best_sqrt_seed(g0,0;backend))
    probe("bury_seeded empty $backend", () -> Q.bury_seeded(g0,0;backend))
end

# Rejecting invalid sampled-search targets should ideally not advance the RNG;
# record the actual effect separately from the valid-input contract.
for target in (-1, 5, 1.5)
    Random.seed!(712)
    expected = rand(UInt64)
    Random.seed!(712)
    err = try
        redirect_stdout(devnull) do
            Q.best_sqrt_seed(path_graph(4),target)
        end
        nothing
    catch e
        e
    end
    println("invalid sampled target ",target,"\t",typeof(err),"\trng_unchanged=",rand(UInt64)==expected)
end

# Safe integer representations and copy/round-trip boundaries.
@testset "Public supported-value boundaries" begin
    for n in (0,1,2,6)
        g=path_graph(n)
        for target in 0:n, T in (Int8,Int16,UInt8,UInt16,BigInt)
            expected=Q.bury(g,target;backend=:legacy)
            for backend in (:heap,:scan)
                @test Q.bury(g,T(target);backend)==expected
                if n>0
                    for seed in unique([1,n])
                        @test Q.bury_seeded(g,T(target),T(seed);backend)==Q.bury_seeded(g,target,seed;backend=:legacy)
                    end
                end
            end
        end
    end
    for labels in (Int[],[1],[1,2],[1,3,1]), extra in (0,1,3)
        k=(isempty(labels) ? 1 : maximum(labels))+extra
        p=Q.GraphPartition(labels;nparts=k)
        reg=Q.registry(p)
        land=Q.quantumLand(p)
        roundtrip=Q.GraphPartition(land)
        @test Q.membership(roundtrip)==labels
        @test Q.partition_sizes(roundtrip)==Q.partition_sizes(p)
        @test sum(Q.partition_sizes(p))==length(labels)
        @test Q.registry(roundtrip)==reg
        if !isempty(labels)
            copied=Q.membership(p);copied[1]=k
            reg[1]=k;land.registry[1]=k
            @test Q.membership(p)==labels
        end
    end
end
println("API_PROBES_COMPLETE")
