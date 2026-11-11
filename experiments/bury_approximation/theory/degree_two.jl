using Graphs, QuantumHamlets, Test
const Q = QuantumHamlets

# Verify the key guarantee used by the restricted-class theorem: every BURY
# output on a graph of maximum degree two has at most two crossing edges.
# Labels, disconnected components, targets, and explicit seeds are exhaustive.
graphs_checked = 0
partitions_checked = 0
@testset "BURY maximum-degree-two crossing-edge bound" begin
    for n in 0:6
        pairs = [(u, v) for u in 1:n for v in (u+1):n]
        for mask in 0:((1 << length(pairs))-1)
            g = SimpleGraph(n)
            for (i, (u, v)) in enumerate(pairs)
                if (mask >> (i-1)) & 1 == 1
                    add_edge!(g, u, v)
                end
            end
            all(Graphs.degree(g) .<= 2) || continue
            global graphs_checked += 1
            for target in 0:n
                _, colors = Q.bury(g, target)
                @test sum(colors) == target
                @test count(e -> colors[src(e)] != colors[dst(e)], edges(g)) <= 2
                global partitions_checked += 1
                for seed in 1:n
                    _, colors = Q.bury_seeded(g, target, seed)
                    @test sum(colors) == target
                    @test count(e -> colors[src(e)] != colors[dst(e)], edges(g)) <= 2
                    global partitions_checked += 1
                end
            end
        end
    end
end
println("graphs=", graphs_checked, "; partitions=", partitions_checked)
