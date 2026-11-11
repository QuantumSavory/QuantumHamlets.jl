# A Boolean entrywise elimination oracle, independent of packed-word arithmetic.
function boolean_rank_oracle(input::AbstractMatrix{Bool})
    work = Matrix{Bool}(input)
    nrows, ncols = size(work)
    pivotrow = 1
    for col in 1:ncols
        found = findfirst(row -> work[row, col], pivotrow:nrows)
        isnothing(found) && continue
        selected = pivotrow + found - 1
        for j in 1:ncols
            work[pivotrow, j], work[selected, j] = work[selected, j], work[pivotrow, j]
        end
        for row in (pivotrow + 1):nrows
            if work[row, col]
                for j in col:ncols
                    work[row, j] = xor(work[row, j], work[pivotrow, j])
                end
            end
        end
        pivotrow += 1
        pivotrow > nrows && break
    end
    return pivotrow - 1
end

function graph_from_biadjacency(matrix::AbstractMatrix{Bool})
    nleft, nright = size(matrix)
    g = SimpleGraph(nleft + nright)
    for i in 1:nleft, j in 1:nright
        matrix[i, j] && add_edge!(g, i, nleft + j)
    end
    return g, Q.GraphPartition(vcat(fill(1, nleft), fill(2, nright)); nparts=3)
end

@testset "Exact cut rank across packed word boundaries" begin
    rng = Xoshiro(0x63757472616e6b)
    for (nleft, nright) in ((0, 0), (0, 65), (1, 129), (63, 67), (64, 65),
                            (65, 64), (65, 65), (67, 129), (129, 67), (130, 131))
        random_matrix = rand(rng, Bool, nleft, nright)
        identity_matrix = [i == j for i in 1:nleft, j in 1:nright]
        repeated_matrix = [isodd(i) && isodd(j) for i in 1:nleft, j in 1:nright]
        for input in (random_matrix, identity_matrix, repeated_matrix, falses(nleft, nright))
            g, p = graph_from_biadjacency(input)
            original = copy(g)
            expected = boolean_rank_oracle(input)
            @test Q.pairwise_cut_rank(g, p) == expected
            @test Q.cut_rank_cost(Q.quantumLand(p), g) == (expected, nothing)
            @test g == original
            @test Q.membership(p) == vcat(fill(1, nleft), fill(2, nright))
        end
    end

    # Nonzero pivots at and after bit positions 64 and 128, with duplicate rows.
    matrix = falses(140, 140)
    for col in (1, 63, 64, 65, 127, 128, 129, 140)
        matrix[col, col] = true
        matrix[2, col] = true
        matrix[3, col] = true
    end
    g, p = graph_from_biadjacency(matrix)
    @test Q.pairwise_cut_rank(g, p) == boolean_rank_oracle(matrix) == 8

    # More than two parts: check the sum independently and retain empty labels.
    g = erdos_renyi(140, 0.18; rng)
    labels = [v <= 66 ? 1 : v <= 136 ? 3 : 4 for v in vertices(g)]
    p = Q.GraphPartition(labels; nparts=5)
    expected = 0
    for i in 1:5, j in (i + 1):5
        left, right = findall(==(i), labels), findall(==(j), labels)
        expected += boolean_rank_oracle(Bool[has_edge(g, u, v) for u in left, v in right])
    end
    @test Q.pairwise_cut_rank(g, p) == expected
    @test Q.cut_rank_cost(Q.quantumLand(p), g) == (expected, nothing)
    @test Q.pairwise_cut_rank(g, Q.GraphPartition(ones(Int, nv(g)))) == 0

    # Only actual cut pairs need work, even with many declared empty parts or
    # many active singleton parts whose pair graph is sparse.
    @test Q.pairwise_cut_rank(path_graph(2), Q.GraphPartition([1, 1_000_000]; nparts=1_000_000)) == 1
    @test Q.pairwise_cut_rank(path_graph(128), Q.GraphPartition(collect(1:128))) == 127
end
