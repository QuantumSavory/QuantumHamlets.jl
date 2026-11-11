# Independent small-graph numerical checks. No production source is modified.
using Graphs, QuantumHamlets, Random, Test
const Q = QuantumHamlets

function matching_oracle(g, left, right)
    function search(i, used)
        i > length(left) && return 0
        best = search(i + 1, used)
        for (j, v) in enumerate(right)
            bit = UInt64(1) << (j - 1)
            if iszero(used & bit) && has_edge(g, left[i], v)
                best = max(best, 1 + search(i + 1, used | bit))
            end
        end
        return best
    end
    return search(1, UInt64(0))
end

function rank_oracle(g, left, right)
    # Gaussian elimination over GF(2) using packed rows and XOR only.
    pivots = zeros(UInt64, length(right))
    rank = 0
    for u in left
        row = UInt64(0)
        for (j, v) in enumerate(right)
            has_edge(g, u, v) && (row |= UInt64(1) << (j-1))
        end
        while !iszero(row)
            pivot = 64 - leading_zeros(row)
            if iszero(pivots[pivot])
                pivots[pivot] = row
                rank += 1
                break
            end
            row ⊻= pivots[pivot]
        end
    end
    return rank
end

function cover_oracle(n, cross_edges)
    best = n
    for mask in UInt64(0):(UInt64(1)<<n)-1
        count_ones(mask) >= best && continue
        if all(e -> !iszero(mask & ((UInt64(1)<<(src(e)-1)) | (UInt64(1)<<(dst(e)-1)))), cross_edges)
            best = count_ones(mask)
        end
    end
    return best
end

function check_fixture(g, labels, k)
    p = GraphPartition(labels; nparts=k)
    original = copy(g)
    original_labels = copy(labels)
    cross = [e for e in edges(g) if labels[src(e)] != labels[dst(e)]]
    pair_matching = Dict{Tuple{Int,Int},Int}()
    rank = 0
    for i in 1:k, j in i+1:k
        left, right = findall(==(i),labels), findall(==(j),labels)
        pair_matching[(i,j)] = matching_oracle(g, left, right)
        rank += rank_oracle(g, left, right)
    end
    expected_matching = sum(values(pair_matching); init=0)
    @test edge_cut_cost(g,p) == length(cross)
    @test pairwise_cut_rank(g,p) == rank
    @test matching_cost(g,p) == expected_matching
    details = matching_details(g,p)
    @test details.cost == expected_matching
    @test length(details.edges) == expected_matching
    @test length(unique(details.edges)) == length(details.edges)
    @test all(e -> has_edge(g,e) && labels[src(e)] != labels[dst(e)],details.edges)
    for ((i,j),expected) in pair_matching
        es = [e for e in details.edges if minmax(labels[src(e)],labels[dst(e)]) == (i,j)]
        @test length(es) == expected
        @test length(Set(v for e in es for v in (src(e),dst(e)))) == 2length(es)
    end
    cover = vertex_cover_costs(g,p)
    @test cover_oracle(nv(g),cross) <= cover.ghz_free <= nv(g)
    @test expected_matching <= cover.ghz_n_minus_one <= (k-1)*cover.ghz_free
    @test 0 <= rank <= expected_matching <= length(cross)
    @test g == original
    @test labels == original_labels == membership(p)
end

function main()
    count = 0
    # Known failure is recorded separately, not rediscovered thousands of times.
    try
        println("ONE_PART_RESULT=",matching_cost(path_graph(3),GraphPartition([1,1,1])))
    catch ex
        println("KNOWN_ONE_PART_ERROR=",sprint(showerror,ex))
    end
    @testset "Independent cost oracles" begin
        @testset "Empty graphs and empty declared parts" begin
            for k in (2,3,4)
                check_fixture(SimpleGraph(0),Int[],k)
                count += 1
            end
        end
        @testset "Every labeled graph and 2/3-part membership through n=4" begin
            for n in 1:4
                possible = [(u,v) for u in 1:n for v in u+1:n]
                for mask in 0:(1<<length(possible))-1
                    g = SimpleGraph(n)
                    for (i,(u,v)) in enumerate(possible)
                        !iszero(mask & (1<<(i-1))) && add_edge!(g,u,v)
                    end
                    for k in (2,3), code in 0:k^n-1
                        labels = digits(code; base=k,pad=n) .+ 1
                        check_fixture(g,labels,k)
                        count += 1
                    end
                end
                println("EXHAUSTIVE_FINISHED n=",n," fixtures=",count); flush(stdout)
            end
        end
        @testset "Random n=5/6 with empty parts and relabeling" begin
            rng = Xoshiro(0x434f535441554449)
            for n in (5,6), _ in 1:150
                g = erdos_renyi(n,rand(rng);rng=rng)
                k = rand(rng,2:5)
                labels = rand(rng,1:k,n)
                check_fixture(g,labels,k)
                permutation = randperm(rng,n)
                relabeled = SimpleGraph(n)
                for e in edges(g)
                    add_edge!(relabeled,permutation[src(e)],permutation[dst(e)])
                end
                label_permutation = randperm(rng,k)
                relabeled_labels = similar(labels)
                for v in vertices(g)
                    relabeled_labels[permutation[v]] = label_permutation[labels[v]]
                end
                check_fixture(relabeled,relabeled_labels,k)
                a,b = GraphPartition(labels;nparts=k),GraphPartition(relabeled_labels;nparts=k)
                @test edge_cut_cost(g,a) == edge_cut_cost(relabeled,b)
                @test matching_cost(g,a) == matching_cost(relabeled,b)
                @test pairwise_cut_rank(g,a) == pairwise_cut_rank(relabeled,b)
                count += 2
            end
        end
    end
    println("COST_ORACLE_COMPLETE fixtures=",count)
end

main()
