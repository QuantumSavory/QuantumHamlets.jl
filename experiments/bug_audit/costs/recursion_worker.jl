# Bounded exploration of the original matching evaluator's recursive DFS.
using Graphs, QuantumHamlets

function copied_dfs(g)
    # Exact function body/signature copied from src/costs.jl, lines 99–107.
    function flip_nodes(v::Int, g::Graphs.Graph, visted::Vector{Bool}, vtypes::Vector{Bool})
        for neighbor in neighbors(g,v)
            if !visted[neighbor]
                visted[neighbor] = true
                vtypes[neighbor] = !vtypes[v]
                flip_nodes(neighbor, g, visted, vtypes)
            end
        end
    end
    visited = fill(false,nv(g))
    vtypes = fill(false,nv(g))
    visited[1] = true
    flip_nodes(1,g,visited,vtypes)
    return count(visited)
end

function full_matching(g)
    p = GraphPartition([isodd(v) ? 1 : 2 for v in vertices(g)])
    return matching_cost(g,p)
end

copied_dfs(path_graph(8))
full_matching(path_graph(8))
println("READY");flush(stdout)
for line in eachline(stdin)
    line == "QUIT" && break
    kind, value = split(line,'\t')
    n = parse(Int,value)
    @assert 1 <= n <= 100000
    println("START\t",kind,'\t',n);flush(stdout)
    g = path_graph(n)
    try
        timed = @timed kind == "copied" ? copied_dfs(g) : full_matching(g)
        println("RESULT\t",kind,'\t',n,"\tOK\t",timed.time,'\t',timed.bytes,'\t',timed.value)
    catch ex
        println("RESULT\t",kind,'\t',n,'\t',typeof(ex),'\t',sprint(showerror,ex))
    end
    flush(stdout)
end
