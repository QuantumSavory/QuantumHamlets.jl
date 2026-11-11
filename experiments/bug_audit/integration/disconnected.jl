using QuantumHamlets, Graphs
const Q=QuantumHamlets

function graph_from_components(sizes)
    g=SimpleGraph(sum(sizes));first=1
    for size in sizes
        for v in first:first+size-2
            add_edge!(g,v,v+1)
        end
        first+=size
    end
    return g
end

function zero_cut_labels(sizes)
    n=sum(sizes)
    for mask in 1:(1<<length(sizes))-2
        sum(sizes[i] for i in eachindex(sizes) if !iszero(mask & (1<<(i-1)));init=0)==n÷2 || continue
        labels=Int[]
        for i in eachindex(sizes)
            append!(labels,fill(iszero(mask & (1<<(i-1))) ? 1 : 2,sizes[i]))
        end
        return labels
    end
    error("fixture lacks zero-cut bipartition")
end

println("READY");flush(stdout)
for sizes in ([3,3,2,2],[4,4,2,2],[2,2,2,2],[3,3,3,3],[3,3,2,2,1,1], [2,3,2,3])
    println("START\t",sizes);flush(stdout)
    g=graph_from_components(sizes)
    zero_labels=zero_cut_labels(sizes)
    known=Q.GraphPartition(zero_labels)
    @assert Q.edge_cut_cost(g,known)==0
    reg,trace=Q.k_partition_saran_vazirani(g,2)
    p=Q.GraphPartition(reg;nparts=2)
    println("RESULT\t",sizes,"\tcut=",Q.edge_cut_cost(g,p),"\tlabels=",Q.membership(p),
            "\tzero_cut_labels=",zero_labels,"\ttrace_length=",length(trace))
    println("END\t",sizes);flush(stdout)
end
println("DISCONNECTED_CHECKS_COMPLETE");flush(stdout)
