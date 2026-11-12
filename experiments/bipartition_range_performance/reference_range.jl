# Exact range-function body from commit a05adb5, isolated from production.
# Existing kernels and costs are shared: this measures only the range changes.
module BeforeRange
using Graphs
using QuantumHamlets: _bury_check_backend, _bury_check_target, bury, bury_seeded,
                      matching_cost_of_partition, quantumLand

function bury_not_balanced_bipartition(g, min_v, max_v; backend::Symbol=:heap,
                                      bury_v=bury, cost_f=matching_cost_of_partition)
    _bury_check_backend(backend)
    min_v = _bury_check_target(g, min_v)
    max_v = _bury_check_target(g, max_v)
    min_v <= max_v || throw(ArgumentError("minimum target must not exceed maximum target"))
    partitioner = if bury_v === bury
        (g, target) -> bury(g, target; backend=backend)
    elseif bury_v === bury_seeded
        (g, target) -> bury_seeded(g, target; backend=backend)
    else
        bury_v
    end
    n = nv(g)
    lowest_cost = nothing
    best_reg = nothing
    best_colors = nothing
    for v in min_v:max_v
        bury_reg, colors = partitioner(g, v)
        buryLand = quantumLand(bury_reg, 2, n)
        cost, _ = cost_f(buryLand, g)
        if best_reg === nothing || cost <= lowest_cost
            best_reg = bury_reg
            best_colors = colors
            lowest_cost = cost
        end
    end

    return best_reg, best_colors
end

end # module BeforeRange
