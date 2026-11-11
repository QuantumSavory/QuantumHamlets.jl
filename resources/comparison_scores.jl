# Shared objective for the grid and sandbox partitioner comparisons.
# Keep these functions independent of external partitioners and full sweeps.
using QuantumHamlets: GraphPartition, pairwise_cut_rank

const COMPARISON_OBJECTIVE_LABEL = "Pairwise cut rank (GF(2))"

comparison_score(g, labels_or_registry, nparts) =
    pairwise_cut_rank(g, GraphPartition(labels_or_registry; nparts=nparts))
