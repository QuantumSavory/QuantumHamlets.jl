"""The vendored min-cutrank annealer used by QuantumHamlets' Julia baselines."""

from math import isfinite
from random import Random

from .cut_rank_annealing import cut_rank_annealing_row_formula
from .graph_partition import GraphPartition


class _BestPartition(GraphPartition):
    """Observe completed swaps without changing the annealing trajectory."""

    def __init__(self, adjacency, flags):
        super().__init__(adjacency, flags)
        self.best_rank = self.cut_rank
        self.best_flags = self.row_flag[:]

    def apply_swap(self, row, column):
        super().apply_swap(row, column)
        if self.cut_rank < self.best_rank:
            self.best_rank = self.cut_rank
            self.best_flags = self.row_flag[:]


def anneal(adjacency, flags, temperatures, seed=None, keep_best=True):
    """Return Boolean membership after fixed-capacity cut-rank annealing.

    The Julia caller supplies the initial partition. ``True`` denotes a row
    vertex; each move swaps a row and column, preserving both part sizes.
    Inputs are copied. A local RNG leaves Python's global random state alone.
    By default retain the first strictly best state, including the initializer;
    ``keep_best=False`` instead returns the upstream algorithm's final state.
    Empty temperatures return a copy of the initializer without any sweeps.
    """
    adjacency = [[int(value) for value in row] for row in adjacency]
    flags = [bool(value) for value in flags]
    temperatures = [float(value) for value in temperatures]
    n = len(adjacency)
    if len(flags) != n or not 0 < sum(flags) < n:
        raise ValueError("initial membership must have two nonempty parts")
    if any(len(row) != n for row in adjacency):
        raise ValueError("adjacency must be square")
    if any(not isfinite(t) or t <= 0 for t in temperatures):
        raise ValueError("temperatures must contain only positive finite values")
    if not temperatures:
        return flags[:]
    partition = _BestPartition(adjacency, flags)
    cut_rank_annealing_row_formula(
        partition, temperatures, log=False, rng=Random(seed)
    )
    return (partition.best_flags if keep_best else partition.row_flag)[:]
