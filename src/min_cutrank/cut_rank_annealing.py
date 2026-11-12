from math import exp
import random
from .matrix_tools import copy_matrix, rank_matrix_positions
from .swap_rank_calculator import row_swap_cut_ranks

from .graph_partition import GraphPartition


def cut_rank_annealing_direct(partition : GraphPartition, temperatures, log: bool, rng=None) -> None:

    rng = random if rng is None else rng
    rows = partition.rows[:]
    cols = partition.columns[:]
    cut_rank = partition.cut_rank
    nmb_rows = len(rows)
    nmb_cols = len(cols)
    row_ranks = [0] * nmb_cols

    if log:
        print(f"Starting with cut-rank {cut_rank}")

    for temp in temperatures:
        limits = [exp(-1.0 / temp), exp(-2.0 / temp)]

        for i in range(nmb_rows):
            for j in range(nmb_cols):

                rows[i], cols[j] = cols[j], rows[i]
                copy_matrix(partition.adjacencies, partition.buffer, rows, cols)
                base_rows, _ = rank_matrix_positions(partition.buffer, rows, cols)
                row_ranks[j] = len(base_rows)
                rows[i], cols[j] = cols[j], rows[i]

            min_cutrank=min(row_ranks)
            if min_cutrank <= cut_rank or rng.random() < limits[min_cutrank - cut_rank - 1]:
                min_idx = [j for j in range(nmb_cols) if row_ranks[j] == min_cutrank]
                j = rng.choice(min_idx)
                partition.apply_swap(rows[i], cols[j])
                rows[i], cols[j] = cols[j], rows[i]
                cut_rank = min_cutrank

        if log:
            print(f"Cut-rank is {cut_rank} after sweep with temperature {temp}")


def cut_rank_annealing_row_formula(partition : GraphPartition, temperatures, log: bool, rng=None) -> None:

    rng = random if rng is None else rng
    rows = partition.rows[:]
    cols = partition.columns[:]
    row_ranks = [-1] * partition.nmb_nodes
    cut_rank = partition.cut_rank
    nmb_rows = len(rows)
    nmb_cols = len(cols)
    if log:
        print(f"Starting with cut-rank {cut_rank}")

    for temp in temperatures:
        limits = [exp(-1.0 / temp), exp(-2.0 / temp)]

        for i in range(nmb_rows):
            row = rows[i]
            for n in partition.nodes:
                row_ranks[n] = -1
            row_swap_cut_ranks(partition, row, row_ranks)
            min_cutrank=min(row_ranks[cols[j]] for j in range(nmb_cols))
            if min_cutrank <= cut_rank or rng.random() < limits[min_cutrank - cut_rank - 1]:
                min_idx = [j for j in range(nmb_cols) if row_ranks[cols[j]] == min_cutrank]
                j = rng.choice(min_idx)
                partition.apply_swap(row, cols[j])
                rows[i], cols[j] = cols[j], rows[i]
                cut_rank = min_cutrank

            if partition.cut_rank != cut_rank:
                raise Exception("Partition cut-rank does not fit with directly calculated rank")

        if log:
            print(f"Cut-rank is {cut_rank} after sweep with temperature {temp}")
