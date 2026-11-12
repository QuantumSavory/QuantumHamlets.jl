from .matrix_tools import create_zero_matrix, insert_zero_matrix, copy_matrix, rank_matrix_positions, matrix_inverse, add_product_matrix

class GraphPartition:

    """
    A representation of a simple graph and a partition of the nodes of the graph into two sets, identified as the rows and columns.
    """

    nmb_nodes : int
    """The number of nodes in the graph."""

    nodes : list[int]
    """List of all nodes, i.e. the range from 0 to nmb_nodes exclusive."""

    row_flag : list[bool]
    """Flag telling which partition set each node belongs to. True for partition set 1 (rows), False for partition set 2 (columns)."""

    rows : list[int]
    """The nodes in the first partition set."""

    columns : list[int]
    """The nodes in the second partition set."""

    base_flag : list[bool]
    """Flag telling if a node represents a row or column in the selected invertible cut-rank submatrix of the adjacency matrix."""

    cut_rank : int
    """The cut-rank from the current partition."""

    base_rows : list[int]
    """The nodes from first partition set that represent rows in the invertible cut-rank submatrix of the adjacency matrix. Same as nodes n where row_flag[n] = True and base_falg[n] is True. Length should equal cut_rank."""

    base_columns : list[int]
    """The nodes from second partition set that represent columns in the invertible cut-rank submatrix of the adjacency matrix. Same as nodes n where row_flag[n] = False and base_falg[n] is True. Length should equal cut_rank."""

    free_rows : list[int]
    """The nodes from first partition set that represent rows outside the invertible cut-rank submatrix of the adjacency matrix. Same as nodes n where row_flag[n] = True and base_falg[n] is False."""

    free_columns : list[int]
    """The nodes from second partition set that represent columns outside the invertible cut-rank submatrix of the adjacency matrix. Same as nodes n where row_flag[n] = True and base_falg[n] is False."""

    adjacencies : list[list[int]]
    """The adjacency matrix of the grap. Should be a square symetric matrix with a row and column for each graph node, 0 on main diagonal, 1 in position (i,j) if (i,j) is an edge in the graph, 0 if not."""

    base_inverse : list[list[int]]
    """A square nmb_nodes x nmb_nodes matrix where the base_columns x base_rows submatrix is the inverse of the base_rows x base_columns submatrix of the adjacency matrix that defines the selected cut-rank sub-matrix of the current partition."""

    adj_b_inverse: list[list[int]]
    """A square nmb_nodes x nmb_nodes matrix where the nodes x base_rows submatrix represents 'D = A^{base_columns} * C^(-1)' used in the cut-rank calculations."""

    b_inverse_adj: list[list[int]]
    """A square nmb_nodes x nmb_nodes matrix where the base_columns x nodes submatrix represents 'E = C^(-1) * A_{base_rows}' used in the cut-rank calculations."""

    adj_b_inv_adj: list[list[int]]
    """The square nmb_nodes x nmb_nodes matrix 'F = A^{base_columns} * C^(-1) * A_{base_rows} + A' used in the cut-rank calculations."""

    buffer : list[list[int]]
    """A square nmb_nodes x nmb_nodes used for caching intermediate calculations when updating the variables after the partition has been changed."""


    def __init__(self, adjacencies : list[list[int]], partition_flags : list[bool]):
        self.adjacencies = adjacencies
        self.nmb_nodes = len(adjacencies)
        self.nodes = list(range(self.nmb_nodes))

        self.row_flag = partition_flags[:]
        self.rows = [n for n in self.nodes if self.row_flag[n]]
        self.columns = [n for n in self.nodes if not self.row_flag[n]]

        self._build_matrices()


    def _empty_matrix(self) -> list[list[int]]:

        return create_zero_matrix(self.nmb_nodes, self.nmb_nodes)


    def _build_matrices(self) -> None:

        self.base_inverse = self._empty_matrix()
        copy_matrix(self.adjacencies, self.base_inverse, self.rows, self.columns)
        (self.base_rows, self.base_columns) = rank_matrix_positions(self.base_inverse, self.rows, self.columns)
        self.cut_rank = len(self.base_rows)
        copy_matrix(self.adjacencies, self.base_inverse, self.base_rows, self.base_columns)
        matrix_inverse(self.base_inverse, self.base_inverse, self.base_rows, self.base_columns)

        self.adj_b_inverse = self._empty_matrix()
        add_product_matrix(self.adjacencies, self.base_inverse, self.adj_b_inverse, self.nodes, self.base_columns, self.base_rows)
        self.b_inverse_adj = self._empty_matrix()
        add_product_matrix(self.base_inverse, self.adjacencies, self.b_inverse_adj, self.base_columns, self.base_rows, self.nodes)
        self.adj_b_inv_adj = self._empty_matrix()
        copy_matrix(self.adjacencies, self.adj_b_inv_adj, self.nodes, self.nodes)
        add_product_matrix(self.adj_b_inverse, self.adjacencies, self.adj_b_inv_adj, self.nodes, self.base_rows, self.nodes)
        self.buffer = self._empty_matrix()

        self.base_flag = [False] * self.nmb_nodes
        for r in self.base_rows:
            self.base_flag[r] = True
        for c in self.base_columns:
            self.base_flag[c] = True
        self._build_free_nodes()


    def _build_free_nodes(self) -> None:

        self.free_rows = [row for row in self.rows if not self.base_flag[row]]
        self.free_columns = [col for col in self.columns if not self.base_flag[col]]


    def _reduce_base(self, removed_rows : list[int], removed_cols : list[int]) -> None:
    
        if len(removed_rows) == 0:
            return
        else:

            # Set base nodes
            for row in removed_rows:
                self.base_flag[row] = False
            for col in removed_cols:
                self.base_flag[col] = False
            self.base_rows = [row for row in self.rows if self.base_flag[row]]
            self.base_columns = [col for col in self.columns if self.base_flag[col]]
            self.cut_rank = len(self.base_rows)

            # Get Z
            copy_matrix(self.base_inverse, self.buffer, removed_cols, removed_rows)
            matrix_inverse(self.buffer, self.buffer, removed_cols, removed_rows)

            # Store D^(Delta X) * Z in D^(Delta Y), update D and F
            insert_zero_matrix(self.adj_b_inverse, self.nodes, removed_cols)
            add_product_matrix(self.adj_b_inverse, self.buffer, self.adj_b_inverse, self.nodes, removed_rows, removed_cols)
            add_product_matrix(self.adj_b_inverse, self.b_inverse_adj, self.adj_b_inv_adj, self.nodes, removed_cols, self.nodes)
            add_product_matrix(self.adj_b_inverse, self.base_inverse, self.adj_b_inverse, self.nodes, removed_cols, self.base_rows)

            # Store (C^-1)_YN^(Delta X) * Z in D^(Delta Y), update C^-1 and E
            insert_zero_matrix(self.adj_b_inverse, self.nodes, removed_cols)
            add_product_matrix(self.base_inverse, self.buffer, self.adj_b_inverse, self.base_columns, removed_rows, removed_cols)
            add_product_matrix(self.adj_b_inverse, self.b_inverse_adj, self.b_inverse_adj, self.base_columns, removed_cols, self.nodes)
            add_product_matrix(self.adj_b_inverse, self.base_inverse, self.base_inverse, self.base_columns, removed_cols, self.base_rows)


    def _extend_base(self, added_rows : list[int], added_cols : list[int]) -> None:
    
        if len(added_rows) == 0:
            return
        else:

            # Determine new base
            for row in added_rows:
                self.base_flag[row] = True
            for col in added_cols:
                self.base_flag[col] = True
            new_base_rows = [row for row in self.rows if self.base_flag[row]]
            new_base_columns = [col for col in self.columns if self.base_flag[col]]

			# Store Z in (C^-1)_(Delta Y)^(Delta X)
            copy_matrix(self.adjacencies, self.buffer, added_rows, added_cols)  # Stores (C_N)_(Delta X)^(Delta Y) in position for Z-inverse
            insert_zero_matrix(self.buffer, added_rows, self.base_rows)
            add_product_matrix(self.adjacencies, self.base_inverse, self.buffer, added_rows, self.base_columns, self.base_rows)
            add_product_matrix(self.buffer, self.adjacencies, self.buffer, added_rows, self.base_rows, added_cols)  # Gives Z-inverse
            insert_zero_matrix(self.base_inverse, added_cols, new_base_rows)
            insert_zero_matrix(self.base_inverse, self.base_columns, added_rows)
            matrix_inverse(self.buffer, self.base_inverse, added_rows, added_cols)

			# Get new C^1
            insert_zero_matrix(self.buffer, self.base_columns, added_cols)
            add_product_matrix(self.base_inverse, self.adjacencies, self.buffer, self.base_columns, self.base_rows, added_cols)
            add_product_matrix(self.base_inverse, self.buffer, self.base_inverse, added_cols, added_rows, self.base_rows)
            add_product_matrix(self.buffer, self.base_inverse, self.base_inverse, self.base_columns, added_cols, new_base_rows)

			# Get new D
            copy_matrix(self.adjacencies, self.adj_b_inverse, self.nodes, added_cols)
            add_product_matrix(self.adj_b_inverse, self.adjacencies, self.adj_b_inverse, self.nodes, self.base_rows, added_cols)  # D_O * C_XO^(Delta Y) + A^(Delta Y) stored in (Delta Y)-column of D
            insert_zero_matrix(self.adj_b_inverse, self.nodes, added_rows)
            add_product_matrix(self.adj_b_inverse, self.base_inverse, self.adj_b_inverse, self.nodes, added_cols, new_base_rows)

			# Get new E
            copy_matrix(self.adjacencies, self.b_inverse_adj, added_rows, self.nodes)
            add_product_matrix(self.adjacencies, self.b_inverse_adj, self.b_inverse_adj, added_rows, self.base_columns, self.nodes)  # C_(Delta X)^YO * E_0 + A_(Delta X) stored in (Delta X)-row of E
            insert_zero_matrix(self.b_inverse_adj, added_cols, self.nodes)
            add_product_matrix(self.base_inverse, self.b_inverse_adj, self.b_inverse_adj, new_base_columns, added_rows, self.nodes)

			# Get new F
            insert_zero_matrix(self.buffer, added_cols, self.nodes)
            add_product_matrix(self.base_inverse, self.b_inverse_adj, self.buffer, added_cols, added_rows, self.nodes)
            add_product_matrix(self.adj_b_inverse, self.buffer, self.adj_b_inv_adj, self.nodes, added_cols, self.nodes)

            self.base_rows = new_base_rows
            self.base_columns = new_base_columns
            self.cut_rank = len(self.base_rows)


    def apply_swap(self, row : int, column : int) -> None:

        remove_rows : list[int]
        remove_columns : list[int]
        add_rows : list[int]
        add_columns : list[int]

        if (not self.base_flag[column]):

            if (not self.base_flag[row]):

                # row in X^D, column in Y^D
                remove_rows = []
                remove_columns = []

                k2 = next((k2 for k2 in self.free_rows if k2 != row and self.adj_b_inv_adj[k2][row] == 1), -1)
                l2 = next((l2 for l2 in self.free_columns if l2 != column and self.adj_b_inv_adj[column][l2] == 1), -1)
                if k2 >= 0:
                    if l2 >= 0:
                        add_rows = [column, k2]
                        add_columns = [row, l2]
                    else:
                        add_rows = [k2]
                        add_columns = [row]
                else:
                    if l2 >= 0:
                        add_rows = [column]
                        add_columns = [l2]
                    else:
                        if self.adj_b_inv_adj[column][row] == 1:
                            add_rows = [column]
                            add_columns = [row]
                        else:
                            add_rows = []
                            add_columns = []

            else:

                # row in X^B, column in Y^D
                alpha = next(a for a in self.base_columns if self.base_inverse[a][row] == 1)
                remove_rows = [row]
                remove_columns = [alpha]

                k1 = next((k1 for k1 in self.free_rows if self.adj_b_inverse[k1][row] == 1), -1)
                if k1 >= 0:
                    if self.adj_b_inv_adj[k1][row] == 1:
                        k2 = next((k2 for k2 in self.free_rows if k2 != k1 and self.adj_b_inv_adj[k2][row] != self.adj_b_inverse[k2][row]), -1)
                    else:
                        k2 = next((k2 for k2 in self.free_rows if k2 != k1 and self.adj_b_inv_adj[k2][row] == 1), -1)
                    l2 = next((l2 for l2 in self.free_columns if l2 != column and self.adj_b_inv_adj[column][l2] == 1), -1)
                    if k2 >= 0:
                        if l2 >= 0:
                            add_rows = [column, k1, k2]
                            add_columns = [row, l2, alpha]
                        else:
                            add_rows = [k1, k2]
                            add_columns = [row, alpha]
                    else:
                        if l2 >= 0:
                            add_rows = [column, k1]
                            add_columns = [l2, alpha]
                        else:
                            if self.adj_b_inv_adj[column][row] != (self.adj_b_inverse[column][row] & self.adj_b_inv_adj[k1][row]):
                                add_rows = [column, k1]
                                add_columns = [row, alpha]
                            else:
                                add_rows = [k1]
                                add_columns = [alpha]
                else:
                    k2 = next((k2 for k2 in self.free_rows if self.adj_b_inv_adj[k2][row] == 1), -1)
                    if self.adj_b_inverse[column][row] == 1:
                        if k2 >= 0:
                            add_rows = [column, k2]
                            add_columns = [row, alpha]
                        else:
                            add_rows = [column]
                            add_columns = [alpha]
                    else:
                        l2 = next((l2 for l2 in self.free_columns if l2 != column and self.adj_b_inv_adj[column][l2] == 1), -1)
                        if k2 >= 0:
                            if l2 >= 0:
                                add_rows = [column, k2]
                                add_columns = [row, l2]
                            else:
                                add_rows = [k2]
                                add_columns = [row]
                        else:
                            if l2 >= 0:
                                add_rows = [column]
                                add_columns = [l2]
                            else:
                                if self.adj_b_inv_adj[column][row] == 1:
                                    add_rows = [column]
                                    add_columns = [row]
                                else:
                                    add_rows = []
                                    add_columns = []

        else:

            if (not self.base_flag[row]):

                # row in X^D, column in Y^B
                beta = next(b for b in self.base_rows if self.base_inverse[column][b] == 1)
                remove_rows = [beta]
                remove_columns = [column]

                l1 = next((l1 for l1 in self.free_columns if self.b_inverse_adj[column][l1] == 1), -1)
                if l1 >= 0:
                    if self.adj_b_inv_adj[column][l1] == 1:
                        l2 = next((l2 for l2 in self.free_columns if l2 != l1 and self.adj_b_inv_adj[column][l2] != self.b_inverse_adj[column][l2]), -1)
                    else:
                        l2 = next((l2 for l2 in self.free_columns if l2 != l1 and self.adj_b_inv_adj[column][l2] == 1), -1)
                    k2 = next((k2 for k2 in self.free_rows if k2 != row and self.adj_b_inv_adj[k2][row] == 1), -1)
                    if l2 >= 0:
                        if k2 >= 0:
                            add_rows = [column, k2, beta]
                            add_columns = [row, l1, l2]
                        else:
                            add_rows = [column, beta]
                            add_columns = [l1, l2]
                    else:
                        if k2 >= 0:
                            add_rows = [k2, beta]
                            add_columns = [row, l1]
                        else:
                            if self.adj_b_inv_adj[column][row] != (self.b_inverse_adj[column][row] & self.adj_b_inv_adj[column][l1]):
                                add_rows = [column, beta]
                                add_columns = [row, l1]
                            else:
                                add_rows = [beta]
                                add_columns = [l1]
                else:
                    l2 = next((l2 for l2 in self.free_columns if self.adj_b_inv_adj[column][l2] == 1), -1)
                    if self.b_inverse_adj[column][row] == 1:
                        if l2 >= 0:
                            add_rows = [column, beta]
                            add_columns = [row, l2]
                        else:
                            add_rows = [beta]
                            add_columns = [row]
                    else:
                        k2 = next((k2 for k2 in self.free_rows if k2 != row and self.adj_b_inv_adj[k2][row] == 1), -1)
                        if l2 >= 0:
                            if k2 >= 0:
                                add_rows = [column, k2]
                                add_columns = [row, l2]
                            else:
                                add_rows = [column]
                                add_columns = [l2]
                        else:
                            if k2 >= 0:
                                add_rows = [k2]
                                add_columns = [row]
                            else:
                                if self.adj_b_inv_adj[column][row] == 1:
                                    add_rows = [column]
                                    add_columns = [row]
                                else:
                                    add_rows = []
                                    add_columns = []

            else:

                # row in X^B, column in Y^B
                k1 = next((k1 for k1 in self.free_rows if self.adj_b_inverse[k1][row] == 1), -1)
                l1 = next((l1 for l1 in self.free_columns if self.b_inverse_adj[column][l1] == 1), -1)
                if (self.base_inverse[column][row] == 1):

                    # Full rank matrix with row and column removed is invertible
                    remove_rows = [row]
                    remove_columns = [column]

                    if k1 >= 0 and l1 >= 0:
                        if self.adj_b_inv_adj[k1][row] == 1:
                            k2 = next((k2 for k2 in self.free_rows if k2 != k1 and self.adj_b_inv_adj[k2][row] != self.adj_b_inverse[k2][row]), -1)
                        else:
                            k2 = next((k2 for k2 in self.free_rows if k2 != k1 and self.adj_b_inv_adj[k2][row] == 1), -1)
                        if self.adj_b_inv_adj[column][l1] == 1:
                            l2 = next((l2 for l2 in self.free_columns if l2 != l1 and self.adj_b_inv_adj[column][l2] != self.b_inverse_adj[column][l2]), -1)
                        else:
                            l2 = next((l2 for l2 in self.free_columns if l2 != l1 and self.adj_b_inv_adj[column][l2] == 1), -1)
                        if k2 >= 0:
                            if l2 >= 0:
                                add_rows = [column, k1, k2]
                                add_columns = [row, l1, l2]
                            else:
                                add_rows = [k1, k2]
                                add_columns = [row, l1]
                        else:
                            if l2 >= 0:
                                add_rows = [column, k1]
                                add_columns = [l1, l2]
                            else:
                                if ((self.adj_b_inv_adj[k1][row] & self.adj_b_inv_adj[column][l1]) ^ (self.adj_b_inv_adj[k1][row] & self.adj_b_inverse[column][row]) ^ (self.adj_b_inv_adj[column][l1] & self.b_inverse_adj[column][row])) != self.adj_b_inv_adj[column][row]:
                                    add_rows = [column, k1]
                                    add_columns = [row, l1]
                                else:
                                    add_rows = [k1]
                                    add_columns = [l1]

                    else:
                        k = next((k for k in self.free_rows if self.adj_b_inv_adj[k][row] != (self.adj_b_inverse[k][row] & self.b_inverse_adj[column][row])), -1)
                        l = next((l for l in self.free_columns if self.adj_b_inv_adj[column][l] != (self.adj_b_inverse[column][row] & self.b_inverse_adj[column][l])), -1)
                        if k >= 0:
                            if l >= 0:
                                add_rows = [column, k]
                                add_columns = [row, l]
                            else:
                                add_rows = [k]
                                add_columns = [row]
                        else:
                            if l >= 0:
                                add_rows = [column]
                                add_columns = [l]
                            else:
                                if self.adj_b_inv_adj[column][row] != (self.adj_b_inverse[column][row] & self.b_inverse_adj[column][row]):
                                    add_rows = [column]
                                    add_columns = [row]
                                else:
                                    add_rows = []
                                    add_columns = []

                else:

                    # Full rank matrix with row and column removed is singular
                    alpha = next(a for a in self.base_columns if self.base_inverse[a][row] == 1)
                    beta = next(b for b in self.base_rows if self.base_inverse[column][b] == 1)
                    remove_rows = [row, beta]
                    remove_columns = [column, alpha]

                    if k1 >= 0:

                        if l1 >= 0:

                            # Case k1 >= 0 and l1 >= 0
                            if self.adj_b_inv_adj[k1][row] == 1:
                                k2 = next((k2 for k2 in self.free_rows if k2 != k1 and self.adj_b_inv_adj[k2][row] != self.adj_b_inverse[k2][row]), -1)
                            else:
                                k2 = next((k2 for k2 in self.free_rows if k2 != k1 and self.adj_b_inv_adj[k2][row] == 1), -1)
                            if self.adj_b_inv_adj[column][l1] == 1:
                                l2 = next((l2 for l2 in self.free_columns if l2 != l1 and self.adj_b_inv_adj[column][l2] != self.b_inverse_adj[column][l2]), -1)
                            else:
                                l2 = next((l2 for l2 in self.free_columns if l2 != l1 and self.adj_b_inv_adj[column][l2] == 1), -1)
                            if k2 >= 0:
                                if l2 >= 0:
                                    add_rows = [column, k1, k2, beta]
                                    add_columns = [row, l1, l2, alpha]
                                else:
                                    add_rows = [k1, k2, beta]
                                    add_columns = [row, l1, alpha]
                            else:
                                if l2 >= 0:
                                    add_rows = [column, k1, beta]
                                    add_columns = [l1, l2, alpha]
                                else:
                                    if ((self.adj_b_inv_adj[k1][row] & self.adj_b_inverse[column][row]) ^ (self.adj_b_inv_adj[column][l1] & self.b_inverse_adj[column][row])) != self.adj_b_inv_adj[column][row]:
                                        add_rows = [column, k1, beta]
                                        add_columns = [row, l1, alpha]
                                    else:
                                        add_rows = [k1, beta]
                                        add_columns = [l1, alpha]

                        else:

                            # Case k1 >= 0 and l1 < 0
                            l2 = next((l2 for l2 in self.free_columns if self.adj_b_inv_adj[column][l2] == 1), -1)
                            if l2 >= 0:
                                if self.b_inverse_adj[column][row] == 1:
                                    add_rows = [column, k1, beta]
                                    add_columns = [row, l2, alpha]
                                else:
                                    if self.adj_b_inv_adj[k1][row] == 1:
                                        k2 = next((k2 for k2 in self.free_rows if k2 != k1 and self.adj_b_inv_adj[k2][row] != self.adj_b_inverse[k2][row]), -1)
                                    else:
                                        k2 = next((k2 for k2 in self.free_rows if k2 != k1 and self.adj_b_inv_adj[k2][row] == 1), -1)
                                    if k2 >= 0:
                                        add_rows = [column, k1, k2]
                                        add_columns = [row, l2, alpha]
                                    else:
                                        add_rows = [column, k1]
                                        add_columns = [l2, alpha]
                            else:
                                if self.b_inverse_adj[column][row] == 1:
                                    add_rows = [k1, beta]
                                    add_columns = [row, alpha]
                                else:
                                    if self.adj_b_inv_adj[k1][row] == 1:
                                        k2 = next((k2 for k2 in self.free_rows if k2 != k1 and self.adj_b_inv_adj[k2][row] != self.adj_b_inverse[k2][row]), -1)
                                    else:
                                        k2 = next((k2 for k2 in self.free_rows if k2 != k1 and self.adj_b_inv_adj[k2][row] == 1), -1)
                                    if k2 >= 0:
                                        add_rows = [k1, k2]
                                        add_columns = [row, alpha]
                                    else:
                                        if (self.adj_b_inv_adj[k1][row] & self.adj_b_inverse[column][row]) != self.adj_b_inv_adj[column][row]:
                                            add_rows = [column, k1]
                                            add_columns = [row, alpha]
                                        else:
                                            add_rows = [k1]
                                            add_columns = [alpha]

                    else:

                        if l1 >= 0:

                            # Case k1 < 0 and l1 >= 0
                            k2 = next((k2 for k2 in self.free_rows if self.adj_b_inv_adj[k2][row] == 1), -1)
                            if k2 >= 0:
                                if self.adj_b_inverse[column][row] == 1:
                                    add_rows = [column, k2, beta]
                                    add_columns = [row, l1, alpha]
                                else:
                                    if self.adj_b_inv_adj[column][l1] == 1:
                                        l2 = next((l2 for l2 in self.free_columns if l2 != l1 and self.adj_b_inv_adj[column][l2] != self.b_inverse_adj[column][l2]), -1)
                                    else:
                                        l2 = next((l2 for l2 in self.free_columns if l2 != l1 and self.adj_b_inv_adj[column][l2] == 1), -1)
                                    if l2 >= 0:
                                        add_rows = [column, k2, beta]
                                        add_columns = [row, l1, l2]
                                    else:
                                        add_rows = [k2, beta]
                                        add_columns = [row, l1]
                            else:
                                if self.adj_b_inverse[column][row] == 1:
                                    add_rows = [column, beta]
                                    add_columns = [l1, alpha]
                                else:
                                    if self.adj_b_inv_adj[column][l1] == 1:
                                        l2 = next((l2 for l2 in self.free_columns if l2 != l1 and self.adj_b_inv_adj[column][l2] != self.b_inverse_adj[column][l2]), -1)
                                    else:
                                        l2 = next((l2 for l2 in self.free_columns if l2 != l1 and self.adj_b_inv_adj[column][l2] == 1), -1)
                                    if l2 >= 0:
                                        add_rows = [column, beta]
                                        add_columns = [l1, l2]
                                    else:
                                        if (self.adj_b_inv_adj[column][l1] & self.b_inverse_adj[column][row]) != self.adj_b_inv_adj[column][row]:
                                            add_rows = [column, beta]
                                            add_columns = [row, l1]
                                        else:
                                            add_rows = [beta]
                                            add_columns = [l1]

                        else:

                            # Case k1 < 0 and l1 < 0
                            if self.adj_b_inverse[column][row] == 1:
                                if self.b_inverse_adj[column][row] == 1:
                                    add_rows = [column, beta]
                                    add_columns = [row, alpha]
                                else:
                                    k2 = next((k2 for k2 in self.free_rows if self.adj_b_inv_adj[k2][row] == 1), -1)
                                    if k2 >= 0:
                                        add_rows = [column, k2]
                                        add_columns = [row, alpha]
                                    else:
                                        add_rows = [column]
                                        add_columns = [alpha]
                            else:
                                if self.b_inverse_adj[column][row] == 1:
                                    l2 = next((l2 for l2 in self.free_columns if self.adj_b_inv_adj[column][l2] == 1), -1)
                                    if l2 >= 0:
                                        add_rows = [column, beta]
                                        add_columns = [row, l2]
                                    else:
                                        add_rows = [beta]
                                        add_columns = [row]
                                else:
                                    k2 = next((k2 for k2 in self.free_rows if self.adj_b_inv_adj[k2][row] == 1), -1)
                                    l2 = next((l2 for l2 in self.free_columns if self.adj_b_inv_adj[column][l2] == 1), -1)
                                    if k2 >= 0:
                                        if l2 >= 0:
                                            add_rows = [column, k2]
                                            add_columns = [row, l2]
                                        else:
                                            add_rows = [k2]
                                            add_columns = [row]
                                    else:
                                        if l2 >= 0:
                                            add_rows = [column]
                                            add_columns = [l2]
                                        else:
                                            if self.adj_b_inv_adj[column][row] == 1:
                                                add_rows = [column]
                                                add_columns = [row]
                                            else:
                                                add_rows = []
                                                add_columns = []

        # Set new set of rows and columns
        self.row_flag[row] = False
        self.row_flag[column] = True
        row_idx = next(i for i in range(len(self.rows)) if self.rows[i] == row)
        col_idx = next(i for i in range(len(self.columns)) if self.columns[i] == column)
        self.rows[row_idx] = column
        self.columns[col_idx] = row

        # Apply reduction and extension
        self._reduce_base(remove_rows, remove_columns)
        self._extend_base(add_rows, add_columns)

        # Update set of fre rows and free columns
        self._build_free_nodes()
