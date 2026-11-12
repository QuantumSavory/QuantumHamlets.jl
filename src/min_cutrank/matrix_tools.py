def create_zero_matrix(nmb_rows : int, nmb_columns : int) -> list[list[int]]:

    return [[0] * nmb_columns for _ in range(nmb_rows)]


def insert_zero_matrix(matrix : list[list[int]], rows : list[int], columns : list[int]) -> None:

    for row in rows:
        for col in columns:
            matrix[row][col] = 0


def is_zero_matrix(matrix : list[list[int]], rows : list[int], columns : list[int]) -> bool:

    for row in rows:
        for col in columns:
            if matrix[row][col] == 1:
                return False
    return True

def is_identity_matrix(matrix: list[list[int]], rows_and_columns : list[int]) -> bool:

    for row in rows_and_columns:
        for col in rows_and_columns:
            if matrix[row][col] != (1 if row == col else 0):
                return False
    return True

def set_common_matrix_value(value : int, matrix : list[list[int]], rows : list[int], columns : list[int]) -> None:

    for row in rows:
        for col in columns:
            matrix[row][col] = value


def copy_matrix(from_mat : list[list[int]], to_mat : list[list[int]], rows : list[int], columns : list[int]) -> None:

    for row in rows:
        for col in columns:
            to_mat[row][col] = from_mat[row][col]


def add_matrix(from_mat : list[list[int]], to_mat : list[list[int]], rows : list[int], columns : list[int]) -> None:

    for row in rows:
        for col in columns:
            to_mat[row][col] += from_mat[row][col]


def add_product_matrix(fac1 : list[list[int]], fac2 : list[list[int]], to_mat : list[list[int]], rows : list[int], common : list[int], columns : list[int]) -> None:

    for row in rows:
        for com in common:
            if fac1[row][com] == 1:
                for col in columns:
                    to_mat[row][col] ^= fac2[com][col]


def rank_matrix_positions(matrix : list[list[int]], rows : list[int], columns : list[int]) -> tuple[list[int], list[int]]:

    pos_selected : list[bool] = [False] * (max(max(rows), max(columns)) + 1)

    for row in rows:
        for col in columns:
            if not pos_selected[row] and not pos_selected[col] and matrix[row][col] == 1:
                pos_selected[row] = True
                pos_selected[col] = True
                for r in rows:
                    if not pos_selected[r] and matrix[r][col] == 1:
                        for c in columns:
                            matrix[r][c] ^= matrix[row][c]
                for c in columns:
                    if not pos_selected[c] and matrix[row][c] == 1:
                        for r in rows:
                            matrix[r][c] ^= matrix[r][col]

    return ([row for row in rows if pos_selected[row]], [col for col in columns if pos_selected[col]])


def matrix_inverse(to_be_inverted : list[list[int]], inverse : list[list[int]], rows : list[int], columns : list[int]) -> None:

    size = len(rows)
    insert_zero_matrix(inverse, columns, rows)
    for n in range(size):
        inverse[columns[n]][rows[n]] = 1

    for n in range(size):
        row = rows[n]
        col = columns[n]

        if (to_be_inverted[row][col] == 0):
            n_swap = next(nn for nn in range(n + 1, size) if to_be_inverted[rows[nn]][col] == 1)
            row_swap = rows[n_swap]
            for c in columns:
                (to_be_inverted[row][c], to_be_inverted[row_swap][c]) = (to_be_inverted[row_swap][c], to_be_inverted[row][c])
            col_swap = columns[n_swap]
            for r in rows:
                (inverse[col][r], inverse[col_swap][r]) = (inverse[col_swap][r], inverse[col][r])

        for n2 in range(size):
            if n2 != n and to_be_inverted[rows[n2]][col] == 1:
                to_row = rows[n2]
                for c in columns:
                    to_be_inverted[to_row][c] ^= to_be_inverted[row][c]
                to_col = columns[n2]
                for r in rows:
                    inverse[to_col][r] ^= inverse[col][r]
