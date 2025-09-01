"""Wrapper for SciPy's linear_sum_assignment."""

from __future__ import annotations

import numpy as np
from scipy.optimize import linear_sum_assignment


class LinearAssignment:
    """Wrapper for SciPy's linear_sum_assignment.

    Args:
        costs: The cost matrix of the problem. cost[i,j] should be the
            cost of matching x[i] to y[j]. The cost matrix may be
            rectangular
        epsilon: Tolerance for determining if solution vector is < 0

    Attributes:
        min_cost: The minimum cost of the matching.
        solution: The matching of the rows to columns. i.e solution = [1, 2, 0]
            would match row 0 to column 1, row 1 to column 2 and row 2
            to column 0. Total cost would be c[0, 1] + c[1, 2] + c[2, 0].
    """

    def __init__(self, cost_matrix: np.ndarray) -> None:
        """Wrapper for SciPy's linear_sum_assignment."""
        self.orig_c = np.asarray(cost_matrix, dtype=float)

        n_rows, n_cols = self.orig_c.shape
        if n_rows > n_cols:
            raise ValueError("cost matrix must have at least as many columns as rows")

        row_ind, col_ind = linear_sum_assignment(self.orig_c)

        self.min_cost = self.orig_c[row_ind, col_ind].sum()
        self.solution = col_ind
