# cython: language_level=3
"""
This module contains the LAPJV algorithm to solve the Linear Assignment Problem.
"""

import numpy as np

cimport cython
cimport numpy as np
from libc.math cimport fabs
from libc.stdlib cimport free, malloc

np.import_array()

__author__ = "Will Richards"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Will Richards"
__email__ = "wrichards@mit.edu"
__date__ = "Jan 28, 2013"

class LinearAssignment:
    """
    This class finds the solution to the Linear Assignment Problem.
    It finds a minimum cost matching between two sets, given a cost
    matrix.

    This class is an implementation of the LAPJV algorithm described in:
    R. Jonker, A. Volgenant. A Shortest Augmenting Path Algorithm for
    Dense and Sparse Linear Assignment Problems. Computing 38, 325-340
    (1987)

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

    def __init__(self, costs: np.ndarray, epsilon: float=1e-13) -> None:
        self.orig_c = np.asarray(costs, dtype=np.float64, order="C")
        self.nx, self.ny = self.orig_c.shape
        self.n = self.ny

        self.epsilon = fabs(epsilon)

        # Check that cost matrix is square
        if self.nx > self.ny:
            raise ValueError("cost matrix must have at least as many columns as rows")

        if self.nx == self.ny:
            self.c = self.orig_c
        else:
            self.c = np.zeros((self.n, self.n), dtype=np.float64)
            self.c[:self.nx] = self.orig_c

        # Initialize solution vectors
        self._x = np.empty(self.n, dtype=np.int64)
        self._y = np.empty(self.n, dtype=np.int64)

        self.min_cost = compute(self.n, self.c, self._x, self._y, self.epsilon)
        self.solution = self._x[:self.nx]


@cython.boundscheck(False)
@cython.wraparound(False)
cdef np.float_t compute(int size, np.float_t[:, :] c, np.int64_t[:] x, np.int64_t[:] y, np.float_t eps) nogil:

    # Augment
    cdef int i, j, k, i1, j1, f, f0, cnt, low, up, z, last, nrr
    cdef int n = size
    cdef bint b
    cdef np.int64_t * col = <np.int64_t *> malloc(n * sizeof(np.int64_t))
    cdef np.int64_t * fre = <np.int64_t *> malloc(n * sizeof(np.int64_t))
    cdef np.int64_t * pred = <np.int64_t *> malloc(n * sizeof(np.int64_t))
    cdef np.float_t * v = <np.float_t *> malloc(n * sizeof(np.float_t))
    cdef np.float_t * d = <np.float_t *> malloc(n * sizeof(np.float_t))
    cdef np.float_t h, m, u1, u2, cost

    for i in range(n):
        x[i] = -1

    # Column reduction
    for j from n > j >= 0:
        col[j] = j
        h = c[0, j]
        i1 = 0
        for i in range(1, n):
            if c[i, j] < h:
                h = c[i, j]
                i1 = i
        v[j] = h
        if x[i1] == -1:
            x[i1] = j
            y[j] = i1
        else:
            # NOTE: in the paper it's x[i], but likely a typo
            if x[i1] > -1:
                x[i1] = -2 - x[i1]
            y[j] = -1

    # Reduction transfer
    f = -1
    for i in range(n):
        if x[i] == -1:
            f += 1
            fre[f] = i
        elif x[i] < -1:
            x[i] = -2 - x[i]
        else:
            j1 = x[i]
            m = 1e300
            for j in range(n):
                if j != j1:
                    if c[i, j] - v[j] < m:
                        m = c[i, j] - v[j]
            v[j1] = v[j1] - m

    # Augmenting row reduction
    for cnt in range(2):
        k = 0
        f0 = f
        f = -1
        # This step isn't strictly necessary, and
        # time is proportional to 1/eps in the worst case,
        # so break early by keeping track of nrr
        nrr = 0
        while k <= f0:
            nrr += 1
            i = fre[k]
            k += 1
            u1 = c[i, 0] - v[0]
            j1 = 0
            u2 = 1e300
            for j in range(1, n):
                h = c[i, j] - v[j]
                if h < u2:
                    if h >= u1:
                        u2 = h
                        j2 = j
                    else:
                        u2 = u1
                        u1 = h
                        j2 = j1
                        j1 = j
            i1 = y[j1]
            if u1 + eps < u2 and nrr < n * k:
                v[j1] = v[j1] - u2 + u1
            elif i1 > -1 and nrr < n * k:
                j1 = j2
                i1 = y[j1]
            if i1 > -1:
                if u1 + eps < u2 and nrr < n * k:
                    k -= 1
                    fre[k] = i1
                else:
                    f += 1
                    fre[f] = i1
            x[i] = j1
            y[j1] = i

    # Augmentation
    f0 = f
    for f in range(f0 + 1):
        i1 = fre[f]
        low = 0
        up = 0
        for j in range(n):
            d[j] = c[i1, j] - v[j]
            pred[j] = i1
        while True:
            # The pascal code ends when a single augmentation is found
            # really we need to get back to the for f in range(f0+1) loop
            b = False
            if up == low:
                last = low-1
                m = d[col[up]]
                up = up + 1
                for k in range(up, n):
                    j = col[k]
                    h = d[j]
                    if h <= m + eps:
                        if h + eps < m:
                            up = low
                            m = h
                        col[k] = col[up]
                        col[up] = j
                        up = up + 1

                for z in range(low, up):
                    j = col[z]
                    if y[j] == -1:
                        # augment
                        for k in range(last+1):
                            j1 = col[k]
                            v[j1] = v[j1] + d[j1] - m
                        while True:
                            i = pred[j]
                            y[j] = i
                            k = j
                            j = x[i]
                            x[i] = k
                            if i == i1:
                                b = True
                                break
                        break
            if b:
                break
            j1 = col[low]
            low = low + 1
            i = y[j1]
            u1 = c[i, j1] - v[j1] - m
            for k in range(up, n):
                j = col[k]
                h = c[i, j] - v[j] - u1
                if h + eps < d[j]:
                    d[j] = h
                    pred[j] = i
                    if fabs(h - m) < eps:
                        if y[j] == -1:
                            # Augment
                            for k in range(last+1):
                                j1 = col[k]
                                v[j1] = v[j1] + d[j1] - m
                            while True:
                                i = pred[j]
                                y[j] = i
                                k = j
                                j = x[i]
                                x[i] = k
                                if i == i1:
                                    b = True
                                    break
                            break
                        else:
                            col[k] = col[up]
                            col[up] = j
                            up = up + 1
            if b:
                break
    cost = 0
    for i in range(n):
        cost += c[i, x[i]]

    free(<void *>col)
    free(<void *>fre)
    free(<void *>pred)
    free(<void *>v)
    free(<void *>d)

    return cost
