# coding: utf-8

#from __future__ import division, unicode_literals

"""
Utilities for manipulating coordinates or list of coordinates, under periodic
boundary conditions or otherwise. Many of these are heavily vectorized in
numpy for performance.
"""

from six.moves import zip

__author__ = "Will Richards"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Will Richards"
__email__ = "wmdrichards@gmail.com"
__date__ = "Nov 27, 2011"

import numpy as np

from libc.stdlib cimport malloc, free
from libc.math cimport round, abs, sqrt, acos, M_PI
cimport numpy as np
cimport cython

#create images, 2d array of all length 3 combinations of [-1,0,1]
r = np.arange(-1, 2, dtype=np.float64)
arange = r[:, None] * np.array([1, 0, 0])[None, :]
brange = r[:, None] * np.array([0, 1, 0])[None, :]
crange = r[:, None] * np.array([0, 0, 1])[None, :]
images_t = arange[:, None, None] + brange[None, :, None] + \
    crange[None, None, :]
images = images_t.reshape((27, 3))

cdef np.float64_t[:, :] images_view = images

@cython.boundscheck(False)
@cython.wraparound(False)
cdef void dot_2d(np.float64_t[:, :] a, np.float64_t[:, :] b, np.float64_t* o) nogil:
    cdef int i, j, k, I, J, K
    I = a.shape[0]
    J = b.shape[1]
    K = a.shape[1]

    for j in range(J):
        for i in range(I):
            o[3 * i + j] = 0
        for k in range(K):
            for i in range(I):
                o[3 * i + j] += a[i, k] * b[k, j]

@cython.boundscheck(False)
@cython.wraparound(False)
cdef void dot_2d_mod(np.float64_t[:, :] a, np.float64_t[:, :] b, np.float64_t* o) nogil:
    cdef int i, j, k, I, J, K
    I = a.shape[0]
    J = b.shape[1]
    K = a.shape[1]

    for j in range(J):
        for i in range(I):
            o[3 * i + j] = 0
        for k in range(K):
            for i in range(I):
                o[3 * i + j] += a[i, k] % 1 * b[k, j]

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
def pbc_shortest_vectors(lattice, fcoords1, fcoords2):
    """
    Returns the shortest vectors between two lists of coordinates taking into
    account periodic boundary conditions and the lattice.

    Args:
        lattice: lattice to use
        fcoords1: First set of fractional coordinates. e.g., [0.5, 0.6, 0.7]
            or [[1.1, 1.2, 4.3], [0.5, 0.6, 0.7]]. It can be a single
            coord or any array of coords.
        fcoords2: Second set of fractional coordinates.

    Returns:
        array of displacement vectors from fcoords1 to fcoords2
        first index is fcoords1 index, second is fcoords2 index
    """

    #ensure correct shape
    #fcoords1, fcoords2 = np.atleast_2d(fcoords1, fcoords2)

    cdef np.ndarray[np.float64_t, ndim=2] lat = lattice._matrix

    cdef int i, j, k, l, I, J, L = 3
    I = len(fcoords1)
    J = len(fcoords2)

    cdef np.float64_t * cart_f1 = <np.float64_t *> malloc(3 * I * sizeof(np.float64_t))
    cdef np.float64_t * cart_f2 = <np.float64_t *> malloc(3 * J * sizeof(np.float64_t))
    cdef np.float64_t * cart_im = <np.float64_t *> malloc(81 * sizeof(np.float64_t))

    dot_2d_mod(fcoords1, lat, cart_f1)
    dot_2d_mod(fcoords2, lat, cart_f2)
    dot_2d(images_view, lat, cart_im)

    output = np.empty((I, J, 3))
    cdef np.float64_t[:, :, :] vectors = output
    cdef np.float64_t best, d

    for i in range(I):
        for j in range(J):
            best = 1e100
            for k in range(27):
                d = 0
                for l in range(3):
                    d += (cart_f2[3 * j + l] + cart_im[3 * k + l] - cart_f1[3 * i + l]) ** 2
                if d < best:
                    best = d
                    for l in range(3):
                        vectors[i, j, l] = cart_f2[3 * j + l] + cart_im[3 * k + l] - cart_f1[3 * i + l]

    free(<void *>cart_f1)
    free(<void *>cart_f2)
    free(<void *>cart_im)

    return output

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
def is_coord_subset_pbc(fcoords1, fcoords2, tol, mask):
    """
    Tests if all fractional coords in subset are contained in superset.
    Allows specification of a mask determining pairs that are not
    allowed to match to each other

    Args:
        subset, superset: List of fractional coords

    Returns:
        True if all of subset is in superset.
    """

    cdef np.float64_t[:, :] fc1 = fcoords1
    cdef np.float64_t[:, :] fc2 = fcoords2
    cdef np.float64_t[:] t = tol
    cdef np.int_t[:, :] m = np.array(mask, dtype=np.int)

    cdef int i, j, k, I, J
    cdef np.float64_t d
    cdef bint ok = False

    I = fc1.shape[0]
    J = fc2.shape[0]

    for i in range(I):
        for j in range(J):
            if m[i, j]:
                continue
            ok = True
            for k in range(3):
                d = fc1[i, k] - fc2[j, k]
                if abs(d - round(d)) > t[k]:
                    ok = False
                    break
            if ok:
                break
        if not ok:
            break
    return bool(ok)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
@cython.cdivision(True)
def lengths_and_angles(matrix):
    angles = np.empty(3, dtype=np.float64)
    lengths = np.empty(3, dtype=np.float64)

    cdef np.float64_t[:] a = angles
    cdef np.float64_t[:] l = lengths
    cdef np.float64_t[:, :] m = matrix
    cdef np.float64_t v

    cdef int i, j, k, n

    for i in range(3):
        l[i] = 0
        for j in range(3):
            l[i] += m[i, j] ** 2
        l[i] = sqrt(l[i])

    for i in range(3):
        j = (i + 1) % 3
        k = (i + 2) % 3
        v = 0
        for n in range(3):
            v += m[j, n] * m[k, n]
        v /= l[j] * l[k]
        a[i] = acos(max(min(v, 1), -1)) * 180 / M_PI

    return lengths, angles



