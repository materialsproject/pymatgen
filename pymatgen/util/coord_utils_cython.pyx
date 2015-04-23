# coding: utf-8

#from __future__ import division, unicode_literals

"""
Utilities for manipulating coordinates or list of coordinates, under periodic
boundary conditions or otherwise. Many of these are heavily vectorized in
numpy for performance.
"""

from six.moves import zip

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Nov 27, 2011"

import numpy as np

from libc.stdlib cimport malloc, free
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

    cdef int i, j, k, l, I, J, K = 27, L = 3
    I = fcoords1.shape[0]
    J = fcoords2.shape[0]

    cdef np.float64_t * cart_f1 = <np.float64_t *> malloc(3 * I * sizeof(np.float64_t*))
    cdef np.float64_t * cart_f2 = <np.float64_t *> malloc(3 * J * sizeof(np.float64_t*))
    cdef np.float64_t * cart_im = <np.float64_t *> malloc(81 * sizeof(np.float64_t*))

    dot_2d_mod(fcoords1, lat, cart_f1)
    dot_2d_mod(fcoords2, lat, cart_f2)
    dot_2d(images, lat, cart_im)

    output = np.empty((I, J, 3))
    cdef np.float64_t[:, :, :] vectors = output
    cdef np.float64_t best, d

    for i in range(I):
        for j in range(J):
            best = 1e100
            for k in range(K):
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
