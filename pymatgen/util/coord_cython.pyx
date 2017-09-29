# coding: utf-8

#from __future__ import division, unicode_literals

"""
Utilities for manipulating coordinates or list of coordinates, under periodic
boundary conditions or otherwise.
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
from libc.math cimport round, fabs, sqrt, acos, M_PI, cos, sin
cimport numpy as np
cimport cython

#create images, 2d array of all length 3 combinations of [-1,0,1]
r = np.arange(-1, 2, dtype=np.float_)
arange = r[:, None] * np.array([1, 0, 0])[None, :]
brange = r[:, None] * np.array([0, 1, 0])[None, :]
crange = r[:, None] * np.array([0, 0, 1])[None, :]
images_t = arange[:, None, None] + brange[None, :, None] + \
    crange[None, None, :]
images = images_t.reshape((27, 3))

cdef np.float_t[:, ::1] images_view = images

@cython.boundscheck(False)
@cython.wraparound(False)
cdef void dot_2d(np.float_t[:, ::1] a, np.float_t[:, ::1] b, np.float_t[:, ::1] o) nogil:
    cdef int i, j, k, I, J, K
    I = a.shape[0]
    J = b.shape[1]
    K = a.shape[1]

    for j in range(J):
        for i in range(I):
            o[i, j] = 0
        for k in range(K):
            for i in range(I):
                o[i, j] += a[i, k] * b[k, j]

@cython.boundscheck(False)
@cython.wraparound(False)
cdef void dot_2d_mod(np.float_t[:, ::1] a, np.float_t[:, ::1] b, np.float_t[:, ::1] o) nogil:
    cdef int i, j, k, I, J, K
    I = a.shape[0]
    J = b.shape[1]
    K = a.shape[1]

    for j in range(J):
        for i in range(I):
            o[i, j] = 0
        for k in range(K):
            for i in range(I):
                o[i, j] += a[i, k] % 1 * b[k, j]

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
def pbc_shortest_vectors(lattice, fcoords1, fcoords2, mask=None, return_d2=False, lll_frac_tol=None):
    """
    Returns the shortest vectors between two lists of coordinates taking into
    account periodic boundary conditions and the lattice.

    Args:
        lattice: lattice to use
        fcoords1: First set of fractional coordinates. e.g., [0.5, 0.6, 0.7]
            or [[1.1, 1.2, 4.3], [0.5, 0.6, 0.7]]. Must be np.float_
        fcoords2: Second set of fractional coordinates.
        mask (int_ array): Mask of matches that are not allowed.
            i.e. if mask[1,2] == True, then subset[1] cannot be matched
            to superset[2]
        lll_frac_tol (float_ array of length 3): Fractional tolerance (per LLL lattice vector) over which
            the calculation of minimum vectors will be skipped.
            Can speed up calculation considerably for large structures.

    Returns:
        array of displacement vectors from fcoords1 to fcoords2
        first index is fcoords1 index, second is fcoords2 index
    """

    #ensure correct shape
    fcoords1, fcoords2 = np.atleast_2d(fcoords1, fcoords2)


    fcoords1 = lattice.get_lll_frac_coords(fcoords1)
    fcoords2 = lattice.get_lll_frac_coords(fcoords2)

    cdef np.float_t[:, ::1] lat = np.array(lattice.lll_matrix, dtype=np.float_, copy=False, order='C')

    cdef int i, j, k, l, I, J, bestK

    I = len(fcoords1)
    J = len(fcoords2)

    cdef np.float_t[:, ::1] fc1 = fcoords1
    cdef np.float_t[:, ::1] fc2 = fcoords2

    cdef np.float_t[:, ::1] cart_f1 = <np.float_t[:I, :3]> malloc(3 * I * sizeof(np.float_t))
    cdef np.float_t[:, ::1] cart_f2 = <np.float_t[:J, :3]> malloc(3 * J * sizeof(np.float_t))
    cdef np.float_t[:, ::1] cart_im = <np.float_t[:27, :3]> malloc(81 * sizeof(np.float_t))

    cdef bint has_mask = mask is not None
    cdef np.int_t[:, :] m
    if has_mask:
        m = np.array(mask, dtype=np.int_, copy=False, order='C')

    cdef bint has_ftol = (lll_frac_tol is not None)
    cdef np.float_t[:] ftol
    if has_ftol:
        ftol = np.array(lll_frac_tol, dtype=np.float_, order='C', copy=False)


    dot_2d_mod(fc1, lat, cart_f1)
    dot_2d_mod(fc2, lat, cart_f2)
    dot_2d(images_view, lat, cart_im)

    vectors = np.empty((I, J, 3))
    d2 = np.empty((I, J))
    cdef np.float_t[:, :, ::1] vs = vectors
    cdef np.float_t[:, ::1] ds = d2
    cdef np.float_t best, d, inc_d, da, db, dc, fdist
    cdef bint within_frac = True
    cdef np.float_t[:] pre_im = <np.float_t[:3]> malloc(3 * sizeof(np.float_t))

    for i in range(I):
        for j in range(J):
            within_frac = False
            if (not has_mask) or (m[i, j] == 0):
                within_frac = True
                if has_ftol:
                    for l in range(3):
                        fdist = fc2[j, l] - fc1[i, l]
                        if fabs(fdist - round(fdist)) > ftol[l]:
                            within_frac = False
                            break
                if within_frac:
                    for l in range(3):
                        pre_im[l] = cart_f2[j, l] - cart_f1[i, l]
                    best = 1e100
                    for k in range(27):
                        # compilers have a hard time unrolling this
                        da = pre_im[0] + cart_im[k, 0]
                        db = pre_im[1] + cart_im[k, 1]
                        dc = pre_im[2] + cart_im[k, 2]
                        d = da * da + db * db + dc * dc
                        if d < best:
                            best = d
                            bestk = k
                    ds[i, j] = best
                    for l in range(3):
                        vs[i, j, l] = pre_im[l] + cart_im[bestk, l]
            if not within_frac:
                ds[i, j] = 1e20
                for l in range(3):
                    vs[i, j, l] = 1e20

    free(&cart_f1[0,0])
    free(&cart_f2[0,0])
    free(&cart_im[0,0])
    free(&pre_im[0])

    if return_d2:
        return vectors, d2
    else:
        return vectors

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
def is_coord_subset_pbc(subset, superset, atol, mask):
    """
    Tests if all fractional coords in subset are contained in superset.
    Allows specification of a mask determining pairs that are not
    allowed to match to each other

    Args:
        subset, superset: List of fractional coords

    Returns:
        True if all of subset is in superset.
    """

    cdef np.float_t[:, :] fc1 = subset
    cdef np.float_t[:, :] fc2 = superset
    cdef np.float_t[:] t = atol
    cdef np.int_t[:, :] m = np.array(mask, dtype=np.int_, copy=False, order='C')

    cdef int i, j, k, I, J
    cdef np.float_t d
    cdef bint ok

    I = fc1.shape[0]
    J = fc2.shape[0]

    for i in range(I):
        ok = False
        for j in range(J):
            if m[i, j]:
                continue
            ok = True
            for k in range(3):
                d = fc1[i, k] - fc2[j, k]
                if fabs(d - round(d)) > t[k]:
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
def coord_list_mapping_pbc(subset, superset, atol=1e-8):
    """
    Gives the index mapping from a subset to a superset.
    Superset cannot contain duplicate matching rows

    Args:
        subset, superset: List of frac_coords

    Returns:
        list of indices such that superset[indices] = subset
    """
    inds = np.zeros(len(subset), dtype=np.int) - 1
    subset = np.atleast_2d(subset)
    superset = np.atleast_2d(superset)

    cdef np.float_t[:, :] fc1 = subset
    cdef np.float_t[:, :] fc2 = superset
    cdef np.float_t[:] t = atol
    cdef np.int_t[:] c_inds = inds
    cdef np.float_t d
    cdef bint ok_inner, ok_outer

    I = fc1.shape[0]
    J = fc2.shape[0]

    for i in range(I):
        ok_outer = False
        for j in range(J):
            ok_inner = True
            for k in range(3):
                d = fc1[i, k] - fc2[j, k]
                if fabs(d - round(d)) > t[k]:
                    ok_inner = False
                    break
            if ok_inner:
                if c_inds[i] >= 0:
                    raise ValueError("Something wrong with the inputs, likely duplicates "
                             "in superset")
                c_inds[i] = j
                ok_outer = True
                # we don't break here so we can check for duplicates in superset
        if not ok_outer:
            break

    if not ok_outer:
        raise ValueError("subset is not a subset of superset")

    return inds
