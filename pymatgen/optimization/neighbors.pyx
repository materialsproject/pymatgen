# cython: boundscheck=False
# cython: wraparound=False
# cython: nonecheck=False
# cython: cdivision=False
# cython: profile=True
# cython: language_level=3
# distutils: language = c
# distutils: define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION

import numpy as np

cimport numpy as np
from libc.math cimport ceil, floor, pi, sqrt
from libc.stdlib cimport free, malloc, realloc
from libc.string cimport memset


cdef void *safe_malloc(size_t size) except? NULL:
    """Raise memory error if malloc fails"""
    if size == 0:
        return NULL
    cdef void *ptr = malloc(size)
    if ptr == NULL:
        raise MemoryError("Memory allocation of %s bytes failed!" % size)
    return ptr


cdef void *safe_realloc(void *ptr_orig, size_t size) except? NULL:
    """Raise memory error if realloc fails"""
    if size == 0:
        free(ptr_orig)
        return NULL
    cdef void *ptr = realloc(ptr_orig, size)
    if ptr == NULL:
        raise MemoryError("Realloc memory of %s bytes failed!" % size)
    return ptr


def find_points_in_spheres(double[:, ::1] all_coords, double[:, ::1] center_coords,
                           float r, long[:] pbc, double[:, ::1] lattice,
                           double tol=1e-8, float min_r=1.0):
    """
    For each point in `center_coords`, get all the neighboring points in `all_coords` that are within the
    cutoff radius `r`. All the coordinates should be in Cartesian.

    Args:
        all_coords: (np.ndarray[double, dim=2]) all available points. When periodic boundary is considered,
            this is all the points in the lattice.
        center_coords: (np.ndarray[double, dim=2]) all centering points
        r: (float) cutoff radius
        pbc: (list of bool) whether to set periodic boundaries
        lattice: (np.ndarray[double, dim=2]) 3x3 lattice matrix
        tol: (float) numerical tolerance
        min_r: (float) minimal cutoff to calculate the neighbor list
            directly. If the cutoff is less than this value, the algorithm
            will calculate neighbor list using min_r as cutoff and discard
            those that have larger distances.
    Returns:
        index1 (n, ), index2 (n, ), offset_vectors (n, 3), distances (n, ). index1 of center_coords, and index2 of all_coords that form the neighbor pair
            offset_vectors are the periodic image offsets for the all_coords.
    """
    if r < min_r:
        findex1, findex2, foffset_vectors, fdistances = find_points_in_spheres(
            all_coords=all_coords, center_coords=center_coords,
            r=min_r + tol, pbc=pbc, lattice=lattice, tol=tol, min_r=min_r)
        mask = fdistances <= r
        return findex1[mask], findex2[mask], foffset_vectors[mask], fdistances[
            mask]

    cdef int i, j, k, l, m, n
    cdef double maxr[3]
    # valid boundary, that is the minimum in center_coords - r
    cdef double valid_min[3]
    cdef double valid_max[3]
    cdef double ledge
    if r < 0.1:
        ledge = 0.1
    else:
        ledge = r
    max_and_min(center_coords, valid_max, valid_min)
    for i in range(3):
        valid_max[i] = valid_max[i] + r + tol
        valid_min[i] = valid_min[i] - r - tol

    # Process pbc
    cdef int n_center = center_coords.shape[0]
    cdef int n_total = all_coords.shape[0]
    cdef double [:, ::1] frac_coords =  <double[:n_center, :3]> safe_malloc(n_center * 3 * sizeof(double))
    cdef long[3] max_bounds = [1, 1, 1]
    cdef long[3] min_bounds = [0, 0, 0]
    cdef long nlattice = 1
    cdef long count = 0
    cdef double[:, ::1] all_fcoords = <double[:n_total, :3]> safe_malloc(n_total * 3 * sizeof(double))
    cdef double[:, ::1] coords_in_cell = <double[:n_total, :3]> safe_malloc(n_total * 3 * sizeof(double))
    cdef double[:, ::1] offset_correction = <double[:n_total, :3]> safe_malloc(n_total * 3 * sizeof(double))

    get_frac_coords(lattice, all_coords, offset_correction)
    for i in range(n_total):
        for j in range(3):
            if pbc[j]:
                # only wrap atoms when this dimension is PBC
                all_fcoords[i, j] = offset_correction[i, j] % 1
                offset_correction[i, j] = offset_correction[i, j] - all_fcoords[i, j]
            else:
                all_fcoords[i, j] = offset_correction[i, j]
                offset_correction[i, j] = 0
    get_max_r(lattice, maxr, r)
    # Get fractional coordinates of center points
    get_frac_coords(lattice, center_coords, frac_coords)
    get_bounds(frac_coords, maxr, pbc, max_bounds, min_bounds)
    for i in range(3):
        nlattice *= (max_bounds[i] - min_bounds[i])
    matmul(all_fcoords, lattice, coords_in_cell)

    # Get translated images, coordinates and indices
    cdef long natoms = n_total
    cdef double *offsets_p_temp = <double*> safe_malloc(natoms * 3 * sizeof(double))
    cdef double *expanded_coords_p_temp = <double *> safe_malloc(natoms * 3 * sizeof(double))
    cdef long *indices_p_temp = <long*> safe_malloc(natoms * sizeof(long))
    cdef double coord_temp[3]

    count = 0
    for i in range(min_bounds[0], max_bounds[0]):
        for j in range(min_bounds[1], max_bounds[1]):
            for k in range(min_bounds[2], max_bounds[2]):
                for l in range(n_total):
                    for m in range(3):
                        coord_temp[m] = <double>i * lattice[0, m] + <double>j * lattice[1, m] + \
                            <double>k * lattice[2, m] + coords_in_cell[l, m]
                    if (coord_temp[0] > valid_min[0]) & (coord_temp[0] < valid_max[0]) & \
                        (coord_temp[1] > valid_min[1]) & (coord_temp[1] < valid_max[1]) & \
                        (coord_temp[2] > valid_min[2]) & (coord_temp[2] < valid_max[2]):
                        offsets_p_temp[3*count] = i
                        offsets_p_temp[3*count+1] = j
                        offsets_p_temp[3*count+2] = k
                        indices_p_temp[count] = l
                        expanded_coords_p_temp[3*count] = coord_temp[0]
                        expanded_coords_p_temp[3*count+1] = coord_temp[1]
                        expanded_coords_p_temp[3*count+2] = coord_temp[2]
                        count += 1
                        if count >= natoms:  # exceeding current memory
                            natoms += natoms
                            offsets_p_temp = <double*>safe_realloc(offsets_p_temp, natoms * 3 * sizeof(double))
                            expanded_coords_p_temp = <double*>safe_realloc(expanded_coords_p_temp, natoms * 3 * sizeof(double))
                            indices_p_temp = <long*>safe_realloc(indices_p_temp, natoms * sizeof(long))

    # if no valid neighbors were found return empty
    if count == 0:
        free(&frac_coords[0, 0])
        free(&all_fcoords[0, 0])
        free(&coords_in_cell[0, 0])
        free(&offset_correction[0, 0])
        free(offsets_p_temp)
        free(expanded_coords_p_temp)
        free(indices_p_temp)
        return (np.array([], dtype=int), np.array([], dtype=int),
            np.array([[], [], []], dtype=float).T, np.array([], dtype=float))

    # Delete those beyond (min_center_coords - r, max_center_coords + r)
    cdef double *offsets_p = <double*> safe_realloc(offsets_p_temp, count * 3 * sizeof(double))
    cdef double *expanded_coords_p = <double*> safe_realloc(expanded_coords_p_temp, count * 3 * sizeof(double))
    cdef long *indices_p = <long*> safe_realloc(indices_p_temp, count * sizeof(long))

    cdef double[:, ::1] offsets = <double[:count, :3]> offsets_p
    cdef double[:, ::1] expanded_coords = <double[:count, :3]> expanded_coords_p
    cdef long[::1] indices = <long[:count]> indices_p

    # Construct linked cell list
    natoms = count
    cdef long ncube[3]
    cdef long[:, ::1] all_indices3 = <long[:natoms, :3]> safe_malloc(natoms * 3 * sizeof(long))
    cdef long[::1] all_indices1 = <long[:natoms]> safe_malloc(natoms * sizeof(long))

    for i in range(3):
        ncube[i] = <long>(ceil((valid_max[i] - valid_min[i]) / ledge))

    compute_cube_index(expanded_coords, valid_min, ledge, all_indices3)
    three_to_one(all_indices3, ncube[1], ncube[2], all_indices1)

    cdef long nb_cubes = ncube[0] * ncube[1] * ncube[2]
    cdef long *head = <long*> safe_malloc(nb_cubes*sizeof(long))
    cdef long *atom_indices = <long*> safe_malloc(natoms*sizeof(long))
    memset(<void*>head, -1, nb_cubes*sizeof(long))
    memset(<void*>atom_indices, -1, natoms*sizeof(long))

    cdef long[:, ::1] neighbor_map = <long[:nb_cubes, :27]> safe_malloc(nb_cubes * 27 * sizeof(long))
    get_cube_neighbors(ncube, neighbor_map)
    for i in range(natoms):
        atom_indices[i] = head[all_indices1[i]]
        head[all_indices1[i]] = i

    # Get center atoms' cube indices
    cdef long[:, ::1] center_indices3 = <long[:n_center, :3]> safe_malloc(n_center*3*sizeof(long))
    cdef long *center_indices1 = <long*> safe_malloc(n_center*sizeof(long))
    compute_cube_index(center_coords, valid_min, ledge, center_indices3)
    three_to_one(center_indices3, ncube[1], ncube[2], <long[:n_center]>center_indices1)

    # it works by allocate incrementally n more memory locations, I found it 3x faster to do so
    # compared to using vectors in cpp
    n = 10000
    cdef long *index_1 = <long*> safe_malloc(n*sizeof(long))
    cdef long *index_2 = <long*> safe_malloc(n*sizeof(long))
    cdef double *offset_final = <double*> safe_malloc(3*n*sizeof(double))
    cdef double *distances = <double*> safe_malloc(n*sizeof(double))
    cdef long[:] ncube_indices_map
    cdef long cube_index_temp
    cdef long link_index
    cdef double d_temp2
    cdef double r2 = r * r

    count = 0
    for i in range(n_center):
        ncube_indices_temp = neighbor_map[center_indices1[i]]
        for j in ncube_indices_temp:
            if j == -1:
                continue
            cube_index_temp = j
            link_index = head[cube_index_temp]
            while link_index != -1:
                d_temp2 = distance2(expanded_coords, center_coords, link_index, i, 3)
                if d_temp2 < r2 + tol:
                    index_1[count] = i
                    index_2[count] = indices[link_index]
                    offset_final[3*count] = offsets[link_index, 0] - offset_correction[indices[link_index], 0]
                    offset_final[3*count + 1] = offsets[link_index, 1] - offset_correction[indices[link_index], 1]
                    offset_final[3*count + 2] = offsets[link_index, 2] - offset_correction[indices[link_index], 2]
                    distances[count] = sqrt(d_temp2)

                    count += 1
                    # increasing the memory size
                    if count >= n:
                        n += n  # double the size
                        index_1 = <long*>safe_realloc(index_1, n * sizeof(long))
                        index_2 = <long*>safe_realloc(index_2, n*sizeof(long))
                        offset_final = <double*>safe_realloc(offset_final, 3*n*sizeof(double))
                        distances = <double*>safe_realloc(distances, n*sizeof(double))
                link_index = atom_indices[link_index]

    if count == 0:
        free(index_1)
        free(index_2)
        free(offset_final)
        free(distances)
        free(&offset_correction[0, 0])
        free(&frac_coords[0, 0])
        free(&all_fcoords[0, 0])
        free(&coords_in_cell[0, 0])
        free(&center_indices3[0, 0])
        free(&neighbor_map[0, 0])
        free(offsets_p)
        free(expanded_coords_p)
        free(indices_p)
        free(&all_indices3[0, 0])
        free(&all_indices1[0])
        free(center_indices1)
        free(head)
        free(atom_indices)

        return np.array([], dtype=int), np.array([], dtype=int), np.array([[], [], []], dtype=float).T, np.array([], dtype=float)

    index_1 = <long*>safe_realloc(index_1, count * sizeof(long))
    index_2 = <long*>safe_realloc(index_2, count*sizeof(long))
    offset_final = <double*>safe_realloc(offset_final, 3*count*sizeof(double))
    distances = <double*>safe_realloc(distances, count*sizeof(double))

    # convert to python objects
    py_index_1 = np.array(<long[:count]>index_1)
    py_index_2 = np.array(<long[:count]>index_2)
    py_offsets = np.array(<double[:count, :3]>offset_final)
    py_distances = np.array(<double[:count]>distances)

    # free allocated memories
    free(index_1)
    free(index_2)
    free(offset_final)
    free(distances)
    free(&offset_correction[0, 0])
    free(&frac_coords[0, 0])
    free(&all_fcoords[0, 0])
    free(&coords_in_cell[0, 0])
    free(&center_indices3[0, 0])
    free(&neighbor_map[0, 0])
    free(offsets_p)
    free(expanded_coords_p)
    free(indices_p)
    free(&all_indices3[0, 0])
    free(&all_indices1[0])
    free(center_indices1)
    free(head)
    free(atom_indices)
    return py_index_1, py_index_2, py_offsets, py_distances

cdef double distance2(double[:, ::1] m1, double[:, ::1] m2, long index1, long index2, long size):
    """
    Faster way to compute the distance squared by not using slice but providing indices in each matrix

    """
    cdef double s = 0
    cdef long i
    for i in range(size):
        s += (m1[index1, i] - m2[index2, i]) * (m1[index1, i] - m2[index2, i])
    return s

cdef void get_cube_neighbors(long [:] ncube, long[:, ::1] neighbor_map):
    """
    Get {cube_index: cube_neighbor_indices} map
    """
    cdef long ncubes = ncube[0] * ncube[1] * ncube[2]
    memset(<void*>&neighbor_map[0, 0], -1, neighbor_map.shape[0] * 27 * sizeof(long))
    cdef long[::1] counts = <long[:ncubes]> safe_malloc(ncubes * sizeof(long))
    cdef long[:, ::1] cube_indices_3d = <long[:ncubes, :3]> safe_malloc(ncubes*3*sizeof(long))
    cdef long[::1] cube_indices_1d = <long[:ncubes]> safe_malloc(ncubes*sizeof(long))
    cdef long count = 0
    cdef int i, j, k

    for i in range(ncubes):
        counts[i] = 0

    for i in range(ncube[0]):
        for j in range(ncube[1]):
            for k in range(ncube[2]):
                cube_indices_3d[count, 0] = i
                cube_indices_3d[count, 1] = j
                cube_indices_3d[count, 2] = k
                count += 1
    three_to_one(cube_indices_3d, ncube[1], ncube[2], cube_indices_1d)
    cdef long[:, ::1] index3 = <long[:1, :3]> safe_malloc(3*sizeof(long))
    cdef long index1[1]
    cdef long[:, ::1] ovectors = compute_offset_vectors(1)

    for i in range(ncubes):
        for j in range(27):
            index3[0, 0] = ovectors[j, 0] + cube_indices_3d[i, 0]
            index3[0, 1] = ovectors[j, 1] + cube_indices_3d[i, 1]
            index3[0, 2] = ovectors[j, 2] + cube_indices_3d[i, 2]
            if (index3[0, 0] < ncube[0]) & (index3[0, 0] >= 0) & (index3[0, 1] < ncube[1]) & \
                (index3[0, 1] >= 0) & (index3[0, 2] < ncube[2]) & (index3[0, 2] >= 0):
                three_to_one(index3, ncube[1], ncube[2], index1)
                neighbor_map[i, counts[i]] = index1[0]
                counts[i] += 1

    free(&cube_indices_3d[0, 0])
    free(&cube_indices_1d[0])
    free(&index3[0, 0])
    free(&counts[0])

cdef void get_bounds(double[:, ::1] frac_coords, double[::1] maxr, long[:] pbc, long[:] max_bounds, long[:] min_bounds):
    """Given the fractional coordinates and the number of repeation needed in each direction, maxr,
    compute the translational bounds in each dimension
    """
    cdef double max_fcoords[3]
    cdef double min_fcoords[3]
    max_and_min(frac_coords, max_fcoords, min_fcoords)
    cdef int i
    for i in range(3):
        min_bounds[i] = 0
        max_bounds[i] = 1

    for i in range(3):
        if pbc[i]:
            min_bounds[i] = <long>(floor(min_fcoords[i] - maxr[i] - 1e-8))
            max_bounds[i] = <long>(ceil(max_fcoords[i] + maxr[i] + 1e-8))

cdef void get_frac_coords(double[:, ::1] lattice, double[:, ::1] cart_coords, double [:, ::1] frac_coords):
    """
    Compute the fractional coordinates
    """
    cdef double[:, ::1] inv_lattice = np.empty((3, 3), dtype=np.float64)
    matrix_inv(lattice, inv_lattice)
    matmul(cart_coords, inv_lattice, frac_coords)

cdef void matmul(double [:, ::1] m1, double [:, ::1] m2, double [:, ::1] out):
    """
    Matrix multiplication
    """
    cdef int m = m1.shape[0], n = m1.shape[1], l = m2.shape[1]
    cdef int i, j, k
    for i in range(m):
        for j in range(l):
            out[i, j] = 0
            for k in range(n):
                out[i, j] += m1[i, k] * m2[k, j]

cdef void matrix_inv(double[:, ::1] matrix, double[:, ::1] inv):
    """
    Matrix inversion
    """
    cdef double det = matrix_det(matrix)
    cdef int i, j
    for i in range(3):
        for j in range(3):
            inv[i, j] = (matrix[(j+1)%3, (i+1)%3] * matrix[(j+2)%3, (i+2)%3] -
                matrix[(j+2)%3, (i+1)%3] * matrix[(j+1)%3, (i+2)%3]) / det

cdef double matrix_det(double[:, ::1] matrix):
    """
    Matrix determinant
    """
    return matrix[0, 0] * (matrix[1, 1] * matrix[2, 2] - matrix[1, 2] * matrix[2, 1]) + \
        matrix[0, 1] * (matrix[1, 2] * matrix[2, 0] - matrix[1, 0] * matrix[2, 2]) + \
            matrix[0, 2] * (matrix[1, 0] * matrix[2, 1] - matrix[1, 1] * matrix[2, 0])

cdef void get_max_r(double[:, ::1] lattice, double [::1] maxr, double r):
    """
    Get maximum repetition in each directions
    """
    cdef double [:, ::1] reciprocal_lattice = np.eye(3, dtype=np.float64)
    cdef int i
    cdef int n = lattice.shape[1]
    cdef double recp_len
    get_reciprocal_lattice(lattice, reciprocal_lattice)
    for i in range(n):
        recp_len = norm(reciprocal_lattice[i, :])
        maxr[i] = ceil((r + 0.15) * recp_len / (2 * pi))

cdef void get_reciprocal_lattice(double[:, ::1] lattice, double[:, ::1] reciprocal):
    """
    Compute the reciprocal lattice
    """
    cdef int i
    for i in range(3):
        recip_component(lattice[i, :], lattice[(i+1)%3, :], lattice[(i+2)%3, :], reciprocal[i, :])

cdef void recip_component(double[::1] a1, double[::1] a2, double[::1] a3, double[::1] out):
    """
    Compute the reciprocal lattice vector
    """
    cdef double ai_cross_aj[3]
    cdef double prod
    cdef int i
    cross(a2, a3, ai_cross_aj)
    prod = inner(a1, ai_cross_aj)
    for i in range(a1.shape[0]):
        out[i] = 2 * pi * ai_cross_aj[i] / prod

cdef double inner(double[:] x, double[:] y):
    """
    Compute inner product
    """
    cdef double sum = 0
    cdef int i, n = x.shape[0]
    for i in range(n):
        sum += x[i] * y[i]
    return sum

cdef void cross(double[:] x, double[:] y, double[:] out) nogil:
    """
    Cross product of vector x and y, output in out
    """
    out[0] = x[1] * y[2] - x[2] * y[1]
    out[1] = x[2] * y[0] - x[0] * y[2]
    out[2] = x[0] * y[1] - x[1] * y[0]

cdef double norm(double[::1] vec) nogil:
    """
    Vector norm
    """
    cdef int i, n = vec.shape[0]
    cdef double sum = 0
    for i in range(n):
        sum += vec[i] * vec[i]
    return sqrt(sum)

cdef max_and_min(double[:, ::1] coords, double[::1] max_coords, double[::1] min_coords):
    """
    Compute the min and max of coords
    """
    cdef int i, j
    cdef int M = coords.shape[0]
    cdef int N = coords.shape[1]
    for i in range(N):
        max_coords[i] = coords[0, i]
        min_coords[i] = coords[0, i]
    for i in range(M):
        for j in range(N):
            if coords[i, j] >= max_coords[j]:
                max_coords[j] = coords[i, j]
            if coords[i, j] <= min_coords[j]:
                min_coords[j] = coords[i, j]

cdef void compute_cube_index(double[:, ::1] coords, double [:] global_min, double radius, long[:, ::1] return_indice) nogil:
    cdef int i, j
    for i in range(coords.shape[0]):
        for j in range(coords.shape[1]):
            return_indice[i, j] = <long>(floor((coords[i, j] - global_min[j] + 1e-8) / radius))


cdef void three_to_one(long[:, ::1] label3d, long ny, long nz, long[::1] label1d) nogil:
    """
    3D vector representation to 1D
    """
    cdef int i
    cdef int n = label3d.shape[0]
    for i in range(n):
        label1d[i] = label3d[i, 0] * ny * nz + label3d[i, 1] * nz + label3d[i, 2]


def compute_offset_vectors(long n):
    cdef long i, j, k
    cdef long v[3]
    cdef double center[8][3] # center vertices coords
    cdef int ind
    cdef bint is_within
    cdef long ntotal = (2*n+1) * (2*n+1) * (2*n+1)
    cdef long *ovectors = <long*> safe_malloc(ntotal*3*sizeof(long))
    cdef long count = 0
    for i in range(2):
        for j in range(2):
            for k in range(2):
                ind = i * 4 + j * 2 + k
                center[ind][0] = i - 0.5
                center[ind][1] = j - 0.5
                center[ind][2] = k - 0.5

    cdef double off[8][3] # offseted vertices
    for i in range(-n, n + 1):
        for j in range(-n, n + 1):
            for k in range(-n, n + 1):
                offset_cube(center, i, j, k, off)
                if distance_vertices(center, off, n):
                    ovectors[3*count] = i
                    ovectors[3*count + 1] = j
                    ovectors[3*count + 2] = k
                    count += 1
    ovectors = <long*>safe_realloc(ovectors, count*3*sizeof(long))
    array = np.array(<long[:count, :3]>ovectors)
    free(ovectors)
    return array

cdef bint distance_vertices(double center[8][3], double off[8][3], double r):
    cdef int i, j
    cdef double d2
    cdef double r2 = r * r
    for i in range(8):
        for j in range(8):
            d2 = (center[i][0] - off[j][0]) * (center[i][0] - off[j][0]) + (center[i][1] - off[j][1]) * (center[i][1] - off[j][1]) + \
                (center[i][2] - off[j][2]) * (center[i][2] - off[j][2])
            if d2 <= r2:
                return 1
    return 0

cdef void offset_cube(double center[8][3], long n, long m, long l, double (&offseted)[8][3]):
    cdef int i, j, k
    for i in range(2):
        for j in range(2):
            for k in range(2):
                ind = i * 4 + j * 2 + k
                offseted[ind][0] = center[ind][0] + n
                offseted[ind][1] = center[ind][1] + m
                offseted[ind][2] = center[ind][2] + l
