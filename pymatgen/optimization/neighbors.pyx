# cython: language_level=3
# cython: initializedcheck=False
# cython: boundscheck=False
# cython: wraparound=False
# cython: nonecheck=False
# cython: cdivision=False
# cython: profile=False
# distutils: language = c

# isort: dont-add-imports

# Setting cdivision=True gives a small speed improvement, but the code is currently
# written based on Python division so using cdivision may result in missing neighbors
# in some off cases. See https://github.com/materialsproject/pymatgen/issues/2226

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
        raise MemoryError(f"Memory allocation of {size} bytes failed!")
    return ptr


cdef void *safe_realloc(void *ptr_orig, size_t size) except? NULL:
    """Raise memory error if realloc fails"""
    if size == 0:
        return NULL
    cdef void *ptr = realloc(ptr_orig, size)
    if ptr == NULL:
        raise MemoryError(f"Realloc memory of {size} bytes failed!")
    return ptr


def find_points_in_spheres(
        const double[:, ::1] all_coords,
        const double[:, ::1] center_coords,
        const double r,
        const long[::1] pbc,
        const double[:, ::1] lattice,
        const double tol=1e-8,
        const double min_r=1.0):
    """
    For each point in `center_coords`, get all the neighboring points in `all_coords`
    that are within the cutoff radius `r`. All the coordinates should be Cartesian.

    Args:
        all_coords: (np.ndarray[double, dim=2]) all available points.
            When periodic boundary is considered, this is all the points in the lattice.
        center_coords: (np.ndarray[double, dim=2]) all centering points
        r: (float) cutoff radius
        pbc: (np.ndarray[long, dim=1]) whether to set periodic boundaries
        lattice: (np.ndarray[double, dim=2]) 3x3 lattice matrix
        tol: (float) numerical tolerance
        min_r: (float) minimal cutoff to calculate the neighbor list
            directly. If the cutoff is less than this value, the algorithm
            will calculate neighbor list using min_r as cutoff and discard
            those that have larger distances.
    Returns:
        index1 (n, ), index2 (n, ), offset_vectors (n, 3), distances (n, ).
        index1 of center_coords, and index2 of all_coords that form the neighbor pair
        offset_vectors are the periodic image offsets for the all_coords.
    """
    if r < min_r:
        findex1, findex2, foffset_vectors, fdistances = find_points_in_spheres(
            all_coords=all_coords, center_coords=center_coords,
            r=min_r + tol, pbc=pbc, lattice=lattice, tol=tol, min_r=min_r)
        mask = fdistances <= r
        return findex1[mask], findex2[mask], foffset_vectors[mask], fdistances[
            mask]

    cdef:
        int i, j, k, l, m
        double[3] maxr
        # valid boundary, that is the minimum in center_coords - r
        double[3] valid_min
        double[3] valid_max
        double ledge

        int n_center = center_coords.shape[0]
        int n_total = all_coords.shape[0]
        long nlattice = 1

        long[3] max_bounds = [1, 1, 1]
        long[3] min_bounds = [0, 0, 0]
        double [:, ::1] frac_coords =  <double[:n_center, :3]> safe_malloc(
            n_center * 3 * sizeof(double)
        )
        double[:, ::1] all_fcoords = <double[:n_total, :3]> safe_malloc(
            n_total * 3 * sizeof(double)
        )
        double[:, ::1] coords_in_cell = <double[:n_total, :3]> safe_malloc(
            n_total * 3 * sizeof(double)
        )
        double[:, ::1] offset_correction = <double[:n_total, :3]> safe_malloc(
            n_total * 3 * sizeof(double)
        )
        double[3][3] inv_lattice_arr
        double[:, ::1] inv_lattice = inv_lattice_arr
        double[3][3] reciprocal_lattice_arr = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        double[:, ::1] reciprocal_lattice = reciprocal_lattice_arr

        int count = 0
        int natoms = n_total
        double *offsets_p_temp = <double*> safe_malloc(natoms * 3 * sizeof(double))
        double *expanded_coords_p_temp = <double*> safe_malloc(
            natoms * 3 * sizeof(double)
        )
        long *indices_p_temp = <long*> safe_malloc(natoms * sizeof(long))
        double coord_temp[3]
        long ncube[3]

        long[:, ::1] center_indices3 = <long[:n_center, :3]> safe_malloc(
            n_center*3*sizeof(long)
        )
        long[::1] center_indices1 = <long[:n_center]> safe_malloc(n_center*sizeof(long))

        int malloc_chunk = 10000  # size of memory chunks to re-allocate dynamically
        int failed_malloc = 0  # flag for failed reallocation within loops
        long *index_1 = <long*> safe_malloc(malloc_chunk*sizeof(long))
        long *index_2 = <long*> safe_malloc(malloc_chunk*sizeof(long))
        double *offset_final = <double*> safe_malloc(3*malloc_chunk*sizeof(double))
        double *distances = <double*> safe_malloc(malloc_chunk*sizeof(double))
        long cube_index_temp
        long link_index
        double d_temp2
        double r2 = r * r

    if r < 0.1:
        ledge = 0.1
    else:
        ledge = r
    max_and_min(center_coords, valid_max, valid_min)
    for i in range(3):
        valid_max[i] = valid_max[i] + r + tol
        valid_min[i] = valid_min[i] - r - tol

    # Process pbc
    get_frac_coords(lattice, inv_lattice, all_coords, offset_correction)
    for i in range(n_total):
        for j in range(3):
            if pbc[j]:
                # only wrap atoms when this dimension is PBC
                all_fcoords[i, j] = offset_correction[i, j] % 1
                offset_correction[i, j] = offset_correction[i, j] - all_fcoords[i, j]
            else:
                all_fcoords[i, j] = offset_correction[i, j]
                offset_correction[i, j] = 0

    # compute the reciprocal lattice in place
    get_reciprocal_lattice(lattice, reciprocal_lattice)
    get_max_r(reciprocal_lattice, maxr, r)

    # Get fractional coordinates of center points in place
    get_frac_coords(lattice, inv_lattice, center_coords, frac_coords)
    get_bounds(frac_coords, maxr, &pbc[0], max_bounds, min_bounds)

    for i in range(3):
        nlattice *= (max_bounds[i] - min_bounds[i])
    matmul(all_fcoords, lattice, coords_in_cell)

    # Get translated images, coordinates and indices
    for i in range(min_bounds[0], max_bounds[0]):
        for j in range(min_bounds[1], max_bounds[1]):
            for k in range(min_bounds[2], max_bounds[2]):
                for l in range(n_total):
                    for m in range(3):
                        coord_temp[m] = <double>i * lattice[0, m] + \
                                        <double>j * lattice[1, m] + \
                                        <double>k * lattice[2, m] + \
                                        coords_in_cell[l, m]
                    if (
                            (coord_temp[0] > valid_min[0]) &
                            (coord_temp[0] < valid_max[0]) &
                            (coord_temp[1] > valid_min[1]) &
                            (coord_temp[1] < valid_max[1]) &
                            (coord_temp[2] > valid_min[2]) &
                            (coord_temp[2] < valid_max[2])
                    ):
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
                            offsets_p_temp = <double*> realloc(
                                offsets_p_temp, natoms * 3 * sizeof(double)
                            )
                            expanded_coords_p_temp = <double*> realloc(
                                expanded_coords_p_temp, natoms * 3 * sizeof(double)
                            )
                            indices_p_temp = <long*> realloc(
                                indices_p_temp, natoms * sizeof(long)
                            )
                        if (
                                offset_final is NULL or
                                expanded_coords_p_temp is NULL or
                                indices_p_temp is NULL
                        ):
                            failed_malloc = 1
                            break
                else:
                    continue
                break
            else:
                continue
            break
        else:
            continue
        break

    if failed_malloc:
        raise MemoryError("A realloc of memory of failed!")
    else:
        failed_malloc = 0

    # if no valid neighbors were found return empty
    if count == 0:
        free(&frac_coords[0, 0])
        free(&all_fcoords[0, 0])
        free(&coords_in_cell[0, 0])
        free(&offset_correction[0, 0])
        free(&center_indices1[0])
        free(&center_indices3[0, 0])

        free(offsets_p_temp)
        free(expanded_coords_p_temp)
        free(indices_p_temp)

        free(index_1)
        free(index_2)
        free(offset_final)
        free(distances)

        return (np.array([], dtype=int), np.array([], dtype=int),
            np.array([[], [], []], dtype=float).T, np.array([], dtype=float))

    natoms = count
    cdef:
        # Delete those beyond (min_center_coords - r, max_center_coords + r)
        double *offsets_p = <double*> safe_realloc(
            offsets_p_temp, count * 3 * sizeof(double)
        )
        double *expanded_coords_p = <double*> safe_realloc(
            expanded_coords_p_temp, count * 3 * sizeof(double)
        )
        long *indices_p = <long*> safe_realloc(
            indices_p_temp, count * sizeof(long)
        )

        double[:, ::1] offsets = <double[:count, :3]> offsets_p
        double[:, ::1] expanded_coords = <double[:count, :3]> expanded_coords_p
        long[::1] indices = <long[:count]> indices_p

        # Construct linked cell list
        long[:, ::1] all_indices3 = <long[:natoms, :3]> safe_malloc(
            natoms * 3 * sizeof(long)
        )
        long[::1] all_indices1 = <long[:natoms]> safe_malloc(
            natoms * sizeof(long)
        )

    for i in range(3):
        ncube[i] = <long>(ceil((valid_max[i] - valid_min[i]) / ledge))

    compute_cube_index(expanded_coords, valid_min, ledge, all_indices3)
    three_to_one(all_indices3, ncube[1], ncube[2], all_indices1)

    cdef:
        long nb_cubes = ncube[0] * ncube[1] * ncube[2]
        long *head = <long*> safe_malloc(nb_cubes*sizeof(long))
        long *atom_indices = <long*> safe_malloc(natoms*sizeof(long))
        long[:, ::1] neighbor_map = <long[:nb_cubes, :27]> safe_malloc(
            nb_cubes * 27 * sizeof(long)
        )

    memset(<void*>head, -1, nb_cubes*sizeof(long))
    memset(<void*>atom_indices, -1, natoms*sizeof(long))

    get_cube_neighbors(ncube, neighbor_map)
    for i in range(natoms):
        atom_indices[i] = head[all_indices1[i]]
        head[all_indices1[i]] = i

    # Get center atoms' cube indices
    compute_cube_index(center_coords, valid_min, ledge, center_indices3)
    three_to_one(center_indices3, ncube[1], ncube[2], center_indices1)

    count = 0
    for i in range(n_center):
        for j in range(27):
            if neighbor_map[center_indices1[i], j] == -1:
                continue
            cube_index_temp = neighbor_map[center_indices1[i], j]
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
                    # increasing the memory size by allocating incrementally
                    # malloc_chunk more memory locations, I found it 3x faster to do so
                    # compared to using vectors in cpp
                    if count >= malloc_chunk:
                        malloc_chunk += malloc_chunk  # double the size
                        index_1 = <long*> realloc(index_1, malloc_chunk * sizeof(long))
                        index_2 = <long*> realloc(index_2, malloc_chunk*sizeof(long))
                        offset_final = <double*> realloc(
                            offset_final, 3*malloc_chunk*sizeof(double)
                        )
                        distances = <double*> realloc(
                            distances, malloc_chunk*sizeof(double)
                        )
                        if (
                                index_1 is NULL or index_2 is NULL or
                                offset_final is NULL or distances is NULL
                        ):
                            failed_malloc = 1
                            break
                link_index = atom_indices[link_index]
            else:
                continue
            break
        else:
            continue
        break

    if failed_malloc:
        raise MemoryError("A realloc of memory of failed!")
    else:
        failed_malloc = 0

    if count == 0:
        py_index_1 = np.array([], dtype=int)
        py_index_2 = np.array([], dtype=int)
        py_offsets = np.array([[], [], []], dtype=float).T
        py_distances = np.array([], dtype=float)
    else:
        # resize to the actual size
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
    free(offsets_p)
    free(expanded_coords_p)
    free(indices_p)
    free(head)
    free(atom_indices)

    free(&offset_correction[0, 0])
    free(&frac_coords[0, 0])
    free(&all_fcoords[0, 0])
    free(&coords_in_cell[0, 0])
    free(&center_indices1[0])
    free(&center_indices3[0, 0])
    free(&neighbor_map[0, 0])

    free(&all_indices3[0, 0])
    free(&all_indices1[0])

    return py_index_1, py_index_2, py_offsets, py_distances


cdef void get_cube_neighbors(long[3] ncube, long[:, ::1] neighbor_map):
    """
    Get {cube_index: cube_neighbor_indices} map
    """
    cdef:
        int i, j, k
        int count = 0
        long ncubes = ncube[0] * ncube[1] * ncube[2]
        long[::1] counts = <long[:ncubes]> safe_malloc(ncubes * sizeof(long))
        long[:, ::1] cube_indices_3d = <long[:ncubes, :3]> safe_malloc(
            ncubes*3*sizeof(long)
        )
        long[::1] cube_indices_1d = <long[:ncubes]> safe_malloc(ncubes*sizeof(long))

        # creating the memviews of c-arrays once substantially improves speed
        # but for some reason it makes the runtime scaling with the number of
        # atoms worse
        long[1][3] index3_arr
        long[:, ::1] index3 = index3_arr
        long[1] index1_arr
        long[::1] index1 = index1_arr

        int n = 1
        long ntotal = (2 * n + 1) * (2 * n + 1) * (2 * n + 1)
        long[:, ::1] ovectors
        long *ovectors_p = <long *> safe_malloc(ntotal * 3 * sizeof(long))
        int n_ovectors = compute_offset_vectors(ovectors_p, n)

    # now resize to the actual size
    ovectors_p = <long *> safe_realloc(ovectors_p, n_ovectors * 3 * sizeof(long))
    ovectors = <long[:n_ovectors, :3]>ovectors_p

    memset(<void*>&neighbor_map[0, 0], -1, neighbor_map.shape[0] * 27 * sizeof(long))

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

    for i in range(ncubes):
        for j in range(27):
            index3[0, 0] = ovectors[j, 0] + cube_indices_3d[i, 0]
            index3[0, 1] = ovectors[j, 1] + cube_indices_3d[i, 1]
            index3[0, 2] = ovectors[j, 2] + cube_indices_3d[i, 2]
            if (
                    (index3[0, 0] < ncube[0]) &
                    (index3[0, 0] >= 0) &
                    (index3[0, 1] < ncube[1]) &
                    (index3[0, 1] >= 0) &
                    (index3[0, 2] < ncube[2]) &
                    (index3[0, 2] >= 0)
            ):
                three_to_one(index3, ncube[1], ncube[2], index1)
                neighbor_map[i, counts[i]] = index1[0]
                counts[i] += 1

    free(&cube_indices_3d[0, 0])
    free(&cube_indices_1d[0])
    free(&counts[0])
    free(ovectors_p)


cdef int compute_offset_vectors(long* ovectors, long n) nogil:
    cdef:
        int i, j, k, ind
        int count = 0
        double center[8][3] # center vertices coords
        double offset[8][3]  # offsetted vertices

    for i in range(2):
        for j in range(2):
            for k in range(2):
                ind = i * 4 + j * 2 + k
                center[ind][0] = i - 0.5
                center[ind][1] = j - 0.5
                center[ind][2] = k - 0.5

    for i in range(-n, n + 1):
        for j in range(-n, n + 1):
            for k in range(-n, n + 1):
                offset_cube(center, i, j, k, offset)
                if distance_vertices(center, offset, n):
                    ovectors[3*count] = i
                    ovectors[3*count + 1] = j
                    ovectors[3*count + 2] = k
                    count += 1

    return count


cdef double distance2(
        const double[:, ::1] m1,
        const double[:, ::1] m2,
        long index1,
        long index2,
        long size
    ) nogil:
    """
    Faster way to compute the distance squared by not using slice but providing indices
    in each matrix
    """
    cdef:
        int i
        double s = 0

    for i in range(size):
        s += (m1[index1, i] - m2[index2, i]) * (m1[index1, i] - m2[index2, i])
    return s


cdef void get_bounds(
        const double[:, ::1] frac_coords,
        const double[3] maxr,
        const long[3] pbc,
        long[3] max_bounds,
        long[3] min_bounds
    ) nogil:
    """
    Given the fractional coordinates and the number of repeation needed in each
    direction, maxr, compute the translational bounds in each dimension
    """
    cdef:
        int i
        double[3] max_fcoords
        double[3] min_fcoords

    max_and_min(frac_coords, max_fcoords, min_fcoords)

    for i in range(3):
        min_bounds[i] = 0
        max_bounds[i] = 1

    for i in range(3):
        if pbc[i]:
            min_bounds[i] = <long>(floor(min_fcoords[i] - maxr[i] - 1e-8))
            max_bounds[i] = <long>(ceil(max_fcoords[i] + maxr[i] + 1e-8))

cdef void get_frac_coords(
        const double[:, ::1] lattice,
        double[:, ::1] inv_lattice,
        const double[:, ::1] cart_coords,
        double[:, ::1] frac_coords
    ) nogil:
    """
    Compute the fractional coordinates
    """
    matrix_inv(lattice, inv_lattice)
    matmul(cart_coords, inv_lattice, frac_coords)

cdef void matmul(
        const double[:, ::1] m1,
        const double[:, ::1] m2,
        double [:, ::1] out
    ) nogil:
    """
    Matrix multiplication
    """
    cdef:
        int i, j, k
        int m = m1.shape[0], n = m1.shape[1], l = m2.shape[1]

    for i in range(m):
        for j in range(l):
            out[i, j] = 0
            for k in range(n):
                out[i, j] += m1[i, k] * m2[k, j]

cdef void matrix_inv(const double[:, ::1] matrix, double[:, ::1] inv) nogil:
    """
    Matrix inversion
    """
    cdef:
        int i, j
        double det = matrix_det(matrix)

    for i in range(3):
        for j in range(3):
            inv[i, j] = (matrix[(j+1)%3, (i+1)%3] * matrix[(j+2)%3, (i+2)%3] -
                matrix[(j+2)%3, (i+1)%3] * matrix[(j+1)%3, (i+2)%3]) / det

cdef double matrix_det(const double[:, ::1] matrix) nogil:
    """
    Matrix determinant
    """
    return (
        matrix[0, 0] * (matrix[1, 1] * matrix[2, 2] - matrix[1, 2] * matrix[2, 1]) +
        matrix[0, 1] * (matrix[1, 2] * matrix[2, 0] - matrix[1, 0] * matrix[2, 2]) +
        matrix[0, 2] * (matrix[1, 0] * matrix[2, 1] - matrix[1, 1] * matrix[2, 0])
    )

cdef void get_max_r(
        const double[:, ::1] reciprocal_lattice,
        double[3] maxr,
        double r
    ) nogil:
    """
    Get maximum repetition in each directions
    """
    cdef:
        int i
        double recp_len

    for i in range(3):  # is it ever not 3x3 for our cases?
        recp_len = norm(reciprocal_lattice[i, :])
        maxr[i] = ceil((r + 0.15) * recp_len / (2 * pi))

cdef void get_reciprocal_lattice(
        const double[:, ::1] lattice,
        double[:, ::1] reciprocal_lattice
    ) nogil:
    """
    Compute the reciprocal lattice
    """
    cdef int i
    for i in range(3):
        recip_component(
            lattice[i, :], lattice[(i+1)%3, :],
            lattice[(i+2)%3, :],
            reciprocal_lattice[i, :]
        )

cdef void recip_component(
        const double[::1] a1,
        const double[::1] a2,
        const double[::1] a3,
        double[::1] out
    ) nogil:
    """
    Compute the reciprocal lattice vector
    """
    cdef:
        int i
        double prod
        double ai_cross_aj[3]

    cross(&a2[0], &a3[0], ai_cross_aj)  # need to pass address of memviews
    prod = inner(&a1[0], ai_cross_aj)
    for i in range(3):
        out[i] = 2 * pi * ai_cross_aj[i] / prod

cdef double inner(const double[3] x, const double[3] y) nogil:
    """
    Compute inner product of 3d vectors
    """
    cdef:
        double sum = 0
        int i

    for i in range(3):
        sum += x[i] * y[i]
    return sum

cdef void cross(const double[3] x, const double[3] y, double[3] out) nogil:
    """
    Cross product of vector x and y, output in out
    """
    out[0] = x[1] * y[2] - x[2] * y[1]
    out[1] = x[2] * y[0] - x[0] * y[2]
    out[2] = x[0] * y[1] - x[1] * y[0]

cdef double norm(const double[::1] vec) nogil:
    """
    Vector norm
    """
    cdef:
        int i
        int n = vec.shape[0]
        double sum = 0

    for i in range(n):
        sum += vec[i] * vec[i]
    return sqrt(sum)

cdef void max_and_min(
        const double[:, ::1] coords,
        double[3] max_coords,
        double[3] min_coords
    ) nogil:
    """
    Compute the min and max of coords
    """
    cdef:
        int i, j
        int M = coords.shape[0]
        int N = coords.shape[1]

    for i in range(N):
        max_coords[i] = coords[0, i]
        min_coords[i] = coords[0, i]
    for i in range(M):
        for j in range(N):
            if coords[i, j] >= max_coords[j]:
                max_coords[j] = coords[i, j]
            if coords[i, j] <= min_coords[j]:
                min_coords[j] = coords[i, j]

cdef void compute_cube_index(
        const double[:, ::1] coords,
        const double[3] global_min,
        double radius, long[:, ::1] return_indices
    ) nogil:
    cdef int i, j
    for i in range(coords.shape[0]):
        for j in range(coords.shape[1]):
            return_indices[i, j] = <long>(
                floor((coords[i, j] - global_min[j] + 1e-8) / radius)
            )


cdef void three_to_one(
        const long[:, ::1] label3d, long ny, long nz, long[::1] label1d
    ) nogil:
    """
    3D vector representation to 1D
    """
    cdef:
        int i
        int n = label3d.shape[0]

    for i in range(n):
        label1d[i] = label3d[i, 0] * ny * nz + label3d[i, 1] * nz + label3d[i, 2]


cdef bint distance_vertices(
        const double[8][3] center, const double[8][3] off, double r
    ) nogil:
    cdef:
        int i, j
        double d2
        double r2 = r * r

    for i in range(8):
        for j in range(8):
            d2 = (center[i][0] - off[j][0]) * (center[i][0] - off[j][0]) + \
                 (center[i][1] - off[j][1]) * (center[i][1] - off[j][1]) + \
                 (center[i][2] - off[j][2]) * (center[i][2] - off[j][2])
            if d2 <= r2:
                return 1
    return 0

cdef void offset_cube(
        const double[8][3] center,
        long n, long m, long l,
        const double[8][3] (&offsetted)
    ) nogil:
    cdef int i, j, k
    for i in range(2):
        for j in range(2):
            for k in range(2):
                ind = i * 4 + j * 2 + k
                offsetted[ind][0] = center[ind][0] + n
                offsetted[ind][1] = center[ind][1] + m
                offsetted[ind][2] = center[ind][2] + l
