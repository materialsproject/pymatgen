---
layout: default
title: pymatgen.core.lattice.md
nav_exclude: true
---

# pymatgen.core.lattice module

Defines the classes relating to 3D lattices.


### _class_ pymatgen.core.lattice.Lattice(matrix: ArrayLike, pbc: tuple[bool, bool, bool] = (True, True, True))
Bases: `MSONable`

A lattice object. Essentially a matrix with conversion matrices. In
general, it is assumed that length units are in Angstroms and angles are in
degrees unless otherwise stated.

Create a lattice from any sequence of 9 numbers. Note that the sequence
is assumed to be read one row at a time. Each row represents one
lattice vector.


* **Parameters**


    * **matrix** – Sequence of numbers in any form. Examples of acceptable
    input.
    i) An actual numpy array.
    ii) [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    iii) [1, 0, 0 , 0, 1, 0, 0, 0, 1]
    iv) (1, 0, 0, 0, 1, 0, 0, 0, 1)
    Each row should correspond to a lattice vector.
    E.g., [[10, 0, 0], [20, 10, 0], [0, 0, 30]] specifies a lattice
    with lattice vectors [10, 0, 0], [20, 10, 0] and [0, 0, 30].


    * **pbc** – a tuple defining the periodic boundary conditions along the three
    axis of the lattice. If None periodic in all directions.



#### _property_ a(_: floa_ )
*a* lattice parameter.


#### _property_ abc(_: Vector3_ )
Lengths of the lattice vectors, i.e. (a, b, c).


#### _property_ alpha(_: floa_ )
Angle alpha of lattice in degrees.


#### _property_ angles(_: Vector3_ )
Lattice angles.


* **Returns**

    The angles (alpha, beta, gamma) of the lattice.



#### as_dict(verbosity: int = 0)
Json-serialization dict representation of the Lattice.


* **Parameters**

    **verbosity** (*int*) – Verbosity level. Default of 0 only includes the
    matrix representation. Set to 1 for more details.



#### _property_ b(_: floa_ )
*b* lattice parameter.


#### _property_ beta(_: floa_ )
Angle beta of lattice in degrees.


#### _property_ c(_: floa_ )
*c* lattice parameter.


#### copy()
Deep copy of self.


#### _static_ cubic(a: float, pbc: tuple[bool, bool, bool] = (True, True, True))
Convenience constructor for a cubic lattice.


* **Parameters**


    * **a** (*float*) – The *a* lattice parameter of the cubic cell.


    * **pbc** (*tuple*) – a tuple defining the periodic boundary conditions along the three
    axis of the lattice. If None periodic in all directions.



* **Returns**

    Cubic lattice of dimensions a x a x a.



#### d_hkl(miller_index: ArrayLike)
Returns the distance between the hkl plane and the origin.


* **Parameters**

    **miller_index** (*[**h**,**k**,**l**]*) – Miller index of plane



* **Returns**

    d_hkl (float)



#### dot(coords_a: ArrayLike, coords_b: ArrayLike, frac_coords: bool = False)
Compute the scalar product of vector(s).


* **Parameters**


    * **coords_a** – Array-like coordinates.


    * **coords_b** – Array-like coordinates.


    * **frac_coords** (*bool*) – Boolean stating whether the vector
    corresponds to fractional or Cartesian coordinates.



* **Returns**

    one-dimensional numpy array.



#### find_all_mappings(other_lattice: Lattice, ltol: float = 1e-05, atol: float = 1, skip_rotation_matrix: bool = False)
Finds all mappings between current lattice and another lattice.


* **Parameters**


    * **other_lattice** (*Lattice*) – Another lattice that is equivalent to
    this one.


    * **ltol** (*float*) – Tolerance for matching lengths. Defaults to 1e-5.


    * **atol** (*float*) – Tolerance for matching angles. Defaults to 1.


    * **skip_rotation_matrix** (*bool*) – Whether to skip calculation of the
    rotation matrix



* **Yields**

    (aligned_lattice, rotation_matrix, scale_matrix) if a mapping is
    found. aligned_lattice is a rotated version of other_lattice that
    has the same lattice parameters, but which is aligned in the
    coordinate system of this lattice so that translational points
    match up in 3D. rotation_matrix is the rotation that has to be
    applied to other_lattice to obtain aligned_lattice, i.e.,
    aligned_matrix = np.inner(other_lattice, rotation_matrix) and
    op = SymmOp.from_rotation_and_translation(rotation_matrix)
    aligned_matrix = op.operate_multi(latt.matrix)
    Finally, scale_matrix is the integer matrix that expresses
    aligned_matrix as a linear combination of this
    lattice, i.e., aligned_matrix = np.dot(scale_matrix, self.matrix)

    None is returned if no matches are found.



#### find_mapping(other_lattice: Lattice, ltol: float = 1e-05, atol: float = 1, skip_rotation_matrix: bool = False)
Finds a mapping between current lattice and another lattice. There
are an infinite number of choices of basis vectors for two entirely
equivalent lattices. This method returns a mapping that maps
other_lattice to this lattice.


* **Parameters**


    * **other_lattice** (*Lattice*) – Another lattice that is equivalent to
    this one.


    * **ltol** (*float*) – Tolerance for matching lengths. Defaults to 1e-5.


    * **atol** (*float*) – Tolerance for matching angles. Defaults to 1.


    * **skip_rotation_matrix** (*bool*) – Whether to skip calculation of the rotation matrix.
    Defaults to False.



* **Returns**

    (aligned_lattice, rotation_matrix, scale_matrix) if a mapping is
    found. aligned_lattice is a rotated version of other_lattice that
    has the same lattice parameters, but which is aligned in the
    coordinate system of this lattice so that translational points
    match up in 3D. rotation_matrix is the rotation that has to be
    applied to other_lattice to obtain aligned_lattice, i.e.,
    aligned_matrix = np.inner(other_lattice, rotation_matrix) and
    op = SymmOp.from_rotation_and_translation(rotation_matrix)
    aligned_matrix = op.operate_multi(latt.matrix)
    Finally, scale_matrix is the integer matrix that expresses
    aligned_matrix as a linear combination of this
    lattice, i.e., aligned_matrix = np.dot(scale_matrix, self.matrix)

    None is returned if no matches are found.




#### _classmethod_ from_dict(d: dict, fmt: str | None = None, \*\*kwargs)
Create a Lattice from a dictionary containing the a, b, c, alpha, beta,
and gamma parameters if fmt is None.

If fmt == “abivars”, the function build a Lattice object from a
dictionary with the Abinit variables acell and rprim in Bohr.
If acell is not given, the Abinit default is used i.e. [1,1,1] Bohr

### Example

Lattice.from_dict(fmt=”abivars”, acell=3\*[10], rprim=np.eye(3))


#### _classmethod_ from_parameters(a: float, b: float, c: float, alpha: float, beta: float, gamma: float, vesta: bool = False, pbc: tuple[bool, bool, bool] = (True, True, True))
Create a Lattice using unit cell lengths (in Angstrom) and angles (in degrees).


* **Parameters**


    * **a** (*float*) – *a* lattice parameter.


    * **b** (*float*) – *b* lattice parameter.


    * **c** (*float*) – *c* lattice parameter.


    * **alpha** (*float*) – *alpha* angle in degrees.


    * **beta** (*float*) – *beta* angle in degrees.


    * **gamma** (*float*) – *gamma* angle in degrees.


    * **vesta** – True if you import Cartesian coordinates from VESTA.


    * **pbc** (*tuple*) – a tuple defining the periodic boundary conditions along the three
    axis of the lattice. If None periodic in all directions.



* **Returns**

    Lattice with the specified lattice parameters.



#### _property_ gamma(_: floa_ )
Angle gamma of lattice in degrees.


#### get_all_distances(fcoords1: ArrayLike, fcoords2: ArrayLike)
Returns the distances between two lists of coordinates taking into
account periodic boundary conditions and the lattice. Note that this
computes an MxN array of distances (i.e. the distance between each
point in fcoords1 and every coordinate in fcoords2). This is
different functionality from pbc_diff.


* **Parameters**


    * **fcoords1** – First set of fractional coordinates. e.g., [0.5, 0.6,
    0.7] or [[1.1, 1.2, 4.3], [0.5, 0.6, 0.7]]. It can be a single
    coord or any array of coords.


    * **fcoords2** – Second set of fractional coordinates.



* **Returns**

    2d array of Cartesian distances. E.g the distance between
    fcoords1[i] and fcoords2[j] is distances[i,j]



#### get_brillouin_zone()
Returns the Wigner-Seitz cell for the reciprocal lattice, aka the
Brillouin Zone.


* **Returns**

    A list of list of coordinates.
    Each element in the list is a “facet” of the boundary of the
    Brillouin Zone. For instance, a list of four coordinates will
    represent a square facet.



#### get_cartesian_coords(fractional_coords: ArrayLike)
Returns the Cartesian coordinates given fractional coordinates.


* **Parameters**

    **fractional_coords** (*3x1 array*) – Fractional coords.



* **Returns**

    Cartesian coordinates



#### get_distance_and_image(frac_coords1: ArrayLike, frac_coords2: ArrayLike, jimage: ArrayLike | None = None)
Gets distance between two frac_coords assuming periodic boundary
conditions. If the index jimage is not specified it selects the j
image nearest to the i atom and returns the distance and jimage
indices in terms of lattice vector translations. If the index jimage
is specified it returns the distance between the frac_coords1 and
the specified jimage of frac_coords2, and the given jimage is also
returned.


* **Parameters**


    * **frac_coords1** (*3x1 array*) – Reference fcoords to get distance from.


    * **frac_coords2** (*3x1 array*) – fcoords to get distance from.


    * **jimage** (*3x1 array*) – Specific periodic image in terms of
    lattice translations, e.g., [1,0,0] implies to take periodic
    image that is one a-lattice vector away. If jimage is None,
    the image that is nearest to the site is found.



* **Returns**

    distance and periodic lattice translations
    of the other site for which the distance applies. This means that
    the distance between frac_coords1 and (jimage + frac_coords2) is
    equal to distance.



* **Return type**

    (distance, jimage)



#### get_frac_coords_from_lll(lll_frac_coords: ArrayLike)
Given fractional coordinates in the lll basis, returns corresponding
fractional coordinates in the lattice basis.


#### get_fractional_coords(cart_coords: ArrayLike)
Returns the fractional coordinates given Cartesian coordinates.


* **Parameters**

    **cart_coords** (*3x1 array*) – Cartesian coords.



* **Returns**

    Fractional coordinates.



#### get_lll_frac_coords(frac_coords: ArrayLike)
Given fractional coordinates in the lattice basis, returns corresponding
fractional coordinates in the lll basis.


#### get_lll_reduced_lattice(delta: float = 0.75)

* **Parameters**

    **delta** – Delta parameter.



* **Returns**

    LLL reduced Lattice.



#### get_miller_index_from_coords(coords: ArrayLike, coords_are_cartesian: bool = True, round_dp: int = 4, verbose: bool = True)
Get the Miller index of a plane from a list of site coordinates.

A minimum of 3 sets of coordinates are required. If more than 3 sets of
coordinates are given, the best plane that minimises the distance to all
points will be calculated.


* **Parameters**


    * **coords** (*iterable*) – A list or numpy array of coordinates. Can be
    Cartesian or fractional coordinates. If more than three sets of
    coordinates are provided, the best plane that minimises the
    distance to all sites will be calculated.


    * **coords_are_cartesian** (*bool**, **optional*) – Whether the coordinates are
    in Cartesian space. If using fractional coordinates set to
    False.


    * **round_dp** (*int**, **optional*) – The number of decimal places to round the
    miller index to.


    * **verbose** (*bool**, **optional*) – Whether to print warnings.



* **Returns**

    The Miller index.



* **Return type**

    (tuple)



#### get_niggli_reduced_lattice(tol: float = 1e-05)
Get the Niggli reduced lattice using the numerically stable algo
proposed by R. W. Grosse-Kunstleve, N. K. Sauter, & P. D. Adams,
Acta Crystallographica Section A Foundations of Crystallography, 2003,
60(1), 1-6. doi:10.1107/S010876730302186X.


* **Parameters**

    **tol** (*float*) – The numerical tolerance. The default of 1e-5 should
    result in stable behavior for most cases.



* **Returns**

    Niggli-reduced lattice.



#### get_points_in_sphere(frac_points: ArrayLike, center: ArrayLike, r: float, zip_results=True)
Find all points within a sphere from the point taking into account
periodic boundary conditions. This includes sites in other periodic
images.

Algorithm:


1. place sphere of radius r in crystal and determine minimum supercell
(parallelpiped) which would contain a sphere of radius r. for this
we need the projection of a_1 on a unit vector perpendicular
to a_2 & a_3 (i.e. the unit vector in the direction b_1) to
determine how many a_1”s it will take to contain the sphere.

Nxmax = r \* length_of_b_1 / (2 Pi)


2. keep points falling within r.


* **Parameters**


    * **frac_points** – All points in the lattice in fractional coordinates.


    * **center** – Cartesian coordinates of center of sphere.


    * **r** – radius of sphere.


    * **zip_results** (*bool*) – Whether to zip the results together to group by
    point, or return the raw fcoord, dist, index arrays



* **Returns**

    [(fcoord, dist, index, supercell_image) …] since most of the time, subsequent

        processing requires the distance, index number of the atom, or index of the image

    else:

        fcoords, dists, inds, image




* **Return type**

    if zip_results



#### get_points_in_sphere_old(\*\*kwargs)

#### get_points_in_sphere_py(frac_points: ArrayLike, center: ArrayLike, r: float, zip_results=True)
Find all points within a sphere from the point taking into account
periodic boundary conditions. This includes sites in other periodic
images.

Algorithm:


1. place sphere of radius r in crystal and determine minimum supercell
(parallelpiped) which would contain a sphere of radius r. for this
we need the projection of a_1 on a unit vector perpendicular
to a_2 & a_3 (i.e. the unit vector in the direction b_1) to
determine how many a_1”s it will take to contain the sphere.

Nxmax = r \* length_of_b_1 / (2 Pi)


2. keep points falling within r.


* **Parameters**


    * **frac_points** – All points in the lattice in fractional coordinates.


    * **center** – Cartesian coordinates of center of sphere.


    * **r** – radius of sphere.


    * **zip_results** (*bool*) – Whether to zip the results together to group by
    point, or return the raw fcoord, dist, index arrays



* **Returns**

    [(fcoord, dist, index, supercell_image) …] since most of the time, subsequent

        processing requires the distance, index number of the atom, or index of the image

    else:

        fcoords, dists, inds, image




* **Return type**

    if zip_results



#### get_recp_symmetry_operation(symprec: float = 0.01)
Find the symmetric operations of the reciprocal lattice,
to be used for hkl transformations
:param symprec: default is 0.001.


#### get_vector_along_lattice_directions(cart_coords: ArrayLike)
Returns the coordinates along lattice directions given Cartesian coordinates.

Note, this is different than a projection of the Cartesian vector along the
lattice parameters. It is simply the fractional coordinates multiplied by the
lattice vector magnitudes.

For example, this method is helpful when analyzing the dipole moment (in
units of electron Angstroms) of a ferroelectric crystal. See the Polarization
class in pymatgen.analysis.ferroelectricity.polarization.


* **Parameters**

    **cart_coords** (*3x1 array*) – Cartesian coords.



* **Returns**

    Lattice coordinates.



#### get_wigner_seitz_cell()
Returns the Wigner-Seitz cell for the given lattice.


* **Returns**

    A list of list of coordinates.
    Each element in the list is a “facet” of the boundary of the
    Wigner Seitz cell. For instance, a list of four coordinates will
    represent a square facet.



#### _static_ hexagonal(a: float, c: float, pbc: tuple[bool, bool, bool] = (True, True, True))
Convenience constructor for a hexagonal lattice.


* **Parameters**


    * **a** (*float*) – *a* lattice parameter of the hexagonal cell.


    * **c** (*float*) – *c* lattice parameter of the hexagonal cell.


    * **pbc** (*tuple*) – a tuple defining the periodic boundary conditions along the three
    axis of the lattice. If None periodic in all directions.



* **Returns**

    Hexagonal lattice of dimensions a x a x c.



#### _property_ inv_matrix(_: ndarra_ )
Inverse of lattice matrix.


#### _property_ is_3d_periodic(_: boo_ )
True if the Lattice is periodic in all directions.


#### is_hexagonal(hex_angle_tol: float = 5, hex_length_tol: float = 0.01)

* **Parameters**


    * **hex_angle_tol** – Angle tolerance


    * **hex_length_tol** – Length tolerance



* **Returns**

    Whether lattice corresponds to hexagonal lattice.



#### _property_ is_orthogonal(_: boo_ )
Whether all angles are 90 degrees.


* **Type**

    return



#### _property_ lengths(_: Vector3_ )
Lattice lengths.


* **Returns**

    The lengths (a, b, c) of the lattice.



#### _property_ lll_inverse(_: ndarra_ )
Inverse of self.lll_mapping.


* **Type**

    return



#### _property_ lll_mapping(_: ndarra_ )
The mapping between the LLL reduced lattice and the original
lattice.


* **Type**

    return



#### _property_ lll_matrix(_: ndarra_ )
The matrix for LLL reduction


* **Type**

    return



#### _property_ matrix(_: ndarra_ )
Copy of matrix representing the Lattice.


#### _property_ metric_tensor(_: ndarra_ )
The metric tensor of the lattice.


#### _static_ monoclinic(a: float, b: float, c: float, beta: float, pbc: tuple[bool, bool, bool] = (True, True, True))
Convenience constructor for a monoclinic lattice.


* **Parameters**


    * **a** (*float*) – *a* lattice parameter of the monoclinc cell.


    * **b** (*float*) – *b* lattice parameter of the monoclinc cell.


    * **c** (*float*) – *c* lattice parameter of the monoclinc cell.


    * **beta** (*float*) – *beta* angle between lattice vectors b and c in
    degrees.


    * **pbc** (*tuple*) – a tuple defining the periodic boundary conditions along the three
    axis of the lattice. If None periodic in all directions.



* **Returns**

    Monoclinic lattice of dimensions a x b x c with non right-angle
    beta between lattice vectors a and c.



#### norm(coords: ArrayLike, frac_coords: bool = True)
Compute the norm of vector(s).


* **Parameters**


    * **coords** – Array-like object with the coordinates.


    * **frac_coords** – Boolean stating whether the vector corresponds to fractional or
    Cartesian coordinates.



* **Returns**

    one-dimensional numpy array.



#### _static_ orthorhombic(a: float, b: float, c: float, pbc: tuple[bool, bool, bool] = (True, True, True))
Convenience constructor for an orthorhombic lattice.


* **Parameters**


    * **a** (*float*) – *a* lattice parameter of the orthorhombic cell.


    * **b** (*float*) – *b* lattice parameter of the orthorhombic cell.


    * **c** (*float*) – *c* lattice parameter of the orthorhombic cell.


    * **pbc** (*tuple*) – a tuple defining the periodic boundary conditions along the three
    axis of the lattice. If None periodic in all directions.



* **Returns**

    Orthorhombic lattice of dimensions a x b x c.



#### _property_ parameters(_: tuple[float, float, float, float, float, float_ )
(a, b, c, alpha, beta, gamma).


* **Type**

    Returns



#### _property_ pbc(_: tuple[bool, bool, bool_ )
Tuple defining the periodicity of the Lattice.


#### _property_ reciprocal_lattice(_: Lattic_ )
Return the reciprocal lattice. Note that this is the standard
reciprocal lattice used for solid state physics with a factor of 2 \*
pi. If you are looking for the crystallographic reciprocal lattice,
use the reciprocal_lattice_crystallographic property.
The property is lazily generated for efficiency.


#### _property_ reciprocal_lattice_crystallographic(_: Lattic_ )
Returns the *crystallographic* reciprocal lattice, i.e., no factor of
2 \* pi.


#### _static_ rhombohedral(a: float, alpha: float, pbc: tuple[bool, bool, bool] = (True, True, True))
Convenience constructor for a rhombohedral lattice.


* **Parameters**


    * **a** (*float*) – *a* lattice parameter of the rhombohedral cell.


    * **alpha** (*float*) – Angle for the rhombohedral lattice in degrees.


    * **pbc** (*tuple*) – a tuple defining the periodic boundary conditions along the three
    axis of the lattice. If None periodic in all directions.



* **Returns**

    Rhombohedral lattice of dimensions a x a x a.



#### scale(new_volume: float)
Return a new Lattice with volume new_volume by performing a
scaling of the lattice vectors so that length proportions and angles
are preserved.


* **Parameters**

    **new_volume** – New volume to scale to.



* **Returns**

    New lattice with desired volume.



#### selling_dist(other)
Returns the minimum Selling distance between two lattices.


#### _property_ selling_vector(_: ndarra_ )
Returns the (1,6) array of Selling Scalars.


#### _static_ tetragonal(a: float, c: float, pbc: tuple[bool, bool, bool] = (True, True, True))
Convenience constructor for a tetragonal lattice.


* **Parameters**


    * **a** (*float*) – *a* lattice parameter of the tetragonal cell.


    * **c** (*float*) – *c* lattice parameter of the tetragonal cell.


    * **pbc** (*tuple*) – a tuple defining the periodic boundary conditions along the three
    axis of the lattice. If None periodic in all directions.



* **Returns**

    Tetragonal lattice of dimensions a x a x c.



#### _property_ volume(_: floa_ )
Volume of the unit cell in Angstrom^3.


### pymatgen.core.lattice.find_neighbors(label: ndarray, nx: int, ny: int, nz: int)
Given a cube index, find the neighbor cube indices.


* **Parameters**


    * **label** – (array) (n,) or (n x 3) indice array


    * **nx** – (int) number of cells in y direction


    * **ny** – (int) number of cells in y direction


    * **nz** – (int) number of cells in z direction



* **Returns**

    Neighbor cell indices.



### pymatgen.core.lattice.get_integer_index(miller_index: Sequence[float], round_dp: int = 4, verbose: bool = True)
Attempt to convert a vector of floats to whole numbers.


* **Parameters**


    * **miller_index** (*list** of **float*) – A list miller indexes.


    * **round_dp** (*int**, **optional*) – The number of decimal places to round the
    miller index to.


    * **verbose** (*bool**, **optional*) – Whether to print warnings.



* **Returns**

    The Miller index.



* **Return type**

    (tuple)



### pymatgen.core.lattice.get_points_in_spheres(all_coords: np.ndarray, center_coords: np.ndarray, r: float, pbc: bool | list[bool] | tuple[bool, bool, bool] = True, numerical_tol: float = 1e-08, lattice: Lattice | None = None, return_fcoords: bool = False)
For each point in center_coords, get all the neighboring points in all_coords that are within the
cutoff radius r.


* **Parameters**


    * **all_coords** – (list of Cartesian coordinates) all available points


    * **center_coords** – (list of Cartesian coordinates) all centering points


    * **r** – (float) cutoff radius


    * **pbc** – (bool or a list of bool) whether to set periodic boundaries


    * **numerical_tol** – (float) numerical tolerance


    * **lattice** – (Lattice) lattice to consider when PBC is enabled


    * **return_fcoords** – (bool) whether to return fractional coords when pbc is set.



* **Returns**

    List[List[Tuple[coords, distance, index, image]]]