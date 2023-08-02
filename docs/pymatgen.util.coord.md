---
layout: default
title: pymatgen.util.coord.md
nav_exclude: true
---

# pymatgen.util.coord module

Utilities for manipulating coordinates or list of coordinates, under periodic
boundary conditions or otherwise. Many of these are heavily vectorized in
numpy for performance.


### _class_ pymatgen.util.coord.Simplex(coords)
Bases: `MSONable`

A generalized simplex object. See [http://en.wikipedia.org/wiki/Simplex](http://en.wikipedia.org/wiki/Simplex).

<!-- attribute: space_dim

Dimension of the space. Usually, this is 1 more than the simplex_dim. -->
<!-- attribute: simplex_dim

Dimension of the simplex coordinate space. -->
Initializes a Simplex from vertex coordinates.


* **Parameters**

    **coords** (*[**[**float**]**]*) – Coords of the vertices of the simplex. E.g.,
    [[1, 2, 3], [2, 4, 5], [6, 7, 8], [8, 9, 10].



#### bary_coords(point)

* **Parameters**

    **(****)** (*point*) – Point coordinates.



* **Returns**

    Barycentric coordinations.



#### _property_ coords()
Returns a copy of the vertex coordinates in the simplex.


#### in_simplex(point, tolerance=1e-08)
Checks if a point is in the simplex using the standard barycentric
coordinate system algorithm.

Taking an arbitrary vertex as an origin, we compute the basis for the
simplex from this origin by subtracting all other vertices from the
origin. We then project the point into this coordinate system and
determine the linear decomposition coefficients in this coordinate
system. If the coeffs satisfy that all coeffs >= 0, the composition
is in the facet.


* **Parameters**


    * **point** (*[**float**]*) – Point to test


    * **tolerance** (*float*) – Tolerance to test if point is in simplex.



#### line_intersection(point1, point2, tolerance=1e-08)
Computes the intersection points of a line with a simplex
:param point1: Points that determine the line
:type point1: [float]
:param point2: Points that determine the line
:type point2: [float]


* **Returns**

    points where the line intersects the simplex (0, 1, or 2).



#### point_from_bary_coords(bary_coords)

* **Parameters**

    **(****)** (*bary_coords*) – Barycentric coordinates.



* **Returns**

    Point coordinates



#### _property_ volume()
Volume of the simplex.


### pymatgen.util.coord.all_distances(coords1, coords2)
Returns the distances between two lists of coordinates.


* **Parameters**


    * **coords1** – First set of Cartesian coordinates.


    * **coords2** – Second set of Cartesian coordinates.



* **Returns**

    2d array of Cartesian distances. E.g the distance between
    coords1[i] and coords2[j] is distances[i,j]



### pymatgen.util.coord.barycentric_coords(coords, simplex)
Converts a list of coordinates to barycentric coordinates, given a
simplex with d+1 points. Only works for d >= 2.


* **Parameters**


    * **coords** – list of n coords to transform, shape should be (n,d)


    * **simplex** – list of coordinates that form the simplex, shape should be
    (d+1, d)



* **Returns**

    a LIST of barycentric coordinates (even if the original input was 1d)



### pymatgen.util.coord.coord_list_mapping(subset: ArrayLike, superset: ArrayLike, atol: float = 1e-08)
Gives the index mapping from a subset to a superset.
Subset and superset cannot contain duplicate rows.


* **Parameters**


    * **subset** (*ArrayLike*) – List of coords


    * **superset** (*ArrayLike*) – List of coords


    * **atol** (*float*) – Absolute tolerance. Defaults to 1e-8.



* **Returns**

    list of indices such that superset[indices] = subset



### pymatgen.util.coord.coord_list_mapping_pbc(subset, superset, atol=1e-08, pbc=(True, True, True))
Gives the index mapping from a subset to a superset.
Superset cannot contain duplicate matching rows.


* **Parameters**


    * **subset** (*ArrayLike*) – List of frac_coords


    * **superset** (*ArrayLike*) – List of frac_coords


    * **atol** (*float*) – Absolute tolerance. Defaults to 1e-8.


    * **pbc** (*tuple*) – A tuple defining the periodic boundary conditions along the three
    axis of the lattice.



* **Returns**

    list of indices such that superset[indices] = subset



### pymatgen.util.coord.find_in_coord_list(coord_list, coord, atol=1e-08)
Find the indices of matches of a particular coord in a coord_list.


* **Parameters**


    * **coord_list** – List of coords to test


    * **coord** – Specific coordinates


    * **atol** – Absolute tolerance. Defaults to 1e-8. Accepts both scalar and
    array.



* **Returns**

    Indices of matches, e.g., [0, 1, 2, 3]. Empty list if not found.



### pymatgen.util.coord.find_in_coord_list_pbc(fcoord_list, fcoord, atol=1e-08, pbc=(True, True, True))
Get the indices of all points in a fractional coord list that are
equal to a fractional coord (with a tolerance), taking into account
periodic boundary conditions.


* **Parameters**


    * **fcoord_list** – List of fractional coords


    * **fcoord** – A specific fractional coord to test.


    * **atol** – Absolute tolerance. Defaults to 1e-8.


    * **pbc** – a tuple defining the periodic boundary conditions along the three
    axis of the lattice.



* **Returns**

    Indices of matches, e.g., [0, 1, 2, 3]. Empty list if not found.



### pymatgen.util.coord.get_angle(v1, v2, units='degrees')
Calculates the angle between two vectors.


* **Parameters**


    * **v1** – Vector 1


    * **v2** – Vector 2


    * **units** – “degrees” or “radians”. Defaults to “degrees”.



* **Returns**

    Angle between them in degrees.



### pymatgen.util.coord.get_linear_interpolated_value(x_values, y_values, x)
Returns an interpolated value by linear interpolation between two values.
This method is written to avoid dependency on scipy, which causes issues on
threading servers.


* **Parameters**


    * **x_values** – Sequence of x values.


    * **y_values** – Corresponding sequence of y values


    * **x** – Get value at particular x



* **Returns**

    Value at x.



### pymatgen.util.coord.in_coord_list(coord_list, coord, atol=1e-08)
Tests if a particular coord is within a coord_list.


* **Parameters**


    * **coord_list** – List of coords to test


    * **coord** – Specific coordinates


    * **atol** – Absolute tolerance. Defaults to 1e-8. Accepts both scalar and
    array.



* **Returns**

    True if coord is in the coord list.



### pymatgen.util.coord.in_coord_list_pbc(fcoord_list, fcoord, atol=1e-08, pbc=(True, True, True))
Tests if a particular fractional coord is within a fractional coord_list.


* **Parameters**


    * **fcoord_list** – List of fractional coords to test


    * **fcoord** – A specific fractional coord to test.


    * **atol** – Absolute tolerance. Defaults to 1e-8.


    * **pbc** – a tuple defining the periodic boundary conditions along the three
    axis of the lattice.



* **Returns**

    True if coord is in the coord list.



### pymatgen.util.coord.is_coord_subset(subset: ArrayLike, superset: ArrayLike, atol: float = 1e-08)
Tests if all coords in subset are contained in superset.
Doesn’t use periodic boundary conditions.


* **Parameters**


    * **subset** (*ArrayLike*) – List of coords


    * **superset** (*ArrayLike*) – List of coords


    * **atol** (*float*) – Absolute tolerance for comparing coordinates. Defaults to 1e-8.



* **Returns**

    True if all of subset is in superset.



### pymatgen.util.coord.is_coord_subset_pbc(subset, superset, atol=1e-08, mask=None, pbc=(True, True, True))
Tests if all fractional coords in subset are contained in superset.


* **Parameters**


    * **subset** (*list*) – List of fractional coords to test


    * **superset** (*list*) – List of fractional coords to test against


    * **atol** (*float** or **size 3 array*) – Tolerance for matching


    * **mask** (*boolean array*) – Mask of matches that are not allowed.
    i.e. if mask[1,2] is True, then subset[1] cannot be matched
    to superset[2]


    * **pbc** (*tuple*) – a tuple defining the periodic boundary conditions along the three
    axis of the lattice.



* **Returns**

    True if all of subset is in superset.



### pymatgen.util.coord.lattice_points_in_supercell(supercell_matrix)
Returns the list of points on the original lattice contained in the
supercell in fractional coordinates (with the supercell basis).
e.g. [[2,0,0],[0,1,0],[0,0,1]] returns [[0,0,0],[0.5,0,0]].


* **Parameters**

    **supercell_matrix** – 3x3 matrix describing the supercell



* **Returns**

    numpy array of the fractional coordinates



### pymatgen.util.coord.pbc_diff(fcoords1: ArrayLike, fcoords2: ArrayLike, pbc: tuple[bool, bool, bool] = (True, True, True))
Returns the ‘fractional distance’ between two coordinates taking into
account periodic boundary conditions.


* **Parameters**


    * **fcoords1** – First set of fractional coordinates. e.g., [0.5, 0.6,
    0.7] or [[1.1, 1.2, 4.3], [0.5, 0.6, 0.7]]. It can be a single
    coord or any array of coords.


    * **fcoords2** – Second set of fractional coordinates.


    * **pbc** – a tuple defining the periodic boundary conditions along the three
    axis of the lattice.



* **Returns**

    Fractional distance. Each coordinate must have the property that
    abs(a) <= 0.5. Examples:
    pbc_diff([0.1, 0.1, 0.1], [0.3, 0.5, 0.9]) = [-0.2, -0.4, 0.2]
    pbc_diff([0.9, 0.1, 1.01], [0.3, 0.5, 0.9]) = [-0.4, -0.4, 0.11]



### pymatgen.util.coord.pbc_shortest_vectors(lattice, fcoords1, fcoords2, mask=None, return_d2=False)
Returns the shortest vectors between two lists of coordinates taking into
account periodic boundary conditions and the lattice.


* **Parameters**


    * **lattice** – lattice to use


    * **fcoords1** – First set of fractional coordinates. e.g., [0.5, 0.6, 0.7]
    or [[1.1, 1.2, 4.3], [0.5, 0.6, 0.7]]. It can be a single
    coord or any array of coords.


    * **fcoords2** – Second set of fractional coordinates.


    * **mask** (*boolean array*) – Mask of matches that are not allowed.
    i.e. if mask[1,2] is True, then subset[1] cannot be matched
    to superset[2]


    * **return_d2** (*bool*) – whether to also return the squared distances



* **Returns**

    array of displacement vectors from fcoords1 to fcoords2
    first index is fcoords1 index, second is fcoords2 index