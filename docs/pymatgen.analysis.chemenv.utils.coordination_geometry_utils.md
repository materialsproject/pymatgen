---
layout: default
title: pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md
nav_exclude: true
---

# pymatgen.analysis.chemenv.utils.coordination_geometry_utils module

This module contains some utility functions and classes that are used in the chemenv package.


### _class_ pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane(coefficients, p1=None, p2=None, p3=None)
Bases: `object`

Class used to describe a plane.

Initializes a plane from the 4 coefficients a, b, c and d of ax + by + cz + d = 0
:param coefficients: abcd coefficients of the plane.


#### TEST_2D_POINTS(_ = (array([0., 0.]), array([1., 0.]), array([0., 1.]), array([-1.,  0.]), array([ 0., -1.]), array([0., 2.]), array([2., 0.]), array([ 0., -2.]), array([-2.,  0.]), array([1., 1.]), array([2., 2.]), array([-1., -1.]), array([-2., -2.]), array([1., 2.]), array([ 1., -2.]), array([-1.,  2.]), array([-1., -2.]), array([2., 1.]), array([ 2., -1.]), array([-2.,  1.]), array([-2., -1.])_ )

#### _property_ a()
Coefficient a of the plane.


#### _property_ abcd()
Return a tuple with the plane coefficients.


* **Returns**

    Tuple with the plane coefficients.



#### _property_ b()
Coefficient b of the plane.


#### _property_ c()
Coefficient c of the plane.


#### _property_ coefficients()
Return a copy of the plane coefficients.


* **Returns**

    Plane coefficients as a numpy array.



#### _property_ crosses_origin()
Whether this plane crosses the origin (i.e. coefficient d is 0.0).


#### _property_ d()
Coefficient d of the plane.


#### _property_ distance_to_origin()
Distance of the plane to the origin.


#### distance_to_point(point)
Computes the absolute distance from the plane to the point
:param point: Point for which distance is computed
:return: Distance between the plane and the point.


#### distances(points)
Computes the distances from the plane to each of the points. Positive distances are on the side of the
normal of the plane while negative distances are on the other side
:param points: Points for which distances are computed
:return: Distances from the plane to the points (positive values on the side of the normal to the plane,

> negative values on the other side).


#### distances_indices_groups(points, delta=None, delta_factor=0.05, sign=False)
Computes the distances from the plane to each of the points. Positive distances are on the side of the
normal of the plane while negative distances are on the other side. Indices sorting the points from closest
to furthest is also computed. Grouped indices are also given, for which indices of the distances that are
separated by less than delta are grouped together. The delta parameter is either set explicitly or taken as
a fraction (using the delta_factor parameter) of the maximal point distance.
:param points: Points for which distances are computed
:param delta: Distance interval for which two points are considered in the same group.
:param delta_factor: If delta is None, the distance interval is taken as delta_factor times the maximal

> point distance.


* **Parameters**

    **sign** – Whether to add sign information in the indices sorting the points distances



* **Returns**

    Distances from the plane to the points (positive values on the side of the normal to the plane,
    negative values on the other side), as well as indices of the points from closest to furthest and
    grouped indices of distances separated by less than delta. For the sorting list and the grouped
    indices, when the sign parameter is True, items are given as tuples of (index, sign).



#### distances_indices_sorted(points, sign=False)
Computes the distances from the plane to each of the points. Positive distances are on the side of the
normal of the plane while negative distances are on the other side. Indices sorting the points from closest
to furthest is also computed.
:param points: Points for which distances are computed
:param sign: Whether to add sign information in the indices sorting the points distances
:return: Distances from the plane to the points (positive values on the side of the normal to the plane,

> negative values on the other side), as well as indices of the points from closest to furthest. For the
> latter, when the sign parameter is True, items of the sorting list are given as tuples of (index, sign).


#### fit_error(points, fit='least_square_distance')
Evaluate the error for a list of points with respect to this plane.


* **Parameters**


    * **points** – List of points.


    * **fit** – Type of fit error.



* **Returns**

    Error for a list of points with respect to this plane.



#### fit_least_square_distance_error(points)
Evaluate the sum of squared distances error for a list of points with respect to this plane.


* **Parameters**

    **points** – List of points.



* **Returns**

    Sum of squared distances error for a list of points with respect to this plane.



#### fit_maximum_distance_error(points)
Evaluate the max distance error for a list of points with respect to this plane.


* **Parameters**

    **points** – List of points.



* **Returns**

    Max distance error for a list of points with respect to this plane.



#### _classmethod_ from_2points_and_origin(p1, p2)
Initializes plane from two points and the origin.


* **Parameters**


    * **p1** – First point.


    * **p2** – Second point.



* **Returns**

    Plane.



#### _classmethod_ from_3points(p1, p2, p3)
Initializes plane from three points.


* **Parameters**


    * **p1** – First point.


    * **p2** – Second point.


    * **p3** – Third point.



* **Returns**

    Plane.



#### _classmethod_ from_coefficients(a, b, c, d)
Initialize plane from its coefficients.


* **Parameters**


    * **a** – a coefficient of the plane.


    * **b** – b coefficient of the plane.


    * **c** – c coefficient of the plane.


    * **d** – d coefficient of the plane.



* **Returns**

    Plane.



#### _classmethod_ from_npoints(points, best_fit='least_square_distance')
Initializes plane from a list of points.

If the number of points is larger than 3, will use a least square fitting or max distance fitting.


* **Parameters**


    * **points** – List of points.


    * **best_fit** – Type of fitting procedure for more than 3 points.



* **Returns**

    Plane



#### _classmethod_ from_npoints_least_square_distance(points)
Initializes plane from a list of points using a least square fitting procedure.


* **Parameters**

    **points** – List of points.



* **Returns**

    Plane.



#### _classmethod_ from_npoints_maximum_distance(points)
Initializes plane from a list of points using a max distance fitting procedure.


* **Parameters**

    **points** – List of points.



* **Returns**

    Plane.



#### indices_separate(points, dist_tolerance)
Returns three lists containing the indices of the points lying on one side of the plane, on the plane
and on the other side of the plane. The dist_tolerance parameter controls the tolerance to which a point
is considered to lie on the plane or not (distance to the plane)
:param points: list of points
:param dist_tolerance: tolerance to which a point is considered to lie on the plane

> or not (distance to the plane)


* **Returns**

    The lists of indices of the points on one side of the plane, on the plane and
    on the other side of the plane.



#### init_3points(non_zeros, zeros)
Initialize three random points on this plane.


* **Parameters**


    * **non_zeros** – Indices of plane coefficients ([a, b, c]) that are not zero.


    * **zeros** – Indices of plane coefficients ([a, b, c]) that are equal to zero.



* **Returns**

    None



#### is_in_list(plane_list)
Checks whether the plane is identical to one of the Planes in the plane_list list of Planes
:param plane_list: List of Planes to be compared to
:return: True if the plane is in the list, False otherwise.


#### is_in_plane(pp, dist_tolerance)
Determines if point pp is in the plane within the tolerance dist_tolerance
:param pp: point to be tested
:param dist_tolerance: tolerance on the distance to the plane within which point pp is considered in the plane
:return: True if pp is in the plane, False otherwise.


#### is_same_plane_as(plane)
Checks whether the plane is identical to another Plane “plane”
:param plane: Plane to be compared to
:return: True if the two facets are identical, False otherwise.


#### orthonormal_vectors()
Returns a list of three orthogonal vectors, the two first being parallel to the plane and the
third one is the normal vector of the plane
:return: List of orthogonal vectors
:raise: ValueError if all the coefficients are zero or if there is some other strange error.


#### _classmethod_ perpendicular_bisector(p1, p2)
Initialize a plane from the perpendicular bisector of two points.

The perpendicular bisector of two points is the plane perpendicular to the vector joining these two points
and passing through the middle of the segment joining the two points.


* **Parameters**


    * **p1** – First point.


    * **p2** – Second point.



* **Returns**

    Plane.



#### project_and_to2dim(pps, plane_center)
Projects the list of points pps to the plane and changes the basis from 3D to the 2D basis of the plane
:param pps: List of points to be projected
:return: :raise:


#### project_and_to2dim_ordered_indices(pps, plane_center='mean')
Projects each points in the point list pps on plane and returns the indices that would sort the
list of projected points in anticlockwise order
:param pps: List of points to project on plane
:return: List of indices that would sort the list of projected points.


#### projectionpoints(pps)
Projects each points in the point list pps on plane and returns the list of projected points
:param pps: List of points to project on plane
:return: List of projected point on plane.


### pymatgen.analysis.chemenv.utils.coordination_geometry_utils.anticlockwise_sort(pps)
Sort a list of 2D points in anticlockwise order
:param pps: List of points to be sorted
:return: Sorted list of points.


### pymatgen.analysis.chemenv.utils.coordination_geometry_utils.anticlockwise_sort_indices(pps)
Returns the indices that would sort a list of 2D points in anticlockwise order
:param pps: List of points to be sorted
:return: Indices of the sorted list of points.


### pymatgen.analysis.chemenv.utils.coordination_geometry_utils.changebasis(uu, vv, nn, pps)
For a list of points given in standard coordinates (in terms of e1, e2 and e3), returns the same list
expressed in the basis (uu, vv, nn), which is supposed to be orthonormal.
:param uu: First vector of the basis
:param vv: Second vector of the basis
:param nn: Third vector of the bais
:param pps: List of points in basis (e1, e2, e3)
:return: List of points in basis (uu, vv, nn).


### pymatgen.analysis.chemenv.utils.coordination_geometry_utils.collinear(p1, p2, p3=None, tolerance=0.25)
Checks if the three points p1, p2 and p3 are collinear or not within a given tolerance. The collinearity is
checked by computing the area of the triangle defined by the three points p1, p2 and p3. If the area of this
triangle is less than (tolerance x largest_triangle), then the three points are considered collinear. The
largest_triangle is defined as the right triangle whose legs are the two smallest distances between the three

> points ie, its area is : 0.5 x (min(

> ```
> |p2-p1|
> ```

> ,|p3-p1|,|p3-p2|) x secondmin(

> ```
> |p2-p1|
> ```

> ,|p3-p1|,|p3-p2|))


* **Parameters**


    * **p1** – First point


    * **p2** – Second point


    * **p3** – Third point (origin [0.0, 0.0, 0.0 if not given])


    * **tolerance** – Area tolerance for the collinearity test (0.25 gives about 0.125 deviation from the line)



* **Returns**

    True if the three points are considered as collinear within the given tolerance, False otherwise.



### pymatgen.analysis.chemenv.utils.coordination_geometry_utils.diamond_functions(xx, yy, y_x0, x_y0)
Method that creates two upper and lower functions based on points xx and yy
as well as intercepts defined by y_x0 and x_y0. The resulting functions
form kind of a distorted diamond-like structure aligned from
point xx to point yy.

Schematically :

xx is symbolized by x, yy is symbolized by y, y_x0 is equal to the distance
from x to a, x_y0 is equal to the distance from x to b, the lines a-p and
b-q are parallel to the line x-y such that points p and q are
obtained automatically.
In case of an increasing diamond the lower function is x-b-q and the upper
function is a-p-y while in case of a
decreasing diamond, the lower function is a-p-y and the upper function is
x-b-q.

> Increasing diamond      |     Decreasing diamond

>     > > > > p–y                    x—-b

>     > > > /  /|

>     > > > ```
>     > > > |
>     > > > ```



>     > > /  / |                    |    q

>     > /  /  |                    a    |

>     a  /   |                       |
>     | /    q                       |


>     ```
>     |
>     ```

>     /    /                         |
>     x—-b                          p–y


* **Parameters**


    * **xx** – First point


    * **yy** – Second point



* **Returns**

    A dictionary with the lower and upper diamond functions.



### pymatgen.analysis.chemenv.utils.coordination_geometry_utils.function_comparison(f1, f2, x1, x2, numpoints_check=500)
Method that compares two functions.


* **Parameters**


    * **f1** – First function to compare


    * **f2** – Second function to compare


    * **x1** – Lower bound of the interval to compare


    * **x2** – Upper bound of the interval to compare


    * **numpoints_check** – Number of points used to compare the functions



* **Returns**

    Whether the function are equal (“=”), f1 is always lower than f2 (“<”), f1 is always larger than f2 (“>”),

        f1 is always lower than or equal to f2 (“<”), f1 is always larger than or equal to f2 (“>”) on the
        interval [x1, x2]. If the two functions cross, a RuntimeError is thrown (i.e. we expect to compare
        functions that do not cross…)




### pymatgen.analysis.chemenv.utils.coordination_geometry_utils.get_lower_and_upper_f(surface_calculation_options)
Get the lower and upper functions defining a surface in the distance-angle space of neighbors.


* **Parameters**

    **surface_calculation_options** – Options for the surface.



* **Returns**

    Dictionary containing the “lower” and “upper” functions for the surface.



### pymatgen.analysis.chemenv.utils.coordination_geometry_utils.is_anion_cation_bond(valences, ii, jj)
Checks if two given sites are an anion and a cation.
:param valences: list of site valences
:param ii: index of a site
:param jj: index of another site
:return: True if one site is an anion and the other is a cation (from the valences).


### pymatgen.analysis.chemenv.utils.coordination_geometry_utils.matrixTimesVector(MM, aa)

* **Parameters**


    * **MM** – A matrix of size 3x3


    * **aa** – A vector of size 3



* **Returns**

    A vector of size 3 which is the product of the matrix by the vector



### pymatgen.analysis.chemenv.utils.coordination_geometry_utils.my_solid_angle(center, coords)
Helper method to calculate the solid angle of a set of coords from the
center.


* **Parameters**


    * **center** – Center to measure solid angle from.


    * **coords** – List of coords to determine solid angle.



* **Returns**

    The solid angle.



### pymatgen.analysis.chemenv.utils.coordination_geometry_utils.quarter_ellipsis_functions(xx, yy)
Method that creates two quarter-ellipse functions based on points xx and yy. The ellipsis is supposed to
be aligned with the axes. The two ellipsis pass through the two points xx and yy.


* **Parameters**


    * **xx** – First point


    * **yy** – Second point



* **Returns**

    A dictionary with the lower and upper quarter ellipsis functions.



### pymatgen.analysis.chemenv.utils.coordination_geometry_utils.rectangle_surface_intersection(rectangle, f_lower, f_upper, bounds_lower=None, bounds_upper=None, check=True, numpoints_check=500)
Method to calculate the surface of the intersection of a rectangle (aligned with axes) and another surface
defined by two functions f_lower and f_upper.


* **Parameters**


    * **rectangle** – Rectangle defined as : ((x1, x2), (y1, y2)).


    * **f_lower** – Function defining the lower bound of the surface.


    * **f_upper** – Function defining the upper bound of the surface.


    * **bounds_lower** – Interval in which the f_lower function is defined.


    * **bounds_upper** – Interval in which the f_upper function is defined.


    * **check** – Whether to check if f_lower is always lower than f_upper.


    * **numpoints_check** – Number of points used to check whether f_lower is always lower than f_upper



* **Returns**

    The surface of the intersection of the rectangle and the surface defined by f_lower and f_upper.



### pymatgen.analysis.chemenv.utils.coordination_geometry_utils.rotateCoords(coords, R)
Rotate the list of points using rotation matrix R
:param coords: List of points to be rotated
:param R: Rotation matrix
:return: List of rotated points.


### pymatgen.analysis.chemenv.utils.coordination_geometry_utils.rotateCoordsOpt(coords, R)
Rotate the list of points using rotation matrix R
:param coords: List of points to be rotated
:param R: Rotation matrix
:return: List of rotated points.


### pymatgen.analysis.chemenv.utils.coordination_geometry_utils.separation_in_list(separation_indices, separation_indices_list)
Checks if the separation indices of a plane are already in the list
:param separation_indices: list of separation indices (three arrays of integers)
:param separation_indices_list: list of the list of separation indices to be compared to
:return: True if the separation indices are already in the list, False otherwise.


### pymatgen.analysis.chemenv.utils.coordination_geometry_utils.sort_separation(separation)
Sort a separation.


* **Parameters**

    **separation** – Initial separation.



* **Returns**

    Sorted list of separation.



### pymatgen.analysis.chemenv.utils.coordination_geometry_utils.sort_separation_tuple(separation)
Sort a separation.


* **Parameters**

    **separation** – Initial separation



* **Returns**

    Sorted tuple of separation



### pymatgen.analysis.chemenv.utils.coordination_geometry_utils.spline_functions(lower_points, upper_points, degree=3)
Method that creates two (upper and lower) spline functions based on points lower_points and upper_points.


* **Parameters**


    * **lower_points** – Points defining the lower function.


    * **upper_points** – Points defining the upper function.


    * **degree** – Degree for the spline function



* **Returns**

    A dictionary with the lower and upper spline functions.



### pymatgen.analysis.chemenv.utils.coordination_geometry_utils.vectorsToMatrix(aa, bb)
Performs the vector multiplication of the elements of two vectors, constructing the 3x3 matrix.
:param aa: One vector of size 3
:param bb: Another vector of size 3
:return: A 3x3 matrix M composed of the products of the elements of aa and bb :

> M_ij = aa_i \* bb_j.