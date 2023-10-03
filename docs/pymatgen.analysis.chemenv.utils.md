---
layout: default
title: pymatgen.analysis.chemenv.utils.md
nav_exclude: true
---

1. TOC
{:toc}

# pymatgen.analysis.chemenv.utils package

Utility package for chemenv.


## pymatgen.analysis.chemenv.utils.chemenv_config module

This module contains the classes for configuration of the chemenv package.


### _class_ ChemEnvConfig(package_options=None)
Bases: `object`

Class used to store the configuration of the chemenv package :


    * Materials project access


    * ICSD database access


    * Default options (strategies, …).


* **Parameters**

    **package_options** –



#### DEFAULT_PACKAGE_OPTIONS(_ = {'default_max_distance_factor': 1.5, 'default_strategy': {'strategy': 'SimplestChemenvStrategy', 'strategy_options': {'additional_condition': 1, 'angle_cutoff': 0.3, 'continuous_symmetry_measure_cutoff': 10, 'distance_cutoff': 1.4}}_ )

#### _classmethod_ auto_load(root_dir=None)
Autoload options.
:param root_dir:


#### _property_ has_materials_project_access()
Whether MP access is enabled.


#### package_options_description()
Describe package options.


#### save(root_dir=None)
Save the options.
:param root_dir:


#### setup()
Setup the class.


#### setup_package_options()
Setup the package options.

## pymatgen.analysis.chemenv.utils.chemenv_errors module

This module contains the error classes for the chemenv package.


### _exception_ AbstractChemenvError(cls, method, msg)
Bases: `Exception`

Abstract class for Chemenv errors.


* **Parameters**


    * **cls** –


    * **method** –


    * **msg** –



### _exception_ ChemenvError(cls, method, msg)
Bases: `Exception`

Chemenv error.


* **Parameters**


    * **cls** –


    * **method** –


    * **msg** –



### _exception_ EquivalentSiteSearchError(site)
Bases: `AbstractChemenvError`

Equivalent site search error.


* **Parameters**

    **site** –



### _exception_ NeighborsNotComputedChemenvError(site)
Bases: `AbstractChemenvError`

Neighbors not computed error.


* **Parameters**

    **site** –



### _exception_ SolidAngleError(cosinus)
Bases: `AbstractChemenvError`

Solid angle error.


* **Parameters**

    **cosinus** –


## pymatgen.analysis.chemenv.utils.coordination_geometry_utils module

This module contains some utility functions and classes that are used in the chemenv package.


### _class_ Plane(coefficients, p1=None, p2=None, p3=None)
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


* **Returns**

    Distance between the plane and the point.



#### distances(points)
Computes the distances from the plane to each of the points. Positive distances are on the side of the
normal of the plane while negative distances are on the other side
:param points: Points for which distances are computed


* **Returns**

    Distances from the plane to the points (positive values on the side of the normal to the plane,
    negative values on the other side).



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


* **Returns**

    Distances from the plane to the points (positive values on the side of the normal to the plane,
    negative values on the other side), as well as indices of the points from closest to furthest. For the
    latter, when the sign parameter is True, items of the sorting list are given as tuples of (index, sign).



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



#### is_in_list(plane_list)
Checks whether the plane is identical to one of the Planes in the plane_list list of Planes
:param plane_list: List of Planes to be compared to


* **Returns**

    True if the plane is in the list.



* **Return type**

    bool



#### is_in_plane(pp, dist_tolerance)
Determines if point pp is in the plane within the tolerance dist_tolerance
:param pp: point to be tested
:param dist_tolerance: tolerance on the distance to the plane within which point pp is considered in the plane


* **Returns**

    True if pp is in the plane.



* **Return type**

    bool



#### is_same_plane_as(plane)
Checks whether the plane is identical to another Plane “plane”
:param plane: Plane to be compared to


* **Returns**

    True if the two facets are identical.



* **Return type**

    bool



#### orthonormal_vectors()
Returns a list of three orthogonal vectors, the two first being parallel to the plane and the
third one is the normal vector of the plane


* **Returns**

    List of orthogonal vectors



* **Raise**

    ValueError if all the coefficients are zero or if there is some other strange error.



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


* **Returns**

    raise:



#### project_and_to2dim_ordered_indices(pps, plane_center='mean')
Projects each points in the point list pps on plane and returns the indices that would sort the
list of projected points in anticlockwise order
:param pps: List of points to project on plane


* **Returns**

    List of indices that would sort the list of projected points.



#### projectionpoints(pps)
Projects each points in the point list pps on plane and returns the list of projected points
:param pps: List of points to project on plane


* **Returns**

    List of projected point on plane.



### anticlockwise_sort(pps)
Sort a list of 2D points in anticlockwise order
:param pps: List of points to be sorted


* **Returns**

    Sorted list of points.



### anticlockwise_sort_indices(pps)
Returns the indices that would sort a list of 2D points in anticlockwise order
:param pps: List of points to be sorted


* **Returns**

    Indices of the sorted list of points.



### changebasis(uu, vv, nn, pps)
For a list of points given in standard coordinates (in terms of e1, e2 and e3), returns the same list
expressed in the basis (uu, vv, nn), which is supposed to be orthonormal.
:param uu: First vector of the basis
:param vv: Second vector of the basis
:param nn: Third vector of the bais
:param pps: List of points in basis (e1, e2, e3)
:returns: List of points in basis (uu, vv, nn).


### collinear(p1, p2, p3=None, tolerance=0.25)
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

    True if the three points are considered as collinear within the given tolerance.



* **Return type**

    bool



### diamond_functions(xx, yy, y_x0, x_y0)
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



### function_comparison(f1, f2, x1, x2, numpoints_check=500)
Method that compares two functions.


* **Parameters**


    * **f1** – First function to compare


    * **f2** – Second function to compare


    * **x1** – Lower bound of the interval to compare


    * **x2** – Upper bound of the interval to compare


    * **numpoints_check** – Number of points used to compare the functions



* **Returns**

    ‘=’ if the functions are equal, ‘<’ if f1 is always lower than f2, ‘>’ if f1 is always larger than f2,

        f1 is always lower than or equal to f2 (“<”), f1 is always larger than or equal to f2 (“>”) on the
        interval [x1, x2]. If the two functions cross, a RuntimeError is thrown (i.e. we expect to compare
        functions that do not cross…)




* **Return type**

    str



### get_lower_and_upper_f(surface_calculation_options)
Get the lower and upper functions defining a surface in the distance-angle space of neighbors.


* **Parameters**

    **surface_calculation_options** – Options for the surface.



* **Returns**

    Dictionary containing the “lower” and “upper” functions for the surface.



### is_anion_cation_bond(valences, ii, jj)
Checks if two given sites are an anion and a cation.
:param valences: list of site valences
:param ii: index of a site
:param jj: index of another site


* **Returns**

    True if one site is an anion and the other is a cation (based on valences).



* **Return type**

    bool



### matrixTimesVector(MM, aa)

* **Parameters**


    * **MM** – A matrix of size 3x3


    * **aa** – A vector of size 3



* **Returns**

    A vector of size 3 which is the product of the matrix by the vector



### quarter_ellipsis_functions(xx: ArrayLike, yy: ArrayLike)
Method that creates two quarter-ellipse functions based on points xx and yy. The ellipsis is supposed to
be aligned with the axes. The two ellipsis pass through the two points xx and yy.


* **Parameters**


    * **xx** – First point


    * **yy** – Second point



* **Returns**

    A dictionary with the lower and upper quarter ellipsis functions.



### rectangle_surface_intersection(rectangle, f_lower, f_upper, bounds_lower=None, bounds_upper=None, check=True, numpoints_check=500)
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



### rotateCoords(coords, R)
Rotate the list of points using rotation matrix R
:param coords: List of points to be rotated
:param R: Rotation matrix


* **Returns**

    List of rotated points.



### rotateCoordsOpt(coords, R)
Rotate the list of points using rotation matrix R
:param coords: List of points to be rotated
:param R: Rotation matrix


* **Returns**

    List of rotated points.



### separation_in_list(separation_indices, separation_indices_list)
Checks if the separation indices of a plane are already in the list
:param separation_indices: list of separation indices (three arrays of integers)
:param separation_indices_list: list of the list of separation indices to be compared to


* **Returns**

    True if the separation indices are already in the list.



* **Return type**

    bool



### solid_angle(center, coords)
Helper method to calculate the solid angle of a set of coords from the center.


* **Parameters**


    * **center** – Center to measure solid angle from.


    * **coords** – List of coords to determine solid angle.



* **Returns**

    The solid angle.



### sort_separation(separation)
Sort a separation.


* **Parameters**

    **separation** – Initial separation.



* **Returns**

    Sorted list of separation.



### sort_separation_tuple(separation)
Sort a separation.


* **Parameters**

    **separation** – Initial separation



* **Returns**

    Sorted tuple of separation



### spline_functions(lower_points, upper_points, degree=3)
Method that creates two (upper and lower) spline functions based on points lower_points and upper_points.


* **Parameters**


    * **lower_points** – Points defining the lower function.


    * **upper_points** – Points defining the upper function.


    * **degree** – Degree for the spline function



* **Returns**

    A dictionary with the lower and upper spline functions.



### vectorsToMatrix(aa, bb)
Performs the vector multiplication of the elements of two vectors, constructing the 3x3 matrix.
:param aa: One vector of size 3
:param bb: Another vector of size 3


* **Returns**

    M_ij = aa_i \* bb_j.



* **Return type**

    A 3x3 matrix M composed of the products of the elements of aa and bb


## pymatgen.analysis.chemenv.utils.defs_utils module

This module contains the definition of some objects used in the chemenv package.


### _class_ AdditionalConditions()
Bases: `object`

Class for additional conditions.


#### ALL(_ = (0, 1, 2, 3, 4_ )

#### CONDITION_DESCRIPTION(_ = {0: 'No additional condition', 1: 'Only anion-cation bonds', 2: 'No element-element bonds (same elements)', 3: 'Only anion-cation bonds and no element-element bonds (same elements)', 4: 'Only element-oxygen bonds'_ )

#### NONE(_ = _ )

#### NO_AC(_ = _ )

#### NO_ADDITIONAL_CONDITION(_ = _ )

#### NO_E2SEB(_ = _ )

#### NO_ELEMENT_TO_SAME_ELEMENT_BONDS(_ = _ )

#### ONLY_ACB(_ = _ )

#### ONLY_ACB_AND_NO_E2SEB(_ = _ )

#### ONLY_ANION_CATION_BONDS(_ = _ )

#### ONLY_ANION_CATION_BONDS_AND_NO_ELEMENT_TO_SAME_ELEMENT_BONDS(_ = _ )

#### ONLY_E2OB(_ = _ )

#### ONLY_ELEMENT_TO_OXYGEN_BONDS(_ = _ )

#### check_condition(condition, structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), parameters)

* **Parameters**


    * **condition** –


    * **structure** –


    * **parameters** –


## pymatgen.analysis.chemenv.utils.func_utils module

This module contains some utility functions and classes that are used in the chemenv package.


### _class_ AbstractRatioFunction(function, options_dict=None)
Bases: `object`

Abstract class for all ratio functions.

Constructor for AbstractRatioFunction.


* **Parameters**


    * **function** – Ration function name.


    * **options_dict** – Dictionary containing the parameters for the ratio function.



#### ALLOWED_FUNCTIONS(_: ClassVar[dict[str, list]_ _ = {_ )

#### evaluate(value)
Evaluate the ratio function for the given value.


* **Parameters**

    **value** – Value for which ratio function has to be evaluated.



* **Returns**

    Ratio function corresponding to the value.



#### _classmethod_ from_dict(dct)
Construct ratio function from dict.


* **Parameters**

    **dct** – Dict representation of the ratio function



* **Returns**

    Ratio function object.



#### setup_parameters(options_dict)
Set up the parameters for this ratio function.


* **Parameters**

    **options_dict** – Dictionary containing the parameters for the ratio function.



### _class_ CSMFiniteRatioFunction(function, options_dict=None)
Bases: `AbstractRatioFunction`

Concrete implementation of a series of ratio functions applied to the continuous symmetry measure (CSM).

Uses “finite” ratio functions.

See the following reference for details:
ChemEnv: a fast and robust coordination environment identification tool,
D. Waroquiers et al., Acta Cryst. B 76, 683 (2020).

Constructor for AbstractRatioFunction.


* **Parameters**


    * **function** – Ration function name.


    * **options_dict** – Dictionary containing the parameters for the ratio function.



#### ALLOWED_FUNCTIONS(_: ClassVar[dict[str, list]_ _ = {'power2_decreasing_exp': ['max_csm', 'alpha'], 'smootherstep': ['lower_csm', 'upper_csm'], 'smoothstep': ['lower_csm', 'upper_csm']_ )

#### fractions(data)
Get the fractions from the CSM ratio function applied to the data.


* **Parameters**

    **data** – List of CSM values to estimate fractions.



* **Returns**

    Corresponding fractions for each CSM.



#### mean_estimator(data)
Get the weighted CSM using this CSM ratio function applied to the data.


* **Parameters**

    **data** – List of CSM values to estimate the weighted CSM.



* **Returns**

    Weighted CSM from this ratio function.



#### power2_decreasing_exp(vals)
Get the evaluation of the ratio function f(x)=exp(-a\*x)\*(x-1)^2.

The CSM values (i.e. “x”), are scaled to the “max_csm” parameter. The “a” constant
correspond to the “alpha” parameter.


* **Parameters**

    **vals** – CSM values for which the ratio function has to be evaluated.



* **Returns**

    Result of the ratio function applied to the CSM values.



#### ratios(data)
Get the fractions from the CSM ratio function applied to the data.


* **Parameters**

    **data** – List of CSM values to estimate fractions.



* **Returns**

    Corresponding fractions for each CSM.



#### smootherstep(vals)
Get the evaluation of the smootherstep ratio function: f(x)=6\*x^5-15\*x^4+10\*x^3.

The CSM values (i.e. “x”), are scaled between the “lower_csm” and “upper_csm” parameters.


* **Parameters**

    **vals** – CSM values for which the ratio function has to be evaluated.



* **Returns**

    Result of the ratio function applied to the CSM values.



#### smoothstep(vals)
Get the evaluation of the smoothstep ratio function: f(x)=3\*x^2-2\*x^3.

The CSM values (i.e. “x”), are scaled between the “lower_csm” and “upper_csm” parameters.


* **Parameters**

    **vals** – CSM values for which the ratio function has to be evaluated.



* **Returns**

    Result of the ratio function applied to the CSM values.



### _class_ CSMInfiniteRatioFunction(function, options_dict=None)
Bases: `AbstractRatioFunction`

Concrete implementation of a series of ratio functions applied to the continuous symmetry measure (CSM).

Uses “infinite” ratio functions.

See the following reference for details:
ChemEnv: a fast and robust coordination environment identification tool,
D. Waroquiers et al., Acta Cryst. B 76, 683 (2020).

Constructor for AbstractRatioFunction.


* **Parameters**


    * **function** – Ration function name.


    * **options_dict** – Dictionary containing the parameters for the ratio function.



#### ALLOWED_FUNCTIONS(_: ClassVar[dict[str, list]_ _ = {'power2_inverse_decreasing': ['max_csm'], 'power2_inverse_power2_decreasing': ['max_csm']_ )

#### fractions(data)
Get the fractions from the CSM ratio function applied to the data.


* **Parameters**

    **data** – List of CSM values to estimate fractions.



* **Returns**

    Corresponding fractions for each CSM.



#### mean_estimator(data)
Get the weighted CSM using this CSM ratio function applied to the data.


* **Parameters**

    **data** – List of CSM values to estimate the weighted CSM.



* **Returns**

    Weighted CSM from this ratio function.



#### power2_inverse_decreasing(vals)
Get the evaluation of the ratio function f(x)=(x-1)^2 / x.

The CSM values (i.e. “x”), are scaled to the “max_csm” parameter. The “a” constant
correspond to the “alpha” parameter.


* **Parameters**

    **vals** – CSM values for which the ratio function has to be evaluated.



* **Returns**

    Result of the ratio function applied to the CSM values.



#### power2_inverse_power2_decreasing(vals)
Get the evaluation of the ratio function f(x)=(x-1)^2 / x^2.

The CSM values (i.e. “x”), are scaled to the “max_csm” parameter. The “a” constant
correspond to the “alpha” parameter.


* **Parameters**

    **vals** – CSM values for which the ratio function has to be evaluated.



* **Returns**

    Result of the ratio function applied to the CSM values.



#### ratios(data)
Get the fractions from the CSM ratio function applied to the data.


* **Parameters**

    **data** – List of CSM values to estimate fractions.



* **Returns**

    Corresponding fractions for each CSM.



### _class_ DeltaCSMRatioFunction(function, options_dict=None)
Bases: `AbstractRatioFunction`

Concrete implementation of a series of ratio functions applied to differences of
continuous symmetry measures (DeltaCSM).

Uses “finite” ratio functions.

See the following reference for details:
ChemEnv: a fast and robust coordination environment identification tool,
D. Waroquiers et al., Acta Cryst. B 76, 683 (2020).

Constructor for AbstractRatioFunction.


* **Parameters**


    * **function** – Ration function name.


    * **options_dict** – Dictionary containing the parameters for the ratio function.



#### ALLOWED_FUNCTIONS(_: ClassVar[dict[str, list]_ _ = {'smootherstep': ['delta_csm_min', 'delta_csm_max']_ )

#### smootherstep(vals)
Get the evaluation of the smootherstep ratio function: f(x)=6\*x^5-15\*x^4+10\*x^3.

The DeltaCSM values (i.e. “x”), are scaled between the “delta_csm_min” and “delta_csm_max” parameters.


* **Parameters**

    **vals** – DeltaCSM values for which the ratio function has to be evaluated.



* **Returns**

    Result of the ratio function applied to the DeltaCSM values.



### _class_ RatioFunction(function, options_dict=None)
Bases: `AbstractRatioFunction`

Concrete implementation of a series of ratio functions.

Constructor for AbstractRatioFunction.


* **Parameters**


    * **function** – Ration function name.


    * **options_dict** – Dictionary containing the parameters for the ratio function.



#### ALLOWED_FUNCTIONS(_: ClassVar[dict[str, list]_ _ = {'inverse_smootherstep': ['lower', 'upper'], 'inverse_smoothstep': ['lower', 'upper'], 'power2_decreasing_exp': ['max', 'alpha'], 'power2_inverse_decreasing': ['max'], 'power2_inverse_power2_decreasing': ['max'], 'smootherstep': ['lower', 'upper'], 'smoothstep': ['lower', 'upper']_ )

#### inverse_smootherstep(vals)
Get the evaluation of the “inverse” smootherstep ratio function: f(x)=1-(6\*x^5-15\*x^4+10\*x^3).

The values (i.e. “x”), are scaled between the “lower” and “upper” parameters.


* **Parameters**

    **vals** – Values for which the ratio function has to be evaluated.



* **Returns**

    Result of the ratio function applied to the values.



#### inverse_smoothstep(vals)
Get the evaluation of the “inverse” smoothstep ratio function: f(x)=1-(3\*x^2-2\*x^3).

The values (i.e. “x”), are scaled between the “lower” and “upper” parameters.


* **Parameters**

    **vals** – Values for which the ratio function has to be evaluated.



* **Returns**

    Result of the ratio function applied to the values.



#### power2_decreasing_exp(vals)
Get the evaluation of the ratio function f(x)=exp(-a\*x)\*(x-1)^2.

The values (i.e. “x”), are scaled to the “max” parameter. The “a” constant
correspond to the “alpha” parameter.


* **Parameters**

    **vals** – Values for which the ratio function has to be evaluated.



* **Returns**

    Result of the ratio function applied to the values.



#### power2_inverse_decreasing(vals)
Get the evaluation of the ratio function f(x)=(x-1)^2 / x.

The values (i.e. “x”), are scaled to the “max” parameter.


* **Parameters**

    **vals** – Values for which the ratio function has to be evaluated.



* **Returns**

    Result of the ratio function applied to the values.



#### power2_inverse_power2_decreasing(vals)
Get the evaluation of the ratio function f(x)=(x-1)^2 / x^2.

The values (i.e. “x”), are scaled to the “max” parameter.


* **Parameters**

    **vals** – Values for which the ratio function has to be evaluated.



* **Returns**

    Result of the ratio function applied to the values.



#### smootherstep(vals)
Get the evaluation of the smootherstep ratio function: f(x)=6\*x^5-15\*x^4+10\*x^3.

The values (i.e. “x”), are scaled between the “lower” and “upper” parameters.


* **Parameters**

    **vals** – Values for which the ratio function has to be evaluated.



* **Returns**

    Result of the ratio function applied to the values.



#### smoothstep(vals)
Get the evaluation of the smoothstep ratio function: f(x)=3\*x^2-2\*x^3.

The values (i.e. “x”), are scaled between the “lower” and “upper” parameters.


* **Parameters**

    **vals** – Values for which the ratio function has to be evaluated.



* **Returns**

    Result of the ratio function applied to the values.


## pymatgen.analysis.chemenv.utils.graph_utils module

This module contains some graph utils that are used in the chemenv package.


### _class_ MultiGraphCycle(nodes, edge_indices, validate=True, ordered=None)
Bases: `MSONable`

Class used to describe a cycle in a multigraph.

nodes are the nodes of the cycle and edge_indices are the indices of the edges in the cycle.
The nth index in edge_indices corresponds to the edge index between the nth node in nodes and
the (n+1)th node in nodes with the exception of the last one being the edge index between
the last node in nodes and the first node in nodes

Example: A cycle

    nodes:          1 - 3 - 4 - 0 - 2 - (1)
    edge_indices:     0 . 1 . 0 . 2 . 0 . (0)


* **Parameters**


    * **nodes** –


    * **edge_indices** –


    * **validate** –


    * **ordered** –



#### _is_valid(check_strict_ordering=False)
Check if a MultiGraphCycle is valid.

This method checks that:
1. there are no duplicate nodes,
2. there are either 1 or more than 2 nodes


* **Returns**

    True if the SimpleGraphCycle is valid.



* **Return type**

    bool



#### order(raise_on_fail=True)
Orders the SimpleGraphCycle.

The ordering is performed such that the first node is the “lowest” one
and the second node is the lowest one of the two neighbor nodes of the
first node. If raise_on_fail is set to True a RuntimeError will be
raised if the ordering fails.


* **Parameters**

    **raise_on_fail** – If set to True, will raise a RuntimeError if the ordering fails.



#### validate(check_strict_ordering=False)

* **Parameters**

    **check_strict_ordering** –



### _class_ SimpleGraphCycle(nodes, validate=True, ordered=None)
Bases: `MSONable`

Class used to describe a cycle in a simple graph (graph without multiple edges).

Note that the convention used here is the networkx convention for which simple graphs allow
to have self-loops in a simple graph.
No simple graph cycle with two nodes is possible in a simple graph. The graph will not
be validated if validate is set to False.
By default, the “ordered” parameter is None, in which case the SimpleGraphCycle will be ordered.
If the user explicitly sets ordered to False, the SimpleGraphCycle will not be ordered.


* **Parameters**


    * **nodes** –


    * **validate** –


    * **ordered** –



#### _is_valid(check_strict_ordering=False)
Check if a SimpleGraphCycle is valid.

This method checks :
- that there are no duplicate nodes,
- that there are either 1 or more than 2 nodes


* **Returns**

    True if the SimpleGraphCycle is valid.



* **Return type**

    bool



#### as_dict()
MSONable dict


#### _classmethod_ from_dict(d, validate=False)
Serialize from dict.
:param d:
:param validate:


#### _classmethod_ from_edges(edges, edges_are_ordered=True)
Constructs SimpleGraphCycle from a list edges.

By default, the edges list is supposed to be ordered as it will be
much faster to construct the cycle. If edges_are_ordered is set to
False, the code will automatically try to find the corresponding edge
order in the list.


#### order(raise_on_fail=True)
Orders the SimpleGraphCycle.

The ordering is performed such that the first node is the “lowest” one and the
second node is the lowest one of the two neighbor nodes of the first node. If
raise_on_fail is set to True a RuntimeError will be raised if the ordering fails.


* **Parameters**

    **raise_on_fail** (*bool*) – If set to True, will raise a RuntimeError if the ordering fails.



#### validate(check_strict_ordering=False)

* **Parameters**

    **check_strict_ordering** –



### _c2index_isreverse(c1, c2)
Private helper function to get the index c2_0_index of the first node of cycle c1
in cycle c2 and whether the cycle c2 should be reversed or not.

Returns None if the first node of cycle c1 is not found in cycle c2.
The reverse value depends on the index c2_1_index of the second node of cycle c1 in
cycle c2 : if it is *just after* the c2_0_index, reverse is False, if it is
*just before* the c2_0_index, reverse is True, otherwise the function returns None).


### get_all_elementary_cycles(graph)

* **Parameters**

    **graph** –



### get_all_simple_paths_edges(graph, source, target, cutoff=None, data=True)
Get all the simple path and edges.
:param graph:
:param source:
:param target:
:param cutoff:
:param data:


### get_delta(node1, node2, edge_data)
Get the delta.
:param node1:
:param node2:
:param edge_data:

## pymatgen.analysis.chemenv.utils.math_utils module

This module contains some math utils that are used in the chemenv package.


### _append_es2sequences(sequences, es)

### _cartesian_product(lists)
given a list of lists,
returns all the possible combinations taking one element from each list
The list does not have to be of equal length.


### _factor_generator(n)
From a given natural integer, returns the prime factors and their multiplicity
:param n: Natural integer


### cosinus_step(xx, edges=None, inverse=False)

* **Parameters**


    * **xx** –


    * **edges** –


    * **inverse** –



### divisors(n)
From a given natural integer, returns the list of divisors in ascending order
:param n: Natural integer


* **Returns**

    List of divisors of n in ascending order.



### get_center_of_arc(p1, p2, radius)

* **Parameters**


    * **p1** –


    * **p2** –


    * **radius** –



### get_linearly_independent_vectors(vectors_list)

* **Parameters**

    **vectors_list** –



### normal_cdf_step(xx, mean, scale)

* **Parameters**


    * **xx** –


    * **mean** –


    * **scale** –



### power2_decreasing_exp(xx, edges=None, alpha=1.0)

* **Parameters**


    * **xx** –


    * **edges** –


    * **alpha** –



### power2_inverse_decreasing(xx, edges=None, prefactor=None)

* **Parameters**


    * **xx** –


    * **edges** –


    * **prefactor** –



### power2_inverse_power2_decreasing(xx, edges=None, prefactor=None)

* **Parameters**


    * **xx** –


    * **edges** –


    * **prefactor** –



### power2_inverse_powern_decreasing(xx, edges=None, prefactor=None, powern=2.0)

* **Parameters**


    * **xx** –


    * **edges** –


    * **prefactor** –


    * **powern** –



### power2_tangent_decreasing(xx, edges=None, prefactor=None)

* **Parameters**


    * **xx** –


    * **edges** –


    * **prefactor** –



### power3_step(xx, edges=None, inverse=False)

* **Parameters**


    * **xx** –


    * **edges** –


    * **inverse** –



### powern_decreasing(xx, edges=None, nn=2)

* **Parameters**


    * **xx** –


    * **edges** –


    * **nn** –



### powern_parts_step(xx, edges=None, inverse=False, nn=2)

* **Parameters**


    * **xx** –


    * **edges** –


    * **inverse** –


    * **nn** –



### prime_factors(n: int)
Lists prime factors of a given natural integer, from greatest to smallest
:param n: Natural integer
:rtype : list of all prime factors of the given natural n.


### scale_and_clamp(xx, edge0, edge1, clamp0, clamp1)

* **Parameters**


    * **xx** –


    * **edge0** –


    * **edge1** –


    * **clamp0** –


    * **clamp1** –



### smootherstep(xx, edges=None, inverse=False)

* **Parameters**


    * **xx** –


    * **edges** –


    * **inverse** –



### smoothstep(xx, edges=None, inverse=False)

* **Parameters**


    * **xx** –


    * **edges** –


    * **inverse** –


## pymatgen.analysis.chemenv.utils.scripts_utils module

This module contains some script utils that are used in the chemenv package.


### compute_environments(chemenv_configuration)
Compute the environments.


* **Parameters**

    **chemenv_configuration** –



### draw_cg(vis, site, neighbors, cg=None, perm=None, perfect2local_map=None, show_perfect=False, csm_info=None, symmetry_measure_type='csm_wcs_ctwcc', perfect_radius=0.1, show_distorted=True, faces_color_override=None)
Draw cg.


* **Parameters**


    * **vis** –


    * **site** –


    * **neighbors** –


    * **cg** –


    * **perm** –


    * **perfect2local_map** –


    * **show_perfect** –


    * **csm_info** –


    * **symmetry_measure_type** –


    * **perfect_radius** –


    * **show_distorted** –


    * **faces_color_override** –



### visualize(cg, zoom=None, vis=None, factor=1.0, view_index=True, faces_color_override=None)
Visualizing a coordination geometry
:param cg:
:param zoom:
:param vis:
:param factor:
:param view_index:
:param faces_color_override: