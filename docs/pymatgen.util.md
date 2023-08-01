---
layout: default
title: pymatgen.util.md
nav_exclude: true
---

# pymatgen.util package

The util package implements various utilities that are commonly used by various
packages.

## Subpackages


* [pymatgen.util.tests package](pymatgen.util.tests.md)




    * [pymatgen.util.tests.test_convergence module](pymatgen.util.tests.md#module-pymatgen.util.tests.test_convergence)


        * [`ConvergenceTest`](pymatgen.util.tests.md#pymatgen.util.tests.test_convergence.ConvergenceTest)


            * [`ConvergenceTest.test_determine_convergence()`](pymatgen.util.tests.md#pymatgen.util.tests.test_convergence.ConvergenceTest.test_determine_convergence)


    * [pymatgen.util.tests.test_coord module](pymatgen.util.tests.md#module-pymatgen.util.tests.test_coord)


        * [`CoordUtilsTest`](pymatgen.util.tests.md#pymatgen.util.tests.test_coord.CoordUtilsTest)


            * [`CoordUtilsTest.test_all_distances()`](pymatgen.util.tests.md#pymatgen.util.tests.test_coord.CoordUtilsTest.test_all_distances)


            * [`CoordUtilsTest.test_barycentric()`](pymatgen.util.tests.md#pymatgen.util.tests.test_coord.CoordUtilsTest.test_barycentric)


            * [`CoordUtilsTest.test_coord_list_mapping()`](pymatgen.util.tests.md#pymatgen.util.tests.test_coord.CoordUtilsTest.test_coord_list_mapping)


            * [`CoordUtilsTest.test_coord_list_mapping_pbc()`](pymatgen.util.tests.md#pymatgen.util.tests.test_coord.CoordUtilsTest.test_coord_list_mapping_pbc)


            * [`CoordUtilsTest.test_find_in_coord_list()`](pymatgen.util.tests.md#pymatgen.util.tests.test_coord.CoordUtilsTest.test_find_in_coord_list)


            * [`CoordUtilsTest.test_find_in_coord_list_pbc()`](pymatgen.util.tests.md#pymatgen.util.tests.test_coord.CoordUtilsTest.test_find_in_coord_list_pbc)


            * [`CoordUtilsTest.test_get_angle()`](pymatgen.util.tests.md#pymatgen.util.tests.test_coord.CoordUtilsTest.test_get_angle)


            * [`CoordUtilsTest.test_get_linear_interpolated_value()`](pymatgen.util.tests.md#pymatgen.util.tests.test_coord.CoordUtilsTest.test_get_linear_interpolated_value)


            * [`CoordUtilsTest.test_in_coord_list()`](pymatgen.util.tests.md#pymatgen.util.tests.test_coord.CoordUtilsTest.test_in_coord_list)


            * [`CoordUtilsTest.test_in_coord_list_pbc()`](pymatgen.util.tests.md#pymatgen.util.tests.test_coord.CoordUtilsTest.test_in_coord_list_pbc)


            * [`CoordUtilsTest.test_is_coord_subset()`](pymatgen.util.tests.md#pymatgen.util.tests.test_coord.CoordUtilsTest.test_is_coord_subset)


            * [`CoordUtilsTest.test_is_coord_subset_pbc()`](pymatgen.util.tests.md#pymatgen.util.tests.test_coord.CoordUtilsTest.test_is_coord_subset_pbc)


            * [`CoordUtilsTest.test_lattice_points_in_supercell()`](pymatgen.util.tests.md#pymatgen.util.tests.test_coord.CoordUtilsTest.test_lattice_points_in_supercell)


            * [`CoordUtilsTest.test_pbc_diff()`](pymatgen.util.tests.md#pymatgen.util.tests.test_coord.CoordUtilsTest.test_pbc_diff)


            * [`CoordUtilsTest.test_pbc_shortest_vectors()`](pymatgen.util.tests.md#pymatgen.util.tests.test_coord.CoordUtilsTest.test_pbc_shortest_vectors)


        * [`SimplexTest`](pymatgen.util.tests.md#pymatgen.util.tests.test_coord.SimplexTest)


            * [`SimplexTest.setUp()`](pymatgen.util.tests.md#pymatgen.util.tests.test_coord.SimplexTest.setUp)


            * [`SimplexTest.test_2dtriangle()`](pymatgen.util.tests.md#pymatgen.util.tests.test_coord.SimplexTest.test_2dtriangle)


            * [`SimplexTest.test_bary_coords()`](pymatgen.util.tests.md#pymatgen.util.tests.test_coord.SimplexTest.test_bary_coords)


            * [`SimplexTest.test_equal()`](pymatgen.util.tests.md#pymatgen.util.tests.test_coord.SimplexTest.test_equal)


            * [`SimplexTest.test_in_simplex()`](pymatgen.util.tests.md#pymatgen.util.tests.test_coord.SimplexTest.test_in_simplex)


            * [`SimplexTest.test_intersection()`](pymatgen.util.tests.md#pymatgen.util.tests.test_coord.SimplexTest.test_intersection)


            * [`SimplexTest.test_str()`](pymatgen.util.tests.md#pymatgen.util.tests.test_coord.SimplexTest.test_str)


            * [`SimplexTest.test_to_json()`](pymatgen.util.tests.md#pymatgen.util.tests.test_coord.SimplexTest.test_to_json)


            * [`SimplexTest.test_volume()`](pymatgen.util.tests.md#pymatgen.util.tests.test_coord.SimplexTest.test_volume)


    * [pymatgen.util.tests.test_graph_hashing module](pymatgen.util.tests.md#module-pymatgen.util.tests.test_graph_hashing)


        * [`test_graph_hash()`](pymatgen.util.tests.md#pymatgen.util.tests.test_graph_hashing.test_graph_hash)


        * [`test_subgraph_hashes()`](pymatgen.util.tests.md#pymatgen.util.tests.test_graph_hashing.test_subgraph_hashes)


    * [pymatgen.util.tests.test_io_utils module](pymatgen.util.tests.md#module-pymatgen.util.tests.test_io_utils)


        * [`FuncTest`](pymatgen.util.tests.md#pymatgen.util.tests.test_io_utils.FuncTest)


            * [`FuncTest.test_micro_pyawk()`](pymatgen.util.tests.md#pymatgen.util.tests.test_io_utils.FuncTest.test_micro_pyawk)


    * [pymatgen.util.tests.test_num_utils module](pymatgen.util.tests.md#module-pymatgen.util.tests.test_num_utils)


        * [`FuncTestCase`](pymatgen.util.tests.md#pymatgen.util.tests.test_num_utils.FuncTestCase)


            * [`FuncTestCase.test_abs_cap()`](pymatgen.util.tests.md#pymatgen.util.tests.test_num_utils.FuncTestCase.test_abs_cap)


            * [`FuncTestCase.test_min_max_indexes()`](pymatgen.util.tests.md#pymatgen.util.tests.test_num_utils.FuncTestCase.test_min_max_indexes)


            * [`FuncTestCase.test_round()`](pymatgen.util.tests.md#pymatgen.util.tests.test_num_utils.FuncTestCase.test_round)


    * [pymatgen.util.tests.test_plotting module](pymatgen.util.tests.md#module-pymatgen.util.tests.test_plotting)


        * [`FuncTestCase`](pymatgen.util.tests.md#pymatgen.util.tests.test_plotting.FuncTestCase)


            * [`FuncTestCase.test_plot_periodic_heatmap()`](pymatgen.util.tests.md#pymatgen.util.tests.test_plotting.FuncTestCase.test_plot_periodic_heatmap)


            * [`FuncTestCase.test_van_arkel_triangle()`](pymatgen.util.tests.md#pymatgen.util.tests.test_plotting.FuncTestCase.test_van_arkel_triangle)


    * [pymatgen.util.tests.test_provenance module](pymatgen.util.tests.md#module-pymatgen.util.tests.test_provenance)


        * [`StructureNLCase`](pymatgen.util.tests.md#pymatgen.util.tests.test_provenance.StructureNLCase)


            * [`StructureNLCase.setUp()`](pymatgen.util.tests.md#pymatgen.util.tests.test_provenance.StructureNLCase.setUp)


            * [`StructureNLCase.test_authors()`](pymatgen.util.tests.md#pymatgen.util.tests.test_provenance.StructureNLCase.test_authors)


            * [`StructureNLCase.test_data()`](pymatgen.util.tests.md#pymatgen.util.tests.test_provenance.StructureNLCase.test_data)


            * [`StructureNLCase.test_eq()`](pymatgen.util.tests.md#pymatgen.util.tests.test_provenance.StructureNLCase.test_eq)


            * [`StructureNLCase.test_from_structures()`](pymatgen.util.tests.md#pymatgen.util.tests.test_provenance.StructureNLCase.test_from_structures)


            * [`StructureNLCase.test_history_nodes()`](pymatgen.util.tests.md#pymatgen.util.tests.test_provenance.StructureNLCase.test_history_nodes)


            * [`StructureNLCase.test_references()`](pymatgen.util.tests.md#pymatgen.util.tests.test_provenance.StructureNLCase.test_references)


            * [`StructureNLCase.test_remarks()`](pymatgen.util.tests.md#pymatgen.util.tests.test_provenance.StructureNLCase.test_remarks)


            * [`StructureNLCase.test_to_from_dict()`](pymatgen.util.tests.md#pymatgen.util.tests.test_provenance.StructureNLCase.test_to_from_dict)


    * [pymatgen.util.tests.test_string module](pymatgen.util.tests.md#module-pymatgen.util.tests.test_string)


        * [`FuncTest`](pymatgen.util.tests.md#pymatgen.util.tests.test_string.FuncTest)


            * [`FuncTest.test_charge_string()`](pymatgen.util.tests.md#pymatgen.util.tests.test_string.FuncTest.test_charge_string)


            * [`FuncTest.test_disordered_formula()`](pymatgen.util.tests.md#pymatgen.util.tests.test_string.FuncTest.test_disordered_formula)


            * [`FuncTest.test_formula_double_format()`](pymatgen.util.tests.md#pymatgen.util.tests.test_string.FuncTest.test_formula_double_format)


            * [`FuncTest.test_htmlify()`](pymatgen.util.tests.md#pymatgen.util.tests.test_string.FuncTest.test_htmlify)


            * [`FuncTest.test_latexify()`](pymatgen.util.tests.md#pymatgen.util.tests.test_string.FuncTest.test_latexify)


            * [`FuncTest.test_latexify_spacegroup()`](pymatgen.util.tests.md#pymatgen.util.tests.test_string.FuncTest.test_latexify_spacegroup)


            * [`FuncTest.test_transformation_to_string()`](pymatgen.util.tests.md#pymatgen.util.tests.test_string.FuncTest.test_transformation_to_string)


            * [`FuncTest.test_unicodeify()`](pymatgen.util.tests.md#pymatgen.util.tests.test_string.FuncTest.test_unicodeify)


        * [`StringifyTest`](pymatgen.util.tests.md#pymatgen.util.tests.test_string.StringifyTest)


            * [`StringifyTest.test_to_html_string()`](pymatgen.util.tests.md#pymatgen.util.tests.test_string.StringifyTest.test_to_html_string)


            * [`StringifyTest.test_to_latex_string()`](pymatgen.util.tests.md#pymatgen.util.tests.test_string.StringifyTest.test_to_latex_string)


            * [`StringifyTest.test_to_unicode_string()`](pymatgen.util.tests.md#pymatgen.util.tests.test_string.StringifyTest.test_to_unicode_string)


        * [`SubStr`](pymatgen.util.tests.md#pymatgen.util.tests.test_string.SubStr)


        * [`SupStr`](pymatgen.util.tests.md#pymatgen.util.tests.test_string.SupStr)


            * [`SupStr.STRING_MODE`](pymatgen.util.tests.md#pymatgen.util.tests.test_string.SupStr.STRING_MODE)


            * [`SupStr.to_pretty_string()`](pymatgen.util.tests.md#pymatgen.util.tests.test_string.SupStr.to_pretty_string)


    * [pymatgen.util.tests.test_typing module](pymatgen.util.tests.md#module-pymatgen.util.tests.test_typing)


        * [`test_composition_like()`](pymatgen.util.tests.md#pymatgen.util.tests.test_typing.test_composition_like)


        * [`test_entry_like()`](pymatgen.util.tests.md#pymatgen.util.tests.test_typing.test_entry_like)


        * [`test_matrix_like()`](pymatgen.util.tests.md#pymatgen.util.tests.test_typing.test_matrix_like)


        * [`test_path_like()`](pymatgen.util.tests.md#pymatgen.util.tests.test_typing.test_path_like)


        * [`test_species_like()`](pymatgen.util.tests.md#pymatgen.util.tests.test_typing.test_species_like)



## pymatgen.util.convergence module

Functions for calculating the convergence of an x, y data set.

Main API:

test_conv(xs, ys, name, tol)

tries to fit multiple functions to the x, y data

calculates which function fits best
for tol < 0
returns the x value for which y is converged within tol of the asymptotic value
for tol > 0
returns the x_value for which dy(x)/dx < tol for all x >= x_value, conv is true is such a x_value exists
for the best fit a gnuplot line is printed plotting the data, the function and the asymptotic value


### _exception_ pymatgen.util.convergence.SplineInputError(msg)
Bases: `Exception`

Error for Spline input.


* **Parameters**

    **msg** (*str*) – Message.



### pymatgen.util.convergence.determine_convergence(xs, ys, name, tol: float = 0.0001, extra='', verbose=False, mode='extra', plots=True)
Test it and at which x_value dy(x)/dx < tol for all x >= x_value, conv is true is such a x_value exists.


### pymatgen.util.convergence.exponential(x, a, b, n)
Exponential function base n to fit convergence data.


### pymatgen.util.convergence.extrapolate_reciprocal(xs, ys, n, noise)
Return the parameters such that a + b / x^n hits the last two data points.


### pymatgen.util.convergence.extrapolate_simple_reciprocal(xs, ys)
Extrapolate simple reciprocal function to fit convergence data.


* **Parameters**


    * **xs** – List of x values.


    * **ys** – List of y values.



* **Returns**

    List of parameters [a, b].



### pymatgen.util.convergence.get_derivatives(xs, ys, fd=False)
return the derivatives of y(x) at the points x
if scipy is available a spline is generated to calculate the derivatives
if scipy is not available the left and right slopes are calculated, if both exist the average is returned
putting fd to zero always returns the finite difference slopes.


### pymatgen.util.convergence.get_weights(xs, ys, mode=2)

* **Parameters**


    * **xs** (*list*) – List of x values.


    * **ys** (*list*) – List of y values.


    * **mode** (*int*) – Mode for calculating weights.



* **Returns**

    List of weights.



* **Return type**

    list



### pymatgen.util.convergence.id_generator(size: int = 8, chars: str = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789')
Generate a random string of specified size and characters.


* **Parameters**


    * **size** (*int*) – The length of the generated string.


    * **chars** (*str*) – The characters to use for generating the string.



* **Returns**

    The generated random string.



* **Return type**

    str



### pymatgen.util.convergence.measure(function, xs, ys, popt, weights)
Measure the quality of a fit.


### pymatgen.util.convergence.multi_curve_fit(xs, ys, verbose)
Fit multiple functions to the x, y data, return the best fit.


### pymatgen.util.convergence.multi_reciprocal_extra(xs, ys, noise=False)
Calculates for a series of powers ns the parameters for which the last two points are at the curve.
With these parameters measure how well the other data points fit.
return the best fit.


### pymatgen.util.convergence.p0_exponential(xs, ys)
Calculate the initial guess parameters for the exponential function.


* **Parameters**


    * **xs** (*list*) – List of x values.


    * **ys** (*list*) – List of y values.



* **Returns**

    List of initial guess parameters [a, b, n].



* **Return type**

    list



### pymatgen.util.convergence.p0_reciprocal(xs, ys)
Predictor for first guess for reciprocal.


### pymatgen.util.convergence.p0_simple_2reciprocal(xs, ys)
Calculate the initial guess parameters for the simple reciprocal function with a power of 2.


* **Parameters**


    * **xs** (*list*) – List of x values.


    * **ys** (*list*) – List of y values.



* **Returns**

    List of initial guess parameters [a, b].



* **Return type**

    list



### pymatgen.util.convergence.p0_simple_4reciprocal(xs, ys)
Calculate the initial guess parameters for the simple reciprocal function with a power of 4.


* **Parameters**


    * **xs** (*list*) – List of x values.


    * **ys** (*list*) – List of y values.



* **Returns**

    The initial guess parameters [a, b].



* **Return type**

    list



### pymatgen.util.convergence.p0_simple_5reciprocal(xs, ys)
Calculate the initial guess parameters for the simple reciprocal function with a power of 0.5.


* **Parameters**


    * **xs** (*list*) – List of x values.


    * **ys** (*list*) – List of y values.



* **Returns**

    List of parameters [a, b].



* **Return type**

    list



### pymatgen.util.convergence.p0_simple_reciprocal(xs, ys)
Calculate the initial guess parameters for the simple reciprocal function.


* **Parameters**


    * **xs** (*list*) – List of x values.


    * **ys** (*list*) – List of y values.



* **Returns**

    List of initial guess parameters [a, b].



* **Return type**

    list



### pymatgen.util.convergence.p0_single_reciprocal(xs, ys)
Calculate the initial guess parameters for the single reciprocal function.


* **Parameters**


    * **xs** (*list*) – List of x values.


    * **ys** (*list*) – List of y values.



* **Returns**

    List of initial guess parameters [a, b, c].



* **Return type**

    list



### pymatgen.util.convergence.print_and_raise_error(xs, ys, name)
Print error message and raise a RuntimeError.


* **Parameters**


    * **xs** (*list*) – List of x values.


    * **ys** (*list*) – List of y values.


    * **name** (*str*) – Name of the function where the error occurred.



### pymatgen.util.convergence.print_plot_line(function, popt, xs, ys, name, tol: float = 0.05, extra='')
Print the gnuplot command line to plot the x, y data with the fitted function using the popt parameters.


### pymatgen.util.convergence.reciprocal(x, a, b, n)
Reciprocal function to the power n to fit convergence data.


### pymatgen.util.convergence.simple_2reciprocal(x, a, b)
Reciprocal function to fit convergence data.


### pymatgen.util.convergence.simple_4reciprocal(x, a, b)
Reciprocal function to fit convergence data.


### pymatgen.util.convergence.simple_5reciprocal(x, a, b)
Reciprocal function to fit convergence data.


### pymatgen.util.convergence.simple_reciprocal(x, a, b)
Reciprocal function to fit convergence data.


### pymatgen.util.convergence.single_reciprocal(x, a, b, c)
Reciprocal function to fit convergence data.

## pymatgen.util.coord module

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


## pymatgen.util.coord_cython module

Utilities for manipulating coordinates or list of coordinates, under periodic
boundary conditions or otherwise.


### pymatgen.util.coord_cython.coord_list_mapping_pbc(subset, superset, atol=1e-08, pbc=(True, True, True))
Gives the index mapping from a subset to a superset.
Superset cannot contain duplicate matching rows


* **Parameters**


    * **subset** – List of frac_coords


    * **superset** – List of frac_coords


    * **pbc** – a tuple defining the periodic boundary conditions along the three
    axis of the lattice.



* **Returns**

    list of indices such that superset[indices] = subset



### pymatgen.util.coord_cython.is_coord_subset_pbc(subset, superset, atol, mask, pbc=(True, True, True))
Tests if all fractional coords in subset are contained in superset.
Allows specification of a mask determining pairs that are not
allowed to match to each other


* **Parameters**


    * **subset** – List of fractional coords


    * **superset** – List of fractional coords


    * **pbc** – a tuple defining the periodic boundary conditions along the three
    axis of the lattice.



* **Returns**

    True if all of subset is in superset.



### pymatgen.util.coord_cython.pbc_shortest_vectors(lattice, fcoords1, fcoords2, mask=None, return_d2=False, lll_frac_tol=None)
Returns the shortest vectors between two lists of coordinates taking into
account periodic boundary conditions and the lattice.


* **Parameters**


    * **lattice** – lattice to use


    * **fcoords1** – First set of fractional coordinates. e.g., [0.5, 0.6, 0.7]
    or [[1.1, 1.2, 4.3], [0.5, 0.6, 0.7]]. Must be

    ```
    np.float_
    ```



    * **fcoords2** – Second set of fractional coordinates.


    * **mask** (

    ```
    int_
    ```

     array) – Mask of matches that are not allowed.
    i.e. if mask[1,2] == True, then subset[1] cannot be matched
    to superset[2]


    * **lll_frac_tol** (

    ```
    float_
    ```

     array of length 3) – Fractional tolerance (per LLL lattice vector) over which
    the calculation of minimum vectors will be skipped.
    Can speed up calculation considerably for large structures.



* **Returns**

    array of displacement vectors from fcoords1 to fcoords2
    first index is fcoords1 index, second is fcoords2 index


## pymatgen.util.due module

Stub file for a guaranteed safe import of duecredit constructs: if duecredit
is not available.

Then use in your code as

> from .due import due, Doi, BibTeX, Text

See  [https://github.com/duecredit/duecredit/blob/master/README.md](https://github.com/duecredit/duecredit/blob/master/README.md) for examples.

Origin:     Originally a part of the duecredit
Copyright:  2015-2021  DueCredit developers
License:    BSD-2


### pymatgen.util.due.BibTeX(\*args, \*\*kwargs)
Perform no good and no bad.


### pymatgen.util.due.Doi(\*args, \*\*kwargs)
Perform no good and no bad.


### _class_ pymatgen.util.due.InactiveDueCreditCollector()
Bases: `object`

Just a stub at the Collector which would not do anything.


#### activate(\*args, \*\*kwargs)
Perform no good and no bad.


#### active(_ = Fals_ )

#### add(\*args, \*\*kwargs)
Perform no good and no bad.


#### cite(\*args, \*\*kwargs)
Perform no good and no bad.


#### dcite(\*args, \*\*kwargs)
If I could cite I would.


#### dump(\*args, \*\*kwargs)
Perform no good and no bad.


#### load(\*args, \*\*kwargs)
Perform no good and no bad.


### pymatgen.util.due.Text(\*args, \*\*kwargs)
Perform no good and no bad.


### pymatgen.util.due.Url(\*args, \*\*kwargs)
Perform no good and no bad.

## pymatgen.util.graph_hashing module

Copyright (C) 2004-2022, NetworkX Developers
Aric Hagberg <[hagberg@lanl.gov](mailto:hagberg@lanl.gov)>
Dan Schult <[dschult@colgate.edu](mailto:dschult@colgate.edu)>
Pieter Swart <[swart@lanl.gov](mailto:swart@lanl.gov)>
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

>
> * Redistributions of source code must retain the above copyright
> notice, this list of conditions and the following disclaimer.


> * Redistributions in binary form must reproduce the above
> copyright notice, this list of conditions and the following
> disclaimer in the documentation and/or other materials provided
> with the distribution.


> * Neither the name of the NetworkX Developers nor the names of its
> contributors may be used to endorse or promote products derived
> from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
“AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


---

Functions for hashing graphs to strings.
Isomorphic graphs should be assigned identical hashes.
For now, only Weisfeiler-Lehman hashing is implemented.


### pymatgen.util.graph_hashing.weisfeiler_lehman_graph_hash(G, edge_attr=None, node_attr=None, iterations=3, digest_size=16)
Return Weisfeiler Lehman (WL) graph hash.

The function iteratively aggregates and hashes neighborhoods of each node.
After each node’s neighbors are hashed to obtain updated node labels,
a hashed histogram of resulting labels is returned as the final hash.

Hashes are identical for isomorphic graphs and strong guarantees that
non-isomorphic graphs will get different hashes. See

```
[1]_
```

 for details.

If no node or edge attributes are provided, the degree of each node
is used as its initial label.
Otherwise, node and/or edge labels are used to compute the hash.


* **Parameters**


    * **G** – graph
    The graph to be hashed.
    Can have node and/or edge attributes. Can also have no attributes.


    * **edge_attr** – string, default=None
    The key in edge attribute dictionary to be used for hashing.
    If None, edge labels are ignored.


    * **node_attr** – string, default=None
    The key in node attribute dictionary to be used for hashing.
    If None, and no edge_attr given, use the degrees of the nodes as labels.


    * **iterations** – int, default=3
    Number of neighbor aggregations to perform.
    Should be larger for larger graphs.


    * **digest_size** – int, default=16
    Size (in bits) of blake2b hash digest to use for hashing node labels.



* **Returns**

    string

        Hexadecimal string corresponding to hash of the input graph.




* **Return type**

    h


### Notes

To return the WL hashes of each subgraph of a graph, use
weisfeiler_lehman_subgraph_hashes

Similarity between hashes does not imply similarity between graphs.

### References

Kurt Mehlhorn, and Karsten M. Borgwardt. Weisfeiler Lehman
Graph Kernels. Journal of Machine Learning Research. 2011.
[http://www.jmlr.org/papers/volume12/shervashidze11a/shervashidze11a.pdf](http://www.jmlr.org/papers/volume12/shervashidze11a/shervashidze11a.pdf)


### pymatgen.util.graph_hashing.weisfeiler_lehman_subgraph_hashes(G, edge_attr=None, node_attr=None, iterations=3, digest_size=16)
Return a dictionary of subgraph hashes by node.

The dictionary is keyed by node to a list of hashes in increasingly
sized induced subgraphs containing the nodes within 2\*k edges
of the key node for increasing integer k until all nodes are included.

The function iteratively aggregates and hashes neighborhoods of each node.
This is achieved for each step by replacing for each node its label from
the previous iteration with its hashed 1-hop neighborhood aggregate.
The new node label is then appended to a list of node labels for each
node.

To aggregate neighborhoods at each step for a node $n$, all labels of
nodes adjacent to $n$ are concatenated. If the edge_attr parameter is set,
labels for each neighboring node are prefixed with the value of this attribute
along the connecting edge from this neighbor to node $n$. The resulting string
is then hashed to compress this information into a fixed digest size.

Thus, at the $i$th iteration nodes within $2i$ distance influence any given
hashed node label. We can therefore say that at depth $i$ for node $n$
we have a hash for a subgraph induced by the $2i$-hop neighborhood of $n$.

Can be used to to create general Weisfeiler-Lehman graph kernels, or
generate features for graphs or nodes, for example to generate ‘words’ in a
graph as seen in the ‘graph2vec’ algorithm.
See

```
[1]_
```

 &

```
[2]_
```

 respectively for details.

Hashes are identical for isomorphic subgraphs and there exist strong
guarantees that non-isomorphic graphs will get different hashes.
See

```
[1]_
```

 for details.

If no node or edge attributes are provided, the degree of each node
is used as its initial label.
Otherwise, node and/or edge labels are used to compute the hash.


* **Parameters**


    * **G** – graph
    The graph to be hashed.
    Can have node and/or edge attributes. Can also have no attributes.


    * **edge_attr** – string, default=None
    The key in edge attribute dictionary to be used for hashing.
    If None, edge labels are ignored.


    * **node_attr** – string, default=None
    The key in node attribute dictionary to be used for hashing.
    If None, and no edge_attr given, use the degrees of the nodes as labels.


    * **iterations** – int, default=3
    Number of neighbor aggregations to perform.
    Should be larger for larger graphs.


    * **digest_size** – int, default=16
    Size (in bits) of blake2b hash digest to use for hashing node labels.
    The default size is 16 bits



* **Returns**

    dict

        A dictionary with each key given by a node in G, and each value given
        by the subgraph hashes in order of depth from the key node.




* **Return type**

    node_subgraph_hashes


### Notes

To hash the full graph when subgraph hashes are not needed, use
weisfeiler_lehman_graph_hash for efficiency.

Similarity between hashes does not imply similarity between graphs.

### References

Kurt Mehlhorn, and Karsten M. Borgwardt. Weisfeiler Lehman
Graph Kernels. Journal of Machine Learning Research. 2011.
[http://www.jmlr.org/papers/volume12/shervashidze11a/shervashidze11a.pdf](http://www.jmlr.org/papers/volume12/shervashidze11a/shervashidze11a.pdf)
.. [2] Annamalai Narayanan, Mahinthan Chandramohan, Rajasekar Venkatesan,
Lihui Chen, Yang Liu and Shantanu Jaiswa. graph2vec: Learning
Distributed Representations of Graphs. arXiv. 2017
[https://arxiv.org/pdf/1707.05005.pdf](https://arxiv.org/pdf/1707.05005.pdf)

## pymatgen.util.io_utils module

This module provides utility classes for io operations.


### pymatgen.util.io_utils.clean_lines(string_list, remove_empty_lines=True)
Strips whitespace, carriage returns and empty lines from a list of strings.


* **Parameters**


    * **string_list** – List of strings


    * **remove_empty_lines** – Set to True to skip lines which are empty after
    stripping.



* **Returns**

    List of clean strings with no whitespaces.



### pymatgen.util.io_utils.micro_pyawk(filename, search, results=None, debug=None, postdebug=None)
Small awk-mimicking search routine.

‘file’ is file to search through.
‘search’ is the “search program”, a list of lists/tuples with 3 elements;
i.e. [[regex,test,run],[regex,test,run],…]
‘results’ is a an object that your search program will have access to for
storing results.

Here regex is either as a Regex object, or a string that we compile into a
Regex. test and run are callable objects.

This function goes through each line in filename, and if regex matches that
line *and* test(results,line)==True (or test is None) we execute
run(results,match),where match is the match object from running
Regex.match.

The default results is an empty dictionary. Passing a results object let
you interact with it in run() and test(). Hence, in many occasions it is
thus clever to use results=self.

Author: Rickard Armiento, Ioannis Petousis


* **Returns**

    results


## pymatgen.util.num module

This module provides utilities for basic math operations.


### pymatgen.util.num.abs_cap(val, max_abs_val=1)
Returns the value with its absolute value capped at max_abs_val.
Particularly useful in passing values to trigonometric functions where
numerical errors may result in an argument > 1 being passed in.


* **Parameters**


    * **val** (*float*) – Input value.


    * **max_abs_val** (*float*) – The maximum absolute value for val. Defaults to 1.



* **Returns**

    val if abs(val) < 1 else sign of val \* max_abs_val.



### pymatgen.util.num.make_symmetric_matrix_from_upper_tri(val)
Given a symmetric matrix in upper triangular matrix form as flat array indexes as:
[A_xx,A_yy,A_zz,A_xy,A_xz,A_yz]
This will generate the full matrix:
[[A_xx,A_xy,A_xz],[A_xy,A_yy,A_yz],[A_xz,A_yz,A_zz].


### pymatgen.util.num.maxloc(seq)
Return the index of the (first) maximum in seq.

```python
>>> assert maxloc([1,3,2,3]) == 1
```


### pymatgen.util.num.min_max_indexes(seq)
Uses enumerate, max, and min to return the indices of the values
in a list with the maximum and minimum value.


* **Parameters**

    **seq** – A sequence of numbers.



### pymatgen.util.num.minloc(seq)
Return the index of the (first) minimum in seq.

```python
>>> assert minloc(range(3)) == 0
```


### pymatgen.util.num.monotonic(values, mode='<', atol=1e-08)
True if values are monotonically (decreasing|increasing).
mode is “<” for a decreasing sequence, “>” for an increasing sequence.
Two numbers are considered equal if they differ less than atol.

<!-- warning:
Not very efficient for large data sets. -->
```python
>>> values = [1.2, 1.3, 1.4]
>>> monotonic(values, mode="<")
False
>>> monotonic(values, mode=">")
True
```


### pymatgen.util.num.non_decreasing(values)
True if values are not decreasing.


### pymatgen.util.num.non_increasing(values)
True if values are not increasing.


### pymatgen.util.num.round_to_sigfigs(num, sig_figs)
Rounds a number rounded to a specific number of significant
figures instead of to a specific precision.


### pymatgen.util.num.strictly_decreasing(values)
True if values are strictly decreasing.


### pymatgen.util.num.strictly_increasing(values)
True if values are strictly increasing.

## pymatgen.util.numba module

This module provides a wrapper for numba such that no functionality
is lost if numba is not available. Numba is a just-in-time compiler
that can significantly accelerate the evaluation of certain functions
if installed.


### pymatgen.util.numba.jit(func)
Replacement for numba.jit when numba is not installed that does nothing.


### pymatgen.util.numba.njit(func)
Replacement for numba.njit when numba is not installed that does nothing.

## pymatgen.util.plotting module

Utilities for generating nicer plots.


### pymatgen.util.plotting.add_fig_kwargs(func)
Decorator that adds keyword arguments for functions returning matplotlib
figures.

The function should return either a matplotlib figure or None to signal
some sort of error/unexpected event.
See doc string below for the list of supported options.


### pymatgen.util.plotting.format_formula(formula)
Converts str of chemical formula into
latex format for labelling purposes.


* **Parameters**

    **formula** (*str*) – Chemical formula



### pymatgen.util.plotting.get_ax3d_fig_plt(ax=None, \*\*kwargs)
Helper function used in plot functions supporting an optional Axes3D
argument. If ax is None, we build the matplotlib figure and create the
Axes3D else we return the current active figure.


* **Parameters**


    * **ax** (*Axes3D**, **optional*) – Axes3D object. Defaults to None.


    * **kwargs** – keyword arguments are passed to plt.figure if ax is not None.



* **Returns**

    `Axes` object
    figure: matplotlib figure
    plt: matplotlib pyplot module.



* **Return type**

    ax



### pymatgen.util.plotting.get_ax_fig_plt(ax=None, \*\*kwargs)
Helper function used in plot functions supporting an optional Axes argument.
If ax is None, we build the matplotlib figure and create the Axes else
we return the current active figure.


* **Parameters**


    * **ax** (*Axes**, **optional*) – Axes object. Defaults to None.


    * **kwargs** – keyword arguments are passed to plt.figure if ax is not None.



* **Returns**

    `Axes` object
    figure: matplotlib figure
    plt: matplotlib pyplot module.



* **Return type**

    ax



### pymatgen.util.plotting.get_axarray_fig_plt(ax_array, nrows=1, ncols=1, sharex=False, sharey=False, squeeze=True, subplot_kw=None, gridspec_kw=None, \*\*fig_kw)
Helper function used in plot functions that accept an optional array of Axes
as argument. If ax_array is None, we build the matplotlib figure and
create the array of Axes by calling plt.subplots else we return the
current active figure.


* **Returns**

    Array of `Axes` objects
    figure: matplotlib figure
    plt: matplotlib pyplot module.



* **Return type**

    ax



### pymatgen.util.plotting.periodic_table_heatmap(elemental_data=None, cbar_label='', cbar_label_size=14, show_plot=False, cmap='YlOrRd', cmap_range=None, blank_color='grey', edge_color='white', value_format=None, value_fontsize=10, symbol_fontsize=14, max_row=9, readable_fontcolor=False, pymatviz: bool = True, \*\*kwargs)
A static method that generates a heat map overlaid on a periodic table.


* **Parameters**


    * **elemental_data** (*dict*) – A dictionary with the element as a key and a
    value assigned to it, e.g. surface energy and frequency, etc.
    Elements missing in the elemental_data will be grey by default
    in the final table elemental_data={“Fe”: 4.2, “O”: 5.0}.


    * **cbar_label** (*str*) – Label of the color bar. Default is “”.


    * **cbar_label_size** (*float*) – Font size for the color bar label. Default is 14.


    * **cmap_range** (*tuple*) – Minimum and maximum value of the color map scale.
    If None, the color map will automatically scale to the range of the
    data.


    * **show_plot** (*bool*) – Whether to show the heatmap. Default is False.


    * **value_format** (*str*) – Formatting string to show values. If None, no value
    is shown. Example: “%.4f” shows float to four decimals.


    * **value_fontsize** (*float*) – Font size for values. Default is 10.


    * **symbol_fontsize** (*float*) – Font size for element symbols. Default is 14.


    * **cmap** (*str*) – Color scheme of the heatmap. Default is ‘YlOrRd’.
    Refer to the matplotlib documentation for other options.


    * **blank_color** (*str*) – Color assigned for the missing elements in
    elemental_data. Default is “grey”.


    * **edge_color** (*str*) – Color assigned for the edge of elements in the
    periodic table. Default is “white”.


    * **max_row** (*int*) – Maximum number of rows of the periodic table to be
    shown. Default is 9, which means the periodic table heat map covers
    the standard 7 rows of the periodic table + 2 rows for the lanthanides
    and actinides. Use a value of max_row = 7 to exclude the lanthanides and
    actinides.


    * **readable_fontcolor** (*bool*) – Whether to use readable font color depending
    on background color. Default is False.


    * **pymatviz** (*bool*) – Whether to use pymatviz to generate the heatmap. Defaults to True.
    See [https://github.com/janosh/pymatviz](https://github.com/janosh/pymatviz).


    * **kwargs** – Passed to pymatviz.ptable_heatmap_plotly



### pymatgen.util.plotting.pretty_plot(width=8, height=None, plt=None, dpi=None, color_cycle=('qualitative', 'Set1_9'))
Provides a publication quality plot, with nice defaults for font sizes etc.


* **Parameters**


    * **width** (*float*) – Width of plot in inches. Defaults to 8in.


    * **height** (*float*) – Height of plot in inches. Defaults to width \* golden
    ratio.


    * **plt** (*matplotlib.pyplot*) – If plt is supplied, changes will be made to an
    existing plot. Otherwise, a new plot will be created.


    * **dpi** (*int*) – Sets dot per inch for figure. Defaults to 300.


    * **color_cycle** (*tuple*) – Set the color cycle for new plots to one of the
    color sets in palettable. Defaults to a qualitative Set1_9.



* **Returns**

    Matplotlib plot object with properly sized fonts.



### pymatgen.util.plotting.pretty_plot_two_axis(x, y1, y2, xlabel=None, y1label=None, y2label=None, width=8, height=None, dpi=300, \*\*plot_kwargs)
Variant of pretty_plot that does a dual axis plot. Adapted from matplotlib
examples. Makes it easier to create plots with different axes.


* **Parameters**


    * **x** (*np.ndarray/list*) – Data for x-axis.


    * **y1** (*dict/np.ndarray/list*) – Data for y1 axis (left). If a dict, it will
    be interpreted as a {label: sequence}.


    * **y2** (*dict/np.ndarray/list*) – Data for y2 axis (right). If a dict, it will
    be interpreted as a {label: sequence}.


    * **xlabel** (*str*) – If not None, this will be the label for the x-axis.


    * **y1label** (*str*) – If not None, this will be the label for the y1-axis.


    * **y2label** (*str*) – If not None, this will be the label for the y2-axis.


    * **width** (*float*) – Width of plot in inches. Defaults to 8in.


    * **height** (*float*) – Height of plot in inches. Defaults to width \* golden
    ratio.


    * **dpi** (*int*) – Sets dot per inch for figure. Defaults to 300.


    * **plot_kwargs** – Passthrough kwargs to matplotlib’s plot method. E.g.,
    linewidth, etc.



* **Returns**

    matplotlib.pyplot



### pymatgen.util.plotting.pretty_polyfit_plot(x, y, deg=1, xlabel=None, ylabel=None, \*\*kwargs)
Convenience method to plot data with trend lines based on polynomial fit.


* **Parameters**


    * **x** – Sequence of x data.


    * **y** – Sequence of y data.


    * **deg** (*int*) – Degree of polynomial. Defaults to 1.


    * **xlabel** (*str*) – Label for x-axis.


    * **ylabel** (*str*) – Label for y-axis.


    * **kwargs** – Keyword args passed to pretty_plot.



* **Returns**

    matplotlib.pyplot object.



### pymatgen.util.plotting.van_arkel_triangle(list_of_materials, annotate=True)
A static method that generates a binary van Arkel-Ketelaar triangle to

    quantify the ionic, metallic and covalent character of a compound
    by plotting the electronegativity difference (y) vs average (x).
    See:

    > A.E. van Arkel, Molecules and Crystals in Inorganic Chemistry,

    >     Interscience, New York (1956)

    and

        J.A.A Ketelaar, Chemical Constitution (2nd edition), An Introduction

            to the Theory of the Chemical Bond, Elsevier, New York (1958).


* **Parameters**


    * **list_of_materials** (*list*) – A list of computed entries of binary
    materials or a list of lists containing two elements (str).


    * **annotate** (*bool*) – Whether or not to label the points on the
    triangle with reduced formula (if list of entries) or pair
    of elements (if list of list of str).


## pymatgen.util.provenance module

Classes and methods related to the Structure Notation Language (SNL).


### _class_ pymatgen.util.provenance.Author(name, email)
Bases: `Author`

An Author contains two fields: name and email. It is meant to represent
the author of a Structure or the author of a code that was applied to a Structure.

Create new instance of Author(name, email)


#### as_dict()
Returns: MSONable dict.


#### _static_ from_dict(d)

* **Parameters**

    **d** (*dict*) – Dict representation.



* **Returns**

    Author



#### _static_ parse_author(author)
Parses an Author object from either a String, dict, or tuple.


* **Parameters**

    **author** – A String formatted as “NAME <[email@domain.com](mailto:email@domain.com)>”,
    (name, email) tuple, or a dict with name and email keys.



* **Returns**

    An Author object.



### _class_ pymatgen.util.provenance.HistoryNode(name, url, description)
Bases: `HistoryNode`

A HistoryNode represents a step in the chain of events that lead to a
Structure. HistoryNodes leave ‘breadcrumbs’ so that you can trace back how
a Structure was created. For example, a HistoryNode might represent pulling
a Structure from an external database such as the ICSD or CSD. Or, it might
represent the application of a code (e.g. pymatgen) to the Structure, with
a custom description of how that code was applied (e.g. a site removal
Transformation was applied).

A HistoryNode contains three fields:


#### name()
The name of a code or resource that this Structure encountered in
its history (String)


#### url()
The URL of that code/resource (String)


#### description()
A free-form description of how the code/resource is related to the
Structure (dict).

Create new instance of HistoryNode(name, url, description)


#### as_dict()
Returns: Dict.


#### _static_ from_dict(dct: dict[str, str])

* **Parameters**

    **dct** (*dict*) – Dict representation.



* **Returns**

    HistoryNode



#### _static_ parse_history_node(h_node)
Parses a History Node object from either a dict or a tuple.


* **Parameters**

    **h_node** – A dict with name/url/description fields or a 3-element
    tuple.



* **Returns**

    History node.



### _class_ pymatgen.util.provenance.StructureNL(struct_or_mol, authors, projects=None, references='', remarks=None, data=None, history=None, created_at=None)
Bases: `object`

The Structure Notation Language (SNL, pronounced ‘snail’) is a container for a pymatgen
Structure/Molecule object with some additional fields for enhanced provenance.

It is meant to be imported/exported in a JSON file format with the following structure:


    * sites


    * lattice (optional)


    * about


            * created_at


            * authors


            * projects


            * references


            * remarks


            * data


            * history


* **Parameters**


    * **struct_or_mol** – A pymatgen.core.structure Structure/Molecule object


    * **authors** – *List* of {“name”:’’, “email”:’’} dicts,
    *list* of Strings as ‘John Doe <[johndoe@gmail.com](mailto:johndoe@gmail.com)>’,
    or a single String with commas separating authors


    * **projects** – List of Strings [‘Project A’, ‘Project B’]


    * **references** – A String in BibTeX format


    * **remarks** – List of Strings [‘Remark A’, ‘Remark B’]


    * **data** – A free form dict. Namespaced at the root level with an
    underscore, e.g. {“_materialsproject”: <custom data>}


    * **history** – List of dicts - [{‘name’:’’, ‘url’:’’, ‘description’:{}}]


    * **created_at** – A datetime object.



#### as_dict()
Returns: MSONable dict.


#### _classmethod_ from_dict(d)

* **Parameters**

    **d** (*dict*) – Dict representation.



* **Returns**

    Class



#### _classmethod_ from_structures(structures: Sequence[[Structure](pymatgen.core.md#pymatgen.core.structure.Structure)], authors: Sequence[dict[str, str]], projects=None, references='', remarks=None, data=None, histories=None, created_at=None)
A convenience method for getting a list of StructureNL objects by
specifying structures and metadata separately. Some of the metadata
is applied to all of the structures for ease of use.


* **Parameters**


    * **structures** – A list of Structure objects


    * **authors** – *List* of {“name”:’’, “email”:’’} dicts,
    *list* of Strings as ‘John Doe <[johndoe@gmail.com](mailto:johndoe@gmail.com)>’,
    or a single String with commas separating authors


    * **projects** – List of Strings [‘Project A’, ‘Project B’]. This
    applies to all structures.


    * **references** – A String in BibTeX format. Again, this applies to all
    structures.


    * **remarks** – List of Strings [‘Remark A’, ‘Remark B’]


    * **data** – A list of free form dict. Namespaced at the root level
    with an underscore, e.g. {“_materialsproject”:<custom data>}
    . The length of data should be the same as the list of
    structures if not None.


    * **histories** – List of list of dicts - [[{‘name’:’’, ‘url’:’’,
    ‘description’:{}}], …] The length of histories should be the
    same as the list of structures if not None.


    * **created_at** – A datetime object



### pymatgen.util.provenance.is_valid_bibtex(reference: str)
Use pybtex to validate that a reference is in proper BibTeX format.


* **Parameters**

    **reference** – A String reference in BibTeX format.



* **Returns**

    Boolean indicating if reference is valid bibtex.


## pymatgen.util.string module

This module provides utility classes for string operations.


### _class_ pymatgen.util.string.Stringify()
Bases: `object`

Mix-in class for string formatting, e.g. superscripting numbers and symbols or superscripting.


#### STRING_MODE(_ = 'SUBSCRIPT_ )

#### to_html_string()
Generates a HTML formatted string. This uses the output from to_latex_string to generate a HTML output.
:return: HTML formatted string.


#### to_latex_string()
Generates a LaTeX formatted string. The mode is set by the class variable STRING_MODE, which defaults to
“SUBSCRIPT”. E.g., Fe2O3 is transformed to Fe$_{2}$O$_{3}$. Setting STRING_MODE to “SUPERSCRIPT” creates
superscript, e.g., Fe2+ becomes Fe^{2+}. The initial string is obtained from the class’s __str__ method.


* **Returns**

    String for display as in LaTeX with proper superscripts and subscripts.



#### to_pretty_string()

* **Returns**

    A pretty string representation. By default, the __str__ output is used, but this method can be
    overridden if a different representation from default is desired.



#### to_unicode_string()

* **Returns**

    Unicode string with proper sub and superscripts. Note that this works only with systems where the sub
    and superscripts are pure integers.



### pymatgen.util.string.charge_string(charge, brackets=True, explicit_one=True)
Returns a string representing the charge of an Ion. By default, the
charge is placed in brackets with the sign preceding the magnitude, e.g.,
‘[+2]’. For uncharged species, the string returned is ‘(aq)’.


* **Parameters**


    * **charge** – the charge of the Ion


    * **brackets** – whether to enclose the charge in brackets, e.g. [+2]. Default: True


    * **explicit_one** – whether to include the number one for monovalent ions, e.g.
    +1 rather than +. Default: True



### pymatgen.util.string.disordered_formula(disordered_struct, symbols=('x', 'y', 'z'), fmt='plain')
Returns a formula of a form like AxB1-x (x=0.5)
for disordered structures. Will only return a
formula for disordered structures with one
kind of disordered site at present.


* **Parameters**


    * **disordered_struct** – a disordered structure


    * **symbols** – a tuple of characters to use for


    * **subscripts** (*'x'**, **'y'**, **'z'*) –


    * **is** (*by default this*) –


    * **disordered** (*but if you have more than three*) –


    * **added** (*species more symbols will need to be*) –


    * **fmt** (*str*) – ‘plain’, ‘HTML’ or ‘LaTeX’


Returns (str): a disordered formula string


### pymatgen.util.string.formula_double_format(afloat, ignore_ones=True, tol: float = 1e-08)
This function is used to make pretty formulas by formatting the amounts.
Instead of Li1.0 Fe1.0 P1.0 O4.0, you get LiFePO4.


* **Parameters**


    * **afloat** (*float*) – a float


    * **ignore_ones** (*bool*) – if true, floats of 1 are ignored.


    * **tol** (*float*) – Tolerance to round to nearest int. i.e. 2.0000000001 -> 2



* **Returns**

    A string representation of the float for formulas.



### pymatgen.util.string.htmlify(formula)
Generates a HTML formatted formula, e.g. Fe2O3 is transformed to
Fe<sub>2</sub>O</sub>3</sub>.

Note that Composition now has a to_html_string() method that may
be used instead.


* **Parameters**

    **formula** –



* **Returns**




### pymatgen.util.string.latexify(formula)
Generates a LaTeX formatted formula. E.g., Fe2O3 is transformed to
Fe$_{2}$O$_{3}$.

Note that Composition now has a to_latex_string() method that may
be used instead.


* **Parameters**

    **formula** (*str*) – Input formula.



* **Returns**

    Formula suitable for display as in LaTeX with proper subscripts.



### pymatgen.util.string.latexify_spacegroup(spacegroup_symbol)
Generates a latex formatted spacegroup. E.g., P2_1/c is converted to
P2$_{1}$/c and P-1 is converted to P$\\overline{1}$.

Note that SymmetryGroup now has a to_latex_string() method that may
be called instead.


* **Parameters**

    **spacegroup_symbol** (*str*) – A spacegroup symbol



* **Returns**

    A latex formatted spacegroup with proper subscripts and overlines.



### pymatgen.util.string.str_delimited(results, header=None, delimiter='\\t')
Given a tuple of tuples, generate a delimited string form.
>>> results = [[“a”,”b”,”c”],[“d”,”e”,”f”],[1,2,3]]
>>> print(str_delimited(results,delimiter=”,”))
a,b,c
d,e,f
1,2,3.


* **Parameters**


    * **results** – 2d sequence of arbitrary types.


    * **header** – optional header


    * **delimiter** – Defaults to “t” for tab-delimited output.



* **Returns**

    Aligned string output in a table-like format.



### pymatgen.util.string.stream_has_colours(stream)
True if stream supports colours. Python cookbook, #475186.


### pymatgen.util.string.transformation_to_string(matrix, translation_vec=(0, 0, 0), components=('x', 'y', 'z'), c='', delim=',')
Convenience method. Given matrix returns string, e.g. x+2y+1/4
:param matrix
:param translation_vec
:param components: either (‘x’, ‘y’, ‘z’) or (‘a’, ‘b’, ‘c’)
:param c: optional additional character to print (used for magmoms)
:param delim: delimiter
:return: xyz string.


### pymatgen.util.string.unicodeify(formula)
Generates a formula with unicode subscripts, e.g. Fe2O3 is transformed
to Fe₂O₃. Does not support formulae with decimal points.

Note that Composition now has a to_unicode_string() method that may
be used instead.


* **Parameters**

    **formula** –



* **Returns**




### pymatgen.util.string.unicodeify_spacegroup(spacegroup_symbol)
Generates a unicode formatted spacegroup. E.g., P2$_{1}$/c is converted to
P2₁/c and P$\\overline{1}$ is converted to P̅1.

Note that SymmetryGroup now has a to_unicode_string() method that
may be called instead.


* **Parameters**

    **spacegroup_symbol** (*str*) – A spacegroup symbol as LaTeX



* **Returns**

    A unicode spacegroup with proper subscripts and overlines.



### pymatgen.util.string.unicodeify_species(specie_string)
Generates a unicode formatted species string, with appropriate
superscripts for oxidation states.

Note that Species now has a to_unicode_string() method that
may be used instead.


* **Parameters**

    **specie_string** (*str*) – Species string, e.g. O2-



* **Returns**

    Species string, e.g. O²⁻


## pymatgen.util.testing module

Common test support for pymatgen test scripts.

This single module should provide all the common functionality for pymatgen
tests in a single location, so that test scripts can just import it and work
right away.


### _class_ pymatgen.util.testing.PymatgenTest(methodName='runTest')
Bases: `TestCase`

Extends unittest.TestCase with several assert methods for array and str comparison.

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### MODULE_DIR(_ = PosixPath('/Users/shyue/repos/pymatgen/pymatgen/util'_ )

#### STRUCTURES_DIR(_ = PosixPath('/Users/shyue/repos/pymatgen/pymatgen/util/structures'_ )

#### TEST_FILES_DIR(_ = PosixPath('/Users/shyue/repos/pymatgen/pymatgen/util/../../test_files'_ )

#### TEST_STRUCTURES(_: ClassVar[dict[str, [pymatgen.core.structure.Structure](pymatgen.core.md#pymatgen.core.structure.Structure)]_ _ = {'BaNiO3': <pymatgen.alchemy.materials.TransformedStructure object>, 'CsCl': Structure Summary Lattice     abc : 4.209 4.209 4.209  angles : 90.0 90.0 90.0  volume : 74.56530132899998       A : 4.209 0.0 0.0       B : 0.0 4.209 0.0       C : 0.0 0.0 4.209     pbc : True True True PeriodicSite: Cs (0.0, 0.0, 0.0) [0.0, 0.0, 0.0] PeriodicSite: Cl (2.104, 2.104, 2.104) [0.5, 0.5, 0.5], 'Graphite': Structure Summary Lattice     abc : 2.4700000000000006 2.47 6.8  angles : 90.0 90.0 120.00000000000001  volume : 35.928033824449685       A : -1.2350000000000005 -2.139082747347564 -3.024877593893963e-16       B : -1.2349999999999997 2.139082747347564 1.5124387969469814e-16       C : 0.0 0.0 -6.8     pbc : True True True PeriodicSite: C0+ (-1.235, -0.7132, -5.1) [0.6667, 0.3333, 0.75] PeriodicSite: C0+ (-1.235, 0.7132, -1.7) [0.3333, 0.6667, 0.25] PeriodicSite: C0+ (0.0, 0.0, -5.1) [0.0, 0.0, 0.75] PeriodicSite: C0+ (0.0, 0.0, -1.7) [0.0, 0.0, 0.25], 'He_BCC': Structure Summary Lattice     abc : 2.737172073807164 2.7371720737557492 2.73717207  angles : 109.47122060669534 109.47122060631484 109.47122066304705  volume : 15.786447515629305       A : 2.58063058 0.0 -0.91239069       B : -1.29031529 2.23489164 -0.91239069       C : 0.0 0.0 2.73717207     pbc : True True True PeriodicSite: He (0.0, 0.0, 0.0) [0.0, 0.0, 0.0], 'K2O2': Structure Summary Lattice     abc : 4.862375062662279 4.862375062662278 6.361  angles : 90.0 90.0 91.89149308502677  volume : 150.30921504900002       A : 3.3810000000000002 -3.4945000000000004 0.0       B : 3.381 3.4944999999999995 4.2795282396204257e-16       C : 0.0 0.0 -6.361     pbc : True True True PeriodicSite: K+ (1.69, -1.086, -1.59) [0.4054, 0.0946, 0.25] PeriodicSite: K+ (1.69, 1.086, -4.771) [0.0946, 0.4054, 0.75] PeriodicSite: K+ (5.071, 1.086, -4.771) [0.5946, 0.9054, 0.75] PeriodicSite: K+ (5.072, -1.086, -1.59) [0.9054, 0.5946, 0.25] PeriodicSite: O- (3.381, -0.6262, -3.642) [0.5896, 0.4104, 0.5725] PeriodicSite: O- (3.381, 0.6262, -2.719) [0.4104, 0.5896, 0.4275] PeriodicSite: O- (3.381, -2.868, -0.4612) [0.9104, 0.0896, 0.0725] PeriodicSite: O- (3.381, 2.868, -5.9) [0.0896, 0.9104, 0.9275], 'La2CoO4F': Structure Summary Lattice     abc : 6.847249 6.847249 5.748369  angles : 90.0 90.0 132.754623  volume : 197.89340779864318       A : 6.847249 0.0 4.192730785407457e-16       B : -4.648323410186749 5.027714027499061 4.192730785407457e-16       C : 0.0 0.0 5.748369     pbc : True True True PeriodicSite: La (La3+) (-3.061, 4.334, 2.771) [0.138, 0.862, 0.482] PeriodicSite: La (La3+) (5.26, 0.6941, 2.978) [0.862, 0.138, 0.518] PeriodicSite: La (La3+) (-0.4875, 3.208, 5.645) [0.362, 0.638, 0.982] PeriodicSite: La (La3+) (2.686, 1.82, 0.1033) [0.638, 0.362, 0.01798] PeriodicSite: Co (Co3+) (1.099, 2.514, 2.874) [0.5, 0.5, 0.5] PeriodicSite: Co (Co3+) (0.0, 0.0, 0.0) [0.0, 0.0, 0.0] PeriodicSite: O (O2-) (-0.9933, 3.429, 3.297) [0.318, 0.6821, 0.5735] PeriodicSite: O (O2-) (3.192, 1.599, 2.452) [0.6821, 0.318, 0.4265] PeriodicSite: O (O2-) (-2.556, 4.112, 0.4224) [0.182, 0.8179, 0.07348] PeriodicSite: O (O2-) (4.754, 0.9153, 5.326) [0.8179, 0.182, 0.9265] PeriodicSite: O (O2-) (1.426, 3.869, 1.437) [0.7305, 0.7695, 0.25] PeriodicSite: O (O2-) (1.873, 3.673, 4.311) [0.7695, 0.7305, 0.75] PeriodicSite: O (O2-) (0.3261, 1.355, 1.437) [0.2305, 0.2695, 0.25] PeriodicSite: O (O2-) (0.7734, 1.159, 4.311) [0.2695, 0.2305, 0.75] PeriodicSite: F (F-) (-2.324, 2.514, 4.311) [0.0, 0.5, 0.75] PeriodicSite: F (F2-) (3.424, 0.0, 1.437) [0.5, 0.0, 0.25], 'Li10GeP2S12': Structure Summary Lattice     abc : 8.69407 8.69407 12.5994  angles : 90.0 90.0 90.0  volume : 952.3489977658411       A : -8.69407 0.0 -5.323582498531514e-16       B : 5.323582498531514e-16 -8.69407 -5.323582498531514e-16       C : 0.0 0.0 -12.5994     pbc : True True True PeriodicSite: Li:0.691 (-6.466, -6.331, -10.29) [0.7437, 0.7282, 0.8168] PeriodicSite: Li:0.691 (-2.228, -2.363, -10.29) [0.2563, 0.2718, 0.8168] PeriodicSite: Li:0.691 (-6.71, -2.119, -3.991) [0.7718, 0.2437, 0.3168] PeriodicSite: Li:0.691 (-1.984, -6.575, -3.991) [0.2282, 0.7563, 0.3168] PeriodicSite: Li:0.691 (-6.575, -1.984, -8.608) [0.7563, 0.2282, 0.6832] PeriodicSite: Li:0.691 (-2.119, -6.71, -8.608) [0.2437, 0.7718, 0.6832] PeriodicSite: Li:0.691 (-6.331, -6.466, -2.308) [0.7282, 0.7437, 0.1832] PeriodicSite: Li:0.691 (-2.363, -2.228, -2.308) [0.2718, 0.2563, 0.1832] PeriodicSite: Li:0.691 (-6.575, -6.71, -8.608) [0.7563, 0.7718, 0.6832] PeriodicSite: Li:0.691 (-2.119, -1.984, -8.608) [0.2437, 0.2282, 0.6832] PeriodicSite: Li:0.691 (-6.331, -2.228, -2.308) [0.7282, 0.2563, 0.1832] PeriodicSite: Li:0.691 (-2.363, -6.466, -2.308) [0.2718, 0.7437, 0.1832] PeriodicSite: Li:0.691 (-6.466, -2.363, -10.29) [0.7437, 0.2718, 0.8168] PeriodicSite: Li:0.691 (-2.228, -6.331, -10.29) [0.2563, 0.7282, 0.8168] PeriodicSite: Li:0.691 (-6.71, -6.575, -3.991) [0.7718, 0.7563, 0.3168] PeriodicSite: Li:0.691 (-1.984, -2.119, -3.991) [0.2282, 0.2437, 0.3168] PeriodicSite: Li (2.662e-16, -4.347, -0.698) [0.0, 0.5, 0.0554] PeriodicSite: Li (2.662e-16, -4.347, -6.998) [0.0, 0.5, 0.5554] PeriodicSite: Li (-4.347, 0.0, -5.602) [0.5, 0.0, 0.4446] PeriodicSite: Li (-4.347, 0.0, -11.9) [0.5, 0.0, 0.9446] PeriodicSite: Li:0.643 (-6.553, -6.553, -8.025e-16) [0.7537, 0.7537, 0.0] PeriodicSite: Li:0.643 (-2.141, -2.141, -2.622e-16) [0.2463, 0.2463, 0.0] PeriodicSite: Li:0.643 (-6.488, -2.206, -6.3) [0.7463, 0.2537, 0.5] PeriodicSite: Li:0.643 (-2.206, -6.488, -6.3) [0.2537, 0.7463, 0.5] PeriodicSite: Li:0.643 (-6.488, -6.488, -6.3) [0.7463, 0.7463, 0.5] PeriodicSite: Li:0.643 (-2.206, -2.206, -6.3) [0.2537, 0.2537, 0.5] PeriodicSite: Li:0.643 (-6.553, -2.141, -5.324e-16) [0.7537, 0.2463, 0.0] PeriodicSite: Li:0.643 (-2.141, -6.553, -5.324e-16) [0.2463, 0.7537, 0.0] PeriodicSite: Ge:0.515, P:0.485 (2.662e-16, -4.347, -3.897) [0.0, 0.5, 0.3093] PeriodicSite: Ge:0.515, P:0.485 (2.662e-16, -4.347, -10.2) [0.0, 0.5, 0.8093] PeriodicSite: Ge:0.515, P:0.485 (-4.347, 0.0, -2.403) [0.5, 0.0, 0.1907] PeriodicSite: Ge:0.515, P:0.485 (-4.347, 0.0, -8.702) [0.5, 0.0, 0.6907] PeriodicSite: P (0.0, 0.0, -6.3) [0.0, 0.0, 0.5] PeriodicSite: P (-4.347, -4.347, -5.324e-16) [0.5, 0.5, 0.0] PeriodicSite: S (-8.694, -7.092, -7.43) [1.0, 0.8157, 0.5897] PeriodicSite: S (9.811e-17, -1.602, -7.43) [0.0, 0.1843, 0.5897] PeriodicSite: S (-5.949, -4.347, -1.13) [0.6843, 0.5, 0.0897] PeriodicSite: S (-2.745, -4.347, -1.13) [0.3157, 0.5, 0.0897] PeriodicSite: S (-4.347, -2.745, -11.47) [0.5, 0.3157, 0.9103] PeriodicSite: S (-4.347, -5.949, -11.47) [0.5, 0.6843, 0.9103] PeriodicSite: S (-7.092, 0.0, -5.17) [0.8157, 0.0, 0.4103] PeriodicSite: S (-1.602, 0.0, -5.17) [0.1843, 0.0, 0.4103] PeriodicSite: S (-8.694, -6.094, -11.4) [1.0, 0.7009, 0.905] PeriodicSite: S (1.592e-16, -2.6, -11.4) [0.0, 0.2991, 0.905] PeriodicSite: S (-6.947, -4.347, -5.103) [0.7991, 0.5, 0.405] PeriodicSite: S (-1.747, -4.347, -5.103) [0.2009, 0.5, 0.405] PeriodicSite: S (-4.347, -1.747, -7.497) [0.5, 0.2009, 0.595] PeriodicSite: S (-4.347, -6.947, -7.497) [0.5, 0.7991, 0.595] PeriodicSite: S (-6.094, 0.0, -1.197) [0.7009, 0.0, 0.095] PeriodicSite: S (-2.6, 0.0, -1.197) [0.2991, 0.0, 0.095] PeriodicSite: S (1.602e-16, -2.617, -2.628) [0.0, 0.301, 0.2086] PeriodicSite: S (-8.694, -6.077, -2.628) [1.0, 0.699, 0.2086] PeriodicSite: S (-1.73, -4.347, -8.928) [0.199, 0.5, 0.7086] PeriodicSite: S (-6.964, -4.347, -8.928) [0.801, 0.5, 0.7086] PeriodicSite: S (-4.347, -6.964, -3.671) [0.5, 0.801, 0.2914] PeriodicSite: S (-4.347, -1.73, -3.671) [0.5, 0.199, 0.2914] PeriodicSite: S (-2.617, 0.0, -9.971) [0.301, 0.0, 0.7914] PeriodicSite: S (-6.077, 0.0, -9.971) [0.699, 0.0, 0.7914], 'Li2O': Structure Summary Lattice     abc : 3.291071792359756 3.291071899625086 3.2910720568557887  angles : 60.12971043288485 60.12970952137675 60.12970313039097  volume : 25.279668381289053       A : 2.91738857 0.09789437 1.52000466       B : 0.96463406 2.75503561 1.52000466       C : 0.13320635 0.09789443 3.28691771     pbc : True True True PeriodicSite: O2- (0.0, 0.0, 0.0) [0.0, 0.0, 0.0] PeriodicSite: Li+ (3.012, 2.214, 4.746) [0.7502, 0.7502, 0.7502] PeriodicSite: Li+ (1.003, 0.7372, 1.581) [0.2498, 0.2498, 0.2498], 'Li2O2': Structure Summary Lattice     abc : 3.1830000000000007 3.1830000000000003 7.7258  angles : 90.0 90.0 120.00000000000001  volume : 67.7871492344378       A : -1.5915000000000006 -2.756558860245869 -3.898050761686025e-16       B : -1.5914999999999992 2.756558860245869 1.9490253808430124e-16       C : 0.0 0.0 -7.7258     pbc : True True True PeriodicSite: Li+ (0.0, 0.0, -3.863) [0.0, 0.0, 0.5] PeriodicSite: Li+ (0.0, 0.0, 0.0) [0.0, 0.0, 0.0] PeriodicSite: Li+ (-1.592, -0.919, -5.794) [0.6667, 0.3333, 0.75] PeriodicSite: Li+ (-1.591, 0.919, -1.931) [0.3333, 0.6667, 0.25] PeriodicSite: O- (-1.592, -0.919, -1.157) [0.6667, 0.3333, 0.1497] PeriodicSite: O- (-1.591, 0.919, -6.569) [0.3333, 0.6667, 0.8503] PeriodicSite: O- (-1.591, 0.919, -5.019) [0.3333, 0.6667, 0.6497] PeriodicSite: O- (-1.592, -0.919, -2.706) [0.6667, 0.3333, 0.3503], 'Li3V2(PO4)3': Structure Summary Lattice     abc : 8.64068 8.67265 8.64068  angles : 61.21317 61.3442 61.21317  volume : 470.7978186322687       A : -2.470939478415996 -7.168424451387926 -4.143609518420984       B : -7.600861301025969 0.0 -4.17633397910965       C : 0.0 0.0 -8.64068     pbc : True True True PeriodicSite: Li0+ (-9.822, -7.12, -8.177) [0.9933, 0.9693, 0.00152] PeriodicSite: Li0+ (-5.265, -3.573, -8.66) [0.4985, 0.5307, 0.5067] PeriodicSite: Li0+ (-6.086, -1.117, -6.752) [0.1558, 0.75, 0.3442] PeriodicSite: Li0+ (-3.999, -6.087, -10.19) [0.8492, 0.25, 0.6508] PeriodicSite: Li0+ (-4.474, -5.342, -5.905) [0.7452, 0.3464, 0.1586] PeriodicSite: Li0+ (-2.011, -2.447, -8.578) [0.3414, 0.1536, 0.7548] PeriodicSite: V0+ (-8.617, -6.127, -14.56) [0.8548, 0.8558, 0.8621] PeriodicSite: V0+ (-6.473, -4.573, -10.91) [0.6379, 0.6442, 0.6452] PeriodicSite: V0+ (-3.574, -2.541, -5.995) [0.3544, 0.355, 0.3522] PeriodicSite: V0+ (-1.467, -1.059, -2.475) [0.1478, 0.145, 0.1456] PeriodicSite: P0+ (-7.85, -1.779, -9.759) [0.2482, 0.9521, 0.5503] PeriodicSite: P0+ (-6.511, -6.808, -8.4) [0.9497, 0.5479, 0.2518] PeriodicSite: P0+ (-6.831, -3.279, -5.396) [0.4574, 0.75, 0.04264] PeriodicSite: P0+ (-3.254, -3.927, -11.54) [0.5478, 0.25, 0.9522] PeriodicSite: P0+ (-3.56, -0.2774, -8.539) [0.0387, 0.4558, 0.7494] PeriodicSite: P0+ (-2.191, -5.381, -7.281) [0.7506, 0.04424, 0.4613] PeriodicSite: O0+ (-8.398, -5.114, -10.91) [0.7133, 0.873, 0.4991] PeriodicSite: O0+ (-6.51, -3.511, -12.57) [0.4899, 0.6973, 0.8827] PeriodicSite: O0+ (-7.426, -0.5548, -10.63) [0.07739, 0.9518, 0.7325] PeriodicSite: O0+ (-7.945, -1.4, -8.259) [0.1953, 0.9818, 0.3876] PeriodicSite: O0+ (-6.793, -2.912, -9.816) [0.4062, 0.7617, 0.5732] PeriodicSite: O0+ (-5.754, -6.299, -11.69) [0.8787, 0.4713, 0.7039] PeriodicSite: O0+ (-7.902, -6.644, -7.735) [0.9268, 0.7384, 0.09382] PeriodicSite: O0+ (-6.063, -5.502, -9.121) [0.7675, 0.5482, 0.4226] PeriodicSite: O0+ (-7.627, -4.425, -5.998) [0.6173, 0.8027, 0.01015] PeriodicSite: O0+ (-7.57, -1.893, -5.467) [0.2641, 0.9101, 0.06612] PeriodicSite: O0+ (-4.485, -4.075, -10.61) [0.5685, 0.4053, 0.7588] PeriodicSite: O0+ (-4.768, -0.006738, -9.42) [0.00094, 0.627, 0.7867] PeriodicSite: O0+ (-5.39, -7.063, -7.378) [0.9852, 0.3888, 0.1934] PeriodicSite: O0+ (-5.556, -3.11, -6.3) [0.4339, 0.5899, 0.2359] PeriodicSite: O0+ (-2.551, -5.313, -11.52) [0.7412, 0.09467, 0.9316] PeriodicSite: O0+ (-2.396, -2.797, -10.93) [0.3902, 0.1884, 0.9869] PeriodicSite: O0+ (-3.831, -1.659, -7.841) [0.2314, 0.4288, 0.5892] PeriodicSite: O0+ (-2.202, -0.4825, -9.291) [0.06731, 0.2678, 0.9135] PeriodicSite: O0+ (-3.214, -4.204, -7.139) [0.5865, 0.2322, 0.4327] PeriodicSite: O0+ (-4.216, -0.8054, -5.262) [0.1124, 0.5182, 0.3047] PeriodicSite: O0+ (-2.185, -5.707, -8.787) [0.7961, 0.02868, 0.6213] PeriodicSite: O0+ (-2.792, -6.529, -6.392) [0.9108, 0.07118, 0.2686] PeriodicSite: O0+ (-3.636, -3.678, -4.376) [0.5131, 0.3116, 0.1098] PeriodicSite: O0+ (-1.603, -2.198, -6.183) [0.3066, 0.1112, 0.5148], 'LiFePO4': Structure Summary Lattice     abc : 4.7448 6.0657700000000006 10.41037  angles : 90.50178999999999 90.00019 90.00362  volume : 299.60796771125047       A : 0.0 0.0 -4.7448       B : -0.053122654367700306 -6.065537365280947 0.00038324092231544317       C : 10.41036999994276 0.0 3.452209424235223e-05     pbc : True True True PeriodicSite: Li (10.36, -6.065, 0.0003703) [1e-05, 1.0, 1.0] PeriodicSite: Li (-0.02646, -3.033, -4.745) [1.0, 0.5, 1e-05] PeriodicSite: Li (5.152, -6.065, -2.372) [0.5, 1.0, 0.5] PeriodicSite: Li (5.179, -3.033, -2.372) [0.5, 0.5, 0.5] PeriodicSite: Fe (2.265, -1.537, -2.491) [0.5251, 0.2534, 0.2188] PeriodicSite: Fe (2.887, -4.528, -0.1187) [0.02508, 0.7465, 0.2812] PeriodicSite: Fe (7.47, -1.537, -4.626) [0.975, 0.2535, 0.7188] PeriodicSite: Fe (8.093, -4.528, -2.253) [0.4749, 0.7465, 0.7812] PeriodicSite: P (0.9432, -4.559, -2.761) [0.5821, 0.7517, 0.09444] PeriodicSite: P (4.209, -1.506, -0.3893) [0.08207, 0.2483, 0.4056] PeriodicSite: P (6.148, -4.56, -4.355) [0.9179, 0.7517, 0.5944] PeriodicSite: P (9.414, -1.506, -1.983) [0.4179, 0.2483, 0.9056] PeriodicSite: O (0.4361, -1.523, -1.383) [0.2916, 0.2511, 0.04317] PeriodicSite: O (0.9618, -4.552, -1.226) [0.2585, 0.7504, 0.09622] PeriodicSite: O (1.675, -5.798, -3.385) [0.7135, 0.9559, 0.1658] PeriodicSite: O (1.695, -3.328, -3.398) [0.7163, 0.5486, 0.1656] PeriodicSite: O (3.457, -2.738, -1.026) [0.2163, 0.4514, 0.3344] PeriodicSite: O (3.477, -0.2672, -1.013) [0.2135, 0.04406, 0.3342] PeriodicSite: O (4.19, -1.514, -3.599) [0.7585, 0.2495, 0.4038] PeriodicSite: O (4.716, -4.542, -3.756) [0.7916, 0.7489, 0.4568] PeriodicSite: O (5.641, -1.523, -0.9889) [0.2084, 0.2511, 0.5432] PeriodicSite: O (6.167, -4.552, -1.146) [0.2415, 0.7505, 0.5962] PeriodicSite: O (6.88, -5.799, -3.732) [0.7866, 0.956, 0.6658] PeriodicSite: O (6.9, -3.328, -3.718) [0.7837, 0.5487, 0.6656] PeriodicSite: O (8.662, -2.738, -1.346) [0.2837, 0.4513, 0.8344] PeriodicSite: O (8.682, -0.2668, -1.36) [0.2866, 0.04399, 0.8342] PeriodicSite: O (9.395, -1.514, -3.518) [0.7415, 0.2496, 0.9038] PeriodicSite: O (9.921, -4.542, -3.361) [0.7084, 0.7489, 0.9568], 'NaFePO4': Structure Summary Lattice     abc : 4.9955 6.28746 10.440590000000002  angles : 90.0 89.97269 90.0  volume : 327.9285211911844       A : 0.0 0.0 -4.9955       B : 3.849958881883509e-16 -6.28746 -3.849958881883509e-16       C : -10.440588813976833 0.0 -0.004976500966151329     pbc : True True True PeriodicSite: Na (-10.44, -6.287, -5.0) [1.0, 0.9999, 1.0] PeriodicSite: Na (-10.44, -3.144, -5.0) [1.0, 0.5001, 1.0] PeriodicSite: Na (-5.22, -6.287, -2.501) [0.5001, 0.9999, 0.4999] PeriodicSite: Na (-5.22, -3.144, -2.501) [0.5001, 0.5001, 0.4999] PeriodicSite: Fe (-8.206, -1.572, -2.619) [0.5235, 0.25, 0.7859] PeriodicSite: Fe (-7.458, -4.716, -0.1238) [0.02408, 0.75, 0.7144] PeriodicSite: Fe (-2.984, -1.572, -4.879) [0.9764, 0.25, 0.2859] PeriodicSite: Fe (-2.237, -4.716, -2.38) [0.4763, 0.75, 0.2142] PeriodicSite: P (-9.337, -4.716, -2.838) [0.5672, 0.75, 0.8943] PeriodicSite: P (-6.325, -1.572, -0.337) [0.06686, 0.25, 0.6058] PeriodicSite: P (-4.116, -4.716, -4.662) [0.9329, 0.75, 0.3943] PeriodicSite: P (-1.104, -1.572, -2.164) [0.4331, 0.25, 0.1057] PeriodicSite: O (-10.08, -1.572, -1.664) [0.3321, 0.25, 0.965] PeriodicSite: O (-9.28, -4.716, -1.299) [0.2591, 0.75, 0.8888] PeriodicSite: O (-8.613, -5.953, -3.488) [0.6974, 0.9468, 0.825] PeriodicSite: O (-8.613, -3.478, -3.488) [0.6974, 0.5532, 0.825] PeriodicSite: O (-7.047, -2.81, -0.9879) [0.1971, 0.447, 0.6749] PeriodicSite: O (-7.047, -0.3333, -0.9879) [0.1971, 0.05301, 0.6749] PeriodicSite: O (-6.383, -1.572, -3.794) [0.7588, 0.25, 0.6114] PeriodicSite: O (-5.584, -4.716, -4.156) [0.8314, 0.75, 0.5348] PeriodicSite: O (-4.855, -1.572, -0.842) [0.1681, 0.25, 0.465] PeriodicSite: O (-4.059, -4.716, -1.206) [0.241, 0.75, 0.3887] PeriodicSite: O (-3.393, -5.954, -4.011) [0.8026, 0.9469, 0.325] PeriodicSite: O (-3.393, -3.478, -4.011) [0.8026, 0.5531, 0.325] PeriodicSite: O (-1.826, -2.81, -1.513) [0.3027, 0.4469, 0.1749] PeriodicSite: O (-1.826, -0.3338, -1.513) [0.3027, 0.05309, 0.1749] PeriodicSite: O (-1.162, -1.572, -3.703) [0.7412, 0.25, 0.1113] PeriodicSite: O (-0.365, -4.716, -3.339) [0.6683, 0.75, 0.03496], 'Pb2TiZrO6': Structure Summary Lattice     abc : 4.505383 5.67298 5.67298  angles : 90.0 90.0 90.0  volume : 144.99539884709878       A : 0.0 0.0 4.505383       B : 5.67298 0.0 0.0       C : 0.0 5.67298 0.0     pbc : True True True PeriodicSite: Zr (2.836, 0.0, 2.438) [0.541, 0.5, -0.0] PeriodicSite: Ti (0.0, 2.836, 2.419) [0.5369, -0.0, 0.5] PeriodicSite: Pb (0.0, 0.0, 4.411) [0.9791, -0.0, -0.0] PeriodicSite: Pb (2.836, 2.836, 4.411) [0.9791, 0.5, 0.5] PeriodicSite: O (1.37, 1.466, 2.875) [0.6381, 0.2416, 0.2584] PeriodicSite: O (1.37, 4.207, 2.875) [0.6381, 0.2416, 0.7416] PeriodicSite: O (4.303, 1.466, 2.875) [0.6381, 0.7584, 0.2584] PeriodicSite: O (4.303, 4.207, 2.875) [0.6381, 0.7584, 0.7416] PeriodicSite: O (0.0, 2.836, 0.6698) [0.1487, -0.0, 0.5] PeriodicSite: O (2.836, 0.0, 0.4587) [0.1018, 0.5, -0.0], 'Si': Structure Summary Lattice     abc : 3.8401979337 3.840198994344244 3.8401979337177736  angles : 119.99999086398421 90.0 60.0000091373222  volume : 40.04479464425159       A : 3.8401979337 0.0 0.0       B : 1.9200989668 3.3257101909 0.0       C : 0.0 -2.2171384943 3.1355090603     pbc : True True True PeriodicSite: Si (0.0, 0.0, 0.0) [0.0, 0.0, 0.0] PeriodicSite: Si (3.84, 1.225e-06, 2.352) [0.75, 0.5, 0.75], 'SiO2': Structure Summary Lattice     abc : 5.0277818 5.0277817873408415 5.51891759  angles : 90.0 90.0 120.00000005747829  volume : 120.81959693044213       A : 5.0277818 6e-08 -0.0       B : -2.51389095 4.35418672 0.0       C : -0.0 0.0 5.51891759     pbc : True True True PeriodicSite: Si (1.314, 2.276, 5.519) [0.5227, 0.5227, 1.0] PeriodicSite: Si (2.4, -1.242e-05, 3.679) [0.4773, -2.86e-06, 0.6667] PeriodicSite: Si (3.828, 2.078, 1.84) [1.0, 0.4773, 0.3333] PeriodicSite: O (0.8313, 3.654, 4.805) [0.585, 0.8393, 0.8707] PeriodicSite: O (-1.066, 3.247, 2.965) [0.1607, 0.7457, 0.5373] PeriodicSite: O (0.2352, 1.807, 1.126) [0.2543, 0.415, 0.204] PeriodicSite: O (2.749, 2.547, 0.7137) [0.8393, 0.585, 0.1293] PeriodicSite: O (3.345, 0.6999, 2.553) [0.7457, 0.1607, 0.4626] PeriodicSite: O (1.447, 1.107, 4.393) [0.415, 0.2543, 0.796], 'Si_SiO2_Interface': Structure Summary Lattice     abc : 10.160210588671637 10.160210588671637 29.156770002399192  angles : 90.0 90.0 59.99999999999999  volume : 2606.6064276838156       A : 10.160210588671637 0.0 6.221334688039882e-16       B : 5.080105294335819 8.799000477589283 6.221334688039882e-16       C : 0.0 0.0 29.156770002399192     pbc : True True True PeriodicSite: Si (10.34, 3.876, 9.35) [0.7976, 0.4405, 0.3207] PeriodicSite: Si (13.97, 5.133, 9.35) [1.083, 0.5833, 0.3207] PeriodicSite: Si (7.439, 6.39, 9.35) [0.369, 0.7262, 0.3207] PeriodicSite: Si (11.07, 7.647, 9.35) [0.6548, 0.869, 0.3207] PeriodicSite: Si (9.616, 0.1048, 9.35) [0.9405, 0.0119, 0.3207] PeriodicSite: Si (3.084, 1.362, 9.35) [0.2262, 0.1548, 0.3207] PeriodicSite: Si (6.713, 2.619, 9.35) [0.5119, 0.2976, 0.3207] PeriodicSite: Si (10.34, 3.876, 11.7) [0.7976, 0.4405, 0.4013] PeriodicSite: Si (13.97, 5.133, 11.7) [1.083, 0.5833, 0.4013] PeriodicSite: Si (7.439, 6.39, 11.7) [0.369, 0.7262, 0.4013] PeriodicSite: Si (11.07, 7.647, 11.7) [0.6548, 0.869, 0.4013] PeriodicSite: Si (9.616, 0.1048, 11.7) [0.9405, 0.0119, 0.4013] PeriodicSite: Si (3.084, 1.362, 11.7) [0.2262, 0.1548, 0.4013] PeriodicSite: Si (6.713, 2.619, 11.7) [0.5119, 0.2976, 0.4013] PeriodicSite: Si (3.868, 2.1, 13.7) [0.2614, 0.2386, 0.4699] PeriodicSite: Si (6.408, 6.499, 13.7) [0.2614, 0.7386, 0.4699] PeriodicSite: Si (8.948, 2.1, 13.7) [0.7614, 0.2386, 0.4699] PeriodicSite: Si (11.49, 6.499, 13.7) [0.7614, 0.7386, 0.4699] PeriodicSite: Si (4.965, 4.4, 15.54) [0.2386, 0.5, 0.533] PeriodicSite: Si (2.425, 0.0, 15.54) [0.2386, 0.0, 0.533] PeriodicSite: Si (10.04, 4.4, 15.54) [0.7386, 0.5, 0.533] PeriodicSite: Si (7.505, 0.0, 15.54) [0.7386, 0.0, 0.533] PeriodicSite: Si (1.328, 2.3, 17.38) [0.0, 0.2614, 0.5961] PeriodicSite: Si (3.868, 6.699, 17.38) [0.0, 0.7614, 0.5961] PeriodicSite: Si (6.408, 2.3, 17.38) [0.5, 0.2614, 0.5961] PeriodicSite: Si (8.948, 6.699, 17.38) [0.5, 0.7614, 0.5961] PeriodicSite: O (3.38, 0.7072, 14.42) [0.2925, 0.08037, 0.4944] PeriodicSite: O (5.92, 5.107, 14.42) [0.2925, 0.5804, 0.4944] PeriodicSite: O (8.46, 0.7072, 14.42) [0.7925, 0.08037, 0.4944] PeriodicSite: O (11.0, 5.107, 14.42) [0.7925, 0.5804, 0.4944] PeriodicSite: O (1.463, 1.119, 16.25) [0.08037, 0.1271, 0.5575] PeriodicSite: O (4.003, 5.518, 16.25) [0.08037, 0.6271, 0.5575] PeriodicSite: O (6.543, 1.119, 16.25) [0.5804, 0.1271, 0.5575] PeriodicSite: O (9.083, 5.518, 16.25) [0.5804, 0.6271, 0.5575] PeriodicSite: O (2.778, 2.574, 18.09) [0.1271, 0.2925, 0.6206] PeriodicSite: O (5.318, 6.973, 18.09) [0.1271, 0.7925, 0.6206] PeriodicSite: O (7.858, 2.574, 18.09) [0.6271, 0.2925, 0.6206] PeriodicSite: O (10.4, 6.973, 18.09) [0.6271, 0.7925, 0.6206] PeriodicSite: O (5.318, 1.826, 18.51) [0.4196, 0.2075, 0.6347] PeriodicSite: O (7.858, 6.225, 18.51) [0.4196, 0.7075, 0.6347] PeriodicSite: O (10.4, 1.826, 18.51) [0.9196, 0.2075, 0.6347] PeriodicSite: O (12.94, 6.225, 18.51) [0.9196, 0.7075, 0.6347] PeriodicSite: O (5.92, 3.692, 16.67) [0.3729, 0.4196, 0.5716] PeriodicSite: O (8.46, 8.092, 16.67) [0.3729, 0.9196, 0.5716] PeriodicSite: O (11.0, 3.692, 16.67) [0.8729, 0.4196, 0.5716] PeriodicSite: O (13.54, 8.092, 16.67) [0.8729, 0.9196, 0.5716] PeriodicSite: O (4.003, 3.281, 14.83) [0.2075, 0.3729, 0.5085] PeriodicSite: O (6.543, 7.68, 14.83) [0.2075, 0.8729, 0.5085] PeriodicSite: O (9.083, 3.281, 14.83) [0.7075, 0.3729, 0.5085] PeriodicSite: O (11.62, 7.68, 14.83) [0.7075, 0.8729, 0.5085], 'Sn': Structure Summary Lattice     abc : 6.65061477 6.65061477 6.65061477  angles : 90.0 90.0 90.0  volume : 294.16119253915326       A : 6.65061477 0.0 0.0       B : 0.0 6.65061477 0.0       C : 0.0 0.0 6.65061477     pbc : True True True PeriodicSite: Sn (2.494, 5.819, 2.494) [0.375, 0.875, 0.375] PeriodicSite: Sn (0.8313, 0.8313, 0.8313) [0.125, 0.125, 0.125] PeriodicSite: Sn (2.494, 2.494, 5.819) [0.375, 0.375, 0.875] PeriodicSite: Sn (0.8313, 4.157, 4.157) [0.125, 0.625, 0.625] PeriodicSite: Sn (5.819, 5.819, 5.819) [0.875, 0.875, 0.875] PeriodicSite: Sn (4.157, 0.8313, 4.157) [0.625, 0.125, 0.625] PeriodicSite: Sn (5.819, 2.494, 2.494) [0.875, 0.375, 0.375] PeriodicSite: Sn (4.157, 4.157, 0.8313) [0.625, 0.625, 0.125], 'SrTiO3': Structure Summary Lattice     abc : 3.905 3.905 3.905  angles : 90.0 90.0 90.0  volume : 59.54744262499999       A : 3.905 0.0 2.391122875335207e-16       B : -2.391122875335207e-16 3.905 2.391122875335207e-16       C : 0.0 0.0 3.905     pbc : True True True PeriodicSite: Sr2+ (1.952, 1.952, 1.953) [0.5, 0.5, 0.5] PeriodicSite: Ti4+ (0.0, 0.0, 0.0) [0.0, 0.0, 0.0] PeriodicSite: O2- (0.0, 0.0, 1.952) [0.0, 0.0, 0.5] PeriodicSite: O2- (-1.196e-16, 1.952, 1.196e-16) [0.0, 0.5, 0.0] PeriodicSite: O2- (1.952, 0.0, 1.196e-16) [0.5, 0.0, 0.0], 'TiO2': <pymatgen.alchemy.materials.TransformedStructure object>, 'TlBiSe2': Structure Summary Lattice     abc : 4.372201544153475 61.149528191254596 59.08304095958131  angles : 3.6737691014605547 58.80413768014368 59.97014152897428  volume : 825.8885392060803       A : 3.74140218 -0.00027862 2.26231209       B : -1.13595131 3.78041753 61.02198666       C : 0.036913 0.04093202 59.08301525     pbc : True True True PeriodicSite: Tl (2.515, 3.674, 63.43) [0.9669, 0.9717, 0.03303] PeriodicSite: Tl (3.604, 0.2035, 6.715) [0.9793, 0.05368, 0.02071] PeriodicSite: Tl (3.573, 0.5326, 10.99) [0.9977, 0.1409, 0.002243] PeriodicSite: Tl (-0.2434, 0.8101, 13.08) [9.01e-06, 0.2143, -6.6e-07] PeriodicSite: Tl (-0.2812, 1.128, 76.51) [0.002322, 0.2876, 0.9978] PeriodicSite: Tl (-0.3122, 1.457, 80.78) [0.02071, 0.3749, 0.9793] PeriodicSite: Tl (-0.3594, 1.767, 85.09) [0.03312, 0.4569, 0.967] PeriodicSite: Bi (-0.4251, 2.055, 89.42) [0.03883, 0.5333, 0.9611] PeriodicSite: Bi (-0.5892, 2.228, 93.91) [0.008427, 0.5786, 0.9916] PeriodicSite: Bi (-0.675, 2.493, 98.28) [0.006704, 0.6486, 0.9933] PeriodicSite: Bi (-0.8113, 2.7, 43.59) [3.182e-05, 0.7143, -3.841e-05] PeriodicSite: Bi (2.831, 2.949, 50.24) [0.9934, 0.7801, 0.006551] PeriodicSite: Bi (2.744, 3.213, 54.61) [0.9914, 0.8498, 0.0086] PeriodicSite: Bi (2.58, 3.386, 59.1) [0.9611, 0.8952, 0.03891] PeriodicSite: Se (1.738, 0.2085, 35.15) [0.4743, 0.0495, 0.5256] PeriodicSite: Se (1.699, 0.5271, 39.44) [0.4897, 0.1339, 0.5104] PeriodicSite: Se (1.646, 0.8311, 43.75) [0.5002, 0.2145, 0.4998] PeriodicSite: Se (1.594, 1.135, 48.06) [0.5106, 0.2949, 0.4893] PeriodicSite: Se (1.553, 1.451, 52.35) [0.5255, 0.3788, 0.4747] PeriodicSite: Se (1.487, 1.741, 56.68) [0.5312, 0.4556, 0.4685] PeriodicSite: Se (1.42, 2.028, 61.02) [0.5364, 0.5314, 0.4634] PeriodicSite: Se (1.244, 2.185, 65.54) [0.5013, 0.5726, 0.4986] PeriodicSite: Se (1.183, 2.479, 69.86) [0.5087, 0.6503, 0.4913] PeriodicSite: Se (1.078, 2.721, 74.26) [0.5, 0.7143, 0.5] PeriodicSite: Se (0.9728, 2.963, 78.66) [0.4913, 0.7782, 0.5088] PeriodicSite: Se (0.9121, 3.256, 82.98) [0.4987, 0.856, 0.5014] PeriodicSite: Se (0.7353, 3.414, 87.5) [0.4636, 0.8973, 0.5365] PeriodicSite: Se (0.668, 3.7, 91.84) [0.4687, 0.9729, 0.5316], 'VO2': <pymatgen.alchemy.materials.TransformedStructure object>_ )

#### _static_ assert_all_close(actual, desired, decimal=7, err_msg='', verbose=True)
Tests if two arrays are almost equal up to some relative or absolute tolerance.


#### assert_msonable(obj, test_is_subclass=True)
Test if obj is MSONable and verify the contract is fulfilled.

By default, the method tests whether obj is an instance of MSONable.
This check can be deactivated by setting test_is_subclass=False.


#### _static_ assert_str_content_equal(actual, expected)
Tests if two strings are equal, ignoring things like trailing spaces, etc.


#### fn(_ = PosixPath('/Users/shyue/repos/pymatgen/pymatgen/util/structures/SrTiO3.json'_ )

#### _classmethod_ get_structure(name: str)
Get a structure from the template directories.


* **Parameters**

    **name** (*str*) – Name of structure file.



* **Returns**

    Structure



#### serialize_with_pickle(objects, protocols=None, test_eq=True)
Test whether the object(s) can be serialized and deserialized with
pickle. This method tries to serialize the objects with pickle and the
protocols specified in input. Then it deserializes the pickle format
and compares the two objects with the __eq__ operator if
test_eq is True.


* **Parameters**


    * **objects** – Object or list of objects.


    * **protocols** – List of pickle protocols to test. If protocols is None,
    HIGHEST_PROTOCOL is tested.


    * **test_eq** – If True, the deserialized object is compared with the
    original object using the __eq__ method.



* **Returns**

    Nested list with the objects deserialized with the specified
    protocols.


## pymatgen.util.typing module

This module defines convenience types for type hinting purposes.
Type hinting is new to pymatgen, so this module is subject to
change until best practices are established.