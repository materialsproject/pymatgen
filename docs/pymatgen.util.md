---
layout: default
title: pymatgen.util.md
nav_exclude: true
---

1. TOC
{:toc}

# pymatgen.util package

The util package implements various utilities that are commonly used by various
packages.


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


### _exception_ SplineInputError(msg)
Bases: `Exception`

Error for Spline input.


* **Parameters**

    **msg** (*str*) – Message.



### determine_convergence(xs, ys, name, tol: float = 0.0001, extra='', verbose=False, mode='extra', plots=True)
Test it and at which x_value dy(x)/dx < tol for all x >= x_value, conv is true is such a x_value exists.


### exponential(x, a, b, n)
Exponential function base n to fit convergence data.


### extrapolate_reciprocal(xs, ys, n, noise)
Return the parameters such that a + b / x^n hits the last two data points.


### extrapolate_simple_reciprocal(xs, ys)
Extrapolate simple reciprocal function to fit convergence data.


* **Parameters**


    * **xs** – List of x values.


    * **ys** – List of y values.



* **Returns**

    List of parameters [a, b].



### get_derivatives(xs, ys, fd=False)
Return the derivatives of y(x) at the points x
if scipy is available a spline is generated to calculate the derivatives
if scipy is not available the left and right slopes are calculated, if both exist the average is returned
putting fd to zero always returns the finite difference slopes.


### get_weights(xs, ys, mode=2)

* **Parameters**


    * **xs** (*list*) – List of x values.


    * **ys** (*list*) – List of y values.


    * **mode** (*int*) – Mode for calculating weights.



* **Returns**

    List of weights.



* **Return type**

    list



### id_generator(size: int = 8, chars: str = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789')
Generate a random string of specified size and characters.


* **Parameters**


    * **size** (*int*) – The length of the generated string.


    * **chars** (*str*) – The characters to use for generating the string.



* **Returns**

    The generated random string.



* **Return type**

    str



### measure(function, xs, ys, popt, weights)
Measure the quality of a fit.


### multi_curve_fit(xs, ys, verbose)
Fit multiple functions to the x, y data, return the best fit.


### multi_reciprocal_extra(xs, ys, noise=False)
Calculates for a series of powers ns the parameters for which the last two points are at the curve.
With these parameters measure how well the other data points fit.
return the best fit.


### p0_exponential(xs, ys)
Calculate the initial guess parameters for the exponential function.


* **Parameters**


    * **xs** (*list*) – List of x values.


    * **ys** (*list*) – List of y values.



* **Returns**

    List of initial guess parameters [a, b, n].



* **Return type**

    list



### p0_reciprocal(xs, ys)
Predictor for first guess for reciprocal.


### p0_simple_2reciprocal(xs, ys)
Calculate the initial guess parameters for the simple reciprocal function with a power of 2.


* **Parameters**


    * **xs** (*list*) – List of x values.


    * **ys** (*list*) – List of y values.



* **Returns**

    List of initial guess parameters [a, b].



* **Return type**

    list



### p0_simple_4reciprocal(xs, ys)
Calculate the initial guess parameters for the simple reciprocal function with a power of 4.


* **Parameters**


    * **xs** (*list*) – List of x values.


    * **ys** (*list*) – List of y values.



* **Returns**

    The initial guess parameters [a, b].



* **Return type**

    list



### p0_simple_5reciprocal(xs, ys)
Calculate the initial guess parameters for the simple reciprocal function with a power of 0.5.


* **Parameters**


    * **xs** (*list*) – List of x values.


    * **ys** (*list*) – List of y values.



* **Returns**

    List of parameters [a, b].



* **Return type**

    list



### p0_simple_reciprocal(xs, ys)
Calculate the initial guess parameters for the simple reciprocal function.


* **Parameters**


    * **xs** (*list*) – List of x values.


    * **ys** (*list*) – List of y values.



* **Returns**

    List of initial guess parameters [a, b].



* **Return type**

    list



### p0_single_reciprocal(xs, ys)
Calculate the initial guess parameters for the single reciprocal function.


* **Parameters**


    * **xs** (*list*) – List of x values.


    * **ys** (*list*) – List of y values.



* **Returns**

    List of initial guess parameters [a, b, c].



* **Return type**

    list



### print_and_raise_error(xs, ys, name)
Print error message and raise a RuntimeError.


* **Parameters**


    * **xs** (*list*) – List of x values.


    * **ys** (*list*) – List of y values.


    * **name** (*str*) – Name of the function where the error occurred.



### print_plot_line(function, popt, xs, ys, name, tol: float = 0.05, extra='')
Print the gnuplot command line to plot the x, y data with the fitted function using the popt parameters.


### reciprocal(x, a, b, n)
Reciprocal function to the power n to fit convergence data.


### simple_2reciprocal(x, a, b)
Reciprocal function to fit convergence data.


### simple_4reciprocal(x, a, b)
Reciprocal function to fit convergence data.


### simple_5reciprocal(x, a, b)
Reciprocal function to fit convergence data.


### simple_reciprocal(x, a, b)
Reciprocal function to fit convergence data.


### single_reciprocal(x, a, b, c)
Reciprocal function to fit convergence data.

## pymatgen.util.coord module

Utilities for manipulating coordinates or list of coordinates, under periodic
boundary conditions or otherwise. Many of these are heavily vectorized in
numpy for performance.


### _class_ Simplex(coords)
Bases: `MSONable`

A generalized simplex object. See [http://en.wikipedia.org/wiki/Simplex](http://en.wikipedia.org/wiki/Simplex).


#### space_dim()
Dimension of the space. Usually, this is 1 more than the simplex_dim.


* **Type**

    int



#### simplex_dim()
Dimension of the simplex coordinate space.


* **Type**

    int


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
Computes the intersection points of a line with a simplex.


* **Parameters**


    * **point1** (*Sequence**[**float**]*) – 1st point to determine the line.


    * **point2** (*Sequence**[**float**]*) – 2nd point to determine the line.


    * **tolerance** (*float*) – Tolerance for checking if an intersection is in the simplex. Defaults to 1e-8.



* **Returns**

    points where the line intersects the simplex (0, 1, or 2).



#### point_from_bary_coords(bary_coords)

* **Parameters**

    **(****)** (*bary_coords*) – Barycentric coordinates.



* **Returns**

    Point coordinates



#### _property_ volume()
Volume of the simplex.


### all_distances(coords1, coords2)
Returns the distances between two lists of coordinates.


* **Parameters**


    * **coords1** – First set of Cartesian coordinates.


    * **coords2** – Second set of Cartesian coordinates.



* **Returns**

    2d array of Cartesian distances. E.g the distance between
    coords1[i] and coords2[j] is distances[i,j]



### barycentric_coords(coords, simplex)
Converts a list of coordinates to barycentric coordinates, given a
simplex with d+1 points. Only works for d >= 2.


* **Parameters**


    * **coords** – list of n coords to transform, shape should be (n,d)


    * **simplex** – list of coordinates that form the simplex, shape should be
    (d+1, d)



* **Returns**

    a LIST of barycentric coordinates (even if the original input was 1d)



### coord_list_mapping(subset: ArrayLike, superset: ArrayLike, atol: float = 1e-08)
Gives the index mapping from a subset to a superset.
Subset and superset cannot contain duplicate rows.


* **Parameters**


    * **subset** (*ArrayLike*) – List of coords


    * **superset** (*ArrayLike*) – List of coords


    * **atol** (*float*) – Absolute tolerance. Defaults to 1e-8.



* **Returns**

    list of indices such that superset[indices] = subset



### coord_list_mapping_pbc(subset, superset, atol=1e-08, pbc=(True, True, True))
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



### find_in_coord_list(coord_list, coord, atol=1e-08)
Find the indices of matches of a particular coord in a coord_list.


* **Parameters**


    * **coord_list** – List of coords to test


    * **coord** – Specific coordinates


    * **atol** – Absolute tolerance. Defaults to 1e-8. Accepts both scalar and
    array.



* **Returns**

    Indices of matches, e.g., [0, 1, 2, 3]. Empty list if not found.



### find_in_coord_list_pbc(fcoord_list, fcoord, atol=1e-08, pbc=(True, True, True))
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



### get_angle(v1, v2, units='degrees')
Calculates the angle between two vectors.


* **Parameters**


    * **v1** – Vector 1


    * **v2** – Vector 2


    * **units** – “degrees” or “radians”. Defaults to “degrees”.



* **Returns**

    Angle between them in degrees.



### get_linear_interpolated_value(x_values, y_values, x)
Returns an interpolated value by linear interpolation between two values.
This method is written to avoid dependency on scipy, which causes issues on
threading servers.


* **Parameters**


    * **x_values** – Sequence of x values.


    * **y_values** – Corresponding sequence of y values


    * **x** – Get value at particular x



* **Returns**

    Value at x.



### in_coord_list(coord_list, coord, atol=1e-08)
Tests if a particular coord is within a coord_list.


* **Parameters**


    * **coord_list** – List of coords to test


    * **coord** – Specific coordinates


    * **atol** – Absolute tolerance. Defaults to 1e-8. Accepts both scalar and
    array.



* **Returns**

    True if coord is in the coord list.



* **Return type**

    bool



### in_coord_list_pbc(fcoord_list, fcoord, atol=1e-08, pbc=(True, True, True))
Tests if a particular fractional coord is within a fractional coord_list.


* **Parameters**


    * **fcoord_list** – List of fractional coords to test


    * **fcoord** – A specific fractional coord to test.


    * **atol** – Absolute tolerance. Defaults to 1e-8.


    * **pbc** – a tuple defining the periodic boundary conditions along the three
    axis of the lattice.



* **Returns**

    True if coord is in the coord list.



* **Return type**

    bool



### is_coord_subset(subset: ArrayLike, superset: ArrayLike, atol: float = 1e-08)
Tests if all coords in subset are contained in superset.
Doesn’t use periodic boundary conditions.


* **Parameters**


    * **subset** (*ArrayLike*) – List of coords


    * **superset** (*ArrayLike*) – List of coords


    * **atol** (*float*) – Absolute tolerance for comparing coordinates. Defaults to 1e-8.



* **Returns**

    True if all of subset is in superset.



* **Return type**

    bool



### is_coord_subset_pbc(subset, superset, atol=1e-08, mask=None, pbc=(True, True, True))
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



* **Return type**

    bool



### lattice_points_in_supercell(supercell_matrix)
Returns the list of points on the original lattice contained in the
supercell in fractional coordinates (with the supercell basis).
e.g. [[2,0,0],[0,1,0],[0,0,1]] returns [[0,0,0],[0.5,0,0]].


* **Parameters**

    **supercell_matrix** – 3x3 matrix describing the supercell



* **Returns**

    numpy array of the fractional coordinates



### pbc_diff(fcoords1: ArrayLike, fcoords2: ArrayLike, pbc: tuple[bool, bool, bool] = (True, True, True))
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



### pbc_shortest_vectors(lattice, fcoords1, fcoords2, mask=None, return_d2=False)
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


### coord_list_mapping_pbc(subset, superset, atol=1e-08, pbc=(True, True, True))
Gives the index mapping from a subset to a superset.
Superset cannot contain duplicate matching rows


* **Parameters**


    * **subset** – List of frac_coords


    * **superset** – List of frac_coords


    * **pbc** – a tuple defining the periodic boundary conditions along the three
    axis of the lattice.



* **Returns**

    list of indices such that superset[indices] = subset



### is_coord_subset_pbc(subset, superset, atol, mask, pbc=(True, True, True))
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



### pbc_shortest_vectors(lattice, fcoords1, fcoords2, mask=None, return_d2=False, lll_frac_tol=None)
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


### BibTeX(\*args, \*\*kwargs)
Perform no good and no bad.


### Doi(\*args, \*\*kwargs)
Perform no good and no bad.


### _class_ InactiveDueCreditCollector()
Bases: `object`

Just a stub at the Collector which would not do anything.


#### _donothing(\*args, \*\*kwargs)
Perform no good and no bad.


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


### Text(\*args, \*\*kwargs)
Perform no good and no bad.


### Url(\*args, \*\*kwargs)
Perform no good and no bad.


### _donothing_func(\*args, \*\*kwargs)
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


### _hash_label(label, digest_size)

### _init_node_labels(G, edge_attr, node_attr)

### _neighborhood_aggregate(G, node, node_labels, edge_attr=None)
Compute new labels for given node by aggregating
the labels of each node’s neighbors.


### weisfeiler_lehman_graph_hash(G, edge_attr=None, node_attr=None, iterations=3, digest_size=16)
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


### weisfeiler_lehman_subgraph_hashes(G, edge_attr=None, node_attr=None, iterations=3, digest_size=16)
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


### clean_lines(string_list, remove_empty_lines=True)
Strips whitespace, carriage returns and empty lines from a list of strings.


* **Parameters**


    * **string_list** – List of strings


    * **remove_empty_lines** – Set to True to skip lines which are empty after
    stripping.



* **Returns**

    List of clean strings with no whitespaces.



### micro_pyawk(filename, search, results=None, debug=None, postdebug=None)
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


### make_symmetric_matrix_from_upper_tri(val)
Given a symmetric matrix in upper triangular matrix form as flat array indexes as:
[A_xx,A_yy,A_zz,A_xy,A_xz,A_yz]
This will generate the full matrix:
[[A_xx,A_xy,A_xz],[A_xy,A_yy,A_yz],[A_xz,A_yz,A_zz].


### round_to_sigfigs(num, sig_figs)
Rounds a number rounded to a specific number of significant
figures instead of to a specific precision.

## pymatgen.util.numba module

This module provides a wrapper for numba such that no functionality
is lost if numba is not available. Numba is a just-in-time compiler
that can significantly accelerate the evaluation of certain functions
if installed.


### jit(func)
Replacement for numba.jit when numba is not installed that does nothing.


### njit(func)
Replacement for numba.njit when numba is not installed that does nothing.

## pymatgen.util.plotting module

Utilities for generating nicer plots.


### _decide_fontcolor(rgba: tuple)

### add_fig_kwargs(func)
Decorator that adds keyword arguments for functions returning matplotlib
figures.

The function should return either a matplotlib figure or None to signal
some sort of error/unexpected event.
See doc string below for the list of supported options.


### format_formula(formula)
Converts str of chemical formula into
latex format for labelling purposes.


* **Parameters**

    **formula** (*str*) – Chemical formula



### get_ax3d_fig(ax: plt.Axes = None, \*\*kwargs)
Helper function used in plot functions supporting an optional Axes3D
argument. If ax is None, we build the matplotlib figure and create the
Axes3D else we return the current active figure.


* **Parameters**


    * **ax** (*Axes3D**, **optional*) – Axes3D object. Defaults to None.


    * **kwargs** – keyword arguments are passed to plt.figure if ax is not None.



* **Returns**

    matplotlib Axes3D and corresponding figure objects



* **Return type**

    tuple[Axes3D, Figure]



### get_ax_fig(ax: Axes | None = None, \*\*kwargs)
Helper function used in plot functions supporting an optional Axes argument.
If ax is None, we build the matplotlib figure and create the Axes else
we return the current active figure.


* **Parameters**


    * **ax** (*Axes**, **optional*) – Axes object. Defaults to None.


    * **kwargs** – keyword arguments are passed to plt.figure if ax is not None.



* **Returns**

    matplotlib Axes object and Figure objects



* **Return type**

    tuple[plt.Axes, plt.Figure]



### get_axarray_fig_plt(ax_array, nrows=1, ncols=1, sharex=False, sharey=False, squeeze=True, subplot_kw=None, gridspec_kw=None, \*\*fig_kw)
Helper function used in plot functions that accept an optional array of Axes
as argument. If ax_array is None, we build the matplotlib figure and
create the array of Axes by calling plt.subplots else we return the
current active figure.


* **Returns**

    Array of Axes objects
    figure: matplotlib figure
    plt: matplotlib pyplot module.



* **Return type**

    ax



### periodic_table_heatmap(elemental_data=None, cbar_label='', cbar_label_size=14, show_plot=False, cmap='YlOrRd', cmap_range=None, blank_color='grey', edge_color='white', value_format=None, value_fontsize=10, symbol_fontsize=14, max_row=9, readable_fontcolor=False, pymatviz: bool = True, \*\*kwargs)
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



### pretty_plot(width: float = 8, height: float | None = None, ax: plt.Axes = None, dpi: float | None = None, color_cycle: tuple[str, str] = ('qualitative', 'Set1_9'))
Provides a publication quality plot, with nice defaults for font sizes etc.


* **Parameters**


    * **width** (*float*) – Width of plot in inches. Defaults to 8in.


    * **height** (*float*) – Height of plot in inches. Defaults to width \* golden
    ratio.


    * **ax** (*plt.Axes*) – If ax is supplied, changes will be made to an
    existing plot. Otherwise, a new plot will be created.


    * **dpi** (*int*) – Sets dot per inch for figure. Defaults to 300.


    * **color_cycle** (*tuple*) – Set the color cycle for new plots to one of the
    color sets in palettable. Defaults to a qualitative Set1_9.



* **Returns**

    matplotlib axes object with properly sized fonts.



* **Return type**

    plt.Axes



### pretty_plot_two_axis(x, y1, y2, xlabel=None, y1label=None, y2label=None, width=8, height=None, dpi=300, \*\*plot_kwargs)
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



### pretty_polyfit_plot(x, y, deg=1, xlabel=None, ylabel=None, \*\*kwargs)
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



### van_arkel_triangle(list_of_materials, annotate=True)
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


### _class_ Author(name, email)
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



### _class_ HistoryNode(name, url, description)
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
The name of a code or resource that this Structure encountered in its history.


* **Type**

    str



#### url()
The URL of that code/resource.


* **Type**

    str



#### description()
A free-form description of how the code/resource is related to the Structure.


* **Type**

    dict


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



### _class_ StructureNL(struct_or_mol, authors, projects=None, references='', remarks=None, data=None, history=None, created_at=None)
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


    * **struct_or_mol** – A pymatgen Structure/Molecule object


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



### is_valid_bibtex(reference: str)
Use pybtex to validate that a reference is in proper BibTeX format.


* **Parameters**

    **reference** – A String reference in BibTeX format.



* **Returns**

    Boolean indicating if reference is valid bibtex.


## pymatgen.util.string module

This module provides utility classes for string operations.


### _class_ Stringify()
Bases: `object`

Mix-in class for string formatting, e.g. superscripting numbers and symbols or superscripting.


#### STRING_MODE(_ = 'SUBSCRIPT_ )

#### to_html_string()
Generates a HTML formatted string. This uses the output from to_latex_string to generate a HTML output.


* **Returns**

    HTML formatted string.



#### to_latex_string()
Generates a LaTeX formatted string. The mode is set by the class variable STRING_MODE, which defaults to
“SUBSCRIPT”. E.g., Fe2O3 is transformed to Fe$_{2}$O$_{3}$. Setting STRING_MODE to “SUPERSCRIPT” creates
superscript, e.g., Fe2+ becomes Fe^{2+}. The initial string is obtained from the class’s __str__ method.


* **Returns**

    String for display as in LaTeX with proper superscripts and subscripts.



#### to_pretty_string()
A pretty string representation. By default, the __str__ output is used, but this method can be
overridden if a different representation from default is desired.


#### to_unicode_string()
Unicode string with proper sub and superscripts. Note that this works only
with systems where the sub and superscripts are pure integers.


### charge_string(charge, brackets=True, explicit_one=True)
Returns a string representing the charge of an Ion. By default, the
charge is placed in brackets with the sign preceding the magnitude, e.g.,
‘[+2]’. For uncharged species, the string returned is ‘(aq)’.


* **Parameters**


    * **charge** – the charge of the Ion


    * **brackets** – whether to enclose the charge in brackets, e.g. [+2]. Default: True


    * **explicit_one** – whether to include the number one for monovalent ions, e.g.
    +1 rather than +. Default: True



### disordered_formula(disordered_struct, symbols=('x', 'y', 'z'), fmt='plain')
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


### formula_double_format(afloat, ignore_ones=True, tol: float = 1e-08)
This function is used to make pretty formulas by formatting the amounts.
Instead of Li1.0 Fe1.0 P1.0 O4.0, you get LiFePO4.


* **Parameters**


    * **afloat** (*float*) – a float


    * **ignore_ones** (*bool*) – if true, floats of 1 are ignored.


    * **tol** (*float*) – Tolerance to round to nearest int. i.e. 2.0000000001 -> 2



* **Returns**

    A string representation of the float for formulas.



### htmlify(formula)
Generates a HTML formatted formula, e.g. Fe2O3 is transformed to
Fe<sub>2</sub>O</sub>3</sub>.

Note that Composition now has a to_html_string() method that may
be used instead.


* **Parameters**

    **formula** –



### latexify(formula)
Generates a LaTeX formatted formula. E.g., Fe2O3 is transformed to
Fe$_{2}$O$_{3}$.

Note that Composition now has a to_latex_string() method that may
be used instead.


* **Parameters**

    **formula** (*str*) – Input formula.



* **Returns**

    Formula suitable for display as in LaTeX with proper subscripts.



### latexify_spacegroup(spacegroup_symbol)
Generates a latex formatted spacegroup. E.g., P2_1/c is converted to
P2$_{1}$/c and P-1 is converted to P$\\overline{1}$.

Note that SymmetryGroup now has a to_latex_string() method that may
be called instead.


* **Parameters**

    **spacegroup_symbol** (*str*) – A spacegroup symbol



* **Returns**

    A latex formatted spacegroup with proper subscripts and overlines.



### str_delimited(results, header=None, delimiter='\\t')
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



### stream_has_colors(stream)
True if stream supports colors. Python cookbook, #475186.


### transformation_to_string(matrix, translation_vec=(0, 0, 0), components=('x', 'y', 'z'), c='', delim=',')
Convenience method. Given matrix returns string, e.g. x+2y+1/4.


* **Parameters**


    * **matrix** – A 3x3 matrix.


    * **translation_vec** – A 3-element tuple representing the translation vector. Defaults to (0, 0, 0).


    * **components** – A tuple of 3 strings representing the components. Either (‘x’, ‘y’, ‘z’) or (‘a’, ‘b’, ‘c’).
    Defaults to (‘x’, ‘y’, ‘z’).


    * **c** – An optional additional character to print (used for magmoms). Defaults to “”.


    * **delim** – A delimiter. Defaults to “,”.



* **Returns**

    xyz string.



### unicodeify(formula)
Generates a formula with unicode subscripts, e.g. Fe2O3 is transformed
to Fe₂O₃. Does not support formulae with decimal points.

Note that Composition now has a to_unicode_string() method that may
be used instead.


* **Parameters**

    **formula** –



### unicodeify_spacegroup(spacegroup_symbol)
Generates a unicode formatted spacegroup. E.g., P2$_{1}$/c is converted to
P2₁/c and P$\\overline{1}$ is converted to P̅1.

Note that SymmetryGroup now has a to_unicode_string() method that
may be called instead.


* **Parameters**

    **spacegroup_symbol** (*str*) – A spacegroup symbol as LaTeX



* **Returns**

    A unicode spacegroup with proper subscripts and overlines.



### unicodeify_species(specie_string)
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


### _class_ PymatgenTest(methodName='runTest')
Bases: `TestCase`

Extends unittest.TestCase with several assert methods for array and str comparison.

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### TEST_STRUCTURES(_: ClassVar[dict[str | Path, [Structure](pymatgen.core.md#pymatgen.core.structure.Structure) | None]_ _ = {PosixPath('/Users/shyue/repos/pymatgen/pymatgen/util/structures/BaNiO3.json'): None, PosixPath('/Users/shyue/repos/pymatgen/pymatgen/util/structures/CsCl.json'): None, PosixPath('/Users/shyue/repos/pymatgen/pymatgen/util/structures/Graphite.json'): None, PosixPath('/Users/shyue/repos/pymatgen/pymatgen/util/structures/He_BCC.json'): None, PosixPath('/Users/shyue/repos/pymatgen/pymatgen/util/structures/K2O2.json'): None, PosixPath('/Users/shyue/repos/pymatgen/pymatgen/util/structures/La2CoO4F.json'): None, PosixPath('/Users/shyue/repos/pymatgen/pymatgen/util/structures/Li10GeP2S12.json'): None, PosixPath('/Users/shyue/repos/pymatgen/pymatgen/util/structures/Li2O.json'): None, PosixPath('/Users/shyue/repos/pymatgen/pymatgen/util/structures/Li2O2.json'): None, PosixPath('/Users/shyue/repos/pymatgen/pymatgen/util/structures/Li3V2(PO4)3.json'): None, PosixPath('/Users/shyue/repos/pymatgen/pymatgen/util/structures/LiFePO4.json'): None, PosixPath('/Users/shyue/repos/pymatgen/pymatgen/util/structures/NaFePO4.json'): None, PosixPath('/Users/shyue/repos/pymatgen/pymatgen/util/structures/Pb2TiZrO6.json'): None, PosixPath('/Users/shyue/repos/pymatgen/pymatgen/util/structures/Si.json'): None, PosixPath('/Users/shyue/repos/pymatgen/pymatgen/util/structures/SiO2.json'): None, PosixPath('/Users/shyue/repos/pymatgen/pymatgen/util/structures/Si_SiO2_Interface.json'): None, PosixPath('/Users/shyue/repos/pymatgen/pymatgen/util/structures/Sn.json'): None, PosixPath('/Users/shyue/repos/pymatgen/pymatgen/util/structures/SrTiO3.json'): None, PosixPath('/Users/shyue/repos/pymatgen/pymatgen/util/structures/TiO2.json'): None, PosixPath('/Users/shyue/repos/pymatgen/pymatgen/util/structures/TlBiSe2.json'): None, PosixPath('/Users/shyue/repos/pymatgen/pymatgen/util/structures/VO2.json'): None_ )

#### _multiprocess_shared_(_ = Tru_ )

#### _tmp_dir(tmp_path: Path, monkeypatch: MonkeyPatch)

#### assert_msonable(obj, test_is_subclass=True)
Test if obj is MSONable and verify the contract is fulfilled.

By default, the method tests whether obj is an instance of MSONable.
This check can be deactivated by setting test_is_subclass=False.


#### _static_ assert_str_content_equal(actual, expected)
Tests if two strings are equal, ignoring things like trailing spaces, etc.


#### _classmethod_ get_structure(name: str)
Lazily load a structure from pymatgen/util/testing/structures.


* **Parameters**

    **name** (*str*) – Name of structure file.



* **Returns**

    Structure



#### serialize_with_pickle(objects: Any, protocols: Sequence[int] | None = None, test_eq: bool = True)
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