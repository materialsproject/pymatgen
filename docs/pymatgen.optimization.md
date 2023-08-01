---
layout: default
title: pymatgen.optimization.md
nav_exclude: true
---

# pymatgen.optimization package

Optimization utilities.

## Subpackages


* [pymatgen.optimization.tests package](pymatgen.optimization.tests.md)




    * [pymatgen.optimization.tests.test_linear_assignment module](pymatgen.optimization.tests.md#module-pymatgen.optimization.tests.test_linear_assignment)


        * [`LinearAssignmentTest`](pymatgen.optimization.tests.md#pymatgen.optimization.tests.test_linear_assignment.LinearAssignmentTest)


            * [`LinearAssignmentTest.another_test_case()`](pymatgen.optimization.tests.md#pymatgen.optimization.tests.test_linear_assignment.LinearAssignmentTest.another_test_case)


            * [`LinearAssignmentTest.test()`](pymatgen.optimization.tests.md#pymatgen.optimization.tests.test_linear_assignment.LinearAssignmentTest.test)


            * [`LinearAssignmentTest.test_boolean_inputs()`](pymatgen.optimization.tests.md#pymatgen.optimization.tests.test_linear_assignment.LinearAssignmentTest.test_boolean_inputs)


            * [`LinearAssignmentTest.test_rectangular()`](pymatgen.optimization.tests.md#pymatgen.optimization.tests.test_linear_assignment.LinearAssignmentTest.test_rectangular)


            * [`LinearAssignmentTest.test_small_range()`](pymatgen.optimization.tests.md#pymatgen.optimization.tests.test_linear_assignment.LinearAssignmentTest.test_small_range)


    * [pymatgen.optimization.tests.test_neighbors module](pymatgen.optimization.tests.md#module-pymatgen.optimization.tests.test_neighbors)


        * [`NeighborsTestCase`](pymatgen.optimization.tests.md#pymatgen.optimization.tests.test_neighbors.NeighborsTestCase)


            * [`NeighborsTestCase.setUp()`](pymatgen.optimization.tests.md#pymatgen.optimization.tests.test_neighbors.NeighborsTestCase.setUp)


            * [`NeighborsTestCase.test_points_in_spheres()`](pymatgen.optimization.tests.md#pymatgen.optimization.tests.test_neighbors.NeighborsTestCase.test_points_in_spheres)



## pymatgen.optimization.linear_assignment module

This module contains an algorithm to solve the Linear Assignment Problem


### _class_ pymatgen.optimization.linear_assignment.LinearAssignment(costs, epsilon=1e-13)
Bases: `object`

This class finds the solution to the Linear Assignment Problem.
It finds a minimum cost matching between two sets, given a cost
matrix.

This class is an implementation of the LAPJV algorithm described in:
R. Jonker, A. Volgenant. A Shortest Augmenting Path Algorithm for
Dense and Sparse Linear Assignment Problems. Computing 38, 325-340
(1987)


* **Parameters**


    * **costs** – The cost matrix of the problem. cost[i,j] should be the
    cost of matching x[i] to y[j]. The cost matrix may be
    rectangular


    * **epsilon** – Tolerance for determining if solution vector is < 0


<!-- attribute: min_cost:

The minimum cost of the matching -->
<!-- attribute: solution:

The matching of the rows to columns. i.e solution = [1, 2, 0]
would match row 0 to column 1, row 1 to column 2 and row 2
to column 0. Total cost would be c[0, 1] + c[1, 2] + c[2, 0] -->
## pymatgen.optimization.linear_assignment_numpy module

This module contains an algorithm to solve the Linear Assignment Problem.
It has the same functionality as linear_assignment.pyx, but is much slower
as it is vectorized in numpy rather than cython.


### _class_ pymatgen.optimization.linear_assignment_numpy.LinearAssignment(costs, epsilon=1e-06)
Bases: `object`

This class finds the solution to the Linear Assignment Problem.
It finds a minimum cost matching between two sets, given a cost
matrix.

This class is an implementation of the LAPJV algorithm described in:
R. Jonker, A. Volgenant. A Shortest Augmenting Path Algorithm for
Dense and Sparse Linear Assignment Problems. Computing 38, 325-340
(1987)

<!-- attribute: min_cost:

The minimum cost of the matching -->
<!-- attribute: solution:

The matching of the rows to columns. i.e solution = [1, 2, 0]
would match row 0 to column 1, row 1 to column 2 and row 2
to column 0. Total cost would be c[0, 1] + c[1, 2] + c[2, 0] -->

* **Parameters**


    * **costs** – The cost matrix of the problem. cost[i,j] should be the
    cost of matching x[i] to y[j]. The cost matrix may be
    rectangular


    * **epsilon** – Tolerance for determining if solution vector is < 0.



#### _property_ min_cost()
Returns the cost of the best assignment.

## pymatgen.optimization.neighbors module


### pymatgen.optimization.neighbors.find_points_in_spheres(all_coords, center_coords, r, pbc, lattice, tol=1e-08, min_r=1.0)
For each point in center_coords, get all the neighboring points in all_coords
that are within the cutoff radius r. All the coordinates should be Cartesian.


* **Parameters**


    * **all_coords** – (np.ndarray[double, dim=2]) all available points.
    When periodic boundary is considered, this is all the points in the lattice.


    * **center_coords** – (np.ndarray[double, dim=2]) all centering points


    * **r** – (float) cutoff radius


    * **pbc** – (np.ndarray[long, dim=1]) whether to set periodic boundaries


    * **lattice** – (np.ndarray[double, dim=2]) 3x3 lattice matrix


    * **tol** – (float) numerical tolerance


    * **min_r** – (float) minimal cutoff to calculate the neighbor list
    directly. If the cutoff is less than this value, the algorithm
    will calculate neighbor list using min_r as cutoff and discard
    those that have larger distances.



* **Returns**

    index1 (n, ), index2 (n, ), offset_vectors (n, 3), distances (n, ).
    index1 of center_coords, and index2 of all_coords that form the neighbor pair
    offset_vectors are the periodic image offsets for the all_coords.