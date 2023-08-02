---
layout: default
title: pymatgen.optimization.linear_assignment_numpy.md
nav_exclude: true
---

# pymatgen.optimization.linear_assignment_numpy module

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