---
layout: default
title: pymatgen.optimization.neighbors.md
nav_exclude: true
---

# pymatgen.optimization.neighbors module


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