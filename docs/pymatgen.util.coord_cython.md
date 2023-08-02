---
layout: default
title: pymatgen.util.coord_cython.md
nav_exclude: true
---

# pymatgen.util.coord_cython module

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