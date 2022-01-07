Change log
==========

v2022.1.7
---------
* ASE io improvements (e.g., magnetic moments and selective dynamics transfer). @arosen93
* New automatic k-point generation scheme, `automatic_density_by_lengths`, which allows the user to specify a density of k-points in each dimension (rather than just for the entire volume). @arosen93 
* Build improvements to dynamically generate C code by running Cython on pyx files rather than having hard-generated .c files.
