Change log
==========

v2017.06.08
-----------
* Switch to date-based version for pymatgen.
* Electronegativities now available for all elements except for He, Ne and
  Ar, which are set to infinity with a warning.
* Bond lengths are now set to sum of atomic radii with warning if not available.
* Bug fixes to boltztrap, symmetry for trigonal-hex systems, etc.
