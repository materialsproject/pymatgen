Change log
==========

v2017.11.30
-----------
* Fix for severe enumlib_caller bug. This causes enumerations not to be carried
  out properly due to bad accounting of symmetry of ordered sites. It results
  in too few orderings.
* New method to extract clusters of atoms from a Molecule based on bonds.
