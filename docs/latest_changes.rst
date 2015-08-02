Change log
==========

v3.1.3
------
* Major refactoring of pymatgen.io. Now, the io suffix is dropped from all io
  classes. i.e., it is just pymatgen.io.vasp, not pymatgen.io.vaspio. Also, all
  input sets have been moved within the relevant package, e.g.,
  pymatgen.io.vasp.sets. All changes are backwards compatible for now. But
  deprecation messages have been included which states that the stubs will be
  removed in pymatgen 4.0. Pls migrate code when you see the deprecation
  messages.
* Make Composition.anonymized_formula truly chemistry independent (No A2B2
  for peroxides or A2 for diatomic gasses) 
* Allowing CIF data_* header to be prefixed with spaces and tabulations.
