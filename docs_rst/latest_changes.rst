Change log
==========

v2022.0.1
---------
* `pymatgen`, `pymatgen.ext`, `pymatgen.io` and `pymatgen.analysis` are now
  namespace packages. Note that this does not affect normal usage of pymatgen
  from v2022.0.0. All imports remain the same. However, it does allow developers
  to write "add-ons" to these subpackages. A full documentation with examples
  and templates is in the works to guide developers on how to write these
  packages.
