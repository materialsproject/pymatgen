Change log
==========

v2023.7.14
----------
- Emergency bug fix release to remove use of sys.path in pymatgen.io.ase package.
- Fix "Incompatible POTCAR" error on ComputedEntries with oxidation states.
- New global config variable `PMG_POTCAR_CHECKS` provides means to disable all POTCAR checking.
