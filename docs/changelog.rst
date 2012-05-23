Older versions
==============

Version 1.9.0
-------------

1. Completely new json encoder and decoder that support serialization of almost
   all pymatgen objects.
2. Simplification to Borg API utilizing the new json API.
3. Bandstructure classes now support spin-polarized runs.
4. Beta classes for battery (insertion and conversion) analysis.

Version 1.8.3
-------------

1. spglib_adaptor now supports disordered structures.
2. Update to support new spglib with angle_tolerance.
3. Changes to Borg API to support both file and directory style paths.
4. Speed up for COMPLETE_ORDERING algo for PartialRemoveSpecieTransformation.


Version 1.8.1
-------------

1. Revamped transmuter classes for better readability and long term support.
2. Much improved speed for PartialRemoveSpecieTransformations.
3. Misc bug fixes.

Version 1.8.0
-------------

1. Support for additional properties on Specie (Spin) and Site (magmom, charge).
2. Molecule class to support molecules without periodicity.
3. Beta io class for XYZ and GaussianInput.

Version 1.7.2
-------------

1. Bug fixes for vaspio_set and compatibility classes.

Version 1.7.0
-------------

1. Complete reorganization of modules for electronic structure.
2. Beta of band structure classes.
3. Misc improvements to vaspio classes.
4. Bug fixes.

Version 1.6.0
-------------

1. Beta of pymatgen.borg package implemented for high-throughput data assimilation.
2. Added ComputedEntry classes for handling calculated data.
3. New method of specifying VASP pseudopotential location using a VASP_PSP_DIR 
   environment variable. 
4. Bug fix for pymatgen.symmetry
5. Ewald sum speed up by factor of 2 or more.
