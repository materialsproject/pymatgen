Change log
==========

v2021.3.3
---------
* **Backwards incompatible*: pymatgen.SETTINGS have been moved to 
  pymatgen.settings.SETTINGS. In general, this should not lead to many breakages
  since most of these settings are used within pymatgen itself.
* **Backwards incompatible*: pymatgen.loadfn and get_structure_from_mp have been
  removed since no one was using them. 
* critic2_caller has been refactored. (@samblau)
* Improved hash for Compositon (@CompRhys)
* Fixes Outcar parsing for VASP 6.2.0. (@MichaelWolloch)
* Allow None for Gaussian functional, bset, charge and multiplicity (@eimrek)
