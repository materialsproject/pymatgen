Change log
==========

v2022.4.26
----------
* Fix dipole units in recent vasp versions (at least 6.3, maybe even before) (@@fraricci)
* Removed complex numbers from the definition of WSWQ (@jmmshn)
* MP database version logging is now no longer logged in the .pmgrc.yaml but rather in the .mprester.log.yaml.
  This avoids the MPRester constantly rewriting a config file and causing users' pymatgen to completely fail.
