Change log
==========

v2018.9.19
----------
* Fix to composition handling in `MolecularOrbitals` (@dyllamt)
* Fix to allow mixed compressed/uncompressed loading of VASP band structures (@ajjackson)
* New features and fixes to `chemenv` analysis module (@davidwaroquiers)
* Fix to include structure predictor data with pip/conda-installed pymatgen (@shyamd)
* Fixes to `Defect` objects, icluding allowing rotational supercell transformations (@dbroberg)
* Fix to `BSDOSPlotter` to correctly fill in parts of DOS (@fraricci)
* Added '@' notation parsing in `Composition` (@tamuhey)
* BibTex reference extraction updated in `CifParser` to support ICSD CIFs (@shyamd)
* Various updates to speed up and fix test suite (@shyuep, @fraricci)
* Improvements to BoltzTraP 2 support (@shyuep, @fraricci)
