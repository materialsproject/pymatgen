Change log
==========

v2020.3.13
----------
* Added angle_tolerance to CifWriter.
* Change default float precision in CifWriter to 8. Adds float_prec kwarg to 
  allow setting of arbitrary precision. 
* Rudimentary pymatgen.io.vasp.help.VaspDoc class for obtaining help from VASP wiki.
* Massive documentation cleanup.
* Reorganization of Entry, ComputedEntry (@ayushsgupta).
* Bug fix for PourbaixDiagram (@montoyjh).
* Read WAVECAR from gamma-point only VASP executable. (@bernstei)
