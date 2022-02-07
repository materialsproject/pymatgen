Change log
==========

v2022.2.7
---------
* Critical bug fix for pmgrc.yaml being overwritten in MPRester in a non-standard way.
* Change in config file for Lobster basis. Removed the 2p orbitals for Be as they led to problems in our computations and probably should be optional during the projection. (@JaGeo) 
* Return None for ISPIN=1 for `Vasprun('vasprun.xml').complete_dos.spin_polarization`.

