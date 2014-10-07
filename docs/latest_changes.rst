Change log
==========

v3.0.5
------
* Completely revamped symmetry package. The finder.SymmetryFinder and
  pointgroup and spacegroup modules are now deprecated. Instead,
  all symmetry analysis is in the :module:`pymatgen.symmetry.analyzer`_
  module. There is also a completely rewritten support for symmetry groups in
  :module:`pymatgen.symmetry.groups`_. Structure now supports a static
  constructor to generate a structure from a spacegroup (see examples).
* BatteryAnalyzer class (Anubhav Jain) to provide for some common analysis of
  intercalation electrodes.
* Minor bug fixes for structure_matcher, lattice, abinitio.
* MOAB qadapter for abinit. (Liam Damewood)
