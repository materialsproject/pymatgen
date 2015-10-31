Change log
==========

v3.2.4
------
* GaussianOutput can now parse frequencies, normal modes and cartesian forces
  (Xin Chen).
* Support for Aiida<->pymatgen conversion by the Aiida development team (Andrius
  Merkys).
* Specialized BSVasprun parser that is ~2-3x faster than Vasprun.
* Refactor the boltztrap package (merge a few methods together) and add several
  new methods (power factor, seebeck...)
* Support of the new PCM format in QChem 4.3
* Local environment analysis to pmg script.
* Deprecate prettytable in favor of tabulate package.
* Improvements to MITNEBVaspInputSet.
* Misc bug fixes.
