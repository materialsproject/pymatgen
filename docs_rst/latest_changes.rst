Change log
==========

v2018.9.30
----------
* Fix: increased cut-off to VoronoiNN to avoid scipy crash (@utf)
* Fix: Outcar parsing issues with certain values of electrostatic potential (@sivonxay)
* Fix: bug in EnumlibAdaptor/EnumerateStructureTransformation involving incorrect
  stoichiometries in some instances (#1286) (@shyuep)
* Fix: fractional co-ordinate finite precision errors in CifParser, now
  also includes additional warnings for implicit hydrogens (@mkhorton)
* New features and improvements to GBGenerator (@ucsdlxg, @shyuep)
* New analysis options in StructureGraph, speed up tests (@mkhorton)
* New utility function to pretty print disordered formulae, along with a
  ordered-to-disordered structure transformation (@mkhorton)
* Ability to use pymatgen's StructureMatcher against AFLOW's library of
  crystallographic prototypes (@mkhorton)
* Make Chgcar serializable to/from dict for database insertion (@jmmshn)
