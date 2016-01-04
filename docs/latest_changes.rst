Change log
==========

v3.3.0
------
* Updated and checked for Python 3.5.* compatibility.
* Element, Spin, Orbital and various other Enum-like classes are now actually
  implemented using Enum (with enum34 dependency for Python < 3.4).
* Speed up Site creation by 20% for ordered sites, with cost in terms of
  slightly slower non-ordered Sites. Since ordered Sites is the far more common
  case, this gives significant boost for large scale manipulations of
  structures.
* Alternative, more pythonic syntax for creating supercells via simply
  Structure * 3 or Structure * (3, 1, 1).
* zeo++ fixes.
* More stable incar settings for MITMDVaspInputSet.
