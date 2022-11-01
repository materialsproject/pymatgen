Change log
==========

v2022.11.1
----------
* Order of kwargs `fmt` and `filename` in `Structure.to()` swapped for ease of use (note: this can break codes that do not use these options as kwargs).
* @yuzie007 Parse "Atomic configuration" in POTCAR (52 and 54). Useful for estimating a reasonable NBANDS value.
* EnumerateStructureTransformation now supports `m3gnet_relax` or `m3gnet_static` options.
