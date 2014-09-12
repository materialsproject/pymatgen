Change log
==========

v3.0.0
------
* Pymatgen is now completely Python 2.7 and Python 3.x compatible!
* Spglib and pyhull have been updated to support Python 3.x.
* Completely rewritten pure python cifio module (courtesy of William Davidson
  Richards) removed dependency on PyCIFRW, which has been causing many issues
  with installation.
* Structure and Molecule now supports a very convenient to() and from_str and
  from_file functionality. Instead of trying to load the appropriate parser,
  you can output and read from the appropriate formats directly. See example
  usage.
* ~50% speedup to LinearAssignment code.
* Continuous integration and contribution guidelines now include Python 3.
* **Backwards incompatible changes**
* matgenie.py has now been renamed simply "pmg" for brevity.
* All deprecated methods in pymatgen 2.x have been removed. E.g.,
  pymatgen.core.structure_modifier is no longer available.
* Pymatgen classes now uses the as_dict() method protocol implemented in the
  Monty package. The to_dict property is now deprecated and will be removed
  in pymatgen v3.1.
* Update main docs page examples with the new Structure to, from formats.
