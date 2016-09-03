Change log
==========

v4.2.2
------
* Global configuration variables such as VASP\_PSP\_DIR and MAPI\_KEY are now
  stored in "~/.pmgrc.yaml". If you are setting these as environmental
  variables right now, you can easily transition to the new system using
  "pmg config --add VASP\_PSP\_DIR $VASP\_PSP\_DIR MAPI\_KEY $MAPI\_KEY".
  This new scheme will provide greater flexibility for user-defined
  global behavior in pymatgen, e.g., tolerances, default input sets for
  transmuters, etc., in future.
* Beta of k-point weight calculator.
* Use default MSONable as and from_dict for all transformations.
