---
layout: default
title: pymatgen.command_line.enumlib_caller.md
nav_exclude: true
---

# pymatgen.command_line.enumlib_caller module

This module implements an interface to enumlib, Gus Hart’s excellent Fortran
code for enumerating derivative structures.

This module depends on a compiled enumlib with the executables enum.x and
makestr.x available in the path. Please download the library at
[https://github.com/msg-byu/enumlib](https://github.com/msg-byu/enumlib) and follow the instructions in the README to
compile these two executables accordingly.

If you use this module, please cite:

Gus L. W. Hart and Rodney W. Forcade, “Algorithm for generating derivative
structures,” Phys. Rev. B 77 224115 (26 June 2008)

Gus L. W. Hart and Rodney W. Forcade, “Generating derivative structures from
multilattices: Application to hcp alloys,” Phys. Rev. B 80 014120 (July 2009)

Gus L. W. Hart, Lance J. Nelson, and Rodney W. Forcade, “Generating
derivative structures at a fixed concentration,” Comp. Mat. Sci. 59
101-107 (March 2012)

Wiley S. Morgan, Gus L. W. Hart, Rodney W. Forcade, “Generating derivative
superstructures for systems with high configurational freedom,” Comp. Mat.
Sci. 136 144-149 (May 2017)


### _exception_ pymatgen.command_line.enumlib_caller.EnumError()
Bases: `BaseException`

Error subclass for enumeration errors.