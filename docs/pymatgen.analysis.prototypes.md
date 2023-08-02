---
layout: default
title: pymatgen.analysis.prototypes.md
nav_exclude: true
---

# pymatgen.analysis.prototypes module

This module is intended to match crystal structures against known crystallographic “prototype”
structures.

In this module, the AflowPrototypeMatcher uses the AFLOW LIBRARY OF CRYSTALLOGRAPHIC PROTOTYPES.
If using this particular class, please cite their publication appropriately:

Mehl, M. J., Hicks, D., Toher, C., Levy, O., Hanson, R. M., Hart, G., & Curtarolo, S. (2017).
The AFLOW library of crystallographic prototypes: part 1.
Computational Materials Science, 136, S1-S828.
[https://doi.org/10.1016/j.commatsci.2017.01.017](https://doi.org/10.1016/j.commatsci.2017.01.017)


### _class_ pymatgen.analysis.prototypes.AflowPrototypeMatcher(initial_ltol=0.2, initial_stol=0.3, initial_angle_tol=5)
Bases: `object`

This class will match structures to their crystal prototypes, and will
attempt to group species together to match structures derived from
prototypes (e.g. an A_xB_1-x_C from a binary prototype), and will
give these the names the “-like” suffix.

This class uses data from the AFLOW LIBRARY OF CRYSTALLOGRAPHIC PROTOTYPES.
If using this class, please cite their publication appropriately:

Mehl, M. J., Hicks, D., Toher, C., Levy, O., Hanson, R. M., Hart, G., & Curtarolo, S. (2017).
The AFLOW library of crystallographic prototypes: part 1.
Computational Materials Science, 136, S1-S828.
[https://doi.org/10.1016/j.commatsci.2017.01.017](https://doi.org/10.1016/j.commatsci.2017.01.017)

Tolerances as defined in StructureMatcher. Tolerances will be
gradually decreased until only a single match is found (if possible).


* **Parameters**


    * **initial_ltol** – fractional length tolerance


    * **initial_stol** – site tolerance


    * **initial_angle_tol** – angle tolerance



#### get_prototypes(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))
Get prototype(s) structures for a given input structure. If you use this method in
your work, please cite the appropriate AFLOW publication:

Mehl, M. J., Hicks, D., Toher, C., Levy, O., Hanson, R. M., Hart, G., & Curtarolo,
S. (2017). The AFLOW library of crystallographic prototypes: part 1. Computational
Materials Science, 136, S1-S828. [https://doi.org/10.1016/j.commatsci.2017.01.017](https://doi.org/10.1016/j.commatsci.2017.01.017)


* **Parameters**

    **structure** – structure to match


Returns (list): A list of dicts with keys ‘snl’ for the matched prototype and

    ‘tags’, a dict of tags (‘mineral’, ‘strukturbericht’ and ‘aflow’) of that
    prototype. This should be a list containing just a single entry, but it is
    possible a material can match multiple prototypes.