---
layout: default
title: pymatgen.analysis.chemenv.utils.defs_utils.md
nav_exclude: true
---

# pymatgen.analysis.chemenv.utils.defs_utils module

This module contains the definition of some objects used in the chemenv package.


### _class_ pymatgen.analysis.chemenv.utils.defs_utils.AdditionalConditions()
Bases: `object`

Class for additional conditions.


#### ALL(_ = (0, 1, 2, 3, 4_ )

#### CONDITION_DESCRIPTION(_ = {0: 'No additional condition', 1: 'Only anion-cation bonds', 2: 'No element-element bonds (same elements)', 3: 'Only anion-cation bonds and no element-element bonds (same elements)', 4: 'Only element-oxygen bonds'_ )

#### NONE(_ = _ )

#### NO_AC(_ = _ )

#### NO_ADDITIONAL_CONDITION(_ = _ )

#### NO_E2SEB(_ = _ )

#### NO_ELEMENT_TO_SAME_ELEMENT_BONDS(_ = _ )

#### ONLY_ACB(_ = _ )

#### ONLY_ACB_AND_NO_E2SEB(_ = _ )

#### ONLY_ANION_CATION_BONDS(_ = _ )

#### ONLY_ANION_CATION_BONDS_AND_NO_ELEMENT_TO_SAME_ELEMENT_BONDS(_ = _ )

#### ONLY_E2OB(_ = _ )

#### ONLY_ELEMENT_TO_OXYGEN_BONDS(_ = _ )

#### check_condition(condition, structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), parameters)

* **Parameters**


    * **condition** –


    * **structure** –


    * **parameters** –



* **Returns**