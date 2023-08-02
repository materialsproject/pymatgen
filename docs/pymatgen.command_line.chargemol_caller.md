---
layout: default
title: pymatgen.command_line.chargemol_caller.md
nav_exclude: true
---

# pymatgen.command_line.chargemol_caller module

This module implements an interface to Thomas Manz’s Chargemol code
[https://sourceforge.net/projects/ddec](https://sourceforge.net/projects/ddec) for calculating DDEC3, DDEC6, and CM5 population analyses.

This module depends on a compiled chargemol executable being available in the path.
If you use this module, please cite the following based on which modules you use:

Chargemol:
(1) T. A. Manz and N. Gabaldon Limas, Chargemol program for performing DDEC analysis,
Version 3.5, 2017, ddec.sourceforge.net.

DDEC6 Charges:
(1) T. A. Manz and N. Gabaldon Limas, “Introducing DDEC6 atomic population analysis:
part 1. Charge partitioning theory and methodology,” RSC Adv., 6 (2016) 47771-47801.
(2) N. Gabaldon Limas and T. A. Manz, “Introducing DDEC6 atomic population analysis:
part 2. Computed results for a wide range of periodic and nonperiodic materials,”
(3) N. Gabaldon Limas and T. A. Manz, “Introducing DDEC6 atomic population analysis:
part 4. Efficient parallel computation of net atomic charges, atomic spin moments,
bond orders, and more,” RSC Adv., 8 (2018) 2678-2707.

CM5 Charges:
(1) A.V. Marenich, S.V. Jerome, C.J. Cramer, D.G. Truhlar, “Charge Model 5: An Extension
of Hirshfeld Population Analysis for the Accurate Description of Molecular Interactions
in Gaseous and Condensed Phases”, J. Chem. Theory. Comput., 8 (2012) 527-541.

Spin Moments:
(1) T. A. Manz and D. S. Sholl, “Methods for Computing Accurate Atomic Spin Moments for
Collinear and Noncollinear Magnetism in Periodic and Nonperiodic Materials,”
J. Chem. Theory Comput. 7 (2011) 4146-4164.

Bond Orders:
(1) “Introducing DDEC6 atomic population analysis: part 3. Comprehensive method to compute
bond orders,” RSC Adv., 7 (2017) 45552-45581.

DDEC3 Charges:
(1) T. A. Manz and D. S. Sholl, “Improved Atoms-in-Molecule Charge Partitioning Functional
for Simultaneously Reproducing the Electrostatic Potential and Chemical States in Periodic
and Non-Periodic Materials,” J. Chem. Theory Comput. 8 (2012) 2844-2867.
(2) T. A. Manz and D. S. Sholl, “Chemically Meaningful Atomic Charges that Reproduce the
Electrostatic Potential in Periodic and Nonperiodic Materials,” J. Chem. Theory Comput. 6
(2010) 2455-2468.


### _class_ pymatgen.command_line.chargemol_caller.ChargemolAnalysis(path=None, atomic_densities_path=None, run_chargemol=True)
Bases: `object`

Chargemol analysis for DDEC3, DDEC6, and/or CM5 population analyses,
including the calculation of partial atomic charges, atomic spin moments,
bond orders, and related properties.

Initializes the Chargemol Analysis.


* **Parameters**


    * **path** (*str*) – Path to the CHGCAR, POTCAR, AECCAR0, and AECCAR files.


    * **not.** (*Note that it doesn't matter if the files gzip'd or*) – Default: None (current working directory).


    * **atomic_densities_path** (*str**|**None*) – Path to the atomic densities directory


    * **None** (*required by Chargemol. If*) –


    * **is** (*Pymatgen assumes that this*) –


    * **variable.** (*defined in a "DDEC6_ATOMIC_DENSITIES_DIR" environment*) –


    * **True.** (*Only used if run_chargemol is*) – Default: None.


    * **run_chargemol** (*bool*) – Whether to run the Chargemol analysis. If False,


    * **path.** (*the existing Chargemol output files will be read from*) – Default: True.



#### get_bond_order(index_from, index_to)
Convenience method to get the bond order between two atoms.


* **Parameters**


    * **index_from** (*int*) – Index of atom to get bond order from.


    * **index_to** (*int*) – Index of atom to get bond order to.



* **Returns**

    bond order between atoms



* **Return type**

    float



#### get_charge(atom_index, nelect=None, charge_type='ddec')
Convenience method to get the charge on a particular atom using the same
sign convention as the BaderAnalysis. Note that this is *not* the partial
atomic charge. This value is nelect (e.g. ZVAL from the POTCAR) + the
charge transferred. If you want the partial atomic charge, use
get_partial_charge().


* **Parameters**


    * **atom_index** (*int*) – Index of atom to get charge for.


    * **nelect** (*int*) – number of electrons associated with an isolated atom at this index.


    * **electrons** (*For most DFT codes this corresponds to the number** of **valence*) –


    * **None** (*associated with the pseudopotential. If*) –


    * **automatically** (*this value will be*) –


    * **POTCAR** (*obtained from the*) – Default: None.


    * **charge_type** (*str*) – Type of charge to use (“ddec” or “cm5”).



* **Returns**

    charge on atom_index



* **Return type**

    float



#### get_charge_transfer(atom_index, charge_type='ddec')
Returns the charge transferred for a particular atom. A positive value means
that the site has gained electron density (i.e. exhibits anionic character)
whereas a negative value means the site has lost electron density (i.e. exhibits
cationic character). This is the same thing as the negative of the partial atomic
charge.


* **Parameters**


    * **atom_index** (*int*) – Index of atom to get charge transfer for.


    * **charge_type** (*str*) – Type of charge to use (“ddec” or “cm5”).



* **Returns**

    charge transferred at atom_index



* **Return type**

    float



#### get_partial_charge(atom_index, charge_type='ddec')
Convenience method to get the partial atomic charge on a particular atom.
This is the value printed in the Chargemol analysis.


* **Parameters**


    * **atom_index** (*int*) – Index of atom to get charge for.


    * **charge_type** (*str*) – Type of charge to use (“ddec” or “cm5”).



#### get_property_decorated_structure()
Takes CHGCAR’s structure object and updates it with properties
from the Chargemol analysis.


* **Returns**

    Pymatgen structure with site properties added



#### _property_ summary()
Returns a dictionary summary of the Chargemol analysis
{

> “ddec”: {

>     > “partial_charges”: List[float],
>     > “spin_moments”: List[float],
>     > “dipoles”: List[float],
>     > “rsquared_moments”: List[float],
>     > “rcubed_moments”: List[float],
>     > “rfourth_moments”: List[float],
>     > “bond_order_dict”: Dict

>     },

> “cm5”: {

>     > “partial_charges”: List[float],

>     }

}.