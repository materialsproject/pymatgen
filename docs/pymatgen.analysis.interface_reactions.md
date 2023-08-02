---
layout: default
title: pymatgen.analysis.interface_reactions.md
nav_exclude: true
---

# pymatgen.analysis.interface_reactions module

This module provides a class to predict and analyze interfacial reactions between two
solids, with or without an open element (e.g., flowing O2).


### _class_ pymatgen.analysis.interface_reactions.GrandPotentialInterfacialReactivity(c1: [Composition](pymatgen.core.composition.md#pymatgen.core.composition.Composition), c2: [Composition](pymatgen.core.composition.md#pymatgen.core.composition.Composition), grand_pd: [GrandPotentialPhaseDiagram](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.GrandPotentialPhaseDiagram), pd_non_grand: [PhaseDiagram](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PhaseDiagram), include_no_mixing_energy: bool = False, norm: bool = True, use_hull_energy: bool = True)
Bases: `InterfacialReactivity`

Extends upon InterfacialReactivity to allow for modelling possible reactions
at the interface between two solids in the presence of an open element. The
thermodynamics of the open system are provided by the user via the
GrandPotentialPhaseDiagram class.


* **Parameters**


    * **c1** – Reactant 1 composition


    * **c2** – Reactant 2 composition


    * **grand_pd** – Grand potential phase diagram object built from all elements in
    composition c1 and c2.


    * **include_no_mixing_energy** – No_mixing_energy for a reactant is the
    opposite number of its energy above grand potential convex hull. In
    cases where reactions involve elements reservoir, this param
    determines whether no_mixing_energy of reactants will be included
    in the final reaction energy calculation. By definition, if pd is
    not a GrandPotentialPhaseDiagram object, this param is False.


    * **pd_non_grand** – PhaseDiagram object but not
    GrandPotentialPhaseDiagram object built from elements in c1 and c2.


    * **norm** – Whether or not the total number of atoms in composition
    of reactant will be normalized to 1.


    * **use_hull_energy** – Whether or not use the convex hull energy for
    a given composition for reaction energy calculation. If false,
    the energy of ground state structure will be used instead.
    Note that in case when ground state can not be found for a
    composition, convex hull energy will be used associated with a
    warning message.



#### get_no_mixing_energy()
Generates the opposite number of energy above grand potential
convex hull for both reactants.


* **Returns**

    [(reactant1, no_mixing_energy1),(reactant2,no_mixing_energy2)].



### _class_ pymatgen.analysis.interface_reactions.InterfacialReactivity(c1: [Composition](pymatgen.core.composition.md#pymatgen.core.composition.Composition), c2: [Composition](pymatgen.core.composition.md#pymatgen.core.composition.Composition), pd: [PhaseDiagram](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PhaseDiagram), norm: bool = True, use_hull_energy: bool = False, \*\*kwargs)
Bases: `MSONable`

Class for modeling an interface between two solids and its possible reactions.
The two reactants are provided as Composition objects (c1 and c2), along with the
relevant compositional PhaseDiagram object. Possible reactions are calculated by
finding all points along a tie-line between c1 and c2 where there is a “kink” in
the phase diagram; i.e. a point or facet of the phase diagram.

Please consider citing one or both of the following papers if you use this code
in your own work.

### References

Richards, W. D., Miara, L. J., Wang, Y., Kim, J. C., &amp; Ceder, G. (2015).
Interface stability in solid-state batteries. Chemistry of Materials, 28(1),
266-273. [https://doi.org/10.1021/acs.chemmater.5b04082](https://doi.org/10.1021/acs.chemmater.5b04082)

Xiao, Y., Wang, Y., Bo, S.-H., Kim, J. C., Miara, L. J., &amp; Ceder, G. (2019).
Understanding interface stability in solid-state batteries.
Nature Reviews Materials, 5(2), 105-126.
[https://doi.org/10.1038/s41578-019-0157-5](https://doi.org/10.1038/s41578-019-0157-5)


* **Parameters**


    * **c1** – Reactant 1 composition


    * **c2** – Reactant 2 composition


    * **pd** – Phase diagram object built from all elements in composition c1 and c2.


    * **norm** – Whether or not the total number of atoms in composition
    of reactant will be normalized to 1.


    * **use_hull_energy** – Whether or not use the convex hull energy for
    a given composition for reaction energy calculation. If false,
    the energy of ground state structure will be used instead.
    Note that in case when ground state can not be found for a
    composition, convex hull energy will be used associated with a
    warning message.



#### EV_TO_KJ_PER_MOL(_ = 96.485_ )

#### _classmethod_ get_chempot_correction(element: str, temp: float, pres: float)
Get the normalized correction term Δμ for chemical potential of a gas
phase consisting of element at given temperature and pressure,
referenced to that in the standard state (T_std = 298.15 K,
T_std = 1 bar). The gas phase is limited to be one of O2, N2, Cl2,
F2, H2. Calculation formula can be found in the documentation of
Materials Project website.


* **Parameters**


    * **element** – The string representing the element.


    * **temp** – The temperature of the gas phase in Kelvin.


    * **pres** – The pressure of the gas phase in Pa.



* **Returns**

    The correction of chemical potential in eV/atom of the gas
    phase at given temperature and pressure.



#### get_critical_original_kink_ratio()
Returns a list of molar mixing ratio for each kink between ORIGINAL
(instead of processed) reactant compositions. This is the
same list as mixing ratio obtained from get_kinks method
if self.norm = False.


* **Returns**

    A list of floats representing molar mixing ratios between
    the original reactant compositions for each kink.



#### get_dataframe()
Returns a pandas DataFrame representation of the data produced by the
get_kinks() method.


#### get_kinks()
Finds all the kinks in mixing ratio where reaction products changes
along the tie-line of composition self.c1 and composition self.c2.


* **Returns**

    (index, mixing ratio, reaction energy in eV/atom, Reaction object, reaction
    energy per mol of formula in kJ/mol).



* **Return type**

    List object of tuples, each of which contains 5 elements



#### _property_ labels()
Returns a dictionary containing kink information:
{index: ‘x= mixing_ratio energy= reaction_energy reaction_equation’}.
E.g., {1: ‘x= 0 energy = 0 Mn -> Mn’,

> 2: ‘x= 0.5 energy = -15 O2 + Mn -> MnO2’,
> 3: ‘x= 1 energy = 0 O2 -> O2’}.


#### _property_ minimum()
Finds the minimum reaction energy E_min and corresponding
mixing ratio x_min.


* **Returns**

    Tuple (x_min, E_min).



#### plot(backend: Literal['plotly', 'matplotlib'] = 'plotly')
Plots reaction energy as a function of mixing ratio x in self.c1 - self.c2
tie line.


* **Parameters**

    **backend** (*"plotly"** | **"matplotlib"*) – Plotting library used to create the plot. Defaults to
    “plotly” but can also be “matplotlib”.



* **Returns**

    Plot of reaction energies as a function of mixing ratio



#### _property_ products()
List of formulas of potential products. E.g., [‘Li’,’O2’,’Mn’].