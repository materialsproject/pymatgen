---
layout: default
title: pymatgen.analysis.pourbaix_diagram.md
nav_exclude: true
---

# pymatgen.analysis.pourbaix_diagram module

This module is intended to be used to compute Pourbaix diagrams of arbitrary compositions
and formation energies.


### _class_ pymatgen.analysis.pourbaix_diagram.IonEntry(ion, energy, name=None, attribute=None)
Bases: [`PDEntry`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PDEntry)

Object similar to PDEntry, but contains an Ion object instead of a
Composition object.


#### name()
A name for the entry. This is the string shown in the phase diagrams.
By default, this is the reduced formula for the composition, but can be
set to some other string for display purposes.


* **Parameters**


    * **ion** – Ion object


    * **energy** – Energy for composition.


    * **name** – Optional parameter to name the entry. Defaults to the
    chemical formula.



#### as_dict()
Creates a dict of composition, energy, and ion name.


#### _classmethod_ from_dict(d)
Returns an IonEntry object from a dict.


### _class_ pymatgen.analysis.pourbaix_diagram.MultiEntry(entry_list, weights=None)
Bases: `PourbaixEntry`

PourbaixEntry-like object for constructing multi-elemental Pourbaix
diagrams.

Initializes a MultiEntry.


* **Parameters**


    * **entry_list** (*[**PourbaixEntry**]*) – List of component PourbaixEntries


    * **weights** (*[**float**]*) – Weights associated with each entry. Default is None



#### as_dict()
Returns: MSONable dict.


#### _classmethod_ from_dict(d)

* **Parameters**

    **(****)** (*d*) – Dict representation.



* **Returns**

    MultiEntry



#### _property_ name()
MultiEntry name, i. e. the name of each entry joined by ‘ + ‘.


### _class_ pymatgen.analysis.pourbaix_diagram.PourbaixDiagram(entries: list[PourbaixEntry] | list[MultiEntry], comp_dict: dict[str, float] | None = None, conc_dict: dict[str, float] | None = None, filter_solids: bool = True, nproc: int | None = None)
Bases: `MSONable`

Class to create a Pourbaix diagram from entries.


* **Parameters**


    * **entries** (*[**PourbaixEntry**] or **[**MultiEntry**]*) – Entries list
    containing Solids and Ions or a list of MultiEntries


    * **comp_dict** (*dict**[**str**, **float**]*) – Dictionary of compositions,
    defaults to equal parts of each elements


    * **conc_dict** (*dict**[**str**, **float**]*) – Dictionary of ion concentrations,
    defaults to 1e-6 for each element


    * **filter_solids** (*bool*) – applying this filter to a Pourbaix
    diagram ensures all included solid phases are filtered by
    stability on the compositional phase diagram. Defaults to True.
    The practical consequence of this is that highly oxidized or reduced
    phases that might show up in experiments due to kinetic limitations
    on oxygen/hydrogen evolution won’t appear in the diagram, but they are
    not actually “stable” (and are frequently overstabilized from DFT errors).
    Hence, including only the stable solid phases generally leads to the
    most accurate Pourbaix diagrams.


    * **nproc** (*int*) – number of processes to generate multientries with
    in parallel. Defaults to None (serial processing).



#### _property_ all_entries()
Return all entries used to generate the Pourbaix diagram.


#### as_dict()

* **Returns**

    MSONable dict.



#### find_stable_entry(pH, V)
Finds stable entry at a pH,V condition
:param pH: pH to find stable entry
:type pH: float
:param V: V to find stable entry.
:type V: float

Returns:


#### _classmethod_ from_dict(d)

* **Parameters**

    **(****)** (*d*) – Dict representation.



* **Returns**

    PourbaixDiagram



#### get_decomposition_energy(entry, pH, V)
Finds decomposition to most stable entries in eV/atom,
supports vectorized inputs for pH and V.


* **Parameters**


    * **entry** (*PourbaixEntry*) – PourbaixEntry corresponding to
    compound to find the decomposition for


    * **pH** (*float**, **[**float**]*) – pH at which to find the decomposition


    * **V** (*float**, **[**float**]*) – voltage at which to find the decomposition



* **Returns**

    Decomposition energy for the entry, i. e. the energy above

        the “Pourbaix hull” in eV/atom at the given conditions




#### get_hull_energy(pH, V)
Gets the minimum energy of the Pourbaix “basin” that is formed
from the stable Pourbaix planes. Vectorized.


* **Parameters**


    * **pH** (*float** or **[**float**]*) – pH at which to find the hull energy


    * **V** (*float** or **[**float**]*) – V at which to find the hull energy



* **Returns**

    (float or [float]) minimum Pourbaix energy at conditions



#### _static_ get_pourbaix_domains(pourbaix_entries, limits=None)
Returns a set of Pourbaix stable domains (i. e. polygons) in
pH-V space from a list of pourbaix_entries.

This function works by using scipy’s HalfspaceIntersection
function to construct all of the 2-D polygons that form the
boundaries of the planes corresponding to individual entry
gibbs free energies as a function of pH and V. Hyperplanes
of the form a\*pH + b\*V + 1 - g(0, 0) are constructed and
supplied to HalfspaceIntersection, which then finds the
boundaries of each Pourbaix region using the intersection
points.


* **Parameters**


    * **pourbaix_entries** (*[**PourbaixEntry**]*) – Pourbaix entries
    with which to construct stable Pourbaix domains


    * **limits** (*[**[**float**]**]*) – limits in which to do the pourbaix
    analysis



* **Returns**

    [boundary_points]}.
    The list of boundary points are the sides of the N-1
    dim polytope bounding the allowable ph-V range of each entry.



* **Return type**

    Returns a dict of the form {entry



#### get_stable_entry(pH, V)
Gets the stable entry at a given pH, V condition.


* **Parameters**


    * **pH** (*float*) – pH at a given condition


    * **V** (*float*) – V at a given condition



* **Returns**

    Pourbaix or multi-entry

        corresponding ot the minimum energy entry at a given
        pH, V condition




* **Return type**

    (PourbaixEntry or MultiEntry)



#### _static_ process_multientry(entry_list, prod_comp, coeff_threshold=0.0001)
Static method for finding a multientry based on
a list of entries and a product composition.
Essentially checks to see if a valid aqueous
reaction exists between the entries and the
product composition and returns a MultiEntry
with weights according to the coefficients if so.


* **Parameters**


    * **entry_list** (*[*[*Entry*](pymatgen.entries.md#pymatgen.entries.Entry)*]*) – list of entries from which to
    create a MultiEntry


    * **prod_comp** ([*Composition*](pymatgen.core.composition.md#pymatgen.core.composition.Composition)) – composition constraint for setting
    weights of MultiEntry


    * **coeff_threshold** (*float*) – threshold of stoichiometric
    coefficients to filter, if weights are lower than
    this value, the entry is not returned



#### _property_ stable_entries()
Returns the stable entries in the Pourbaix diagram.


#### _property_ unprocessed_entries()
Return unprocessed entries.


#### _property_ unstable_entries()
Returns all unstable entries in the Pourbaix diagram.


### _class_ pymatgen.analysis.pourbaix_diagram.PourbaixEntry(entry, entry_id=None, concentration=1e-06)
Bases: `MSONable`, [`Stringify`](pymatgen.util.string.md#pymatgen.util.string.Stringify)

An object encompassing all data relevant to a solid or ion
in a Pourbaix diagram. Each bulk solid/ion has an energy
g of the form: e = e0 + 0.0591 log10(conc) - nO mu_H2O
+ (nH - 2nO) pH + phi (-nH + 2nO + q).

Note that the energies corresponding to the input entries
should be formation energies with respect to hydrogen and
oxygen gas in order for the Pourbaix diagram formalism to
work. This may be changed to be more flexible in the future.


* **Parameters**


    * **entry** (*ComputedEntry/ComputedStructureEntry/PDEntry/IonEntry*) – An
    entry object


    * **(****)** (*concentration*) –


    * **(****)** –



#### as_dict()
Returns dict which contains Pourbaix Entry data.
Note that the pH, voltage, H2O factors are always calculated when
constructing a PourbaixEntry object.


#### _property_ composition()
Returns composition.


#### _property_ conc_term()
Returns the concentration contribution to the free energy,
and should only be present when there are ions in the entry.


#### _property_ energy()
total energy of the Pourbaix
entry (at pH, V = 0 vs. SHE).


* **Type**

    Returns (float)



#### energy_at_conditions(pH, V)
Get free energy for a given pH and V.


* **Parameters**


    * **pH** (*float*) – pH at which to evaluate free energy


    * **V** (*float*) – voltage at which to evaluate free energy



* **Returns**

    free energy at conditions



#### _property_ energy_per_atom()
energy per atom of the Pourbaix entry.

Returns (float): energy per atom


#### _classmethod_ from_dict(d)
Invokes a PourbaixEntry from a dictionary.


#### get_element_fraction(element)
Gets the elemental fraction of a given non-OH element.


* **Parameters**

    **element** ([*Element*](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element)* or **str*) – string or element corresponding
    to element to get from composition



* **Returns**

    fraction of element / sum(all non-OH elements)



#### _property_ nH2O()
Get the number of H2O.


#### _property_ nPhi()
Get the number of electrons.


#### _property_ name()
Get the name for entry.


#### _property_ normalization_factor()
Sum of number of atoms minus the number of H and O in composition.


#### _property_ normalized_energy()
Returns:
energy normalized by number of non H or O atoms, e. g.
for Zn2O6, energy / 2 or for AgTe3(OH)3, energy / 4.


#### normalized_energy_at_conditions(pH, V)
Energy at an electrochemical condition, compatible with
numpy arrays for pH/V input.


* **Parameters**


    * **pH** (*float*) – pH at condition


    * **V** (*float*) – applied potential at condition



* **Returns**

    energy normalized by number of non-O/H atoms at condition



#### _property_ npH()
Get the number of H.


#### _property_ num_atoms()
Return number of atoms in current formula. Useful for normalization.


#### to_pretty_string()

* **Returns**

    A pretty string representation.



### _class_ pymatgen.analysis.pourbaix_diagram.PourbaixPlotter(pourbaix_diagram)
Bases: `object`

A plotter class for phase diagrams.


* **Parameters**

    **pourbaix_diagram** (*PourbaixDiagram*) – A PourbaixDiagram object.



#### domain_vertices(entry)
Returns the vertices of the Pourbaix domain.


* **Parameters**

    **entry** – Entry for which domain vertices are desired



* **Returns**

    list of vertices



#### get_pourbaix_plot(limits=None, title='', label_domains=True, label_fontsize=20, show_water_lines=True, show_neutral_axes=True, plt=None)
Plot Pourbaix diagram.


* **Parameters**


    * **limits** – 2D list containing limits of the Pourbaix diagram
    of the form [[xlo, xhi], [ylo, yhi]]


    * **title** (*str*) – Title to display on plot


    * **label_domains** (*bool*) – whether to label Pourbaix domains


    * **label_fontsize** – font size for domain labels


    * **show_water_lines** – whether to show dashed lines indicating the region
    of water stability.


    * **lines** (*show_neutral_axes; whether to show dashed horizontal and vertical*) – at 0 V and pH 7, respectively.


    * **plt** (*pyplot*) – Pyplot instance for plotting



* **Returns**

    plt (pyplot) - matplotlib plot object with Pourbaix diagram



#### plot_entry_stability(entry, pH_range=None, pH_resolution=100, V_range=None, V_resolution=100, e_hull_max=1, cmap='RdYlBu_r', \*\*kwargs)

* **Parameters**


    * **(****)** (*\*\*kwargs*) –


    * **(****)** –


    * **(****)** –


    * **(****)** –


    * **(****)** –


    * **(****)** –


    * **(****)** –


    * **(****)** –


Returns:


#### show(\*args, \*\*kwargs)
Shows the Pourbaix plot.


* **Parameters**


    * **\*args** – args to get_pourbaix_plot


    * **\*\*kwargs** – kwargs to get_pourbaix_plot



* **Returns**

    None



### pymatgen.analysis.pourbaix_diagram.generate_entry_label(entry)
Generates a label for the Pourbaix plotter.


* **Parameters**

    **entry** (*PourbaixEntry** or **MultiEntry*) – entry to get a label for



### pymatgen.analysis.pourbaix_diagram.ion_or_solid_comp_object(formula)
Returns either an ion object or composition object given
a formula.


* **Parameters**

    **formula** – String formula. Eg. of ion: NaOH(aq), Na[+];
    Eg. of solid: Fe2O3(s), Fe(s), Na2O



* **Returns**

    Composition/Ion object