---
layout: default
title: pymatgen.analysis.bond_valence.md
nav_exclude: true
---

# pymatgen.analysis.bond_valence module

This module implements classes to perform bond valence analyses.


### _class_ pymatgen.analysis.bond_valence.BVAnalyzer(symm_tol=0.1, max_radius=4, max_permutations=100000, distance_scale_factor=1.015, charge_neutrality_tolerance=1e-05, forbidden_species=None)
Bases: `object`

This class implements a maximum a posteriori (MAP) estimation method to
determine oxidation states in a structure. The algorithm is as follows:
1) The bond valence sum of all symmetrically distinct sites in a structure
is calculated using the element-based parameters in M. O’Keefe, & N. Brese,
JACS, 1991, 113(9), 3226-3229. doi:10.1021/ja00009a002.
2) The posterior probabilities of all oxidation states is then calculated
using: P(oxi_state/BV) = K \* P(BV/oxi_state) \* P(oxi_state), where K is
a constant factor for each element. P(BV/oxi_state) is calculated as a
Gaussian with mean and std deviation determined from an analysis of
the ICSD. The posterior P(oxi_state) is determined from a frequency
analysis of the ICSD.
3) The oxidation states are then ranked in order of decreasing probability
and the oxidation state combination that result in a charge neutral cell
is selected.

Initializes the BV analyzer, with useful defaults.


* **Parameters**


    * **symm_tol** – Symmetry tolerance used to determine which sites are
    symmetrically equivalent. Set to 0 to turn off symmetry.


    * **max_radius** – Maximum radius in Angstrom used to find nearest neighbors.


    * **max_permutations** – The maximum number of permutations of oxidation states to test.


    * **distance_scale_factor** – A scale factor to be applied. This is useful for scaling
    distances, esp in the case of calculation-relaxed structures
    which may tend to under (GGA) or over bind (LDA). The default
    of 1.015 works for GGA. For experimental structure, set this to
    1.


    * **charge_neutrality_tolerance** – Tolerance on the charge neutrality when unordered structures
    are at stake.


    * **forbidden_species** – List of species that are forbidden (example : [“O-”] cannot be
    used) It is used when e.g. someone knows that some oxidation
    state cannot occur for some atom in a structure or list of
    structures.



#### CHARGE_NEUTRALITY_TOLERANCE(_ = 1e-0_ )

#### get_oxi_state_decorated_structure(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))
Get an oxidation state decorated structure. This currently works only
for ordered structures only.


* **Parameters**

    **structure** – Structure to analyze



* **Returns**

    A modified structure that is oxidation state decorated.



* **Raises**

    **ValueError if the valences cannot be determined.** –



#### get_valences(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))
Returns a list of valences for each site in the structure.


* **Parameters**

    **structure** – Structure to analyze



* **Returns**

    A list of valences for each site in the structure (for an ordered structure),
    e.g., [1, 1, -2] or a list of lists with the valences for each fractional
    element of each site in the structure (for an unordered structure), e.g., [[2,
    4], [3], [-2], [-2], [-2]]



* **Raises**

    **A ValueError if the valences cannot be determined.** –



### pymatgen.analysis.bond_valence.add_oxidation_state_by_site_fraction(structure, oxidation_states)
Add oxidation states to a structure by fractional site.


* **Parameters**

    **oxidation_states** (*list*) – List of list of oxidation states for each
    site fraction for each site.
    E.g., [[2, 4], [3], [-2], [-2], [-2]]



### pymatgen.analysis.bond_valence.calculate_bv_sum(site, nn_list, scale_factor=1.0)
Calculates the BV sum of a site.


* **Parameters**


    * **site** ([*PeriodicSite*](pymatgen.core.sites.md#pymatgen.core.sites.PeriodicSite)) – The central site to calculate the bond valence


    * **nn_list** (*[*[*Neighbor*](pymatgen.core.structure.md#pymatgen.core.structure.Neighbor)*]*) – A list of namedtuple Neighbors having “distance”
    and “site” attributes


    * **scale_factor** (*float*) – A scale factor to be applied. This is useful for
    scaling distance, esp in the case of calculation-relaxed structures
    which may tend to under (GGA) or over bind (LDA).



### pymatgen.analysis.bond_valence.calculate_bv_sum_unordered(site, nn_list, scale_factor=1)
Calculates the BV sum of a site for unordered structures.


* **Parameters**


    * **site** ([*PeriodicSite*](pymatgen.core.sites.md#pymatgen.core.sites.PeriodicSite)) – The central site to calculate the bond valence


    * **nn_list** (*[*[*Neighbor*](pymatgen.core.structure.md#pymatgen.core.structure.Neighbor)*]*) – A list of namedtuple Neighbors having “distance”
    and “site” attributes


    * **scale_factor** (*float*) – A scale factor to be applied. This is useful for
    scaling distance, esp in the case of calculation-relaxed structures
    which may tend to under (GGA) or over bind (LDA).



### pymatgen.analysis.bond_valence.get_z_ordered_elmap(comp)
Arbitrary ordered element map on the elements/species of a composition of a
given site in an unordered structure. Returns a list of tuples (
element_or_specie: occupation) in the arbitrary order.

The arbitrary order is based on the Z of the element and the smallest
fractional occupations first.
Example : {“Ni3+”: 0.2, “Ni4+”: 0.2, “Cr3+”: 0.15, “Zn2+”: 0.34,
“Cr4+”: 0.11} will yield the species in the following order :
Cr4+, Cr3+, Ni3+, Ni4+, Zn2+ … or
Cr4+, Cr3+, Ni4+, Ni3+, Zn2+