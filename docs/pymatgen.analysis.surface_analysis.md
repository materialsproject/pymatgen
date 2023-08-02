---
layout: default
title: pymatgen.analysis.surface_analysis.md
nav_exclude: true
---

# pymatgen.analysis.surface_analysis module

This module defines tools to analyze surface and adsorption related
quantities as well as related plots. If you use this module, please
consider citing the following works:

```default
R. Tran, Z. Xu, B. Radhakrishnan, D. Winston, W. Sun, K. A. Persson,
S. P. Ong, "Surface Energies of Elemental Crystals", Scientific
Data, 2016, 3:160080, doi: 10.1038/sdata.2016.80.

and

Kang, S., Mo, Y., Ong, S. P., & Ceder, G. (2014). Nanoscale
stabilization of sodium oxides: Implications for Na-O2 batteries.
Nano Letters, 14(2), 1016-1020. https://doi.org/10.1021/nl404557w

and

Montoya, J. H., & Persson, K. A. (2017). A high-throughput framework
    for determining adsorption energies on solid surfaces. Npj
    Computational Materials, 3(1), 14.
    https://doi.org/10.1038/s41524-017-0017-z
```

Todo:
- Still assumes individual elements have their own chempots

> in a molecular adsorbate instead of considering a single
> chempot for a single molecular adsorbate. E.g. for an OH
> adsorbate, the surface energy is a function of delu_O and
> delu_H instead of delu_OH


* Need a method to automatically get chempot range when

    dealing with non-stoichiometric slabs


* Simplify the input for SurfaceEnergyPlotter such that the

    user does not need to generate a dict


### _class_ pymatgen.analysis.surface_analysis.NanoscaleStability(se_analyzers, symprec=1e-05)
Bases: `object`

A class for analyzing the stability of nanoparticles of different

    polymorphs with respect to size. The Wulff shape will be the
    model for the nanoparticle. Stability will be determined by
    an energetic competition between the weighted surface energy
    (surface energy of the Wulff shape) and the bulk energy. A
    future release will include a 2D phase diagram (e.g. wrt size
    vs chempot for adsorbed or non-stoichiometric surfaces). Based
    on the following work:

    Kang, S., Mo, Y., Ong, S. P., & Ceder, G. (2014). Nanoscale

        stabilization of sodium oxides: Implications for Na-O2
        batteries. Nano Letters, 14(2), 1016-1020.
        [https://doi.org/10.1021/nl404557w](https://doi.org/10.1021/nl404557w)


#### se_analyzers()
List of SurfaceEnergyPlotter objects. Each item corresponds to a

    different polymorph.


#### symprec()
See WulffShape.

Analyzes the nanoscale stability of different polymorphs.


#### _static_ bulk_gform(bulk_entry)
Returns the formation energy of the bulk.


* **Parameters**

    **bulk_entry** ([*ComputedStructureEntry*](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedStructureEntry)) – Entry of the corresponding bulk.



* **Returns**

    bulk formation energy (in eV)



* **Return type**

    float



#### plot_all_stability_map(max_r, increments=50, delu_dict=None, delu_default=0, plt=None, labels=None, from_sphere_area=False, e_units='keV', r_units='nanometers', normalize=False, scale_per_atom=False)
Returns the plot of the formation energy of a particles

    of different polymorphs against its effect radius.


* **Parameters**


    * **max_r** (*float*) – The maximum radius of the particle to plot up to.


    * **increments** (*int*) – Number of plot points


    * **delu_dict** (*dict*) – Dictionary of the chemical potentials to be set as
    constant. Note the key should be a sympy Symbol object of the
    format: Symbol(“delu_el”) where el is the name of the element.


    * **delu_default** (*float*) – Default value for all unset chemical potentials


    * **plt** (*pyplot*) – Plot


    * **labels** (*list*) – List of labels for each plot, corresponds to the
    list of se_analyzers


    * **from_sphere_area** (*bool*) – There are two ways to calculate the bulk
    formation energy. Either by treating the volume and thus surface
    area of the particle as a perfect sphere, or as a Wulff shape.



#### plot_one_stability_map(analyzer, max_r, delu_dict=None, label='', increments=50, delu_default=0, plt=None, from_sphere_area=False, e_units='keV', r_units='nanometers', normalize=False, scale_per_atom=False)
Returns the plot of the formation energy of a particle against its

    effect radius.


* **Parameters**


    * **analyzer** (*SurfaceEnergyPlotter*) – Analyzer associated with the
    first polymorph


    * **max_r** (*float*) – The maximum radius of the particle to plot up to.


    * **delu_dict** (*dict*) – Dictionary of the chemical potentials to be set as
    constant. Note the key should be a sympy Symbol object of the
    format: Symbol(“delu_el”) where el is the name of the element.


    * **label** (*str*) – Label of the plot for legend


    * **increments** (*int*) – Number of plot points


    * **delu_default** (*float*) – Default value for all unset chemical potentials


    * **plt** (*pyplot*) – Plot


    * **from_sphere_area** (*bool*) – There are two ways to calculate the bulk
    formation energy. Either by treating the volume and thus surface
    area of the particle as a perfect sphere, or as a Wulff shape.


    * **r_units** (*str*) – Can be nanometers or Angstrom


    * **e_units** (*str*) – Can be keV or eV


    * **normalize** (*str*) – Whether or not to normalize energy by volume



#### scaled_wulff(wulffshape, r)
Scales the Wulff shape with an effective radius r. Note that the resulting

    Wulff does not necessarily have the same effective radius as the one
    provided. The Wulff shape is scaled by its surface energies where first
    the surface energies are scale by the minimum surface energy and then
    multiplied by the given effective radius.


* **Parameters**


    * **wulffshape** ([*WulffShape*](pymatgen.analysis.wulff.md#pymatgen.analysis.wulff.WulffShape)) – Initial, unscaled WulffShape


    * **r** (*float*) – Arbitrary effective radius of the WulffShape



* **Returns**

    WulffShape (scaled by r)



#### solve_equilibrium_point(analyzer1, analyzer2, delu_dict=None, delu_default=0, units='nanometers')
Gives the radial size of two particles where equilibrium is reached

    between both particles. NOTE: the solution here is not the same
    as the solution visualized in the plot because solving for r
    requires that both the total surface area and volume of the
    particles are functions of r.


* **Parameters**


    * **analyzer1** (*SurfaceEnergyPlotter*) – Analyzer associated with the
    first polymorph


    * **analyzer2** (*SurfaceEnergyPlotter*) – Analyzer associated with the
    second polymorph


    * **delu_dict** (*dict*) – Dictionary of the chemical potentials to be set as
    constant. Note the key should be a sympy Symbol object of the
    format: Symbol(“delu_el”) where el is the name of the element.


    * **delu_default** (*float*) – Default value for all unset chemical potentials


    * **units** (*str*) – Can be nanometers or Angstrom



* **Returns**

    Particle radius in nm



#### wulff_gform_and_r(wulffshape, bulk_entry, r, from_sphere_area=False, r_units='nanometers', e_units='keV', normalize=False, scale_per_atom=False)
Calculates the formation energy of the particle with arbitrary radius r.


* **Parameters**


    * **wulffshape** ([*WulffShape*](pymatgen.analysis.wulff.md#pymatgen.analysis.wulff.WulffShape)) – Initial, unscaled WulffShape


    * **bulk_entry** ([*ComputedStructureEntry*](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedStructureEntry)) – Entry of the corresponding bulk.


    * **r** (*float** (**Ang**)*) – Arbitrary effective radius of the WulffShape


    * **from_sphere_area** (*bool*) – There are two ways to calculate the bulk
    formation energy. Either by treating the volume and thus surface
    area of the particle as a perfect sphere, or as a Wulff shape.


    * **r_units** (*str*) – Can be nanometers or Angstrom


    * **e_units** (*str*) – Can be keV or eV


    * **normalize** (*bool*) – Whether or not to normalize energy by volume


    * **scale_per_atom** (*True*) – Whether or not to normalize by number of
    atoms in the particle



* **Returns**

    particle formation energy (float in keV), effective radius



### _class_ pymatgen.analysis.surface_analysis.SlabEntry(structure, energy, miller_index, correction=0.0, parameters=None, data=None, entry_id=None, label=None, adsorbates=None, clean_entry=None, marker=None, color=None)
Bases: [`ComputedStructureEntry`](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedStructureEntry)

A ComputedStructureEntry object encompassing all data relevant to a

    slab for analyzing surface thermodynamics.


#### miller_index()
Miller index of plane parallel to surface.


#### label()
Brief description for this slab.


#### adsorbates()
List of ComputedStructureEntry for the types of adsorbates

..attribute:: clean_entry

> SlabEntry for the corresponding clean slab for an adsorbed slab

..attribute:: ads_entries_dict

> Dictionary where the key is the reduced composition of the

>     adsorbate entry and value is the entry itself

Make a SlabEntry containing all relevant surface thermodynamics data.


* **Parameters**


    * **structure** ([*Slab*](pymatgen.core.surface.md#pymatgen.core.surface.Slab)) – The primary slab associated with this entry.


    * **energy** (*float*) – Energy from total energy calculation


    * **miller_index** (*tuple**(**h**, **k**, **l**)*) – Miller index of plane parallel
    to surface


    * **correction** (*float*) – See ComputedSlabEntry


    * **parameters** (*dict*) – See ComputedSlabEntry


    * **data** (*dict*) – See ComputedSlabEntry


    * **entry_id** (*str*) – See ComputedSlabEntry


    * **data** – See ComputedSlabEntry


    * **entry_id** – See ComputedSlabEntry


    * **label** (*str*) – Any particular label for this slab, e.g. “Tasker 2”,
    “non-stoichiometric”, “reconstructed”


    * **adsorbates** (*[*[*ComputedStructureEntry*](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedStructureEntry)*]*) – List of reference entries
    for the adsorbates on the slab, can be an isolated molecule
    (e.g. O2 for O or O2 adsorption), a bulk structure (eg. fcc
    Cu for Cu adsorption) or anything.


    * **clean_entry** ([*ComputedStructureEntry*](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedStructureEntry)) – If the SlabEntry is for an
    adsorbed slab, this is the corresponding SlabEntry for the
    clean slab


    * **marker** (*str*) – Custom marker for gamma plots (”–” and “-” are typical)


    * **color** (*str** or **rgba*) – Custom color for gamma plots



#### _property_ Nads_in_slab()
Returns the TOTAL number of adsorbates in the slab on BOTH sides.


#### _property_ Nsurfs_ads_in_slab()
Returns the TOTAL number of adsorbed surfaces in the slab.


#### as_dict()
Returns dict which contains Slab Entry data.


#### _property_ cleaned_up_slab()
Returns a slab with the adsorbates removed.


#### _property_ create_slab_label()
Returns a label (str) for this particular slab based on composition, coverage and Miller index.


#### _static_ from_computed_structure_entry(entry, miller_index, label=None, adsorbates=None, clean_entry=None, \*\*kwargs)
Returns SlabEntry from a ComputedStructureEntry.


#### _classmethod_ from_dict(dct)
Returns a SlabEntry by reading in an dictionary.


#### _property_ get_monolayer()
Returns the primitive unit surface area density of the
adsorbate.


#### _property_ get_unit_primitive_area()
Returns the surface area of the adsorbed system per
unit area of the primitive slab system.


#### gibbs_binding_energy(eads=False)
Returns the adsorption energy or Gibbs binding energy

    of an adsorbate on a surface


* **Parameters**

    **eads** (*bool*) – Whether to calculate the adsorption energy
    (True) or the binding energy (False) which is just
    adsorption energy normalized by number of adsorbates.



#### _property_ surface_area()
Calculates the surface area of the slab.


#### surface_energy(ucell_entry, ref_entries=None)
Calculates the surface energy of this SlabEntry.


* **Parameters**


    * **ucell_entry** (*entry*) – An entry object for the bulk


    * **(****list** (*ref_entries*) – [entry]): A list of entries for each type
    of element to be used as a reservoir for non-stoichiometric
    systems. The length of this list MUST be n-1 where n is the
    number of different elements in the bulk entry. The chempot
    of the element ref_entry that is not in the list will be
    treated as a variable.


Returns (Add (Sympy class)): Surface energy


### _class_ pymatgen.analysis.surface_analysis.SurfaceEnergyPlotter(all_slab_entries, ucell_entry, ref_entries=None)
Bases: `object`

A class used for generating plots to analyze the thermodynamics of surfaces

    of a material. Produces stability maps of different slab configurations,
    phases diagrams of two parameters to determine stability of configurations
    (future release), and Wulff shapes.


#### all_slab_entries()
Either a list of SlabEntry objects (note for a list, the SlabEntry must

    have the adsorbates and clean_entry parameter pulgged in) or a Nested
    dictionary containing a list of entries for slab calculations as
    items and the corresponding Miller index of the slab as the key.
    To account for adsorption, each value is a sub-dictionary with the
    entry of a clean slab calculation as the sub-key and a list of
    entries for adsorption calculations as the sub-value. The sub-value
    can contain different adsorption configurations such as a different
    site or a different coverage, however, ordinarily only the most stable
    configuration for a particular coverage will be considered as the
    function of the adsorbed surface energy has an intercept dependent on
    the adsorption energy (ie an adsorption site with a higher adsorption
    energy will always provide a higher surface energy than a site with a
    lower adsorption energy). An example parameter is provided:
    {(h1,k1,l1): {clean_entry1: [ads_entry1, ads_entry2, …],

    > clean_entry2: […], …}, (h2,k2,l2): {…}}

    where clean_entry1 can be a pristine surface and clean_entry2 can be a
    reconstructed surface while ads_entry1 can be adsorption at site 1 with
    a 2x2 coverage while ads_entry2 can have a 3x3 coverage. If adsorption
    entries are present (i.e. if all_slab_entries[(h,k,l)][clean_entry1]), we
    consider adsorption in all plots and analysis for this particular facet.

..attribute:: color_dict

> Dictionary of colors (r,g,b,a) when plotting surface energy stability. The

>     keys are individual surface entries where clean surfaces have a solid
>     color while the corresponding adsorbed surface will be transparent.


#### ucell_entry()
ComputedStructureEntry of the bulk reference for this particular material.


#### ref_entries()
List of ComputedStructureEntries to be used for calculating chemical potential.


#### color_dict()
Randomly generated dictionary of colors associated with each facet.

Object for plotting surface energy in different ways for clean and

    adsorbed surfaces.


* **Parameters**


    * **all_slab_entries** (*dict** or **list*) – Dictionary or list containing
    all entries for slab calculations. See attributes.


    * **ucell_entry** ([*ComputedStructureEntry*](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedStructureEntry)) – ComputedStructureEntry
    of the bulk reference for this particular material.


    * **ref_entries** (*[**ComputedStructureEntries**]*) – A list of entries for
    each type of element to be used as a reservoir for
    non-stoichiometric systems. The length of this list MUST be
    n-1 where n is the number of different elements in the bulk
    entry. The bulk energy term in the grand surface potential can
    be defined by a summation of the chemical potentials for each
    element in the system. As the bulk energy is already provided,
    one can solve for one of the chemical potentials as a function
    of the other chemical potetinals and bulk energy. i.e. there
    are n-1 variables (chempots). e.g. if your ucell_entry is for
    LiFePO4 than your ref_entries should have an entry for Li, Fe,
    and P if you want to use the chempot of O as the variable.



#### BE_vs_clean_SE(delu_dict, delu_default=0, plot_eads=False, annotate_monolayer=True, JPERM2=False)
For each facet, plot the clean surface energy against the most

    stable binding energy.


* **Parameters**


    * **delu_dict** (*dict*) – Dictionary of the chemical potentials to be set as
    constant. Note the key should be a sympy Symbol object of the
    format: Symbol(“delu_el”) where el is the name of the element.


    * **delu_default** (*float*) – Default value for all unset chemical potentials


    * **plot_eads** (*bool*) – Option to plot the adsorption energy (binding
    energy multiplied by number of adsorbates) instead.


    * **annotate_monolayer** (*bool*) – Whether or not to label each data point
    with its monolayer (adsorbate density per unit primiitve area)


    * **JPERM2** (*bool*) – Whether to plot surface energy in /m^2 (True) or
    eV/A^2 (False)



* **Returns**

    Plot of clean surface energy vs binding energy for

        all facets.




* **Return type**

    (Plot)



#### area_frac_vs_chempot_plot(ref_delu, chempot_range, delu_dict=None, delu_default=0, increments=10, no_clean=False, no_doped=False)
1D plot. Plots the change in the area contribution
of each facet as a function of chemical potential.


* **Parameters**


    * **ref_delu** (*sympy Symbol*) – The free variable chempot with the format:
    Symbol(“delu_el”) where el is the name of the element.


    * **chempot_range** (*list*) – Min/max range of chemical potential to plot along


    * **delu_dict** (*dict*) – Dictionary of the chemical potentials to be set as
    constant. Note the key should be a sympy Symbol object of the
    format: Symbol(“delu_el”) where el is the name of the element.


    * **delu_default** (*float*) – Default value for all unset chemical potentials


    * **increments** (*int*) – Number of data points between min/max or point
    of intersection. Defaults to 10 points.



* **Returns**

    Plot of area frac on the Wulff shape

        for each facet vs chemical potential.




* **Return type**

    (pyplot)



#### _static_ chempot_plot_addons(plt, xrange, ref_el, axes, pad=2.4, rect=None, ylim=None)
Helper function to a chempot plot look nicer.


* **Parameters**


    * **plt** (*Plot*) –


    * **xrange** (*list*) – xlim parameter


    * **ref_el** (*str*) – Element of the referenced chempot.


    * **axes** (*axes*) –


    * **pad** (*float*) –


    * **rect** (*list*) – For tight layout


    * **ylim** (*ylim parameter*) –


return (Plot): Modified plot with addons.
return (Plot): Modified plot with addons.


#### chempot_vs_gamma(ref_delu, chempot_range, miller_index=(), delu_dict=None, delu_default=0, JPERM2=False, show_unstable=False, ylim=None, plt=None, no_clean=False, no_doped=False, use_entry_labels=False, no_label=False)
Plots the surface energy as a function of chemical potential.

    Each facet will be associated with its own distinct colors.
    Dashed lines will represent stoichiometries different from that
    of the mpid’s compound. Transparent lines indicates adsorption.


* **Parameters**


    * **ref_delu** (*sympy Symbol*) – The range stability of each slab is based
    on the chempot range of this chempot. Should be a sympy Symbol
    object of the format: Symbol(“delu_el”) where el is the name of
    the element


    * **chempot_range** (*[**max_chempot**, **min_chempot**]*) – Range to consider the
    stability of the slabs.


    * **miller_index** (*list*) – Miller index for a specific facet to get a
    dictionary for.


    * **delu_dict** (*dict*) – Dictionary of the chemical potentials to be set as
    constant. Note the key should be a sympy Symbol object of the
    format: Symbol(“delu_el”) where el is the name of the element.


    * **delu_default** (*float*) – Default value for all unset chemical potentials


    * **JPERM2** (*bool*) – Whether to plot surface energy in /m^2 (True) or
    eV/A^2 (False)


    * **show_unstable** (*bool*) – Whether or not to show parts of the surface
    energy plot outside the region of stability.


    * **ylim** (*[**ymax**, **ymin**]*) – Range of y axis


    * **no_doped** (*bool*) – Whether to plot for the clean slabs only.


    * **no_clean** (*bool*) – Whether to plot for the doped slabs only.


    * **use_entry_labels** (*bool*) – If True, will label each slab configuration
    according to their given label in the SlabEntry object.


    * **no_label** (*bool*) – Option to turn off labels.



* **Returns**

    Plot of surface energy vs chempot for all entries.



* **Return type**

    (Plot)



#### chempot_vs_gamma_plot_one(plt, entry, ref_delu, chempot_range, delu_dict=None, delu_default=0, label='', JPERM2=False)
Helper function to  help plot the surface energy of a
single SlabEntry as a function of chemical potential.


* **Parameters**


    * **plt** (*Plot*) – A plot.


    * **entry** (*SlabEntry*) – Entry of the slab whose surface energy we want
    to plot


    * **ref_delu** (*sympy Symbol*) – The range stability of each slab is based
    on the chempot range of this chempot. Should be a sympy Symbol
    object of the format: Symbol(“delu_el”) where el is the name of
    the element


    * **chempot_range** (*[**max_chempot**, **min_chempot**]*) – Range to consider the
    stability of the slabs.


    * **delu_dict** (*dict*) – Dictionary of the chemical potentials to be set as
    constant. Note the key should be a sympy Symbol object of the
    format: Symbol(“delu_el”) where el is the name of the element.


    * **delu_default** (*float*) – Default value for all unset chemical potentials


    * **label** (*str*) – Label of the slab for the legend.


    * **JPERM2** (*bool*) – Whether to plot surface energy in /m^2 (True) or
    eV/A^2 (False)



* **Returns**

    Plot of surface energy vs chemical potential for one entry.



* **Return type**

    (Plot)



#### color_palette_dict(alpha=0.35)
Helper function to assign each facet a unique color using a dictionary.


* **Parameters**

    **alpha** (*float*) – Degree of transparency


return (dict): Dictionary of colors (r,g,b,a) when plotting surface

    energy stability. The keys are individual surface entries where
    clean surfaces have a solid color while the corresponding adsorbed
    surface will be transparent.


#### get_stable_entry_at_u(miller_index, delu_dict=None, delu_default=0, no_doped=False, no_clean=False)
Returns the entry corresponding to the most stable slab for a particular

    facet at a specific chempot. We assume that surface energy is constant
    so all free variables must be set with delu_dict, otherwise they are
    assumed to be equal to delu_default.


* **Parameters**


    * **miller_index** (*(**h**,**k**,**l**)*) – The facet to find the most stable slab in


    * **delu_dict** (*dict*) – Dictionary of the chemical potentials to be set as
    constant. Note the key should be a sympy Symbol object of the
    format: Symbol(“delu_el”) where el is the name of the element.


    * **delu_default** (*float*) – Default value for all unset chemical potentials


    * **no_doped** (*bool*) – Consider stability of clean slabs only.


    * **no_clean** (*bool*) – Consider stability of doped slabs only.



* **Returns**

    SlabEntry, surface_energy (float)



#### get_surface_equilibrium(slab_entries, delu_dict=None)
Takes in a list of SlabEntries and calculates the chemical potentials

    at which all slabs in the list coexists simultaneously. Useful for
    building surface phase diagrams. Note that to solve for x equations
    (x slab_entries), there must be x free variables (chemical potentials).
    Adjust delu_dict as need be to get the correct number of free variables.


* **Parameters**


    * **slab_entries** (*array*) – The coefficients of the first equation


    * **delu_dict** (*dict*) – Dictionary of the chemical potentials to be set as
    constant. Note the key should be a sympy Symbol object of the
    format: Symbol(“delu_el”) where el is the name of the element.



* **Returns**

    Array containing a solution to x equations with x

        variables (x-1 chemical potential and 1 surface energy)




* **Return type**

    (array)



#### monolayer_vs_BE(plot_eads=False)
Plots the binding energy as a function of monolayers (ML), i.e.

    the fractional area adsorbate density for all facets. For each
    facet at a specific monlayer, only plot the lowest binding energy.


* **Parameters**

    **plot_eads** (*bool*) – Option to plot the adsorption energy (binding
    energy multiplied by number of adsorbates) instead.



* **Returns**

    Plot of binding energy vs monolayer for all facets.



* **Return type**

    (Plot)



#### set_all_variables(delu_dict, delu_default)
Sets all chemical potential values and returns a dictionary where

    the key is a sympy Symbol and the value is a float (chempot).


* **Parameters**


    * **entry** (*SlabEntry*) – Computed structure entry of the slab


    * **delu_dict** (*dict*) – Dictionary of the chemical potentials to be set as
    constant. Note the key should be a sympy Symbol object of the
    format: Symbol(“delu_el”) where el is the name of the element.


    * **delu_default** (*float*) – Default value for all unset chemical potentials



* **Returns**

    Dictionary of set chemical potential values



#### stable_u_range_dict(chempot_range, ref_delu, no_doped=True, no_clean=False, delu_dict=None, miller_index=(), dmu_at_0=False, return_se_dict=False)
Creates a dictionary where each entry is a key pointing to a
chemical potential range where the surface of that entry is stable.
Does so by enumerating through all possible solutions (intersect)
for surface energies of a specific facet.


* **Parameters**


    * **chempot_range** (*[**max_chempot**, **min_chempot**]*) – Range to consider the
    stability of the slabs.


    * **ref_delu** (*sympy Symbol*) – The range stability of each slab is based
    on the chempot range of this chempot. Should be a sympy Symbol
    object of the format: Symbol(“delu_el”) where el is the name of
    the element


    * **no_doped** (*bool*) – Consider stability of clean slabs only.


    * **no_clean** (*bool*) – Consider stability of doped slabs only.


    * **delu_dict** (*dict*) – Dictionary of the chemical potentials to be set as
    constant. Note the key should be a sympy Symbol object of the
    format: Symbol(“delu_el”) where el is the name of the element.


    * **miller_index** (*list*) – Miller index for a specific facet to get a
    dictionary for.


    * **dmu_at_0** (*bool*) – If True, if the surface energies corresponding to
    the chemical potential range is between a negative and positive
    value, the value is a list of three chemical potentials with the
    one in the center corresponding a surface energy of 0. Uselful
    in identifying unphysical ranges of surface energies and their
    chemical potential range.


    * **return_se_dict** (*bool*) – Whether or not to return the corresponding
    dictionary of surface energies



#### surface_chempot_range_map(elements, miller_index, ranges, incr=50, no_doped=False, no_clean=False, delu_dict=None, plt=None, annotate=True, show_unphyiscal_only=False, fontsize=10)
Adapted from the get_chempot_range_map() method in the PhaseDiagram

    class. Plot the chemical potential range map based on surface
    energy stability. Currently works only for 2-component PDs. At
    the moment uses a brute force method by enumerating through the
    range of the first element chempot with a specified increment
    and determines the chempot rangeo fht e second element for each
    SlabEntry. Future implementation will determine the chempot range
    map first by solving systems of equations up to 3 instead of 2.


* **Parameters**


    * **elements** (*list*) – Sequence of elements to be considered as independent
    variables. E.g., if you want to show the stability ranges of
    all Li-Co-O phases wrt to duLi and duO, you will supply
    [Element(“Li”), Element(“O”)]


    * **miller_index** (*[**h**, **k**, **l**]*) – Miller index of the surface we are interested in


    * **ranges** (*[**[**range1**]**, **[**range2**]**]*) – List of chempot ranges (max and min values)
    for the first and second element.


    * **incr** (*int*) – Number of points to sample along the range of the first chempot


    * **no_doped** (*bool*) – Whether or not to include doped systems.


    * **no_clean** (*bool*) – Whether or not to include clean systems.


    * **delu_dict** (*dict*) – Dictionary of the chemical potentials to be set as
    constant. Note the key should be a sympy Symbol object of the
    format: Symbol(“delu_el”) where el is the name of the element.


    * **plt** (*Plot*) – Plot object to plot on. If None, will create a new plot.


    * **annotate** (*bool*) – Whether to annotate each “phase” with the label of
    the entry. If no label, uses the reduced formula


    * **show_unphyiscal_only** (*bool*) – Whether to only show the shaded region where
    surface energy is negative. Useful for drawing other chempot range maps.


    * **fontsize** (*int*) – Font size of the annotation



#### wulff_from_chempot(delu_dict=None, delu_default=0, symprec=1e-05, no_clean=False, no_doped=False)
Method to get the Wulff shape at a specific chemical potential.


* **Parameters**


    * **delu_dict** (*dict*) – Dictionary of the chemical potentials to be set as
    constant. Note the key should be a sympy Symbol object of the
    format: Symbol(“delu_el”) where el is the name of the element.


    * **delu_default** (*float*) – Default value for all unset chemical potentials


    * **symprec** (*float*) – See WulffShape.


    * **no_doped** (*bool*) – Consider stability of clean slabs only.


    * **no_clean** (*bool*) – Consider stability of doped slabs only.



* **Returns**

    The WulffShape at u_ref and u_ads.



* **Return type**

    ([WulffShape](pymatgen.analysis.wulff.md#pymatgen.analysis.wulff.WulffShape))



### _class_ pymatgen.analysis.surface_analysis.WorkFunctionAnalyzer(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), locpot_along_c, efermi, shift=0, blength=3.5)
Bases: `object`

A class used for calculating the work function

    from a slab model and visualizing the behavior
    of the local potential along the slab.


#### efermi()
The Fermi energy


#### locpot_along_c()
Local potential in eV along points along the  axis


#### vacuum_locpot()
The maximum local potential along the c direction for

    the slab model, ie the potential at the vacuum


#### work_function()
The minimum energy needed to move an electron from the

    surface to infinity. Defined as the difference between
    the potential at the vacuum and the Fermi energy.


#### slab()
The slab structure model


#### along_c()
Points along the c direction with same

    increments as the locpot in the c axis


#### ave_locpot()
Mean of the minimum and maximmum (vacuum) locpot along c


#### sorted_sites()
List of sites from the slab sorted along the c direction


#### ave_bulk_p()
The average locpot of the slab region along the c direction

Initializes the WorkFunctionAnalyzer class.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Structure object modelling the surface


    * **locpot_along_c** (*list*) – Local potential along the c direction


    * **outcar** (*MSONable*) – Outcar vasp output object


    * **shift** (*float*) – Parameter to translate the slab (and
    therefore the vacuum) of the slab structure, thereby
    translating the plot along the x axis.


    * **blength** (*float** (**Ang**)*) – The longest bond length in the material.
    Used to handle pbc for noncontiguous slab layers



#### _static_ from_files(poscar_filename, locpot_filename, outcar_filename, shift=0, blength=3.5)

* **Parameters**


    * **poscar_filename** – POSCAR file


    * **locpot_filename** – LOCPOT file


    * **outcar_filename** – OUTCAR file


    * **shift** – shift


    * **blength** – The longest bond length in the material.
    Used to handle pbc for noncontiguous slab layers



* **Returns**

    WorkFunctionAnalyzer



#### get_labels(plt, label_fontsize=10)
Handles the optional labelling of the plot with relevant quantities
:param plt: Plot of the locpot vs c axis
:type plt: plt
:param label_fontsize: Fontsize of labels
:type label_fontsize: float

Returns Labelled plt.


#### get_locpot_along_slab_plot(label_energies=True, plt=None, label_fontsize=10)
Returns a plot of the local potential (eV) vs the

    position along the c axis of the slab model (Ang).


* **Parameters**


    * **label_energies** (*bool*) – Whether to label relevant energy
    quantities such as the work function, Fermi energy,
    vacuum locpot, bulk-like locpot


    * **plt** (*plt*) – Matplotlib pyplot object


    * **label_fontsize** (*float*) – Fontsize of labels


Returns plt of the locpot vs c axis


#### is_converged(min_points_frac=0.015, tol: float = 0.0025)
A well converged work function should have a flat electrostatic

    potential within some distance (min_point) about where the peak
    electrostatic potential is found along the c direction of the
    slab. This is dependent on the size of the slab.


* **Parameters**


    * **min_point** (*fractional coordinates*) – The number of data points
    +/- the point of where the electrostatic potential is at
    its peak along the c direction.


    * **tol** (*float*) – If the electrostatic potential stays the same
    within this tolerance, within the min_points, it is converged.


Returns a bool (whether or not the work function is converged)


### pymatgen.analysis.surface_analysis.entry_dict_from_list(all_slab_entries)
Converts a list of SlabEntry to an appropriate dictionary. It is
assumed that if there is no adsorbate, then it is a clean SlabEntry
and that adsorbed SlabEntry has the clean_entry parameter set.


* **Parameters**

    **all_slab_entries** (*list*) – List of SlabEntry objects



* **Returns**

    Dictionary of SlabEntry with the Miller index as the main

        key to a dictionary with a clean SlabEntry as the key to a
        list of adsorbed SlabEntry.




* **Return type**

    (dict)



### pymatgen.analysis.surface_analysis.sub_chempots(gamma_dict, chempots)
Uses dot product of numpy array to sub chemical potentials

    into the surface grand potential. This is much faster
    than using the subs function in sympy.


* **Parameters**


    * **gamma_dict** (*dict*) – Surface grand potential equation
    as a coefficient dictionary


    * **chempots** (*dict*) – Dictionary assigning each chemical
    potential (key) in gamma a value



* **Returns**

    Surface energy as a float