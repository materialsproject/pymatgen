---
layout: default
title: pymatgen.analysis.phase_diagram.md
nav_exclude: true
---

# pymatgen.analysis.phase_diagram module

This module defines tools to generate and analyze phase diagrams.


### _class_ pymatgen.analysis.phase_diagram.CompoundPhaseDiagram(entries, terminal_compositions, normalize_terminal_compositions=True)
Bases: `PhaseDiagram`

Generates phase diagrams from compounds as terminations instead of
elements.

Initializes a CompoundPhaseDiagram.


* **Parameters**


    * **entries** (*[**PDEntry**]*) – Sequence of input entries. For example,
    if you want a Li2O-P2O5 phase diagram, you might have all
    Li-P-O entries as an input.


    * **terminal_compositions** (*[*[*Composition*](pymatgen.core.composition.md#pymatgen.core.composition.Composition)*]*) – Terminal compositions of
    phase space. In the Li2O-P2O5 example, these will be the
    Li2O and P2O5 compositions.


    * **normalize_terminal_compositions** (*bool*) – Whether to normalize the
    terminal compositions to a per atom basis. If normalized,
    the energy above hulls will be consistent
    for comparison across systems. Non-normalized terminals are
    more intuitive in terms of compositional breakdowns.



#### amount_tol(_ = 1e-0_ )

#### as_dict()

* **Returns**

    MSONable dictionary representation of CompoundPhaseDiagram.



#### _classmethod_ from_dict(d)

* **Parameters**

    **d** (*dict*) – dictionary representation of CompoundPhaseDiagram.



* **Returns**

    CompoundPhaseDiagram



#### transform_entries(entries, terminal_compositions)
Method to transform all entries to the composition coordinate in the
terminal compositions. If the entry does not fall within the space
defined by the terminal compositions, they are excluded. For example,
Li3PO4 is mapped into a Li2O:1.5, P2O5:0.5 composition. The terminal
compositions are represented by DummySpecies.


* **Parameters**


    * **entries** – Sequence of all input entries


    * **terminal_compositions** – Terminal compositions of phase space.



* **Returns**

    Sequence of TransformedPDEntries falling within the phase space.



### _class_ pymatgen.analysis.phase_diagram.GrandPotPDEntry(entry, chempots, name=None)
Bases: `PDEntry`

A grand potential pd entry object encompassing all relevant data for phase
diagrams. Chemical potentials are given as a element-chemical potential
dict.


* **Parameters**


    * **entry** – A PDEntry-like object.


    * **chempots** – Chemical potential specification as {Element: float}.


    * **name** – Optional parameter to name the entry. Defaults to the reduced
    chemical formula of the original entry.



#### as_dict()

* **Returns**

    MSONable dictionary representation of GrandPotPDEntry.



#### _property_ chemical_energy()
The chemical energy term mu\*N in the grand potential.


* **Returns**

    The chemical energy term mu\*N in the grand potential



#### _property_ composition(_: [Composition](pymatgen.core.composition.md#pymatgen.core.composition.Composition_ )
The composition after removing free species.


* **Returns**

    Composition



#### _property_ energy()
Returns:
The grand potential energy.


#### _classmethod_ from_dict(d)

* **Parameters**

    **d** (*dict*) – dictionary representation of GrandPotPDEntry.



* **Returns**

    GrandPotPDEntry



### _class_ pymatgen.analysis.phase_diagram.GrandPotentialPhaseDiagram(entries, chempots, elements=None, \*, computed_data=None)
Bases: `PhaseDiagram`

A class representing a Grand potential phase diagram. Grand potential phase
diagrams are essentially phase diagrams that are open to one or more
components. To construct such phase diagrams, the relevant free energy is
the grand potential, which can be written as the Legendre transform of the
Gibbs free energy as follows.

Grand potential = G - u_X N_X

The algorithm is based on the work in the following papers:


1. S. P. Ong, L. Wang, B. Kang, and G. Ceder, Li-Fe-P-O2 Phase Diagram from
First Principles Calculations. Chem. Mater., 2008, 20(5), 1798-1807.
doi:10.1021/cm702327g


2. S. P. Ong, A. Jain, G. Hautier, B. Kang, G. Ceder, Thermal stabilities
of delithiated olivine MPO4 (M=Fe, Mn) cathodes investigated using first
principles calculations. Electrochem. Comm., 2010, 12(3), 427-430.
doi:10.1016/j.elecom.2010.01.010

Standard constructor for grand potential phase diagram.


* **Parameters**


    * **entries** (*[**PDEntry**]*) – A list of PDEntry-like objects having an
    energy, energy_per_atom and composition.


    * **(****{Element** (*chempots*) – float}): Specify the chemical potentials
    of the open elements.


    * **elements** (*[*[*Element*](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element)*]*) – Optional list of elements in the phase
    diagram. If set to None, the elements are determined from
    the entries themselves.


    * **computed_data** (*dict*) – A dict containing pre-computed data. This allows
    PhaseDiagram object to be reconstituted without performing the
    expensive convex hull computation. The dict is the output from the
    PhaseDiagram._compute() method and is stored in PhaseDiagram.computed_data
    when generated for the first time.



#### as_dict()

* **Returns**

    MSONable dictionary representation of GrandPotentialPhaseDiagram.



#### _classmethod_ from_dict(d)

* **Parameters**

    **d** (*dict*) – dictionary representation of GrandPotentialPhaseDiagram.



* **Returns**

    GrandPotentialPhaseDiagram



### _class_ pymatgen.analysis.phase_diagram.PDEntry(composition: [Composition](pymatgen.core.composition.md#pymatgen.core.composition.Composition), energy: float, name: str | None = None, attribute: object = None)
Bases: [`Entry`](pymatgen.entries.md#pymatgen.entries.Entry)

An object encompassing all relevant data for phase diagrams.


#### composition()
The composition associated with the PDEntry.


* **Type**

    [Composition](pymatgen.core.composition.md#pymatgen.core.composition.Composition)



#### energy()
The energy associated with the entry.


* **Type**

    float



#### name()
A name for the entry. This is the string shown in the phase diagrams.
By default, this is the reduced formula for the composition, but can be
set to some other string for display purposes.


* **Type**

    str



#### attribute()
A arbitrary attribute. Can be used to specify that the
entry is a newly found compound, or to specify a particular label for
the entry, etc. An attribute can be anything but must be MSONable.


* **Type**

    MSONable



* **Parameters**


    * **composition** ([*Composition*](pymatgen.core.composition.md#pymatgen.core.composition.Composition)) – Composition


    * **energy** (*float*) – Energy for composition.


    * **name** (*str*) – Optional parameter to name the entry. Defaults
    to the reduced chemical formula.


    * **attribute** – Optional attribute of the entry. Must be MSONable.



#### as_dict()

* **Returns**

    MSONable dictionary representation of PDEntry.



#### _property_ energy(_: floa_ )
Returns:
the energy of the entry.


#### _classmethod_ from_dict(dct)

* **Parameters**

    **dct** (*dict*) – dictionary representation of PDEntry.



* **Returns**

    PDEntry



### _class_ pymatgen.analysis.phase_diagram.PDPlotter(phasediagram: PhaseDiagram, show_unstable: float = 0.2, backend: Literal['plotly', 'matplotlib'] = 'plotly', ternary_style: Literal['2d', '3d'] = '2d', \*\*plotkwargs)
Bases: `object`

A plotting class for compositional phase diagrams.

To use, initialize this class with a PhaseDiagram object containing 1-4 components
and call get_plot() or show().


* **Parameters**


    * **phasediagram** (*PhaseDiagram*) – PhaseDiagram object (must be 1-4 components).


    * **show_unstable** (*float*) – Whether unstable (above the hull) phases will be
    plotted. If a number > 0 is entered, all phases with
    e_hull < show_unstable (eV/atom) will be shown.


    * **backend** (*"plotly"** | **"matplotlib"*) – Python package to use for plotting.
    Defaults to “plotly”.


    * **ternary_style** (*"2d"** | **"3d"*) – Ternary phase diagrams are typically plotted in
    two-dimensions (2d), but can be plotted in three dimensions (3d) to visualize
    the depth of the hull. This argument only applies when backend=”plotly”.
    Defaults to “2d”.


    * **\*\*plotkwargs** (*dict*) – Keyword args passed to matplotlib.pyplot.plot (only
    applies when backend=”matplotlib”). Can be used to customize markers
    etc. If not set, the default is:

    > {

    >     “markerfacecolor”: “#4daf4a”,
    >     “markersize”: 10,
    >     “linewidth”: 3

    > }.




#### get_chempot_range_map_plot(elements, referenced=True)
Returns a plot of the chemical potential range _map. Currently works
only for 3-component PDs.

Note: this functionality is now included in the ChemicalPotentialDiagram
class (pymatgen.analysis.chempot_diagram).


* **Parameters**


    * **elements** – Sequence of elements to be considered as independent
    variables. E.g., if you want to show the stability ranges of
    all Li-Co-O phases wrt to uLi and uO, you will supply
    [Element(“Li”), Element(“O”)]


    * **referenced** – if True, gives the results with a reference being the
    energy of the elemental phase. If False, gives absolute values.



* **Returns**

    A matplotlib plot object.



#### get_contour_pd_plot()
Plot a contour phase diagram plot, where phase triangles are colored
according to degree of instability by interpolation. Currently only
works for 3-component phase diagrams.


* **Returns**

    A matplotlib plot object.



#### get_plot(label_stable: bool = True, label_unstable: bool = True, ordering: Sequence[str] | None = None, energy_colormap=None, process_attributes: bool = False, plt=None, label_uncertainties: bool = False, fill: bool = True, highlight_entries: Collection[PDEntry] | None = None)

* **Parameters**


    * **label_stable** – Whether to label stable compounds.


    * **label_unstable** – Whether to label unstable compounds.


    * **ordering** – Ordering of vertices, given as a list [‘Up’,
    ‘Left’,’Right’] (matplotlib only).


    * **energy_colormap** – Colormap for coloring energy (matplotlib only).


    * **process_attributes** – Whether to process the attributes (matplotlib only).


    * **plt** – Existing matplotlib.pyplot object if plotting multiple phase diagrams
    (matplotlib only).


    * **label_uncertainties** – Whether to add error bars to the hull.
    For binaries, this also shades the hull with the uncertainty window.
    (plotly only).


    * **fill** – Whether to shade the hull. For ternary_2d and quaternary plots, this
    colors facets arbitrarily for visual clarity. For ternary_3d plots, this
    shades the hull by formation energy (plotly only).


    * **highlight_entries** – Entries to highlight in the plot (plotly only). This will
    create a new marker trace that is separate from the other entries.



* **Returns**

    go.Figure (backend=”plotly”) or matplotlib.pyplot (backend=”matplotlib”)



#### _property_ pd_plot_data()
Plotting data for phase diagram. Cached for repetitive calls.

2-comp - Full hull with energies
3/4-comp - Projection into 2D or 3D Gibbs triangles


* **Returns**


    * lines is a list of list of coordinates for lines in the PD.


    * stable_entries is a dict of {coordinates

        in the phase diagram. (Each coordinate can only have one
        stable phase)


    * unstable_entries is a dict of {entry: coordinates} for all unstable

        nodes in the phase diagram.




* **Return type**

    A tuple containing three objects (lines, stable_entries, unstable_entries)



#### plot_chempot_range_map(elements, referenced=True)
Plot the chemical potential range _map using matplotlib. Currently works only for
3-component PDs. This shows the plot but does not return it.

Note: this functionality is now included in the ChemicalPotentialDiagram
class (pymatgen.analysis.chempot_diagram).


* **Parameters**


    * **elements** – Sequence of elements to be considered as independent
    variables. E.g., if you want to show the stability ranges of
    all Li-Co-O phases wrt to uLi and uO, you will supply
    [Element(“Li”), Element(“O”)]


    * **referenced** – if True, gives the results with a reference being the
    energy of the elemental phase. If False, gives absolute values.



#### plot_element_profile(element, comp, show_label_index=None, xlim=5)
Draw the element profile plot for a composition varying different
chemical potential of an element.

X value is the negative value of the chemical potential reference to
elemental chemical potential. For example, if choose Element(“Li”),
X= -(µLi-µLi0), which corresponds to the voltage versus metal anode.
Y values represent for the number of element uptake in this composition
(unit: per atom). All reactions are printed to help choosing the
profile steps you want to show label in the plot.


* **Parameters**


    * **element** ([*Element*](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element)) – An element of which the chemical potential is
    considered. It also must be in the phase diagram.


    * **comp** ([*Composition*](pymatgen.core.composition.md#pymatgen.core.composition.Composition)) – A composition.


    * **show_label_index** (*list** of **integers*) – The labels for reaction products
    you want to show in the plot. Default to None (not showing any
    annotation for reaction products). For the profile steps you want
    to show the labels, just add it to the show_label_index. The
    profile step counts from zero. For example, you can set
    show_label_index=[0, 2, 5] to label profile step 0,2,5.


    * **xlim** (*float*) – The max x value. x value is from 0 to xlim. Default to
    5 eV.



* **Returns**

    Plot of element profile evolution by varying the chemical potential
    of an element.



#### show(\*args, \*\*kwargs)
Draw the phase diagram with the provided arguments and display it. This shows
the figure but does not return it.


* **Parameters**


    * **\*args** – Passed to get_plot.


    * **\*\*kwargs** – Passed to get_plot.



#### write_image(stream: str | StringIO, image_format: str = 'svg', \*\*kwargs)
Directly save the plot to a file. This is a wrapper for calling plt.savefig() or
fig.write_image(), depending on the backend. For more customization, it is
recommended to call those methods directly.


* **Parameters**


    * **stream** (*str** | **StringIO*) – Filename or StringIO stream.


    * **image_format** (*str*) – Can be any supported image format for the plotting backend.
    Defaults to ‘svg’ (vector graphics).


    * **\*\*kwargs** – Optinoal kwargs passed to the get_plot function.



### _class_ pymatgen.analysis.phase_diagram.PatchedPhaseDiagram(entries: Sequence[PDEntry] | set[PDEntry], elements: Sequence[[Element](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element)] | None = None, keep_all_spaces: bool = False, verbose: bool = False)
Bases: `PhaseDiagram`

Computing the Convex Hull of a large set of data in multiple dimensions is
highly expensive. This class acts to breakdown large chemical spaces into
smaller chemical spaces which can be computed much more quickly due to having
both reduced dimensionality and data set sizes.


### subspaces ({str()
{Element, }}): Dictionary of the sets of elements for each of the
PhaseDiagrams within the PatchedPhaseDiagram.


### pds ({str()
PhaseDiagram}): Dictionary of PhaseDiagrams within the
PatchedPhaseDiagram.


#### all_entries()
All entries provided for Phase Diagram construction.
Note that this does not mean that all these entries are actually used in
the phase diagram. For example, this includes the positive formation energy
entries that are filtered out before Phase Diagram construction.


* **Type**

    list[PDEntry]



#### min_entries()
List of the  lowest energy entries for each composition
in the data provided for Phase Diagram construction.


* **Type**

    list[PDEntry]



#### el_refs()
List of elemental references for the phase diagrams.
These are entries corresponding to the lowest energy element entries for
simple compositional phase diagrams.


* **Type**

    list[PDEntry]



#### elements()
List of elements in the phase diagram.


* **Type**

    list[[Element](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element)]



* **Parameters**


    * **entries** (*list**[**PDEntry**]*) – A list of PDEntry-like objects having an
    energy, energy_per_atom and composition.


    * **elements** (*list**[*[*Element*](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element)*]**, **optional*) – Optional list of elements in the phase
    diagram. If set to None, the elements are determined from
    the entries themselves and are sorted alphabetically.
    If specified, element ordering (e.g. for pd coordinates)
    is preserved.


    * **keep_all_spaces** (*bool*) – Boolean control on whether to keep chemical spaces
    that are subspaces of other spaces.


    * **verbose** (*bool*) – Whether to show progress bar during convex hull construction.



#### as_dict()

* **Returns**

    MSONable dictionary representation of PatchedPhaseDiagram.



* **Return type**

    dict[str, Any]



#### _classmethod_ from_dict(dct)

* **Parameters**

    **d** (*dict*) – dictionary representation of PatchedPhaseDiagram.



* **Returns**

    PatchedPhaseDiagram



#### get_all_chempots()
Not Implemented - See PhaseDiagram.


#### get_chempot_range_map()
Not Implemented - See PhaseDiagram.


#### get_chempot_range_stability_phase()
Not Implemented - See PhaseDiagram.


#### get_composition_chempots()
Not Implemented - See PhaseDiagram.


#### get_critical_compositions()
Not Implemented - See PhaseDiagram.


#### get_decomp_and_e_above_hull(entry: PDEntry, allow_negative: bool = False, check_stable: bool = False, on_error: Literal['raise', 'warn', 'ignore'] = 'raise')
Same as method on parent class PhaseDiagram except check_stable defaults to False
for speed. See [https://github.com/materialsproject/pymatgen/issues/2840](https://github.com/materialsproject/pymatgen/issues/2840) for details.


#### get_decomposition(comp: [Composition](pymatgen.core.composition.md#pymatgen.core.composition.Composition))
See PhaseDiagram.


* **Parameters**

    **comp** ([*Composition*](pymatgen.core.composition.md#pymatgen.core.composition.Composition)) – A composition



* **Returns**

    amount} where amount
    is the amount of the fractional composition.



* **Return type**

    Decomposition as a dict of {PDEntry



#### get_element_profile()
Not Implemented - See PhaseDiagram.


#### get_equilibrium_reaction_energy(entry: [Entry](pymatgen.entries.md#pymatgen.entries.Entry))
See PhaseDiagram.

NOTE this is only approximately the same as the what we would get
from PhaseDiagram as we make use of the slsqp approach inside
get_phase_separation_energy().


* **Parameters**

    **entry** (*PDEntry*) – A PDEntry like object



* **Returns**

    Equilibrium reaction energy of entry. Stable entries should have
    equilibrium reaction energy <= 0. The energy is given per atom.



#### get_pd_for_entry(entry: [Entry](pymatgen.entries.md#pymatgen.entries.Entry) | [Composition](pymatgen.core.composition.md#pymatgen.core.composition.Composition))
Get the possible phase diagrams for an entry.


* **Parameters**

    **entry** (*PDEntry** | *[*Composition*](pymatgen.core.composition.md#pymatgen.core.composition.Composition)) – A PDEntry or Composition-like object



* **Returns**

    phase diagram that the entry is part of



* **Return type**

    PhaseDiagram



#### get_transition_chempots()
Not Implemented - See PhaseDiagram.


#### getmu_vertices_stability_phase()
Not Implemented - See PhaseDiagram.


### _class_ pymatgen.analysis.phase_diagram.PhaseDiagram(entries: Sequence[PDEntry] | set[PDEntry], elements: Sequence[[Element](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element)] = (), \*, computed_data: dict[str, Any] | None = None)
Bases: `MSONable`

Simple phase diagram class taking in elements and entries as inputs.
The algorithm is based on the work in the following papers:


1.
    1.
        1. Ong, L. Wang, B. Kang, and G. Ceder, Li-Fe-P-O2 Phase Diagram from

> First Principles Calculations. Chem. Mater., 2008, 20(5), 1798-1807.
> doi:10.1021/cm702327g


2.
    1.
        1. Ong, A. Jain, G. Hautier, B. Kang, G. Ceder, Thermal stabilities

> of delithiated olivine MPO4 (M=Fe, Mn) cathodes investigated using first
> principles calculations. Electrochem. Comm., 2010, 12(3), 427-430.
> doi:10.1016/j.elecom.2010.01.010


#### dim()
The dimensionality of the phase diagram.


* **Type**

    int



#### elements()
Elements in the phase diagram.


#### el_refs()
List of elemental references for the phase diagrams. These are
entries corresponding to the lowest energy element entries for simple
compositional phase diagrams.


#### all_entries()
All entries provided for Phase Diagram construction. Note that this
does not mean that all these entries are actually used in the phase
diagram. For example, this includes the positive formation energy
entries that are filtered out before Phase Diagram construction.


#### qhull_entries()
Actual entries used in convex hull. Excludes all positive formation
energy entries.


#### qhull_data()
Data used in the convex hull operation. This is essentially a matrix of
composition data and energy per atom values created from qhull_entries.


#### facets()
Facets of the phase diagram in the form of  [[1,2,3],[4,5,6]…].
For a ternary, it is the indices (references to qhull_entries and
qhull_data) for the vertices of the phase triangles. Similarly
extended to higher D simplices for higher dimensions.


#### simplices()
The simplices of the phase diagram as a list of np.ndarray, i.e.,
the list of stable compositional coordinates in the phase diagram.


* **Parameters**


    * **entries** (*list**[**PDEntry**]*) – A list of PDEntry-like objects having an
    energy, energy_per_atom and composition.


    * **elements** (*list**[*[*Element*](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element)*]*) – Optional list of elements in the phase
    diagram. If set to None, the elements are determined from
    the entries themselves and are sorted alphabetically.
    If specified, element ordering (e.g. for pd coordinates)
    is preserved.


    * **computed_data** (*dict*) – A dict containing pre-computed data. This allows
    PhaseDiagram object to be reconstituted without performing the
    expensive convex hull computation. The dict is the output from the
    PhaseDiagram._compute() method and is stored in PhaseDiagram.computed_data
    when generated for the first time.



#### _property_ all_entries_hulldata()
Returns:
The actual ndarray used to construct the convex hull.


#### as_dict()

* **Returns**

    MSONable dictionary representation of PhaseDiagram.



#### formation_energy_tol(_ = 1e-1_ )

#### _classmethod_ from_dict(dct: dict[str, Any])

* **Parameters**

    **d** (*dict*) – dictionary representation of PhaseDiagram.



* **Returns**

    PhaseDiagram



#### get_all_chempots(comp)
Get chemical potentials at a given composition.


* **Parameters**

    **comp** ([*Composition*](pymatgen.core.composition.md#pymatgen.core.composition.Composition)) – Composition



* **Returns**

    Chemical potentials.



#### get_chempot_range_map(elements: Sequence[[Element](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element)], referenced: bool = True, joggle: bool = True)
Returns a chemical potential range map for each stable entry.


* **Parameters**


    * **elements** – Sequence of elements to be considered as independent variables.
    E.g., if you want to show the stability ranges
    of all Li-Co-O phases with respect to mu_Li and mu_O, you will supply
    [Element(“Li”), Element(“O”)]


    * **referenced** – If True, gives the results with a reference being the
    energy of the elemental phase. If False, gives absolute values.


    * **joggle** (*bool*) – Whether to joggle the input to avoid precision
    errors.



* **Returns**

    [simplices]}. The list of
    simplices are the sides of the N-1 dim polytope bounding the
    allowable chemical potential range of each entry.



* **Return type**

    Returns a dict of the form {entry



#### get_chempot_range_stability_phase(target_comp, open_elt)
Returns a set of chemical potentials corresponding to the max and min
chemical potential of the open element for a given composition. It is
quite common to have for instance a ternary oxide (e.g., ABO3) for
which you want to know what are the A and B chemical potential leading
to the highest and lowest oxygen chemical potential (reducing and
oxidizing conditions). This is useful for defect computations.


* **Parameters**


    * **target_comp** – A Composition object


    * **open_elt** – Element that you want to constrain to be max or min



* **Returns**

    (mu_min, mu_max)}: Chemical potentials are given in
    “absolute” values (i.e., not referenced to 0)



* **Return type**

    {Element



#### get_composition_chempots(comp)
Get the chemical potentials for all elements at a given composition.


* **Parameters**

    **comp** ([*Composition*](pymatgen.core.composition.md#pymatgen.core.composition.Composition)) – Composition



* **Returns**

    Dictionary of chemical potentials.



#### get_critical_compositions(comp1, comp2)
Get the critical compositions along the tieline between two
compositions. I.e. where the decomposition products change.
The endpoints are also returned.


* **Parameters**


    * **comp1** ([*Composition*](pymatgen.core.composition.md#pymatgen.core.composition.Composition)) – First composition to define the tieline


    * **comp2** ([*Composition*](pymatgen.core.composition.md#pymatgen.core.composition.Composition)) – Second composition to define the tieline



* **Returns**

    list of critical compositions. All are of

        the form x \* comp1 + (1-x) \* comp2




* **Return type**

    [([Composition](pymatgen.core.composition.md#pymatgen.core.composition.Composition))]



#### get_decomp_and_e_above_hull(entry: PDEntry, allow_negative: bool = False, check_stable: bool = True, on_error: Literal['raise', 'warn', 'ignore'] = 'raise')
Provides the decomposition and energy above convex hull for an entry.
Due to caching, can be much faster if entries with the same composition
are processed together.


* **Parameters**


    * **entry** (*PDEntry*) – A PDEntry like object


    * **allow_negative** (*bool*) – Whether to allow negative e_above_hulls. Used to
    calculate equilibrium reaction energies. Defaults to False.


    * **check_stable** (*bool*) – Whether to first check whether an entry is stable.
    In normal circumstances, this is the faster option since checking for
    stable entries is relatively fast. However, if you have a huge proportion
    of unstable entries, then this check can slow things down. You should then
    set this to False.


    * **on_error** (*'raise'** | **'warn'** | **'ignore'*) – What to do if no valid decomposition was
    found. ‘raise’ will throw ValueError. ‘warn’ will print return (None, None).
    ‘ignore’ just returns (None, None). Defaults to ‘raise’.



* **Raises**

    **ValueError** – If no valid decomposition exists in this phase diagram for given entry.



* **Returns**

    (decomp, energy_above_hull). The decomposition is provided

        as a dict of {PDEntry: amount} where amount is the amount of the
        fractional composition. Stable entries should have energy above
        convex hull of 0. The energy is given per atom.




#### get_decomp_and_hull_energy_per_atom(comp: [Composition](pymatgen.core.composition.md#pymatgen.core.composition.Composition))

* **Parameters**

    **comp** ([*Composition*](pymatgen.core.composition.md#pymatgen.core.composition.Composition)) – Input composition.



* **Returns**

    Energy of lowest energy equilibrium at desired composition per atom



#### get_decomp_and_phase_separation_energy(entry: PDEntry, space_limit: int = 200, stable_only: bool = False, tols: Sequence[float] = (1e-08,), maxiter: int = 1000, \*\*kwargs: Any)
Provides the combination of entries in the PhaseDiagram that gives the
lowest formation enthalpy with the same composition as the given entry
excluding entries with the same composition and the energy difference
per atom between the given entry and the energy of the combination found.

For unstable entries that are not polymorphs of stable entries (or completely
novel entries) this is simply the energy above (or below) the convex hull.

For entries with the same composition as one of the stable entries in the
phase diagram setting stable_only to False (Default) allows for entries
not previously on the convex hull to be considered in the combination.
In this case the energy returned is what is referred to as the decomposition
enthalpy in:


1. Bartel, C., Trewartha, A., Wang, Q., Dunn, A., Jain, A., Ceder, G.,

    A critical examination of compound stability predictions from
    machine-learned formation energies, npj Computational Materials 6, 97 (2020)

For stable entries setting stable_only to True returns the same energy
as get_equilibrium_reaction_energy. This function is based on a constrained
optimization rather than recalculation of the convex hull making it
algorithmically cheaper. However, if tol is too loose there is potential
for this algorithm to converge to a different solution.


* **Parameters**


    * **entry** (*PDEntry*) – A PDEntry like object.


    * **space_limit** (*int*) – The maximum number of competing entries to consider
    before calculating a second convex hull to reducing the complexity
    of the optimization.


    * **stable_only** (*bool*) – Only use stable materials as competing entries.


    * **tols** (*list**[**float**]*) – Tolerances for convergence of the SLSQP optimization
    when finding the equilibrium reaction. Tighter tolerances tested first.


    * **maxiter** (*int*) – The maximum number of iterations of the SLSQP optimizer
    when finding the equilibrium reaction.


    * **\*\*kwargs** – Passed to get_decomp_and_e_above_hull.



* **Returns**

    (decomp, energy). The decomposition  is given as a dict of {PDEntry, amount}
    for all entries in the decomp reaction where amount is the amount of the
    fractional composition. The phase separation energy is given per atom.



#### get_decomposition(comp: [Composition](pymatgen.core.composition.md#pymatgen.core.composition.Composition))
Provides the decomposition at a particular composition.


* **Parameters**

    **comp** ([*Composition*](pymatgen.core.composition.md#pymatgen.core.composition.Composition)) – A composition



* **Returns**

    amount} where amount
    is the amount of the fractional composition.



* **Return type**

    Decomposition as a dict of {PDEntry



#### get_e_above_hull(entry: PDEntry, \*\*kwargs: Any)
Provides the energy above convex hull for an entry.


* **Parameters**


    * **entry** (*PDEntry*) – A PDEntry like object.


    * **\*\*kwargs** – Passed to get_decomp_and_e_above_hull().



* **Returns**

    Energy above convex hull of entry. Stable entries should have

        energy above hull of 0. The energy is given per atom.




* **Return type**

    float | None



#### get_element_profile(element, comp, comp_tol=1e-05)
Provides the element evolution data for a composition. For example, can be used
to analyze Li conversion voltages by varying mu_Li and looking at the phases
formed. Also can be used to analyze O2 evolution by varying mu_O2.


* **Parameters**


    * **element** – An element. Must be in the phase diagram.


    * **comp** – A Composition


    * **comp_tol** – The tolerance to use when calculating decompositions.
    Phases with amounts less than this tolerance are excluded.
    Defaults to 1e-5.



* **Returns**

    [ {‘chempot’: -10.487582010000001, ‘evolution’: -2.0,
    ‘reaction’: Reaction Object], …]



* **Return type**

    Evolution data as a list of dictionaries of the following format



#### get_equilibrium_reaction_energy(entry: PDEntry)
Provides the reaction energy of a stable entry from the neighboring
equilibrium stable entries (also known as the inverse distance to
hull).


* **Parameters**

    **entry** (*PDEntry*) – A PDEntry like object



* **Returns**

    Equilibrium reaction energy of entry. Stable entries should have

        equilibrium reaction energy <= 0. The energy is given per atom.




* **Return type**

    float | None



#### get_form_energy(entry: PDEntry)
Returns the formation energy for an entry (NOT normalized) from the
elemental references.


* **Parameters**

    **entry** (*PDEntry*) – A PDEntry-like object.



* **Returns**

    Formation energy from the elemental references.



* **Return type**

    float



#### get_form_energy_per_atom(entry: PDEntry)
Returns the formation energy per atom for an entry from the
elemental references.


* **Parameters**

    **entry** (*PDEntry*) – An PDEntry-like object



* **Returns**

    Formation energy **per atom** from the elemental references.



#### get_hull_energy(comp: [Composition](pymatgen.core.composition.md#pymatgen.core.composition.Composition))

* **Parameters**

    **comp** ([*Composition*](pymatgen.core.composition.md#pymatgen.core.composition.Composition)) – Input composition.



* **Returns**

    Energy of lowest energy equilibrium at desired composition. Not

        normalized by atoms, i.e. E(Li4O2) = 2 \* E(Li2O)




#### get_hull_energy_per_atom(comp: [Composition](pymatgen.core.composition.md#pymatgen.core.composition.Composition), \*\*kwargs)

* **Parameters**

    **comp** ([*Composition*](pymatgen.core.composition.md#pymatgen.core.composition.Composition)) – Input composition.



* **Returns**

    Energy of lowest energy equilibrium at desired composition.



#### get_phase_separation_energy(entry, \*\*kwargs)
Provides the energy to the convex hull for the given entry. For stable entries
already in the phase diagram the algorithm provides the phase separation energy
which is referred to as the decomposition enthalpy in:


1. Bartel, C., Trewartha, A., Wang, Q., Dunn, A., Jain, A., Ceder, G.,

    A critical examination of compound stability predictions from
    machine-learned formation energies, npj Computational Materials 6, 97 (2020)


* **Parameters**


    * **entry** (*PDEntry*) – A PDEntry like object


    * **\*\*kwargs** – Keyword args passed to get_decomp_and_decomp_energy
    space_limit (int): The maximum number of competing entries to consider.
    stable_only (bool): Only use stable materials as competing entries
    tol (float): The tolerance for convergence of the SLSQP optimization

    > when finding the equilibrium reaction.

    maxiter (int): The maximum number of iterations of the SLSQP optimizer

        when finding the equilibrium reaction.




* **Returns**

    phase separation energy per atom of entry. Stable entries should have
    energies <= 0, Stable elemental entries should have energies = 0 and
    unstable entries should have energies > 0. Entries that have the same
    composition as a stable energy may have positive or negative phase
    separation energies depending on their own energy.



#### get_plot(show_unstable: float = 0.2, backend: Literal['plotly', 'matplotlib'] = 'plotly', ternary_style: Literal['2d', '3d'] = '2d', label_stable: bool = True, label_unstable: bool = True, ordering: Sequence[str] | None = None, energy_colormap=None, process_attributes: bool = False, plt=None, label_uncertainties: bool = False, fill: bool = True, \*\*plotkwargs)
Convenient wrapper for PDPlotter. Initializes a PDPlotter object and calls
get_plot() with provided combined arguments.

Plotting is only supported for phase diagrams with <=4 elements (unary,
binary, ternary, or quaternary systems).


* **Parameters**


    * **show_unstable** (*float*) – Whether unstable (above the hull) phases will be
    plotted. If a number > 0 is entered, all phases with
    e_hull < show_unstable (eV/atom) will be shown.


    * **backend** (*"plotly"** | **"matplotlib"*) – Python package to use for plotting.
    Defaults to “plotly”.


    * **ternary_style** (*"2d"** | **"3d"*) – Ternary phase diagrams are typically plotted in
    two-dimensions (2d), but can be plotted in three dimensions (3d) to visualize
    the depth of the hull. This argument only applies when backend=”plotly”.
    Defaults to “2d”.


    * **label_stable** – Whether to label stable compounds.


    * **label_unstable** – Whether to label unstable compounds.


    * **ordering** – Ordering of vertices (matplotlib backend only).


    * **energy_colormap** – Colormap for coloring energy (matplotlib backend only).


    * **process_attributes** – Whether to process the attributes (matplotlib
    backend only).


    * **plt** – Existing plt object if plotting multiple phase diagrams (
    matplotlib backend only).


    * **label_uncertainties** – Whether to add error bars to the hull (plotly
    backend only). For binaries, this also shades the hull with the
    uncertainty window.


    * **fill** – Whether to shade the hull. For ternary_2d and quaternary plots, this
    colors facets arbitrarily for visual clarity. For ternary_3d plots, this
    shades the hull by formation energy (plotly backend only).


    * **\*\*plotkwargs** (*dict*) – Keyword args passed to matplotlib.pyplot.plot (only
    applies when backend=”matplotlib”). Can be used to customize markers
    etc. If not set, the default is:

    > {

    >     “markerfacecolor”: “#4daf4a”,
    >     “markersize”: 10,
    >     “linewidth”: 3

    > }




#### get_reference_energy_per_atom(comp: [Composition](pymatgen.core.composition.md#pymatgen.core.composition.Composition))

* **Parameters**

    **comp** ([*Composition*](pymatgen.core.composition.md#pymatgen.core.composition.Composition)) – Input composition.



* **Returns**

    Reference energy of the terminal species at a given composition.



#### get_transition_chempots(element)
Get the critical chemical potentials for an element in the Phase
Diagram.


* **Parameters**

    **element** – An element. Has to be in the PD in the first place.



* **Returns**

    A sorted sequence of critical chemical potentials, from less
    negative to more negative.



#### getmu_vertices_stability_phase(target_comp, dep_elt, tol_en=0.01)
Returns a set of chemical potentials corresponding to the vertices of
the simplex in the chemical potential phase diagram.
The simplex is built using all elements in the target_composition
except dep_elt.
The chemical potential of dep_elt is computed from the target
composition energy.
This method is useful to get the limiting conditions for
defects computations for instance.


* **Parameters**


    * **target_comp** – A Composition object


    * **dep_elt** – the element for which the chemical potential is computed
    from the energy of the stable phase at the target composition


    * **tol_en** – a tolerance on the energy to set



* **Returns**

    mu}]: An array of conditions on simplex vertices for
    which each element has a chemical potential set to a given
    value. “absolute” values (i.e., not referenced to element energies)



* **Return type**

    [{Element



#### numerical_tol(_ = 1e-0_ )

#### pd_coords(comp: [Composition](pymatgen.core.composition.md#pymatgen.core.composition.Composition))
The phase diagram is generated in a reduced dimensional space
(n_elements - 1). This function returns the coordinates in that space.
These coordinates are compatible with the stored simplex objects.


* **Parameters**

    **comp** ([*Composition*](pymatgen.core.composition.md#pymatgen.core.composition.Composition)) – A composition



* **Returns**

    The coordinates for a given composition in the PhaseDiagram’s basis



#### _property_ stable_entries(_: set[[pymatgen.entries.Entry](pymatgen.entries.md#pymatgen.entries.Entry)_ )
Returns:
set[Entry]: of stable entries in the phase diagram.


#### _property_ unstable_entries(_: set[[pymatgen.entries.Entry](pymatgen.entries.md#pymatgen.entries.Entry)_ )
Returns:
set[Entry]: unstable entries in the phase diagram. Includes positive formation energy entries.


### _exception_ pymatgen.analysis.phase_diagram.PhaseDiagramError()
Bases: `Exception`

An exception class for Phase Diagram generation.


### _class_ pymatgen.analysis.phase_diagram.ReactionDiagram(entry1, entry2, all_entries, tol: float = 0.0001, float_fmt='%.4f')
Bases: `object`

Analyzes the possible reactions between a pair of compounds, e.g.,
an electrolyte and an electrode.


* **Parameters**


    * **entry1** ([*ComputedEntry*](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedEntry)) – Entry for 1st component. Note that
    corrections, if any, must already be pre-applied. This is to
    give flexibility for different kinds of corrections, e.g.,
    if a particular entry is fitted to an experimental data (such
    as EC molecule).


    * **entry2** ([*ComputedEntry*](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedEntry)) – Entry for 2nd component. Note that
    corrections must already be pre-applied. This is to
    give flexibility for different kinds of corrections, e.g.,
    if a particular entry is fitted to an experimental data (such
    as EC molecule).


    * **all_entries** (*[*[*ComputedEntry*](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedEntry)*]*) – All other entries to be
    considered in the analysis. Note that corrections, if any,
    must already be pre-applied.


    * **tol** (*float*) – Tolerance to be used to determine validity of reaction.


    * **float_fmt** (*str*) – Formatting string to be applied to all floats.
    Determines number of decimal places in reaction string.



#### get_compound_pd()
Get the CompoundPhaseDiagram object, which can then be used for
plotting.


* **Returns**

    CompoundPhaseDiagram



### _class_ pymatgen.analysis.phase_diagram.TransformedPDEntry(entry, sp_mapping, name=None)
Bases: `PDEntry`

This class represents a TransformedPDEntry, which allows for a PDEntry to be
transformed to a different composition coordinate space. It is used in the
construction of phase diagrams that do not have elements as the terminal
compositions.


* **Parameters**


    * **entry** (*PDEntry*) – Original entry to be transformed.


    * **(****{Composition** (*sp_mapping*) – DummySpecies}): dictionary mapping Terminal Compositions to Dummy Species.



#### amount_tol(_ = 1e-0_ )

#### as_dict()

* **Returns**

    MSONable dictionary representation of TransformedPDEntry.



#### _property_ composition(_: [Composition](pymatgen.core.composition.md#pymatgen.core.composition.Composition_ )
The composition in the dummy species space.


* **Returns**

    Composition



#### _classmethod_ from_dict(d)

* **Parameters**

    **d** (*dict*) – dictionary representation of TransformedPDEntry.



* **Returns**

    TransformedPDEntry



### _exception_ pymatgen.analysis.phase_diagram.TransformedPDEntryError()
Bases: `Exception`

An exception class for TransformedPDEntry.


### pymatgen.analysis.phase_diagram.get_facets(qhull_data: ArrayLike, joggle: bool = False)
Get the simplex facets for the Convex hull.


* **Parameters**


    * **qhull_data** (*np.ndarray*) – The data from which to construct the convex
    hull as a Nxd array (N being number of data points and d being the
    dimension)


    * **joggle** (*bool*) – Whether to joggle the input to avoid precision
    errors.



* **Returns**

    List of simplices of the Convex Hull.



### pymatgen.analysis.phase_diagram.order_phase_diagram(lines, stable_entries, unstable_entries, ordering)
Orders the entries (their coordinates) in a phase diagram plot according
to the user specified ordering.
Ordering should be given as [‘Up’, ‘Left’, ‘Right’], where Up,
Left and Right are the names of the entries in the upper, left and right
corners of the triangle respectively.


* **Parameters**


    * **lines** – list of list of coordinates for lines in the PD.


    * **stable_entries** – {coordinate : entry} for each stable node in the
    phase diagram. (Each coordinate can only have one stable phase)


    * **unstable_entries** – {entry: coordinates} for all unstable nodes in the
    phase diagram.


    * **ordering** – Ordering of the phase diagram, given as a list [‘Up’,
    ‘Left’,’Right’]



* **Returns**


    * newlines is a list of list of coordinates for lines in the PD.


    * newstable_entries is a {coordinate : entry} for each stable node

    in the phase diagram. (Each coordinate can only have one
    stable phase)
    - newunstable_entries is a {entry: coordinates} for all unstable
    nodes in the phase diagram.




* **Return type**

    (newlines, newstable_entries, newunstable_entries)



### pymatgen.analysis.phase_diagram.tet_coord(coord)
Convert a 3D coordinate into a tetrahedron based coordinate system for a
prettier phase diagram.


* **Parameters**

    **coord** – coordinate used in the convex hull computation.



* **Returns**

    coordinates in a tetrahedron-based coordinate system.



### pymatgen.analysis.phase_diagram.triangular_coord(coord)
Convert a 2D coordinate into a triangle-based coordinate system for a
prettier phase diagram.


* **Parameters**

    **coord** – coordinate used in the convex hull computation.



* **Returns**

    coordinates in a triangular-based coordinate system.



### pymatgen.analysis.phase_diagram.uniquelines(q)
Given all the facets, convert it into a set of unique lines. Specifically
used for converting convex hull facets into line pairs of coordinates.


* **Parameters**

    **q** – A 2-dim sequence, where each row represents a facet. E.g.,
    [[1,2,3],[3,6,7],…]



* **Returns**

    A set of tuple of lines. E.g., ((1,2), (1,3), (2,3), ….)



* **Return type**

    setoflines