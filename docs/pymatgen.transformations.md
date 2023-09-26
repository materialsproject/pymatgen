---
layout: default
title: pymatgen.transformations.md
nav_exclude: true
---

1. TOC
{:toc}

# pymatgen.transformations package

The transformations package defines various transformations that can be applied
on structures, i.e., converting one structure to another.


## pymatgen.transformations.advanced_transformations module

This module implements more advanced transformations.


### _class_ AddAdsorbateTransformation(adsorbate, selective_dynamics=False, height=0.9, mi_vec=None, repeat=None, min_lw=5.0, translate=True, reorient=True, find_args=None)
Bases: `AbstractTransformation`

Create absorbate structures.

Use AdsorbateSiteFinder to add an absorbate to a slab.


* **Parameters**


    * **adsorbate** ([*Molecule*](pymatgen.core.md#pymatgen.core.structure.Molecule)) – molecule to add as adsorbate


    * **selective_dynamics** (*bool*) – flag for whether to assign
    non-surface sites as fixed for selective dynamics


    * **height** (*float*) – height criteria for selection of surface sites


    * **mi_vec** – vector corresponding to the vector
    concurrent with the miller index, this enables use with
    slabs that have been reoriented, but the miller vector
    must be supplied manually


    * **repeat** (*3-tuple** or **list*) – repeat argument for supercell generation


    * **min_lw** (*float*) – minimum length and width of the slab, only used
    if repeat is None


    * **translate** (*bool*) – flag on whether to translate the molecule so
    that its CoM is at the origin prior to adding it to the surface


    * **reorient** (*bool*) – flag on whether or not to reorient adsorbate
    along the miller index


    * **find_args** (*dict*) – dictionary of arguments to be passed to the
    call to self.find_adsorption_sites, e.g. {“distance”:2.0}



#### _abc_impl(_ = <_abc._abc_data object_ )

#### apply_transformation(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), return_ranked_list: bool | int = False)

* **Parameters**


    * **structure** – Must be a Slab structure


    * **return_ranked_list** (*bool** | **int**, **optional*) – If return_ranked_list is int, that number of structures.

    is returned. If False, only the single lowest energy structure is returned. Defaults to False.




* **Returns**

    with adsorbate



* **Return type**

    [Slab](pymatgen.core.md#pymatgen.core.surface.Slab)



#### _property_ inverse()
None.


* **Type**

    Returns



#### _property_ is_one_to_many(_: boo_ )
True.


* **Type**

    Returns



### _class_ ChargeBalanceTransformation(charge_balance_sp)
Bases: `AbstractTransformation`

This is a transformation that disorders a structure to make it charge
balanced, given an oxidation state-decorated structure.


* **Parameters**

    **charge_balance_sp** – specie to add or remove. Currently only removal
    is supported.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### apply_transformation(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))
Applies the transformation.


* **Parameters**

    **structure** – Input Structure



* **Returns**

    Charge balanced structure.



#### _property_ inverse()
None.


* **Type**

    Returns



#### _property_ is_one_to_many(_: boo_ )
False.


* **Type**

    Returns



### _class_ CubicSupercellTransformation(min_atoms: int | None = None, max_atoms: int | None = None, min_length: float = 15.0, force_diagonal: bool = False, force_90_degrees: bool = False, angle_tolerance: float = 0.001)
Bases: `AbstractTransformation`

A transformation that aims to generate a nearly cubic supercell structure
from a structure.

The algorithm solves for a transformation matrix that makes the supercell
cubic. The matrix must have integer entries, so entries are rounded (in such
a way that forces the matrix to be non-singular). From the supercell
resulting from this transformation matrix, vector projections are used to
determine the side length of the largest cube that can fit inside the
supercell. The algorithm will iteratively increase the size of the supercell
until the largest inscribed cube’s side length is at least ‘min_length’
and the number of atoms in the supercell falls in the range
`min_atoms < n < max_atoms`.


* **Parameters**


    * **max_atoms** – Maximum number of atoms allowed in the supercell.


    * **min_atoms** – Minimum number of atoms allowed in the supercell.


    * **min_length** – Minimum length of the smallest supercell lattice vector.


    * **force_diagonal** – If True, return a transformation with a diagonal
    transformation matrix.


    * **force_90_degrees** – If True, return a transformation for a supercell
    with 90 degree angles (if possible). To avoid long run times,
    please use max_atoms


    * **angle_tolerance** – tolerance to determine the 90 degree angles.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### apply_transformation(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))
The algorithm solves for a transformation matrix that makes the
supercell cubic. The matrix must have integer entries, so entries are
rounded (in such a way that forces the matrix to be non-singular). From
the supercell resulting from this transformation matrix, vector
projections are used to determine the side length of the largest cube
that can fit inside the supercell. The algorithm will iteratively
increase the size of the supercell until the largest inscribed cube’s
side length is at least ‘num_nn_dists’ times the nearest neighbor
distance and the number of atoms in the supercell falls in the range
defined by min_atoms and max_atoms.


* **Returns**

    Transformed supercell.



* **Return type**

    supercell



#### _property_ inverse()
Returns None.


#### _property_ is_one_to_many(_: boo_ )
Returns False.


### _class_ DisorderOrderedTransformation(max_sites_to_merge=2)
Bases: `AbstractTransformation`

Not to be confused with OrderDisorderedTransformation,
this transformation attempts to obtain a
*disordered* structure from an input ordered structure.
This may or may not be physically plausible, further
inspection of the returned structures is advised.
The main purpose for this transformation is for structure
matching to crystal prototypes for structures that have
been derived from a parent prototype structure by
substitutions or alloying additions.


* **Parameters**

    **max_sites_to_merge** – only merge this number of sites together.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### _static_ _get_disorder_mappings(composition, partitions)
Private method to obtain the mapping to create
a disordered structure from a given partition.


#### _static_ _partition_species(composition, max_components=2)
Private method to split a list of species into
various partitions.


#### apply_transformation(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), return_ranked_list: bool | int = False)

* **Parameters**


    * **structure** – ordered structure


    * **return_ranked_list** (*bool** | **int**, **optional*) – If return_ranked_list is int, that number of structures.

    is returned. If False, only the single lowest energy structure is returned. Defaults to False.




* **Returns**

    Transformed disordered structure(s)



#### _property_ inverse()
None.


* **Type**

    Returns



#### _property_ is_one_to_many(_: boo_ )
True.


* **Type**

    Returns



### _class_ DopingTransformation(dopant, ionic_radius_tol=inf, min_length=10, alio_tol=0, codopant=False, max_structures_per_enum=100, allowed_doping_species=None, \*\*kwargs)
Bases: `AbstractTransformation`

A transformation that performs doping of a structure.


* **Parameters**


    * **dopant** (*Species-like*) – E.g., Al3+. Must have oxidation state.


    * **ionic_radius_tol** (*float*) – E.g., Fractional allowable ionic radii
    mismatch for dopant to fit into a site. Default of inf means
    that any dopant with the right oxidation state is allowed.


    * **min_length** (*float*) – Min. lattice parameter between periodic
    images of dopant. Defaults to 10A for now.


    * **alio_tol** (*int*) – If this is not 0, attempt will be made to dope
    sites with oxidation_states +- alio_tol of the dopant. E.g.,
    1 means that the ions like Ca2+ and Ti4+ are considered as
    potential doping sites for Al3+.


    * **codopant** (*bool*) – If True, doping will be carried out with a
    codopant to maintain charge neutrality. Otherwise, vacancies
    will be used.


    * **max_structures_per_enum** (*float*) – Maximum number of structures to
    return per enumeration. Note that there can be more than one
    candidate doping site, and each site enumeration will return at
    max max_structures_per_enum structures. Defaults to 100.


    * **allowed_doping_species** (*list*) – Species that are allowed to be
    doping sites. This is an inclusionary list. If specified,
    any sites which are not


    * **\*\*kwargs** – Same keyword args as EnumerateStructureTransformation,
    i.e., min_cell_size, etc.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### apply_transformation(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), return_ranked_list: bool | int = False)

* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Input structure to dope


    * **return_ranked_list** (*bool** | **int**, **optional*) – If return_ranked_list is int, that number of structures.

    is returned. If False, only the single lowest energy structure is returned. Defaults to False.




* **Returns**

    Structure, “energy”: float}]



* **Return type**

    [{“structure”



#### _property_ inverse()
None.


* **Type**

    Returns



#### _property_ is_one_to_many(_: boo_ )
True.


* **Type**

    Returns



### _class_ EnumerateStructureTransformation(min_cell_size: int = 1, max_cell_size: int = 1, symm_prec: float = 0.1, refine_structure: bool = False, enum_precision_parameter: float = 0.001, check_ordered_symmetry: bool = True, max_disordered_sites: int | None = None, sort_criteria: str | Callable = 'ewald', timeout: float | None = None, n_jobs: int = -1)
Bases: `AbstractTransformation`

Order a disordered structure using enumlib. For complete orderings, this
generally produces fewer structures that the OrderDisorderedStructure
transformation, and at a much faster speed.


* **Parameters**


    * **min_cell_size** – The minimum cell size wanted. Must be an int. Defaults to 1.


    * **max_cell_size** – The maximum cell size wanted. Must be an int. Defaults to 1.


    * **symm_prec** – Tolerance to use for symmetry.


    * **refine_structure** – This parameter has the same meaning as in enumlib_caller.
    If you are starting from a structure that has been relaxed via
    some electronic structure code, it is usually much better to
    start with symmetry determination and then obtain a refined
    structure. The refined structure have cell parameters and
    atomic positions shifted to the expected symmetry positions,
    which makes it much less sensitive precision issues in enumlib.
    If you are already starting from an experimental cif, refinement
    should have already been done and it is not necessary. Defaults
    to False.


    * **enum_precision_parameter** (*float*) – Finite precision parameter for
    enumlib. Default of 0.001 is usually ok, but you might need to
    tweak it for certain cells.


    * **check_ordered_symmetry** (*bool*) – Whether to check the symmetry of
    the ordered sites. If the symmetry of the ordered sites is
    lower, the lowest symmetry ordered sites is included in the
    enumeration. This is important if the ordered sites break
    symmetry in a way that is important getting possible
    structures. But sometimes including ordered sites
    slows down enumeration to the point that it cannot be
    completed. Switch to False in those cases. Defaults to True.


    * **max_disordered_sites** (*int*) – An alternate parameter to max_cell size. Will sequentially try
    larger and larger cell sizes until (i) getting a result or (ii)
    the number of disordered sites in the cell exceeds
    max_disordered_sites. Must set max_cell_size to None when using
    this parameter.


    * **sort_criteria** (*str** or **callable*) – Sort by Ewald energy (“ewald”, must have oxidation states and slow) or
    M3GNet relaxed energy (“m3gnet_relax”, which is the most accurate but most expensive and provides
    pre-relaxed structures - needs m3gnet package installed) or by M3GNet static energy (“m3gnet_static”)
    or by number of sites (“nsites”, the fastest, the default). The expense of m3gnet_relax or m3gnet_static
    can be worth it if it significantly reduces the number of structures to be considered. m3gnet_relax
    speeds up the subsequent DFT calculations. Alternatively, a callable can be supplied that returns a
    (Structure, energy) tuple.


    * **timeout** (*float*) – timeout in minutes to pass to EnumlibAdaptor.


    * **n_jobs** (*int*) – Number of parallel jobs used to compute energy criteria. This is used only when the Ewald
    or m3gnet or callable sort_criteria is used. Default is -1, which uses all available CPUs.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### apply_transformation(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), return_ranked_list: bool | int = False)
Returns either a single ordered structure or a sequence of all ordered
structures.


* **Parameters**


    * **structure** – Structure to order.


    * **return_ranked_list** (*bool** | **int**, **optional*) – If return_ranked_list is int, that number of structures

    is returned. If False, only the single lowest energy structure is returned. Defaults to False.




* **Returns**

    Depending on returned_ranked list, either a transformed structure
    or a list of dictionaries, where each dictionary is of the form
    {“structure” = …. , “other_arguments”}

    The list of ordered structures is ranked by Ewald energy / atom, if
    the input structure is an oxidation state decorated structure.
    Otherwise, it is ranked by number of sites, with smallest number of
    sites first.




#### _property_ inverse()
None.


* **Type**

    Returns



#### _property_ is_one_to_many(_: boo_ )
True.


* **Type**

    Returns



### _class_ GrainBoundaryTransformation(rotation_axis, rotation_angle, expand_times=4, vacuum_thickness=0.0, ab_shift=None, normal=False, ratio=True, plane=None, max_search=20, tol_coi=1e-08, rm_ratio=0.7, quick_gen=False)
Bases: `AbstractTransformation`

A transformation that creates a gb from a bulk structure.


* **Parameters**


    * **rotation_axis** (*list*) – Rotation axis of GB in the form of a list of integer
    e.g.: [1, 1, 0]


    * **rotation_angle** (*float**, **in unit** of **degree*) – rotation angle used to generate GB.
    Make sure the angle is accurate enough. You can use the enum\* functions
    in this class to extract the accurate angle.
    e.g.: The rotation angle of sigma 3 twist GB with the rotation axis
    [1, 1, 1] and GB plane (1, 1, 1) can be 60.000000000 degree.
    If you do not know the rotation angle, but know the sigma value, we have
    provide the function get_rotation_angle_from_sigma which is able to return
    all the rotation angles of sigma value you provided.


    * **expand_times** (*int*) – The multiple times used to expand one unit grain to larger grain.
    This is used to tune the grain length of GB to warrant that the two GBs in one
    cell do not interact with each other. Default set to 4.


    * **vacuum_thickness** (*float*) – The thickness of vacuum that you want to insert between
    two grains of the GB. Default to 0.


    * **ab_shift** (*list** of **float**, **in unit** of **a**, **b vectors** of **Gb*) – in plane shift of two grains


    * **normal** (*logic*) – determine if need to require the c axis of top grain (first transformation matrix)
    perpendicular to the surface or not.
    default to false.


    * **ratio** (*list** of **integers*) – lattice axial ratio.
    If True, will try to determine automatically from structure.
    For cubic system, ratio is not needed and can be set to None.
    For tetragonal system, ratio = [mu, mv], list of two integers,
    that is, mu/mv = c2/a2. If it is irrational, set it to None.
    For orthorhombic system, ratio = [mu, lam, mv], list of three integers,

    > that is, mu:lam:mv = c2:b2:a2. If irrational for one axis, set it to None.

    e.g. mu:lam:mv = c2,None,a2, means b2 is irrational.
    For rhombohedral system, ratio = [mu, mv], list of two integers,
    that is, mu/mv is the ratio of (1+2\*cos(alpha))/cos(alpha).
    If irrational, set it to None.
    For hexagonal system, ratio = [mu, mv], list of two integers,
    that is, mu/mv = c2/a2. If it is irrational, set it to none.



    * **plane** (*list*) – Grain boundary plane in the form of a list of integers
    e.g.: [1, 2, 3]. If none, we set it as twist GB. The plane will be perpendicular
    to the rotation axis.


    * **max_search** (*int*) – max search for the GB lattice vectors that give the smallest GB
    lattice. If normal is true, also max search the GB c vector that perpendicular
    to the plane. For complex GB, if you want to speed up, you can reduce this value.
    But too small of this value may lead to error.


    * **tol_coi** (*float*) – tolerance to find the coincidence sites. When making approximations to
    the ratio needed to generate the GB, you probably need to increase this tolerance to
    obtain the correct number of coincidence sites. To check the number of coincidence
    sites are correct or not, you can compare the generated Gb object’s sigma with enum\*
    sigma values (what user expected by input).


    * **rm_ratio** (*float*) – the criteria to remove the atoms which are too close with each other.
    rm_ratio \* bond_length of bulk system is the criteria of bond length, below which the atom
    will be removed. Default to 0.7.


    * **quick_gen** (*bool*) – whether to quickly generate a supercell, if set to true, no need to
    find the smallest cell.



* **Returns**

    Grain boundary structure (gb (Structure) object).



#### _abc_impl(_ = <_abc._abc_data object_ )

#### apply_transformation(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))
Applies the transformation.


* **Parameters**


    * **structure** – Input Structure


    * **return_ranked_list** (*bool** | **int**, **optional*) – If return_ranked_list is int, that number of structures
    is returned. If False, only the single lowest energy structure is returned. Defaults to False.



* **Returns**

    Grain boundary Structures.



#### _property_ inverse()
None.


* **Type**

    Returns



#### _property_ is_one_to_many(_: boo_ )
False.


* **Type**

    Returns



### _class_ MagOrderParameterConstraint(order_parameter, species_constraints=None, site_constraint_name=None, site_constraints=None)
Bases: `MSONable`

This class can be used to supply MagOrderingTransformation
to just a specific subset of species or sites that satisfy the
provided constraints. This can be useful for setting an order
parameters for, for example, ferrimagnetic structures which
might order on certain motifs, with the global order parameter
dependent on how many sites satisfy that motif.


* **Parameters**


    * **(****float****)** (*order_parameter*) – any number from 0.0 to 1.0,
    typically 0.5 (antiferromagnetic) or 1.0 (ferromagnetic)


    * **(****list****)** (*site_constraints*) – str or list of strings
    of Species symbols that the constraint should apply to


    * **(****str****)** (*site_constraint_name*) – name of the site property
    that the constraint should apply to, e.g. “coordination_no”


    * **(****list****)** – list of values of the site
    property that the constraints should apply to



#### satisfies_constraint(site)
Checks if a periodic site satisfies the constraint.


### _class_ MagOrderingTransformation(mag_species_spin, order_parameter=0.5, energy_model=None, \*\*kwargs)
Bases: `AbstractTransformation`

This transformation takes a structure and returns a list of collinear
magnetic orderings. For disordered structures, make an ordered
approximation first.


* **Parameters**


    * **mag_species_spin** – A mapping of elements/species to their
    spin magnitudes, e.g. {“Fe3+”: 5, “Mn3+”: 4}


    * **list****)** (*order_parameter** (**float or*) – if float, a specifies a
    global order parameter and can take values from 0.0 to 1.0
    (e.g. 0.5 for antiferromagnetic or 1.0 for ferromagnetic), if
    list has to be a list of
    pymatgen.transformations.advanced_transformations.MagOrderParameterConstraint
    to specify more complicated orderings, see documentation for
    MagOrderParameterConstraint more details on usage


    * **energy_model** – Energy model to rank the returned structures,
    see :mod: pymatgen.analysis.energy_models for more information (note
    that this is not necessarily a physical energy). By default, returned
    structures use SymmetryModel() which ranks structures from most
    symmetric to least.


    * **kwargs** – Additional kwargs that are passed to


EnumerateStructureTransformation such as min_cell_size etc.


#### _abc_impl(_ = <_abc._abc_data object_ )

#### _static_ _add_dummy_species(structure, order_parameters)

* **Parameters**


    * **structure** – ordered Structure


    * **order_parameters** – list of MagOrderParameterConstraints



* **Returns**

    A structure decorated with disordered
    DummySpecies on which to perform the enumeration.
    Note that the DummySpecies are super-imposed on
    to the original sites, to make it easier to
    retrieve the original site after enumeration is
    performed (this approach is preferred over a simple
    mapping since multiple species may have the same
    DummySpecies, depending on the constraints specified).
    This approach can also preserve site properties even after
    enumeration.



#### _add_spin_magnitudes(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))
Replaces Spin.up/Spin.down with spin magnitudes specified by mag_species_spin.


* **Parameters**

    **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Structure to modify.



* **Returns**

    Structure with spin magnitudes added.



* **Return type**

    [Structure](pymatgen.core.md#pymatgen.core.structure.Structure)



#### _static_ _remove_dummy_species(structure)
Structure with dummy species removed, but their corresponding spin properties
merged with the original sites. Used after performing enumeration.


#### apply_transformation(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), return_ranked_list: bool | int = False)
Apply MagOrderTransformation to an input structure.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Any ordered structure.


    * **return_ranked_list** (*bool** | **int**, **optional*) – If return_ranked_list is int, that number of structures
    is returned. If False, only the single lowest energy structure is returned. Defaults to False.



* **Raises**

    **ValueError** – On disordered structures.



* **Returns**

    Structure(s) after MagOrderTransformation.



* **Return type**

    [Structure](pymatgen.core.md#pymatgen.core.structure.Structure) | list[[Structure](pymatgen.core.md#pymatgen.core.structure.Structure)]



#### _static_ determine_min_cell(disordered_structure)
Determine the smallest supercell that is able to enumerate
the provided structure with the given order parameter.


#### _property_ inverse(_: Non_ )
None.


* **Type**

    Returns



#### _property_ is_one_to_many(_: boo_ )
True.


* **Type**

    Returns



### _class_ MonteCarloRattleTransformation(rattle_std: float, min_distance: float, seed: int | None = None, \*\*kwargs)
Bases: `AbstractTransformation`

Uses a Monte Carlo rattle procedure to randomly perturb the sites in a
structure.

This class requires the hiPhive package to be installed.

Rattling atom i is carried out as a Monte Carlo move that is accepted with
a probability determined from the minimum interatomic distance
$d_{ij}$. If $\\min(d_{ij})$ is smaller than $d_{min}$
the move is only accepted with a low probability.

This process is repeated for each atom a number of times meaning
the magnitude of the final displacements is not *directly*
connected to rattle_std.


* **Parameters**


    * **rattle_std** – Rattle amplitude (standard deviation in normal
    distribution). Note: this value is not *directly* connected to the
    final average displacement for the structures


    * **min_distance** – Interatomic distance used for computing the probability
    for each rattle move.


    * **seed** – Seed for setting up NumPy random state from which random numbers
    are generated. If `None`, a random seed will be generated
    (default). This option allows the output of this transformation
    to be deterministic.


    * **\*\*kwargs** – Additional keyword arguments to be passed to the hiPhive
    mc_rattle function.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### apply_transformation(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))
Apply the transformation.


* **Parameters**

    **structure** – Input Structure



* **Returns**

    Structure with sites perturbed.



#### _property_ inverse()
None.


* **Type**

    Returns



#### _property_ is_one_to_many(_: boo_ )
False.


* **Type**

    Returns



### _class_ MultipleSubstitutionTransformation(sp_to_replace, r_fraction, substitution_dict, charge_balance_species=None, order=True)
Bases: `object`

Performs multiple substitutions on a structure. For example, can do a
fractional replacement of Ge in LiGePS with a list of species, creating one
structure for each substitution. Ordering is done using a dummy element so
only one ordering must be done per substitution oxidation state. Charge
balancing of the structure is optionally performed.

**NOTE**: There are no checks to make sure that removal fractions are possible and rounding
may occur. Currently charge balancing only works for removal of species.

Performs multiple fractional substitutions on a transmuter.


* **Parameters**


    * **sp_to_replace** – species to be replaced


    * **r_fraction** – fraction of that specie to replace


    * **substitution_dict** – dictionary of the format
    {2: [“Mg”, “Ti”, “V”, “As”, “Cr”, “Ta”, “N”, “Nb”],
    3: [“Ru”, “Fe”, “Co”, “Ce”, “As”, “Cr”, “Ta”, “N”, “Nb”],
    4: [“Ru”, “V”, “Cr”, “Ta”, “N”, “Nb”],
    5: [“Ru”, “W”, “Mn”]
    }
    The number is the charge used for each of the list of elements
    (an element can be present in multiple lists)


    * **charge_balance_species** – If specified, will balance the charge on
    the structure using that specie.


    * **order** – Whether to order the structures.



#### apply_transformation(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), return_ranked_list: bool | int = False)
Applies the transformation.


* **Parameters**


    * **structure** – Input Structure


    * **return_ranked_list** (*bool** | **int**, **optional*) – If return_ranked_list is int, that number of structures

    is returned. If False, only the single lowest energy structure is returned. Defaults to False.




* **Returns**

    Structures with all substitutions applied.



#### _property_ inverse()
None.


* **Type**

    Returns



#### _property_ is_one_to_many(_: boo_ )
True.


* **Type**

    Returns



### _class_ SQSTransformation(scaling, cluster_size_and_shell=None, search_time=60, directory=None, instances=None, temperature=1, wr=1, wn=1, wd=0.5, tol=0.001, best_only=True, remove_duplicate_structures=True, reduction_algo='LLL')
Bases: `AbstractTransformation`

A transformation that creates a special quasirandom structure (SQS) from a structure with partial occupancies.


* **Parameters**


    * **scaling** (*int** or **list*) – Scaling factor to determine supercell. Two options are possible:
    a. (preferred) Scales number of atoms, e.g., for a structure with 8 atoms,

    > scaling=4 would lead to a 32 atom supercell


        1. A sequence of three scaling factors, e.g., [2, 1, 1], which
    specifies that the supercell should have dimensions 2a x b x c



    * **cluster_size_and_shell** (*Optional**[**Dict**[**int**, **int**]**]*) – Dictionary of cluster interactions with entries in
    the form number of atoms: nearest neighbor shell



* **Keyword Arguments**


    * **search_time** (*float*) – Time spent looking for the ideal SQS in minutes (default: 60)


    * **directory** (*str*) – Directory to run mcsqs calculation and store files (default: None
    runs calculations in a temp directory)


    * **instances** (*int*) – Specifies the number of parallel instances of mcsqs to run
    (default: number of cpu cores detected by Python)


    * **temperature** (*int** or **float*) – Monte Carlo temperature (default: 1), “T” in atat code


    * **wr** (*int** or **float*) – Weight assigned to range of perfect correlation match in objective
    function (default = 1)


    * **wn** (*int** or **float*) – Multiplicative decrease in weight per additional point in cluster (default: 1)


    * **wd** (*int** or **float*) – Exponent of decay in weight as function of cluster diameter (default: 0)


    * **tol** (*int** or **float*) – Tolerance for matching correlations (default: 1e-3)


    * **best_only** (*bool*) – only return structures with lowest objective function


    * **remove_duplicate_structures** (*bool*) – only return unique structures


    * **reduction_algo** (*str*) – The lattice reduction algorithm to use.
    Currently supported options are “niggli” or “LLL”.
    “False” does not reduce structure.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### _static_ _get_disordered_substructure(struc_disordered)
Converts disordered structure into a substructure consisting of only disordered sites.


* **Parameters**

    **struc_disordered** – pymatgen disordered Structure object.



* **Returns**

    pymatgen Structure object representing a substructure of disordered sites.



#### _static_ _get_max_neighbor_distance(struct, shell)
Calculate maximum nearest neighbor distance.


* **Parameters**


    * **struct** – pymatgen Structure object


    * **shell** – nearest neighbor shell, such that shell=1 is the first nearest
    neighbor, etc.



* **Returns**

    maximum nearest neighbor distance, in angstroms



#### _static_ _get_unique_bestsqs_strucs(sqs, best_only, return_ranked_list, remove_duplicate_structures, reduction_algo)
Gets unique sqs structures with lowest objective function. Requires an mcsqs output that has been run

    in parallel, otherwise returns Sqs.bestsqs.


* **Parameters**


    * **sqs** ([*Sqs*](pymatgen.command_line.md#pymatgen.command_line.mcsqs_caller.Sqs)) – Sqs class object.


    * **best_only** (*bool*) – only return structures with lowest objective function.


    * **return_ranked_list** (*bool** | **int**, **optional*) – If return_ranked_list is int, that number of structures.

    is returned. If False, only the single lowest energy structure is returned. Defaults to False.



    * **remove_duplicate_structures** (*bool*) – only return unique structures.


    * **reduction_algo** (*str*) – The lattice reduction algorithm to use.
    Currently supported options are “niggli” or “LLL”.
    “False” does not reduce structure.



* **Returns**

    Structure, ‘objective_function’: …}, unless run in serial

        (returns a single structure Sqs.bestsqs)




* **Return type**

    list of dicts of the form {‘structure’



#### _static_ _sqs_cluster_estimate(struc_disordered, cluster_size_and_shell: dict[int, int] | None = None)
Set up an ATAT cluster.out file for a given structure and set of constraints.


* **Parameters**


    * **struc_disordered** – disordered pymatgen Structure object


    * **cluster_size_and_shell** – dict of integers {cluster: shell}.



* **Returns**

    distance in angstroms} for mcsqs calculation



* **Return type**

    dict of {cluster size



#### apply_transformation(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), return_ranked_list: bool | int = False)
Applies SQS transformation.


* **Parameters**


    * **structure** (*pymatgen Structure*) – pymatgen Structure with partial occupancies


    * **return_ranked_list** (*bool** | **int**, **optional*) – If return_ranked_list is int, that number of structures

    is returned. If False, only the single lowest energy structure is returned. Defaults to False.




* **Returns**

    pymatgen Structure which is an SQS of the input structure



#### _property_ inverse()
None.


* **Type**

    Returns



#### _property_ is_one_to_many(_: boo_ )
True.


* **Type**

    Returns



### _class_ SlabTransformation(miller_index, min_slab_size, min_vacuum_size, lll_reduce=False, center_slab=False, in_unit_planes=False, primitive=True, max_normal_search=None, shift=0, tol=0.1)
Bases: `AbstractTransformation`

A transformation that creates a slab from a structure.


* **Parameters**


    * **miller_index** (*3-tuple** or **list*) – miller index of slab


    * **min_slab_size** (*float*) – minimum slab size in angstroms


    * **min_vacuum_size** (*float*) – minimum size of vacuum


    * **lll_reduce** (*bool*) – whether to apply LLL reduction


    * **center_slab** (*bool*) – whether to center the slab


    * **primitive** (*bool*) – whether to reduce slabs to most primitive cell


    * **in_unit_planes** (*bool*) – Whether to set min_slab_size and min_vac_size
    in units of hkl planes (True) or Angstrom (False, the default). Setting in
    units of planes is useful for ensuring some slabs have a certain n_layer of
    atoms. e.g. for Cs (100), a 10 Ang slab will result in a slab with only 2
    layer of atoms, whereas Fe (100) will have more layer of atoms. By using units
    of hkl planes instead, we ensure both slabs have the same number of atoms. The
    slab thickness will be in min_slab_size/math.ceil(self._proj_height/dhkl)
    multiples of oriented unit cells.


    * **max_normal_search** (*int*) – maximum index to include in linear
    combinations of indices to find c lattice vector orthogonal
    to slab surface


    * **shift** (*float*) – shift to get termination


    * **tol** (*float*) – tolerance for primitive cell finding.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### apply_transformation(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))
Applies the transformation.


* **Parameters**

    **structure** – Input Structure



* **Returns**

    Slab Structures.



#### _property_ inverse()
None.


* **Type**

    Returns



#### _property_ is_one_to_many(_: boo_ )
False.


* **Type**

    Returns



### _class_ SubstituteSurfaceSiteTransformation(atom, selective_dynamics=False, height=0.9, mi_vec=None, target_species=None, sub_both_sides=False, range_tol=0.01, dist_from_surf=0)
Bases: `AbstractTransformation`

Use AdsorptionSiteFinder to perform substitution-type doping on the surface
and returns all possible configurations where one dopant is substituted
per surface. Can substitute one surface or both.


* **Parameters**


    * **atom** (*str*) – atom corresponding to substitutional dopant


    * **selective_dynamics** (*bool*) – flag for whether to assign
    non-surface sites as fixed for selective dynamics


    * **height** (*float*) – height criteria for selection of surface sites


    * **mi_vec** – vector corresponding to the vector
    concurrent with the miller index, this enables use with
    slabs that have been reoriented, but the miller vector
    must be supplied manually


    * **target_species** – List of specific species to substitute


    * **sub_both_sides** (*bool*) – If true, substitute an equivalent
    site on the other surface


    * **range_tol** (*float*) – Find viable substitution sites at a specific
    distance from the surface +- this tolerance


    * **dist_from_surf** (*float*) – Distance from the surface to find viable
    substitution sites, defaults to 0 to substitute at the surface.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### apply_transformation(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), return_ranked_list: bool | int = False)

* **Parameters**


    * **structure** – Must be a Slab structure


    * **return_ranked_list** (*bool** | **int**, **optional*) – If return_ranked_list is int, that number of structures.

    is returned. If False, only the single lowest energy structure is returned. Defaults to False.




* **Returns**

    each dict has key ‘structure’ which is a Slab with sites substituted



* **Return type**

    list[dict]



#### _property_ inverse()
None.


* **Type**

    Returns



#### _property_ is_one_to_many(_: boo_ )
True.


* **Type**

    Returns



### _class_ SubstitutionPredictorTransformation(threshold=0.01, scale_volumes=True, \*\*kwargs)
Bases: `AbstractTransformation`

This transformation takes a structure and uses the structure
prediction module to find likely site substitutions.


* **Parameters**


    * **threshold** – Threshold for substitution.


    * **scale_volumes** – Whether to scale volumes after substitution.


    * **\*\*kwargs** – Args for SubstitutionProbability class lambda_table, alpha.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### apply_transformation(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), return_ranked_list: bool | int = False)
Applies the transformation.


* **Parameters**


    * **structure** – Input Structure


    * **return_ranked_list** (*bool** | **int**, **optional*) – If return_ranked_list is int, that number of structures

    is returned. If False, only the single lowest energy structure is returned. Defaults to False.




* **Returns**

    Predicted Structures.



#### _property_ inverse()
None.


* **Type**

    Returns



#### _property_ is_one_to_many(_: boo_ )
True.


* **Type**

    Returns



### _class_ SuperTransformation(transformations, nstructures_per_trans=1)
Bases: `AbstractTransformation`

This is a transformation that is inherently one-to-many. It is constructed
from a list of transformations and returns one structure for each
transformation. The primary use for this class is extending a transmuter
object.


* **Parameters**


    * **transformations** (*[**transformations**]*) – List of transformations to apply
    to a structure. One transformation is applied to each output
    structure.


    * **nstructures_per_trans** (*int*) – If the transformations are one-to-many and,
    nstructures_per_trans structures from each transformation are
    added to the full list. Defaults to 1, i.e., only best structure.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### apply_transformation(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), return_ranked_list: bool | int = False)
Applies the transformation.


* **Parameters**


    * **structure** – Input Structure


    * **return_ranked_list** (*bool** | **int**, **optional*) – If return_ranked_list is int, that number of structures

    is returned. If False, only the single lowest energy structure is returned. Defaults to False.




* **Returns**

    Structures with all transformations applied.



#### _property_ inverse()
None.


* **Type**

    Returns



#### _property_ is_one_to_many(_: boo_ )
True.


* **Type**

    Returns



### _find_codopant(target, oxidation_state, allowed_elements=None)
Finds the element from “allowed elements” that (i) possesses the desired
“oxidation state” and (ii) is closest in ionic radius to the target specie.


* **Parameters**


    * **target** – (Species) provides target ionic radius.


    * **oxidation_state** – (float) codopant oxidation state.


    * **allowed_elements** – ([str]) List of allowed elements. If None,
    all elements are tried.



* **Returns**

    (Species) with oxidation_state that has ionic radius closest to
    target.



### _proj(b, a)
Returns vector projection (np.ndarray) of vector b (np.ndarray)
onto vector a (np.ndarray).


### _round_and_make_arr_singular(arr: ndarray)
This function rounds all elements of a matrix to the nearest integer,
unless the rounding scheme causes the matrix to be singular, in which
case elements of zero rows or columns in the rounded matrix with the
largest absolute valued magnitude in the unrounded matrix will be
rounded to the next integer away from zero rather than to the
nearest integer.

The transformation is as follows. First, all entries in ‘arr’ will be
rounded to the nearest integer to yield ‘arr_rounded’. If ‘arr_rounded’
has any zero rows, then one element in each zero row of ‘arr_rounded’
corresponding to the element in ‘arr’ of that row with the largest
absolute valued magnitude will be rounded to the next integer away from
zero (see the ‘_round_away_from_zero(x)’ function) rather than the
nearest integer. This process is then repeated for zero columns. Also
note that if ‘arr’ already has zero rows or columns, then this function
will not change those rows/columns.


* **Parameters**

    **arr** – Input matrix



* **Returns**

    Transformed matrix.


## pymatgen.transformations.site_transformations module

This module defines site transformations which transforms a structure into
another structure. Site transformations differ from standard transformations
in that they operate in a site-specific manner.
All transformations should inherit the AbstractTransformation ABC.


### _class_ AddSitePropertyTransformation(site_properties)
Bases: `AbstractTransformation`

Simple transformation to add site properties to a given structure.


* **Parameters**

    **site_properties** (*dict*) – site properties to be added to a structure.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### apply_transformation(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))
Apply the transformation.

Arg:

    structure (Structure): A structurally similar structure in

        regards to crystal and site positions.


* **Returns**

    Returns a copy of structure with sites properties added.



#### _property_ inverse()
None.


* **Type**

    Return



#### _property_ is_one_to_many(_: boo_ )
False.


* **Type**

    Return



### _class_ InsertSitesTransformation(species, coords, coords_are_cartesian=False, validate_proximity=True)
Bases: `AbstractTransformation`

This transformation substitutes certain sites with certain species.


* **Parameters**


    * **species** – A list of species. e.g., [“Li”, “Fe”]


    * **coords** – A list of coords corresponding to those species. e.g.,
    [[0,0,0],[0.5,0.5,0.5]].


    * **coords_are_cartesian** (*bool*) – Set to True if coords are given in
    Cartesian coords. Defaults to False.


    * **validate_proximity** (*bool*) – Set to False if you do not wish to ensure
    that added sites are not too close to other sites. Defaults to True.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### apply_transformation(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))
Apply the transformation.

Arg:

    structure (Structure): A structurally similar structure in

        regards to crystal and site positions.


* **Returns**

    Returns a copy of structure with sites inserted.



#### _property_ inverse()
None.


* **Type**

    Return



#### _property_ is_one_to_many(_: boo_ )
False.


* **Type**

    Return



### _class_ PartialRemoveSitesTransformation(indices, fractions, algo=1)
Bases: `AbstractTransformation`

Remove fraction of specie from a structure.
Requires an oxidation state decorated structure for Ewald sum to be
computed.

Given that the solution to selecting the right removals is NP-hard, there
are several algorithms provided with varying degrees of accuracy and speed.
The options are as follows:

ALGO_FAST:

    This is a highly optimized algorithm to quickly go through the search
    tree. It is guaranteed to find the optimal solution, but will return
    only a single lowest energy structure. Typically, you will want to use
    this.

ALGO_COMPLETE:

    The complete algo ensures that you get all symmetrically distinct
    orderings, ranked by the estimated Ewald energy. But this can be an
    extremely time-consuming process if the number of possible orderings is
    very large. Use this if you really want all possible orderings. If you
    want just the lowest energy ordering, ALGO_FAST is accurate and faster.

ALGO_BEST_FIRST:

    This algorithm is for ordering the really large cells that defeats even
    ALGO_FAST. For example, if you have 48 sites of which you want to
    remove 16 of them, the number of possible orderings is around
    2 x 10^12. ALGO_BEST_FIRST shortcircuits the entire search tree by
    removing the highest energy site first, then followed by the next
    highest energy site, and so on. It is guaranteed to find a solution
    in a reasonable time, but it is also likely to be highly inaccurate.

ALGO_ENUMERATE:

    This algorithm uses the EnumerateStructureTransformation to perform
    ordering. This algo returns *complete* orderings up to a single unit
    cell size. It is more robust than the ALGO_COMPLETE, but requires
    Gus Hart’s enumlib to be installed.


* **Parameters**


    * **indices** – A list of list of indices.
    e.g. [[0, 1], [2, 3, 4, 5]]


    * **fractions** – The corresponding fractions to remove. Must be same length as
    indices. e.g., [0.5, 0.25]


    * **algo** – This parameter allows you to choose the algorithm to perform
    ordering. Use one of PartialRemoveSpecieTransformation.ALGO_\*
    variables to set the algo.



#### ALGO_BEST_FIRST(_ = _ )

#### ALGO_COMPLETE(_ = _ )

#### ALGO_ENUMERATE(_ = _ )

#### ALGO_FAST(_ = _ )

#### _abc_impl(_ = <_abc._abc_data object_ )

#### _best_first_ordering(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), num_remove_dict)

#### _complete_ordering(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), num_remove_dict)

#### _enumerate_ordering(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))

#### _fast_ordering(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), num_remove_dict, num_to_return=1)
This method uses the matrix form of ewaldsum to calculate the ewald
sums of the potential structures. This is on the order of 4 orders of
magnitude faster when there are large numbers of permutations to
consider. There are further optimizations possible (doing a smarter
search of permutations for example), but this won’t make a difference
until the number of permutations is on the order of 30,000.


#### apply_transformation(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), return_ranked_list: bool | int = False)
Apply the transformation.


* **Parameters**


    * **structure** – input structure


    * **return_ranked_list** (*bool** | **int*) – Whether or not multiple structures are returned.
    If return_ranked_list is int, that number of structures is returned.



* **Returns**

    Depending on returned_ranked list, either a transformed structure
    or a list of dictionaries, where each dictionary is of the form
    {“structure” = …. , “other_arguments”}
    the key “transformation” is reserved for the transformation that
    was actually applied to the structure.
    This transformation is parsed by the alchemy classes for generating
    a more specific transformation history. Any other information will
    be stored in the transformation_parameters dictionary in the
    transmuted structure class.



#### _property_ inverse()
None.


* **Type**

    Return



#### _property_ is_one_to_many(_: boo_ )
True.


* **Type**

    Return



### _class_ RadialSiteDistortionTransformation(site_index, displacement=0.1, nn_only=False)
Bases: `AbstractTransformation`

Radially perturbs atoms around a site. Can be used to create spherical distortion due to a
point defect.


* **Parameters**


    * **site_index** (*int*) – index of the site in structure to place at the center of the distortion (will
    not be distorted). This index must be provided before the structure is provided in
    apply_transformation in order to keep in line with the base class.


    * **displacement** (*float*) – distance to perturb the atoms around the objective site


    * **nn_only** (*bool*) – Whether or not to perturb beyond the nearest neighbors. If True, then only the
    nearest neighbors will be perturbed, leaving the other sites undisturbed. If False, then
    the nearest neighbors will receive the full displacement, and then subsequent sites will receive
    a displacement=0.1 / r, where r is the distance each site to the origin site. For small displacements,
    atoms beyond the NN environment will receive very small displacements, and these are almost equal.
    For large displacements, this difference is noticeable.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### apply_transformation(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))
Apply the transformation.


* **Parameters**

    **structure** – Structure or Molecule to apply the transformation to



* **Returns**

    the transformed structure



#### _property_ inverse()
Returns the inverse transformation if available.
Otherwise, should return None.


#### _property_ is_one_to_many(_: boo_ )
Determines if a Transformation is a one-to-many transformation. If a
Transformation is a one-to-many transformation, the
apply_transformation method should have a keyword arg
“return_ranked_list” which allows for the transformed structures to be
returned as a ranked list.


#### _property_ use_multiprocessing()
Indicates whether the transformation can be applied by a
subprocessing pool. This should be overridden to return True for
transformations that the transmuter can parallelize.


### _class_ RemoveSitesTransformation(indices_to_remove)
Bases: `AbstractTransformation`

Remove certain sites in a structure.


* **Parameters**

    **indices_to_remove** – List of indices to remove. E.g., [0, 1, 2].



#### _abc_impl(_ = <_abc._abc_data object_ )

#### apply_transformation(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))
Apply the transformation.

Arg:

    structure (Structure): A structurally similar structure in

        regards to crystal and site positions.


* **Returns**

    Returns a copy of structure with sites removed.



#### _property_ inverse()
None.


* **Type**

    Return



#### _property_ is_one_to_many(_: boo_ )
False.


* **Type**

    Return



### _class_ ReplaceSiteSpeciesTransformation(indices_species_map)
Bases: `AbstractTransformation`

This transformation substitutes certain sites with certain species.


* **Parameters**

    **indices_species_map** – A dict containing the species mapping in
    int-string pairs. E.g., { 1:”Na”} or {2:”Mn2+”}. Multiple
    substitutions can be done. Overloaded to accept sp_and_occu
    dictionary. E.g. {1: {“Ge”:0.75, “C”:0.25} }, which
    substitutes a single species with multiple species to generate a
    disordered structure.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### apply_transformation(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))
Apply the transformation.

Arg:

    structure (Structure): A structurally similar structure in

        regards to crystal and site positions.


* **Returns**

    Returns a copy of structure with sites replaced.



#### _property_ inverse()
None.


* **Type**

    Return



#### _property_ is_one_to_many(_: boo_ )
False.


* **Type**

    Return



### _class_ TranslateSitesTransformation(indices_to_move, translation_vector, vector_in_frac_coords=True)
Bases: `AbstractTransformation`

This class translates a set of sites by a certain vector.


* **Parameters**


    * **indices_to_move** – The indices of the sites to move


    * **translation_vector** – Vector to move the sites. If a list of list or numpy
    array of shape, (len(indices_to_move), 3), is provided then each
    translation vector is applied to the corresponding site in the
    indices_to_move.


    * **vector_in_frac_coords** – Set to True if the translation vector is in
    fractional coordinates, and False if it is in cartesian
    coordinations. Defaults to True.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### apply_transformation(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))
Apply the transformation.

Arg:

    structure (Structure): A structurally similar structure in

        regards to crystal and site positions.


* **Returns**

    Returns a copy of structure with sites translated.



#### as_dict()
JSON-serializable dict representation.


#### _property_ inverse()
Returns:
TranslateSitesTransformation with the reverse translation.


#### _property_ is_one_to_many(_: boo_ )
False.


* **Type**

    Return


## pymatgen.transformations.standard_transformations module

This module defines standard transformations which transforms a structure into
another structure. Standard transformations operate in a structure-wide manner,
rather than site-specific manner.
All transformations should inherit the AbstractTransformation ABC.


### _class_ AutoOxiStateDecorationTransformation(symm_tol=0.1, max_radius=4, max_permutations=100000, distance_scale_factor=1.015)
Bases: `AbstractTransformation`

This transformation automatically decorates a structure with oxidation
states using a bond valence approach.


* **Parameters**


    * **symm_tol** (*float*) – Symmetry tolerance used to determine which sites are
    symmetrically equivalent. Set to 0 to turn off symmetry.


    * **max_radius** (*float*) – Maximum radius in Angstrom used to find nearest
    neighbors.


    * **max_permutations** (*int*) – Maximum number of permutations of oxidation
    states to test.


    * **distance_scale_factor** (*float*) – A scale factor to be applied. This is
    useful for scaling distances, esp in the case of
    calculation-relaxed structures, which may tend to under (GGA) or
    over bind (LDA). The default of 1.015 works for GGA. For
    experimental structure, set this to 1.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### apply_transformation(structure)
Apply the transformation.


* **Parameters**

    **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Input Structure



* **Returns**

    Oxidation state decorated Structure.



#### _property_ inverse()
None.


* **Type**

    Returns



#### _property_ is_one_to_many(_: boo_ )
False.


* **Type**

    Returns



### _class_ ChargedCellTransformation(charge=0)
Bases: `AbstractTransformation`

The ChargedCellTransformation applies a charge to a structure (or defect
object).


* **Parameters**

    **charge** – A integer charge to apply to the structure.
    Defaults to zero. Has to be a single integer. e.g. 2.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### apply_transformation(structure)
Apply the transformation.


* **Parameters**

    **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Input Structure



* **Returns**

    Charged Structure.



#### _property_ inverse()
NotImplementedError.


* **Type**

    Raises



#### _property_ is_one_to_many(_: boo_ )
False.


* **Type**

    Returns



### _class_ ConventionalCellTransformation(symprec: float = 0.01, angle_tolerance=5, international_monoclinic=True)
Bases: `AbstractTransformation`

This class finds the conventional cell of the input structure.


* **Parameters**


    * **symprec** (*float*) – tolerance as in SpacegroupAnalyzer


    * **angle_tolerance** (*float*) – angle tolerance as in SpacegroupAnalyzer


    * **international_monoclinic** (*bool*) – whether to use beta (True) or alpha (False)


as the non-right-angle in the unit cell.


#### _abc_impl(_ = <_abc._abc_data object_ )

#### apply_transformation(structure)
Returns most primitive cell for structure.


* **Parameters**

    **structure** – A structure



* **Returns**

    The same structure in a conventional standard setting



#### _property_ inverse()
None.


* **Type**

    Returns



#### _property_ is_one_to_many(_: boo_ )
False.


* **Type**

    Returns



### _class_ DeformStructureTransformation(deformation=((1, 0, 0), (0, 1, 0), (0, 0, 1)))
Bases: `AbstractTransformation`

This transformation deforms a structure by a deformation gradient matrix.


* **Parameters**

    **deformation** (*array*) – deformation gradient for the transformation.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### apply_transformation(structure)
Apply the transformation.


* **Parameters**

    **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Input Structure



* **Returns**

    Deformed Structure.



#### _property_ inverse()
Returns inverse Transformation.


#### _property_ is_one_to_many(_: boo_ )
False.


* **Type**

    Returns



### _class_ DiscretizeOccupanciesTransformation(max_denominator=5, tol: float | None = None, fix_denominator=False)
Bases: `AbstractTransformation`

Discretizes the site occupancies in a disordered structure; useful for
grouping similar structures or as a pre-processing step for order-disorder
transformations.


* **Parameters**


    * **max_denominator** – An integer maximum denominator for discretization. A higher
    denominator allows for finer resolution in the site occupancies.


    * **tol** – A float that sets the maximum difference between the original and
    discretized occupancies before throwing an error. If None, it is
    set to 1 / (4 \* max_denominator).


    * **fix_denominator** (*bool*) – If True, will enforce a common denominator for all species.
    This prevents a mix of denominators (for example, 1/3, 1/4)
    that might require large cell sizes to perform an enumeration.
    ‘tol’ needs to be > 1.0 in some cases.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### apply_transformation(structure)
Discretizes the site occupancies in the structure.


* **Parameters**

    **structure** – disordered Structure to discretize occupancies



* **Returns**

    A new disordered Structure with occupancies discretized



#### _property_ inverse()
None.


* **Type**

    Returns



#### _property_ is_one_to_many(_: boo_ )
False.


* **Type**

    Returns



### _class_ OrderDisorderedStructureTransformation(algo=0, symmetrized_structures=False, no_oxi_states=False)
Bases: `AbstractTransformation`

Order a disordered structure. The disordered structure must be oxidation
state decorated for Ewald sum to be computed. No attempt is made to perform
symmetry determination to reduce the number of combinations.

Hence, attempting to performing ordering on a large number of disordered
sites may be extremely expensive. The time scales approximately with the
number of possible combinations. The algorithm can currently compute
approximately 5,000,000 permutations per minute.

Also, simple rounding of the occupancies are performed, with no attempt
made to achieve a target composition. This is usually not a problem for
most ordering problems, but there can be times where rounding errors may
result in structures that do not have the desired composition.
This second step will be implemented in the next iteration of the code.

If multiple fractions for a single species are found for different sites,
these will be treated separately if the difference is above a threshold
tolerance. currently this is .1

For example, if a fraction of .25 Li is on sites 0,1,2,3  and .5 on sites
4, 5, 6, 7 then 1 site from [0,1,2,3] will be filled and 2 sites from [4,5,6,7]
will be filled, even though a lower energy combination might be found by
putting all lithium in sites [4,5,6,7].

USE WITH CARE.


* **Parameters**


    * **algo** (*int*) – Algorithm to use.


    * **symmetrized_structures** (*bool*) – Whether the input structures are
    instances of SymmetrizedStructure, and that their symmetry
    should be used for the grouping of sites.


    * **no_oxi_states** (*bool*) – Whether to remove oxidation states prior to
    ordering.



#### ALGO_BEST_FIRST(_ = _ )

#### ALGO_COMPLETE(_ = _ )

#### ALGO_FAST(_ = _ )

#### _abc_impl(_ = <_abc._abc_data object_ )

#### apply_transformation(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), return_ranked_list: bool | int = False)
For this transformation, the apply_transformation method will return
only the ordered structure with the lowest Ewald energy, to be
consistent with the method signature of the other transformations.
However, all structures are stored in the  all_structures attribute in
the transformation object for easy access.


* **Parameters**


    * **structure** – Oxidation state decorated disordered structure to order


    * **return_ranked_list** (*bool** | **int**, **optional*) – If return_ranked_list is int, that number of structures

    is returned. If False, only the single lowest energy structure is returned. Defaults to False.




* **Returns**

    Depending on returned_ranked list, either a transformed structure
    or a list of dictionaries, where each dictionary is of the form
    {“structure” = …. , “other_arguments”}
    the key “transformation” is reserved for the transformation that
    was actually applied to the structure.
    This transformation is parsed by the alchemy classes for generating
    a more specific transformation history. Any other information will
    be stored in the transformation_parameters dictionary in the
    transmuted structure class.



#### _property_ inverse()
None.


* **Type**

    Returns



#### _property_ is_one_to_many(_: boo_ )
True.


* **Type**

    Returns



#### _property_ lowest_energy_structure()
Lowest energy structure found.


### _class_ OxidationStateDecorationTransformation(oxidation_states)
Bases: `AbstractTransformation`

This transformation decorates a structure with oxidation states.


* **Parameters**


    * **oxidation_states** (*dict*) – Oxidation states supplied as a dict,


    * **e.g.** – 1, “O”:-2}.


    * **{"Li"** – 1, “O”:-2}.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### apply_transformation(structure)
Apply the transformation.


* **Parameters**

    **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Input Structure



* **Returns**

    Oxidation state decorated Structure.



#### _property_ inverse()
None.


* **Type**

    Returns



#### _property_ is_one_to_many(_: boo_ )
False.


* **Type**

    Returns



### _class_ OxidationStateRemovalTransformation()
Bases: `AbstractTransformation`

This transformation removes oxidation states from a structure.

No arg needed.


#### _abc_impl(_ = <_abc._abc_data object_ )

#### apply_transformation(structure)
Apply the transformation.


* **Parameters**

    **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Input Structure



* **Returns**

    Non-oxidation state decorated Structure.



#### _property_ inverse()
None.


* **Type**

    Returns



#### _property_ is_one_to_many(_: boo_ )
False.


* **Type**

    Returns



### _class_ PartialRemoveSpecieTransformation(specie_to_remove, fraction_to_remove, algo=0)
Bases: `AbstractTransformation`

Remove fraction of specie from a structure.

Requires an oxidation state decorated structure for Ewald sum to be
computed.

Given that the solution to selecting the right removals is NP-hard, there
are several algorithms provided with varying degrees of accuracy and speed.
Please see
pymatgen.transformations.site_transformations.PartialRemoveSitesTransformation.


* **Parameters**


    * **specie_to_remove** – Species to remove. Must have oxidation state E.g.,
    “Li+”


    * **fraction_to_remove** – Fraction of specie to remove. E.g., 0.5


    * **algo** – This parameter allows you to choose the algorithm to perform
    ordering. Use one of PartialRemoveSpecieTransformation.ALGO_\*
    variables to set the algo.



#### ALGO_BEST_FIRST(_ = _ )

#### ALGO_COMPLETE(_ = _ )

#### ALGO_ENUMERATE(_ = _ )

#### ALGO_FAST(_ = _ )

#### _abc_impl(_ = <_abc._abc_data object_ )

#### apply_transformation(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), return_ranked_list: bool | int = False)
Apply the transformation.


* **Parameters**


    * **structure** – input structure


    * **return_ranked_list** (*bool** | **int**, **optional*) – If return_ranked_list is int, that number of structures

    is returned. If False, only the single lowest energy structure is returned. Defaults to False.




* **Returns**

    Depending on returned_ranked list, either a transformed structure
    or a list of dictionaries, where each dictionary is of the form
    {“structure” = …. , “other_arguments”}
    the key “transformation” is reserved for the transformation that
    was actually applied to the structure.
    This transformation is parsed by the alchemy classes for generating
    a more specific transformation history. Any other information will
    be stored in the transformation_parameters dictionary in the
    transmuted structure class.



#### _property_ inverse()
None.


* **Type**

    Returns



#### _property_ is_one_to_many(_: boo_ )
True.


* **Type**

    Returns



### _class_ PerturbStructureTransformation(distance: float = 0.01, min_distance: float | None = None)
Bases: `AbstractTransformation`

This transformation perturbs a structure by a specified distance in random
directions. Used for breaking symmetries.


* **Parameters**


    * **distance** – Distance of perturbation in angstroms. All sites
    will be perturbed by exactly that distance in a random
    direction.


    * **min_distance** – if None, all displacements will be equidistant. If int
    or float, perturb each site a distance drawn from the uniform
    distribution between ‘min_distance’ and ‘distance’.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### apply_transformation(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))
Apply the transformation.


* **Parameters**

    **structure** – Input Structure



* **Returns**

    Structure with sites perturbed.



#### _property_ inverse()
None.


* **Type**

    Returns



#### _property_ is_one_to_many(_: boo_ )
False.


* **Type**

    Returns



### _class_ PrimitiveCellTransformation(tolerance=0.5)
Bases: `AbstractTransformation`

This class finds the primitive cell of the input structure.
It returns a structure that is not necessarily orthogonalized
Author: Will Richards.


* **Parameters**

    **tolerance** (*float*) – Tolerance for each coordinate of a particular
    site. For example, [0.5, 0, 0.5] in Cartesian coordinates will be
    considered to be on the same coordinates as [0, 0, 0] for a
    tolerance of 0.5. Defaults to 0.5.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### apply_transformation(structure)
Returns most primitive cell for structure.


* **Parameters**

    **structure** – A structure



* **Returns**

    The most primitive structure found. The returned structure is
    guaranteed to have len(new structure) <= len(structure).



#### _property_ inverse()
None.


* **Type**

    Returns



#### _property_ is_one_to_many(_: boo_ )
False.


* **Type**

    Returns



### _class_ RemoveSpeciesTransformation(species_to_remove)
Bases: `AbstractTransformation`

Remove all occurrences of some species from a structure.


* **Parameters**

    **species_to_remove** – List of species to remove. E.g., [“Li”, “Mn”].



#### _abc_impl(_ = <_abc._abc_data object_ )

#### apply_transformation(structure)
Apply the transformation.


* **Parameters**

    **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Input Structure



* **Returns**

    Structure with species removed.



#### _property_ inverse()
None.


* **Type**

    Returns



#### _property_ is_one_to_many(_: boo_ )
False.


* **Type**

    Returns



### _class_ RotationTransformation(axis, angle, angle_in_radians=False)
Bases: `AbstractTransformation`

The RotationTransformation applies a rotation to a structure.


* **Parameters**


    * **axis** (*3x1 array*) – Axis of rotation, e.g., [1, 0, 0]


    * **angle** (*float*) – Angle to rotate


    * **angle_in_radians** (*bool*) – Set to True if angle is supplied in radians.
    Else degrees are assumed.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### apply_transformation(structure)
Apply the transformation.


* **Parameters**

    **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Input Structure



* **Returns**

    Rotated Structure.



#### _property_ inverse()
Returns inverse Transformation.


#### _property_ is_one_to_many(_: boo_ )
False.


* **Type**

    Returns



### _class_ ScaleToRelaxedTransformation(unrelaxed_structure, relaxed_structure, species_map=None)
Bases: `AbstractTransformation`

Takes the unrelaxed and relaxed structure and applies its site and volume
relaxation to a structurally similar structures (e.g. bulk: NaCl and PbTe
(rock-salt), slab: Sc(10-10) and Mg(10-10) (hcp), GB: Mo(001) sigma 5 GB,
Fe(001) sigma 5). Useful for finding an initial guess of a set of similar
structures closer to its most relaxed state.


* **Parameters**


    * **unrelaxed_structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Initial, unrelaxed structure


    * **relaxed_structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Relaxed structure


    * **species_map** (*dict*) – A dict or list of tuples containing the species mapping in
    string-string pairs. The first species corresponds to the relaxed
    structure while the second corresponds to the species in the
    structure to be scaled. E.g., {“Li”:”Na”} or [(“Fe2+”,”Mn2+”)].
    Multiple substitutions can be done. Overloaded to accept
    sp_and_occu dictionary E.g. {“Si: {“Ge”:0.75, “C”:0.25}},
    which substitutes a single species with multiple species to
    generate a disordered structure.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### apply_transformation(structure)
Returns a copy of structure with lattice parameters
and sites scaled to the same degree as the relaxed_structure.

Arg:

    structure (Structure): A structurally similar structure in

        regards to crystal and site positions.


#### _property_ inverse()
None.


* **Type**

    Returns



#### _property_ is_one_to_many(_: boo_ )
False.


* **Type**

    Returns



### _class_ SubstitutionTransformation(species_map: dict[SpeciesLike, SpeciesLike | dict[SpeciesLike, float]] | list[tuple[SpeciesLike, SpeciesLike]])
Bases: `AbstractTransformation`

This transformation substitutes species for one another.


* **Parameters**

    **species_map** – A dict or list of tuples containing the species mapping in
    string-string pairs. E.g., {“Li”: “Na”} or [(“Fe2+”,”Mn2+”)].
    Multiple substitutions can be done. Overloaded to accept
    sp_and_occu dictionary E.g. {“Si: {“Ge”:0.75, “C”:0.25}},
    which substitutes a single species with multiple species to
    generate a disordered structure.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### apply_transformation(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))
Apply the transformation.


* **Parameters**

    **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Input Structure



* **Returns**

    Substituted Structure.



#### _property_ inverse()
Returns inverse Transformation.


#### _property_ is_one_to_many(_: boo_ )
False.


* **Type**

    Returns



### _class_ SupercellTransformation(scaling_matrix=((1, 0, 0), (0, 1, 0), (0, 0, 1)))
Bases: `AbstractTransformation`

The SupercellTransformation replicates an unitcell to a supercell.


* **Parameters**

    **scaling_matrix** – A matrix of transforming the lattice vectors.
    Defaults to the identity matrix. Has to be all integers. e.g.,
    [[2,1,0],[0,3,0],[0,0,1]] generates a new structure with
    lattice vectors a” = 2a + b, b” = 3b, c” = c where a, b, and c
    are the lattice vectors of the original structure.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### apply_transformation(structure)
Apply the transformation.


* **Parameters**

    **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Input Structure



* **Returns**

    Supercell Structure.



#### _static_ from_boundary_distance(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), min_boundary_dist: float = 6, allow_rotation: bool = False, max_atoms: float = -1)
Get a SupercellTransformation according to the desired minimum distance between periodic
boundaries of the resulting supercell.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Input structure.


    * **min_boundary_dist** (*float*) – Desired minimum distance between all periodic boundaries. Defaults to 6.


    * **allow_rotation** (*bool*) – Whether allowing lattice angles to change. Only useful when
    at least two of the three lattice vectors are required to expand. Defaults to False.
    If True, a SupercellTransformation satisfying min_boundary_dist but with smaller
    number of atoms than the SupercellTransformation with unchanged lattice angles
    can possibly be found. If such a SupercellTransformation cannot be found easily,
    the SupercellTransformation with unchanged lattice angles will be returned.


    * **max_atoms** (*int*) – Maximum number of atoms allowed in the supercell. Defaults to infinity.



* **Returns**

    SupercellTransformation.



#### _static_ from_scaling_factors(scale_a=1, scale_b=1, scale_c=1)
Convenience method to get a SupercellTransformation from a simple
series of three numbers for scaling each lattice vector. Equivalent to
calling the normal with [[scale_a, 0, 0], [0, scale_b, 0],
[0, 0, scale_c]].


* **Parameters**


    * **scale_a** – Scaling factor for lattice direction a. Defaults to 1.


    * **scale_b** – Scaling factor for lattice direction b. Defaults to 1.


    * **scale_c** – Scaling factor for lattice direction c. Defaults to 1.



* **Returns**

    SupercellTransformation.



#### _property_ inverse()
NotImplementedError.


* **Type**

    Raises



#### _property_ is_one_to_many(_: boo_ )
False.


* **Type**

    Returns


## pymatgen.transformations.transformation_abc module

Defines an abstract base class contract for Transformation object.


### _class_ AbstractTransformation()
Bases: `MSONable`

Abstract transformation class.


#### _abc_impl(_ = <_abc._abc_data object_ )

#### _abstract_ apply_transformation(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))
Applies the transformation to a structure. Depending on whether a
transformation is one-to-many, there may be an option to return a
ranked list of structures.


* **Parameters**


    * **structure** – input structure


    * **return_ranked_list** (*bool** | **int**, **optional*) – If return_ranked_list is int, that number of structures

    is returned. If False, only the single lowest energy structure is returned. Defaults to False.




* **Returns**

    depending on returned_ranked list, either a transformed structure
    or
    a list of dictionaries, where each dictionary is of the form
    {‘structure’ = …. , ‘other_arguments’}
    the key ‘transformation’ is reserved for the transformation that
    was actually applied to the structure.
    This transformation is parsed by the alchemy classes for generating
    a more specific transformation history. Any other information will
    be stored in the transformation_parameters dictionary in the
    transmuted structure class.



#### _abstract property_ inverse(_: AbstractTransformation | Non_ )
Returns the inverse transformation if available.
Otherwise, should return None.


#### _abstract property_ is_one_to_many(_: boo_ )
Determines if a Transformation is a one-to-many transformation. If a
Transformation is a one-to-many transformation, the
apply_transformation method should have a keyword arg
“return_ranked_list” which allows for the transformed structures to be
returned as a ranked list.


#### _property_ use_multiprocessing(_: boo_ )
Indicates whether the transformation can be applied by a
subprocessing pool. This should be overridden to return True for
transformations that the transmuter can parallelize.