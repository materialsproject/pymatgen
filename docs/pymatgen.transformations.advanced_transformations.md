---
layout: default
title: pymatgen.transformations.advanced_transformations.md
nav_exclude: true
---

# pymatgen.transformations.advanced_transformations module

This module implements more advanced transformations.


### _class_ pymatgen.transformations.advanced_transformations.AddAdsorbateTransformation(adsorbate, selective_dynamics=False, height=0.9, mi_vec=None, repeat=None, min_lw=5.0, translate=True, reorient=True, find_args=None)
Bases: [`AbstractTransformation`](pymatgen.transformations.transformation_abc.md#pymatgen.transformations.transformation_abc.AbstractTransformation)

Create absorbate structures.

Use AdsorbateSiteFinder to add an absorbate to a slab.


* **Parameters**


    * **adsorbate** ([*Molecule*](pymatgen.core.structure.md#pymatgen.core.structure.Molecule)) – molecule to add as adsorbate


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



#### apply_transformation(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), return_ranked_list: bool | int = False)

* **Parameters**


    * **structure** – Must be a Slab structure


    * **return_ranked_list** (*bool** | **int**, **optional*) – If return_ranked_list is int, that number of structures.

    is returned. If False, only the single lowest energy structure is returned. Defaults to False.



Returns: Slab with adsorbate


#### _property_ inverse()
None.


* **Type**

    Returns



#### _property_ is_one_to_many(_: boo_ )
True.


* **Type**

    Returns



### _class_ pymatgen.transformations.advanced_transformations.ChargeBalanceTransformation(charge_balance_sp)
Bases: [`AbstractTransformation`](pymatgen.transformations.transformation_abc.md#pymatgen.transformations.transformation_abc.AbstractTransformation)

This is a transformation that disorders a structure to make it charge
balanced, given an oxidation state-decorated structure.


* **Parameters**

    **charge_balance_sp** – specie to add or remove. Currently only removal
    is supported.



#### apply_transformation(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))
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



### _class_ pymatgen.transformations.advanced_transformations.CubicSupercellTransformation(min_atoms: int | None = None, max_atoms: int | None = None, min_length: float = 15.0, force_diagonal: bool = False, force_90_degrees: bool = False, angle_tolerance: float = 0.001)
Bases: [`AbstractTransformation`](pymatgen.transformations.transformation_abc.md#pymatgen.transformations.transformation_abc.AbstractTransformation)

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



#### apply_transformation(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))
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
Returns:
None.


#### _property_ is_one_to_many(_: boo_ )
Returns:
False.


### _class_ pymatgen.transformations.advanced_transformations.DisorderOrderedTransformation(max_sites_to_merge=2)
Bases: [`AbstractTransformation`](pymatgen.transformations.transformation_abc.md#pymatgen.transformations.transformation_abc.AbstractTransformation)

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



#### apply_transformation(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), return_ranked_list: bool | int = False)

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



### _class_ pymatgen.transformations.advanced_transformations.DopingTransformation(dopant, ionic_radius_tol=inf, min_length=10, alio_tol=0, codopant=False, max_structures_per_enum=100, allowed_doping_species=None, \*\*kwargs)
Bases: [`AbstractTransformation`](pymatgen.transformations.transformation_abc.md#pymatgen.transformations.transformation_abc.AbstractTransformation)

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


    * **\*\*kwargs** – Same keyword args as `EnumerateStructureTransformation`,
    i.e., min_cell_size, etc.



#### apply_transformation(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), return_ranked_list: bool | int = False)

* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Input structure to dope


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



### _class_ pymatgen.transformations.advanced_transformations.EnumerateStructureTransformation(min_cell_size: int = 1, max_cell_size: int = 1, symm_prec: float = 0.1, refine_structure: bool = False, enum_precision_parameter: float = 0.001, check_ordered_symmetry: bool = True, max_disordered_sites: int | None = None, sort_criteria: str | Callable = 'ewald', timeout: float | None = None, n_jobs: int = -1)
Bases: [`AbstractTransformation`](pymatgen.transformations.transformation_abc.md#pymatgen.transformations.transformation_abc.AbstractTransformation)

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



#### apply_transformation(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), return_ranked_list: bool | int = False)
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



### _class_ pymatgen.transformations.advanced_transformations.GrainBoundaryTransformation(rotation_axis, rotation_angle, expand_times=4, vacuum_thickness=0.0, ab_shift=None, normal=False, ratio=True, plane=None, max_search=20, tol_coi=1e-08, rm_ratio=0.7, quick_gen=False)
Bases: [`AbstractTransformation`](pymatgen.transformations.transformation_abc.md#pymatgen.transformations.transformation_abc.AbstractTransformation)

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



#### apply_transformation(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))
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



### _class_ pymatgen.transformations.advanced_transformations.MagOrderParameterConstraint(order_parameter, species_constraints=None, site_constraint_name=None, site_constraints=None)
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


### _class_ pymatgen.transformations.advanced_transformations.MagOrderingTransformation(mag_species_spin, order_parameter=0.5, energy_model=None, \*\*kwargs)
Bases: [`AbstractTransformation`](pymatgen.transformations.transformation_abc.md#pymatgen.transformations.transformation_abc.AbstractTransformation)

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
    `pymatgen.transformations.advanced_transformations.MagOrderParameterConstraint`
    to specify more complicated orderings, see documentation for
    MagOrderParameterConstraint more details on usage


    * **energy_model** – Energy model to rank the returned structures,
    see :mod: pymatgen.analysis.energy_models for more information (note
    that this is not necessarily a physical energy). By default, returned
    structures use SymmetryModel() which ranks structures from most
    symmetric to least.


    * **kwargs** – Additional kwargs that are passed to


`EnumerateStructureTransformation` such as min_cell_size etc.


#### apply_transformation(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), return_ranked_list: bool | int = False)
Apply MagOrderTransformation to an input structure.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Any ordered structure.


    * **return_ranked_list** (*bool** | **int**, **optional*) – If return_ranked_list is int, that number of structures
    is returned. If False, only the single lowest energy structure is returned. Defaults to False.



* **Raises**

    **ValueError** – On disordered structures.



* **Returns**

    Structure(s) after MagOrderTransformation.



* **Return type**

    [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure) | list[[Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure)]



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



### _class_ pymatgen.transformations.advanced_transformations.MonteCarloRattleTransformation(rattle_std: float, min_distance: float, seed: int | None = None, \*\*kwargs)
Bases: [`AbstractTransformation`](pymatgen.transformations.transformation_abc.md#pymatgen.transformations.transformation_abc.AbstractTransformation)

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



#### apply_transformation(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))
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



### _class_ pymatgen.transformations.advanced_transformations.MultipleSubstitutionTransformation(sp_to_replace, r_fraction, substitution_dict, charge_balance_species=None, order=True)
Bases: `object`

Performs multiple substitutions on a structure. For example, can do a
fractional replacement of Ge in LiGePS with a list of species, creating one
structure for each substitution. Ordering is done using a dummy element so
only one ordering must be done per substitution oxidation state. Charge
balancing of the structure is optionally performed.

**NOTE**: There are no checks to make sure that removal fractions are possible
and rounding may occur. Currently charge balancing only works for
removal of species.

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



#### apply_transformation(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), return_ranked_list: bool | int = False)
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



### _class_ pymatgen.transformations.advanced_transformations.SQSTransformation(scaling, cluster_size_and_shell=None, search_time=60, directory=None, instances=None, temperature=1, wr=1, wn=1, wd=0.5, tol=0.001, best_only=True, remove_duplicate_structures=True, reduction_algo='LLL')
Bases: [`AbstractTransformation`](pymatgen.transformations.transformation_abc.md#pymatgen.transformations.transformation_abc.AbstractTransformation)

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



#### apply_transformation(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), return_ranked_list: bool | int = False)
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



### _class_ pymatgen.transformations.advanced_transformations.SlabTransformation(miller_index, min_slab_size, min_vacuum_size, lll_reduce=False, center_slab=False, in_unit_planes=False, primitive=True, max_normal_search=None, shift=0, tol=0.1)
Bases: [`AbstractTransformation`](pymatgen.transformations.transformation_abc.md#pymatgen.transformations.transformation_abc.AbstractTransformation)

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



#### apply_transformation(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))
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



### _class_ pymatgen.transformations.advanced_transformations.SubstituteSurfaceSiteTransformation(atom, selective_dynamics=False, height=0.9, mi_vec=None, target_species=None, sub_both_sides=False, range_tol=0.01, dist_from_surf=0)
Bases: [`AbstractTransformation`](pymatgen.transformations.transformation_abc.md#pymatgen.transformations.transformation_abc.AbstractTransformation)

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



#### apply_transformation(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), return_ranked_list: bool | int = False)

* **Parameters**


    * **structure** – Must be a Slab structure


    * **return_ranked_list** (*bool** | **int**, **optional*) – If return_ranked_list is int, that number of structures.

    is returned. If False, only the single lowest energy structure is returned. Defaults to False.



Returns: Slab with sites substituted


#### _property_ inverse()
None.


* **Type**

    Returns



#### _property_ is_one_to_many(_: boo_ )
True.


* **Type**

    Returns



### _class_ pymatgen.transformations.advanced_transformations.SubstitutionPredictorTransformation(threshold=0.01, scale_volumes=True, \*\*kwargs)
Bases: [`AbstractTransformation`](pymatgen.transformations.transformation_abc.md#pymatgen.transformations.transformation_abc.AbstractTransformation)

This transformation takes a structure and uses the structure
prediction module to find likely site substitutions.


* **Parameters**


    * **threshold** – Threshold for substitution.


    * **scale_volumes** – Whether to scale volumes after substitution.


    * **\*\*kwargs** – Args for SubstitutionProbability class lambda_table, alpha.



#### apply_transformation(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), return_ranked_list: bool | int = False)
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



### _class_ pymatgen.transformations.advanced_transformations.SuperTransformation(transformations, nstructures_per_trans=1)
Bases: [`AbstractTransformation`](pymatgen.transformations.transformation_abc.md#pymatgen.transformations.transformation_abc.AbstractTransformation)

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



#### apply_transformation(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), return_ranked_list: bool | int = False)
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