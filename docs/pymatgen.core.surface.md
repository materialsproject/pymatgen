---
layout: default
title: pymatgen.core.surface.md
nav_exclude: true
---

# pymatgen.core.surface module

This module implements representations of slabs and surfaces + algorithms for generating them.

If you use this module, please consider citing the following work:

>
> 1. Tran, Z. Xu, B. Radhakrishnan, D. Winston, W. Sun, K. A. Persson,

> S. P. Ong, “Surface Energies of Elemental Crystals”, Scientific Data,
> 2016, 3:160080, doi: 10.1038/sdata.2016.80.

> Sun, W.; Ceder, G. Efficient creation and convergence of surface slabs,
> Surface Science, 2013, 617, 53-59, doi:10.1016/j.susc.2013.05.016.


### _class_ pymatgen.core.surface.ReconstructionGenerator(initial_structure, min_slab_size, min_vacuum_size, reconstruction_name)
Bases: `object`

This class takes in a pre-defined dictionary specifying the parameters
need to build a reconstructed slab such as the SlabGenerator parameters,
transformation matrix, sites to remove/add and slab/vacuum size. It will
then use the formatted instructions provided by the dictionary to build
the desired reconstructed slab from the initial structure.


#### slabgen_params()
Parameters for the SlabGenerator

Todo:
- Right now there is no way to specify what atom is being

> added. In the future, use basis sets?

Generates reconstructed slabs from a set of instructions

    specified by a dictionary or json file.


* **Parameters**


    * **initial_structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Initial input structure. Note
    that to ensure that the miller indices correspond to usual
    crystallographic definitions, you should supply a conventional
    unit cell structure.


    * **min_slab_size** (*float*) – In Angstroms


    * **min_vacuum_size** (*float*) – In Angstroms


    * **reconstruction_name** (*str*) – Name of the dict containing the instructions
    for building a reconstructed slab. The dictionary can contain
    any item the creator deems relevant, however any instructions
    archived in pymatgen for public use needs to contain the
    following keys and items to ensure compatibility with the
    ReconstructionGenerator:

    > ”name” (str): A descriptive name for the type of

    >     reconstruction. Typically the name will have the type
    >     of structure the reconstruction is for, the Miller
    >     index, and Wood’s notation along with anything to
    >     describe the reconstruction: e.g.:
    >     “fcc_110_missing_row_1x2”

    > ”description” (str): A longer description of your

    >     reconstruction. This is to help future contributors who
    >     want to add other types of reconstructions to the
    >     archive on pymatgen to check if the reconstruction
    >     already exists. Please read the descriptions carefully
    >     before adding a new type of reconstruction to ensure it
    >     is not in the archive yet.

    > ”reference” (str): Optional reference to where the

    >     reconstruction was taken from or first observed.

    > ”spacegroup” (dict): e.g. {“symbol”: “Fm-3m”, “number”: 225}

    >     Indicates what kind of structure is this reconstruction.

    > ”miller_index” ([h,k,l]): Miller index of your reconstruction
    > “Woods_notation” (str): For a reconstruction, the a and b

    > > lattice may change to accommodate the symmetry of the
    > > reconstruction. This notation indicates the change in
    > > the vectors relative to the primitive (p) or
    > > conventional (c) slab cell. E.g. p(2x1):

    > > Wood, E. A. (1964). Vocabulary of surface
    > > crystallography. Journal of Applied Physics, 35(4),
    > > 1306-1312.

    > ”transformation_matrix” (numpy array): A 3x3 matrix to

    >     transform the slab. Only the a and b lattice vectors
    >     should change while the c vector remains the same.

    > ”SlabGenerator_parameters” (dict): A dictionary containing

    >     the parameters for the SlabGenerator class excluding the
    >     miller_index, min_slab_size and min_vac_size as the
    >     Miller index is already specified and the min_slab_size
    >     and min_vac_size can be changed regardless of what type
    >     of reconstruction is used. Having a consistent set of
    >     SlabGenerator parameters allows for the instructions to
    >     be reused to consistently build a reconstructed slab.

    > ”points_to_remove” (list of coords): A list of sites to

    >     remove where the first two indices are fraction (in a
    >     and b) and the third index is in units of 1/d (in c).

    > ”points_to_add” (list of frac_coords): A list of sites to add

    >     where the first two indices are fraction (in a an b) and
    >     the third index is in units of 1/d (in c).

    > ”base_reconstruction” (dict): Option to base a reconstruction on

    >     an existing reconstruction model also exists to easily build
    >     the instructions without repeating previous work. E.g. the
    >     alpha reconstruction of halites is based on the octopolar
    >     reconstruction but with the topmost atom removed. The dictionary
    >     for the alpha reconstruction would therefore contain the item
    >     “reconstruction_base”: “halite_111_octopolar_2x2”, and
    >     additional sites for “points_to_remove” and “points_to_add”
    >     can be added to modify this reconstruction.

    > For “points_to_remove” and “points_to_add”, the third index for

    >     the c vector is in units of 1/d where d is the spacing
    >     between atoms along hkl (the c vector) and is relative to
    >     the topmost site in the unreconstructed slab. e.g. a point
    >     of [0.5, 0.25, 1] corresponds to the 0.5 frac_coord of a,
    >     0.25 frac_coord of b and a distance of 1 atomic layer above
    >     the topmost site. [0.5, 0.25, -0.5] where the third index
    >     corresponds to a point half a atomic layer below the topmost
    >     site. [0.5, 0.25, 0] corresponds to a point in the same
    >     position along c as the topmost site. This is done because
    >     while the primitive units of a and b will remain constant,
    >     the user can vary the length of the c direction by changing
    >     the slab layer or the vacuum layer.



    * **NOTE** – THE DICTIONARY SHOULD ONLY CONTAIN “points_to_remove” AND


    * **ReconstructionGenerator** (*"points_to_add" FOR THE TOP SURFACE. THE*) –


    * **WITH** (*WILL MODIFY THE BOTTOM SURFACE ACCORDINGLY TO RETURN A SLAB*) –


    * **SURFACES.** (*EQUIVALENT*) –



#### build_slabs()
Builds the reconstructed slab by:


    1. Obtaining the unreconstructed slab using the specified
    parameters for the SlabGenerator.


    2. Applying the appropriate lattice transformation in the
    a and b lattice vectors.


    3. Remove any specified sites from both surfaces.


    4. Add any specified sites to both surfaces.


* **Returns**

    The reconstructed slab.



* **Return type**

    (Slab)



#### get_unreconstructed_slabs()
Generates the unreconstructed or pristine super slab.


### _class_ pymatgen.core.surface.Slab(lattice, species, coords, miller_index, oriented_unit_cell, shift, scale_factor, reorient_lattice=True, validate_proximity=False, to_unit_cell=False, reconstruction=None, coords_are_cartesian=False, site_properties=None, energy=None)
Bases: [`Structure`](pymatgen.core.structure.md#pymatgen.core.structure.Structure)

Subclass of Structure representing a Slab. Implements additional
attributes pertaining to slabs, but the init method does not
actually implement any algorithm that creates a slab. This is a
DUMMY class who’s init method only holds information about the
slab. Also has additional methods that returns other information
about a slab such as the surface area, normal, and atom adsorption.

Note that all Slabs have the surface normal oriented perpendicular to the a
and b lattice vectors. This means the lattice vectors a and b are in the
surface plane and the c vector is out of the surface plane (though not
necessarily perpendicular to the surface).


#### miller_index()
Miller index of plane parallel to surface.


#### scale_factor()
Final computed scale factor that brings the parent cell to the
surface cell.


#### shift()
The shift value in Angstrom that indicates how much this
slab has been shifted.

Makes a Slab structure, a structure object with additional information
and methods pertaining to slabs.


* **Parameters**


    * **lattice** (*Lattice/3x3 array*) – The lattice, either as a
    [`pymatgen.core.lattice.Lattice`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice) or
    simply as any 2D array. Each row should correspond to a lattice
    vector. E.g., [[10,0,0], [20,10,0], [0,0,30]] specifies a
    lattice with lattice vectors [10,0,0], [20,10,0] and [0,0,30].


    * **species** (*[*[*Species*](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Species)*]*) – Sequence of species on each site. Can take in
    flexible input, including:


        1. A sequence of element / species specified either as string
    symbols, e.g. [“Li”, “Fe2+”, “P”, …] or atomic numbers,
    e.g., (3, 56, …) or actual Element or Species objects.


        2. List of dict of elements/species and occupancies, e.g.,
    [{“Fe” : 0.5, “Mn”:0.5}, …]. This allows the setup of
    disordered structures.



    * **coords** (*Nx3 array*) – list of fractional/cartesian coordinates of
    each species.


    * **miller_index** (*[**h**, **k**, **l**]*) – Miller index of plane parallel to
    surface. Note that this is referenced to the input structure. If
    you need this to be based on the conventional cell,
    you should supply the conventional structure.


    * **oriented_unit_cell** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – The oriented_unit_cell from which
    this Slab is created (by scaling in the c-direction).


    * **shift** (*float*) – The shift in the c-direction applied to get the
    termination.


    * **scale_factor** (*np.ndarray*) – scale_factor Final computed scale factor
    that brings the parent cell to the surface cell.


    * **reorient_lattice** (*bool*) – reorients the lattice parameters such that
    the c direction is along the z axis.


    * **validate_proximity** (*bool*) – Whether to check if there are sites
    that are less than 0.01 Ang apart. Defaults to False.


    * **reconstruction** (*str*) – Type of reconstruction. Defaults to None if
    the slab is not reconstructed.


    * **to_unit_cell** (*bool*) – Translates fractional coordinates into the unit cell. Defaults to False.


    * **coords_are_cartesian** (*bool*) – Set to True if you are providing
    coordinates in Cartesian coordinates. Defaults to False.


    * **site_properties** (*dict*) – Properties associated with the sites as a
    dict of sequences, e.g., {“magmom”:[5,5,5,5]}. The sequences
    have to be the same length as the atomic species and
    fractional_coords. Defaults to None for no properties.


    * **energy** (*float*) – A value for the energy.



#### add_adsorbate_atom(indices, specie, distance)
Gets the structure of single atom adsorption.
slab structure from the Slab class(in [0, 0, 1]).


* **Parameters**


    * **indices** (*[**int**]*) – Indices of sites on which to put the absorbate.
    Absorbed atom will be displaced relative to the center of
    these sites.


    * **specie** (*Species/Element/str*) – adsorbed atom species


    * **distance** (*float*) – between centers of the adsorbed atom and the
    given site in Angstroms.



#### as_dict()

* **Returns**

    MSONable dict



#### _property_ center_of_mass()
Calculates the center of mass of the slab.


#### copy(site_properties=None, sanitize=False)
Convenience method to get a copy of the structure, with options to add
site properties.


* **Parameters**


    * **site_properties** (*dict*) – Properties to add or override. The
    properties are specified in the same way as the constructor,
    i.e., as a dict of the form {property: [values]}. The
    properties should be in the order of the *original* structure
    if you are performing sanitization.


    * **sanitize** (*bool*) – If True, this method will return a sanitized
    structure. Sanitization performs a few things: (i) The sites are
    sorted by electronegativity, (ii) a LLL lattice reduction is
    carried out to obtain a relatively orthogonalized cell,
    (iii) all fractional coords for sites are mapped into the
    unit cell.



* **Returns**

    A copy of the Structure, with optionally new site_properties and
    optionally sanitized.



#### _property_ dipole()
Calculates the dipole of the Slab in the direction of the surface
normal. Note that the Slab must be oxidation state-decorated for this
to work properly. Otherwise, the Slab will always have a dipole of 0.


#### _classmethod_ from_dict(d)

* **Parameters**

    **d** – dict



* **Returns**

    Creates slab from dict.



#### get_orthogonal_c_slab()
This method returns a Slab where the normal (c lattice vector) is
“forced” to be exactly orthogonal to the surface a and b lattice
vectors. **Note that this breaks inherent symmetries in the slab.**
It should be pointed out that orthogonality is not required to get good
surface energies, but it can be useful in cases where the slabs are
subsequently used for postprocessing of some kind, e.g. generating
GBs or interfaces.


#### get_sorted_structure(key=None, reverse=False)
Get a sorted copy of the structure. The parameters have the same
meaning as in list.sort. By default, sites are sorted by the
electronegativity of the species. Note that Slab has to override this
because of the different __init__ args.


* **Parameters**


    * **key** – Specifies a function of one argument that is used to extract
    a comparison key from each list element: key=str.lower. The
    default value is None (compare the elements directly).


    * **reverse** (*bool*) – If set to True, then the list elements are sorted
    as if each comparison were reversed.



#### get_surface_sites(tag=False)
Returns the surface sites and their indices in a dictionary. The
oriented unit cell of the slab will determine the coordination number
of a typical site. We use VoronoiNN to determine the
coordination number of bulk sites and slab sites. Due to the
pathological error resulting from some surface sites in the
VoronoiNN, we assume any site that has this error is a surface
site as well. This will work for elemental systems only for now. Useful
for analysis involving broken bonds and for finding adsorption sites.


* **Parameters**

    **tag** (*bool*) – Option to adds site attribute “is_surfsite” (bool)
    to all sites of slab. Defaults to False



* **Returns**

    A dictionary grouping sites on top and bottom of the slab
    together.
    {“top”: [sites with indices], “bottom”: [sites with indices}



#### get_symmetric_site(point, cartesian=False)
This method uses symmetry operations to find equivalent sites on

    both sides of the slab. Works mainly for slabs with Laue
    symmetry. This is useful for retaining the non-polar and
    symmetric properties of a slab when creating adsorbed
    structures or symmetric reconstructions.

Arg:

    point: Fractional coordinate.


* **Returns**

    Fractional coordinate. A point equivalent to the

        parameter point, but on the other side of the slab




* **Return type**

    point



#### get_tasker2_slabs(tol: float = 0.01, same_species_only=True)
Get a list of slabs that have been Tasker 2 corrected.


* **Parameters**


    * **tol** (*float*) – Tolerance to determine if atoms are within same plane.
    This is a fractional tolerance, not an absolute one.


    * **same_species_only** (*bool*) – If True, only that are of the exact same
    species as the atom at the outermost surface are considered for
    moving. Otherwise, all atoms regardless of species that is
    within tol are considered for moving. Default is True (usually
    the desired behavior).



* **Returns**

    ([Slab]) List of tasker 2 corrected slabs.



#### is_polar(tol_dipole_per_unit_area=0.001)
Checks whether the surface is polar by computing the dipole per unit
area. Note that the Slab must be oxidation state-decorated for this
to work properly. Otherwise, the Slab will always be non-polar.


* **Parameters**

    **tol_dipole_per_unit_area** (*float*) – A tolerance. If the dipole
    magnitude per unit area is less than this value, the Slab is
    considered non-polar. Defaults to 1e-3, which is usually
    pretty good. Normalized dipole per unit area is used as it is
    more reliable than using the total, which tends to be larger for
    slabs with larger surface areas.



#### is_symmetric(symprec: float = 0.1)
Checks if surfaces are symmetric, i.e., contains inversion, mirror on (hkl) plane,

    or screw axis (rotation and translation) about [hkl].


* **Parameters**

    **symprec** (*float*) – Symmetry precision used for SpaceGroup analyzer.



* **Returns**

    Whether surfaces are symmetric.



* **Return type**

    bool



#### _property_ normal()
Calculates the surface normal vector of the slab.


#### _property_ surface_area()
Calculates the surface area of the slab.


#### symmetrically_add_atom(specie, point, coords_are_cartesian=False)
Class method for adding a site at a specified point in a slab.

    Will add the corresponding site on the other side of the
    slab to maintain equivalent surfaces.

Arg:

    specie (str): The specie to add
    point (coords): The coordinate of the site in the slab to add.
    coords_are_cartesian (bool): Is the point in Cartesian coordinates


* **Returns**

    The modified slab



* **Return type**

    (Slab)



#### symmetrically_remove_atoms(indices)
Class method for removing sites corresponding to a list of indices.

    Will remove the corresponding site on the other side of the
    slab to maintain equivalent surfaces.

Arg:

    indices ([indices]): The indices of the sites

        in the slab to remove.


### _class_ pymatgen.core.surface.SlabGenerator(initial_structure, miller_index, min_slab_size, min_vacuum_size, lll_reduce=False, center_slab=False, in_unit_planes=False, primitive=True, max_normal_search=None, reorient_lattice=True)
Bases: `object`

This class generates different slabs using shift values determined by where
a unique termination can be found along with other criteria such as where a
termination doesn’t break a polyhedral bond. The shift value then indicates
where the slab layer will begin and terminate in the slab-vacuum system.


#### oriented_unit_cell()
A unit cell of the parent structure with the miller
index of plane parallel to surface


#### parent()
Parent structure from which Slab was derived.


#### lll_reduce()
Whether or not the slabs will be orthogonalized


#### center_slab()
Whether or not the slabs will be centered between
the vacuum layer


#### slab_scale_factor()
Final computed scale factor that brings the parent cell to the
surface cell.


#### miller_index()
Miller index of plane parallel to surface.


#### min_slab_size()
Minimum size in angstroms of layers containing atoms


#### min_vac_size()
Minimize size in angstroms of layers containing vacuum

Calculates the slab scale factor and uses it to generate a unit cell
of the initial structure that has been oriented by its miller index.
Also stores the initial information needed later on to generate a slab.


* **Parameters**


    * **initial_structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Initial input structure. Note that to
    ensure that the miller indices correspond to usual
    crystallographic definitions, you should supply a conventional
    unit cell structure.


    * **miller_index** (*[**h**, **k**, **l**]*) – Miller index of plane parallel to
    surface. Note that this is referenced to the input structure. If
    you need this to be based on the conventional cell,
    you should supply the conventional structure.


    * **min_slab_size** (*float*) – In Angstroms or number of hkl planes


    * **min_vacuum_size** (*float*) – In Angstroms or number of hkl planes


    * **lll_reduce** (*bool*) – Whether to perform an LLL reduction on the
    eventual structure.


    * **center_slab** (*bool*) – Whether to center the slab in the cell with
    equal vacuum spacing from the top and bottom.


    * **in_unit_planes** (*bool*) – Whether to set min_slab_size and min_vac_size
    in units of hkl planes (True) or Angstrom (False/default).
    Setting in units of planes is useful for ensuring some slabs
    have a certain n_layer of atoms. e.g. for Cs (100), a 10 Ang
    slab will result in a slab with only 2 layer of atoms, whereas
    Fe (100) will have more layer of atoms. By using units of hkl
    planes instead, we ensure both slabs
    have the same number of atoms. The slab thickness will be in
    min_slab_size/math.ceil(self._proj_height/dhkl)
    multiples of oriented unit cells.


    * **primitive** (*bool*) – Whether to reduce any generated slabs to a
    primitive cell (this does **not** mean the slab is generated
    from a primitive cell, it simply means that after slab
    generation, we attempt to find shorter lattice vectors,
    which lead to less surface area and smaller cells).


    * **max_normal_search** (*int*) – If set to a positive integer, the code will
    conduct a search for a normal lattice vector that is as
    perpendicular to the surface as possible by considering
    multiples linear combinations of lattice vectors up to
    max_normal_search. This has no bearing on surface energies,
    but may be useful as a preliminary step to generating slabs
    for absorption and other sizes. It is typical that this will
    not be the smallest possible cell for simulation. Normality
    is not guaranteed, but the oriented cell will have the c
    vector as normal as possible (within the search range) to the
    surface. A value of up to the max absolute Miller index is
    usually sufficient.


    * **reorient_lattice** (*bool*) – reorients the lattice parameters such that
    the c direction is the third vector of the lattice matrix



#### get_slab(shift=0, tol: float = 0.1, energy=None)
This method takes in shift value for the c lattice direction and
generates a slab based on the given shift. You should rarely use this
method. Instead, it is used by other generation algorithms to obtain
all slabs.

Arg:

    shift (float): A shift value in Angstrom that determines how much a

        slab should be shifted.

    tol (float): Tolerance to determine primitive cell.
    energy (float): An energy to assign to the slab.


* **Returns**

    (Slab) A Slab object with a particular shifted oriented unit cell.



#### get_slabs(bonds=None, ftol=0.1, tol=0.1, max_broken_bonds=0, symmetrize=False, repair=False)
This method returns a list of slabs that are generated using the list of
shift values from the method, _calculate_possible_shifts(). Before the
shifts are used to create the slabs however, if the user decides to take
into account whether or not a termination will break any polyhedral
structure (bonds is not None), this method will filter out any shift
values that do so.


* **Parameters**


    * **bonds** (*{**(**specie1**, **specie2*) – max_bond_dist}: bonds are
    specified as a dict of tuples: float of specie1, specie2
    and the max bonding distance. For example, PO4 groups may be
    defined as {(“P”, “O”): 3}.


    * **tol** (*float*) – General tolerance parameter for getting primitive
    cells and matching structures


    * **ftol** (*float*) – Threshold parameter in fcluster in order to check
    if two atoms are lying on the same plane. Default thresh set
    to 0.1 Angstrom in the direction of the surface normal.


    * **max_broken_bonds** (*int*) – Maximum number of allowable broken bonds
    for the slab. Use this to limit # of slabs (some structures
    may have a lot of slabs). Defaults to zero, which means no
    defined bonds must be broken.


    * **symmetrize** (*bool*) – Whether or not to ensure the surfaces of the
    slabs are equivalent.


    * **repair** (*bool*) – Whether to repair terminations with broken bonds
    or just omit them. Set to False as repairing terminations can
    lead to many possible slabs as oppose to just omitting them.



* **Returns**

    ([Slab]) List of all possible terminations of a particular surface.
    Slabs are sorted by the # of bonds broken.



#### move_to_other_side(init_slab, index_of_sites)
This method will Move a set of sites to the
other side of the slab (opposite surface).

Arg:

    init_slab (structure): A structure object representing a slab.
    index_of_sites (list of ints): The list of indices representing

    > the sites we want to move to the other side.


* **Returns**

    (Slab) A Slab object with a particular shifted oriented unit cell.



#### nonstoichiometric_symmetrized_slab(init_slab)
This method checks whether or not the two surfaces of the slab are
equivalent. If the point group of the slab has an inversion symmetry (
ie. belong to one of the Laue groups), then it is assumed that the
surfaces should be equivalent. Otherwise, sites at the bottom of the
slab will be removed until the slab is symmetric. Note the removal of sites
can destroy the stoichiometry of the slab. For non-elemental
structures, the chemical potential will be needed to calculate surface energy.

Arg:

    init_slab (Structure): A single slab structure


* **Returns**

    A symmetrized Slab object.



* **Return type**

    Slab (structure)



#### repair_broken_bonds(slab, bonds)
This method will find undercoordinated atoms due to slab
cleaving specified by the bonds parameter and move them
to the other surface to make sure the bond is kept intact.
In a future release of surface.py, the ghost_sites will be
used to tell us how the repair bonds should look like.

Arg:

    slab (structure): A structure object representing a slab.
    bonds ({(specie1, specie2): max_bond_dist}: bonds are

    > specified as a dict of tuples: float of specie1, specie2
    > and the max bonding distance. For example, PO4 groups may be
    > defined as {(“P”, “O”): 3}.


* **Returns**

    (Slab) A Slab object with a particular shifted oriented unit cell.



### pymatgen.core.surface.center_slab(slab)
The goal here is to ensure the center of the slab region

    is centered close to c=0.5. This makes it easier to
    find the surface sites and apply operations like doping.

There are three cases where the slab in not centered:

1. The slab region is completely between two vacuums in the
box but not necessarily centered. We simply shift the
slab by the difference in its center of mass and 0.5
along the c direction.

2. The slab completely spills outside the box from the bottom
and into the top. This makes it incredibly difficult to
locate surface sites. We iterate through all sites that
spill over (z>c) and shift all sites such that this specific
site is now on the other side. Repeat for all sites with z>c.

3. This is a simpler case of scenario 2. Either the top or bottom
slab sites are at c=0 or c=1. Treat as scenario 2.


* **Parameters**

    **slab** (*Slab*) – Slab structure to center



* **Returns**

    Returns a centered slab structure



### pymatgen.core.surface.generate_all_slabs(structure, max_index, min_slab_size, min_vacuum_size, bonds=None, tol=0.1, ftol=0.1, max_broken_bonds=0, lll_reduce=False, center_slab=False, primitive=True, max_normal_search=None, symmetrize=False, repair=False, include_reconstructions=False, in_unit_planes=False)
A function that finds all different slabs up to a certain miller index.
Slabs oriented under certain Miller indices that are equivalent to other
slabs in other Miller indices are filtered out using symmetry operations
to get rid of any repetitive slabs. For example, under symmetry operations,
CsCl has equivalent slabs in the (0,0,1), (0,1,0), and (1,0,0) direction.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Initial input structure. Note that to
    ensure that the miller indices correspond to usual
    crystallographic definitions, you should supply a conventional
    unit cell structure.


    * **max_index** (*int*) – The maximum Miller index to go up to.


    * **min_slab_size** (*float*) – In Angstroms


    * **min_vacuum_size** (*float*) – In Angstroms


    * **bonds** (*{**(**specie1**, **specie2*) – max_bond_dist}: bonds are
    specified as a dict of tuples: float of specie1, specie2
    and the max bonding distance. For example, PO4 groups may be
    defined as {(“P”, “O”): 3}.


    * **tol** (*float*) – General tolerance parameter for getting primitive
    cells and matching structures


    * **ftol** (*float*) – Threshold parameter in fcluster in order to check
    if two atoms are lying on the same plane. Default thresh set
    to 0.1 Angstrom in the direction of the surface normal.


    * **max_broken_bonds** (*int*) – Maximum number of allowable broken bonds
    for the slab. Use this to limit # of slabs (some structures
    may have a lot of slabs). Defaults to zero, which means no
    defined bonds must be broken.


    * **lll_reduce** (*bool*) – Whether to perform an LLL reduction on the
    eventual structure.


    * **center_slab** (*bool*) – Whether to center the slab in the cell with
    equal vacuum spacing from the top and bottom.


    * **primitive** (*bool*) – Whether to reduce any generated slabs to a
    primitive cell (this does **not** mean the slab is generated
    from a primitive cell, it simply means that after slab
    generation, we attempt to find shorter lattice vectors,
    which lead to less surface area and smaller cells).


    * **max_normal_search** (*int*) – If set to a positive integer, the code will
    conduct a search for a normal lattice vector that is as
    perpendicular to the surface as possible by considering
    multiples linear combinations of lattice vectors up to
    max_normal_search. This has no bearing on surface energies,
    but may be useful as a preliminary step to generating slabs
    for absorption and other sizes. It is typical that this will
    not be the smallest possible cell for simulation. Normality
    is not guaranteed, but the oriented cell will have the c
    vector as normal as possible (within the search range) to the
    surface. A value of up to the max absolute Miller index is
    usually sufficient.


    * **symmetrize** (*bool*) – Whether or not to ensure the surfaces of the
    slabs are equivalent.


    * **repair** (*bool*) – Whether to repair terminations with broken bonds
    or just omit them


    * **include_reconstructions** (*bool*) – Whether to include reconstructed
    slabs available in the reconstructions_archive.json file. Defaults to False.


    * **in_unit_planes** (*bool*) – Whether to generate slabs in units of the primitive
    cell’s c lattice vector. This is useful for generating slabs with
    a specific number of layers, as the number of layers will be
    independent of the Miller index. Defaults to False.


    * **in_unit_planes** – Whether to set min_slab_size and min_vac_size
    in units of hkl planes (True) or Angstrom (False, the default). Setting in
    units of planes is useful for ensuring some slabs have a certain n_layer of
    atoms. e.g. for Cs (100), a 10 Ang slab will result in a slab with only 2
    layer of atoms, whereas Fe (100) will have more layer of atoms. By using units
    of hkl planes instead, we ensure both slabs have the same number of atoms. The
    slab thickness will be in min_slab_size/math.ceil(self._proj_height/dhkl)
    multiples of oriented unit cells.



### pymatgen.core.surface.get_d(slab)
Determine the distance of space between
each layer of atoms along c.


### pymatgen.core.surface.get_slab_regions(slab, blength=3.5)
Function to get the ranges of the slab regions. Useful for discerning where
the slab ends and vacuum begins if the slab is not fully within the cell
:param slab: Structure object modelling the surface
:type slab: Structure
:param blength: The bondlength between atoms. You generally

> want this value to be larger than the actual bondlengths in
> order to find atoms that are part of the slab.


### pymatgen.core.surface.get_symmetrically_distinct_miller_indices(structure, max_index, return_hkil=False)
Returns all symmetrically distinct indices below a certain max-index for
a given structure. Analysis is based on the symmetry of the reciprocal
lattice of the structure.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – input structure.


    * **max_index** (*int*) – The maximum index. For example, a max_index of 1
    means that (100), (110), and (111) are returned for the cubic
    structure. All other indices are equivalent to one of these.


    * **return_hkil** (*bool*) – If true, return hkil form of Miller
    index for hexagonal systems, otherwise return hkl



### pymatgen.core.surface.get_symmetrically_equivalent_miller_indices(structure, miller_index, return_hkil=True, system: Literal['triclinic', 'monoclinic', 'orthorhombic', 'tetragonal', 'trigonal', 'hexagonal', 'cubic'] | None = None)
Returns all symmetrically equivalent indices for a given structure. Analysis
is based on the symmetry of the reciprocal lattice of the structure.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Structure to analyze


    * **miller_index** (*tuple*) – Designates the family of Miller indices
    to find. Can be hkl or hkil for hexagonal systems


    * **return_hkil** (*bool*) – If true, return hkil form of Miller
    index for hexagonal systems, otherwise return hkl


    * **system** – If known, specify the crystal system of the structure
    so that it does not need to be re-calculated.



### pymatgen.core.surface.hkl_transformation(transf, miller_index)
Returns the Miller index from setting
A to B using a transformation matrix
:param transf: The transformation matrix

> that transforms a lattice of A to B


* **Parameters**

    **miller_index** (*[**h**, **k**, **l**]*) – Miller index to transform to setting B.



### pymatgen.core.surface.is_already_analyzed(miller_index: tuple, miller_list: list, symm_ops: list)
Helper function to check if a given Miller index is
part of the family of indices of any index in a list.


* **Parameters**


    * **miller_index** (*tuple*) – The Miller index to analyze


    * **miller_list** (*list*) – List of Miller indices. If the given
    Miller index belongs in the same family as any of the
    indices in this list, return True, else return False


    * **symm_ops** (*list*) – Symmetry operations of a
    lattice, used to define family of indices



### pymatgen.core.surface.miller_index_from_sites(lattice, coords, coords_are_cartesian=True, round_dp=4, verbose=True)
Get the Miller index of a plane from a list of site coordinates.

A minimum of 3 sets of coordinates are required. If more than 3 sets of
coordinates are given, the best plane that minimises the distance to all
points will be calculated.


* **Parameters**


    * **lattice** (*list** or *[*Lattice*](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice)) – A 3x3 lattice matrix or Lattice object (for
    example obtained from Structure.lattice).


    * **coords** (*iterable*) – A list or numpy array of coordinates. Can be
    Cartesian or fractional coordinates. If more than three sets of
    coordinates are provided, the best plane that minimises the
    distance to all sites will be calculated.


    * **coords_are_cartesian** (*bool**, **optional*) – Whether the coordinates are
    in Cartesian space. If using fractional coordinates set to False.


    * **round_dp** (*int**, **optional*) – The number of decimal places to round the
    miller index to.


    * **verbose** (*bool**, **optional*) – Whether to print warnings.



* **Returns**

    The Miller index.



* **Return type**

    (tuple)