---
layout: default
title: pymatgen.analysis.ferroelectricity.md
nav_exclude: true
---

1. TOC
{:toc}

# pymatgen.analysis.ferroelectricity package

Package for analyzing ferroelectricity.


## pymatgen.analysis.ferroelectricity.polarization module

This module contains classes useful for analyzing ferroelectric candidates.
The Polarization class can recover the spontaneous polarization using
multiple calculations along a nonpolar to polar ferroelectric distortion.
The EnergyTrend class is useful for assessing the trend in energy across
the distortion.

See Nicola Spaldin’s “A beginner’s guide to the modern theory of polarization”
([https://arxiv.org/abs/1202.1831](https://arxiv.org/abs/1202.1831)) for an introduction to crystal polarization.

VASP reports dipole moment values (used to derive polarization) along Cartesian
directions (see pead.F around line 970 in the VASP source to confirm this).
However, it is most convenient to perform the adjustments necessary to recover
a same branch polarization by expressing the polarization along lattice directions.
For this reason, calc_ionic calculates ionic contributions to the polarization
along lattice directions. We provide the means to convert Cartesian direction
polarizations to lattice direction polarizations in the Polarization class.

We recommend using our calc_ionic function for calculating the ionic
polarization rather than the values from OUTCAR. We find that the ionic
dipole moment reported in OUTCAR differ from the naive calculation of
\\sum_i Z_i r_i where i is the index of the atom, Z_i is the ZVAL from the
pseudopotential file, and r is the distance in Angstroms along the lattice vectors.
Note, this difference is not simply due to VASP using Cartesian directions and
calc_ionic using lattice direction but rather how the ionic polarization is
computed. Compare calc_ionic to VASP SUBROUTINE POINT_CHARGE_DIPOL in dipol.F in
the VASP source to see the differences. We are able to recover a smooth same
branch polarization more frequently using the naive calculation in calc_ionic
than using the ionic dipole moment reported in the OUTCAR.

Some definitions of terms used in the comments below:

A polar structure belongs to a polar space group. A polar space group has a
one of the 10 polar point group:

> (1, 2, m, mm2, 4, 4mm, 3, 3m, 6, 6m)

Being nonpolar is not equivalent to being centrosymmetric (having inversion
symmetry). For example, any space group with point group 222 is nonpolar but
not centrosymmetric.

By symmetry the polarization of a nonpolar material modulo the quantum
of polarization can only be zero or 1/2. We use a nonpolar structure to help
determine the spontaneous polarization because it serves as a reference point.


### _class_ EnergyTrend(energies)
Bases: `object`

Class for fitting trends to energies.


* **Parameters**

    **energies** – Energies



#### endpoints_minima(slope_cutoff=0.005)
Test if spline endpoints are at minima for a given slope cutoff.


#### max_spline_jump()
Get maximum difference between spline and energy trend.


#### smoothness()
Get rms average difference between spline and energy trend.


#### spline()
Fit spline to energy trend data.


### _class_ Polarization(p_elecs, p_ions, structures, p_elecs_in_cartesian=True, p_ions_in_cartesian=False)
Bases: `object`

Class for recovering the same branch polarization for a set of
polarization calculations along the nonpolar - polar distortion
path of a ferroelectric.

p_elecs, p_ions, and structures lists should be given in order
of nonpolar to polar! For example, the structures returned from:

> nonpolar.interpolate(polar,interpolate_lattices=True)

if nonpolar is the nonpolar Structure and polar is the polar structure.

It is assumed that the electronic and ionic dipole moment values
are given in electron Angstroms along the three lattice directions
(a,b,c).

p_elecs: np.array of electronic contribution to the polarization with shape [N, 3]
p_ions: np.array of ionic contribution to the polarization with shape [N, 3]
p_elecs_in_cartesian: whether p_elecs is along Cartesian directions (rather than lattice directions).

> Default is True because that is the convention for VASP.

p_ions_in_cartesian: whether p_ions is along Cartesian directions (rather than lattice directions).

    Default is False because calc_ionic (which we recommend using for calculating the ionic
    contribution to the polarization) uses lattice directions.


#### _classmethod_ from_outcars_and_structures(outcars, structures, calc_ionic_from_zval=False)
Create Polarization object from list of Outcars and Structures in order
of nonpolar to polar.

Note, we recommend calculating the ionic dipole moment using calc_ionic
than using the values in Outcar (see module comments). To do this set
calc_ionic_from_zval = True


#### get_lattice_quanta(convert_to_muC_per_cm2=True, all_in_polar=True)
Returns the dipole / polarization quanta along a, b, and c for
all structures.


#### get_pelecs_and_pions(convert_to_muC_per_cm2=False)
Get the electronic and ionic dipole moments / polarizations.

convert_to_muC_per_cm2: Convert from electron \* Angstroms to microCoulomb

    per centimeter\*\*2


#### get_polarization_change(convert_to_muC_per_cm2=True, all_in_polar=True)
Get difference between nonpolar and polar same branch polarization.


#### get_polarization_change_norm(convert_to_muC_per_cm2=True, all_in_polar=True)
Get magnitude of difference between nonpolar and polar same branch
polarization.


#### get_same_branch_polarization_data(convert_to_muC_per_cm2=True, all_in_polar=True)
Get same branch dipole moment (convert_to_muC_per_cm2=False)
or polarization for given polarization data (convert_to_muC_per_cm2=True).

Polarization is a lattice vector, meaning it is only defined modulo the
quantum of polarization:

> P = P_0 + \\sum_i \\frac{n_i e R_i}{\\Omega}

where n_i is an integer, e is the charge of the electron in microCoulombs,
R_i is a lattice vector, and \\Omega is the unit cell volume in cm\*\*3
(giving polarization units of microCoulomb per centimeter\*\*2).

The quantum of the dipole moment in electron Angstroms (as given by VASP) is:

> \\sum_i n_i e R_i

where e, the electron charge, is 1 and R_i is a lattice vector, and n_i is an integer.

Given N polarization calculations in order from nonpolar to polar, this algorithm
minimizes the distance between adjacent polarization images. To do this, it
constructs a polarization lattice for each polarization calculation using the
pymatgen.core.structure class and calls the get_nearest_site method to find the
image of a given polarization lattice vector that is closest to the previous polarization
lattice vector image.

Note, using convert_to_muC_per_cm2=True and all_in_polar=True calculates the “proper
polarization” (meaning the change in polarization does not depend on the choice of
polarization branch) while convert_to_muC_per_cm2=True and all_in_polar=False calculates
the “improper polarization” (meaning the change in polarization does depend on the choice
of branch). As one might guess from the names. We recommend calculating the “proper
polarization”.

convert_to_muC_per_cm2: convert polarization from electron \* Angstroms to

    microCoulomb per centimeter\*\*2

all_in_polar: convert polarization to be in polar (final structure) polarization lattice


#### max_spline_jumps(convert_to_muC_per_cm2=True, all_in_polar=True)
Get maximum difference between spline and same branch polarization data.


#### same_branch_splines(convert_to_muC_per_cm2=True, all_in_polar=True)
Fit splines to same branch polarization. This is used to assess any jumps
in the same branch polarizaiton.


#### smoothness(convert_to_muC_per_cm2=True, all_in_polar=True)
Get rms average difference between spline and same branch polarization data.


### _class_ PolarizationLattice(lattice: ArrayLike | [Lattice](pymatgen.core.md#pymatgen.core.lattice.Lattice), species: Sequence[CompositionLike], coords: Sequence[ArrayLike], charge: float | None = None, validate_proximity: bool = False, to_unit_cell: bool = False, coords_are_cartesian: bool = False, site_properties: dict | None = None, labels: Sequence[str | None] | None = None, properties: dict | None = None)
Bases: [`Structure`](pymatgen.core.md#pymatgen.core.structure.Structure)

Why is a Lattice inheriting a structure? This is ridiculous.

Create a periodic structure.


* **Parameters**


    * **lattice** – The lattice, either as a pymatgen.core.Lattice or
    simply as any 2D array. Each row should correspond to a lattice
    vector. E.g., [[10,0,0], [20,10,0], [0,0,30]] specifies a
    lattice with lattice vectors [10,0,0], [20,10,0] and [0,0,30].


    * **species** – List of species on each site. Can take in flexible input,
    including:


        1. A sequence of element / species specified either as string
    symbols, e.g. [“Li”, “Fe2+”, “P”, …] or atomic numbers,
    e.g., (3, 56, …) or actual Element or Species objects.


        2. List of dict of elements/species and occupancies, e.g.,
    [{“Fe” : 0.5, “Mn”:0.5}, …]. This allows the setup of
    disordered structures.



    * **coords** (*Nx3 array*) – list of fractional/cartesian coordinates of
    each species.


    * **charge** (*int*) – overall charge of the structure. Defaults to behavior
    in SiteCollection where total charge is the sum of the oxidation
    states.


    * **validate_proximity** (*bool*) – Whether to check if there are sites
    that are less than 0.01 Ang apart. Defaults to False.


    * **to_unit_cell** (*bool*) – Whether to map all sites into the unit cell,
    i.e., fractional coords between 0 and 1. Defaults to False.


    * **coords_are_cartesian** (*bool*) – Set to True if you are providing
    coordinates in Cartesian coordinates. Defaults to False.


    * **site_properties** (*dict*) – Properties associated with the sites as a
    dict of sequences, e.g., {“magmom”:[5,5,5,5]}. The sequences
    have to be the same length as the atomic species and
    fractional_coords. Defaults to None for no properties.


    * **labels** (*list**[**str**]*) – Labels associated with the sites as a
    list of strings, e.g. [‘Li1’, ‘Li2’]. Must have the same
    length as the species and fractional coords. Defaults to
    None for no labels.


    * **properties** (*dict*) – Properties associated with the whole structure.
    Will be serialized when writing the structure to JSON or YAML but is
    lost when converting to other formats.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### _properties(_: dic_ )

#### get_nearest_site(coords, site, r=None)
Given coords and a site, find closet site to coords.


* **Parameters**


    * **coords** (*3x1 array*) – Cartesian coords of center of sphere


    * **site** – site to find closest to coords


    * **r** – radius of sphere. Defaults to diagonal of unit cell



* **Returns**

    Closest site and distance.



### calc_ionic(site, structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), zval)
Calculate the ionic dipole moment using ZVAL from pseudopotential.

site: PeriodicSite
structure: Structure
zval: Charge value for ion (ZVAL for VASP pseudopotential)

Returns polarization in electron Angstroms.


### get_total_ionic_dipole(structure, zval_dict)
Get the total ionic dipole moment for a structure.

structure: pymatgen Structure
zval_dict: specie, zval dictionary pairs
center (np.array with shape [3,1]) : dipole center used by VASP
tiny (float) : tolerance for determining boundary of calculation.


### zval_dict_from_potcar(potcar)
Creates zval_dictionary for calculating the ionic polarization from
Potcar object.

potcar: Potcar object