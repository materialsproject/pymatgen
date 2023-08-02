---
layout: default
title: pymatgen.core.periodic_table.md
nav_exclude: true
---

# pymatgen.core.periodic_table module

Module contains classes presenting Element, Species (Element + oxidation state) and PeriodicTable.

It should be noted that Element and Species are meant to be immutable objects.


### _class_ pymatgen.core.periodic_table.DummySpecie(symbol: str = 'X', oxidation_state: float | None = 0, properties: dict | None = None, spin: float | None = None)
Bases: `DummySpecies`

This maps the historical grammatically inaccurate DummySpecie to DummySpecies
to maintain backwards compatibility.


* **Parameters**


    * **symbol** (*str*) – An assigned symbol for the dummy specie. Strict
    rules are applied to the choice of the symbol. The dummy
    symbol cannot have any part of first two letters that will
    constitute an Element symbol. Otherwise, a composition may
    be parsed wrongly. E.g., “X” is fine, but “Vac” is not
    because Vac contains V, a valid Element.


    * **oxidation_state** (*float*) – Oxidation state for dummy specie. Defaults to 0.


    * **properties** – Properties associated with the Species, e.g. {“spin”: 5}. Defaults to None. This is now
    deprecated and retained purely for backward compatibility.


    * **spin** – Spin associated with Species. Defaults to None.



### _class_ pymatgen.core.periodic_table.DummySpecies(symbol: str = 'X', oxidation_state: float | None = 0, properties: dict | None = None, spin: float | None = None)
Bases: `Species`

A special specie for representing non-traditional elements or species. For
example, representation of vacancies (charged or otherwise), or special
sites, etc.


#### oxi_state()
Oxidation state associated with Species.


#### Z()
DummySpecies is always assigned an atomic number equal to the hash
number of the symbol. Obviously, it makes no sense whatsoever to use
the atomic number of a Dummy specie for anything scientific. The purpose
of this is to ensure that for most use cases, a DummySpecies behaves no
differently from an Element or Species.


#### X()
DummySpecies is always assigned a Pauling electronegativity of 0.


* **Parameters**


    * **symbol** (*str*) – An assigned symbol for the dummy specie. Strict
    rules are applied to the choice of the symbol. The dummy
    symbol cannot have any part of first two letters that will
    constitute an Element symbol. Otherwise, a composition may
    be parsed wrongly. E.g., “X” is fine, but “Vac” is not
    because Vac contains V, a valid Element.


    * **oxidation_state** (*float*) – Oxidation state for dummy specie. Defaults to 0.


    * **properties** – Properties associated with the Species, e.g. {“spin”: 5}. Defaults to None. This is now
    deprecated and retained purely for backward compatibility.


    * **spin** – Spin associated with Species. Defaults to None.



#### _property_ X(_: floa_ )
DummySpecies is always assigned a Pauling electronegativity of 0. The effect of
this is that DummySpecies are always sorted in front of actual Species.


#### _property_ Z(_: in_ )
DummySpecies is always assigned an atomic number equal to the hash of
the symbol. The expectation is that someone would be an actual dummy
to use atomic numbers for a Dummy specie.


#### as_dict()

* **Returns**

    MSONable dict representation.



#### _classmethod_ from_dict(d)

* **Parameters**

    **d** – Dict representation



* **Returns**

    DummySpecies



#### _static_ from_str(species_string: str)
Returns a Dummy from a string representation.


* **Parameters**

    **species_string** (*str*) – A string representation of a dummy
    species, e.g., “X2+”, “X3+”.



* **Returns**

    A DummySpecies object.



* **Raises**

    **ValueError if species_string cannot be interpreted.** –



#### _property_ oxi_state(_: float | Non_ )
Oxidation state associated with DummySpecies.


#### _property_ symbol(_: st_ )
Symbol for DummySpecies.


* **Type**

    return



### _class_ pymatgen.core.periodic_table.Element(value)
Bases: `ElementBase`

Enum representing an element in the periodic table.

Basic immutable element object with all relevant properties.

Only one instance of Element for each symbol is stored after creation,
ensuring that a particular element behaves like a singleton. For all
attributes, missing data (i.e., data for which is not available) is
represented by a None unless otherwise stated.


* **Parameters**

    **symbol** (*str*) – Element symbol, e.g., “H”, “Fe”



#### Z()
Atomic number


#### symbol()
Element symbol


#### long_name()
Long name for element. E.g., “Hydrogen”.


#### atomic_radius_calculated()
Calculated atomic radius for the element. This is the empirical value.
Data is obtained from
[http://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page](http://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)).


#### van_der_waals_radius()
Van der Waals radius for the element. This is the empirical
value determined from critical reviews of X-ray diffraction, gas kinetic
collision cross-section, and other experimental data by Bondi and later
workers. The uncertainty in these values is on the order of 0.1 Å.

Data are obtained from

“Atomic Radii of the Elements” in CRC Handbook of Chemistry and Physics,

    91st Ed.; Haynes, W.M., Ed.; CRC Press: Boca Raton, FL, 2010.


#### mendeleev_no()
Mendeleev number from definition given by Pettifor, D. G. (1984).
A chemical scale for crystal-structure maps. Solid State Communications,
51 (1), 31-34


#### electrical_resistivity()
Electrical resistivity


#### velocity_of_sound()
Velocity of sound


#### reflectivity()
Reflectivity


#### refractive_index()
Refractice index


#### poissons_ratio()
Poisson’s ratio


#### molar_volume()
Molar volume


#### electronic_structure()
Electronic structure.
E.g., The electronic structure for Fe is represented as
[Ar].3d6.4s2


#### atomic_orbitals()
Atomic Orbitals. Energy of the atomic orbitals as a dict.
E.g., The orbitals energies in eV are represented as
{‘1s’: -1.0, ‘2s’: -0.1}
Data is obtained from
[https://www.nist.gov/pml/data/atomic-reference-data-electronic-structure-calculations](https://www.nist.gov/pml/data/atomic-reference-data-electronic-structure-calculations)
The LDA values for neutral atoms are used


#### thermal_conductivity()
Thermal conductivity


#### boiling_point()
Boiling point


#### melting_point()
Melting point


#### critical_temperature()
Critical temperature


#### superconduction_temperature()
Superconduction temperature


#### liquid_range()
Liquid range


#### bulk_modulus()
Bulk modulus


#### youngs_modulus()
Young’s modulus


#### brinell_hardness()
Brinell hardness


#### rigidity_modulus()
Rigidity modulus


#### mineral_hardness()
Mineral hardness


#### vickers_hardness()
Vicker’s hardness


#### density_of_solid()
Density of solid phase


#### coefficient_of_linear_thermal_expansion()
Coefficient of linear thermal expansion


#### ground_level()
Ground level for element


#### ionization_energies()
List of ionization energies. First value is the first ionization energy, second is the second ionization
energy, etc. Note that this is zero-based indexing! So Element.ionization_energies[0] refer to the 1st
ionization energy. Values are from the NIST Atomic Spectra Database. Missing values are None.


#### Ac(_ = 'Ac_ )

#### Ag(_ = 'Ag_ )

#### Al(_ = 'Al_ )

#### Am(_ = 'Am_ )

#### Ar(_ = 'Ar_ )

#### As(_ = 'As_ )

#### At(_ = 'At_ )

#### Au(_ = 'Au_ )

#### B(_ = 'B_ )

#### Ba(_ = 'Ba_ )

#### Be(_ = 'Be_ )

#### Bh(_ = 'Bh_ )

#### Bi(_ = 'Bi_ )

#### Bk(_ = 'Bk_ )

#### Br(_ = 'Br_ )

#### C(_ = 'C_ )

#### Ca(_ = 'Ca_ )

#### Cd(_ = 'Cd_ )

#### Ce(_ = 'Ce_ )

#### Cf(_ = 'Cf_ )

#### Cl(_ = 'Cl_ )

#### Cm(_ = 'Cm_ )

#### Cn(_ = 'Cn_ )

#### Co(_ = 'Co_ )

#### Cr(_ = 'Cr_ )

#### Cs(_ = 'Cs_ )

#### Cu(_ = 'Cu_ )

#### Db(_ = 'Db_ )

#### Ds(_ = 'Ds_ )

#### Dy(_ = 'Dy_ )

#### Er(_ = 'Er_ )

#### Es(_ = 'Es_ )

#### Eu(_ = 'Eu_ )

#### F(_ = 'F_ )

#### Fe(_ = 'Fe_ )

#### Fl(_ = 'Fl_ )

#### Fm(_ = 'Fm_ )

#### Fr(_ = 'Fr_ )

#### Ga(_ = 'Ga_ )

#### Gd(_ = 'Gd_ )

#### Ge(_ = 'Ge_ )

#### H(_ = 'H_ )

#### He(_ = 'He_ )

#### Hf(_ = 'Hf_ )

#### Hg(_ = 'Hg_ )

#### Ho(_ = 'Ho_ )

#### Hs(_ = 'Hs_ )

#### I(_ = 'I_ )

#### In(_ = 'In_ )

#### Ir(_ = 'Ir_ )

#### K(_ = 'K_ )

#### Kr(_ = 'Kr_ )

#### La(_ = 'La_ )

#### Li(_ = 'Li_ )

#### Lr(_ = 'Lr_ )

#### Lu(_ = 'Lu_ )

#### Lv(_ = 'Lv_ )

#### Mc(_ = 'Mc_ )

#### Md(_ = 'Md_ )

#### Mg(_ = 'Mg_ )

#### Mn(_ = 'Mn_ )

#### Mo(_ = 'Mo_ )

#### Mt(_ = 'Mt_ )

#### N(_ = 'N_ )

#### Na(_ = 'Na_ )

#### Nb(_ = 'Nb_ )

#### Nd(_ = 'Nd_ )

#### Ne(_ = 'Ne_ )

#### Nh(_ = 'Nh_ )

#### Ni(_ = 'Ni_ )

#### No(_ = 'No_ )

#### Np(_ = 'Np_ )

#### O(_ = 'O_ )

#### Og(_ = 'Og_ )

#### Os(_ = 'Os_ )

#### P(_ = 'P_ )

#### Pa(_ = 'Pa_ )

#### Pb(_ = 'Pb_ )

#### Pd(_ = 'Pd_ )

#### Pm(_ = 'Pm_ )

#### Po(_ = 'Po_ )

#### Pr(_ = 'Pr_ )

#### Pt(_ = 'Pt_ )

#### Pu(_ = 'Pu_ )

#### Ra(_ = 'Ra_ )

#### Rb(_ = 'Rb_ )

#### Re(_ = 'Re_ )

#### Rf(_ = 'Rf_ )

#### Rg(_ = 'Rg_ )

#### Rh(_ = 'Rh_ )

#### Rn(_ = 'Rn_ )

#### Ru(_ = 'Ru_ )

#### S(_ = 'S_ )

#### Sb(_ = 'Sb_ )

#### Sc(_ = 'Sc_ )

#### Se(_ = 'Se_ )

#### Sg(_ = 'Sg_ )

#### Si(_ = 'Si_ )

#### Sm(_ = 'Sm_ )

#### Sn(_ = 'Sn_ )

#### Sr(_ = 'Sr_ )

#### Ta(_ = 'Ta_ )

#### Tb(_ = 'Tb_ )

#### Tc(_ = 'Tc_ )

#### Te(_ = 'Te_ )

#### Th(_ = 'Th_ )

#### Ti(_ = 'Ti_ )

#### Tl(_ = 'Tl_ )

#### Tm(_ = 'Tm_ )

#### Ts(_ = 'Ts_ )

#### U(_ = 'U_ )

#### V(_ = 'V_ )

#### W(_ = 'W_ )

#### Xe(_ = 'Xe_ )

#### Y(_ = 'Y_ )

#### Yb(_ = 'Yb_ )

#### Zn(_ = 'Zn_ )

#### Zr(_ = 'Zr_ )

### _class_ pymatgen.core.periodic_table.ElementBase(value)
Bases: `Enum`

Element class defined without any enum values so it can be subclassed.

This class is needed to get nested (as|from)_dict to work properly. All emmet classes that had
Element classes required custom construction whereas this definition behaves more like dataclasses
so serialization is less troublesome. There were many times where objects in as_dict serialized
only when they were top level. See [https://github.com/materialsproject/pymatgen/issues/2999](https://github.com/materialsproject/pymatgen/issues/2999).

Basic immutable element object with all relevant properties.

Only one instance of Element for each symbol is stored after creation,
ensuring that a particular element behaves like a singleton. For all
attributes, missing data (i.e., data for which is not available) is
represented by a None unless otherwise stated.


* **Parameters**

    **symbol** (*str*) – Element symbol, e.g., “H”, “Fe”



#### Z()
Atomic number


#### symbol()
Element symbol


#### long_name()
Long name for element. E.g., “Hydrogen”.


#### atomic_radius_calculated()
Calculated atomic radius for the element. This is the empirical value.
Data is obtained from
[http://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page](http://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)).


#### van_der_waals_radius()
Van der Waals radius for the element. This is the empirical
value determined from critical reviews of X-ray diffraction, gas kinetic
collision cross-section, and other experimental data by Bondi and later
workers. The uncertainty in these values is on the order of 0.1 Å.

Data are obtained from

“Atomic Radii of the Elements” in CRC Handbook of Chemistry and Physics,

    91st Ed.; Haynes, W.M., Ed.; CRC Press: Boca Raton, FL, 2010.


#### mendeleev_no()
Mendeleev number from definition given by Pettifor, D. G. (1984).
A chemical scale for crystal-structure maps. Solid State Communications,
51 (1), 31-34


#### electrical_resistivity()
Electrical resistivity


#### velocity_of_sound()
Velocity of sound


#### reflectivity()
Reflectivity


#### refractive_index()
Refractice index


#### poissons_ratio()
Poisson’s ratio


#### molar_volume()
Molar volume


#### electronic_structure()
Electronic structure.
E.g., The electronic structure for Fe is represented as
[Ar].3d6.4s2


#### atomic_orbitals()
Atomic Orbitals. Energy of the atomic orbitals as a dict.
E.g., The orbitals energies in eV are represented as
{‘1s’: -1.0, ‘2s’: -0.1}
Data is obtained from
[https://www.nist.gov/pml/data/atomic-reference-data-electronic-structure-calculations](https://www.nist.gov/pml/data/atomic-reference-data-electronic-structure-calculations)
The LDA values for neutral atoms are used


#### thermal_conductivity()
Thermal conductivity


#### boiling_point()
Boiling point


#### melting_point()
Melting point


#### critical_temperature()
Critical temperature


#### superconduction_temperature()
Superconduction temperature


#### liquid_range()
Liquid range


#### bulk_modulus()
Bulk modulus


#### youngs_modulus()
Young’s modulus


#### brinell_hardness()
Brinell hardness


#### rigidity_modulus()
Rigidity modulus


#### mineral_hardness()
Mineral hardness


#### vickers_hardness()
Vicker’s hardness


#### density_of_solid()
Density of solid phase


#### coefficient_of_linear_thermal_expansion()
Coefficient of linear thermal expansion


#### ground_level()
Ground level for element


#### ionization_energies()
List of ionization energies. First value is the first ionization energy, second is the second ionization
energy, etc. Note that this is zero-based indexing! So Element.ionization_energies[0] refer to the 1st
ionization energy. Values are from the NIST Atomic Spectra Database. Missing values are None.


#### _property_ X(_: floa_ )
Pauling electronegativity of element. Note that if an element does not
have an Pauling electronegativity, a NaN float is returned.


* **Type**

    return



#### as_dict()
Makes Element obey the general json interface used in pymatgen for
easier serialization.


#### _property_ atomic_mass(_: [FloatWithUnit](pymatgen.core.units.md#pymatgen.core.units.FloatWithUnit_ )
Returns:
float: The atomic mass of the element in amu.


#### _property_ atomic_radius(_: [FloatWithUnit](pymatgen.core.units.md#pymatgen.core.units.FloatWithUnit) | Non_ )
Returns:
float | None: The atomic radius of the element in Ångstroms. Can be None for

> some elements like noble gases.


#### _property_ average_anionic_radius(_: floa_ )
Average anionic radius for element (with units). The average is
taken over all negative oxidation states of the element for which
data is present.


#### _property_ average_cationic_radius(_: [FloatWithUnit](pymatgen.core.units.md#pymatgen.core.units.FloatWithUnit_ )
Average cationic radius for element (with units). The average is
taken over all positive oxidation states of the element for which
data is present.


#### _property_ average_ionic_radius(_: [FloatWithUnit](pymatgen.core.units.md#pymatgen.core.units.FloatWithUnit_ )
Average ionic radius for element (with units). The average is taken
over all oxidation states of the element for which data is present.


#### _property_ block(_: st_ )
Return the block character “s,p,d,f”.


#### _property_ common_oxidation_states(_: tuple[int, ..._ )
Tuple of common oxidation states.


#### _property_ data(_: dict[str, Any_ )
Returns dict of data for element.


#### _property_ electron_affinity(_: floa_ )
The amount of energy released when an electron is attached to a neutral atom.


#### _property_ electronic_structure(_: st_ )
Electronic structure as string, with only valence electrons.
E.g., The electronic structure for Fe is represented as ‘[Ar].3d6.4s2’.


#### _static_ from_Z(Z: int)
Get an element from an atomic number.


* **Parameters**

    **Z** (*int*) – Atomic number



* **Returns**

    Element with atomic number Z.



#### _static_ from_dict(d)
Makes Element obey the general json interface used in pymatgen for
easier serialization.


#### _static_ from_name(name: str)
Get an element from its long name.


* **Parameters**

    **name** – Long name of the element, e.g. ‘Hydrogen’ or
    ‘Iron’. Not case-sensitive.



* **Returns**

    Element with the name ‘name’



#### _static_ from_row_and_group(row: int, group: int)
Returns an element from a row and group number.
Important Note: For lanthanoids and actinoids, the row number must
be 8 and 9, respectively, and the group number must be
between 3 (La, Ac) and 17 (Lu, Lr). This is different than the
value for Element(symbol).row and Element(symbol).group for these
elements.


* **Parameters**


    * **row** (*int*) – (pseudo) row number. This is the
    standard row number except for the lanthanoids
    and actinoids for which it is 8 or 9, respectively.


    * **group** (*int*) – (pseudo) group number. This is the
    standard group number except for the lanthanoids
    and actinoids for which it is 3 (La, Ac) to 17 (Lu, Lr).


**NOTE**: The 18 group number system is used, i.e., Noble gases are group 18.


#### _property_ full_electronic_structure(_: list[tuple[int, str, int]_ )
Full electronic structure as tuple.
E.g., The electronic structure for Fe is represented as:
[(1, “s”, 2), (2, “s”, 2), (2, “p”, 6), (3, “s”, 2), (3, “p”, 6),
(3, “d”, 6), (4, “s”, 2)].


#### _property_ ground_state_term_symbol()
Ground state term symbol
Selected based on Hund’s Rule.


#### _property_ group(_: in_ )
Returns the periodic table group of the element.
Note: For lanthanoids and actinoids, the group is always 3.


#### _property_ icsd_oxidation_states(_: tuple[int, ..._ )
Tuple of all oxidation states with at least 10 instances in
ICSD database AND at least 1% of entries for that element.


#### _property_ ionic_radii(_: dict[int, float_ )
All ionic radii of the element as a dict of
{oxidation state: ionic radii}. Radii are given in angstrom.


#### _property_ ionization_energy(_: floa_ )
First ionization energy of element.


#### _property_ is_actinoid(_: boo_ )
True if element is a actinoid.


#### _property_ is_alkali(_: boo_ )
True if element is an alkali metal.


#### _property_ is_alkaline(_: boo_ )
True if element is an alkaline earth metal (group II).


#### _property_ is_chalcogen(_: boo_ )
True if element is a chalcogen.


#### _property_ is_halogen(_: boo_ )
True if element is a halogen.


#### _property_ is_lanthanoid(_: boo_ )
True if element is a lanthanoid.


#### _property_ is_metal(_: boo_ )
True if is a metal.


#### _property_ is_metalloid(_: boo_ )
True if element is a metalloid.


#### _property_ is_noble_gas(_: boo_ )
True if element is noble gas.


#### _property_ is_post_transition_metal(_: boo_ )
True if element is a post-transition or poor metal.


#### _property_ is_quadrupolar(_: boo_ )
Checks if this element can be quadrupolar.


#### _property_ is_rare_earth_metal(_: boo_ )
True if element is a rare earth metal.


#### _property_ is_transition_metal(_: boo_ )
True if element is a transition metal.


#### _static_ is_valid_symbol(symbol: str)
Returns true if symbol is a valid element symbol.


* **Parameters**

    **symbol** (*str*) – Element symbol



* **Returns**

    True if symbol is a valid element (e.g., “H”). False otherwise
    (e.g., “Zebra”).



#### _property_ iupac_ordering()
Ordering according to Table VI of “Nomenclature of Inorganic Chemistry
(IUPAC Recommendations 2005)”. This ordering effectively follows the
groups and rows of the periodic table, except the Lanthanides, Actinides
and hydrogen.


#### _property_ max_oxidation_state(_: floa_ )
Maximum oxidation state for element.


#### _property_ metallic_radius(_: floa_ )
Metallic radius of the element. Radius is given in ang.


#### _property_ min_oxidation_state(_: floa_ )
Minimum oxidation state for element.


#### _property_ nmr_quadrupole_moment(_: dict[str, [pymatgen.core.units.FloatWithUnit](pymatgen.core.units.md#pymatgen.core.units.FloatWithUnit)_ )
Get a dictionary the nuclear electric quadrupole moment in units of
e\*millibarns for various isotopes.


#### _property_ number(_: in_ )
Alternative attribute for atomic number Z.


#### _property_ oxidation_states(_: tuple[int, ..._ )
Tuple of all known oxidation states.


#### _static_ print_periodic_table(filter_function: Callable | None = None)
A pretty ASCII printer for the periodic table, based on some
filter_function.


* **Parameters**

    **filter_function** – A filtering function taking an Element as input
    and returning a boolean. For example, setting
    filter_function = lambda el: el.X > 2 will print a periodic
    table containing only elements with Pauling electronegativity > 2.



#### _property_ row(_: in_ )
Returns the periodic table row of the element.
Note: For lanthanoids and actinoids, the row is always 6 or 7,
respectively.


#### _property_ term_symbols(_: list[list[str]_ )
All possible  Russell-Saunders term symbol of the Element.
eg. L = 1, n_e = 2 (s2) returns [[‘1D2’], [‘3P0’, ‘3P1’, ‘3P2’], [‘1S0’]].


#### _property_ valence()
From full electron config obtain valence subshell angular moment (L) and number of valence e- (v_e).


### _class_ pymatgen.core.periodic_table.Specie(symbol: SpeciesLike, oxidation_state: float | None = None, properties: dict | None = None, spin: float | None = None)
Bases: `Species`

This maps the historical grammatically inaccurate Specie to Species
to maintain backwards compatibility.


* **Parameters**


    * **symbol** (*str*) – Element symbol optionally incl. oxidation state. E.g. Fe, Fe2+, O2-.


    * **oxidation_state** (*float*) – Explicit oxidation state of element, e.g. -2, -1, 0, 1, 2, …
    If oxidation state is present in symbol, this argument is ignored.


    * **properties** – Properties associated with the Species, e.g., {“spin”: 5}. Defaults to None. This is now
    deprecated and retained purely for backward compatibility.


    * **spin** – Spin associated with Species. Defaults to None.



* **Raises**

    **ValueError** – If oxidation state passed both in symbol string and via
        oxidation_state kwarg.



### _class_ pymatgen.core.periodic_table.Species(symbol: SpeciesLike, oxidation_state: float | None = None, properties: dict | None = None, spin: float | None = None)
Bases: `MSONable`, [`Stringify`](pymatgen.util.string.md#pymatgen.util.string.Stringify)

An extension of Element with an oxidation state and other optional
properties. Properties associated with Species should be “idealized”
values, not calculated values. For example, high-spin Fe2+ may be
assigned an idealized spin of +5, but an actual Fe2+ site may be
calculated to have a magmom of +4.5. Calculated properties should be
assigned to Site objects, and not Species.


* **Parameters**


    * **symbol** (*str*) – Element symbol optionally incl. oxidation state. E.g. Fe, Fe2+, O2-.


    * **oxidation_state** (*float*) – Explicit oxidation state of element, e.g. -2, -1, 0, 1, 2, …
    If oxidation state is present in symbol, this argument is ignored.


    * **properties** – Properties associated with the Species, e.g., {“spin”: 5}. Defaults to None. This is now
    deprecated and retained purely for backward compatibility.


    * **spin** – Spin associated with Species. Defaults to None.



* **Raises**

    **ValueError** – If oxidation state passed both in symbol string and via
        oxidation_state kwarg.



#### STRING_MODE(_ = 'SUPERSCRIPT_ )

#### as_dict()

* **Returns**

    Json-able dictionary representation.



#### _property_ element(_: Elemen_ )
Underlying element object.


#### _classmethod_ from_dict(d)

* **Parameters**

    **d** – Dict representation.



* **Returns**

    Species.



#### _static_ from_str(species_string: str)
Returns a Species from a string representation.


* **Parameters**

    **species_string** (*str*) – A typical string representation of a
    species, e.g., “Mn2+”, “Fe3+”, “O2-“.



* **Returns**

    A Species object.



* **Raises**

    **ValueError if species_string cannot be interpreted.** –



#### _classmethod_ from_string(\*args, \*\*kwds)
from_string is deprecated!
Use from_str instead

Use from_str instead.


#### get_crystal_field_spin(coordination: Literal['oct', 'tet'] = 'oct', spin_config: Literal['low', 'high'] = 'high')
Calculate the crystal field spin based on coordination and spin
configuration. Only works for transition metal species.


* **Parameters**


    * **coordination** (*"oct"** | **"tet"*) – Tetrahedron or octahedron crystal site coordination


    * **spin_config** (*"low"** | **"high"*) – Whether the species is in a high or low spin state



* **Returns**

    Crystal field spin in Bohr magneton.



* **Raises**


    * **AttributeError if species is not a valid transition metal**** or ****has** – an invalid oxidation state.


    * **ValueError if invalid coordination**** or ****spin_config.** –



#### get_nmr_quadrupole_moment(isotope: str | None = None)
Gets the nuclear electric quadrupole moment in units of e \* millibarns.


* **Parameters**

    **isotope** (*str*) – the isotope to get the quadrupole moment for
    default is None, which gets the lowest mass isotope



#### get_shannon_radius(cn: str, spin: Literal['', 'Low Spin', 'High Spin'] = '', radius_type: Literal['ionic', 'crystal'] = 'ionic')
Get the local environment specific ionic radius for species.


* **Parameters**


    * **cn** (*str*) – Coordination using roman letters. Supported values are
    I-IX, as well as IIIPY, IVPY and IVSQ.


    * **spin** (*str*) – Some species have different radii for different
    spins. You can get specific values using “High Spin” or
    “Low Spin”. Leave it as “” if not available. If only one spin
    data is available, it is returned and this spin parameter is
    ignored.


    * **radius_type** (*str*) – Either “crystal” or “ionic” (default).



* **Returns**

    Shannon radius for specie in the specified environment.



#### _property_ ionic_radius(_: float | Non_ )
Ionic radius of specie. Returns None if data is not present.


#### _property_ oxi_state(_: float | Non_ )
Oxidation state of Species.


#### _property_ properties(_: dic_ )
Retained for backwards incompatibility.


#### _property_ spin(_: float | Non_ )
Spin of Species.


#### to_pretty_string()

* **Returns**

    String without properties.



### pymatgen.core.periodic_table.get_el_sp(obj: int | SpeciesLike)
Utility method to get an Element, Species or DummySpecies from any input.

If obj is in itself an element or a specie, it is returned automatically.
If obj is an int or a string representing an integer, the Element with the
atomic number obj is returned.
If obj is a string, Species parsing will be attempted (e.g. Mn2+). Failing that
Element parsing will be attempted (e.g. Mn). Failing that DummyElement parsing
will be attempted.


* **Parameters**

    **obj** (*Element/Species/str/int*) – An arbitrary object. Supported objects
    are actual Element/Species objects, integers (representing atomic
    numbers) or strings (element symbols or species strings).



* **Returns**

    Species or Element, with a bias for the maximum number of properties
    that can be determined.



* **Raises**

    **ValueError** – if obj cannot be converted into an Element or Species.