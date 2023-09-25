---
layout: default
title: pymatgen.ext.md
nav_exclude: true
---

1. TOC
{:toc}

# pymatgen.ext namespace


## pymatgen.ext.cod module

This module provides classes to interface with the Crystallography Open
Database. If you use data from the COD, please cite the following works (as
stipulated by the COD developers).

> Merkys, A., Vaitkus, A., Butkus, J., Okulič-Kazarinas, M., Kairys, V. &
> Gražulis, S. (2016) “COD::CIF::Parser: an error-correcting CIF parser for
> the Perl language”. Journal of Applied Crystallography 49.

> Gražulis, S., Merkys, A., Vaitkus, A. & Okulič-Kazarinas, M. (2015)
> “Computing stoichiometric molecular composition from crystal structures”.
> Journal of Applied Crystallography 48, 85-91.

> Gražulis, S., Daškevič, A., Merkys, A., Chateigner, D., Lutterotti, L.,
> Quirós, M., Serebryanaya, N. R., Moeck, P., Downs, R. T. & LeBail, A.
> (2012) “Crystallography Open Database (COD): an open-access collection of
> crystal structures and platform for world-wide collaboration”. Nucleic
> Acids Research 40, D420-D427.

> Grazulis, S., Chateigner, D., Downs, R. T., Yokochi, A. T., Quiros, M.,
> Lutterotti, L., Manakova, E., Butkus, J., Moeck, P. & Le Bail, A. (2009)
> “Crystallography Open Database - an open-access collection of crystal
> structures”. J. Appl. Cryst. 42, 726-729.

> Downs, R. T. & Hall-Wallace, M. (2003) “The American Mineralogist Crystal
> Structure Database”. American Mineralogist 88, 247-250.


### _class_ COD()
Bases: `object`

An interface to the Crystallography Open Database.


#### get_cod_ids(formula)
Queries the COD for all cod ids associated with a formula. Requires
mysql executable to be in the path.


* **Parameters**

    **formula** (*str*) – Formula.



* **Returns**

    List of cod ids.



#### get_structure_by_formula(formula: str, \*\*kwargs)
Queries the COD for structures by formula. Requires mysql executable to
be in the path.


* **Parameters**


    * **formula** (*str*) – Chemical formula.


    * **kwargs** – All kwargs supported by
    `pymatgen.core.structure.Structure.from_str()`.



* **Returns**

    Structure, “cod_id”: int, “sg”: “P n m a”}]



* **Return type**

    A list of dict of the format [{“structure”



#### get_structure_by_id(cod_id, \*\*kwargs)
Queries the COD for a structure by id.


* **Parameters**


    * **cod_id** (*int*) – COD id.


    * **kwargs** – All kwargs supported by
    `pymatgen.core.structure.Structure.from_str()`.



* **Returns**

    A Structure.



#### query(sql: str)
Perform a query.


* **Parameters**

    **sql** – SQL string



* **Returns**

    Response from SQL query.



#### url(_ = 'www.crystallography.net_ )
## pymatgen.ext.matproj module

This module provides classes to interface with the Materials Project REST
API v2 to enable the creation of data structures and pymatgen objects using
Materials Project data.

To make use of the Materials API, you need to be a registered user of the
Materials Project, and obtain an API key by going to your dashboard at
[https://materialsproject.org/dashboard](https://materialsproject.org/dashboard).


### _exception_ MPRestError()
Bases: `Exception`

Exception class for legacy MPRestAdaptor. Raised when query is malformed.


### _class_ MPRester(\*args, \*\*kwargs)
Bases: `object`

A class to conveniently interface with the new and legacy Materials Project REST interface.

The recommended way to use MPRester is as a context manager to ensure
that sessions are properly closed after usage:

> with MPRester(“API_KEY”) as mpr:

>     docs = mpr.call_some_method()

MPRester uses the “requests” package, which provides HTTP connection
pooling. All connections are made via https for security.

For more advanced uses of the Materials API, please consult the API
documentation at [https://materialsproject.org/api](https://materialsproject.org/api) and [https://docs.materialsproject.org](https://docs.materialsproject.org).

This class handles the transition between old and new MP API, making it easy to switch between them
by passing a new (length 32) or old (15 <= length <= 17) API key. See [https://docs.materialsproject.org](https://docs.materialsproject.org)
for which API to use.


* **Parameters**


    * **\*args** – Pass through to either legacy or new MPRester.


    * **\*\*kwargs** – Pass through to either legacy or new MPRester.



### _class_ TaskType(value)
Bases: `Enum`

task types available in legacy MP data.


#### GGAU_DEF(_ = 'GGA+U Deformation_ )

#### GGAU_LINE(_ = 'GGA+U NSCF Line_ )

#### GGAU_OPT(_ = 'GGA+U Structure Optimization_ )

#### GGAU_STATIC(_ = 'GGA+U Static_ )

#### GGAU_STATIC_DIEL(_ = 'GGA+U Static Dielectric_ )

#### GGAU_UNIFORM(_ = 'GGA+U NSCF Uniform_ )

#### GGA_DEF(_ = 'GGA Deformation_ )

#### GGA_LINE(_ = 'GGA NSCF Line_ )

#### GGA_OPT(_ = 'GGA Structure Optimization_ )

#### GGA_STATIC(_ = 'GGA Static_ )

#### GGA_STATIC_DIEL(_ = 'GGA Static Dielectric_ )

#### GGA_UNIFORM(_ = 'GGA NSCF Uniform_ )

#### LDA_STATIC_DIEL(_ = 'LDA Static Dielectric_ )

#### SCAN_OPT(_ = 'SCAN Structure Optimization_ )

### _class_ _MPResterLegacy(api_key: str | None = None, endpoint: str | None = None, notify_db_version: bool = True, include_user_agent: bool = True)
Bases: `object`

A class to conveniently interface with the Materials Project REST interface.
The recommended way to use MPRester is with the “with” context manager to ensure
sessions are properly closed after usage.

> with MPRester(“API_KEY”) as mpr:

>     mpr.some_method()

MPRester uses the “requests” package, which provides for HTTP connection
pooling. All connections are made via https for security.

For more advanced uses of the legacy Materials API, please consult the API
documentation at [https://github.com/materialsproject/mapidoc](https://github.com/materialsproject/mapidoc).

Note that this class is for the *legacy* API. Upcoming changes to the
Materials Project api are described at [https://materialsproject.org/api](https://materialsproject.org/api).


* **Parameters**


    * **api_key** (*str*) – A String API key for accessing the MaterialsProject
    REST interface. Please obtain your API key at
    [https://materialsproject.org/dashboard](https://materialsproject.org/dashboard). If this is None,
    the code will check if there is a “PMG_MAPI_KEY” setting.
    If so, it will use that environment variable. This makes
    easier for heavy users to simply add this environment variable to
    their setups and MPRester can then be called without any arguments.


    * **endpoint** (*str*) – Url of endpoint to access the MaterialsProject REST
    interface. Defaults to the standard Materials Project REST
    address at “[https://legacy.materialsproject.org/rest/v2](https://legacy.materialsproject.org/rest/v2)”, but
    can be changed to other urls implementing a similar interface.


    * **notify_db_version** (*bool*) – If True, the current MP database version will
    be retrieved and logged locally in the ~/.pmgrc.yaml. If the database
    version changes, you will be notified. The current database version is
    also printed on instantiation. These local logs are not sent to
    materialsproject.org and are not associated with your API key, so be
    aware that a notification may not be presented if you run MPRester
    from multiple computing environments.


    * **include_user_agent** (*bool*) – If True, will include a user agent with the
    HTTP request including information on pymatgen and system version
    making the API request. This helps MP support pymatgen users, and
    is similar to what most web browsers send with each page request.
    Set to False to disable the user agent.



#### _check_get_download_info_url_by_task_id(prefix, task_ids)

#### _static_ _check_nomad_exist(url)

#### _make_request(sub_url: str, payload: Any | None = None, method: Literal['GET', 'POST', 'PUT', 'DELETE'] = 'GET', mp_decode: bool = True)

#### _static_ _print_help_message(nomad_exist_task_ids, task_ids, file_patterns, task_types)

#### delete_snl(snl_ids)
Delete earlier submitted SNLs.

**NOTE**: As of now, this MP REST feature is open only to a select group of
users. Opening up submissions to all users is being planned for the future.


* **Parameters**

    **snl_ids** – List of SNL ids.



* **Raises**

    **MPRestError** –



#### find_structure(filename_or_structure)
Finds matching structures on the Materials Project site.


* **Parameters**

    **filename_or_structure** – filename or Structure object



* **Returns**

    A list of matching materials project ids for structure.



* **Raises**

    **MPRestError** –



#### get_all_substrates()
Gets the list of all possible substrates considered in the
Materials Project substrate database.


* **Returns**

    list of material_ids corresponding to possible substrates



#### get_bandstructure_by_material_id(material_id, line_mode=True)
Get a BandStructure corresponding to a material_id.

REST Endpoint: [https://materialsproject.org/rest/v2/materials](https://materialsproject.org/rest/v2/materials)/<mp-id>/vasp/bandstructure or
[https://materialsproject.org/rest/v2/materials](https://materialsproject.org/rest/v2/materials)/<mp-id>/vasp/bandstructure_uniform


* **Parameters**


    * **material_id** (*str*) – Materials Project material_id.


    * **line_mode** (*bool*) – If True, fetch a BandStructureSymmLine object
    (default). If False, return the uniform band structure.



* **Returns**

    A BandStructure object.



#### get_cohesive_energy(material_id, per_atom=False)
Gets the cohesive for a material (eV per formula unit). Cohesive energy

    is defined as the difference between the bulk energy and the sum of
    total DFT energy of isolated atoms for atom elements in the bulk.


* **Parameters**


    * **material_id** (*str*) – Materials Project material_id, e.g. ‘mp-123’.


    * **per_atom** (*bool*) – Whether or not to return cohesive energy per atom



* **Returns**

    Cohesive energy (eV).



#### get_data(chemsys_formula_id, data_type='vasp', prop='')
Flexible method to get any data using the Materials Project REST
interface. Generally used by other methods for more specific queries.

Format of REST return is *always* a list of dict (regardless of the
number of pieces of data returned. The general format is as follows:

[{“material_id”: material_id, “property_name” : value}, …]

This is generally a call to
[https://materialsproject.org/rest/v2/materials/vasp](https://materialsproject.org/rest/v2/materials/vasp)/<prop>.
See [https://github.com/materialsproject/mapidoc](https://github.com/materialsproject/mapidoc) for details.


* **Parameters**


    * **chemsys_formula_id** (*str*) – A chemical system (e.g., Li-Fe-O),
    or formula (e.g., Fe2O3) or materials_id (e.g., mp-1234).


    * **data_type** (*str*) – Type of data to return. Currently can either be
    “vasp” or “exp”.


    * **prop** (*str*) – Property to be obtained. Should be one of the
    MPRester.supported_task_properties. Leave as empty string for a
    general list of useful properties.



#### get_database_version()
The Materials Project database is periodically updated and has a
database version associated with it. When the database is updated,
consolidated data (information about “a material”) may and does
change, while calculation data about a specific calculation task
remains unchanged and available for querying via its task_id.

The database version is set as a date in the format YYYY-MM-DD,
where “-DD” may be optional. An additional numerical suffix
might be added if multiple releases happen on the same day.


* **Returns**

    database version



* **Return type**

    str



#### get_doc(materials_id)
Get the entire data document for one materials id. Use this judiciously.

REST Endpoint: [https://materialsproject.org/materials](https://materialsproject.org/materials)/<mp-id>/doc.


* **Parameters**

    **materials_id** (*str*) – E.g., mp-1143 for Al2O3



* **Returns**

    Dict of json document of all data that is displayed on a materials
    details page.



#### get_dos_by_material_id(material_id)
Get a Dos corresponding to a material_id.

REST Endpoint: [https://materialsproject.org/rest/v2/materials](https://materialsproject.org/rest/v2/materials)/<mp-id>/vasp/dos


* **Parameters**

    **material_id** (*str*) – Materials Project material_id (a string,
    e.g., mp-1234).



* **Returns**

    A Dos object.



#### get_download_info(material_ids, task_types=None, file_patterns=None)
Get a list of URLs to retrieve raw VASP output files from the NoMaD repository.


* **Parameters**


    * **material_ids** (*list*) – list of material identifiers (mp-id’s)


    * **task_types** (*list*) – list of task types to include in download (see TaskType Enum class)


    * **file_patterns** (*list*) – list of wildcard file names to include for each task



* **Returns**

    a tuple of 1) a dictionary mapping material_ids to task_ids and
    task_types, and 2) a list of URLs to download zip archives from
    NoMaD repository. Each zip archive will contain a manifest.json with
    metadata info, e.g. the task/external_ids that belong to a directory



#### get_entries(chemsys_formula_id_criteria: str | dict[str, Any], compatible_only: bool = True, inc_structure: bool | Literal['initial'] | None = None, property_data: list[str] | None = None, conventional_unit_cell: bool = False, sort_by_e_above_hull: bool = False)
Get a list of ComputedEntries or ComputedStructureEntries corresponding
to a chemical system, formula, or materials_id or full criteria.


* **Parameters**


    * **chemsys_formula_id_criteria** (*str/dict*) – A chemical system
    (e.g., Li-Fe-O), or formula (e.g., Fe2O3) or materials_id
    (e.g., mp-1234) or full Mongo-style dict criteria.


    * **compatible_only** (*bool*) – Whether to return only “compatible”
    entries. Compatible entries are entries that have been
    processed using the MaterialsProject2020Compatibility class,
    which performs adjustments to allow mixing of GGA and GGA+U
    calculations for more accurate phase diagrams and reaction
    energies.


    * **inc_structure** (*str*) – If None, entries returned are
    ComputedEntries. If inc_structure=”initial”,
    ComputedStructureEntries with initial structures are returned.
    Otherwise, ComputedStructureEntries with final structures
    are returned.


    * **property_data** (*list*) – Specify additional properties to include in
    entry.data. If None, no data. Should be a subset of
    supported_properties.


    * **conventional_unit_cell** (*bool*) – Whether to get the standard
    conventional unit cell


    * **sort_by_e_above_hull** (*bool*) – Whether to sort the list of entries by
    e_above_hull (will query e_above_hull as a property_data if True).



* **Returns**

    List of ComputedEntry or ComputedStructureEntry objects.



#### get_entries_in_chemsys(elements, compatible_only=True, inc_structure=None, property_data=None, conventional_unit_cell=False, additional_criteria=None)
Helper method to get a list of ComputedEntries in a chemical system.

For example, elements = [“Li”, “Fe”, “O”] will return a list of all entries in the
Li-Fe-O chemical system, i.e., all LixOy, FexOy, LixFey, LixFeyOz, Li, Fe and O
phases. Extremely useful for creating phase diagrams of entire chemical systems.


* **Parameters**


    * **elements** (*str** or **[**str**]*) – Chemical system string comprising element
    symbols separated by dashes, e.g., “Li-Fe-O” or List of element
    symbols, e.g., [“Li”, “Fe”, “O”].


    * **compatible_only** (*bool*) – Whether to return only “compatible”
    entries. Compatible entries are entries that have been
    processed using the MaterialsProject2020Compatibility class,
    which performs adjustments to allow mixing of GGA and GGA+U
    calculations for more accurate phase diagrams and reaction
    energies.


    * **inc_structure** (*str*) – If None, entries returned are
    ComputedEntries. If inc_structure=”initial”,
    ComputedStructureEntries with initial structures are returned.
    Otherwise, ComputedStructureEntries with final structures
    are returned.


    * **property_data** (*list*) – Specify additional properties to include in
    entry.data. If None, no data. Should be a subset of
    supported_properties.


    * **conventional_unit_cell** (*bool*) – Whether to get the standard
    conventional unit cell


    * **additional_criteria** (*dict*) – Any additional criteria to pass. For instance, if you are only interested in
    stable entries, you can pass {“e_above_hull”: {“$lte”: 0.001}}.



* **Returns**

    List of ComputedEntries.



#### get_entry_by_material_id(material_id: str, compatible_only: bool = True, inc_structure: bool | Literal['initial'] | None = None, property_data: list[str] | None = None, conventional_unit_cell: bool = False)
Get a ComputedEntry corresponding to a material_id.


* **Parameters**


    * **material_id** (*str*) – Materials Project material_id (a string,
    e.g., mp-1234).


    * **compatible_only** (*bool*) – Whether to return only “compatible”
    entries. Compatible entries are entries that have been
    processed using the MaterialsProject2020Compatibility class,
    which performs adjustments to allow mixing of GGA and GGA+U
    calculations for more accurate phase diagrams and reaction
    energies.


    * **inc_structure** (*str*) – If None, entries returned are
    ComputedEntries. If inc_structure=”initial”,
    ComputedStructureEntries with initial structures are returned.
    Otherwise, ComputedStructureEntries with final structures
    are returned.


    * **property_data** (*list*) – Specify additional properties to include in
    entry.data. If None, no data. Should be a subset of
    supported_properties.


    * **conventional_unit_cell** (*bool*) – Whether to get the standard
    conventional unit cell



* **Raises**

    **MPRestError if no data for given material_id is found.** –



* **Returns**

    ComputedEntry or ComputedStructureEntry object.



#### get_exp_entry(formula)
Returns an ExpEntry object, which is the experimental equivalent of a
ComputedEntry and can be used for analyses using experimental data.


* **Parameters**

    **formula** (*str*) – A formula to search for.



* **Returns**

    An ExpEntry object.



#### get_exp_thermo_data(formula)
Get a list of ThermoData objects associated with a formula using the
Materials Project REST interface.


* **Parameters**

    **formula** (*str*) – A formula to search for.



* **Returns**

    List of ThermoData objects.



#### get_gb_data(material_id=None, pretty_formula=None, chemsys=None, sigma=None, gb_plane=None, rotation_axis=None, include_work_of_separation=False)
Gets grain boundary data for a material.


* **Parameters**


    * **material_id** (*str*) – Materials Project material_id, e.g., ‘mp-129’.


    * **pretty_formula** (*str*) – The formula of metals. e.g., ‘Fe’


    * **chemsys** (*str*) – The chemical system. e.g., ‘Fe-O’


    * **sigma** (*int*) – The sigma value of a certain type of grain boundary


    * **gb_plane** (*list** of **integer*) – The Miller index of grain boundary plane. e.g., [1, 1, 1]


    * **rotation_axis** (*list** of **integer*) – The Miller index of rotation axis. e.g.,
    [1, 0, 0], [1, 1, 0], and [1, 1, 1] Sigma value is determined by the combination of
    rotation axis and rotation angle. The five degrees of freedom (DOF) of one grain boundary
    include: rotation axis (2 DOFs), rotation angle (1 DOF), and grain boundary plane (2 DOFs).


    * **include_work_of_separation** (*bool*) – whether to include the work of separation
    (in unit of (J/m^2)). If you want to query the work of separation, please
    specify the material_id.



* **Returns**

    A list of grain boundaries that satisfy the query conditions (sigma, gb_plane).
    Energies are given in SI units (J/m^2).



#### get_interface_reactions(reactant1, reactant2, open_el=None, relative_mu=None, use_hull_energy=False)
Gets critical reactions between two reactants.

Get critical reactions (“kinks” in the mixing ratio where
reaction products change) between two reactants. See the
pymatgen.analysis.interface_reactions module for more info.


* **Parameters**


    * **reactant1** (*str*) – Chemical formula for reactant


    * **reactant2** (*str*) – Chemical formula for reactant


    * **open_el** (*str*) – Element in reservoir available to system


    * **relative_mu** (*float*) – Relative chemical potential of element in
    reservoir with respect to pure substance. Must be non-positive.


    * **use_hull_energy** (*bool*) – Whether to use the convex hull energy for a


    * **false** (*given composition for the reaction energy calculation. If*) –


:param :
:param the energy of the ground state structure will be preferred; if a:
:param ground state can not be found for a composition:
:param the convex hull:
:param energy will be used with a warning message.:


* **Returns**

    list of dicts of form {ratio,energy,rxn} where ratio is the

        reactant mixing ratio, energy is the reaction energy
        in eV/atom, and rxn is a
        pymatgen.analysis.reaction_calculator.Reaction.




* **Return type**

    list



#### get_material_id(chemsys_formula)
Get all materials ids for a formula or chemsys.


* **Parameters**

    **chemsys_formula** (*str*) – A chemical system (e.g., Li-Fe-O),
    or formula (e.g., Fe2O3).



* **Returns**

    ([str]) List of all materials ids.



#### get_materials_id_from_task_id(task_id)
Returns a new MP materials id from a task id (which can be
equivalent to an old materials id).


* **Parameters**

    **task_id** (*str*) – A task id.



* **Returns**

    materials_id (str)



#### get_materials_id_references(material_id)
Returns all references for a materials id.


* **Parameters**

    **material_id** (*str*) – A material id.



* **Returns**

    BibTeX (str)



#### get_materials_ids(chemsys_formula)
Get all materials ids for a formula or chemsys.


* **Parameters**

    **chemsys_formula** (*str*) – A chemical system (e.g., Li-Fe-O),
    or formula (e.g., Fe2O3).



* **Returns**

    ([str]) List of all materials ids.



#### get_phonon_bandstructure_by_material_id(material_id: str)
Get phonon dispersion data corresponding to a material_id.


* **Parameters**

    **material_id** (*str*) – Materials Project material_id.



* **Returns**

    A phonon band structure.



* **Return type**

    [PhononBandStructureSymmLine](pymatgen.phonon.md#pymatgen.phonon.bandstructure.PhononBandStructureSymmLine)



#### get_phonon_ddb_by_material_id(material_id: str)
Get ABINIT Derivative Data Base (DDB) output for phonon calculations.


* **Parameters**

    **material_id** (*str*) – Materials Project material_id.



* **Returns**

    ABINIT DDB file as a string.



* **Return type**

    str



#### get_phonon_dos_by_material_id(material_id: str)
Get phonon density of states data corresponding to a material_id.


* **Parameters**

    **material_id** (*str*) – Materials Project material_id.



* **Returns**

    A phonon DOS object.



* **Return type**

    [CompletePhononDos](pymatgen.phonon.md#pymatgen.phonon.dos.CompletePhononDos)



#### get_pourbaix_entries(chemsys, solid_compat='MaterialsProject2020Compatibility')
A helper function to get all entries necessary to generate
a Pourbaix diagram from the rest interface.


* **Parameters**


    * **chemsys** (*str** or **[**str**]*) – Chemical system string comprising element
    symbols separated by dashes, e.g., “Li-Fe-O” or List of element
    symbols, e.g., [“Li”, “Fe”, “O”].


    * **solid_compat** – Compatibility scheme used to pre-process solid DFT energies prior to applying aqueous
    energy adjustments. May be passed as a class (e.g. MaterialsProject2020Compatibility) or an instance
    (e.g., MaterialsProject2020Compatibility()). If None, solid DFT energies are used as-is.
    Default: MaterialsProject2020Compatibility



#### get_reaction(reactants, products)
Gets a reaction from the Materials Project.


* **Parameters**


    * **reactants** (*[**str**]*) – List of formulas


    * **products** (*[**str**]*) – List of formulas



* **Returns**

    rxn



#### get_stability(entries)
Returns the stability of all entries.


#### get_structure_by_material_id(material_id: str, final: bool = True, conventional_unit_cell: bool = False)
Get a Structure corresponding to a material_id.


* **Parameters**


    * **material_id** (*str*) – Materials Project ID (e.g. mp-1234).


    * **final** (*bool*) – Whether to get the final structure, or the initial
    (pre-relaxation) structure. Defaults to True.


    * **conventional_unit_cell** (*bool*) – Whether to get the standard conventional unit cell



* **Returns**

    Structure object.



#### get_structures(chemsys_formula_id, final=True)
Get a list of Structures corresponding to a chemical system, formula,
or materials_id.


* **Parameters**


    * **chemsys_formula_id** (*str*) – A chemical system (e.g., Li-Fe-O),
    or formula (e.g., Fe2O3) or materials_id (e.g., mp-1234).


    * **final** (*bool*) – Whether to get the final structure, or the initial
    (pre-relaxation) structure. Defaults to True.



* **Returns**

    List of Structure objects.



#### get_substrates(material_id, number=50, orient=None)
Get a substrate list for a material id. The list is in order of
increasing elastic energy if a elastic tensor is available for
the material_id. Otherwise the list is in order of increasing
matching area.


* **Parameters**


    * **material_id** (*str*) – Materials Project material_id, e.g. ‘mp-123’.


    * **orient** (*list*) – substrate orientation to look for


    * **number** (*int*) – number of substrates to return
    n=0 returns all available matches



* **Returns**

    list of dicts with substrate matches



#### get_surface_data(material_id, miller_index=None, inc_structures=False)
Gets surface data for a material. Useful for Wulff shapes.

Reference for surface data:

Tran, R., Xu, Z., Radhakrishnan, B., Winston, D., Sun, W., Persson, K.
A., & Ong, S. P. (2016). Data Descriptor: Surface energies of elemental
crystals. Scientific Data, 3(160080), 1-13.
[https://doi.org/10.1038/sdata.2016.80](https://doi.org/10.1038/sdata.2016.80)


* **Parameters**


    * **material_id** (*str*) – Materials Project material_id, e.g. ‘mp-123’.


    * **miller_index** (*list** of **integer*) – The miller index of the surface.


    * **e.g.** –


    * **[****3** –


    * **2** –


    * **provided** (*1**]**. If miller_index is*) –


    * **dictionary** (*only one*) –


    * **returned.** (*of this specific plane will be*) –


    * **inc_structures** (*bool*) – Include final surface slab structures.
    These are unnecessary for Wulff shape construction.



* **Returns**

    Surface data for material. Energies are given in SI units (J/m^2).



#### get_task_data(chemsys_formula_id, prop='')
Flexible method to get any data using the Materials Project REST
interface. Generally used by other methods for more specific queries.
Unlike the

```
:func:`get_data`_
```

, this method queries the task collection
for specific run information.

Format of REST return is *always* a list of dict (regardless of the
number of pieces of data returned. The general format is as follows:

[{“material_id”: material_id, “property_name” : value}, …]


* **Parameters**


    * **chemsys_formula_id** (*str*) – A chemical system (e.g., Li-Fe-O),
    or formula (e.g., Fe2O3) or materials_id (e.g., mp-1234).


    * **prop** (*str*) – Property to be obtained. Should be one of the
    MPRester.supported_properties. Leave as empty string for a
    general list of useful properties.



#### get_wulff_shape(material_id)
Constructs a Wulff shape for a material.


* **Parameters**

    **material_id** (*str*) – Materials Project material_id, e.g. ‘mp-123’.



* **Returns**

    pymatgen.analysis.wulff.WulffShape



#### get_xas_data(material_id, absorbing_element)
Get X-ray absorption spectroscopy data for absorbing element in the
structure corresponding to a material_id. Only X-ray Absorption Near Edge
Structure (XANES) for K-edge is supported.

REST Endpoint:
[https://materialsproject.org/materials](https://materialsproject.org/materials)/<mp-id>/xas/<absorbing_element>.


* **Parameters**


    * **material_id** (*str*) – E.g., mp-1143 for Al2O3


    * **absorbing_element** (*str*) – The absorbing element in the corresponding
    structure. E.g., Al in Al2O3



#### _static_ parse_criteria(criteria_string)
Parses a powerful and simple string criteria and generates a proper
mongo syntax criteria.


* **Parameters**

    **criteria_string** (*str*) – A string representing a search criteria.
    Also supports wild cards. E.g.,
    something like “

    ```
    *
    ```

    2O” gets converted to
    {‘pretty_formula’: {‘$in’: [u’B2O’, u’Xe2O’, u”Li2O”, …]}}

    Other syntax examples:

        mp-1234: Interpreted as a Materials ID.
        Fe2O3 or

        ```
        *
        ```

        2O3: Interpreted as reduced formulas.
        Li-Fe-O or

        ```
        *
        ```

        -Fe-O: Interpreted as chemical systems.

    You can mix and match with spaces, which are interpreted as
    “OR”. E.g., “mp-1234 FeO” means query for all compounds with
    reduced formula FeO or with materials_id mp-1234.




* **Returns**

    A mongo query dict.



#### query(criteria, properties, chunk_size: int = 500, max_tries_per_chunk: int = 5, mp_decode: bool = True, show_progress_bar: bool = True)
Performs an advanced query using MongoDB-like syntax for directly
querying the Materials Project database. This allows one to perform
queries which are otherwise too cumbersome to perform using the standard
convenience methods.

Please consult the Materials API documentation at
[https://github.com/materialsproject/mapidoc](https://github.com/materialsproject/mapidoc), which provides a
comprehensive explanation of the document schema used in the Materials
Project (supported criteria and properties) and guidance on how best to
query for the relevant information you need.

For queries that request data on more than CHUNK_SIZE materials at once,
this method will chunk a query by first retrieving a list of material
IDs that satisfy CRITERIA, and then merging the criteria with a
restriction to one chunk of materials at a time of size CHUNK_SIZE. You
can opt out of this behavior by setting CHUNK_SIZE=0. To guard against
intermittent server errors in the case of many chunks per query,
possibly-transient server errors will result in re-trying a give chunk
up to MAX_TRIES_PER_CHUNK times.


* **Parameters**


    * **criteria** (*str/dict*) – Criteria of the query as a string or
    mongo-style dict.

    If string, it supports a powerful but simple string criteria.
    E.g., “Fe2O3” means search for materials with reduced_formula
    Fe2O3. Wild cards are also supported. E.g., “\\\*2O” means get
    all materials whose formula can be formed as \\\*2O, e.g.,
    Li2O, K2O, etc.

    Other syntax examples:
    mp-1234: Interpreted as a Materials ID.
    Fe2O3 or

    ```
    *
    ```

    2O3: Interpreted as reduced formulas.
    Li-Fe-O or

    ```
    *
    ```

    -Fe-O: Interpreted as chemical systems.

    You can mix and match with spaces, which are interpreted as
    “OR”. E.g. “mp-1234 FeO” means query for all compounds with
    reduced formula FeO or with materials_id mp-1234.

    Using a full dict syntax, even more powerful queries can be
    constructed. For example, {“elements”:{“$in”:[“Li”,
    “Na”, “K”], “$all”: [“O”]}, “nelements”:2} selects all Li, Na
    and K oxides. {“band_gap”: {“$gt”: 1}} selects all materials
    with band gaps greater than 1 eV.



    * **properties** (*list*) – Properties to request for as a list. For
    example, [“formula”, “formation_energy_per_atom”] returns
    the formula and formation energy per atom.


    * **chunk_size** (*int*) – Number of materials for which to fetch data at a
    time. More data-intensive properties may require smaller chunk
    sizes. Use chunk_size=0 to force no chunking – this is useful
    when fetching only properties such as ‘material_id’.


    * **max_tries_per_chunk** (*int*) – How many times to re-try fetching a given
    chunk when the server gives a 5xx error (e.g. a timeout error).


    * **mp_decode** (*bool*) – Whether to do a decoding to a Pymatgen object
    where possible. In some cases, it might be useful to just get
    the raw python dict, i.e., set to False.


    * **show_progress_bar** (*bool*) – Whether to show a progress bar for large queries.
    Defaults to True. Set to False to reduce visual noise.



* **Returns**

    List of results. E.g.,
    [{u’formula’: {u’O’: 1, u’Li’: 2.0}},
    {u’formula’: {u’Na’: 2.0, u’O’: 2.0}},
    {u’formula’: {u’K’: 1, u’O’: 3.0}},
    …]



#### query_snl(criteria)
Query for submitted SNLs.

**NOTE**: As of now, this MP REST feature is open only to a select group of
users. Opening up submissions to all users is being planned for the future.


* **Parameters**

    **criteria** (*dict*) – Query criteria.



* **Returns**

    A dict, with a list of submitted SNLs in the “response” key.



* **Raises**

    **MPRestError** –



#### submit_snl(snl)
Submits a list of StructureNL to the Materials Project site.

**NOTE**: As of now, this MP REST feature is open only to a select group of
users. Opening up submissions to all users is being planned for the future.


* **Parameters**


    * **snl** (*StructureNL/**[*[*StructureNL*](pymatgen.util.md#pymatgen.util.provenance.StructureNL)*]*) – A single StructureNL, or a list


    * **objects** (*of StructureNL*) –



* **Returns**

    A list of inserted submission ids.



* **Raises**

    **MPRestError** –



#### submit_structures(structures, authors, projects=None, references='', remarks=None, data=None, histories=None, created_at=None)
Submits a list of structures to the Materials Project as SNL files.
The argument list mirrors the arguments for the StructureNL object,
except that a list of structures with the same metadata is used as an
input.

**NOTE**: As of now, this MP REST feature is open only to a select group of
users. Opening up submissions to all users is being planned for the future.


* **Parameters**


    * **structures** – A list of Structure objects


    * **authors** (*list*) – List of {“name”:’’, “email”:’’} dicts,
    *list* of Strings as ‘John Doe <[johndoe@gmail.com](mailto:johndoe@gmail.com)>’,
    or a single String with commas separating authors


    * **projects** (*[**str**]*) – List of Strings [‘Project A’, ‘Project B’].
    This applies to all structures.


    * **references** (*str*) – A String in BibTeX format. Again, this applies to
    all structures.


    * **remarks** (*[**str**]*) – List of Strings [‘Remark A’, ‘Remark B’]


    * **data** (*[**dict**]*) – A list of free form dict. Namespaced at the root
    level with an underscore, e.g. {“_materialsproject”:<custom
    data>}. The length of data should be the same as the list of
    structures if not None.


    * **histories** – List of list of dicts - [[{‘name’:’’, ‘url’:’’,
    ‘description’:{}}], …] The length of histories should be the
    same as the list of structures if not None.


    * **created_at** (*datetime*) – A datetime object



* **Returns**

    A list of inserted submission ids.



#### submit_vasp_directory(rootdir, authors, projects=None, references='', remarks=None, master_data=None, master_history=None, created_at=None, ncpus=None)
Assimilates all vasp run directories beneath a particular
directory using BorgQueen to obtain structures, and then submits thhem
to the Materials Project as SNL files. VASP related meta data like
initial structure and final energies are automatically incorporated.

**NOTE**: As of now, this MP REST feature is open only to a select group of
users. Opening up submissions to all users is being planned for the future.


* **Parameters**


    * **rootdir** (*str*) – Rootdir to start assimilating VASP runs from.


    * **authors** – *List* of {“name”:’’, “email”:’’} dicts,
    *list* of Strings as ‘John Doe <[johndoe@gmail.com](mailto:johndoe@gmail.com)>’,
    or a single String with commas separating authors. The same
    list of authors should apply to all runs.


    * **projects** (*[**str**]*) – List of Strings [‘Project A’, ‘Project B’].
    This applies to all structures.


    * **references** (*str*) – A String in BibTeX format. Again, this applies to
    all structures.


    * **remarks** (*[**str**]*) – List of Strings [‘Remark A’, ‘Remark B’]


    * **master_data** (*dict*) – A free form dict. Namespaced at the root
    level with an underscore, e.g. {“_materialsproject”:<custom
    data>}. This data is added to all structures detected in the
    directory, in addition to other vasp data on a per structure
    basis.


    * **master_history** – A master history to be added to all entries.


    * **created_at** (*datetime*) – A datetime object


    * **ncpus** (*int*) – Number of cpus to use in using BorgQueen to
    assimilate. Defaults to None, which means serial.



#### supported_properties(_ = ('energy', 'energy_per_atom', 'volume', 'formation_energy_per_atom', 'nsites', 'unit_cell_formula', 'pretty_formula', 'is_hubbard', 'elements', 'nelements', 'e_above_hull', 'hubbards', 'is_compatible', 'spacegroup', 'task_ids', 'band_gap', 'density', 'icsd_id', 'icsd_ids', 'cif', 'total_magnetization', 'material_id', 'oxide_type', 'tags', 'elasticity'_ )

#### supported_task_properties(_ = ('energy', 'energy_per_atom', 'volume', 'formation_energy_per_atom', 'nsites', 'unit_cell_formula', 'pretty_formula', 'is_hubbard', 'elements', 'nelements', 'e_above_hull', 'hubbards', 'is_compatible', 'spacegroup', 'band_gap', 'density', 'icsd_id', 'cif'_ )

### _class_ _MPResterNewBasic(api_key: str | None = None, include_user_agent: bool = True)
Bases: `object`

A new MPRester that supports the new MP API. If you are getting your API key from the new dashboard of MP, you will
need to use this instead of the original MPRester because the new API keys do not work with the old MP API (???!).
This is a basic implementation for now and features will be added soon. The current implementation is to enable
users with simple requirements use the API without having to install additional packages.

If you are a power user who needs the full functionality, please get the mp-api package instead.


* **Parameters**


    * **api_key** (*str*) – A String API key for accessing the MaterialsProject
    REST interface. Please obtain your API key at
    [https://www.materialsproject.org/dashboard](https://www.materialsproject.org/dashboard). If this is None,
    the code will check if there is a “PMG_MAPI_KEY” setting.
    If so, it will use that environment variable. This makes
    easier for heavy users to simply add this environment variable to
    their setups and MPRester can then be called without any arguments.


    * **include_user_agent** (*bool*) – If True, will include a user agent with the
    HTTP request including information on pymatgen and system version
    making the API request. This helps MP support pymatgen users, and
    is similar to what most web browsers send with each page request.
    Set to False to disable the user agent.



#### get_doc(material_id: str, fields: list | None = None)
Get a data corresponding to a material_id.


* **Parameters**


    * **material_id** (*str*) – Materials Project ID (e.g. mp-1234).


    * **fields** (*list*) – Fields to query for. If None (the default), all fields are returned.



* **Returns**

    Dict



#### get_entries(criteria, compatible_only=True, inc_structure=None, property_data=None, conventional_unit_cell=False, sort_by_e_above_hull=False)
Get a list of ComputedEntries or ComputedStructureEntries corresponding
to a chemical system, formula, or materials_id or full criteria.


* **Parameters**


    * **criteria** – Chemsys, formula, or mp-id.


    * **compatible_only** (*bool*) – Whether to return only “compatible”
    entries. Compatible entries are entries that have been
    processed using the MaterialsProject2020Compatibility class,
    which performs adjustments to allow mixing of GGA and GGA+U
    calculations for more accurate phase diagrams and reaction
    energies.


    * **inc_structure** (*str*) – If None, entries returned are
    ComputedEntries. If inc_structure=”initial”,
    ComputedStructureEntries with initial structures are returned.
    Otherwise, ComputedStructureEntries with final structures
    are returned.


    * **property_data** (*list*) – Specify additional properties to include in
    entry.data. If None, no data. Should be a subset of
    supported_properties.


    * **conventional_unit_cell** (*bool*) – Whether to get the standard
    conventional unit cell


    * **sort_by_e_above_hull** (*bool*) – Whether to sort the list of entries by
    e_above_hull (will query e_above_hull as a property_data if True).



* **Returns**

    List of ComputedStructureEntry objects.



#### get_entries_in_chemsys(elements, \*args, \*\*kwargs)
Helper method to get a list of ComputedEntries in a chemical system. For example, elements = [“Li”, “Fe”, “O”]
will return a list of all entries in the Li-Fe-O chemical system, i.e., all LixOy, FexOy, LixFey, LixFeyOz,
Li, Fe and O phases. Extremely useful for creating phase diagrams of entire chemical systems.


* **Parameters**


    * **elements** (*str** or **[**str**]*) – Chemical system string comprising element
    symbols separated by dashes, e.g., “Li-Fe-O” or List of element
    symbols, e.g., [“Li”, “Fe”, “O”].


    * **\*args** – Pass-through to get_entries.


    * **\*\*kwargs** – Pass-through to get_entries.



* **Returns**

    List of ComputedEntries.



#### get_entry_by_material_id(material_id: str, \*args, \*\*kwargs)
Get a ComputedEntry corresponding to a material_id.


* **Parameters**


    * **material_id** (*str*) – Materials Project material_id (a string,
    e.g., mp-1234).


    * **\*args** – Pass-through to get_entries.


    * **\*\*kwargs** – Pass-through to get_entries.



* **Returns**

    ComputedStructureEntry object.



#### get_initial_structures_by_material_id(material_id: str, conventional_unit_cell: bool = False)
Get a Structure corresponding to a material_id.


* **Parameters**


    * **material_id** (*str*) – Materials Project ID (e.g. mp-1234).


    * **final** (*bool*) – Whether to get the final structure, or the initial
    (pre-relaxation) structures. Defaults to True.


    * **conventional_unit_cell** (*bool*) – Whether to get the standard conventional unit cell



* **Returns**

    Structure object.



#### get_material_ids(formula)
Get all materials ids for a formula.


* **Parameters**

    **formula** (*str*) – A formula (e.g., Fe2O3).



* **Returns**

    ([str]) List of all materials ids.



#### get_materials_ids(formula)
Get all materials ids for a formula.


* **Parameters**

    **formula** (*str*) – A formula (e.g., Fe2O3).



* **Returns**

    ([str]) List of all materials ids.



#### get_structure_by_material_id(material_id: str, conventional_unit_cell: bool = False)
Get a Structure corresponding to a material_id.


* **Parameters**


    * **material_id** (*str*) – Materials Project ID (e.g. mp-1234).


    * **final** (*bool*) – Whether to get the final structure, or the initial
    (pre-relaxation) structures. Defaults to True.


    * **conventional_unit_cell** (*bool*) – Whether to get the standard conventional unit cell



* **Returns**

    Structure object.



#### get_structures(chemsys_formula: str, final=True)
Get a list of Structures corresponding to a chemical system or formula.


* **Parameters**


    * **chemsys_formula** (*str*) – A chemical system, list of chemical systems
    (e.g., Li-Fe-O, Si-*), or single formula (e.g., Fe2O3, Si*).


    * **final** (*bool*) – Whether to get the final structure, or the list of initial
    (pre-relaxation) structures. Defaults to True.



* **Returns**

    List of Structure objects. ([Structure])



#### get_summary(criteria: dict, fields: list | None = None)
Get a data corresponding to a criteria.


* **Parameters**


    * **criteria** (*dict*) – Materials Project ID (e.g. mp-1234), e.g., {“formula”: “Fe2O3,FeO”}


    * **fields** (*list*) – Fields to query for. If None (the default), all fields are returned.



* **Returns**

    List of dict of summary docs.



#### get_summary_by_material_id(material_id: str, fields: list | None = None)
Get a data corresponding to a material_id.


* **Parameters**


    * **material_id** (*str*) – Materials Project ID (e.g. mp-1234).


    * **fields** (*list*) – Fields to query for. If None (the default), all fields are returned.



* **Returns**

    Dict



#### request(sub_url, payload=None, method='GET', mp_decode=True)
Helper method to make the requests and perform decoding based on MSONable protocol.


### get_chunks(sequence: Sequence[Any], size=1)

* **Parameters**


    * **sequence** (*Sequence**[**Any**]*) – Any sequence.


    * **size** (*int*) – Chunk length. Defaults to 1.



* **Returns**

    input sequence in chunks of length size.



* **Return type**

    list[Sequence[Any]]


## pymatgen.ext.optimade module

Optimade support.


### _class_ OptimadeRester(aliases_or_resource_urls: str | list[str] | None = None, refresh_aliases: bool = False, timeout: int = 5)
Bases: `object`

Class to call OPTIMADE-compliant APIs, see [https://optimade.org](https://optimade.org) and [1].

This class is ready to use but considered in-development and subject to change.

[1] Andersen, C.W., *et al*.

    OPTIMADE, an API for exchanging materials data.
    Sci Data 8, 217 (2021). [https://doi.org/10.1038/s41597-021-00974-z](https://doi.org/10.1038/s41597-021-00974-z)

OPTIMADE is an effort to provide a standardized interface to retrieve information
from many different materials science databases.

This is a client to retrieve structures from OPTIMADE v1 compliant endpoints. It
does not yet support all features of the OPTIMADE v1 specification but is intended
as a way to quickly search an endpoint in a way familiar to users of pymatgen without
needing to know the full OPTIMADE specification.

For advanced usage, please see the OPTIMADE documentation at optimade.org and
consider calling the APIs directly.

For convenience, known OPTIMADE endpoints have been given aliases in pymatgen to save
typing the full URL.

To get an up-to-date list aliases, generated from the current list of OPTIMADE providers
at optimade.org, call the refresh_aliases() method or pass refresh_aliases=True when
creating instances of this class.


* **Parameters**


    * **aliases_or_resource_urls** – the alias or structure resource URL or a list of


    * **URLs** (*aliases** or **resource*) –


    * **not** (*if providing the resource URL directly it should*) –


    * **index** (*be an*) –


    * **"v1/structures"** (*this interface can only currently access the*) –


    * **URL** (*information from the specified resource*) –


    * **refresh_aliases** – if True, use an up-to-date list of providers/aliases from the live


    * **https** (*list** of **OPTIMADE providers hosted at*) – //providers.optimade.org.


    * **timeout** – number of seconds before an attempted request is abandoned, a good


    * **providers** (*timeout is useful when querying many*) –


    * **offline** (*some** of **which may be*) –



#### _static_ _build_filter(elements: str | list[str] | None = None, nelements: int | None = None, nsites: int | None = None, chemical_formula_anonymous: str | None = None, chemical_formula_hill: str | None = None)
Convenience method to build an OPTIMADE filter.


#### _get_json(url)
Retrieves JSON, will attempt to (politely) try again on failure subject to a
random delay and a maximum number of attempts.


#### _static_ _get_snls_from_resource(json, url, identifier)

#### _handle_response_fields(additional_response_fields: str | list[str] | set[str] | None = None)
Used internally to handle the mandatory and additional response fields.


* **Parameters**

    **additional_response_fields** – A set of additional fields to request.



* **Returns**

    A string of comma-separated OPTIMADE response fields.



#### _parse_provider(provider, provider_url)
Used internally to update the list of providers or to
check a given URL is valid.

It does not raise exceptions but will instead _logger.warning and provide
an empty dictionary in the case of invalid data.

In future, when the specification  is sufficiently well adopted,
we might be more strict here.


* **Parameters**


    * **provider** – the provider prefix


    * **provider_url** – An OPTIMADE provider URL



* **Returns**

    A dictionary of keys (in format of “provider.database”) to
    Provider objects.



#### _validate_provider(provider_url)
Checks that a given URL is indeed an OPTIMADE provider,
returning None if it is not a provider, or the provider
prefix if it is.

TODO: careful reading of OPTIMADE specification required
TODO: add better exception handling, intentionally permissive currently


#### aliases(_ = {'aflow': 'http://aflow.org/API/optimade/', 'cod': 'https://www.crystallography.net/cod/optimade', 'jarvis': 'https://jarvis.nist.gov/optimade/jarvisdft', 'mcloud.2dtopo': 'https://aiida.materialscloud.org/2dtopo/optimade', 'mcloud.curated-cofs': 'https://aiida.materialscloud.org/curated-cofs/optimade', 'mcloud.mc2d': 'https://aiida.materialscloud.org/mc2d/optimade', 'mcloud.mc3d': 'https://aiida.materialscloud.org/mc3d/optimade', 'mcloud.optimade-sample': 'https://aiida.materialscloud.org/optimade-sample/optimade', 'mcloud.pyrene-mofs': 'https://aiida.materialscloud.org/pyrene-mofs/optimade', 'mcloud.scdm': 'https://aiida.materialscloud.org/autowannier/optimade', 'mcloud.stoceriaitf': 'https://aiida.materialscloud.org/stoceriaitf/optimade', 'mcloud.tc-applicability': 'https://aiida.materialscloud.org/tc-applicability/optimade', 'mcloud.tin-antimony-sulfoiodide': 'https://aiida.materialscloud.org/tin-antimony-sulfoiodide/optimade', 'mp': 'https://optimade.materialsproject.org', 'mpds': 'https://api.mpds.io', 'nmd': 'https://nomad-lab.eu/prod/rae/optimade/', 'odbx': 'https://optimade.odbx.science', 'odbx.odbx_misc': 'https://optimade-misc.odbx.science', 'omdb.omdb_production': 'http://optimade.openmaterialsdb.se', 'oqmd': 'http://oqmd.org/optimade/', 'tcod': 'https://www.crystallography.net/tcod/optimade', 'twodmatpedia': 'http://optimade.2dmatpedia.org'_ )

#### describe()
Provides human-readable information about the resources being searched by the OptimadeRester.


#### get_snls(elements: list[str] | str | None = None, nelements: int | None = None, nsites: int | None = None, chemical_formula_anonymous: str | None = None, chemical_formula_hill: str | None = None, additional_response_fields: str | list[str] | set[str] | None = None)
Retrieve StructureNL from OPTIMADE providers.

A StructureNL is an object provided by pymatgen which combines Structure with
associated metadata, such as the URL is was downloaded from and any additional namespaced
data.

Not all functionality of OPTIMADE is currently exposed in this convenience method. To
use a custom filter, call get_structures_with_filter().


* **Parameters**


    * **elements** – List of elements


    * **nelements** – Number of elements, e.g. 4 or [2, 5] for the range >=2 and <=5


    * **nsites** – Number of sites, e.g. 4 or [2, 5] for the range >=2 and <=5


    * **chemical_formula_anonymous** – Anonymous chemical formula


    * **chemical_formula_hill** – Chemical formula following Hill convention


    * **additional_response_fields** – Any additional fields desired from the OPTIMADE API,


    * **dictionary.** (*these will be stored under the '_optimade' key in each StructureNL.data*) –



* **Returns**

    keyed by that database provider’s id system



* **Return type**

    dict[str, [StructureNL](pymatgen.util.md#pymatgen.util.provenance.StructureNL)]



#### get_snls_with_filter(optimade_filter: str, additional_response_fields: str | list[str] | set[str] | None = None)
Get structures satisfying a given OPTIMADE filter.


* **Parameters**


    * **optimade_filter** – An OPTIMADE-compliant filter


    * **additional_response_fields** – Any additional fields desired from the OPTIMADE API,



* **Returns**

    keyed by that database provider’s id system



* **Return type**

    dict[str, [Structure](pymatgen.core.md#pymatgen.core.structure.Structure)]



#### get_structures(elements: list[str] | str | None = None, nelements: int | None = None, nsites: int | None = None, chemical_formula_anonymous: str | None = None, chemical_formula_hill: str | None = None)
Retrieve Structures from OPTIMADE providers.

Not all functionality of OPTIMADE is currently exposed in this convenience method. To
use a custom filter, call get_structures_with_filter().


* **Parameters**


    * **elements** – List of elements


    * **nelements** – Number of elements, e.g. 4 or [2, 5] for the range >=2 and <=5


    * **nsites** – Number of sites, e.g. 4 or [2, 5] for the range >=2 and <=5


    * **chemical_formula_anonymous** – Anonymous chemical formula


    * **chemical_formula_hill** – Chemical formula following Hill convention



* **Returns**

    keyed by that database provider’s id system



* **Return type**

    dict[str, [Structure](pymatgen.core.md#pymatgen.core.structure.Structure)]



#### get_structures_with_filter(optimade_filter: str)
Get structures satisfying a given OPTIMADE filter.


* **Parameters**

    **optimade_filter** – An OPTIMADE-compliant filter



* **Returns**

    keyed by that database provider’s id system



* **Return type**

    dict[str, [Structure](pymatgen.core.md#pymatgen.core.structure.Structure)]



#### mandatory_response_fields(_ = ('lattice_vectors', 'cartesian_site_positions', 'species', 'species_at_sites'_ )

#### refresh_aliases(providers_url='https://providers.optimade.org/providers.json')
Updates available OPTIMADE structure resources based on the current list of OPTIMADE
providers.


### _class_ Provider(name, base_url, description, homepage, prefix)
Bases: `tuple`

Create new instance of Provider(name, base_url, description, homepage, prefix)


#### _asdict()
Return a new dict which maps field names to their values.


#### _field_defaults(_ = {_ )

#### _fields(_ = ('name', 'base_url', 'description', 'homepage', 'prefix'_ )

#### _classmethod_ _make(iterable)
Make a new Provider object from a sequence or iterable


#### _replace(\*\*kwds)
Return a new Provider object replacing specified fields with new values


#### base_url()
Alias for field number 1


#### description()
Alias for field number 2


#### homepage()
Alias for field number 3


#### name()
Alias for field number 0


#### prefix()
Alias for field number 4