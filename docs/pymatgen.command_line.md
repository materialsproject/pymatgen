---
layout: default
title: pymatgen.command_line.md
nav_exclude: true
---

# pymatgen.command_line package

This package contains various command line wrappers to programs used in
pymatgen that do not have Python equivalents.

## Subpackages


* [pymatgen.command_line.tests package](pymatgen.command_line.tests.md)




    * [pymatgen.command_line.tests.test_bader_caller module](pymatgen.command_line.tests.md#module-pymatgen.command_line.tests.test_bader_caller)


        * [`BaderAnalysisTest`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_bader_caller.BaderAnalysisTest)


            * [`BaderAnalysisTest.setUp()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_bader_caller.BaderAnalysisTest.setUp)


            * [`BaderAnalysisTest.tearDown()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_bader_caller.BaderAnalysisTest.tearDown)


            * [`BaderAnalysisTest.test_atom_parsing()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_bader_caller.BaderAnalysisTest.test_atom_parsing)


            * [`BaderAnalysisTest.test_automatic_runner()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_bader_caller.BaderAnalysisTest.test_automatic_runner)


            * [`BaderAnalysisTest.test_from_path()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_bader_caller.BaderAnalysisTest.test_from_path)


            * [`BaderAnalysisTest.test_init()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_bader_caller.BaderAnalysisTest.test_init)


            * [`BaderAnalysisTest.test_missing_file_bader_exe_path()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_bader_caller.BaderAnalysisTest.test_missing_file_bader_exe_path)


    * [pymatgen.command_line.tests.test_chargemol_caller module](pymatgen.command_line.tests.md#module-pymatgen.command_line.tests.test_chargemol_caller)


        * [`ChargemolAnalysisTest`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_chargemol_caller.ChargemolAnalysisTest)


            * [`ChargemolAnalysisTest.test_parse_chargemol()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_chargemol_caller.ChargemolAnalysisTest.test_parse_chargemol)


            * [`ChargemolAnalysisTest.test_parse_chargemol2()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_chargemol_caller.ChargemolAnalysisTest.test_parse_chargemol2)


    * [pymatgen.command_line.tests.test_critic2_caller module](pymatgen.command_line.tests.md#module-pymatgen.command_line.tests.test_critic2_caller)


        * [`Critic2AnalysisTest`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_critic2_caller.Critic2AnalysisTest)


            * [`Critic2AnalysisTest.setUp()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_critic2_caller.Critic2AnalysisTest.setUp)


            * [`Critic2AnalysisTest.test_graph_output()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_critic2_caller.Critic2AnalysisTest.test_graph_output)


            * [`Critic2AnalysisTest.test_properties_to_from_dict()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_critic2_caller.Critic2AnalysisTest.test_properties_to_from_dict)


        * [`Critic2CallerTest`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_critic2_caller.Critic2CallerTest)


            * [`Critic2CallerTest.test_from_path()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_critic2_caller.Critic2CallerTest.test_from_path)


            * [`Critic2CallerTest.test_from_structure()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_critic2_caller.Critic2CallerTest.test_from_structure)


    * [pymatgen.command_line.tests.test_enumlib_caller module](pymatgen.command_line.tests.md#module-pymatgen.command_line.tests.test_enumlib_caller)


        * [`EnumlibAdaptorTest`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_enumlib_caller.EnumlibAdaptorTest)


            * [`EnumlibAdaptorTest.test_init()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_enumlib_caller.EnumlibAdaptorTest.test_init)


            * [`EnumlibAdaptorTest.test_partial_disorder()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_enumlib_caller.EnumlibAdaptorTest.test_partial_disorder)


            * [`EnumlibAdaptorTest.test_rounding_errors()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_enumlib_caller.EnumlibAdaptorTest.test_rounding_errors)


            * [`EnumlibAdaptorTest.test_timeout()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_enumlib_caller.EnumlibAdaptorTest.test_timeout)


    * [pymatgen.command_line.tests.test_gulp_caller module](pymatgen.command_line.tests.md#module-pymatgen.command_line.tests.test_gulp_caller)


        * [`BuckinghamPotentialBushTest`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_gulp_caller.BuckinghamPotentialBushTest)


            * [`BuckinghamPotentialBushTest.setUp()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_gulp_caller.BuckinghamPotentialBushTest.setUp)


            * [`BuckinghamPotentialBushTest.test_element_different_valence()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_gulp_caller.BuckinghamPotentialBushTest.test_element_different_valence)


            * [`BuckinghamPotentialBushTest.test_existing_element()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_gulp_caller.BuckinghamPotentialBushTest.test_existing_element)


            * [`BuckinghamPotentialBushTest.test_non_exisitng_element()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_gulp_caller.BuckinghamPotentialBushTest.test_non_exisitng_element)


            * [`BuckinghamPotentialBushTest.test_spring()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_gulp_caller.BuckinghamPotentialBushTest.test_spring)


        * [`BuckinghamPotentialLewisTest`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_gulp_caller.BuckinghamPotentialLewisTest)


            * [`BuckinghamPotentialLewisTest.setUp()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_gulp_caller.BuckinghamPotentialLewisTest.setUp)


            * [`BuckinghamPotentialLewisTest.test_element_different_valence()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_gulp_caller.BuckinghamPotentialLewisTest.test_element_different_valence)


            * [`BuckinghamPotentialLewisTest.test_existing_element()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_gulp_caller.BuckinghamPotentialLewisTest.test_existing_element)


            * [`BuckinghamPotentialLewisTest.test_non_exisitng_element()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_gulp_caller.BuckinghamPotentialLewisTest.test_non_exisitng_element)


            * [`BuckinghamPotentialLewisTest.test_spring()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_gulp_caller.BuckinghamPotentialLewisTest.test_spring)


            * [`BuckinghamPotentialLewisTest.test_values()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_gulp_caller.BuckinghamPotentialLewisTest.test_values)


        * [`GlobalFunctionsTest`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_gulp_caller.GlobalFunctionsTest)


            * [`GlobalFunctionsTest.setUp()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_gulp_caller.GlobalFunctionsTest.setUp)


            * [`GlobalFunctionsTest.test_get_energy_buckingham()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_gulp_caller.GlobalFunctionsTest.test_get_energy_buckingham)


            * [`GlobalFunctionsTest.test_get_energy_relax_structure_buckingham()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_gulp_caller.GlobalFunctionsTest.test_get_energy_relax_structure_buckingham)


            * [`GlobalFunctionsTest.test_get_energy_tersoff()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_gulp_caller.GlobalFunctionsTest.test_get_energy_tersoff)


        * [`GulpCallerTest`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_gulp_caller.GulpCallerTest)


            * [`GulpCallerTest.test_decimal()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_gulp_caller.GulpCallerTest.test_decimal)


            * [`GulpCallerTest.test_run()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_gulp_caller.GulpCallerTest.test_run)


        * [`GulpIOTest`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_gulp_caller.GulpIOTest)


            * [`GulpIOTest.setUp()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_gulp_caller.GulpIOTest.setUp)


            * [`GulpIOTest.test_buckingham_input()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_gulp_caller.GulpIOTest.test_buckingham_input)


            * [`GulpIOTest.test_buckingham_potential()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_gulp_caller.GulpIOTest.test_buckingham_potential)


            * [`GulpIOTest.test_get_energy()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_gulp_caller.GulpIOTest.test_get_energy)


            * [`GulpIOTest.test_get_relaxed_structure()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_gulp_caller.GulpIOTest.test_get_relaxed_structure)


            * [`GulpIOTest.test_keyword_line_with_correct_keywords()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_gulp_caller.GulpIOTest.test_keyword_line_with_correct_keywords)


            * [`GulpIOTest.test_library_line_explicit_path()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_gulp_caller.GulpIOTest.test_library_line_explicit_path)


            * [`GulpIOTest.test_library_line_wrong_file()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_gulp_caller.GulpIOTest.test_library_line_wrong_file)


            * [`GulpIOTest.test_specie_potential()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_gulp_caller.GulpIOTest.test_specie_potential)


            * [`GulpIOTest.test_structure_lines_default_options()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_gulp_caller.GulpIOTest.test_structure_lines_default_options)


            * [`GulpIOTest.test_structure_lines_no_frac_coords()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_gulp_caller.GulpIOTest.test_structure_lines_no_frac_coords)


            * [`GulpIOTest.test_structure_lines_no_unitcell()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_gulp_caller.GulpIOTest.test_structure_lines_no_unitcell)


            * [`GulpIOTest.test_tersoff_inpt()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_gulp_caller.GulpIOTest.test_tersoff_inpt)


            * [`GulpIOTest.test_tersoff_potential()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_gulp_caller.GulpIOTest.test_tersoff_potential)


    * [pymatgen.command_line.tests.test_mcsqs_caller module](pymatgen.command_line.tests.md#module-pymatgen.command_line.tests.test_mcsqs_caller)


        * [`McsqsCallerTest`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_mcsqs_caller.McsqsCallerTest)


            * [`McsqsCallerTest.setUp()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_mcsqs_caller.McsqsCallerTest.setUp)


            * [`McsqsCallerTest.test_mcsqs_caller_parallel()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_mcsqs_caller.McsqsCallerTest.test_mcsqs_caller_parallel)


            * [`McsqsCallerTest.test_mcsqs_caller_runtime_error()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_mcsqs_caller.McsqsCallerTest.test_mcsqs_caller_runtime_error)


            * [`McsqsCallerTest.test_mcsqs_caller_supercell()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_mcsqs_caller.McsqsCallerTest.test_mcsqs_caller_supercell)


            * [`McsqsCallerTest.test_mcsqs_caller_total_atoms()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_mcsqs_caller.McsqsCallerTest.test_mcsqs_caller_total_atoms)


            * [`McsqsCallerTest.test_mcsqs_caller_total_atoms_auto_instances()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_mcsqs_caller.McsqsCallerTest.test_mcsqs_caller_total_atoms_auto_instances)


            * [`McsqsCallerTest.test_mcsqs_perfect_match_error()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_mcsqs_caller.McsqsCallerTest.test_mcsqs_perfect_match_error)


            * [`McsqsCallerTest.test_mcsqs_perfect_match_error_parallel()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_mcsqs_caller.McsqsCallerTest.test_mcsqs_perfect_match_error_parallel)


    * [pymatgen.command_line.tests.test_vampire_caller module](pymatgen.command_line.tests.md#module-pymatgen.command_line.tests.test_vampire_caller)


        * [`VampireCallerTest`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_vampire_caller.VampireCallerTest)


            * [`VampireCallerTest.setUp()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_vampire_caller.VampireCallerTest.setUp)


            * [`VampireCallerTest.setUpClass()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_vampire_caller.VampireCallerTest.setUpClass)


            * [`VampireCallerTest.tearDown()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_vampire_caller.VampireCallerTest.tearDown)


            * [`VampireCallerTest.test_vampire()`](pymatgen.command_line.tests.md#pymatgen.command_line.tests.test_vampire_caller.VampireCallerTest.test_vampire)



## pymatgen.command_line.bader_caller module

This module implements an interface to the Henkelmann et al.’s excellent
Fortran code for calculating a Bader charge analysis.

This module depends on a compiled bader executable available in the path.
Please download the library at [http://theory.cm.utexas.edu/henkelman/code/bader/](http://theory.cm.utexas.edu/henkelman/code/bader/)
and follow the instructions to compile the executable.

If you use this module, please cite:

G. Henkelman, A. Arnaldsson, and H. Jonsson, “A fast and robust algorithm for
Bader decomposition of charge density”, Comput. Mater. Sci. 36, 254-360 (2006).


### _class_ pymatgen.command_line.bader_caller.BaderAnalysis(chgcar_filename=None, potcar_filename=None, chgref_filename=None, parse_atomic_densities=False, cube_filename=None, bader_exe_path: str | None = None)
Bases: `object`

Performs Bader analysis for Cube files and VASP outputs.


#### data()
Atomic data parsed from bader analysis. Each dictionary in the list has the keys:
“atomic_vol”, “min_dist”, “charge”, “x”, “y”, “z”.


* **Type**

    list[dict]



#### vacuum_volume()
Vacuum volume of the Bader analysis.


* **Type**

    float



#### vacuum_charge()
Vacuum charge of the Bader analysis.


* **Type**

    float



#### nelectrons()
Number of electrons of the Bader analysis.


* **Type**

    int



#### chgcar()
Chgcar object associated with input CHGCAR file.


* **Type**

    [Chgcar](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Chgcar)



#### atomic_densities()
List of charge densities for each atom centered on the atom.
Excess 0’s are removed from the array to reduce its size. Each dictionary has the keys:
“data”, “shift”, “dim”, where “data” is the charge density array,
“shift” is the shift used to center the atomic charge density, and
“dim” is the dimension of the original charge density map.


* **Type**

    list[dict]


Initializes the Bader caller.


* **Parameters**


    * **chgcar_filename** (*str*) – The filename of the CHGCAR.


    * **potcar_filename** (*str*) – The filename of the POTCAR.


    * **chgref_filename** (*str*) – The filename of the reference charge density.


    * **parse_atomic_densities** (*bool**, **optional*) – turns on atomic partition of the charge density
    charge densities are atom centered


    * **cube_filename** (*str**, **optional*) – The filename of the cube file.


    * **bader_exe_path** (*str**, **optional*) – The path to the bader executable.



#### _classmethod_ from_path(path, suffix='')
Convenient constructor that takes in the path name of VASP run
to perform Bader analysis.


* **Parameters**


    * **path** (*str*) – Name of directory where VASP output files are
    stored.


    * **suffix** (*str*) – specific suffix to look for (e.g. ‘.relax1’
    for ‘CHGCAR.relax1.gz’).



#### get_charge(atom_index)
Convenience method to get the charge on a particular atom. This is the “raw”
charge generated by the Bader program, not a partial atomic charge. If the cube file
is a spin-density file, then this will return the spin density per atom with
positive being spin up and negative being spin down.


* **Parameters**

    **atom_index** – Index of atom.



* **Returns**

    Charge associated with atom from the Bader analysis.



#### get_charge_decorated_structure()
Returns a charge decorated structure.

Note, this assumes that the Bader analysis was correctly performed on a file
with electron densities


#### get_charge_transfer(atom_index, nelect=None)
Returns the charge transferred for a particular atom. A positive value means
that the site has gained electron density (i.e. exhibits anionic character)
whereas a negative value means the site has lost electron density (i.e. exhibits
cationic character). If the arg nelect is not supplied, then POTCAR must be
supplied to determine nelect.


* **Parameters**


    * **atom_index** – Index of atom.


    * **nelect** – number of electrons associated with an isolated atom at this index.
    For most DFT codes this corresponds to the number of valence electrons
    associated with the pseudopotential (e.g. ZVAL for VASP).



* **Returns**

    Charge transfer associated with atom from the Bader analysis.
    Given by bader charge on atom - nelect for associated atom.



#### get_decorated_structure(property_name, average=False)
Get a property-decorated structure from the Bader analysis.

This is distinct from getting charge decorated structure, which assumes
the “standard” Bader analysis of electron densities followed by converting
electron count to charge. The expected way to use this is to call Bader on
a non-charge density file such as a spin density file, electrostatic potential
file, etc., while using the charge density file as the reference (chgref_filename)
so that the partitioning is determined via the charge, but averaging or integrating
is done for another property.

User warning: Bader analysis cannot automatically determine what property is
inside of the file. So if you want to use this for a non-conventional property
like spin, you must ensure that you have the file is for the appropriate
property and you have an appropriate reference file.


* **Parameters**


    * **property_name** – name of the property to assign to the structure, note that
    if name is “spin” this is handled as a special case, and the appropriate
    spin properties are set on the species in the structure


    * **average** – whether or not to return the average of this property, rather
    than the total, by dividing by the atomic volume.



* **Returns**

    structure with site properties assigned via Bader Analysis



#### get_oxidation_state_decorated_structure(nelects=None)
Returns an oxidation state decorated structure based on bader analysis results.
Each site is assigned a charge based on the computed partial atomic charge from bader.

Note, this assumes that the Bader analysis was correctly performed on a file
with electron densities


#### get_partial_charge(atom_index, nelect=None)
Convenience method to get the partial charge on a particular atom. This is
simply the negative value of the charge transferred. A positive value indicates
that the atom has cationic character, whereas a negative value indicates the
site has anionic character.


* **Parameters**


    * **atom_index** – Index of atom.


    * **nelect** – number of electrons associated with an isolated atom at this index.
    For most DFT codes this corresponds to the number of valence electrons
    associated with the pseudopotential (e.g. ZVAL for VASP).



* **Returns**

    Charge associated with atom from the Bader analysis.



#### _property_ summary()
Dict summary of key analysis, e.g., atomic volume, charge, etc.


* **Type**

    return



### pymatgen.command_line.bader_caller.bader_analysis_from_objects(chgcar, potcar=None, aeccar0=None, aeccar2=None)
Convenience method to run Bader analysis from a set
of pymatgen Chgcar and Potcar objects.

This method will:

1. If aeccar objects are present, constructs a temporary reference
file as AECCAR0 + AECCAR2
2. Runs Bader analysis twice: once for charge, and a second time
for the charge difference (magnetization density).


* **Parameters**


    * **chgcar** – Chgcar object


    * **potcar** – (optional) Potcar object


    * **aeccar0** – (optional) Chgcar object from aeccar0 file


    * **aeccar2** – (optional) Chgcar object from aeccar2 file



* **Returns**

    summary dict



### pymatgen.command_line.bader_caller.bader_analysis_from_path(path, suffix='')
Convenience method to run Bader analysis on a folder containing
typical VASP output files.

This method will:

1. Look for files CHGCAR, AECCAR0, AECCAR2, POTCAR or their gzipped
counterparts.
2. If AECCAR\* files are present, constructs a temporary reference
file as AECCAR0 + AECCAR2
3. Runs Bader analysis twice: once for charge, and a second time
for the charge difference (magnetization density).


* **Parameters**


    * **path** – path to folder to search in


    * **suffix** – specific suffix to look for (e.g. ‘.relax1’ for ‘CHGCAR.relax1.gz’



* **Returns**

    summary dict


## pymatgen.command_line.chargemol_caller module

This module implements an interface to Thomas Manz’s Chargemol code
[https://sourceforge.net/projects/ddec](https://sourceforge.net/projects/ddec) for calculating DDEC3, DDEC6, and CM5 population analyses.

This module depends on a compiled chargemol executable being available in the path.
If you use this module, please cite the following based on which modules you use:

Chargemol:
(1) T. A. Manz and N. Gabaldon Limas, Chargemol program for performing DDEC analysis,
Version 3.5, 2017, ddec.sourceforge.net.

DDEC6 Charges:
(1) T. A. Manz and N. Gabaldon Limas, “Introducing DDEC6 atomic population analysis:
part 1. Charge partitioning theory and methodology,” RSC Adv., 6 (2016) 47771-47801.
(2) N. Gabaldon Limas and T. A. Manz, “Introducing DDEC6 atomic population analysis:
part 2. Computed results for a wide range of periodic and nonperiodic materials,”
(3) N. Gabaldon Limas and T. A. Manz, “Introducing DDEC6 atomic population analysis:
part 4. Efficient parallel computation of net atomic charges, atomic spin moments,
bond orders, and more,” RSC Adv., 8 (2018) 2678-2707.

CM5 Charges:
(1) A.V. Marenich, S.V. Jerome, C.J. Cramer, D.G. Truhlar, “Charge Model 5: An Extension
of Hirshfeld Population Analysis for the Accurate Description of Molecular Interactions
in Gaseous and Condensed Phases”, J. Chem. Theory. Comput., 8 (2012) 527-541.

Spin Moments:
(1) T. A. Manz and D. S. Sholl, “Methods for Computing Accurate Atomic Spin Moments for
Collinear and Noncollinear Magnetism in Periodic and Nonperiodic Materials,”
J. Chem. Theory Comput. 7 (2011) 4146-4164.

Bond Orders:
(1) “Introducing DDEC6 atomic population analysis: part 3. Comprehensive method to compute
bond orders,” RSC Adv., 7 (2017) 45552-45581.

DDEC3 Charges:
(1) T. A. Manz and D. S. Sholl, “Improved Atoms-in-Molecule Charge Partitioning Functional
for Simultaneously Reproducing the Electrostatic Potential and Chemical States in Periodic
and Non-Periodic Materials,” J. Chem. Theory Comput. 8 (2012) 2844-2867.
(2) T. A. Manz and D. S. Sholl, “Chemically Meaningful Atomic Charges that Reproduce the
Electrostatic Potential in Periodic and Nonperiodic Materials,” J. Chem. Theory Comput. 6
(2010) 2455-2468.


### _class_ pymatgen.command_line.chargemol_caller.ChargemolAnalysis(path=None, atomic_densities_path=None, run_chargemol=True)
Bases: `object`

Chargemol analysis for DDEC3, DDEC6, and/or CM5 population analyses,
including the calculation of partial atomic charges, atomic spin moments,
bond orders, and related properties.

Initializes the Chargemol Analysis.


* **Parameters**


    * **path** (*str*) – Path to the CHGCAR, POTCAR, AECCAR0, and AECCAR files.


    * **not.** (*Note that it doesn't matter if the files gzip'd or*) – Default: None (current working directory).


    * **atomic_densities_path** (*str**|**None*) – Path to the atomic densities directory


    * **None** (*required by Chargemol. If*) –


    * **is** (*Pymatgen assumes that this*) –


    * **variable.** (*defined in a "DDEC6_ATOMIC_DENSITIES_DIR" environment*) –


    * **True.** (*Only used if run_chargemol is*) – Default: None.


    * **run_chargemol** (*bool*) – Whether to run the Chargemol analysis. If False,


    * **path.** (*the existing Chargemol output files will be read from*) – Default: True.



#### get_bond_order(index_from, index_to)
Convenience method to get the bond order between two atoms.


* **Parameters**


    * **index_from** (*int*) – Index of atom to get bond order from.


    * **index_to** (*int*) – Index of atom to get bond order to.



* **Returns**

    bond order between atoms



* **Return type**

    float



#### get_charge(atom_index, nelect=None, charge_type='ddec')
Convenience method to get the charge on a particular atom using the same
sign convention as the BaderAnalysis. Note that this is *not* the partial
atomic charge. This value is nelect (e.g. ZVAL from the POTCAR) + the
charge transferred. If you want the partial atomic charge, use
get_partial_charge().


* **Parameters**


    * **atom_index** (*int*) – Index of atom to get charge for.


    * **nelect** (*int*) – number of electrons associated with an isolated atom at this index.


    * **electrons** (*For most DFT codes this corresponds to the number** of **valence*) –


    * **None** (*associated with the pseudopotential. If*) –


    * **automatically** (*this value will be*) –


    * **POTCAR** (*obtained from the*) – Default: None.


    * **charge_type** (*str*) – Type of charge to use (“ddec” or “cm5”).



* **Returns**

    charge on atom_index



* **Return type**

    float



#### get_charge_transfer(atom_index, charge_type='ddec')
Returns the charge transferred for a particular atom. A positive value means
that the site has gained electron density (i.e. exhibits anionic character)
whereas a negative value means the site has lost electron density (i.e. exhibits
cationic character). This is the same thing as the negative of the partial atomic
charge.


* **Parameters**


    * **atom_index** (*int*) – Index of atom to get charge transfer for.


    * **charge_type** (*str*) – Type of charge to use (“ddec” or “cm5”).



* **Returns**

    charge transferred at atom_index



* **Return type**

    float



#### get_partial_charge(atom_index, charge_type='ddec')
Convenience method to get the partial atomic charge on a particular atom.
This is the value printed in the Chargemol analysis.


* **Parameters**


    * **atom_index** (*int*) – Index of atom to get charge for.


    * **charge_type** (*str*) – Type of charge to use (“ddec” or “cm5”).



#### get_property_decorated_structure()
Takes CHGCAR’s structure object and updates it with properties
from the Chargemol analysis.


* **Returns**

    Pymatgen structure with site properties added



#### _property_ summary()
Returns a dictionary summary of the Chargemol analysis
{

> “ddec”: {

>     > “partial_charges”: List[float],
>     > “spin_moments”: List[float],
>     > “dipoles”: List[float],
>     > “rsquared_moments”: List[float],
>     > “rcubed_moments”: List[float],
>     > “rfourth_moments”: List[float],
>     > “bond_order_dict”: Dict

>     },

> “cm5”: {

>     > “partial_charges”: List[float],

>     }

}.

## pymatgen.command_line.critic2_caller module

This module implements an interface to the critic2 Bader analysis code.

For most Bader analysis purposes, users are referred to
pymatgen.command_line.bader_caller instead, this module is for advanced
usage requiring identification of critical points in the charge density.

This module depends on a compiled critic2 executable available in the path.
Please follow the instructions at [https://github.com/aoterodelaroza/critic2](https://github.com/aoterodelaroza/critic2)
to compile.

New users are *strongly* encouraged to read the critic2 manual first.

In brief,
\* critic2 searches for critical points in charge density
\* a critical point can be one of four types: nucleus, bond, ring
or cage
\* it does this by seeding locations for likely critical points
and then searching in these regions
\* there are two lists of critical points in the output, a list
of non-equivalent points (with in-depth information about the
field at those points), and a full list of points generated
by the appropriate symmetry operations
\* connectivity between these points is also provided when
appropriate (e.g. the two nucleus critical points linked to

> a bond critical point)


* critic2 can do many other things besides

If you use this module, please cite:

A. Otero-de-la-Roza, E. R. Johnson and V. Luaña,
Comput. Phys. Communications 185, 1007-1018 (2014)
([https://doi.org/10.1016/j.cpc.2013.10.026](https://doi.org/10.1016/j.cpc.2013.10.026))

A. Otero-de-la-Roza, M. A. Blanco, A. Martín Pendás and
V. Luaña, Comput. Phys. Communications 180, 157-166 (2009)
([https://doi.org/10.1016/j.cpc.2008.07.018](https://doi.org/10.1016/j.cpc.2008.07.018))


### _class_ pymatgen.command_line.critic2_caller.Critic2Analysis(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), stdout=None, stderr=None, cpreport=None, yt=None, zpsp=None)
Bases: `MSONable`

Class to process the standard output from critic2 into pymatgen-compatible objects.

This class is used to store results from the Critic2Caller.

To explore the bond graph, use the “structure_graph”
method, which returns a user-friendly StructureGraph
class with bonding information. By default, this returns
a StructureGraph with edge weights as bond lengths, but
can optionally return a graph with edge weights as any
property supported by the CriticalPoint class, such as
bond ellipticity.

This class also provides an interface to explore just the
non-symmetrically-equivalent critical points via the
critical_points attribute, and also all critical
points (via nodes dict) and connections between them
(via edges dict). The user should be familiar with critic2
before trying to understand these.

Indexes of nucleus critical points in the nodes dict are the
same as the corresponding sites in structure, with indices of
other critical points arbitrarily assigned.

Only one of (stdout, cpreport) required, with cpreport preferred
since this is a new, native JSON output from critic2.


* **Parameters**


    * **structure** – associated Structure


    * **stdout** – stdout from running critic2 in automatic
    mode


    * **stderr** – stderr from running critic2 in automatic
    mode


    * **cpreport** – json output from CPREPORT command


    * **yt** – json output from YT command


    * **(****dict****)** (*zpsp*) – Dict of element/symbol name to number of electrons


(ZVAL in VASP pseudopotential), with which to calculate charge transfer.
Optional.


#### get_critical_point_for_site(n: int)

* **Parameters**

    **n** (*int*) – Site index.


Returns: A CriticalPoint instance


#### get_volume_and_charge_for_site(n)

* **Parameters**

    **n** – Site index n.


Returns: A dict containing “volume” and “charge” keys,
or None if YT integration not performed


#### structure_graph(include_critical_points=('bond', 'ring', 'cage'))
A StructureGraph object describing bonding information
in the crystal.


* **Parameters**


    * **include_critical_points** – add DummySpecies for


    * **themselves** (*the critical points*) –


    * **of** (*a list*) –


    * **"nucleus"** –


    * **"bond"** –


    * **"ring"** –


    * **"cage"** –


    * **None** (*set to*) –


    * **disable** (*to*) –


Returns: a StructureGraph


### _class_ pymatgen.command_line.critic2_caller.Critic2Caller(input_script)
Bases: `object`

Class to call critic2 and store standard output for further processing.

Run Critic2 on a given input script.


* **Parameters**

    **input_script** – string defining the critic2 input



#### _classmethod_ from_chgcar(structure, chgcar=None, chgcar_ref=None, user_input_settings=None, write_cml=False, write_json=True, zpsp=None)
Run Critic2 in automatic mode on a supplied structure, charge
density (chgcar) and reference charge density (chgcar_ref).

The reason for a separate reference field is that in
VASP, the CHGCAR charge density only contains valence
electrons and may be missing substantial charge at
nuclei leading to misleading results. Thus, a reference
field is commonly constructed from the sum of AECCAR0
and AECCAR2 which is the total charge density, but then
the valence charge density is used for the final analysis.

If chgcar_ref is not supplied, chgcar will be used as the
reference field. If chgcar is not supplied, the promolecular
charge density will be used as the reference field – this can
often still give useful results if only topological information
is wanted.

User settings is a dictionary that can contain:
\* GRADEPS, float (field units), gradient norm threshold
\* CPEPS, float (Bohr units in crystals), minimum distance between

> critical points for them to be equivalent


* NUCEPS, same as CPEPS but specifically for nucleus critical
points (critic2 default is dependent on grid dimensions)


* NUCEPSH, same as NUCEPS but specifically for hydrogen nuclei
since associated charge density can be significantly displaced
from hydrogen nucleus


* EPSDEGEN, float (field units), discard critical point if any
element of the diagonal of the Hessian is below this value,
useful for discarding points in vacuum regions


* DISCARD, float (field units), discard critical points with field
value below this value, useful for discarding points in vacuum
regions


* SEED, list of strings, strategies for seeding points, default
is [‘WS 1’, ‘PAIR 10’] which seeds critical points by
sub-dividing the Wigner-Seitz cell and between every atom pair
closer than 10 Bohr, see critic2 manual for more options


* **Parameters**


    * **structure** – Structure to analyze


    * **chgcar** – Charge density to use for analysis. If None, will
    use promolecular density. Should be a Chgcar object or path (string).


    * **chgcar_ref** – Reference charge density. If None, will use
    chgcar as reference. Should be a Chgcar object or path (string).


    * **(****dict****)** (*user_input_settings*) – as explained above


    * **(****bool****)** (*write_json*) – Useful for debug, if True will write all
    critical points to a file ‘table.cml’ in the working directory
    useful for visualization


    * **(****bool****)** – Whether to write out critical points


and YT json. YT integration will be performed with this setting.
:param zpsp (dict): Dict of element/symbol name to number of electrons
(ZVAL in VASP pseudopotential), with which to properly augment core regions
and calculate charge transfer. Optional.


#### _classmethod_ from_path(path, suffix='', zpsp=None)
Convenience method to run critic2 analysis on a folder with typical VASP output files.

This method will:

1. Look for files CHGCAR, AECAR0, AECAR2, POTCAR or their gzipped
counterparts.

2. If AECCAR\* files are present, constructs a temporary reference
file as AECCAR0 + AECCAR2.

3. Runs critic2 analysis twice: once for charge, and a second time
for the charge difference (magnetization density).


* **Parameters**


    * **path** – path to folder to search in


    * **suffix** – specific suffix to look for (e.g. ‘.relax1’ for
    ‘CHGCAR.relax1.gz’)


    * **zpsp** – manually specify ZPSP if POTCAR not present



* **Returns**




### _class_ pymatgen.command_line.critic2_caller.CriticalPoint(index, type, frac_coords, point_group, multiplicity, field, field_gradient, coords=None, field_hessian=None)
Bases: `MSONable`

Access information about a critical point and the field values at that point.

Class to characterise a critical point from a topological
analysis of electron charge density.

Note this class is usually associated with a Structure, so
has information on multiplicity/point group symmetry.


* **Parameters**


    * **index** – index of point


    * **type** – type of point, given as a string


    * **coords** – Cartesian coordinates in Angstroms


    * **frac_coords** – fractional coordinates


    * **point_group** – point group associated with critical point


    * **multiplicity** – number of equivalent critical points


    * **field** – value of field at point (f)


    * **field_gradient** – gradient of field at point (grad f)


    * **field_hessian** – hessian of field at point (del^2 f)



#### _property_ ellipticity()
Most meaningful for bond critical points,
can be physically interpreted as e.g. degree
of pi-bonding in organic molecules. Consult
literature for more information.
Returns: The ellpiticity of the field at the critical point.


#### _property_ laplacian()
The Laplacian of the field at the critical point.


* **Type**

    Returns



#### _property_ type()
Instance of CriticalPointType.


* **Type**

    Returns



### _class_ pymatgen.command_line.critic2_caller.CriticalPointType(value)
Bases: `Enum`

Enum type for the different varieties of critical point.


#### bond(_ = 'bond_ )

#### cage(_ = 'cage_ )

#### nnattr(_ = 'nnattr_ )

#### nucleus(_ = 'nucleus_ )

#### ring(_ = 'ring_ )

### pymatgen.command_line.critic2_caller.get_filepath(filename, warning, path, suffix)

* **Parameters**


    * **filename** – Filename


    * **warning** – Warning message


    * **path** – Path to search


    * **suffix** – Suffixes to search.


## pymatgen.command_line.enumlib_caller module

This module implements an interface to enumlib, Gus Hart’s excellent Fortran
code for enumerating derivative structures.

This module depends on a compiled enumlib with the executables enum.x and
makestr.x available in the path. Please download the library at
[https://github.com/msg-byu/enumlib](https://github.com/msg-byu/enumlib) and follow the instructions in the README to
compile these two executables accordingly.

If you use this module, please cite:

Gus L. W. Hart and Rodney W. Forcade, “Algorithm for generating derivative
structures,” Phys. Rev. B 77 224115 (26 June 2008)

Gus L. W. Hart and Rodney W. Forcade, “Generating derivative structures from
multilattices: Application to hcp alloys,” Phys. Rev. B 80 014120 (July 2009)

Gus L. W. Hart, Lance J. Nelson, and Rodney W. Forcade, “Generating
derivative structures at a fixed concentration,” Comp. Mat. Sci. 59
101-107 (March 2012)

Wiley S. Morgan, Gus L. W. Hart, Rodney W. Forcade, “Generating derivative
superstructures for systems with high configurational freedom,” Comp. Mat.
Sci. 136 144-149 (May 2017)


### _exception_ pymatgen.command_line.enumlib_caller.EnumError()
Bases: `BaseException`

Error subclass for enumeration errors.

## pymatgen.command_line.gulp_caller module

Interface with command line GULP.
[http://projects.ivec.org](http://projects.ivec.org)
WARNING: you need to have GULP installed on your system.


### _class_ pymatgen.command_line.gulp_caller.BuckinghamPotential(bush_lewis_flag)
Bases: `object`

Generate the Buckingham Potential Table from the bush.lib and lewis.lib.

Ref:
T.S.Bush, J.D.Gale, C.R.A.Catlow and P.D. Battle,  J. Mater Chem.,
4, 831-837 (1994).
G.V. Lewis and C.R.A. Catlow, J. Phys. C: Solid State Phys., 18,
1149-1161 (1985)


* **Parameters**

    **bush_lewis_flag** (*str*) – Flag for using Bush or Lewis potential.



### _class_ pymatgen.command_line.gulp_caller.GulpCaller(cmd='gulp')
Bases: `object`

Class to run gulp from commandline.

Initialize with the executable if not in the standard path.


* **Parameters**

    **cmd** – Command. Defaults to gulp.



#### run(gin)
Run GULP using the gin as input.


* **Parameters**

    **gin** – GULP input string



* **Returns**

    GULP output string



* **Return type**

    gout



### _exception_ pymatgen.command_line.gulp_caller.GulpConvergenceError(msg='')
Bases: `Exception`

Exception class for GULP.
Raised when proper convergence is not reached in Mott-Littleton
defect energy optimization procedure in GULP.


* **Parameters**

    **msg** (*str*) – Message.



### _exception_ pymatgen.command_line.gulp_caller.GulpError(msg)
Bases: `Exception`

Exception class for GULP.
Raised when the GULP gives an error.


* **Parameters**

    **msg** (*str*) – Message.



### _class_ pymatgen.command_line.gulp_caller.GulpIO()
Bases: `object`

To generate GULP input and process output.


#### buckingham_input(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), keywords, library=None, uc=True, valence_dict=None)
Gets a GULP input for an oxide structure and buckingham potential
from library.


* **Parameters**


    * **structure** – pymatgen.core.structure.Structure


    * **keywords** – GULP first line keywords.


    * **library** (*Default=None*) – File containing the species and potential.


    * **uc** (*Default=True*) – Unit Cell Flag.


    * **valence_dict** – {El: valence}



#### _static_ buckingham_potential(structure, val_dict=None)
Generate species, buckingham, and spring options for an oxide structure
using the parameters in default libraries.

Ref:


    1. G.V. Lewis and C.R.A. Catlow, J. Phys. C: Solid State Phys.,
    18, 1149-1161 (1985)


    2. T.S.Bush, J.D.Gale, C.R.A.Catlow and P.D. Battle,
    J. Mater Chem., 4, 831-837 (1994)


* **Parameters**


    * **structure** – pymatgen.core.structure.Structure


    * **val_dict** (*Needed if structure is not charge neutral*) – {El:valence}
    dict, where El is element.



#### _static_ get_energy(gout: str)

* **Parameters**

    **gout** (*str*) – GULP output string.



* **Returns**

    Energy



#### _static_ get_relaxed_structure(gout: str)

* **Parameters**

    **gout** (*str*) – GULP output string.



* **Returns**

    (Structure) relaxed structure.



#### _static_ keyword_line(\*args)
Checks if the input args are proper gulp keywords and
generates the 1st line of gulp input. Full keywords are expected.


* **Parameters**

    **args** – 1st line keywords



#### _static_ library_line(file_name)
Specifies GULP library file to read species and potential parameters.
If using library don’t specify species and potential
in the input file and vice versa. Make sure the elements of
structure are in the library file.


* **Parameters**

    **file_name** – Name of GULP library file



* **Returns**

    GULP input string specifying library option



#### _static_ specie_potential_lines(structure, potential, \*\*kwargs)
Generates GULP input specie and potential string for pymatgen
structure.


* **Parameters**


    * **structure** – pymatgen.core.structure.Structure object


    * **potential** – String specifying the type of potential used


    * **kwargs** – Additional parameters related to potential. For
    potential == “buckingham”,
    anion_shell_flg (default = False):
    If True, anions are considered polarizable.
    anion_core_chrg=float
    anion_shell_chrg=float
    cation_shell_flg (default = False):
    If True, cations are considered polarizable.
    cation_core_chrg=float
    cation_shell_chrg=float



* **Returns**

    string containing specie and potential specification for gulp
    input.



#### _static_ structure_lines(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), cell_flg: bool = True, frac_flg: bool = True, anion_shell_flg: bool = True, cation_shell_flg: bool = False, symm_flg: bool = True)
Generates GULP input string corresponding to pymatgen structure.


* **Parameters**


    * **structure** – pymatgen Structure object


    * **cell_flg** (*default = True*) – Option to use lattice parameters.


    * **frac_flg** (*default = True*) – If True, fractional coordinates
    are used. Else, Cartesian coordinates in Angstroms are used.
    **\*\***
    GULP convention is to use fractional coordinates for periodic
    structures and Cartesian coordinates for non-periodic
    structures.
    **\*\***


    * **anion_shell_flg** (*default = True*) – If True, anions are considered
    polarizable.


    * **cation_shell_flg** (*default = False*) – If True, cations are
    considered polarizable.


    * **symm_flg** (*default = True*) – If True, symmetry information is also
    written.



* **Returns**

    string containing structure for GULP input



#### tersoff_input(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), periodic=False, uc=True, \*keywords)
Gets a GULP input with Tersoff potential for an oxide structure.


* **Parameters**


    * **structure** – pymatgen.core.structure.Structure


    * **periodic** (*Default=False*) – Flag denoting whether periodic
    boundary conditions are used


    * **library** (*Default=None*) – File containing the species and potential.


    * **uc** (*Default=True*) – Unit Cell Flag.


    * **keywords** – GULP first line keywords.



#### _static_ tersoff_potential(structure)
Generate the species, Tersoff potential lines for an oxide structure.


* **Parameters**

    **structure** – pymatgen.core.structure.Structure



### _class_ pymatgen.command_line.gulp_caller.TersoffPotential()
Bases: `object`

Generate Tersoff Potential Table from “OxideTersoffPotentialentials” file.

Init TersoffPotential.


### pymatgen.command_line.gulp_caller.get_energy_buckingham(structure, gulp_cmd='gulp', keywords=('optimise', 'conp', 'qok'), valence_dict=None)
Compute the energy of a structure using Buckingham potential.


* **Parameters**


    * **structure** – pymatgen.core.structure.Structure


    * **gulp_cmd** – GULP command if not in standard place


    * **keywords** – GULP first line keywords


    * **valence_dict** – {El: valence}. Needed if the structure is not charge
    neutral.



### pymatgen.command_line.gulp_caller.get_energy_relax_structure_buckingham(structure, gulp_cmd='gulp', keywords=('optimise', 'conp'), valence_dict=None)
Relax a structure and compute the energy using Buckingham potential.


* **Parameters**


    * **structure** – pymatgen.core.structure.Structure


    * **gulp_cmd** – GULP command if not in standard place


    * **keywords** – GULP first line keywords


    * **valence_dict** – {El: valence}. Needed if the structure is not charge
    neutral.



### pymatgen.command_line.gulp_caller.get_energy_tersoff(structure, gulp_cmd='gulp')
Compute the energy of a structure using Tersoff potential.


* **Parameters**


    * **structure** – pymatgen.core.structure.Structure


    * **gulp_cmd** – GULP command if not in standard place


## pymatgen.command_line.mcsqs_caller module

Module to call mcsqs, distributed with AT-AT
[https://www.brown.edu/Departments/Engineering/Labs/avdw/atat/](https://www.brown.edu/Departments/Engineering/Labs/avdw/atat/).


### _class_ pymatgen.command_line.mcsqs_caller.Sqs(bestsqs, objective_function, allsqs, clusters, directory)
Bases: `tuple`

Return type for run_mcsqs.
bestsqs: Structure
objective_function: float | str
allsqs: List
clusters: List
directory: str


#### allsqs()
Alias for field number 2


#### bestsqs()
Alias for field number 0


#### clusters()
Alias for field number 3


#### directory()
Alias for field number 4


#### objective_function()
Alias for field number 1


### pymatgen.command_line.mcsqs_caller.run_mcsqs(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), clusters: dict[int, float], scaling: int | list[int] = 1, search_time: float = 60, directory: str | None = None, instances: int | None = None, temperature: int | float = 1, wr: float = 1, wn: float = 1, wd: float = 0.5, tol: float = 0.001)
Helper function for calling mcsqs with different arguments
:param structure: Disordered pymatgen Structure object
:type structure: Structure
:param clusters: Dictionary of cluster interactions with entries in the form

> number of atoms: cutoff in angstroms


* **Parameters**


    * **scaling** (*int** or **list*) – Scaling factor to determine supercell. Two options are possible:


            1. (preferred) Scales number of atoms, e.g., for a structure with 8 atoms,
        scaling=4 would lead to a 32 atom supercell


            2. A sequence of three scaling factors, e.g., [2, 1, 1], which
        specifies that the supercell should have dimensions 2a x b x c

    Defaults to 1.



    * **search_time** (*float*) – Time spent looking for the ideal SQS in minutes (default: 60)


    * **directory** (*str*) – Directory to run mcsqs calculation and store files (default: None
    runs calculations in a temp directory)


    * **instances** (*int*) – Specifies the number of parallel instances of mcsqs to run
    (default: number of cpu cores detected by Python)


    * **temperature** (*int** or **float*) – Monte Carlo temperature (default: 1), “T” in atat code


    * **wr** (*int** or **float*) – Weight assigned to range of perfect correlation match in objective
    function (default = 1)


    * **wn** (*int** or **float*) – Multiplicative decrease in weight per additional point in cluster (default: 1)


    * **wd** (*int** or **float*) – Exponent of decay in weight as function of cluster diameter (default: 0.5)


    * **tol** (*int** or **float*) – Tolerance for matching correlations (default: 1e-3).



* **Returns**

    Tuple of Pymatgen structure SQS of the input structure, the mcsqs objective function,

        list of all SQS structures, and the directory where calculations are run



## pymatgen.command_line.vampire_caller module

This module implements an interface to the VAMPIRE code for atomistic
simulations of magnetic materials.

This module depends on a compiled vampire executable available in the path.
Please download at [https://vampire.york.ac.uk/download/](https://vampire.york.ac.uk/download/) and
follow the instructions to compile the executable.

If you use this module, please cite:

“Atomistic spin model simulations of magnetic nanomaterials.”
R. F. L. Evans, W. J. Fan, P. Chureemart, T. A. Ostler, M. O. A. Ellis
and R. W. Chantrell. J. Phys.: Condens. Matter 26, 103202 (2014)


### _class_ pymatgen.command_line.vampire_caller.VampireCaller(ordered_structures=None, energies=None, mc_box_size=4.0, equil_timesteps=2000, mc_timesteps=4000, save_inputs=False, hm=None, avg=True, user_input_settings=None)
Bases: `object`

Run Vampire on a material with magnetic ordering and exchange parameter information to compute the critical
temperature with classical Monte Carlo.

user_input_settings is a dictionary that can contain:
\* start_t (int): Start MC sim at this temp, defaults to 0 K.
\* end_t (int): End MC sim at this temp, defaults to 1500 K.
\* temp_increment (int): Temp step size, defaults to 25 K.


* **Parameters**


    * **ordered_structures** (*list*) – Structure objects with magmoms.


    * **energies** (*list*) – Energies of each relaxed magnetic structure.


    * **mc_box_size** (*float*) – x=y=z dimensions (nm) of MC simulation box


    * **equil_timesteps** (*int*) – number of MC steps for equilibrating


    * **mc_timesteps** (*int*) – number of MC steps for averaging


    * **save_inputs** (*bool*) – if True, save scratch dir of vampire input files


    * **hm** ([*HeisenbergModel*](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.heisenberg.HeisenbergModel)) – object already fit to low energy
    magnetic orderings.


    * **avg** (*bool*) – If True, simply use <J> exchange parameter estimate.
    If False, attempt to use NN, NNN, etc. interactions.


    * **user_input_settings** (*dict*) – optional commands for VAMPIRE Monte Carlo


    * **sgraph** ([*StructureGraph*](pymatgen.analysis.md#pymatgen.analysis.graphs.StructureGraph)) – Ground state graph.


    * **unique_site_ids** (*dict*) – Maps each site to its unique identifier


    * **nn_interactions** (*dict*) – {i: j} pairs of NN interactions
    between unique sites.


    * **ex_params** (*dict*) – Exchange parameter values (meV/atom)


    * **mft_t** (*float*) – Mean field theory estimate of critical T


    * **mat_name** (*str*) – Formula unit label for input files


    * **mat_id_dict** (*dict*) – Maps sites to material id # for vampire
    indexing.



#### _static_ parse_stdout(vamp_stdout, n_mats: int)
Parse stdout from Vampire.


* **Parameters**


    * **vamp_stdout** (*txt file*) – Vampire ‘output’ file.


    * **n_mats** (*int*) – Number of materials in Vampire simulation.



* **Returns**

    MSONable vampire output.
    critical_temp (float): Calculated critical temp.



* **Return type**

    parsed_out (DataFrame)



### _class_ pymatgen.command_line.vampire_caller.VampireOutput(parsed_out=None, nmats=None, critical_temp=None)
Bases: `MSONable`

This class processes results from a Vampire Monte Carlo simulation
and returns the critical temperature.


* **Parameters**


    * **parsed_out** (*json*) – json rep of parsed stdout DataFrame.


    * **nmats** (*int*) – Number of distinct materials (1 for each specie and up/down spin).


    * **critical_temp** (*float*) – Monte Carlo Tc result.