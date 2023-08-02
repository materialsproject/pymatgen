---
layout: default
title: pymatgen.io.cp2k.outputs.md
nav_exclude: true
---

# pymatgen.io.cp2k.outputs module

This module defines the Cp2k output parser along with a few other functions for parsing cp2k-related
outputs.


### _class_ pymatgen.io.cp2k.outputs.Cp2kOutput(filename, verbose=False, auto_load=False)
Bases: `object`

Class for parsing output file from CP2K. The CP2K output file is very flexible in the way that
it is returned. This class will automatically parse parameters that should always be present,
but other parsing features may be called depending on the run type.

Initialize the Cp2kOutput object.


* **Parameters**


    * **filename** – (str) Name of the CP2K output file to parse


    * **verbose** – (bool) Whether or not to parse with verbosity (will parse lots of data that
    may not be useful)


    * **auto_load** (*bool*) – Whether or not to automatically load basic info like energies
    and structures.



#### as_dict()
Return dictionary representation of the output.


#### _property_ band_structure(_: [BandStructure](pymatgen.electronic_structure.bandstructure.md#pymatgen.electronic_structure.bandstructure.BandStructure_ )
Returns band structure object if it has been parsed.


#### _property_ calculation_type()
Returns the calculation type (what io.vasp.outputs calls run_type).


#### _property_ charge(_: floa_ )
Get charge from the input file.


#### _property_ complete_dos(_: [CompleteDos](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.CompleteDos_ )
Returns complete dos object if it has been parsed.


#### _property_ completed()
Did the calculation complete.


#### convergence()
Check whether or not the SCF and geometry optimization cycles converged.


#### _property_ cp2k_version()
The cp2k version used in the calculation.


#### _property_ is_hubbard(_: boo_ )
Returns True if hubbard +U correction was used.


#### _property_ is_metal(_: boo_ )
Was a band gap found? i.e. is it a metal.


#### _property_ is_molecule(_: boo_ )
Returns True if the cp2k output was generated for a molecule (i.e.
no periodicity in the cell). Returns false otherwise.


#### _property_ multiplicity(_: in_ )
Get the spin multiplicity from input file.


#### _property_ num_warnings()
How many warnings showed up during the run.


#### parse_atomic_kind_info()
Parse info on what atomic kinds are present and what basis/pseudopotential is describing
each of them.


#### parse_bandstructure(bandstructure_filename=None)
Parse a CP2K bandstructure file.


* **Parameters**


    * **bandstructure_filename** – Filename containing bandstructure info. If


    * **provided** (*not*) –


    * **by** (*then the pmg name** of **"BAND.bs" will be assumed*) –


    * **parser.** (*the filename*) –



#### parse_cell_params()
Parse the lattice parameters (initial) from the output file.


#### parse_chi_tensor(chi_filename=None)
Parse the magnetic susceptibility tensor.


#### parse_cp2k_params()
Parse the CP2K general parameters from CP2K output file into a dictionary.


#### parse_dft_params()
Parse the DFT parameters (as well as functional, HF, vdW params).


#### parse_dos(dos_file=None, pdos_files=None, ldos_files=None)
Parse the dos files produced by cp2k calculation. CP2K produces different files based
on the input file rather than assimilating them all into one file.

One file type is the overall DOS file, which is used for k-point calculations. For
non-kpoint calculation, the overall DOS is generally not calculated, but the
element-projected pDOS is. Separate files are created for each spin channel and each
atom kind. If requested, cp2k can also do site/local projected dos (ldos). Each site
requested will have a separate file for each spin channel (if spin polarized calculation
is performed).

If possible, this function will assimilate the ldos files into a CompleteDos object.
Either provide a list of PDOS file paths, or use glob to find the .pdos_ALPHA extension
in the calculation directory.


* **Parameters**


    * **dos_file** (*str*) – Name of the dos file, otherwise will be inferred


    * **pdos_files** (*list*) – list of pdos file paths, otherwise they will be inferred


    * **ldos_files** (*list*) – list of ldos file paths, otherwise they will be inferred



#### parse_energies()
Get the total energy from a CP2K calculation. Presently, the energy reported in the
trajectory (pos.xyz) file takes presidence over the energy reported in the main output
file. This is because the trajectory file keeps track of energies in between restarts,
while the main output file may or may not depending on whether a particular machine
overwrites or appends it.


#### parse_files()
Identify files present in the directory with the cp2k output file. Looks for trajectories,
dos, and cubes.


#### parse_forces()
Get the forces from the forces file, or from the main output file.


#### parse_global_params()
Parse the GLOBAL section parameters from CP2K output file into a dictionary.


#### parse_gtensor(gtensor_filename=None)
Parse a file containing g tensor.


#### parse_hirshfeld()
Parse the hirshfeld population analysis for each step.


#### parse_homo_lumo()
Find the HOMO - LUMO gap in [eV]. Returns the last value. For gaps/eigenvalues decomposed
by spin up/spin down channel and over many ionic steps, see parse_mo_eigenvalues().


#### parse_hyperfine(hyperfine_filename=None)
Parse a file containing hyperfine coupling tensors for each atomic site.


#### parse_initial_structure()
Parse the initial structure from the main cp2k output file.


#### parse_input()
Load in the input set from the input file (if it can be found).


#### parse_ionic_steps()
Parse the ionic step info. If already parsed, this will just assimilate.


#### parse_mo_eigenvalues()
Parse the MO eigenvalues from the cp2k output file. Will get the eigenvalues (and band gap)
at each ionic step (if more than one exist).

Everything is decomposed by spin channel. If calculation was performed without spin
polarization, then only Spin.up will be present, which represents the average of up and
down.


#### parse_mulliken()
Parse the mulliken population analysis info for each step
:return:


#### parse_nmr_shift()
Parse NMR calculation.


#### parse_opt_steps()
Parse the geometry optimization information.


#### parse_overlap_condition()
Retrieve the overlap condition number
:return:


#### parse_plus_u_params()
Parse the DFT+U params.


#### parse_qs_params()
Parse the DFT parameters (as well as functional, HF, vdW params).


#### parse_raman()
Parse raman calculation.


#### parse_scf_opt()
Parse the SCF cycles (not usually important).


#### parse_scf_params()
Retrieve the most import SCF parameters: the max number of scf cycles (max_scf),
the convergence cutoff for scf (eps_scf),.


#### parse_stresses()
Get the stresses from stress file, or from the main output file.


#### parse_structures(trajectory_file=None, lattice_file=None)
Parses the structures from a cp2k calculation. Static calculations simply use the initial
structure. For calculations with ionic motion, the function will look for the appropriate
trajectory and lattice files based on naming convention. If no file is given, and no file
is found, it is assumed that the lattice/structure remained constant, and the initial
lattice/structure is used. Cp2k does not output the trajectory in the main output file by
default, so non static calculations have to reference the trajectory file.


#### parse_tddfpt()
Parse TDDFPT calculation.


#### parse_timing()
Parse the timing info (how long did the run take).


#### parse_total_numbers()
Parse total numbers (not usually important).


#### _property_ project_name(_: st_ )
What project name was used for this calculation.


#### ran_successfully()
Sanity checks that the program ran successfully. Looks at the bottom of the CP2K output
file for the “PROGRAM ENDED” line, which is printed when successfully ran. Also grabs
the number of warnings issued.


#### read_pattern(patterns, reverse=False, terminate_on_match=False, postprocess=<class 'str'>)
This function originally comes from pymatgen.io.vasp.outputs Outcar class.

General pattern reading. Uses monty’s regrep method. Takes the same
arguments.


* **Parameters**


    * **patterns** (*dict*) – A dict of patterns, e.g.,
    {“energy”: r”energy\\(sigma->0\\)\\s+=\\s+([\\d\\-.]+)”}.


    * **reverse** (*bool*) – Read files in reverse. Defaults to false. Useful for
    large files, esp OUTCARs, especially when used with
    terminate_on_match.


    * **terminate_on_match** (*bool*) – Whether to terminate when there is at
    least one match in each key in pattern.


    * **postprocess** (*callable*) – A post processing function to convert all
    matches. Defaults to str, i.e., no change.


Renders accessible:

    Any attribute in patterns. For example,
    {“energy”: r”energy\\(sigma->0\\)\\s+=\\s+([\\d\\-.]+)”} will set the
    value of self.data[“energy”] = [[-1234], [-3453], …], to the
    results from regex and postprocess. Note that the returned values
    are lists of lists, because you can grep multiple items on one line.


#### read_table_pattern(header_pattern, row_pattern, footer_pattern, postprocess=<class 'str'>, attribute_name=None, last_one_only=True, strip=None)
This function originally comes from pymatgen.io.vasp.outputs Outcar class.

Parse table-like data. A table composes of three parts: header,
main body, footer. All the data matches “row pattern” in the main body
will be returned.


* **Parameters**


    * **header_pattern** (*str*) – The regular expression pattern matches the
    table header. This pattern should match all the text
    immediately before the main body of the table. For multiple
    sections table match the text until the section of
    interest. MULTILINE and DOTALL options are enforced, as a
    result, the “.” meta-character will also match “n” in this
    section.


    * **row_pattern** (*str*) – The regular expression matches a single line in
    the table. Capture interested field using regular expression
    groups.


    * **footer_pattern** (*str*) – The regular expression matches the end of the
    table. E.g. a long dash line.


    * **postprocess** (*callable*) – A post processing function to convert all
    matches. Defaults to str, i.e., no change.


    * **attribute_name** (*str*) – Name of this table. If present the parsed data
    will be attached to “data. e.g. self.data[“efg”] = […]


    * **last_one_only** (*bool*) – All the tables will be parsed, if this option
    is set to True, only the last table will be returned. The
    enclosing list will be removed. i.e. Only a single table will
    be returned. Default to be True.


    * **strip** (*list*) – Whether or not to strip contents out of the file before
    reading for a table pattern. This is mainly used by parse_scf_opt(),
    to strip HFX info out of the SCF loop start or DFT+U warnings out
    of the SCF loop iterations.



* **Returns**

    List of tables. 1) A table is a list of rows. 2) A row if either a list of
    attribute values in case the the capturing group is defined without name in
    row_pattern, or a dict in case that named capturing groups are defined by
    row_pattern.



#### _property_ run_type()
What type of run (Energy, MD, etc.) was performed.


#### _property_ spin_polarized(_: boo_ )
Was the calculation spin polarized.


### pymatgen.io.cp2k.outputs.parse_dos(dos_file=None)
Parse a dos file. This format is different from the pdos files.


### pymatgen.io.cp2k.outputs.parse_energy_file(energy_file)
Parses energy file for calculations with multiple ionic steps.


### pymatgen.io.cp2k.outputs.parse_pdos(dos_file=None, spin_channel=None, total=False)
Parse a single DOS file created by cp2k. Must contain one PDOS snapshot. i.e. you cannot
use this cannot deal with multiple concatenated dos files.


* **Parameters**


    * **dos_file** (*list*) – list of pdos_ALPHA file paths


    * **spin_channel** (*int*) – Which spin channel the file corresponds to. By default, CP2K will
    write the file with ALPHA or BETA in the filename (for spin up or down), but
    you can specify this here, in case you have a manual file name.
    spin_channel == 1 –> spin up, spin_channel == -1 –> spin down.


    * **total** (*bool*) – Whether to grab the total occupations, or the orbital decomposed ones.


    * **sigma** (*float*) – width for gaussian smearing, if desired



* **Returns**

    >
    > 1. orbital decomposed DOS dict:
    > i.e. pdoss = {specie: {orbital.s: {Spin.up: … }, orbital.px: {Spin.up: … } …}}


    > 2. energy levels of this dos file


    > 3. fermi energy (in eV).

    DOS object is not created here




* **Return type**

    Everything necessary to create a dos object, in dict format