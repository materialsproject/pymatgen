---
layout: default
title: pymatgen.io.vasp.outputs.md
nav_exclude: true
---

# pymatgen.io.vasp.outputs module

Classes for reading/manipulating/writing VASP output files.


### _class_ pymatgen.io.vasp.outputs.BSVasprun(filename: str, parse_projected_eigen: bool | str = False, parse_potcar_file: bool | str = False, occu_tol: float = 1e-08, separate_spins: bool = False)
Bases: `Vasprun`

A highly optimized version of Vasprun that parses only eigenvalues for
bandstructures. All other properties like structures, parameters,
etc. are ignored.


* **Parameters**


    * **filename** – Filename to parse


    * **parse_projected_eigen** – Whether to parse the projected
    eigenvalues. Defaults to False. Set to True to obtain projected
    eigenvalues. **Note that this can take an extreme amount of time
    and memory.** So use this wisely.


    * **parse_potcar_file** – Whether to parse the potcar file to read
    the potcar hashes for the potcar_spec attribute. Defaults to True,
    where no hashes will be determined and the potcar_spec dictionaries
    will read {“symbol”: ElSymbol, “hash”: None}. By Default, looks in
    the same directory as the vasprun.xml, with same extensions as

    > Vasprun.xml. If a string is provided, looks at that filepath.



    * **occu_tol** – Sets the minimum tol for the determination of the
    vbm and cbm. Usually the default of 1e-8 works well enough,
    but there may be pathological cases.


    * **separate_spins** (*bool*) – Whether the band gap, CBM, and VBM should be
    reported for each individual spin channel. Defaults to False,
    which computes the eigenvalue band properties independent of
    the spin orientation. If True, the calculation must be spin-polarized.



#### as_dict()
JSON-serializable dict representation.


### _class_ pymatgen.io.vasp.outputs.Chgcar(poscar, data, data_aug=None)
Bases: `VolumetricData`

Simple object for reading a CHGCAR file.


* **Parameters**


    * **poscar** ([*Poscar*](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Poscar)* or *[*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Object containing structure.


    * **data** – Actual data.


    * **data_aug** – Augmentation charge data.



#### _static_ from_file(filename: str)
Read a CHGCAR file.


* **Parameters**

    **filename** – Filename



* **Returns**

    Chgcar



#### _property_ net_magnetization()
Net magnetization from Chgcar


* **Type**

    return



### _class_ pymatgen.io.vasp.outputs.Dynmat(filename)
Bases: `object`

Object for reading a DYNMAT file.


#### data()
A nested dict containing the DYNMAT data of the form::
[atom <int>][disp <int>][‘dispvec’] =

> displacement vector (part of first line in dynmat block, e.g. “0.01 0 0”)

[atom <int>][disp <int>][‘dynmat’] =

    <list> list of dynmat lines for this atom and this displacement

Authors: Patrick Huck


* **Parameters**

    **filename** – Name of file containing DYNMAT.



#### get_phonon_frequencies()
Calculate phonon frequencies.


#### _property_ masses()
Returns the list of atomic masses.


#### _property_ natoms()
Returns the number of atoms.


#### _property_ ndisps()
Returns the number of displacements.


#### _property_ nspecs()
Returns the number of species.


### _class_ pymatgen.io.vasp.outputs.Eigenval(filename, occu_tol=1e-08, separate_spins=False)
Bases: `object`

Object for reading EIGENVAL file.


#### filename()
string containing input filename


#### occu_tol()
tolerance for determining occupation in band properties


#### ispin()
spin polarization tag (int)


#### nelect()
number of electrons


#### nkpt()
number of kpoints


#### nbands()
number of bands


#### kpoints()
list of kpoints


#### kpoints_weights()
weights of each kpoint in the BZ, should sum to 1.


#### eigenvalues()
Eigenvalues as a dict of {(spin): np.ndarray(shape=(nkpt, nbands, 2))}.
This representation is based on actual ordering in VASP and is meant as
an intermediate representation to be converted into proper objects. The
kpoint index is 0-based (unlike the 1-based indexing in VASP).

Reads input from filename to construct Eigenval object.


* **Parameters**


    * **filename** (*str*) – filename of EIGENVAL to read in


    * **occu_tol** (*float*) – tolerance for determining band gap


    * **separate_spins** (*bool*) – whether the band gap, CBM, and VBM should be
    reported for each individual spin channel. Defaults to False,
    which computes the eigenvalue band properties independent of
    the spin orientation. If True, the calculation must be spin-polarized.



* **Returns**

    a pymatgen.io.vasp.outputs.Eigenval object



#### _property_ eigenvalue_band_properties()
Band properties from the eigenvalues as a tuple,
(band gap, cbm, vbm, is_band_gap_direct). In the case of separate_spins=True,
the band gap, cbm, vbm, and is_band_gap_direct are each lists of length 2,
with index 0 representing the spin-up channel and index 1 representing
the spin-down channel.


### _class_ pymatgen.io.vasp.outputs.Elfcar(poscar, data)
Bases: `VolumetricData`

Read an ELFCAR file which contains the Electron Localization Function (ELF)
as calculated by VASP.

For ELF, “total” key refers to Spin.up, and “diff” refers to Spin.down.

This also contains information on the kinetic energy density.


* **Parameters**


    * **poscar** ([*Poscar*](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Poscar)* or *[*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Object containing structure.


    * **data** – Actual data.



#### _classmethod_ from_file(filename)
Reads a ELFCAR file.


* **Parameters**

    **filename** – Filename



* **Returns**

    Elfcar



#### get_alpha()
Get the parameter alpha where ELF = 1/(1+alpha^2).


### _class_ pymatgen.io.vasp.outputs.Locpot(poscar, data)
Bases: `VolumetricData`

Simple object for reading a LOCPOT file.


* **Parameters**


    * **poscar** ([*Poscar*](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Poscar)) – Poscar object containing structure.


    * **data** – Actual data.



#### _classmethod_ from_file(filename, \*\*kwargs)
Reads a LOCPOT file.


* **Parameters**

    **filename** – Filename



* **Returns**

    Locpot



### _class_ pymatgen.io.vasp.outputs.Oszicar(filename)
Bases: `object`

A basic parser for an OSZICAR output from VASP. In general, while the
OSZICAR is useful for a quick look at the output from a VASP run, we
recommend that you use the Vasprun parser instead, which gives far richer
information about a run.


#### electronic_steps()
All electronic steps as a list of list of dict. e.g.,
[[{“rms”: 160.0, “E”: 4507.24605593, “dE”: 4507.2, “N”: 1,
“deps”: -17777.0, “ncg”: 16576}, …], [….]
where electronic_steps[index] refers the list of electronic steps
in one ionic_step, electronic_steps[index][subindex] refers to a
particular electronic step at subindex in ionic step at index. The
dict of properties depends on the type of VASP run, but in general,
“E”, “dE” and “rms” should be present in almost all runs.


### ionic_steps:()
All ionic_steps as a list of dict, e.g.,
[{“dE”: -526.36, “E0”: -526.36024, “mag”: 0.0, “F”: -526.36024},
…]
This is the typical output from VASP at the end of each ionic step.


* **Parameters**

    **filename** (*str*) – Filename of file to parse.



#### _property_ all_energies()
Compilation of all energies from all electronic steps and ionic steps
as a tuple of list of energies, e.g.,
((4507.24605593, 143.824705755, -512.073149912, …), …).


#### as_dict()

* **Returns**

    MSONable dict



#### _property_ final_energy()

### _class_ pymatgen.io.vasp.outputs.Outcar(filename)
Bases: `object`

Parser for data in OUTCAR that is not available in Vasprun.xml.

Note, this class works a bit differently than most of the other
VaspObjects, since the OUTCAR can be very different depending on which
“type of run” performed.

Creating the OUTCAR class with a filename reads “regular parameters” that
are always present.


#### magnetization()
Magnetization on each ion as a tuple of dict, e.g.,
({“d”: 0.0, “p”: 0.003, “s”: 0.002, “tot”: 0.005}, … )
Note that this data is not always present. LORBIT must be set to some
other value than the default.


#### chemical_shielding()
chemical shielding on each ion as a dictionary with core and valence contributions


#### unsym_cs_tensor()
Unsymmetrized chemical shielding tensor matrixes on each ion as a list.
e.g.,
[[[sigma11, sigma12, sigma13],

> > [sigma21, sigma22, sigma23],
> > [sigma31, sigma32, sigma33]],
> > …

> [[sigma11, sigma12, sigma13],

>     [sigma21, sigma22, sigma23],
>     [sigma31, sigma32, sigma33]]]


#### cs_g0_contribution()
G=0 contribution to chemical shielding. 2D rank 3 matrix


#### cs_core_contribution()
Core contribution to chemical shielding. dict. e.g.,
{‘Mg’: -412.8, ‘C’: -200.5, ‘O’: -271.1}


#### efg()
Electric Field Gradient (EFG) tensor on each ion as a tuple of dict, e.g.,
({“cq”: 0.1, “eta”, 0.2, “nuclear_quadrupole_moment”: 0.3},

> {“cq”: 0.7, “eta”, 0.8, “nuclear_quadrupole_moment”: 0.9},
> …)


#### charge()
Charge on each ion as a tuple of dict, e.g.,
({“p”: 0.154, “s”: 0.078, “d”: 0.0, “tot”: 0.232}, …)
Note that this data is not always present. LORBIT must be set to some
other value than the default.


#### is_stopped()
True if OUTCAR is from a stopped run (using STOPCAR, see VASP Manual).


#### run_stats()
Various useful run stats as a dict including “System time (sec)”,
“Total CPU time used (sec)”, “Elapsed time (sec)”,
“Maximum memory used (kb)”, “Average memory used (kb)”,
“User time (sec)”, “cores”


#### elastic_tensor()
Total elastic moduli (Kbar) is given in a 6x6 array matrix.


#### drift()
Total drift for each step in eV/Atom


#### ngf()
Dimensions for the Augementation grid

<!-- attribute: sampling_radii

Size of the sampling radii in VASP for the test charges for
the electrostatic potential at each atom. Total array size is the number
of elements present in the calculation -->
<!-- attribute: electrostatic_potential

Average electrostatic potential at each atomic position in order
of the atoms in POSCAR. -->
..attribute: final_energy_contribs

> Individual contributions to the total final energy as a dictionary.
> Include contirbutions from keys, e.g.:
> {‘DENC’: -505778.5184347, ‘EATOM’: 15561.06492564, ‘EBANDS’: -804.53201231,
> ‘EENTRO’: -0.08932659, ‘EXHF’: 0.0, ‘Ediel_sol’: 0.0,
> ‘PAW double counting’: 664.6726974100002, ‘PSCENC’: 742.48691646,
> ‘TEWEN’: 489742.86847338, ‘XCENC’: -169.64189814}


#### efermi()
Fermi energy


#### filename()
> Filename


#### final_energy()
Final energy after extrapolation of sigma back to 0, i.e. energy(sigma->0).


#### final_energy_wo_entrp()
Final energy before extrapolation of sigma, i.e. energy without entropy.


#### final_fr_energy()
Final “free energy”, i.e. free energy TOTEN


#### has_onsite_density_matrices()
Boolean for if onsite density matrices have been set


#### lcalcpol()
If LCALCPOL has been set


#### lepsilon()
If LEPSILON has been set


#### nelect()
Returns the number of electrons in the calculation


#### spin()
If spin-polarization was enabled via ISPIN


#### total_mag()
Total magnetization (in terms of the number of unpaired electrons)

One can then call a specific reader depending on the type of run being
performed. These are currently: read_igpar(), read_lepsilon() and
read_lcalcpol(), read_core_state_eign(), read_avg_core_pot().

See the documentation of those methods for more documentation.

Authors: Rickard Armiento, Shyue Ping Ong


* **Parameters**

    **filename** (*str*) – OUTCAR filename to parse.



#### as_dict()

* **Returns**

    MSONable dict.



#### read_avg_core_poten()
Read the core potential at each ionic step.


* **Returns**

    A list for each ionic step containing a list of the average core
    potentials for each atom: [[avg core pot]].


### Example

The average core potential of the 2nd atom of the structure at the
last ionic step is: [-1][1]


#### read_chemical_shielding()
Parse the NMR chemical shieldings data. Only the second part “absolute, valence and core”
will be parsed. And only the three right most field (ISO_SHIELDING, SPAN, SKEW) will be retrieved.


* **Returns**

    List of chemical shieldings in the order of atoms from the OUTCAR. Maryland notation is adopted.



#### read_core_state_eigen()
Read the core state eigenenergies at each ionic step.


* **Returns**

    [core state eig]}].
    The core state eigenenergie list for each AO is over all ionic
    step.



* **Return type**

    A list of dict over the atom such as [{“AO”


### Example

The core state eigenenergie of the 2s AO of the 6th atom of the
structure at the last ionic step is [5][“2s”][-1]


#### read_corrections(reverse=True, terminate_on_match=True)
Reads the dipol qudropol corrections into the
Outcar.data[“dipol_quadrupol_correction”].


* **Parameters**


    * **reverse** – Whether to start from end of OUTCAR.


    * **terminate_on_match** – Whether to terminate once match is found.



#### read_cs_core_contribution()
Parse the core contribution of NMR chemical shielding.

Returns:
G0 contribution matrix as list of list.


#### read_cs_g0_contribution()
Parse the  G0 contribution of NMR chemical shielding.


* **Returns**

    G0 contribution matrix as list of list.



#### read_cs_raw_symmetrized_tensors()
Parse the matrix form of NMR tensor before corrected to table.


* **Returns**

    nsymmetrized tensors list in the order of atoms.



#### read_elastic_tensor()
Parse the elastic tensor data.


* **Returns**

    6x6 array corresponding to the elastic tensor from the OUTCAR.



#### read_electrostatic_potential()
Parses the eletrostatic potential for the last ionic step.


#### read_fermi_contact_shift()
Output example:
Fermi contact (isotropic) hyperfine coupling parameter (MHz)
————————————————————-
ion      A_pw      A_1PS     A_1AE     A_1c      A_tot
————————————————————-

> 1      -0.002    -0.002    -0.051     0.000    -0.052
> 2      -0.002    -0.002    -0.051     0.000    -0.052
> 3       0.056     0.056     0.321    -0.048     0.321

> [-0.002, -0.002, -0.051, 0.0, -0.052],
> [0.056, 0.056, 0.321, -0.048, 0.321]] from ‘fch’ data.


#### read_freq_dielectric()
Parses the frequency dependent dielectric function (obtained with
LOPTICS). Frequencies (in eV) are in self.frequencies, and dielectric
tensor function is given as self.dielectric_tensor_function.


#### read_igpar()
Renders accessible:

    er_ev = e<r>_ev (dictionary with Spin.up/Spin.down as keys)
    er_bp = e<r>_bp (dictionary with Spin.up/Spin.down as keys)
    er_ev_tot = spin up + spin down summed
    er_bp_tot = spin up + spin down summed
    p_elc = spin up + spin down summed
    p_ion = spin up + spin down summed.

(See VASP section “LBERRY,  IGPAR,  NPPSTR,  DIPOL tags” for info on
what these are).


#### read_internal_strain_tensor()
Reads the internal strain tensor and populates self.internal_strain_tensor with an array of voigt notation

    tensors for each site.


#### read_lcalcpol()
Reads the lcalpol.

# TODO: Document the actual variables.


#### read_lepsilon()
Reads an LEPSILON run.

# TODO: Document the actual variables.


#### read_lepsilon_ionic()
Reads an LEPSILON run, the ionic component.

# TODO: Document the actual variables.


#### read_neb(reverse=True, terminate_on_match=True)
Reads NEB data. This only works with OUTCARs from both normal
VASP NEB calculations or from the CI NEB method implemented by
Henkelman et al.


* **Parameters**


    * **reverse** (*bool*) – Read files in reverse. Defaults to false. Useful for
    large files, esp OUTCARs, especially when used with
    terminate_on_match. Defaults to True here since we usually
    want only the final value.


    * **terminate_on_match** (*bool*) – Whether to terminate when there is at
    least one match in each key in pattern. Defaults to True here
    since we usually want only the final value.


Renders accessible:

    tangent_force - Final tangent force.
    energy - Final energy.
    These can be accessed under Outcar.data[key]


#### read_nmr_efg()
Parse the NMR Electric Field Gradient interpreted values.


* **Returns**

    Electric Field Gradient tensors as a list of dict in the order of atoms from OUTCAR.
    Each dict key/value pair corresponds to a component of the tensors.



#### read_nmr_efg_tensor()
Parses the NMR Electric Field Gradient Raw Tensors.


* **Returns**

    A list of Electric Field Gradient Tensors in the order of Atoms from OUTCAR



#### read_onsite_density_matrices()
Parse the onsite density matrices, returns list with index corresponding
to atom index in Structure.


#### read_pattern(patterns, reverse=False, terminate_on_match=False, postprocess=<class 'str'>)
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


#### read_piezo_tensor()
Parse the piezo tensor data.


#### read_pseudo_zval()
Create pseudopotential ZVAL dictionary.


#### read_table_pattern(header_pattern, row_pattern, footer_pattern, postprocess=<class 'str'>, attribute_name=None, last_one_only=True, first_one_only=False)
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
    be returned. Default to be True. Incompatible with first_one_only.


    * **first_one_only** (*bool*) – Only the first occurrence of the table will be
    parsed and the parsing procedure will stop. The enclosing list
    will be removed. i.e. Only a single table will be returned.
    Incompatible with last_one_only.



* **Returns**

    List of tables. 1) A table is a list of rows. 2) A row if either a list of
    attribute values in case the capturing group is defined without name in
    row_pattern, or a dict in case that named capturing groups are defined by
    row_pattern.



### _class_ pymatgen.io.vasp.outputs.Procar(filename)
Bases: `object`

Object for reading a PROCAR file.


#### data()
The PROCAR data of the form below. It should VASP uses 1-based indexing,
but all indices are converted to 0-based here.:

```default
{
    spin: nd.array accessed with (k-point index, band index,
                                  ion index, orbital index)
}
```


#### weights()
The weights associated with each k-point as an nd.array of length
nkpoints.

..attribute:: phase_factors

> Phase factors, where present (e.g. LORBIT = 12). A dict of the form:
> {

> > spin: complex nd.array accessed with (k-point index, band index,

> >     ion index, orbital index)

> }

..attribute:: nbands

> Number of bands

..attribute:: nkpoints

> Number of k-points

..attribute:: nions

> Number of ions


* **Parameters**

    **filename** – Name of file containing PROCAR.



#### get_occupation(atom_index, orbital)
Returns the occupation for a particular orbital of a particular atom.


* **Parameters**


    * **atom_num** (*int*) – Index of atom in the PROCAR. It should be noted
    that VASP uses 1-based indexing for atoms, but this is
    converted to 0-based indexing in this parser to be
    consistent with representation of structures in pymatgen.


    * **orbital** (*str*) – An orbital. If it is a single character, e.g., s,
    p, d or f, the sum of all s-type, p-type, d-type or f-type
    orbitals occupations are returned respectively. If it is a
    specific orbital, e.g., px, dxy, etc., only the occupation
    of that orbital is returned.



* **Returns**

    Sum occupation of orbital of atom.



#### get_projection_on_elements(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))
Method returning a dictionary of projections on elements.


* **Parameters**

    **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Input structure.



* **Returns**

    [k index][b index][{Element:values}]]



* **Return type**

    a dictionary in the {Spin.up



### _exception_ pymatgen.io.vasp.outputs.UnconvergedVASPWarning()
Bases: `Warning`

Warning for unconverged vasp run.


### _exception_ pymatgen.io.vasp.outputs.VaspParseError()
Bases: [`ParseError`](pymatgen.io.core.md#pymatgen.io.core.ParseError)

Exception class for VASP parsing.


### _class_ pymatgen.io.vasp.outputs.Vasprun(filename, ionic_step_skip=None, ionic_step_offset=0, parse_dos=True, parse_eigen=True, parse_projected_eigen=False, parse_potcar_file=True, occu_tol=1e-08, separate_spins=False, exception_on_bad_xml=True)
Bases: `MSONable`

Vastly improved cElementTree-based parser for vasprun.xml files. Uses
iterparse to support incremental parsing of large files.
Speedup over Dom is at least 2x for smallish files (~1Mb) to orders of
magnitude for larger files (~10Mb).

**VASP results**


#### ionic_steps()
All ionic steps in the run as a list of
{“structure”: structure at end of run,
“electronic_steps”: {All electronic step data in vasprun file},
“stresses”: stress matrix}


#### tdos()
Total dos calculated at the end of run.


#### idos()
Integrated dos calculated at the end of run.


#### pdos()
List of list of PDos objects. Access as pdos[atomindex][orbitalindex]


#### efermi()
Fermi energy


#### eigenvalues()
Available only if parse_eigen=True. Final eigenvalues as a dict of
{(spin, kpoint index):[[eigenvalue, occu]]}.
This representation is based on actual ordering in VASP and is meant as
an intermediate representation to be converted into proper objects. The
kpoint index is 0-based (unlike the 1-based indexing in VASP).


#### projected_eigenvalues()
Final projected eigenvalues as a dict of {spin: nd-array}. To access
a particular value, you need to do
Vasprun.projected_eigenvalues[spin][kpoint index][band index][atom index][orbital_index]
This representation is based on actual ordering in VASP and is meant as
an intermediate representation to be converted into proper objects. The
kpoint, band and atom indices are 0-based (unlike the 1-based indexing
in VASP).


#### projected_magnetisation()
Final projected magnetisation as a numpy array with the shape (nkpoints, nbands,
natoms, norbitals, 3). Where the last axis is the contribution in the 3
Cartesian directions. This attribute is only set if spin-orbit coupling
(LSORBIT = True) or non-collinear magnetism (LNONCOLLINEAR = True) is turned
on in the INCAR.


#### other_dielectric()
Dictionary, with the tag comment as key, containing other variants of
the real and imaginary part of the dielectric constant (e.g., computed
by RPA) in function of the energy (frequency). Optical properties (e.g.
absorption coefficient) can be obtained through this.
The data is given as a tuple of 3 values containing each of them
the energy, the real part tensor, and the imaginary part tensor
([energies],[[real_partxx,real_partyy,real_partzz,real_partxy,
real_partyz,real_partxz]],[[imag_partxx,imag_partyy,imag_partzz,
imag_partxy, imag_partyz, imag_partxz]])


#### nionic_steps()
The total number of ionic steps. This number is always equal
to the total number of steps in the actual run even if
ionic_step_skip is used.


#### force_constants()
Force constants computed in phonon DFPT run(IBRION = 8).
The data is a 4D numpy array of shape (natoms, natoms, 3, 3).


#### normalmode_eigenvals()
Normal mode frequencies.
1D numpy array of size 3\*natoms.


#### normalmode_eigenvecs()
Normal mode eigen vectors.
3D numpy array of shape (3\*natoms, natoms, 3).


#### md_data()
Available only for ML MD runs, i.e., INCAR with ML_LMLFF = .TRUE.
md_data is a list of dict with the following format:

[

    {

        ‘energy’: {

            ‘e_0_energy’: -525.07195568,
            ‘e_fr_energy’: -525.07195568,
            ‘e_wo_entrp’: -525.07195568,
            ‘kinetic’: 3.17809233,
            ‘lattice kinetic’: 0.0,
            ‘nosekinetic’: 1.323e-05,
            ‘nosepot’: 0.0,
            ‘total’: -521.89385012
            },

        ‘forces’: [[0.17677989, 0.48309874, 1.85806696], …],
        ‘structure’: Structure object

    }

]

**VASP inputs**


#### incar()
Incar object for parameters specified in INCAR file.


#### parameters()
Incar object with parameters that vasp actually used, including all
defaults.


#### kpoints()
Kpoints object for KPOINTS specified in run.


#### actual_kpoints()
List of actual kpoints, e.g.,
[[0.25, 0.125, 0.08333333], [-0.25, 0.125, 0.08333333],
[0.25, 0.375, 0.08333333], ….]


#### actual_kpoints_weights()
List of kpoint weights, E.g.,
[0.04166667, 0.04166667, 0.04166667, 0.04166667, 0.04166667, ….]


#### atomic_symbols()
List of atomic symbols, e.g., [“Li”, “Fe”, “Fe”, “P”, “P”, “P”]


#### potcar_symbols()
List of POTCAR symbols. e.g.,
[“PAW_PBE Li 17Jan2003”, “PAW_PBE Fe 06Sep2000”, ..]

Author: Shyue Ping Ong


* **Parameters**


    * **filename** (*str*) – Filename to parse


    * **ionic_step_skip** (*int*) – If ionic_step_skip is a number > 1,
    only every ionic_step_skip ionic steps will be read for
    structure and energies. This is very useful if you are parsing
    very large vasprun.xml files and you are not interested in every
    single ionic step. Note that the final energies may not be the
    actual final energy in the vasprun.


    * **ionic_step_offset** (*int*) – Used together with ionic_step_skip. If set,
    the first ionic step read will be offset by the amount of
    ionic_step_offset. For example, if you want to start reading
    every 10th structure but only from the 3rd structure onwards,
    set ionic_step_skip to 10 and ionic_step_offset to 3. Main use
    case is when doing statistical structure analysis with
    extremely long time scale multiple VASP calculations of
    varying numbers of steps.


    * **parse_dos** (*bool*) – Whether to parse the dos. Defaults to True. Set
    to False to shave off significant time from the parsing if you
    are not interested in getting those data.


    * **parse_eigen** (*bool*) – Whether to parse the eigenvalues. Defaults to
    True. Set to False to shave off significant time from the
    parsing if you are not interested in getting those data.


    * **parse_projected_eigen** (*bool*) – Whether to parse the projected
    eigenvalues and magnetisation. Defaults to False. Set to True to obtain
    projected eigenvalues and magnetisation. **Note that this can take an
    extreme amount of time and memory.** So use this wisely.


    * **parse_potcar_file** (*bool/str*) – Whether to parse the potcar file to read
    the potcar hashes for the potcar_spec attribute. Defaults to True,
    where no hashes will be determined and the potcar_spec dictionaries
    will read {“symbol”: ElSymbol, “hash”: None}. By Default, looks in
    the same directory as the vasprun.xml, with same extensions as

    > Vasprun.xml. If a string is provided, looks at that filepath.



    * **occu_tol** (*float*) – Sets the minimum tol for the determination of the
    vbm and cbm. Usually the default of 1e-8 works well enough,
    but there may be pathological cases.


    * **separate_spins** (*bool*) – Whether the band gap, CBM, and VBM should be
    reported for each individual spin channel. Defaults to False,
    which computes the eigenvalue band properties independent of
    the spin orientation. If True, the calculation must be spin-polarized.


    * **exception_on_bad_xml** (*bool*) – Whether to throw a ParseException if a
    malformed XML is detected. Default to True, which ensures only
    proper vasprun.xml are parsed. You can set to False if you want
    partial results (e.g., if you are monitoring a calculation during a
    run), but use the results with care. A warning is issued.



#### as_dict()
JSON-serializable dict representation.


#### calculate_efermi(tol: float = 0.001)
Calculate the Fermi level using a robust algorithm.

Sometimes VASP can put the Fermi level just inside of a band due to issues in
the way band occupancies are handled. This algorithm tries to detect and correct
for this bug.

Slightly more details are provided here: [https://www.vasp.at/forum/viewtopic.php?f=4&t=17981](https://www.vasp.at/forum/viewtopic.php?f=4&t=17981)


#### _property_ complete_dos()
A complete dos object which incorporates the total dos and all
projected dos.


#### _property_ complete_dos_normalized(_: [CompleteDos](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.CompleteDos_ )
A CompleteDos object which incorporates the total DOS and all
projected DOS. Normalized by the volume of the unit cell with
units of states/eV/unit cell volume.


#### _property_ converged()
Returns:
True if a relaxation run is converged both ionically and
electronically.


#### _property_ converged_electronic()
Returns:
True if electronic step convergence has been reached in the final
ionic step.


#### _property_ converged_ionic()
Returns:
True if ionic step convergence has been reached, i.e. that vasp
exited before reaching the max ionic steps for a relaxation run.


#### _property_ dielectric()
Returns:
The real and imaginary part of the dielectric constant (e.g., computed
by RPA) in function of the energy (frequency). Optical properties (e.g.
absorption coefficient) can be obtained through this.
The data is given as a tuple of 3 values containing each of them
the energy, the real part tensor, and the imaginary part tensor
([energies],[[real_partxx,real_partyy,real_partzz,real_partxy,
real_partyz,real_partxz]],[[imag_partxx,imag_partyy,imag_partzz,
imag_partxy, imag_partyz, imag_partxz]]).


#### _property_ eigenvalue_band_properties()
Band properties from the eigenvalues as a tuple,
(band gap, cbm, vbm, is_band_gap_direct). In the case of separate_spins=True,
the band gap, cbm, vbm, and is_band_gap_direct are each lists of length 2,
with index 0 representing the spin-up channel and index 1 representing
the spin-down channel.


#### _property_ epsilon_ionic()
Property only available for DFPT calculations and when IBRION=5, 6, 7 or 8.


* **Returns**

    The ionic part of the static dielectric constant. Present when it’s a
    DFPT run (LEPSILON=TRUE) and IBRION=5, 6, 7 or 8



#### _property_ epsilon_static()
Property only available for DFPT calculations.


* **Returns**

    The static part of the dielectric constant. Present when it’s a DFPT run
    (LEPSILON=TRUE)



#### _property_ epsilon_static_wolfe()
Property only available for DFPT calculations.


* **Returns**

    The static part of the dielectric constant without any local field
    effects. Present when it’s a DFPT run (LEPSILON=TRUE)



#### _property_ final_energy()

#### get_band_structure(kpoints_filename: str | None = None, efermi: float | Literal['smart'] | None = None, line_mode: bool = False, force_hybrid_mode: bool = False)
Get the band structure as a BandStructure object.


* **Parameters**


    * **kpoints_filename** – Full path of the KPOINTS file from which
    the band structure is generated.
    If none is provided, the code will try to intelligently
    determine the appropriate KPOINTS file by substituting the
    filename of the vasprun.xml with KPOINTS.
    The latter is the default behavior.


    * **efermi** – The Fermi energy associated with the bandstructure, in eV. By
    default (None), uses the value reported by VASP in vasprun.xml. To
    manually set the Fermi energy, pass a float. Pass ‘smart’ to use the
    calculate_efermi() method, which calculates the Fermi level by first
    checking whether it lies within a small tolerance (by default 0.001 eV)
    of a band edge) If it does, the Fermi level is placed in the center of
    the bandgap. Otherwise, the value is identical to the value reported by
    VASP.


    * **line_mode** – Force the band structure to be considered as
    a run along symmetry lines. (Default: False)


    * **force_hybrid_mode** – Makes it possible to read in self-consistent band
    structure calculations for every type of functional. (Default: False)



* **Returns**

    a BandStructure object (or more specifically a
    BandStructureSymmLine object if the run is detected to be a run
    along symmetry lines)

    Two types of runs along symmetry lines are accepted: non-sc with
    Line-Mode in the KPOINT file or hybrid, self-consistent with a
    uniform grid+a few kpoints along symmetry lines (explicit KPOINTS
    file) (it’s not possible to run a non-sc band structure with hybrid
    functionals). The explicit KPOINTS file needs to have data on the
    kpoint label as commentary.




#### get_computed_entry(inc_structure=True, parameters=None, data=None, entry_id: str | None = None)
Returns a ComputedEntry or ComputedStructureEntry from the Vasprun.


* **Parameters**


    * **inc_structure** (*bool*) – Set to True if you want
    ComputedStructureEntries to be returned instead of
    ComputedEntries.


    * **parameters** (*list*) – Input parameters to include. It has to be one of
    the properties supported by the Vasprun object. If
    parameters is None, a default set of parameters that are
    necessary for typical post-processing will be set.


    * **data** (*list*) – Output data to include. Has to be one of the properties
    supported by the Vasprun object.


    * **entry_id** (*str*) – Specify an entry id for the ComputedEntry. Defaults to
    “vasprun-{current datetime}”



* **Returns**

    ComputedStructureEntry/ComputedEntry



#### get_potcars(path)

* **Parameters**

    **path** – Path to search for POTCARs



* **Returns**

    Potcar from path.



#### get_trajectory()
This method returns a Trajectory object, which is an alternative
representation of self.structures into a single object. Forces are
added to the Trajectory as site properties.

Returns: a Trajectory


#### _property_ hubbards()
Hubbard U values used if a vasprun is a GGA+U run. {} otherwise.


#### _property_ is_hubbard(_: boo_ )
True if run is a DFT+U run.


#### _property_ is_spin(_: boo_ )
True if run is spin-polarized.


#### _property_ optical_absorption_coeff()
Calculate the optical absorption coefficient
from the dielectric constants. Note that this method is only
implemented for optical properties calculated with GGA and BSE.


* **Returns**

    optical absorption coefficient in list



#### _property_ run_type()
Returns the run type. Currently detects GGA, metaGGA, HF, HSE, B3LYP,
and hybrid functionals based on relevant INCAR tags. LDA is assigned if
PAW POTCARs are used and no other functional is detected.

Hubbard U terms and vdW corrections are detected automatically as well.


#### _property_ structures()
Returns:
List of Structure objects for the structure at each ionic step.


#### update_charge_from_potcar(path)
Sets the charge of a structure based on the POTCARs found.


* **Parameters**

    **path** – Path to search for POTCARs



#### update_potcar_spec(path)

* **Parameters**

    **path** – Path to search for POTCARs



* **Returns**

    Potcar spec from path.



### _class_ pymatgen.io.vasp.outputs.VolumetricData(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), data, distance_matrix=None, data_aug=None)
Bases: [`VolumetricData`](pymatgen.io.common.md#pymatgen.io.common.VolumetricData)

Container for volumetric data that allows
for reading/writing with Poscar-type data.

Typically, this constructor is not used directly and the static
from_file constructor is used. This constructor is designed to allow
summation and other operations between VolumetricData objects.


* **Parameters**


    * **structure** – Structure associated with the volumetric data


    * **data** – Actual volumetric data. If the data is provided as in list format,
    it will be converted into an np.array automatically


    * **data_aug** – Any extra information associated with volumetric data
    (typically augmentation charges)


    * **distance_matrix** – A pre-computed distance matrix if available.
    Useful so pass distance_matrices between sums,
    short-circuiting an otherwise expensive operation.



#### _static_ parse_file(filename)
Convenience method to parse a generic volumetric data file in the vasp
like format. Used by subclasses for parsing file.


* **Parameters**

    **filename** (*str*) – Path of file to parse



* **Returns**

    (poscar, data)



#### write_file(file_name, vasp4_compatible=False)
Write the VolumetricData object to a vasp compatible file.


* **Parameters**


    * **file_name** (*str*) – Path to a file


    * **vasp4_compatible** (*bool*) – True if the format is vasp4 compatible



### _class_ pymatgen.io.vasp.outputs.WSWQ(nspin: int, nkpoints: int, nbands: int, me_real: ndarray, me_imag: ndarray)
Bases: `MSONable`

Class for reading a WSWQ file.
The WSWQ file is used to calculation the wave function overlaps between

>
> * W: Wavefunctions in the currenct directory’s WAVECAR file


> * WQ: Wavefunctions stored in a filed named the WAVECAR.qqq.

The overlap is computed using the overlap operator S
which make the PAW wavefunctions orthogonormal:

> <W_k,m| S | W_k,n> = delta_{mn}

The WSWQ file contains matrix elements of the overlap operator S evaluated
between the planewave wavefunctions W and WQ:

> COVL_k,mn = < W_s,k,m | S | WQ_s,k,n >

The indices of WSWQ.data are:

    [spin][kpoint][band_i][band_j]


#### nspin()
Number of spin channels


* **Type**

    int



#### nkpoints()
Number of k-points


* **Type**

    int



#### nbands()
Number of bands


* **Type**

    int



#### me_real()
Real part of the overlap matrix elements


* **Type**

    numpy.ndarray



#### me_imag()
Imaginary part of the overlap matrix elements


* **Type**

    numpy.ndarray



#### _property_ data()
Complex overlap matrix.


#### _classmethod_ from_file(filename)
Constructs a WSWQ object from a file.


* **Parameters**

    **filename** – Name of WSWQ file.



* **Returns**

    WSWQ object.



#### me_imag(_: ndarra_ )

#### me_real(_: ndarra_ )

#### nbands(_: in_ )

#### nkpoints(_: in_ )

#### nspin(_: in_ )

### _class_ pymatgen.io.vasp.outputs.Wavecar(filename='WAVECAR', verbose=False, precision='normal', vasp_type=None)
Bases: `object`

This is a class that contains the (pseudo-) wavefunctions from VASP.

Coefficients are read from the given WAVECAR file and the corresponding
G-vectors are generated using the algorithm developed in WaveTrans (see
acknowledgments below). To understand how the wavefunctions are evaluated,
please see the evaluate_wavefunc docstring.

It should be noted that the pseudopotential augmentation is not included in
the WAVECAR file. As a result, some caution should be exercised when
deriving value from this information.

The usefulness of this class is to allow the user to do projections or band
unfolding style manipulations of the wavefunction. An example of this can
be seen in the work of Shen et al. 2017
([https://doi.org/10.1103/PhysRevMaterials.1.065001](https://doi.org/10.1103/PhysRevMaterials.1.065001)).


#### filename()
String of the input file (usually WAVECAR)


#### vasp_type()
String that determines VASP type the WAVECAR was generated with (either
‘std’, ‘gam’, or ‘ncl’)


#### nk()
Number of k-points from the WAVECAR


#### nb()
Number of bands per k-point


#### encut()
Energy cutoff (used to define G_{cut})


#### efermi()
Fermi energy


#### a()
Primitive lattice vectors of the cell (e.g. a_1 = self.a[0, :])


#### b()
Reciprocal lattice vectors of the cell (e.g. b_1 = self.b[0, :])


#### vol()
The volume of the unit cell in real space


#### kpoints()
The list of k-points read from the WAVECAR file


#### band_energy()
The list of band eigenenergies (and corresponding occupancies) for
each kpoint, where the first index corresponds to the index of the
k-point (e.g. self.band_energy[kp])


#### Gpoints()
The list of generated G-points for each k-point (a double list), which
are used with the coefficients for each k-point and band to recreate
the wavefunction (e.g. self.Gpoints[kp] is the list of G-points for
k-point kp). The G-points depend on the k-point and reciprocal lattice
and therefore are identical for each band at the same k-point. Each
G-point is represented by integer multipliers (e.g. assuming
Gpoints[kp][n] == [n_1, n_2, n_3], then
G_n = n_1\*b_1 + n_2\*b_2 + n_3\*b_3)


#### coeffs()
The list of coefficients for each k-point and band for reconstructing
the wavefunction. For non-spin-polarized, the first index corresponds
to the kpoint and the second corresponds to the band (e.g.
self.coeffs[kp][b] corresponds to k-point kp and band b). For
spin-polarized calculations, the first index is for the spin.
If the calculation was non-collinear, then self.coeffs[kp][b] will have
two columns (one for each component of the spinor).

Acknowledgments:

    This code is based upon the Fortran program, WaveTrans, written by
    R. M. Feenstra and M. Widom from the Dept. of Physics at Carnegie
    Mellon University. To see the original work, please visit:
    [https://www.andrew.cmu.edu/user/feenstra/wavetrans/](https://www.andrew.cmu.edu/user/feenstra/wavetrans/)

Author: Mark Turiansky

Information is extracted from the given WAVECAR.


* **Parameters**


    * **filename** (*str*) – input file (default: WAVECAR)


    * **verbose** (*bool*) – determines whether processing information is shown


    * **precision** (*str*) – determines how fine the fft mesh is (normal or
    accurate), only the first letter matters


    * **vasp_type** (*str*) – determines the VASP type that is used, allowed
    values are [‘std’, ‘gam’, ‘ncl’]
    (only first letter is required)



#### evaluate_wavefunc(kpoint: int, band: int, r: ndarray, spin: int = 0, spinor: int = 0)
Evaluates the wavefunction for a given position, r.

The wavefunction is given by the k-point and band. It is evaluated
at the given position by summing over the components. Formally,

psi_n^k (r) = sum_{i=1}^N c_i^{n,k} exp (i (k + G_i^{n,k}) cdot r)

where psi_n^k is the wavefunction for the nth band at k-point k, N is
the number of plane waves, c_i^{n,k} is the ith coefficient that
corresponds to the nth band and k-point k, and G_i^{n,k} is the ith
G-point corresponding to k-point k.

NOTE: This function is very slow; a discrete fourier transform is the
preferred method of evaluation (see Wavecar.fft_mesh).


* **Parameters**


    * **kpoint** (*int*) – the index of the kpoint where the wavefunction
    will be evaluated


    * **band** (*int*) – the index of the band where the wavefunction will be
    evaluated


    * **r** (*np.array*) – the position where the wavefunction will be evaluated


    * **spin** (*int*) – spin index for the desired wavefunction (only for
    ISPIN = 2, default = 0)


    * **spinor** (*int*) – component of the spinor that is evaluated (only used
    if vasp_type == ‘ncl’)



* **Returns**

    a complex value corresponding to the evaluation of the wavefunction



#### fft_mesh(kpoint: int, band: int, spin: int = 0, spinor: int = 0, shift: bool = True)
Places the coefficients of a wavefunction onto an fft mesh.

Once the mesh has been obtained, a discrete fourier transform can be
used to obtain real-space evaluation of the wavefunction. The output
of this function can be passed directly to numpy’s fft function. For
.. rubric:: Example

mesh = Wavecar(‘WAVECAR’).fft_mesh(kpoint, band)
evals = np.fft.ifftn(mesh)


* **Parameters**


    * **kpoint** (*int*) – the index of the kpoint where the wavefunction
    will be evaluated


    * **band** (*int*) – the index of the band where the wavefunction will be
    evaluated


    * **spin** (*int*) – the spin of the wavefunction for the desired
    wavefunction (only for ISPIN = 2, default = 0)


    * **spinor** (*int*) – component of the spinor that is evaluated (only used
    if vasp_type == ‘ncl’)


    * **shift** (*bool*) – determines if the zero frequency coefficient is
    placed at index (0, 0, 0) or centered



* **Returns**

    a numpy ndarray representing the 3D mesh of coefficients



#### get_parchg(poscar: [Poscar](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Poscar), kpoint: int, band: int, spin: int | None = None, spinor: int | None = None, phase: bool = False, scale: int = 2)
Generates a Chgcar object, which is the charge density of the specified
wavefunction.

This function generates a Chgcar object with the charge density of the
wavefunction specified by band and kpoint (and spin, if the WAVECAR
corresponds to a spin-polarized calculation). The phase tag is a
feature that is not present in VASP. For a real wavefunction, the phase
tag being turned on means that the charge density is multiplied by the
sign of the wavefunction at that point in space. A warning is generated
if the phase tag is on and the chosen kpoint is not Gamma.

Note: Augmentation from the PAWs is NOT included in this function. The
maximal charge density will differ from the PARCHG from VASP, but the
qualitative shape of the charge density will match.


* **Parameters**


    * **poscar** ([*pymatgen.io.vasp.inputs.Poscar*](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Poscar)) – Poscar object that has the
    structure associated with the WAVECAR file


    * **kpoint** (*int*) – the index of the kpoint for the wavefunction


    * **band** (*int*) – the index of the band for the wavefunction


    * **spin** (*int*) – optional argument to specify the spin. If the Wavecar
    has ISPIN = 2, spin is None generates a Chgcar with total spin
    and magnetization, and spin == {0, 1} specifies just the spin
    up or down component.


    * **spinor** (*int*) – optional argument to specify the spinor component
    for noncollinear data wavefunctions (allowed values of None,
    0, or 1)


    * **phase** (*bool*) – flag to determine if the charge density is multiplied
    by the sign of the wavefunction. Only valid for real
    wavefunctions.


    * **scale** (*int*) – scaling for the FFT grid. The default value of 2 is at
    least as fine as the VASP default.



* **Returns**

    a pymatgen.io.vasp.outputs.Chgcar object



#### write_unks(directory: str)
Write the UNK files to the given directory.

Writes the cell-periodic part of the bloch wavefunctions from the
WAVECAR file to each of the UNK files. There will be one UNK file for
each of the kpoints in the WAVECAR file.

**NOTE**: wannier90 expects the full kpoint grid instead of the symmetry-
reduced one that VASP stores the wavefunctions on. You should run
a nscf calculation with ISYM=0 to obtain the correct grid.


* **Parameters**

    **directory** (*str*) – directory where the UNK files are written



### _class_ pymatgen.io.vasp.outputs.Waveder(cder_real: ndarray, cder_imag: ndarray)
Bases: `MSONable`

Representation of the WAVEDER file.

The LOPTICS tag produces a WAVEDER file which contains the derivative of the orbitals with respect to k.
Since the data is complex, we need to split it into the real and imaginary parts for JSON serialization.

**NOTE**: The way that VASP writes the WAVEDER and WAVEDERF has slightly different logic when indexing the bands.
This results in the formatted WAVDERF only indexing between filled bands. (i.e. all the matrix elements
are between the states i=1:8 and j=1:8 in a two atom Si calculation, which is likely a VASP bug).
As such, it is recommended to used the hidden `LVEL=.True.` flag in VASP which will force indexing over
all bands.

The order of the indices of the data are:
[

> band index1,
> band index2,
> kpoint index,
> spin index,
> cartesian direction,

]


#### cder_real()
Real part of the derivative of the orbitals with respect to k.


* **Type**

    numpy.ndarray



#### cder_imag()
Imaginary part of the derivative of the orbitals with respect to k.


* **Type**

    numpy.ndarray


Author: Miguel Dias Costa, Kamal Choudhary, Jimmy-Xuan Shen


#### _property_ cder()
Return the complex derivative of the orbitals with respect to k.


#### cder_imag(_: ndarra_ )

#### cder_real(_: ndarra_ )

#### _classmethod_ from_binary(filename, data_type='complex64')
Read the WAVEDER file and returns a Waveder object.


* **Parameters**


    * **filename** – Name of file containing WAVEDER.


    * **data_type** – Data type of the WAVEDER file. Default is complex64.
    If the file was generated with the “gamma” version of VASP,
    the data type can be either “float64” or “float32”.



* **Returns**

    Waveder object.



#### _classmethod_ from_formatted(filename)
Reads the WAVEDERF file and returns a Waveder object.

Note: This file is only produced when LOPTICS is true AND vasp has been
recompiled after uncommenting the line that calls
WRT_CDER_BETWEEN_STATES_FORMATTED in linear_optics.F
It is recommended to use from_binary instead since the binary file is
much smaller and contains the same information.


* **Parameters**

    **filename** (*str*) – The name of the WAVEDER file.



* **Returns**

    A Waveder object.



#### get_orbital_derivative_between_states(band_i, band_j, kpoint, spin, cart_dir)
Method returning a value
between bands band_i and band_j for k-point index, spin-channel and Cartesian direction.


* **Parameters**


    * **band_i** (*int*) – Index of band i


    * **band_j** (*int*) – Index of band j


    * **kpoint** (*int*) – Index of k-point


    * **spin** (*int*) – Index of spin-channel (0 or 1)


    * **cart_dir** (*int*) – Index of Cartesian direction (0,1,2)



* **Returns**

    a float value



#### _property_ nbands()
Returns the number of bands.


#### _property_ nkpoints()
Returns the number of k-points.


#### _property_ nspin()
Returns the number of spin channels.


### _class_ pymatgen.io.vasp.outputs.Xdatcar(filename, ionicstep_start=1, ionicstep_end=None, comment=None)
Bases: `object`

Class representing an XDATCAR file. Only tested with VASP 5.x files.


#### structures()
List of structures parsed from XDATCAR.


#### comment()
Optional comment string.

Authors: Ram Balachandran

Init a Xdatcar.


* **Parameters**


    * **filename** (*str*) – Filename of input XDATCAR file.


    * **ionicstep_start** (*int*) – Starting number of ionic step.


    * **ionicstep_end** (*int*) – Ending number of ionic step.


    * **comment** (*str*) – Optional comment attached to this set of structures.



#### concatenate(filename, ionicstep_start=1, ionicstep_end=None)
Concatenate structures in file to Xdatcar.


* **Parameters**


    * **filename** (*str*) – Filename of XDATCAR file to be concatenated.


    * **ionicstep_start** (*int*) – Starting number of ionic step.


    * **ionicstep_end** (*int*) – Ending number of ionic step.


TODO(rambalachandran):

    Requires a check to ensure if the new concatenating file has the
    same lattice structure and atoms as the Xdatcar class.


#### get_string(ionicstep_start=1, ionicstep_end=None, significant_figures=8)
Write  Xdatcar class to a string.


* **Parameters**


    * **ionicstep_start** (*int*) – Starting number of ionic step.


    * **ionicstep_end** (*int*) – Ending number of ionic step.


    * **significant_figures** (*int*) – Number of significant figures.



#### _property_ natoms()
Sequence of number of sites of each type associated with the Poscar.
Similar to 7th line in vasp 5+ Xdatcar.


#### _property_ site_symbols()
Sequence of symbols associated with the Xdatcar. Similar to 6th line in
vasp 5+ Xdatcar.


#### write_file(filename, \*\*kwargs)
Write Xdatcar class into a file.


* **Parameters**


    * **filename** (*str*) – Filename of output XDATCAR file.


    * **\*\*kwargs** – Supported kwargs are the same as those for the
    Xdatcar.get_string method and are passed through directly.



### pymatgen.io.vasp.outputs.get_adjusted_fermi_level(efermi, cbm, band_structure)
When running a band structure computations the Fermi level needs to be
take from the static run that gave the charge density used for the non-self
consistent band structure run. Sometimes this Fermi level is however a
little too low because of the mismatch between the uniform grid used in
the static run and the band structure k-points (e.g., the VBM is on Gamma
and the Gamma point is not in the uniform mesh). Here we use a procedure
consisting in looking for energy levels higher than the static Fermi level
(but lower than the LUMO) if any of these levels make the band structure
appears insulating and not metallic anymore, we keep this adjusted fermi
level. This procedure has shown to detect correctly most insulators.


* **Parameters**


    * **efermi** (*float*) – The Fermi energy of the static run.


    * **cbm** (*float*) – The conduction band minimum of the static run.


    * **band_structure** ([*BandStructureSymmLine*](pymatgen.electronic_structure.bandstructure.md#pymatgen.electronic_structure.bandstructure.BandStructureSymmLine)) – A band structure object.



* **Returns**

    A new adjusted Fermi level.



* **Return type**

    float



### pymatgen.io.vasp.outputs.get_band_structure_from_vasp_multiple_branches(dir_name, efermi=None, projections=False)
This method is used to get band structure info from a VASP directory. It
takes into account that the run can be divided in several branches named
“branch_x”. If the run has not been divided in branches the method will
turn to parsing vasprun.xml directly.

The method returns None is there’s a parsing error


* **Parameters**


    * **dir_name** – Directory containing all bandstructure runs.


    * **efermi** – Efermi for bandstructure.


    * **projections** – True if you want to get the data on site projections if
    any. Note that this is sometimes very large



* **Returns**

    A BandStructure Object