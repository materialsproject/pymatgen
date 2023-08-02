---
layout: default
title: pymatgen.io.vasp.sets.md
nav_exclude: true
---

# pymatgen.io.vasp.sets module

This module defines the VaspInputSet abstract base class and a concrete
implementation for the parameters developed and tested by the core team
of pymatgen, including the Materials Virtual Lab, Materials Project and the MIT
high throughput project. The basic concept behind an input set is to specify
a scheme to generate a consistent set of VASP inputs from a structure
without further user intervention. This ensures comparability across
runs.

Read the following carefully before implementing new input sets:


1. 99% of what needs to be done can be done by specifying user_incar_settings
to override some of the defaults of various input sets. Unless there is an
extremely good reason to add a new set, DO NOT add one. E.g., if you want
to turn the Hubbard U off, just set “LDAU”: False as a user_incar_setting.


2. All derivative input sets should inherit from one of the usual MPRelaxSet or
MITRelaxSet, and proper superclass delegation should be used where possible.
In particular, you are not supposed to implement your own as_dict or
from_dict for derivative sets unless you know what you are doing.
Improper overriding the as_dict and from_dict protocols is the major
cause of implementation headaches. If you need an example, look at how the
MPStaticSet or MPNonSCFSets are constructed.

The above are recommendations. The following are UNBREAKABLE rules:


1. All input sets must take in a structure or list of structures as the first
argument.


2. user_incar_settings, user_kpoints_settings and user_<whatever>_settings are
ABSOLUTE. Any new sets you implement must obey this. If a user wants to
override your settings, you assume he knows what he is doing. Do not
magically override user supplied settings. You can issue a warning if you
think the user is wrong.


3. All input sets must save all supplied args and kwargs as instance variables.
E.g., self.my_arg = my_arg and self.kwargs = kwargs in the __init__. This
ensures the as_dict and from_dict work correctly.


### _exception_ pymatgen.io.vasp.sets.BadInputSetWarning()
Bases: `UserWarning`

Warning class for bad but legal inputs.


### _class_ pymatgen.io.vasp.sets.DictSet(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), config_dict: dict[str, Any], files_to_transfer=None, user_incar_settings=None, user_kpoints_settings=None, user_potcar_settings=None, constrain_total_magmom: bool = False, sort_structure: bool = True, user_potcar_functional: UserPotcarFunctional | None = None, force_gamma: bool = False, reduce_structure=None, vdw=None, use_structure_charge: bool = False, standardize: bool = False, sym_prec=0.1, international_monoclinic: bool = True, validate_magmom: bool = True)
Bases: `VaspInputSet`

Concrete implementation of VaspInputSet that is initialized from a dict
settings. This allows arbitrary settings to be input. In general,
this is rarely used directly unless there is a source of settings in yaml
format (e.g., from a REST interface). It is typically used by other
VaspInputSets for initialization.

Special consideration should be paid to the way the MAGMOM initialization
for the INCAR is done. The initialization differs depending on the type of
structure and the configuration settings. The order in which the magmom is
determined is as follows:


1. If the site itself has a magmom setting (i.e. site.properties[“magmom”] = float),

    that is used. This can be set with structure.add_site_property().


2. If the species of the site has a spin setting, that is used. This can be set

    with structure.add_spin_by_element().


3. If the species itself has a particular setting in the config file, that
is used, e.g., Mn3+ may have a different magmom than Mn4+.


4. Lastly, the element symbol itself is checked in the config file. If
there are no settings, a default value of 0.6 is used.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – The Structure to create inputs for.


    * **config_dict** (*dict*) – The config dictionary to use.


    * **files_to_transfer** (*dict*) – A dictionary of {filename: filepath}. This
    allows the transfer of files from a previous calculation.


    * **user_incar_settings** (*dict*) – User INCAR settings. This allows a user
    to override INCAR settings, e.g., setting a different MAGMOM for
    various elements or species. Note that in the new scheme,
    ediff_per_atom and hubbard_u are no longer args. Instead, the
    config_dict supports EDIFF_PER_ATOM and EDIFF keys. The former
    scales with # of atoms, the latter does not. If both are
    present, EDIFF is preferred. To force such settings, just supply
    user_incar_settings={“EDIFF”: 1e-5, “LDAU”: False} for example.
    The keys ‘LDAUU’, ‘LDAUJ’, ‘LDAUL’ are special cases since
    pymatgen defines different values depending on what anions are
    present in the structure, so these keys can be defined in one
    of two ways, e.g. either {“LDAUU”:{“O”:{“Fe”:5}}} to set LDAUU
    for Fe to 5 in an oxide, or {“LDAUU”:{“Fe”:5}} to set LDAUU to
    5 regardless of the input structure.

    If a None value is given, that key is unset. For example,
    {“ENCUT”: None} will remove ENCUT from the incar settings.



    * **user_kpoints_settings** (*dict** or *[*Kpoints*](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Kpoints)) – Allow user to override kpoints
    setting by supplying a dict E.g., {“reciprocal_density”: 1000}.
    User can also supply Kpoints object. Default is None.


    * **(****dict** (*user_potcar_settings*) – Allow user to override POTCARs. E.g.,
    {“Gd”: “Gd_3”}. This is generally not recommended. Default is None.


    * **constrain_total_magmom** (*bool*) – Whether to constrain the total magmom
    (NUPDOWN in INCAR) to be the sum of the expected MAGMOM for all
    species. Defaults to False.


    * **sort_structure** (*bool*) – Whether to sort the structure (using the
    default sort order of electronegativity) before generating input
    files. Defaults to True, the behavior you would want most of the
    time. This ensures that similar atomic species are grouped
    together.


    * **user_potcar_functional** (*str*) – Functional to use. Default (None) is to use
    the functional in the config dictionary. Valid values:
    “PBE”, “PBE_52”, “PBE_54”, “LDA”, “LDA_52”, “LDA_54”, “PW91”,
    “LDA_US”, “PW91_US”.


    * **force_gamma** (*bool*) – Force gamma centered kpoint generation. Default
    (False) is to use the Automatic Density kpoint scheme, which
    will use the Gamma centered generation scheme for hexagonal
    cells, and Monkhorst-Pack otherwise.


    * **reduce_structure** (*None/str*) – Before generating the input files,
    generate the reduced structure. Default (None), does not
    alter the structure. Valid values: None, “niggli”, “LLL”.


    * **vdw** – Adds default parameters for van-der-Waals functionals supported
    by VASP to INCAR. Supported functionals are: DFT-D2, undamped
    DFT-D3, DFT-D3 with Becke-Jonson damping, Tkatchenko-Scheffler,
    Tkatchenko-Scheffler with iterative Hirshfeld partitioning,
    [MBD@rSC](mailto:MBD@rSC), dDsC, Dion’s vdW-DF, DF2, optPBE, optB88, optB86b and
    rVV10.


    * **use_structure_charge** (*bool*) – If set to True, then the public
    variable used for setting the overall charge of the
    structure (structure.charge) is used to set the NELECT
    variable in the INCAR
    Default is False (structure’s overall charge is not used)


    * **standardize** (*float*) – Whether to standardize to a primitive standard
    cell. Defaults to False.


    * **sym_prec** (*float*) – Tolerance for symmetry finding.


    * **international_monoclinic** (*bool*) – Whether to use international convention
    (vs Curtarolo) for monoclinic. Defaults True.


    * **validate_magmom** (*bool*) – Ensure that the missing magmom values are filled
    in with the VASP default value of 1.0



#### calculate_ng(max_prime_factor: int = 7, must_inc_2: bool = True)
Calculates the NGX, NGY, and NGZ values using the information available in the INCAR and POTCAR
This is meant to help with making initial guess for the FFT grid so we can interact with the Charge density API.


* **Parameters**


    * **max_prime_factor** (*int*) – the valid prime factors of the grid size in each direction
    VASP has many different setting for this to handle many compiling options.
    For typical MPI options all prime factors up to 7 are allowed


    * **must_inc_2** (*bool*) – Whether 2 must be a prime factor of the result. Defaults to True.



#### estimate_nbands()
Estimate the number of bands that VASP will initialize a
calculation with by default. Note that in practice this
can depend on # of cores (if not set explicitly).
Note that this formula is slightly different than the formula on the VASP wiki
(as of July 2023). This is because the formula in the source code (main.F) is
slightly different than what is on the wiki.


#### _property_ incar(_: [Incar](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Incar_ )
Incar


* **Type**

    return



#### _property_ kpoints(_: [Kpoints](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Kpoints) | Non_ )
Returns a KPOINTS file using the fully automated grid method. Uses
Gamma centered meshes for hexagonal cells and Monk grids otherwise.

If KSPACING is set in user_incar_settings (or the INCAR file), no
file is created because VASP will automatically generate the kpoints.

Algorithm:

    Uses a simple approach scaling the number of divisions along each
    reciprocal lattice vector proportional to its length.


#### _property_ nelect(_: floa_ )
Gets the default number of electrons for a given structure.


#### _property_ poscar(_: [Poscar](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Poscar_ )
Poscar


* **Type**

    return



#### _property_ potcar_functional(_: Literal['PBE', 'PBE_52', 'PBE_54', 'LDA', 'LDA_52', 'LDA_54', 'PW91', 'LDA_US', 'PW91_US'] | Non_ )
Returns the functional used for POTCAR generation.


#### _property_ structure(_: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure_ )
Structure


* **Type**

    return



#### write_input(output_dir: str, make_dir_if_not_present: bool = True, include_cif: bool = False, potcar_spec: bool = False, zip_output: bool = False)
Writes out all input to a directory.


* **Parameters**


    * **output_dir** (*str*) – Directory to output the VASP input files


    * **make_dir_if_not_present** (*bool*) – Set to True if you want the
    directory (and the whole path) to be created if it is not
    present.


    * **include_cif** (*bool*) – Whether to write a CIF file in the output
    directory for easier opening by VESTA.


    * **potcar_spec** (*bool*) – Instead of writing the POTCAR, write a “POTCAR.spec”.
    This is intended to help sharing an input set with people who might
    not have a license to specific Potcar files. Given a “POTCAR.spec”,
    the specific POTCAR file can be re-generated using pymatgen with the
    “generate_potcar” function in the pymatgen CLI.


    * **zip_output** (*bool*) – Whether to zip each VASP input file written to the output directory.



### _class_ pymatgen.io.vasp.sets.LobsterSet(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), isym: int = 0, ismear: int = -5, reciprocal_density: int | None = None, address_basis_file: str | None = None, user_supplied_basis: dict | None = None, \*\*kwargs)
Bases: `MPRelaxSet`

Input set to prepare VASP runs that can be digested by Lobster (See cohp.de).


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – input structure.


    * **isym** (*int*) – ISYM entry for INCAR, only isym=-1 and isym=0 are allowed


    * **ismear** (*int*) – ISMEAR entry for INCAR, only ismear=-5 and ismear=0 are allowed


    * **reciprocal_density** (*int*) – density of k-mesh by reciprocal volume


    * **user_supplied_basis** (*dict*) – dict including basis functions for all elements in structure,
    e.g. {“Fe”: “3d 3p 4s”, “O”: “2s 2p”}; if not supplied, a standard basis is used


    * **address_basis_file** (*str*) – address to a file similar to “BASIS_PBE_54_standaard.yaml”
    in pymatgen.io.lobster.lobster_basis


    * **user_potcar_settings** (*dict*) – dict including potcar settings for all elements in structure,
    e.g. {“Fe”: “Fe_pv”, “O”: “O”}; if not supplied, a standard basis is used.


    * **\*\*kwargs** – Other kwargs supported by `DictSet`.



#### CONFIG(_ = {'INCAR': {'ALGO': 'FAST', 'EDIFF_PER_ATOM': 5e-05, 'ENCUT': 520, 'IBRION': 2, 'ISIF': 3, 'ISMEAR': -5, 'ISPIN': 2, 'LASPH': True, 'LDAU': True, 'LDAUJ': {'F': {'Co': 0, 'Cr': 0, 'Fe': 0, 'Mn': 0, 'Mo': 0, 'Ni': 0, 'V': 0, 'W': 0}, 'O': {'Co': 0, 'Cr': 0, 'Fe': 0, 'Mn': 0, 'Mo': 0, 'Ni': 0, 'V': 0, 'W': 0}}, 'LDAUL': {'F': {'Co': 2, 'Cr': 2, 'Fe': 2, 'Mn': 2, 'Mo': 2, 'Ni': 2, 'V': 2, 'W': 2}, 'O': {'Co': 2, 'Cr': 2, 'Fe': 2, 'Mn': 2, 'Mo': 2, 'Ni': 2, 'V': 2, 'W': 2}}, 'LDAUPRINT': 1, 'LDAUTYPE': 2, 'LDAUU': {'F': {'Co': 3.32, 'Cr': 3.7, 'Fe': 5.3, 'Mn': 3.9, 'Mo': 4.38, 'Ni': 6.2, 'V': 3.25, 'W': 6.2}, 'O': {'Co': 3.32, 'Cr': 3.7, 'Fe': 5.3, 'Mn': 3.9, 'Mo': 4.38, 'Ni': 6.2, 'V': 3.25, 'W': 6.2}}, 'LORBIT': 11, 'LREAL': 'AUTO', 'LWAVE': False, 'MAGMOM': {'Ce': 5, 'Ce3+': 1, 'Co': 0.6, 'Co3+': 0.6, 'Co4+': 1, 'Cr': 5, 'Dy3+': 5, 'Er3+': 3, 'Eu': 10, 'Eu2+': 7, 'Eu3+': 6, 'Fe': 5, 'Gd3+': 7, 'Ho3+': 4, 'La3+': 0.6, 'Lu3+': 0.6, 'Mn': 5, 'Mn3+': 4, 'Mn4+': 3, 'Mo': 5, 'Nd3+': 3, 'Ni': 5, 'Pm3+': 4, 'Pr3+': 2, 'Sm3+': 5, 'Tb3+': 6, 'Tm3+': 2, 'V': 5, 'W': 5, 'Yb3+': 1}, 'NELM': 100, 'NSW': 99, 'PREC': 'Accurate', 'SIGMA': 0.05}, 'KPOINTS': {'reciprocal_density': 64}, 'PARENT': 'VASPIncarBase', 'POTCAR': {'Ac': 'Ac', 'Ag': 'Ag', 'Al': 'Al', 'Ar': 'Ar', 'As': 'As', 'Au': 'Au', 'B': 'B', 'Ba': 'Ba_sv', 'Be': 'Be_sv', 'Bi': 'Bi', 'Br': 'Br', 'C': 'C', 'Ca': 'Ca_sv', 'Cd': 'Cd', 'Ce': 'Ce', 'Cl': 'Cl', 'Co': 'Co', 'Cr': 'Cr_pv', 'Cs': 'Cs_sv', 'Cu': 'Cu_pv', 'Dy': 'Dy_3', 'Er': 'Er_3', 'Eu': 'Eu', 'F': 'F', 'Fe': 'Fe_pv', 'Ga': 'Ga_d', 'Gd': 'Gd', 'Ge': 'Ge_d', 'H': 'H', 'He': 'He', 'Hf': 'Hf_pv', 'Hg': 'Hg', 'Ho': 'Ho_3', 'I': 'I', 'In': 'In_d', 'Ir': 'Ir', 'K': 'K_sv', 'Kr': 'Kr', 'La': 'La', 'Li': 'Li_sv', 'Lu': 'Lu_3', 'Mg': 'Mg_pv', 'Mn': 'Mn_pv', 'Mo': 'Mo_pv', 'N': 'N', 'Na': 'Na_pv', 'Nb': 'Nb_pv', 'Nd': 'Nd_3', 'Ne': 'Ne', 'Ni': 'Ni_pv', 'Np': 'Np', 'O': 'O', 'Os': 'Os_pv', 'P': 'P', 'Pa': 'Pa', 'Pb': 'Pb_d', 'Pd': 'Pd', 'Pm': 'Pm_3', 'Pr': 'Pr_3', 'Pt': 'Pt', 'Pu': 'Pu', 'Rb': 'Rb_sv', 'Re': 'Re_pv', 'Rh': 'Rh_pv', 'Ru': 'Ru_pv', 'S': 'S', 'Sb': 'Sb', 'Sc': 'Sc_sv', 'Se': 'Se', 'Si': 'Si', 'Sm': 'Sm_3', 'Sn': 'Sn_d', 'Sr': 'Sr_sv', 'Ta': 'Ta_pv', 'Tb': 'Tb_3', 'Tc': 'Tc_pv', 'Te': 'Te', 'Th': 'Th', 'Ti': 'Ti_pv', 'Tl': 'Tl_d', 'Tm': 'Tm_3', 'U': 'U', 'V': 'V_pv', 'W': 'W_pv', 'Xe': 'Xe', 'Y': 'Y_sv', 'Yb': 'Yb_2', 'Zn': 'Zn', 'Zr': 'Zr_sv'}, 'POTCAR_FUNCTIONAL': 'PBE'_ )

### _class_ pymatgen.io.vasp.sets.MITMDSet(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), start_temp, end_temp, nsteps, time_step=2, spin_polarized=False, \*\*kwargs)
Bases: `MITRelaxSet`

Class for writing a vasp md run. This DOES NOT do multiple stage
runs.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Input structure.


    * **start_temp** (*int*) – Starting temperature.


    * **end_temp** (*int*) – Final temperature.


    * **nsteps** (*int*) – Number of time steps for simulations. NSW parameter.


    * **time_step** (*int*) – The time step for the simulation. The POTIM
    parameter. Defaults to 2fs.


    * **spin_polarized** (*bool*) – Whether to do spin polarized calculations.
    The ISPIN parameter. Defaults to False.


    * **\*\*kwargs** – Other kwargs supported by `DictSet`.



#### _property_ kpoints()
Kpoints


* **Type**

    return



### _class_ pymatgen.io.vasp.sets.MITNEBSet(structures, unset_encut=False, \*\*kwargs)
Bases: `MITRelaxSet`

Class for writing NEB inputs. Note that EDIFF is not on a per atom
basis for this input set.


* **Parameters**


    * **structures** – List of Structure objects.


    * **unset_encut** (*bool*) – Whether to unset ENCUT.


    * **\*\*kwargs** – Other kwargs supported by `DictSet`.



#### _property_ poscar()
Poscar for structure of first end point.


* **Type**

    return



#### _property_ poscars()
List of Poscars.


* **Type**

    return



#### write_input(output_dir, make_dir_if_not_present=True, write_cif=False, write_path_cif=False, write_endpoint_inputs=False)
NEB inputs has a special directory structure where inputs are in 00,
01, 02, ….


* **Parameters**


    * **output_dir** (*str*) – Directory to output the VASP input files


    * **make_dir_if_not_present** (*bool*) – Set to True if you want the
    directory (and the whole path) to be created if it is not
    present.


    * **write_cif** (*bool*) – If true, writes a cif along with each POSCAR.


    * **write_path_cif** (*bool*) – If true, writes a cif for each image.


    * **write_endpoint_inputs** (*bool*) – If true, writes input files for
    running endpoint calculations.



### _class_ pymatgen.io.vasp.sets.MITRelaxSet(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), \*\*kwargs)
Bases: `DictSet`

Standard implementation of VaspInputSet utilizing parameters in the MIT
High-throughput project.
The parameters are chosen specifically for a high-throughput project,
which means in general pseudopotentials with fewer electrons were chosen.

Please refer:

```default
A Jain, G. Hautier, C. Moore, S. P. Ong, C. Fischer, T. Mueller,
K. A. Persson, G. Ceder. A high-throughput infrastructure for density
functional theory calculations. Computational Materials Science,
2011, 50(8), 2295-2310. doi:10.1016/j.commatsci.2011.02.023
```


* **Parameters**


    * **structure** – Structure


    * **kwargs** – Same as those supported by DictSet.



#### CONFIG(_ = {'INCAR': {'ALGO': 'FAST', 'EDIFF': 1e-05, 'ENCUT': 520, 'IBRION': 2, 'ICHARG': 1, 'ISIF': 3, 'ISMEAR': -5, 'ISPIN': 2, 'ISYM': 0, 'LDAU': True, 'LDAUJ': {'F': {'Ag': 0, 'Co': 0, 'Cr': 0, 'Cu': 0, 'Fe': 0, 'Mn': 0, 'Mo': 0, 'Nb': 0, 'Ni': 0, 'Re': 0, 'Ta': 0, 'V': 0, 'W': 0}, 'O': {'Ag': 0, 'Co': 0, 'Cr': 0, 'Cu': 0, 'Fe': 0, 'Mn': 0, 'Mo': 0, 'Nb': 0, 'Ni': 0, 'Re': 0, 'Ta': 0, 'V': 0, 'W': 0}, 'S': {'Fe': 0, 'Mn': 0}}, 'LDAUL': {'F': {'Ag': 2, 'Co': 2, 'Cr': 2, 'Cu': 2, 'Fe': 2, 'Mn': 2, 'Mo': 2, 'Nb': 2, 'Ni': 2, 'Re': 2, 'Ta': 2, 'V': 2, 'W': 2}, 'O': {'Ag': 2, 'Co': 2, 'Cr': 2, 'Cu': 2, 'Fe': 2, 'Mn': 2, 'Mo': 2, 'Nb': 2, 'Ni': 2, 'Re': 2, 'Ta': 2, 'V': 2, 'W': 2}, 'S': {'Fe': 2, 'Mn': 2.5}}, 'LDAUPRINT': 1, 'LDAUTYPE': 2, 'LDAUU': {'F': {'Ag': 1.5, 'Co': 3.4, 'Cr': 3.5, 'Cu': 4, 'Fe': 4.0, 'Mn': 3.9, 'Mo': 4.38, 'Nb': 1.5, 'Ni': 6, 'Re': 2, 'Ta': 2, 'V': 3.1, 'W': 4.0}, 'O': {'Ag': 1.5, 'Co': 3.4, 'Cr': 3.5, 'Cu': 4, 'Fe': 4.0, 'Mn': 3.9, 'Mo': 4.38, 'Nb': 1.5, 'Ni': 6, 'Re': 2, 'Ta': 2, 'V': 3.1, 'W': 4.0}, 'S': {'Fe': 1.9, 'Mn': 2.5}}, 'LORBIT': '11', 'LREAL': 'AUTO', 'LWAVE': False, 'MAGMOM': {'Ce': 5, 'Ce3+': 1, 'Co': 0.6, 'Co3+': 0.6, 'Co4+': 1, 'Cr': 5, 'Dy3+': 5, 'Er3+': 3, 'Eu': 10, 'Eu2+': 7, 'Eu3+': 6, 'Fe': 5, 'Gd3+': 7, 'Ho3+': 4, 'La3+': 0.6, 'Lu3+': 0.6, 'Mn': 5, 'Mn3+': 4, 'Mn4+': 3, 'Mo': 5, 'Nd3+': 3, 'Ni': 5, 'Pm3+': 4, 'Pr3+': 2, 'Sm3+': 5, 'Tb3+': 6, 'Tm3+': 2, 'V': 5, 'W': 5, 'Yb3+': 1}, 'NELM': 200, 'NELMIN': 6, 'NSW': 99, 'PREC': 'Accurate', 'SIGMA': 0.05}, 'KPOINTS': {'length': 25}, 'PARENT': 'VASPIncarBase', 'POTCAR': {'Ac': {'hash': 'd6854224d20e3de6e6fd7399503791d1', 'symbol': 'Ac'}, 'Ag': {'hash': 'e8ffa02fe3f3a51338ac1ac91ae968b9', 'symbol': 'Ag'}, 'Al': {'hash': 'a6fd9a46aec185f4ad2acd0cbe4ae2fa', 'symbol': 'Al'}, 'Ar': {'hash': 'e782fc6292623b396091bf8b871c272f', 'symbol': 'Ar'}, 'As': {'hash': '8005364db225a254e52cba350bedd032', 'symbol': 'As'}, 'Au': {'hash': 'a9182d436a13194b744640ac940ab9b0', 'symbol': 'Au'}, 'B': {'hash': '18ed2875dfa6305324cec3d7d59273ae', 'symbol': 'B'}, 'Ba': {'hash': 'c0477913afb63dfae3439f3534fbf0ed', 'symbol': 'Ba_sv'}, 'Be': {'hash': 'fb974e44d56a8c62c6bbd1a1eb70c3a7', 'symbol': 'Be'}, 'Bi': {'hash': 'e29661c79d59abae3b3ba69eae24b1a5', 'symbol': 'Bi'}, 'Br': {'hash': '40f9594b4506684a69158c8975cfb9d6', 'symbol': 'Br'}, 'C': {'hash': 'c0a8167dbb174fe492a3db7f5006c0f8', 'symbol': 'C'}, 'Ca': {'hash': 'eb006721e214c04b3c13146e81b3a27d', 'symbol': 'Ca_sv'}, 'Cd': {'hash': '0506b2d0ac28d5fe2b5ced77a701aa86', 'symbol': 'Cd'}, 'Ce': {'hash': 'ff3a09f2ff91798e58eb4b9854e9be4a', 'symbol': 'Ce'}, 'Cl': {'hash': '779b9901046c78fe51c5d80224642aeb', 'symbol': 'Cl'}, 'Co': {'hash': 'b169bca4e137294d2ab3df8cbdd09083', 'symbol': 'Co'}, 'Cr': {'hash': '82c14307937c7509fda4e9bc023d243d', 'symbol': 'Cr'}, 'Cs': {'hash': '096b53a7d80cc0086976bcda50d536e5', 'symbol': 'Cs_sv'}, 'Cu': {'hash': '8ca4e43a30de0c397e51f16bbb20d678', 'symbol': 'Cu'}, 'Dy': {'hash': 'd4a05220ab0a2d4c03a76872ea724a1e', 'symbol': 'Dy_3'}, 'Er': {'hash': 'daa65a04877317f8c3c593ddeaa8a132', 'symbol': 'Er_3'}, 'Eu': {'hash': 'd466d046adf21f6146ee9644049ea268', 'symbol': 'Eu'}, 'F': {'hash': '180141c33d032bfbfff30b3bea9d23dd', 'symbol': 'F'}, 'Fe': {'hash': '9530da8244e4dac17580869b4adab115', 'symbol': 'Fe'}, 'Ga': {'hash': '6e0b9d58412b1bfcd7252aff13d476c2', 'symbol': 'Ga'}, 'Gd': {'hash': '1f0d42b1e5f6769d319d3f247992aeb9', 'symbol': 'Gd'}, 'Ge': {'hash': '79e788788c31e196a460553010512d3f', 'symbol': 'Ge'}, 'H': {'hash': 'bb43c666e3d36577264afe07669e9582', 'symbol': 'H'}, 'He': {'hash': '47f9434aa3db96c85d7c4b3e4c2df09b', 'symbol': 'He'}, 'Hf': {'hash': 'b113f150cbf9c736f8244a6c25b0482e', 'symbol': 'Hf'}, 'Hg': {'hash': 'c2f15dfb5fd53396c5427635e5019160', 'symbol': 'Hg'}, 'Ho': {'hash': '661891464a27e87cf7e1324dd1893b77', 'symbol': 'Ho_3'}, 'I': {'hash': 'f4ff16a495dd361ff5824ee61b418bb0', 'symbol': 'I'}, 'In': {'hash': '7df38c0cdb4e6d9a9b93f09d690bb3ae', 'symbol': 'In'}, 'Ir': {'hash': 'dbcf7dcc6f4fb40df7b3d26904f60a66', 'symbol': 'Ir'}, 'K': {'hash': '3e84f86d37f203a4fb01de36af57e430', 'symbol': 'K_sv'}, 'Kr': {'hash': '39b9b85ae3982e6c012fb549b2840ce5', 'symbol': 'Kr'}, 'La': {'hash': '9b3ce03d18f7c0b40471a817ff91b287', 'symbol': 'La'}, 'Li': {'hash': '65e83282d1707ec078c1012afbd05be8', 'symbol': 'Li'}, 'Lu': {'hash': 'd40a90babf1224b88ffb4c3273ac3848', 'symbol': 'Lu_3'}, 'Mg': {'hash': '1771eb72adbbfa6310d66e7517e49930', 'symbol': 'Mg'}, 'Mn': {'hash': 'd082dba29b57ab59b3165e605dbf71b8', 'symbol': 'Mn'}, 'Mo': {'hash': '84e18fd84a98e3d7fa8f055952410df0', 'symbol': 'Mo_pv'}, 'N': {'hash': 'b98fd027ddebc67da4063ff2cabbc04b', 'symbol': 'N'}, 'Na': {'hash': '1a89e79f7e21d99e8cf5788979f6a987', 'symbol': 'Na'}, 'Nb': {'hash': '7bcee99a4dc3094be0f9fd7961c02966', 'symbol': 'Nb_pv'}, 'Nd': {'hash': '0c64e63070cee837c967283fffa001df', 'symbol': 'Nd'}, 'Ne': {'hash': '52064eee378b9e37a295a674f1c278f0', 'symbol': 'Ne'}, 'Ni': {'hash': '653f5772e68b2c7fd87ffd1086c0d710', 'symbol': 'Ni'}, 'Np': {'hash': '20cb30b714200c4db870550b288ac4cd', 'symbol': 'Np'}, 'O': {'hash': '7a25bc5b9a5393f46600a4939d357982', 'symbol': 'O'}, 'Os': {'hash': '35c2cb48d48a9c38c40fb82bbe70626d', 'symbol': 'Os'}, 'P': {'hash': '7dc3393307131ae67785a0cdacb61d5f', 'symbol': 'P'}, 'Pa': {'hash': 'a1fdb1089d0727f415416ec8082246ba', 'symbol': 'Pa'}, 'Pb': {'hash': '704c2c967247d7f84090d2536c91877d', 'symbol': 'Pb'}, 'Pd': {'hash': 'a395eb3aaf2fcab12fac3030a1146f61', 'symbol': 'Pd'}, 'Pm': {'hash': 'a2c9485ea86b2a7cf175077e6e5c7b3e', 'symbol': 'Pm'}, 'Pr': {'hash': '92f191499bf5346ea652bb806350ad87', 'symbol': 'Pr'}, 'Pt': {'hash': 'a604ea3c6a9cc23c739b762f625cf449', 'symbol': 'Pt'}, 'Pu': {'hash': 'f1d01e845dccc52d448679911f301a73', 'symbol': 'Pu'}, 'Rb': {'hash': 'e447c648d870b066b3514e6b800727ab', 'symbol': 'Rb_pv'}, 'Re': {'hash': '72385e193c92a8acfe17ea49004c2be1', 'symbol': 'Re'}, 'Rh': {'hash': '2c3dba3fcc6058ca1b1cfa75e45084bc', 'symbol': 'Rh'}, 'Ru': {'hash': '7925f4d4b68076d70af7cd86eef9ba8d', 'symbol': 'Ru_pv'}, 'S': {'hash': 'd368db6899d8839859bbee4811a42a88', 'symbol': 'S'}, 'Sb': {'hash': 'd82c022b02fc5344e85bd1909f9ee3e7', 'symbol': 'Sb'}, 'Sc': {'hash': 'dc386f505ad0c43385a7715b4111cb75', 'symbol': 'Sc_sv'}, 'Se': {'hash': '67a8804ede9f1112726e3d136978ef19', 'symbol': 'Se'}, 'Si': {'hash': 'b2b0ea6feb62e7cde209616683b8f7f5', 'symbol': 'Si'}, 'Sm': {'hash': 'e5e274e7cd99602ca81d146155abdf88', 'symbol': 'Sm_3'}, 'Sn': {'hash': '849b0795e148f93113a06be8fd5f5001', 'symbol': 'Sn_d'}, 'Sr': {'hash': 'ca6a5429c120a0ab705824386a76fe5b', 'symbol': 'Sr_sv'}, 'Ta': {'hash': 'd4e2cfe9338ef80da592d5bb9dc782c7', 'symbol': 'Ta'}, 'Tb': {'hash': '0790955c547003956c0fd4f080f7f508', 'symbol': 'Tb_3'}, 'Tc': {'hash': '9592642886319309a39d55c5717c6f48', 'symbol': 'Tc'}, 'Te': {'hash': '72719856e22fb1d3032df6f96d98a0f2', 'symbol': 'Te'}, 'Th': {'hash': 'aea79f322180fa6f0bfa74cb2a156dcf', 'symbol': 'Th'}, 'Ti': {'hash': 'c617e8b539c3f44a0ab6e8da2a92d318', 'symbol': 'Ti'}, 'Tl': {'hash': '2aa0d5406aaab7ebfbc761da382f1352', 'symbol': 'Tl'}, 'Tm': {'hash': '94a07cb7949b01305cb161da0cbfb492', 'symbol': 'Tm_3'}, 'U': {'hash': '72702eabbb1bc02b4167590dc848ed5d', 'symbol': 'U'}, 'V': {'hash': '7f1297a2e1d963e2a4d81b61f85e4ded', 'symbol': 'V_pv'}, 'W': {'hash': '2a33e0d5c700640535f60ac0a12177ab', 'symbol': 'W_pv'}, 'Xe': {'hash': '338472e581f58b41d37c002a5e22353b', 'symbol': 'Xe'}, 'Y': {'hash': '4ed187e77cd54f198bb88020278b143d', 'symbol': 'Y_sv'}, 'Yb': {'hash': '9f472bd422f640710f7d93e2d9ce89f4', 'symbol': 'Yb'}, 'Zn': {'hash': 'e35ee27f8483a63bb68dbc236a343af3', 'symbol': 'Zn'}, 'Zr': {'hash': 'd221d2c0bac4f8e81af2f5c42a314274', 'symbol': 'Zr'}}, 'POTCAR_FUNCTIONAL': 'PBE'_ )

### _class_ pymatgen.io.vasp.sets.MPAbsorptionSet(structure, mode='IPA', copy_wavecar=True, nbands=None, nbands_factor=2, reciprocal_density=400, nkred=None, nedos=2001, prev_incar=None, \*\*kwargs)
Bases: `MPRelaxSet`

MP input set for generating frequency dependent dielectrics.
Two modes are supported: “IPA” or “RPA”.
A typical sequence is mode=”STATIC” -> mode=”IPA” -> mode=”RPA”(optional)
For all steps other than the first one (static), the
recommendation is to use from_prev_calculation on the preceding run in
the series. It is important to ensure Gamma centred kpoints for the RPA step.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Input structure.


    * **prev_incar** (*Incar/string*) – Incar file from previous run.


    * **mode** (*str*) – Supported modes are “IPA”, “RPA”


    * **copy_wavecar** (*bool*) – Whether to copy the WAVECAR from a previous run. Defaults to True.


    * **nbands** (*int*) – For subsequent calculations, it is generally
    recommended to perform NBANDS convergence starting from the
    NBANDS of the previous run for DIAG, and to use the exact same
    NBANDS for RPA. This parameter is used by
    from_previous_calculation to set nband.


    * **nbands_factor** (*int*) – Multiplicative factor for NBANDS when starting
    from a previous calculation. Only applies if mode==”IPA”.
    Need to be tested for convergence.


    * **reciprocal_density** – the k-points density


    * **nkred** – the reduced number of kpoints to calculate, equal to the k-mesh. Only applies in “RPA” mode
    because of the q->0 limit.


    * **nedos** – the density of DOS, default: 2001.


    * **\*\*kwargs** – All kwargs supported by DictSet. Typically, user_incar_settings is a commonly used option.



#### SUPPORTED_MODES(_ = ('IPA', 'RPA'_ )

#### _classmethod_ from_prev_calc(prev_calc_dir, mode, \*\*kwargs)
Generate a set of VASP input files for absorption calculation
:param prev_calc_dir: The directory contains the outputs(

> vasprun.xml of previous vasp run.


* **Parameters**


    * **mode** (*str*) – Supported modes are “IPA”, “RPA” (default)


    * **\*\*kwargs** – All kwargs supported by MPAbsorptionsSet, other than structure.



#### _property_ incar()
Incar


* **Type**

    return



#### _property_ kpoints()
Generate gamma center k-points mesh grid for optical calculation. It is not mandatory for ‘ALGO = Exact’,
but is requested by ‘ALGO = CHI’ calculation.


#### override_from_prev_calc(prev_calc_dir='.', \*\*kwargs)
Update the input set to include settings from a previous calculation.


* **Parameters**


    * **prev_calc_dir** (*str*) – The path to the previous calculation directory.


    * **\*\*kwargs** – unused



* **Returns**

    The input set with the settings (structure, k-points, incar, etc)
    updated using the previous VASP run.



### _class_ pymatgen.io.vasp.sets.MPHSEBSSet(structure, user_incar_settings=None, added_kpoints=None, mode='Gap', reciprocal_density=None, copy_chgcar=True, kpoints_line_density=20, \*\*kwargs)
Bases: `MPHSERelaxSet`

Implementation of a VaspInputSet for HSE band structure computations.
Remember that HSE band structures must be self-consistent in VASP. A
band structure along symmetry lines for instance needs BOTH a uniform
grid with appropriate weights AND a path along the lines with weight 0.

Thus, the “Uniform” mode is just like regular static SCF but allows
adding custom kpoints (e.g., corresponding to known VBM/CBM) to the
uniform grid that have zero weight (e.g., for better gap estimate).

The “Gap” mode behaves just like the “Uniform” mode, however, if starting
from a previous calculation, the VBM and CBM k-points will automatically
be added to `added_kpoints`.

The “Line” mode is just like Uniform mode, but additionally adds
k-points along symmetry lines with zero weight.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Structure to compute


    * **user_incar_settings** (*dict*) – A dict specifying additional incar
    settings


    * **added_kpoints** (*list*) – a list of kpoints (list of 3 number list)
    added to the run. The k-points are in fractional coordinates


    * **mode** (*str*) – “Line” - generate k-points along symmetry lines for
    bandstructure. “Uniform” - generate uniform k-points grid.


    * **reciprocal_density** (*int*) – k-point density to use for uniform mesh.


    * **copy_chgcar** (*bool*) – Whether to copy the CHGCAR of a previous run.


    * **kpoints_line_density** (*int*) – k-point density for high symmetry lines


    * **\*\*kwargs** (*dict*) – Any other parameters to pass into DictSet.



#### _classmethod_ from_prev_calc(prev_calc_dir, \*\*kwargs)
Generate a set of VASP input files for HSE calculations from a
directory of previous VASP run.


* **Parameters**


    * **prev_calc_dir** (*str*) – Directory containing the outputs
    (vasprun.xml and OUTCAR) of previous vasp run.


    * **\*\*kwargs** – All kwargs supported by MPHSEBSStaticSet, other than
    prev_structure which is determined from the previous calc dir.



#### _property_ kpoints(_: [Kpoints](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Kpoints_ )
Kpoints


* **Type**

    return



#### override_from_prev_calc(prev_calc_dir='.')
Update the input set to include settings from a previous calculation.


* **Parameters**

    **prev_calc_dir** (*str*) – The path to the previous calculation directory.



* **Returns**

    The input set with the settings (structure, k-points, incar, etc)
    updated using the previous VASP run.



### _class_ pymatgen.io.vasp.sets.MPHSERelaxSet(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), \*\*kwargs)
Bases: `DictSet`

Same as the MPRelaxSet, but with HSE parameters.


* **Parameters**


    * **structure** – Structure


    * **kwargs** – Same as those supported by DictSet.



#### CONFIG(_ = {'INCAR': {'ALGO': 'All', 'EDIFF_PER_ATOM': 5e-05, 'ENCUT': 520, 'HFSCREEN': 0.2, 'IBRION': 2, 'ICHARG': 1, 'ISIF': 3, 'ISMEAR': 0, 'ISPIN': 2, 'LHFCALC': True, 'LORBIT': 11, 'LREAL': 'AUTO', 'LWAVE': False, 'MAGMOM': {'Ce': 5, 'Ce3+': 1, 'Co': 0.6, 'Co3+': 0.6, 'Co4+': 1, 'Cr': 5, 'Dy3+': 5, 'Er3+': 3, 'Eu': 10, 'Eu2+': 7, 'Eu3+': 6, 'Fe': 5, 'Gd3+': 7, 'Ho3+': 4, 'La3+': 0.6, 'Lu3+': 0.6, 'Mn': 5, 'Mn3+': 4, 'Mn4+': 3, 'Mo': 5, 'Nd3+': 3, 'Ni': 5, 'Pm3+': 4, 'Pr3+': 2, 'Sm3+': 5, 'Tb3+': 6, 'Tm3+': 2, 'V': 5, 'W': 5, 'Yb3+': 1}, 'NELM': 100, 'NSW': 99, 'PREC': 'Accurate', 'PRECFOCK': 'Fast', 'SIGMA': 0.05}, 'KPOINTS': {'reciprocal_density': 50}, 'PARENT': 'VASPIncarBase', 'POTCAR': {'Ac': 'Ac', 'Ag': 'Ag', 'Al': 'Al', 'Ar': 'Ar', 'As': 'As', 'Au': 'Au', 'B': 'B', 'Ba': 'Ba_sv', 'Be': 'Be_sv', 'Bi': 'Bi', 'Br': 'Br', 'C': 'C', 'Ca': 'Ca_sv', 'Cd': 'Cd', 'Ce': 'Ce', 'Cl': 'Cl', 'Co': 'Co', 'Cr': 'Cr_pv', 'Cs': 'Cs_sv', 'Cu': 'Cu_pv', 'Dy': 'Dy_3', 'Er': 'Er_3', 'Eu': 'Eu', 'F': 'F', 'Fe': 'Fe_pv', 'Ga': 'Ga_d', 'Gd': 'Gd', 'Ge': 'Ge_d', 'H': 'H', 'He': 'He', 'Hf': 'Hf_pv', 'Hg': 'Hg', 'Ho': 'Ho_3', 'I': 'I', 'In': 'In_d', 'Ir': 'Ir', 'K': 'K_sv', 'Kr': 'Kr', 'La': 'La', 'Li': 'Li_sv', 'Lu': 'Lu_3', 'Mg': 'Mg_pv', 'Mn': 'Mn_pv', 'Mo': 'Mo_pv', 'N': 'N', 'Na': 'Na_pv', 'Nb': 'Nb_pv', 'Nd': 'Nd_3', 'Ne': 'Ne', 'Ni': 'Ni_pv', 'Np': 'Np', 'O': 'O', 'Os': 'Os_pv', 'P': 'P', 'Pa': 'Pa', 'Pb': 'Pb_d', 'Pd': 'Pd', 'Pm': 'Pm_3', 'Pr': 'Pr_3', 'Pt': 'Pt', 'Pu': 'Pu', 'Rb': 'Rb_sv', 'Re': 'Re_pv', 'Rh': 'Rh_pv', 'Ru': 'Ru_pv', 'S': 'S', 'Sb': 'Sb', 'Sc': 'Sc_sv', 'Se': 'Se', 'Si': 'Si', 'Sm': 'Sm_3', 'Sn': 'Sn_d', 'Sr': 'Sr_sv', 'Ta': 'Ta_pv', 'Tb': 'Tb_3', 'Tc': 'Tc_pv', 'Te': 'Te', 'Th': 'Th', 'Ti': 'Ti_pv', 'Tl': 'Tl_d', 'Tm': 'Tm_3', 'U': 'U', 'V': 'V_pv', 'W': 'W_pv', 'Xe': 'Xe', 'Y': 'Y_sv', 'Yb': 'Yb_2', 'Zn': 'Zn', 'Zr': 'Zr_sv'}, 'POTCAR_FUNCTIONAL': 'PBE_52'_ )

### _class_ pymatgen.io.vasp.sets.MPMDSet(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), start_temp, end_temp, nsteps, spin_polarized=False, \*\*kwargs)
Bases: `MPRelaxSet`

This a modified version of the old MITMDSet pre 2018/03/12.

This set serves as the basis for the amorphous skyline paper.


1. Aykol, M.; Dwaraknath, S. S.; Sun, W.; Persson, K. A. Thermodynamic
Limit for Synthesis of Metastable Inorganic Materials. Sci. Adv. 2018,
4 (4).

Class for writing a vasp md run. This DOES NOT do multiple stage runs.
Precision remains normal, to increase accuracy of stress tensor.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Input structure.


    * **start_temp** (*int*) – Starting temperature.


    * **end_temp** (*int*) – Final temperature.


    * **nsteps** (*int*) – Number of time steps for simulations. NSW parameter.


    * **time_step** (*int*) – The time step for the simulation. The POTIM
    parameter. Defaults to 2fs.


    * **spin_polarized** (*bool*) – Whether to do spin polarized calculations.
    The ISPIN parameter. Defaults to False.


    * **\*\*kwargs** – Other kwargs supported by `DictSet`.



#### _property_ kpoints()
Kpoints


* **Type**

    return



### _class_ pymatgen.io.vasp.sets.MPMetalRelaxSet(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), \*\*kwargs)
Bases: `MPRelaxSet`

Implementation of VaspInputSet utilizing parameters in the public
Materials Project, but with tuning for metals. Key things are a denser
k point density, and a.


* **Parameters**


    * **structure** – Structure


    * **kwargs** – Same as those supported by DictSet.



#### CONFIG(_ = {'INCAR': {'ALGO': 'FAST', 'EDIFF_PER_ATOM': 5e-05, 'ENCUT': 520, 'IBRION': 2, 'ISIF': 3, 'ISMEAR': -5, 'ISPIN': 2, 'LASPH': True, 'LDAU': True, 'LDAUJ': {'F': {'Co': 0, 'Cr': 0, 'Fe': 0, 'Mn': 0, 'Mo': 0, 'Ni': 0, 'V': 0, 'W': 0}, 'O': {'Co': 0, 'Cr': 0, 'Fe': 0, 'Mn': 0, 'Mo': 0, 'Ni': 0, 'V': 0, 'W': 0}}, 'LDAUL': {'F': {'Co': 2, 'Cr': 2, 'Fe': 2, 'Mn': 2, 'Mo': 2, 'Ni': 2, 'V': 2, 'W': 2}, 'O': {'Co': 2, 'Cr': 2, 'Fe': 2, 'Mn': 2, 'Mo': 2, 'Ni': 2, 'V': 2, 'W': 2}}, 'LDAUPRINT': 1, 'LDAUTYPE': 2, 'LDAUU': {'F': {'Co': 3.32, 'Cr': 3.7, 'Fe': 5.3, 'Mn': 3.9, 'Mo': 4.38, 'Ni': 6.2, 'V': 3.25, 'W': 6.2}, 'O': {'Co': 3.32, 'Cr': 3.7, 'Fe': 5.3, 'Mn': 3.9, 'Mo': 4.38, 'Ni': 6.2, 'V': 3.25, 'W': 6.2}}, 'LORBIT': 11, 'LREAL': 'AUTO', 'LWAVE': False, 'MAGMOM': {'Ce': 5, 'Ce3+': 1, 'Co': 0.6, 'Co3+': 0.6, 'Co4+': 1, 'Cr': 5, 'Dy3+': 5, 'Er3+': 3, 'Eu': 10, 'Eu2+': 7, 'Eu3+': 6, 'Fe': 5, 'Gd3+': 7, 'Ho3+': 4, 'La3+': 0.6, 'Lu3+': 0.6, 'Mn': 5, 'Mn3+': 4, 'Mn4+': 3, 'Mo': 5, 'Nd3+': 3, 'Ni': 5, 'Pm3+': 4, 'Pr3+': 2, 'Sm3+': 5, 'Tb3+': 6, 'Tm3+': 2, 'V': 5, 'W': 5, 'Yb3+': 1}, 'NELM': 100, 'NSW': 99, 'PREC': 'Accurate', 'SIGMA': 0.05}, 'KPOINTS': {'reciprocal_density': 64}, 'PARENT': 'VASPIncarBase', 'POTCAR': {'Ac': 'Ac', 'Ag': 'Ag', 'Al': 'Al', 'Ar': 'Ar', 'As': 'As', 'Au': 'Au', 'B': 'B', 'Ba': 'Ba_sv', 'Be': 'Be_sv', 'Bi': 'Bi', 'Br': 'Br', 'C': 'C', 'Ca': 'Ca_sv', 'Cd': 'Cd', 'Ce': 'Ce', 'Cl': 'Cl', 'Co': 'Co', 'Cr': 'Cr_pv', 'Cs': 'Cs_sv', 'Cu': 'Cu_pv', 'Dy': 'Dy_3', 'Er': 'Er_3', 'Eu': 'Eu', 'F': 'F', 'Fe': 'Fe_pv', 'Ga': 'Ga_d', 'Gd': 'Gd', 'Ge': 'Ge_d', 'H': 'H', 'He': 'He', 'Hf': 'Hf_pv', 'Hg': 'Hg', 'Ho': 'Ho_3', 'I': 'I', 'In': 'In_d', 'Ir': 'Ir', 'K': 'K_sv', 'Kr': 'Kr', 'La': 'La', 'Li': 'Li_sv', 'Lu': 'Lu_3', 'Mg': 'Mg_pv', 'Mn': 'Mn_pv', 'Mo': 'Mo_pv', 'N': 'N', 'Na': 'Na_pv', 'Nb': 'Nb_pv', 'Nd': 'Nd_3', 'Ne': 'Ne', 'Ni': 'Ni_pv', 'Np': 'Np', 'O': 'O', 'Os': 'Os_pv', 'P': 'P', 'Pa': 'Pa', 'Pb': 'Pb_d', 'Pd': 'Pd', 'Pm': 'Pm_3', 'Pr': 'Pr_3', 'Pt': 'Pt', 'Pu': 'Pu', 'Rb': 'Rb_sv', 'Re': 'Re_pv', 'Rh': 'Rh_pv', 'Ru': 'Ru_pv', 'S': 'S', 'Sb': 'Sb', 'Sc': 'Sc_sv', 'Se': 'Se', 'Si': 'Si', 'Sm': 'Sm_3', 'Sn': 'Sn_d', 'Sr': 'Sr_sv', 'Ta': 'Ta_pv', 'Tb': 'Tb_3', 'Tc': 'Tc_pv', 'Te': 'Te', 'Th': 'Th', 'Ti': 'Ti_pv', 'Tl': 'Tl_d', 'Tm': 'Tm_3', 'U': 'U', 'V': 'V_pv', 'W': 'W_pv', 'Xe': 'Xe', 'Y': 'Y_sv', 'Yb': 'Yb_2', 'Zn': 'Zn', 'Zr': 'Zr_sv'}, 'POTCAR_FUNCTIONAL': 'PBE'_ )

### _class_ pymatgen.io.vasp.sets.MPNMRSet(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), mode: Literal['cs', 'efg'] = 'cs', isotopes: list | None = None, prev_incar: [Incar](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Incar) = None, reciprocal_density: int = 100, \*\*kwargs)
Bases: `MPStaticSet`

Init a MPNMRSet.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Structure to compute


    * **mode** (*str*) – The NMR calculation to run
    “cs”: for Chemical Shift
    “efg” for Electric Field Gradient


    * **isotopes** (*list*) – list of Isotopes for quadrupole moments


    * **prev_incar** ([*Incar*](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Incar)) – Incar file from previous run.


    * **reciprocal_density** (*int*) – density of k-mesh by reciprocal volume. Defaults to 100.


    * **\*\*kwargs** – kwargs supported by MPStaticSet.



#### _property_ incar()
Incar


* **Type**

    return



### _class_ pymatgen.io.vasp.sets.MPNonSCFSet(structure, prev_incar=None, mode='line', nedos=2001, dedos=0.005, reciprocal_density=100, sym_prec=0.1, kpoints_line_density=20, optics=False, copy_chgcar=True, nbands_factor=1.2, small_gap_multiply=None, \*\*kwargs)
Bases: `MPRelaxSet`

Init a MPNonSCFSet. Typically, you would use the classmethod
from_prev_calc to initialize from a previous SCF run.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Structure to compute


    * **prev_incar** (*Incar/string*) – Incar file from previous run.


    * **mode** (*str*) – Line, Uniform or Boltztrap mode supported.


    * **nedos** (*int*) – nedos parameter. Default to 2001.


    * **dedos** (*float*) – setting nedos=0 and uniform mode in from_prev_calc,
    an automatic nedos will be calculated using the total energy range
    divided by the energy step dedos


    * **reciprocal_density** (*int*) – density of k-mesh by reciprocal
    volume (defaults to 100)


    * **sym_prec** (*float*) – Symmetry precision (for Uniform mode).


    * **kpoints_line_density** (*int*) – Line density for Line mode.


    * **optics** (*bool*) – whether to add dielectric function


    * **copy_chgcar** – Whether to copy the old CHGCAR when starting from a
    previous calculation.


    * **nbands_factor** (*float*) – Multiplicative factor for NBANDS when starting
    from a previous calculation. Choose a higher number if you are
    doing an LOPTICS calculation.


    * **small_gap_multiply** (*[**float**, **float**]*) – When starting from a previous
    calculation, if the gap is less than 1st index, multiply the default
    reciprocal_density by the 2nd index.


    * **\*\*kwargs** – kwargs supported by MPRelaxSet.



#### _classmethod_ from_prev_calc(prev_calc_dir, \*\*kwargs)
Generate a set of VASP input files for NonSCF calculations from a
directory of previous static VASP run.


* **Parameters**


    * **prev_calc_dir** (*str*) – The directory contains the outputs(
    vasprun.xml and OUTCAR) of previous vasp run.


    * **\*\*kwargs** – All kwargs supported by MPNonSCFSet, other than structure,
    prev_incar and prev_chgcar which are determined from the
    prev_calc_dir.



#### _property_ incar(_: [Incar](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Incar_ )
Incar


* **Type**

    return



#### _property_ kpoints(_: [Kpoints](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Kpoints) | Non_ )
Kpoints


* **Type**

    return



#### override_from_prev_calc(prev_calc_dir='.')
Update the input set to include settings from a previous calculation.


* **Parameters**

    **prev_calc_dir** (*str*) – The path to the previous calculation directory.



* **Returns**

    The input set with the settings (structure, k-points, incar, etc)
    updated using the previous VASP run.



### _class_ pymatgen.io.vasp.sets.MPRelaxSet(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), \*\*kwargs)
Bases: `DictSet`

Implementation of VaspInputSet utilizing parameters in the public
Materials Project. Typically, the pseudopotentials chosen contain more
electrons than the MIT parameters, and the k-point grid is ~50% more dense.
The LDAUU parameters are also different due to the different psps used,
which result in different fitted values.


* **Parameters**


    * **structure** – Structure


    * **kwargs** – Same as those supported by DictSet.



#### CONFIG(_ = {'INCAR': {'ALGO': 'FAST', 'EDIFF_PER_ATOM': 5e-05, 'ENCUT': 520, 'IBRION': 2, 'ISIF': 3, 'ISMEAR': -5, 'ISPIN': 2, 'LASPH': True, 'LDAU': True, 'LDAUJ': {'F': {'Co': 0, 'Cr': 0, 'Fe': 0, 'Mn': 0, 'Mo': 0, 'Ni': 0, 'V': 0, 'W': 0}, 'O': {'Co': 0, 'Cr': 0, 'Fe': 0, 'Mn': 0, 'Mo': 0, 'Ni': 0, 'V': 0, 'W': 0}}, 'LDAUL': {'F': {'Co': 2, 'Cr': 2, 'Fe': 2, 'Mn': 2, 'Mo': 2, 'Ni': 2, 'V': 2, 'W': 2}, 'O': {'Co': 2, 'Cr': 2, 'Fe': 2, 'Mn': 2, 'Mo': 2, 'Ni': 2, 'V': 2, 'W': 2}}, 'LDAUPRINT': 1, 'LDAUTYPE': 2, 'LDAUU': {'F': {'Co': 3.32, 'Cr': 3.7, 'Fe': 5.3, 'Mn': 3.9, 'Mo': 4.38, 'Ni': 6.2, 'V': 3.25, 'W': 6.2}, 'O': {'Co': 3.32, 'Cr': 3.7, 'Fe': 5.3, 'Mn': 3.9, 'Mo': 4.38, 'Ni': 6.2, 'V': 3.25, 'W': 6.2}}, 'LORBIT': 11, 'LREAL': 'AUTO', 'LWAVE': False, 'MAGMOM': {'Ce': 5, 'Ce3+': 1, 'Co': 0.6, 'Co3+': 0.6, 'Co4+': 1, 'Cr': 5, 'Dy3+': 5, 'Er3+': 3, 'Eu': 10, 'Eu2+': 7, 'Eu3+': 6, 'Fe': 5, 'Gd3+': 7, 'Ho3+': 4, 'La3+': 0.6, 'Lu3+': 0.6, 'Mn': 5, 'Mn3+': 4, 'Mn4+': 3, 'Mo': 5, 'Nd3+': 3, 'Ni': 5, 'Pm3+': 4, 'Pr3+': 2, 'Sm3+': 5, 'Tb3+': 6, 'Tm3+': 2, 'V': 5, 'W': 5, 'Yb3+': 1}, 'NELM': 100, 'NSW': 99, 'PREC': 'Accurate', 'SIGMA': 0.05}, 'KPOINTS': {'reciprocal_density': 64}, 'PARENT': 'VASPIncarBase', 'POTCAR': {'Ac': 'Ac', 'Ag': 'Ag', 'Al': 'Al', 'Ar': 'Ar', 'As': 'As', 'Au': 'Au', 'B': 'B', 'Ba': 'Ba_sv', 'Be': 'Be_sv', 'Bi': 'Bi', 'Br': 'Br', 'C': 'C', 'Ca': 'Ca_sv', 'Cd': 'Cd', 'Ce': 'Ce', 'Cl': 'Cl', 'Co': 'Co', 'Cr': 'Cr_pv', 'Cs': 'Cs_sv', 'Cu': 'Cu_pv', 'Dy': 'Dy_3', 'Er': 'Er_3', 'Eu': 'Eu', 'F': 'F', 'Fe': 'Fe_pv', 'Ga': 'Ga_d', 'Gd': 'Gd', 'Ge': 'Ge_d', 'H': 'H', 'He': 'He', 'Hf': 'Hf_pv', 'Hg': 'Hg', 'Ho': 'Ho_3', 'I': 'I', 'In': 'In_d', 'Ir': 'Ir', 'K': 'K_sv', 'Kr': 'Kr', 'La': 'La', 'Li': 'Li_sv', 'Lu': 'Lu_3', 'Mg': 'Mg_pv', 'Mn': 'Mn_pv', 'Mo': 'Mo_pv', 'N': 'N', 'Na': 'Na_pv', 'Nb': 'Nb_pv', 'Nd': 'Nd_3', 'Ne': 'Ne', 'Ni': 'Ni_pv', 'Np': 'Np', 'O': 'O', 'Os': 'Os_pv', 'P': 'P', 'Pa': 'Pa', 'Pb': 'Pb_d', 'Pd': 'Pd', 'Pm': 'Pm_3', 'Pr': 'Pr_3', 'Pt': 'Pt', 'Pu': 'Pu', 'Rb': 'Rb_sv', 'Re': 'Re_pv', 'Rh': 'Rh_pv', 'Ru': 'Ru_pv', 'S': 'S', 'Sb': 'Sb', 'Sc': 'Sc_sv', 'Se': 'Se', 'Si': 'Si', 'Sm': 'Sm_3', 'Sn': 'Sn_d', 'Sr': 'Sr_sv', 'Ta': 'Ta_pv', 'Tb': 'Tb_3', 'Tc': 'Tc_pv', 'Te': 'Te', 'Th': 'Th', 'Ti': 'Ti_pv', 'Tl': 'Tl_d', 'Tm': 'Tm_3', 'U': 'U', 'V': 'V_pv', 'W': 'W_pv', 'Xe': 'Xe', 'Y': 'Y_sv', 'Yb': 'Yb_2', 'Zn': 'Zn', 'Zr': 'Zr_sv'}, 'POTCAR_FUNCTIONAL': 'PBE'_ )

### _class_ pymatgen.io.vasp.sets.MPSOCSet(structure, saxis=(0, 0, 1), copy_chgcar=True, nbands_factor=1.2, reciprocal_density=100, small_gap_multiply=None, magmom=None, \*\*kwargs)
Bases: `MPStaticSet`

An input set for running spin-orbit coupling (SOC) calculations.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – the structure must have the ‘magmom’ site
    property and each magnetic moment value must have 3
    components. eg: `magmom = [[0,0,2], ...]`


    * **saxis** (*tuple*) – magnetic moment orientation


    * **copy_chgcar** – Whether to copy the old CHGCAR. Defaults to True.


    * **nbands_factor** (*float*) – Multiplicative factor for NBANDS. Choose a
    higher number if you are doing an LOPTICS calculation.


    * **reciprocal_density** (*int*) – density of k-mesh by reciprocal volume.


    * **small_gap_multiply** (*[**float**, **float**]*) – If the gap is less than
    1st index, multiply the default reciprocal_density by the 2nd
    index.


    * **magmom** (*list**[**list**[**float**]**]*) – Override for the structure magmoms.


    * **\*\*kwargs** – kwargs supported by MPStaticSet.



#### _classmethod_ from_prev_calc(prev_calc_dir, \*\*kwargs)
Generate a set of VASP input files for SOC calculations from a
directory of previous static VASP run. SOC calc requires all 3
components for MAGMOM for each atom in the structure.


* **Parameters**


    * **prev_calc_dir** (*str*) – The directory contains the outputs(
    vasprun.xml and OUTCAR) of previous vasp run.


    * **\*\*kwargs** – All kwargs supported by MPSOCSet, other than structure,
    prev_incar and prev_chgcar which are determined from the
    prev_calc_dir.



#### _property_ incar(_: [Incar](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Incar_ )
Incar


* **Type**

    return



#### override_from_prev_calc(prev_calc_dir='.')
Update the input set to include settings from a previous calculation.


* **Parameters**

    **prev_calc_dir** (*str*) – The path to the previous calculation directory.



* **Returns**

    The input set with the settings (structure, k-points, incar, etc)
    updated using the previous VASP run.



#### user_potcar_functional(_: UserPotcarFunctiona_ )

### _class_ pymatgen.io.vasp.sets.MPScanRelaxSet(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), bandgap=0, \*\*kwargs)
Bases: `DictSet`

Class for writing a relaxation input set using the accurate and numerically
efficient r2SCAN variant of the Strongly Constrained and Appropriately Normed
(SCAN) metaGGA density functional.

### Notes

1. This functional is officially supported in VASP 6.0.0 and above. On older version,
source code may be obtained by contacting the authors of the referenced manuscript.
The original SCAN functional, available from VASP 5.4.3 onwards, maybe used instead
by passing user_incar_settings={“METAGGA”: “SCAN”} when instantiating this InputSet.
r2SCAN and SCAN are expected to yield very similar results.

2. Meta-GGA calculations require POTCAR files that include
information on the kinetic energy density of the core-electrons,
i.e. “PBE_52” or “PBE_54”. Make sure the POTCARs include the
following lines (see VASP wiki for more details):

> $ grep kinetic POTCAR
> kinetic energy-density
> mkinetic energy-density pseudized
> kinetic energy density (partial)

### References

James W. Furness, Aaron D. Kaplan, Jinliang Ning, John P. Perdew, and Jianwei Sun.
Accurate and Numerically Efficient r2SCAN Meta-Generalized Gradient Approximation.
The Journal of Physical Chemistry Letters 0, 11 DOI: 10.1021/acs.jpclett.0c02405


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Input structure.


    * **bandgap** (*int*) – Bandgap of the structure in eV. The bandgap is used to
    compute the appropriate k-point density and determine the
    smearing settings.

    Metallic systems (default, bandgap = 0) use a KSPACING value of 0.22
    and Methfessel-Paxton order 2 smearing (ISMEAR=2, SIGMA=0.2).

    Non-metallic systems (bandgap > 0) use the tetrahedron smearing
    method (ISMEAR=-5, SIGMA=0.05). The KSPACING value is
    calculated from the bandgap via Eqs. 25 and 29 of Wisesa, McGill,
    and Mueller [1] (see References). Note that if ‘user_incar_settings’
    or ‘user_kpoints_settings’ override KSPACING, the calculation from
    bandgap is not performed.



    * **vdw** (*str*) – set “rVV10” to enable SCAN+rVV10, which is a versatile
    van der Waals density functional by combing the SCAN functional
    with the rVV10 non-local correlation functional. rvv10 is the only
    dispersion correction available for SCAN at this time.


    * **\*\*kwargs** – Same as those supported by DictSet.


### References

[1] P. Wisesa, K.A. McGill, T. Mueller, Efficient generation of
generalized Monkhorst-Pack grids through the use of informatics,
Phys. Rev. B. 93 (2016) 1-10. doi:10.1103/PhysRevB.93.155109.


#### CONFIG(_ = {'INCAR': {'ALGO': 'ALL', 'EDIFF': 1e-05, 'EDIFFG': -0.02, 'ENAUG': 1360, 'ENCUT': 680, 'IBRION': 2, 'ISIF': 3, 'ISPIN': 2, 'LAECHG': True, 'LASPH': True, 'LCHARG': True, 'LELF': True, 'LMIXTAU': True, 'LORBIT': 11, 'LREAL': 'Auto', 'LVTOT': True, 'LWAVE': False, 'MAGMOM': {'Ce': 5, 'Ce3+': 1, 'Co': 0.6, 'Co3+': 0.6, 'Co4+': 1, 'Cr': 5, 'Dy3+': 5, 'Er3+': 3, 'Eu': 10, 'Eu2+': 7, 'Eu3+': 6, 'Fe': 5, 'Gd3+': 7, 'Ho3+': 4, 'La3+': 0.6, 'Lu3+': 0.6, 'Mn': 5, 'Mn3+': 4, 'Mn4+': 3, 'Mo': 5, 'Nd3+': 3, 'Ni': 5, 'Pm3+': 4, 'Pr3+': 2, 'Sm3+': 5, 'Tb3+': 6, 'Tm3+': 2, 'V': 5, 'W': 5, 'Yb3+': 1}, 'METAGGA': 'R2SCAN', 'NELM': 200, 'NSW': 99, 'PREC': 'Accurate'}, 'PARENT': 'VASPIncarBase', 'POTCAR': {'Ac': 'Ac', 'Ag': 'Ag', 'Al': 'Al', 'Am': 'Am', 'Ar': 'Ar', 'As': 'As', 'At': 'At', 'Au': 'Au', 'B': 'B', 'Ba': 'Ba_sv', 'Be': 'Be_sv', 'Bi': 'Bi', 'Br': 'Br', 'C': 'C', 'Ca': 'Ca_sv', 'Cd': 'Cd', 'Ce': 'Ce', 'Cf': 'Cf', 'Cl': 'Cl', 'Cm': 'Cm', 'Co': 'Co', 'Cr': 'Cr_pv', 'Cs': 'Cs_sv', 'Cu': 'Cu_pv', 'Dy': 'Dy_3', 'Er': 'Er_3', 'Eu': 'Eu', 'F': 'F', 'Fe': 'Fe_pv', 'Fr': 'Fr_sv', 'Ga': 'Ga_d', 'Gd': 'Gd', 'Ge': 'Ge_d', 'H': 'H', 'He': 'He', 'Hf': 'Hf_pv', 'Hg': 'Hg', 'Ho': 'Ho_3', 'I': 'I', 'In': 'In_d', 'Ir': 'Ir', 'K': 'K_sv', 'Kr': 'Kr', 'La': 'La', 'Li': 'Li_sv', 'Lu': 'Lu_3', 'Mg': 'Mg_pv', 'Mn': 'Mn_pv', 'Mo': 'Mo_pv', 'N': 'N', 'Na': 'Na_pv', 'Nb': 'Nb_pv', 'Nd': 'Nd_3', 'Ne': 'Ne', 'Ni': 'Ni_pv', 'Np': 'Np', 'O': 'O', 'Os': 'Os_pv', 'P': 'P', 'Pa': 'Pa', 'Pb': 'Pb_d', 'Pd': 'Pd', 'Pm': 'Pm_3', 'Po': 'Po_d', 'Pr': 'Pr_3', 'Pt': 'Pt', 'Pu': 'Pu', 'Ra': 'Ra_sv', 'Rb': 'Rb_sv', 'Re': 'Re_pv', 'Rh': 'Rh_pv', 'Rn': 'Rn', 'Ru': 'Ru_pv', 'S': 'S', 'Sb': 'Sb', 'Sc': 'Sc_sv', 'Se': 'Se', 'Si': 'Si', 'Sm': 'Sm_3', 'Sn': 'Sn_d', 'Sr': 'Sr_sv', 'Ta': 'Ta_pv', 'Tb': 'Tb_3', 'Tc': 'Tc_pv', 'Te': 'Te', 'Th': 'Th', 'Ti': 'Ti_pv', 'Tl': 'Tl_d', 'Tm': 'Tm_3', 'U': 'U', 'V': 'V_pv', 'W': 'W_sv', 'Xe': 'Xe', 'Y': 'Y_sv', 'Yb': 'Yb_3', 'Zn': 'Zn', 'Zr': 'Zr_sv'}, 'POTCAR_FUNCTIONAL': 'PBE_54'_ )

### _class_ pymatgen.io.vasp.sets.MPScanStaticSet(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), bandgap=0, prev_incar=None, lepsilon=False, lcalcpol=False, \*\*kwargs)
Bases: `MPScanRelaxSet`

Creates input files for a static calculation using the accurate and numerically
efficient r2SCAN variant of the Strongly Constrained and Appropriately Normed
(SCAN) metaGGA functional.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Structure from previous run.


    * **bandgap** (*float*) – Bandgap of the structure in eV. The bandgap is used to
    compute the appropriate k-point density and determine the
    smearing settings.


    * **prev_incar** ([*Incar*](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Incar)) – Incar file from previous run.


    * **lepsilon** (*bool*) – Whether to add static dielectric calculation


    * **lcalcpol** (*bool*) – Whether to turn on evaluation of the Berry phase approximations
    for electronic polarization.


    * **\*\*kwargs** – kwargs supported by MPScanRelaxSet.



#### _classmethod_ from_prev_calc(prev_calc_dir, \*\*kwargs)
Generate a set of VASP input files for static calculations from a
directory of previous VASP run.


* **Parameters**


    * **prev_calc_dir** (*str*) – Directory containing the outputs(
    vasprun.xml and OUTCAR) of previous vasp run.


    * **\*\*kwargs** – All kwargs supported by MPScanStaticSet, other than prev_incar
    which is determined from the prev_calc_dir.



#### _property_ incar()
Incar


* **Type**

    return



#### override_from_prev_calc(prev_calc_dir='.')
Update the input set to include settings from a previous calculation.


* **Parameters**

    **prev_calc_dir** (*str*) – The path to the previous calculation directory.



* **Returns**

    The input set with the settings (structure, k-points, incar, etc)
    updated using the previous VASP run.



#### user_potcar_functional(_: UserPotcarFunctiona_ )

### _class_ pymatgen.io.vasp.sets.MPStaticSet(structure, prev_incar=None, prev_kpoints=None, lepsilon=False, lcalcpol=False, reciprocal_density=100, small_gap_multiply=None, \*\*kwargs)
Bases: `MPRelaxSet`

Creates input files for a static calculation.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Structure from previous run.


    * **prev_incar** ([*Incar*](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Incar)) – Incar file from previous run.


    * **prev_kpoints** ([*Kpoints*](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Kpoints)) – Kpoints from previous run.


    * **lepsilon** (*bool*) – Whether to add static dielectric calculation


    * **lcalcpol** (*bool*) – Whether to turn on evaluation of the Berry phase approximations
    for electronic polarization


    * **reciprocal_density** (*int*) – For static calculations, we usually set the
    reciprocal density by volume. This is a convenience arg to change
    that, rather than using user_kpoints_settings. Defaults to 100,
    which is ~50% more than that of standard relaxation calculations.


    * **small_gap_multiply** (*[**float**, **float**]*) – If the gap is less than
    1st index, multiply the default reciprocal_density by the 2nd
    index.


    * **\*\*kwargs** – kwargs supported by MPRelaxSet.



#### _classmethod_ from_prev_calc(prev_calc_dir, \*\*kwargs)
Generate a set of VASP input files for static calculations from a
directory of previous VASP run.


* **Parameters**


    * **prev_calc_dir** (*str*) – Directory containing the outputs(
    vasprun.xml and OUTCAR) of previous vasp run.


    * **\*\*kwargs** – All kwargs supported by MPStaticSet, other than prev_incar
    and prev_structure and prev_kpoints which are determined from
    the prev_calc_dir.



#### _property_ incar()
Incar


* **Type**

    return



#### _property_ kpoints(_: [Kpoints](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Kpoints) | Non_ )
Kpoints


* **Type**

    return



#### override_from_prev_calc(prev_calc_dir='.')
Update the input set to include settings from a previous calculation.


* **Parameters**

    **prev_calc_dir** (*str*) – The path to the previous calculation directory.



* **Returns**

    The input set with the settings (structure, k-points, incar, etc)
    updated using the previous VASP run.



#### user_potcar_functional(_: UserPotcarFunctiona_ )

### _class_ pymatgen.io.vasp.sets.MVLElasticSet(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), potim: float = 0.015, \*\*kwargs)
Bases: `MPRelaxSet`

MVL denotes VASP input sets that are implemented by the Materials Virtual
Lab ([http://materialsvirtuallab.org](http://materialsvirtuallab.org)) for various research.

This input set is used to calculate elastic constants in VASP. It is used
in the following work:

```default
Z. Deng, Z. Wang, I.-H. Chu, J. Luo, S. P. Ong.
“Elastic Properties of Alkali Superionic Conductor Electrolytes
from First Principles Calculations”, J. Electrochem. Soc.
2016, 163(2), A67-A74. doi: 10.1149/2.0061602jes
```

To read the elastic constants, you may use the Outcar class which parses the
elastic constants.


* **Parameters**


    * **structure** (*pymatgen.Structure*) – Input structure.


    * **potim** (*float*) – POTIM parameter. The default of 0.015 is usually fine,
    but some structures may require a smaller step.


    * **kwargs** – Parameters supported by MPRelaxSet.



#### user_potcar_functional(_: UserPotcarFunctiona_ )

### _class_ pymatgen.io.vasp.sets.MVLGBSet(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), k_product=40, slab_mode=False, is_metal=True, \*\*kwargs)
Bases: `MPRelaxSet`

Class for writing a vasp input files for grain boundary calculations, slab
or bulk.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – provide the structure


    * **k_product** – Kpoint number \* length for a & b directions, also for c
    direction in bulk calculations. Default to 40.


    * **slab_mode** (*bool*) – Defaults to False. Use default (False) for a
    bulk supercell. Use True if you are performing calculations on a
    slab-like (i.e., surface) of the GB, for example, when you are
    calculating the work of separation.


    * **is_metal** (*bool*) – Defaults to True. This determines whether an ISMEAR of
    1 is used (for metals) or not (for insulators and semiconductors)
    by default. Note that it does *not* override user_incar_settings,
    which can be set by the user to be anything desired.


    * **\*\*kwargs** – Other kwargs supported by `MPRelaxSet`.



#### _property_ incar()
Incar


* **Type**

    return



#### _property_ kpoints()
k_product, default to 40, is kpoint number \* length for a & b
directions, also for c direction in bulk calculations
Automatic mesh & Gamma is the default setting.


#### user_potcar_functional(_: UserPotcarFunctiona_ )

### _class_ pymatgen.io.vasp.sets.MVLGWSet(structure, prev_incar=None, nbands=None, reciprocal_density=100, mode='STATIC', copy_wavecar=True, nbands_factor=5, ncores=16, \*\*kwargs)
Bases: `DictSet`

MVL denotes VASP input sets that are implemented by the Materials Virtual
Lab ([http://materialsvirtuallab.org](http://materialsvirtuallab.org)) for various research. This is a
flexible input set for GW calculations.

Note that unlike all other input sets in this module, the PBE_54 series of
functional is set as the default. These have much improved performance for
GW calculations.

A typical sequence is mode=”STATIC” -> mode=”DIAG” -> mode=”GW” ->
mode=”BSE”. For all steps other than the first one (static), the
recommendation is to use from_prev_calculation on the preceding run in
the series.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Input structure.


    * **prev_incar** (*Incar/string*) – Incar file from previous run.


    * **mode** (*str*) – Supported modes are “STATIC” (default), “DIAG”, “GW”,
    and “BSE”.


    * **nbands** (*int*) – For subsequent calculations, it is generally
    recommended to perform NBANDS convergence starting from the
    NBANDS of the previous run for DIAG, and to use the exact same
    NBANDS for GW and BSE. This parameter is used by
    from_previous_calculation to set nband.


    * **copy_wavecar** – Whether to copy the old WAVECAR, WAVEDER and associated
    files when starting from a previous calculation.


    * **nbands_factor** (*int*) – Multiplicative factor for NBANDS when starting
    from a previous calculation. Only applies if mode==”DIAG”.
    Need to be tested for convergence.


    * **reciprocal_density** (*int*) – Density of k-mesh by reciprocal atom. Only
    applies if mode==”STATIC”. Defaults to 100.


    * **ncores** (*int*) – Numbers of cores used for the calculation. VASP will alter
    NBANDS if it was not dividable by ncores. Only applies if
    mode==”DIAG”.


    * **\*\*kwargs** – All kwargs supported by DictSet. Typically,
    user_incar_settings is a commonly used option.



#### CONFIG(_ = {'INCAR': {'ALGO': 'Normal', 'EDIFF': 1e-08, 'IBRION': -1, 'ICHARG': 1, 'ISMEAR': 0, 'ISPIN': 2, 'LORBIT': 11, 'LREAL': 'AUTO', 'LWAVE': True, 'MAGMOM': {'Ce': 5, 'Ce3+': 1, 'Co': 0.6, 'Co3+': 0.6, 'Co4+': 1, 'Cr': 5, 'Dy3+': 5, 'Er3+': 3, 'Eu': 10, 'Eu2+': 7, 'Eu3+': 6, 'Fe': 5, 'Gd3+': 7, 'Ho3+': 4, 'La3+': 0.6, 'Lu3+': 0.6, 'Mn': 5, 'Mn3+': 4, 'Mn4+': 3, 'Mo': 5, 'Nd3+': 3, 'Ni': 5, 'Pm3+': 4, 'Pr3+': 2, 'Sm3+': 5, 'Tb3+': 6, 'Tm3+': 2, 'V': 5, 'W': 5, 'Yb3+': 1}, 'NELM': 100, 'PREC': 'Accurate', 'SIGMA': 0.01}, 'KPOINTS': {'reciprocal_density': 100}, 'PARENT': 'VASPIncarBase', 'POTCAR': {'Ac': 'Ac', 'Ag': 'Ag_sv_GW', 'Al': 'Al_GW', 'Ar': 'Ar_GW', 'As': 'As_GW', 'At': 'At_d_GW', 'Au': 'Au_sv_GW', 'B': 'B_GW', 'Ba': 'Ba_sv_GW', 'Be': 'Be_sv_GW', 'Bi': 'Bi_d_GW', 'Br': 'Br_GW', 'C': 'C_GW', 'Ca': 'Ca_sv_GW', 'Cd': 'Cd_sv_GW', 'Ce': 'Ce_GW', 'Cl': 'Cl_GW', 'Co': 'Co_sv_GW', 'Cr': 'Cr_sv_GW', 'Cs': 'Cs_sv_GW', 'Cu': 'Cu_sv_GW', 'Dy': 'Dy_3', 'Er': 'Er_3', 'Eu': 'Eu', 'F': 'F_GW', 'Fe': 'Fe_sv_GW', 'Ga': 'Ga_d_GW', 'Gd': 'Gd', 'Ge': 'Ge_d_GW', 'H': 'H_GW', 'He': 'He_GW', 'Hf': 'Hf_sv_GW', 'Hg': 'Hg_sv_GW', 'Ho': 'Ho_3', 'I': 'I_GW', 'In': 'In_d_GW', 'Ir': 'Ir_sv_GW', 'K': 'K_sv_GW', 'Kr': 'Kr_GW', 'La': 'La_GW', 'Li': 'Li_sv_GW', 'Lu': 'Lu_3', 'Mg': 'Mg_sv_GW', 'Mn': 'Mn_sv_GW', 'Mo': 'Mo_sv_GW', 'N': 'N_GW', 'Na': 'Na_sv_GW', 'Nb': 'Nb_sv_GW', 'Nd': 'Nd_3', 'Ne': 'Ne_GW', 'Ni': 'Ni_sv_GW', 'Np': 'Np', 'O': 'O_GW', 'Os': 'Os_sv_GW', 'P': 'P_GW', 'Pa': 'Pa', 'Pb': 'Pb_d_GW', 'Pd': 'Pd_sv_GW', 'Pm': 'Pm_3', 'Po': 'Po_d_GW', 'Pr': 'Pr_3', 'Pt': 'Pt_sv_GW', 'Pu': 'Pu', 'Rb': 'Rb_sv_GW', 'Re': 'Re_sv_GW', 'Rh': 'Rh_sv_GW', 'Rn': 'Rn_d_GW', 'Ru': 'Ru_sv_GW', 'S': 'S_GW', 'Sb': 'Sb_d_GW', 'Sc': 'Sc_sv_GW', 'Se': 'Se_GW', 'Si': 'Si_GW', 'Sm': 'Sm_3', 'Sn': 'Sn_d_GW', 'Sr': 'Sr_sv_GW', 'Ta': 'Ta_sv_GW', 'Tb': 'Tb_3', 'Tc': 'Tc_sv_GW', 'Te': 'Te_GW', 'Th': 'Th', 'Ti': 'Ti_sv_GW', 'Tl': 'Tl_d_GW', 'Tm': 'Tm_3', 'U': 'U', 'V': 'V_sv_GW', 'W': 'W_sv_GW', 'Xe': 'Xe_GW', 'Y': 'Y_sv_GW', 'Yb': 'Yb_3', 'Zn': 'Zn_sv_GW', 'Zr': 'Zr_sv_GW'}, 'POTCAR_FUNCTIONAL': 'PBE_54'_ )

#### SUPPORTED_MODES(_ = ('DIAG', 'GW', 'STATIC', 'BSE'_ )

#### _classmethod_ from_prev_calc(prev_calc_dir, mode='DIAG', \*\*kwargs)
Generate a set of VASP input files for GW or BSE calculations from a
directory of previous Exact Diag VASP run.


* **Parameters**


    * **prev_calc_dir** (*str*) – The directory contains the outputs(
    vasprun.xml of previous vasp run.


    * **mode** (*str*) – Supported modes are “STATIC”, “DIAG” (default), “GW”,
    and “BSE”.


    * **\*\*kwargs** – All kwargs supported by MVLGWSet, other than structure,
    prev_incar and mode, which are determined from the
    prev_calc_dir.



#### _property_ incar()
Incar


* **Type**

    return



#### _property_ kpoints()
Generate gamma center k-points mesh grid for GW calc,
which is requested by GW calculation.


#### override_from_prev_calc(prev_calc_dir='.')
Update the input set to include settings from a previous calculation.


* **Parameters**

    **prev_calc_dir** (*str*) – The path to the previous calculation directory.



* **Returns**

    The input set with the settings (structure, k-points, incar, etc)
    updated using the previous VASP run.



### _class_ pymatgen.io.vasp.sets.MVLNPTMDSet(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), start_temp, end_temp, nsteps, time_step=2, spin_polarized=False, \*\*kwargs)
Bases: `MITMDSet`

Class for writing a vasp md run in NPT ensemble.

### Notes

To eliminate Pulay stress, the default ENCUT is set to a rather large
value of ENCUT, which is 1.5 \* ENMAX.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – input structure.


    * **start_temp** (*int*) – Starting temperature.


    * **end_temp** (*int*) – Final temperature.


    * **nsteps** (*int*) – Number of time steps for simulations. NSW parameter.


    * **time_step** (*int*) – The time step for the simulation. The POTIM
    parameter. Defaults to 2fs.


    * **spin_polarized** (*bool*) – Whether to do spin polarized calculations.
    The ISPIN parameter. Defaults to False.


    * **\*\*kwargs** – Other kwargs supported by `DictSet`.



#### user_potcar_functional(_: UserPotcarFunctiona_ )

### _class_ pymatgen.io.vasp.sets.MVLRelax52Set(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), \*\*kwargs)
Bases: `DictSet`

Implementation of VaspInputSet utilizing the public Materials Project
parameters for INCAR & KPOINTS and VASP’s recommended PAW potentials for
POTCAR.

Keynotes from VASP manual:


    1. Recommended potentials for calculations using vasp.5.2+


    2. If dimers with short bonds are present in the compound (O2, CO,

        N2, F2, P2, S2, Cl2), it is recommended to use the h potentials.
        Specifically, C_h, O_h, N_h, F_h, P_h, S_h, Cl_h


    3. Released on Oct 28, 2018 by VASP. Please refer to VASP

        Manual 1.2, 1.3 & 10.2.1 for more details.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – input structure.


    * **user_potcar_functional** (*str*) – choose from “PBE_52” and “PBE_54”.


    * **\*\*kwargs** – Other kwargs supported by `DictSet`.



#### CONFIG(_ = {'INCAR': {'ALGO': 'FAST', 'EDIFF_PER_ATOM': 5e-05, 'ENCUT': 520, 'IBRION': 2, 'ICHARG': 1, 'ISIF': 3, 'ISMEAR': -5, 'ISPIN': 2, 'LDAU': True, 'LDAUJ': {'F': {'Co': 0, 'Cr': 0, 'Fe': 0, 'Mn': 0, 'Mo': 0, 'Ni': 0, 'V': 0, 'W': 0}, 'O': {'Co': 0, 'Cr': 0, 'Fe': 0, 'Mn': 0, 'Mo': 0, 'Ni': 0, 'V': 0, 'W': 0}}, 'LDAUL': {'F': {'Co': 2, 'Cr': 2, 'Fe': 2, 'Mn': 2, 'Mo': 2, 'Ni': 2, 'V': 2, 'W': 2}, 'O': {'Co': 2, 'Cr': 2, 'Fe': 2, 'Mn': 2, 'Mo': 2, 'Ni': 2, 'V': 2, 'W': 2}}, 'LDAUPRINT': 1, 'LDAUTYPE': 2, 'LDAUU': {'F': {'Co': 3.32, 'Cr': 3.7, 'Fe': 5.3, 'Mn': 3.9, 'Mo': 4.38, 'Ni': 6.2, 'V': 3.25, 'W': 6.2}, 'O': {'Co': 3.32, 'Cr': 3.7, 'Fe': 5.3, 'Mn': 3.9, 'Mo': 4.38, 'Ni': 6.2, 'V': 3.25, 'W': 6.2}}, 'LORBIT': 11, 'LREAL': 'AUTO', 'LWAVE': False, 'NELM': 100, 'NSW': 99, 'PREC': 'Accurate', 'SIGMA': 0.05}, 'KPOINTS': {'reciprocal_density': 64}, 'POTCAR': {'Ac': 'Ac', 'Ag': 'Ag', 'Al': 'Al', 'Am': 'Am', 'Ar': 'Ar', 'As': 'As', 'At': 'At_d', 'Au': 'Au', 'B': 'B', 'Ba': 'Ba_sv', 'Be': 'Be', 'Bi': 'Bi_d', 'Br': 'Br', 'C': 'C', 'Ca': 'Ca_sv', 'Cd': 'Cd', 'Ce': 'Ce', 'Cl': 'Cl', 'Cm': 'Cm', 'Co': 'Co', 'Cr': 'Cr_pv', 'Cs': 'Cs_sv', 'Cu': 'Cu', 'Dy': 'Dy_3', 'Er': 'Er_3', 'Eu': 'Eu_2', 'F': 'F', 'Fe': 'Fe', 'Fr': 'Fr_sv', 'Ga': 'Ga_d', 'Gd': 'Gd_3', 'Ge': 'Ge_d', 'H': 'H', 'He': 'He', 'Hf': 'Hf_pv', 'Hg': 'Hg', 'Ho': 'Ho_3', 'I': 'I', 'In': 'In_d', 'Ir': 'Ir', 'K': 'K_sv', 'Kr': 'Kr', 'La': 'La', 'Li': 'Li_sv', 'Lu': 'Lu_3', 'Mg': 'Mg', 'Mn': 'Mn_pv', 'Mo': 'Mo_sv', 'N': 'N', 'Na': 'Na_pv', 'Nb': 'Nb_sv', 'Nd': 'Nd_3', 'Ne': 'Ne', 'Ni': 'Ni', 'Np': 'Np', 'O': 'O', 'Os': 'Os', 'P': 'P', 'Pa': 'Pa', 'Pb': 'Pb_d', 'Pd': 'Pd', 'Pm': 'Pm_3', 'Po': 'Po_d', 'Pr': 'Pr_3', 'Pt': 'Pt', 'Pu': 'Pu', 'Ra': 'Ra_sv', 'Rb': 'Rb_sv', 'Re': 'Re', 'Rh': 'Rh_pv', 'Rn': 'Rn', 'Ru': 'Ru_pv', 'S': 'S', 'Sb': 'Sb', 'Sc': 'Sc_sv', 'Se': 'Se', 'Si': 'Si', 'Sm': 'Sm_3', 'Sn': 'Sn_d', 'Sr': 'Sr_sv', 'Ta': 'Ta_pv', 'Tb': 'Tb_3', 'Tc': 'Tc_pv', 'Te': 'Te', 'Th': 'Th', 'Ti': 'Ti_sv', 'Tl': 'Tl_d', 'Tm': 'Tm_3', 'U': 'U', 'V': 'V_sv', 'W': 'W_pv', 'Xe': 'Xe', 'Y': 'Y_sv', 'Yb': 'Yb_2', 'Zn': 'Zn', 'Zr': 'Zr_sv'}, 'POTCAR_FUNCTIONAL': 'PBE_52'_ )

### _class_ pymatgen.io.vasp.sets.MVLScanRelaxSet(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), \*\*kwargs)
Bases: `MPRelaxSet`

Class for writing a relax input set using Strongly Constrained and
Appropriately Normed (SCAN) semilocal density functional.

### Notes


1. This functional is only available from VASP.5.4.3 upwards.

2. Meta-GGA calculations require POTCAR files that include
information on the kinetic energy density of the core-electrons,
i.e. “PBE_52” or “PBE_54”. Make sure the POTCAR including the
following lines (see VASP wiki for more details):

> $ grep kinetic POTCAR
> kinetic energy-density
> mkinetic energy-density pseudized
> kinetic energy density (partial)


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – input structure.


    * **vdw** (*str*) – set “rVV10” to enable SCAN+rVV10, which is a versatile
    van der Waals density functional by combing the SCAN functional
    with the rVV10 non-local correlation functional.


    * **\*\*kwargs** – Other kwargs supported by `DictSet`.



#### user_potcar_functional(_: UserPotcarFunctiona_ )

### _class_ pymatgen.io.vasp.sets.MVLSlabSet(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), k_product=50, bulk=False, auto_dipole=False, set_mix=True, sort_structure=True, \*\*kwargs)
Bases: `MPRelaxSet`

Class for writing a set of slab vasp runs,
including both slabs (along the c direction) and orient unit cells (bulk),
to ensure the same KPOINTS, POTCAR and INCAR criterion.


* **Parameters**


    * **structure** – Structure


    * **k_product** – default to 50, kpoint number \* length for a & b
    directions, also for c direction in bulk calculations


    * **bulk** –


    * **auto_dipole** –


    * **set_mix** –


    * **sort_structure** –


    * **kwargs** – Other kwargs supported by `DictSet`.



#### as_dict(verbosity=2)

* **Parameters**

    **verbosity** – Verbosity of dict. E.g., whether to include Structure.



* **Returns**

    MSONable dict



#### _property_ kpoints()
k_product, default to 50, is kpoint number \* length for a & b

    directions, also for c direction in bulk calculations

Automatic mesh & Gamma is the default setting.


#### user_potcar_functional(_: UserPotcarFunctiona_ )

### _class_ pymatgen.io.vasp.sets.VaspInputSet()
Bases: `MSONable`

Base class representing a set of VASP input parameters with a structure
supplied as init parameters. Typically, you should not inherit from this
class. Start from DictSet or MPRelaxSet or MITRelaxSet.


#### as_dict(verbosity=2)

* **Parameters**


    * **verbosity** – Verbosity for generated dict. If 1, structure is


    * **excluded.** –



* **Returns**

    MSONable dict



#### get_vasp_input()

* **Returns**

    VaspInput.



#### _abstract property_ incar()
Incar object.


#### _abstract property_ kpoints()
Kpoints object.


#### _abstract property_ poscar()
Poscar object.


#### _property_ potcar(_: [Potcar](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Potcar_ )
Potcar object.


#### _property_ potcar_symbols()
List of POTCAR symbols.


#### write_input(output_dir: str, make_dir_if_not_present: bool = True, include_cif: bool = False, potcar_spec: bool = False, zip_output: bool = False)
Writes a set of VASP input to a directory.


* **Parameters**


    * **output_dir** (*str*) – Directory to output the VASP input files


    * **make_dir_if_not_present** (*bool*) – Set to True if you want the
    directory (and the whole path) to be created if it is not
    present.


    * **include_cif** (*bool*) – Whether to write a CIF file in the output
    directory for easier opening by VESTA.


    * **potcar_spec** (*bool*) – Instead of writing the POTCAR, write a “POTCAR.spec”.
    This is intended to help sharing an input set with people who might
    not have a license to specific Potcar files. Given a “POTCAR.spec”,
    the specific POTCAR file can be re-generated using pymatgen with the
    “generate_potcar” function in the pymatgen CLI.


    * **zip_output** (*bool*) – If True, output will be zipped into a file with the
    same name as the InputSet (e.g., MPStaticSet.zip)



### pymatgen.io.vasp.sets.batch_write_input(structures, vasp_input_set=<class 'pymatgen.io.vasp.sets.MPRelaxSet'>, output_dir='.', make_dir_if_not_present=True, subfolder=None, sanitize=False, include_cif=False, potcar_spec=False, zip_output=False, \*\*kwargs)
Batch write vasp input for a sequence of structures to
output_dir, following the format output_dir/{group}/{formula}_{number}.


* **Parameters**


    * **structures** (*[*[*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)*]*) – Sequence of Structures.


    * **vasp_input_set** (*VaspInputSet*) – VaspInputSet class that creates
    vasp input files from structures. Note that a class should be
    supplied. Defaults to MPRelaxSet.


    * **output_dir** (*str*) – Directory to output files. Defaults to current
    directory “.”.


    * **make_dir_if_not_present** (*bool*) – Create the directory if not present.
    Defaults to True.


    * **subfolder** (*callable*) – Function to create subdirectory name from
    structure. Defaults to simply “formula_count”.


    * **sanitize** (*bool*) – Boolean indicating whether to sanitize the
    structure before writing the VASP input files. Sanitized output
    are generally easier for viewing and certain forms of analysis.
    Defaults to False.


    * **include_cif** (*bool*) – Whether to output a CIF as well. CIF files are
    generally better supported in visualization programs.


    * **potcar_spec** (*bool*) – Instead of writing the POTCAR, write a “POTCAR.spec”.
    This is intended to help sharing an input set with people who might
    not have a license to specific Potcar files. Given a “POTCAR.spec”,
    the specific POTCAR file can be re-generated using pymatgen with the
    “generate_potcar” function in the pymatgen CLI.


    * **zip_output** (*bool*) – If True, output will be zipped into a file with the
    same name as the InputSet (e.g., MPStaticSet.zip)


    * **\*\*kwargs** – Additional kwargs are passed to the vasp_input_set class
    in addition to structure.



### pymatgen.io.vasp.sets.get_structure_from_prev_run(vasprun, outcar=None)
Process structure from previous run.


* **Parameters**


    * **vasprun** ([*Vasprun*](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Vasprun)) – Vasprun that contains the final structure
    from previous run.


    * **outcar** ([*Outcar*](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Outcar)) – Outcar that contains the magnetization info from
    previous run.



* **Returns**

    Returns the magmom-decorated structure that can be passed to get
    VASP input files, e.g. get_kpoints.



### pymatgen.io.vasp.sets.get_valid_magmom_struct(structure, inplace=True, spin_mode='auto')
Make sure that the structure has valid magmoms based on the kind of calculation
Fill in missing Magmom values.


* **Parameters**


    * **structure** – The input structure


    * **inplace** – True - edit the magmom of the input structurel; False - return new structure


    * **spin_mode** – “scalar”/”vector”/”none”/”auto” only first letter (s/v/n) is needed.
    dictates how the spin configuration will be determined.


        * auto: read the existing magmom values and decide


        * scalar: use a single scalar value (for spin up/down)


        * vector: use a vector value for spin-orbit systems


        * none: Remove all the magmom information




* **Returns**

    New structure if inplace is False



### pymatgen.io.vasp.sets.get_vasprun_outcar(path, parse_dos=True, parse_eigen=True)

* **Parameters**


    * **path** – Path to get the vasprun.xml and OUTCAR.


    * **parse_dos** – Whether to parse dos. Defaults to True.


    * **parse_eigen** – Whether to parse eigenvalue. Defaults to True.



* **Returns**




### pymatgen.io.vasp.sets.next_num_with_prime_factors(n: int, max_prime_factor: int, must_inc_2: bool = True)
Return the next number greater than or equal to n that only has the desired prime factors.


* **Parameters**


    * **n** (*int*) – Initial guess at the grid density


    * **max_prime_factor** (*int*) – the maximum prime factor


    * **must_inc_2** (*bool*) – 2 must be a prime factor of the result



* **Returns**

    first product of the prime_factors that is >= n



* **Return type**

    int



### pymatgen.io.vasp.sets.primes_less_than(max_val: int)
Get the primes less than or equal to the max value.


### pymatgen.io.vasp.sets.standardize_structure(structure, sym_prec=0.1, international_monoclinic=True)
Get the symmetrically standardized structure.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – The structure.


    * **sym_prec** (*float*) – Tolerance for symmetry finding for standardization.


    * **international_monoclinic** (*bool*) – Whether to use international
    convention (vs Curtarolo) for monoclinic. Defaults True.



* **Returns**

    The symmetrized structure.