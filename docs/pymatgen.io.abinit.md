---
layout: default
title: pymatgen.io.abinit.md
nav_exclude: true
---

1. TOC
{:toc}

# pymatgen.io.abinit package

This package implements basic input and output capabilities for Abinit.


## pymatgen.io.abinit.abiobjects module

Low-level objects providing an abstraction for the objects involved in the calculation.


### _class_ AbivarAble()
Bases: `object`

An AbivarAble object provides a method to_abivars
that returns a dictionary with the abinit variables.


#### _abc_impl(_ = <_abc._abc_data object_ )

#### _abstract_ to_abivars()
Returns a dictionary with the abinit variables.


### _class_ Constraints()
Bases: `AbivarAble`

This object defines the constraints for structural relaxation.


#### _abc_impl(_ = <_abc._abc_data object_ )

#### to_abivars()
Dictionary with Abinit variables.


### _class_ Electrons(spin_mode='polarized', smearing='fermi_dirac:0.1 eV', algorithm=None, nband=None, fband=None, charge=0.0, comment=None)
Bases: `AbivarAble`, `MSONable`

The electronic degrees of freedom.

Constructor for Electrons object.


* **Parameters**


    * **comment** – String comment for Electrons


    * **charge** – Total charge of the system. Default is 0.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### as_dict()
Json friendly dict representation.


#### _classmethod_ from_dict(d)
Build object from dictionary.


#### _property_ nspden()
Number of independent density components.


#### _property_ nspinor()
Number of independent spinor components.


#### _property_ nsppol()
Number of independent spin polarizations.


#### to_abivars()
Return dictionary with Abinit variables.


### _class_ ElectronsAlgorithm(\*args, \*\*kwargs)
Bases: `dict`, `AbivarAble`, `MSONable`

Variables controlling the SCF/NSCF algorithm.

Initialize object.


#### _DEFAULT(_ = {'diecut': None, 'diegap': None, 'dielam': None, 'dielng': None, 'diemac': None, 'diemix': None, 'diemixmag': None, 'iprcell': None, 'iscf': None, 'nstep': 50_ )

#### _abc_impl(_ = <_abc._abc_data object_ )

#### as_dict()
Convert object to dict.


#### _classmethod_ from_dict(d)
Build object from dict.


#### to_abivars()
Dictionary with Abinit input variables.


### _class_ ExcHamiltonian(bs_loband, nband, mbpt_sciss, coulomb_mode, ecuteps, spin_mode='polarized', mdf_epsinf=None, exc_type='TDA', algo='haydock', with_lf=True, bs_freq_mesh=None, zcut=None, \*\*kwargs)
Bases: `AbivarAble`

This object contains the parameters for the solution of the Bethe-Salpeter equation.


* **Parameters**


    * **bs_loband** – Lowest band index (Fortran convention) used in the e-h  basis set.
    Can be scalar or array of shape (nsppol,). Must be >= 1 and <= nband


    * **nband** – Max band index used in the e-h  basis set.


    * **mbpt_sciss** – Scissors energy in Hartree.


    * **coulomb_mode** – Treatment of the Coulomb term.


    * **ecuteps** – Cutoff energy for W in Hartree.


    * **mdf_epsinf** – Macroscopic dielectric function $\\epsilon_\\inf$ used in
    the model dielectric function.


    * **exc_type** – Approximation used for the BSE Hamiltonian


    * **with_lf** – True if local field effects are included <==> exchange term is included


    * **bs_freq_mesh** – Frequency mesh for the macroscopic dielectric function (start, stop, step) in Ha.


    * **zcut** – Broadening parameter in Ha.


    * **\*\*kwargs** – Extra keywords.



#### _ALGO2VAR(_ = {'cg': 3, 'direct_diago': 1, 'haydock': 2_ )

#### _COULOMB_MODES(_ = ('diago', 'full', 'model_df'_ )

#### _EXC_TYPES(_ = {'TDA': 0, 'coupling': 1_ )

#### _abc_impl(_ = <_abc._abc_data object_ )

#### _property_ inclvkb()
Treatment of the dipole matrix element (NC pseudos, default is 2).


#### to_abivars()
Returns a dictionary with the abinit variables.


#### _property_ use_cg()
True if we are using the conjugate gradient method.


#### _property_ use_direct_diago()
True if we are performing the direct diagonalization of the BSE Hamiltonian.


#### _property_ use_haydock()
True if we are using the Haydock iterative technique.


### _class_ HilbertTransform(nomegasf, domegasf=None, spmeth=1, nfreqre=None, freqremax=None, nfreqim=None, freqremin=None)
Bases: `AbivarAble`

Parameters for the Hilbert-transform method (Screening code)
i.e. the parameters defining the frequency mesh used for the spectral function
and the frequency mesh used for the polarizability.


* **Parameters**


    * **nomegasf** – Number of points for sampling the spectral function along the real axis.


    * **domegasf** – Step in Ha for the linear mesh used for the spectral function.


    * **spmeth** – Algorithm for the representation of the delta function.


    * **nfreqre** – Number of points along the real axis (linear mesh).


    * **freqremax** – Maximum frequency for W along the real axis (in hartree).


    * **nfreqim** – Number of point along the imaginary axis (Gauss-Legendre mesh).


    * **freqremin** – Minimum frequency for W along the real axis (in hartree).



#### _abc_impl(_ = <_abc._abc_data object_ )

#### to_abivars()
Returns a dictionary with the abinit variables.


### _class_ KSampling(mode=KSamplingModes.monkhorst, num_kpts=0, kpts=((1, 1, 1),), kpt_shifts=(0.5, 0.5, 0.5), kpts_weights=None, use_symmetries=True, use_time_reversal=True, chksymbreak=None, comment=None)
Bases: `AbivarAble`, `MSONable`

Input variables defining the K-point sampling.

Highly flexible constructor for KSampling objects. The flexibility comes
at the cost of usability and in general, it is recommended that you use
the default constructor only if you know exactly what you are doing and
requires the flexibility. For most usage cases, the object be constructed
far more easily using the convenience static constructors:

>
> 1. gamma_only


> 2. gamma_centered


> 3. monkhorst


> 4. monkhorst_automatic


> 5. path

and it is recommended that you use those.


* **Parameters**


    * **mode** – Mode for generating k-poits. Use one of the KSamplingModes enum types.


    * **num_kpts** – Number of kpoints if mode is “automatic”
    Number of division for the sampling of the smallest segment if mode is “path”.
    Not used for the other modes


    * **kpts** – Number of divisions. Even when only a single specification is
    required, e.g. in the automatic scheme, the kpts should still
    be specified as a 2D array. e.g., [[20]] or [[2,2,2]].


    * **kpt_shifts** – Shifts for Kpoints.


    * **use_symmetries** – False if spatial symmetries should not be used
    to reduce the number of independent k-points.


    * **use_time_reversal** – False if time-reversal symmetry should not be used
    to reduce the number of independent k-points.


    * **kpts_weights** – Optional weights for kpoints. For explicit kpoints.


    * **chksymbreak** – Abinit input variable: check whether the BZ sampling preserves the symmetry of the crystal.


    * **comment** – String comment for Kpoints


**NOTE**: The default behavior of the constructor is monkhorst.


#### _abc_impl(_ = <_abc._abc_data object_ )

#### _classmethod_ _path(ndivsm, structure=None, kpath_bounds=None, comment=None)
Static constructor for path in k-space.


* **Parameters**


    * **structure** – Structure object.


    * **kpath_bounds** – List with the reduced coordinates of the k-points defining the path.


    * **ndivsm** – Number of division for the smallest segment.


    * **comment** – Comment string.



* **Returns**

    KSampling object.



#### as_dict()
Convert object to dict.


#### _classmethod_ automatic_density(structure, kppa, chksymbreak=None, use_symmetries=True, use_time_reversal=True, shifts=(0.5, 0.5, 0.5))
Returns an automatic Kpoint object based on a structure and a kpoint
density. Uses Gamma centered meshes for hexagonal cells and Monkhorst-Pack grids otherwise.

Algorithm:

    Uses a simple approach scaling the number of divisions along each
    reciprocal lattice vector proportional to its length.


* **Parameters**


    * **structure** – Input structure


    * **kppa** – Grid density



#### _classmethod_ explicit_path(ndivsm, kpath_bounds)
See _path for the meaning of the variables.


#### _classmethod_ from_dict(d)
Build object from dict.


#### _classmethod_ gamma_centered(kpts=(1, 1, 1), use_symmetries=True, use_time_reversal=True)
Convenient static constructor for an automatic Gamma centered Kpoint grid.


* **Parameters**


    * **kpts** – Subdivisions N_1, N_2 and N_3 along reciprocal lattice vectors.


    * **use_symmetries** – False if spatial symmetries should not be used
    to reduce the number of independent k-points.


    * **use_time_reversal** – False if time-reversal symmetry should not be used
    to reduce the number of independent k-points.



* **Returns**

    KSampling object.



#### _classmethod_ gamma_only()
Gamma-only sampling.


#### _property_ is_homogeneous(_: boo_ )
Homogeneous sampling.


#### _classmethod_ monkhorst(ngkpt, shiftk=(0.5, 0.5, 0.5), chksymbreak=None, use_symmetries=True, use_time_reversal=True, comment=None)
Convenient static constructor for a Monkhorst-Pack mesh.


* **Parameters**


    * **ngkpt** – Subdivisions N_1, N_2 and N_3 along reciprocal lattice vectors.


    * **shiftk** – Shift to be applied to the kpoints.


    * **use_symmetries** – Use spatial symmetries to reduce the number of k-points.


    * **use_time_reversal** – Use time-reversal symmetry to reduce the number of k-points.



* **Returns**

    KSampling object.



#### _classmethod_ monkhorst_automatic(structure, ngkpt, use_symmetries=True, use_time_reversal=True, chksymbreak=None, comment=None)
Convenient static constructor for an automatic Monkhorst-Pack mesh.


* **Parameters**


    * **structure** – Structure object.


    * **ngkpt** – Subdivisions N_1, N_2 and N_3 along reciprocal lattice vectors.


    * **use_symmetries** – Use spatial symmetries to reduce the number of k-points.


    * **use_time_reversal** – Use time-reversal symmetry to reduce the number of k-points.



* **Returns**

    KSampling object.



#### _classmethod_ path_from_structure(ndivsm, structure)
See _path for the meaning of the variables.


#### to_abivars()
Dictionary with Abinit variables.


### _class_ KSamplingModes(value)
Bases: `Enum`

Enum if the different samplings of the BZ.


#### automatic(_ = _ )

#### monkhorst(_ = _ )

#### path(_ = _ )

### _class_ ModelDielectricFunction(mdf_epsinf)
Bases: `AbivarAble`

Model dielectric function used for BSE calculation.


* **Parameters**

    **mdf_epsinf** – Value of epsilon_infinity.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### to_abivars()
Return dictionary with abinit variables.


### _class_ PPModel(mode='godby', plasmon_freq=None)
Bases: `AbivarAble`, `MSONable`

Parameters defining the plasmon-pole technique.
The common way to instantiate a PPModel object is via the class method PPModel.as_ppmodel(string).


* **Parameters**


    * **mode** – ppmodel type


    * **plasmon_freq** – Plasmon frequency in Ha.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### as_dict()
Convert object to dictionary.


#### _classmethod_ as_ppmodel(obj)
Constructs an instance of PPModel from obj.

Accepts obj in the form:


    * PPmodel instance


    * string. e.g “godby:12.3 eV”, “linden”.


#### _static_ from_dict(d)
Build object from dictionary.


#### _classmethod_ get_noppmodel()
Calculation without plasmon-pole model.


#### to_abivars()
Return dictionary with Abinit variables.


### _class_ PPModelModes(value)
Bases: `Enum`

Different kind of plasmon-pole models.


#### farid(_ = _ )

#### godby(_ = _ )

#### hybersten(_ = _ )

#### linden(_ = _ )

#### noppmodel(_ = _ )

### _class_ RelaxationMethod(\*args, \*\*kwargs)
Bases: `AbivarAble`, `MSONable`

This object stores the variables for the (constrained) structural optimization
ionmov and optcell specify the type of relaxation.
The other variables are optional and their use depend on ionmov and optcell.
A None value indicates that we use abinit default. Default values can
be modified by passing them to the constructor.
The set of variables are constructed in to_abivars depending on ionmov and optcell.

Initialize object.


#### IONMOV_DEFAULT(_ = _ )

#### OPTCELL_DEFAULT(_ = _ )

#### _abc_impl(_ = <_abc._abc_data object_ )

#### _default_vars(_ = {'atoms_constraints': {}, 'dilatmx': 1.05, 'ecutsm': 0.5, 'ionmov': <pymatgen.io.abinit.abiobjects.MandatoryVariable object>, 'ntime': 80, 'optcell': <pymatgen.io.abinit.abiobjects.MandatoryVariable object>, 'strfact': None, 'strtarget': None, 'tolmxf': None_ )

#### as_dict()
Convert object to dict.


#### _classmethod_ atoms_and_cell(atoms_constraints=None)
Relax atomic positions as well as unit cell.


#### _classmethod_ atoms_only(atoms_constraints=None)
Relax atomic positions, keep unit cell fixed.


#### _classmethod_ from_dict(d)
Build object from dictionary.


#### _property_ move_atoms()
True if atoms must be moved.


#### _property_ move_cell()
True if lattice parameters must be optimized.


#### to_abivars()
Returns a dictionary with the abinit variables.


### _class_ Screening(ecuteps, nband, w_type='RPA', sc_mode='one_shot', hilbert=None, ecutwfn=None, inclvkb=2)
Bases: `AbivarAble`

This object defines the parameters used for the
computation of the screening function.


* **Parameters**


    * **ecuteps** – Cutoff energy for the screening (Ha units).


    * **function** (*nband Number** of **bands for the Green's*) –


    * **w_type** – Screening type


    * **sc_mode** – Self-consistency mode.


    * **hilbert** – Instance of HilbertTransform defining the parameters for the Hilber transform method.


    * **ecutwfn** – Cutoff energy for the wavefunctions (Default: ecutwfn == ecut).


    * **inclvkb** – Option for the treatment of the dipole matrix elements (NC pseudos).



#### _SC_MODES(_ = {'energy_only': 1, 'one_shot': 0, 'wavefunctions': 2_ )

#### _WTYPES(_ = {'RPA': 0_ )

#### _abc_impl(_ = <_abc._abc_data object_ )

#### to_abivars()
Returns a dictionary with the abinit variables.


#### _property_ use_hilbert()
True if we are using the Hilbert transform method.


### _class_ SelfEnergy(se_type, sc_mode, nband, ecutsigx, screening, gw_qprange=1, ppmodel=None, ecuteps=None, ecutwfn=None, gwpara=2)
Bases: `AbivarAble`

This object defines the parameters used for the computation of the self-energy.


* **Parameters**


    * **se_type** – Type of self-energy (str)


    * **sc_mode** – Self-consistency mode.


    * **nband** – Number of bands for the Green’s function


    * **ecutsigx** – Cutoff energy for the exchange part of the self-energy (Ha units).


    * **screening** – Screening instance.


    * **gw_qprange** – Option for the automatic selection of k-points and bands for GW corrections.
    See Abinit docs for more detail. The default value makes the code computie the
    QP energies for all the point in the IBZ and one band above and one band below the Fermi level.


    * **ppmodel** – PPModel instance with the parameters used for the plasmon-pole technique.


    * **ecuteps** – Cutoff energy for the screening (Ha units).


    * **ecutwfn** – Cutoff energy for the wavefunctions (Default: ecutwfn == ecut).



#### _SC_MODES(_ = {'energy_only': 1, 'one_shot': 0, 'wavefunctions': 2_ )

#### _SIGMA_TYPES(_ = {'cohsex': 7, 'gw': 0, 'hartree_fock': 5, 'model_gw_cd': 9, 'model_gw_ppm': 8, 'sex': 6_ )

#### _abc_impl(_ = <_abc._abc_data object_ )

#### _property_ gwcalctyp()
Returns the value of the gwcalctyp input variable.


#### _property_ symsigma()
1 if symmetries can be used to reduce the number of q-points.


#### to_abivars()
Returns a dictionary with the abinit variables.


#### _property_ use_ppmodel()
True if we are using the plasmon-pole approximation.


### _class_ Smearing(occopt, tsmear)
Bases: `AbivarAble`, `MSONable`

Variables defining the smearing technique. The preferred way to instantiate
a Smearing object is via the class method Smearing.as_smearing(string).

Build object with occopt and tsmear.


#### _abc_impl(_ = <_abc._abc_data object_ )

#### _mode2occopt(_ = {'fermi_dirac': 3, 'gaussian': 7, 'marzari4': 4, 'marzari5': 5, 'methfessel': 6, 'nosmearing': 1_ )
Mapping string_mode –> occopt


#### as_dict()
JSON-friendly dict representation of Smearing.


#### _classmethod_ as_smearing(obj)
Constructs an instance of Smearing from obj. Accepts obj in the form:

>
> * Smearing instance


> * “name:tsmear”  e.g. “gaussian:0.004”  (Hartree units)


> * “name:tsmear units” e.g. “gaussian:0.1 eV”


> * None –> no smearing


#### _static_ from_dict(d)
Build object from dict.


#### _property_ mode()
String with smearing technique.


#### _static_ nosmearing()
Build object for calculations without smearing.


#### to_abivars()
Return dictionary with Abinit variables.


### _class_ SpinMode(mode, nsppol, nspinor, nspden)
Bases: `SpinMode`, `AbivarAble`, `MSONable`

Different configurations of the electron density as implemented in abinit:
One can use as_spinmode to construct the object via SpinMode.as_spinmode
(string) where string can assume the values:

>
> * polarized


> * unpolarized


> * afm (anti-ferromagnetic)


> * spinor (non-collinear magnetism)


> * spinor_nomag (non-collinear, no magnetism)

Create new instance of SpinMode(mode, nsppol, nspinor, nspden)


#### _abc_impl(_ = <_abc._abc_data object_ )

#### as_dict()
Convert object to dict.


#### _classmethod_ as_spinmode(obj)
Converts obj into a SpinMode instance.


#### _classmethod_ from_dict(d)
Build object from dict.


#### to_abivars()
Dictionary with Abinit input variables.


### contract(s)
assert contract(“1 1 1 2 2 3”) == “3\*1 2\*2 1\*3”
assert contract(“1 1 3 2 3”) == “2\*1 1\*3 1\*2 1\*3”


### lattice_from_abivars(cls=None, \*args, \*\*kwargs)
Returns a Lattice object from a dictionary
with the Abinit variables acell and either rprim in Bohr or angdeg
If acell is not given, the Abinit default is used i.e. [1,1,1] Bohr.


* **Parameters**

    **cls** – Lattice class to be instantiated. Defaults to pymatgen.core.Lattice.


### Example

lattice_from_abivars(acell=3\*[10], rprim=np.eye(3))


### species_by_znucl(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))
Return list of unique specie found in structure **ordered according to sites**.

### Example

Site0: 0.5 0 0 O
Site1: 0   0 0 Si

produces [Specie_O, Specie_Si] and not set([Specie_O, Specie_Si]) as in types_of_specie


### structure_from_abivars(cls=None, \*args, \*\*kwargs)
Build a Structure object from a dictionary with ABINIT variables.


* **Parameters**

    **cls** – Structure class to be instantiated. pymatgen.core.structure.Structure if cls is None


### Example

al_structure = structure_from_abivars(

    acell=3\*[7.5],
    rprim=[0.0, 0.5, 0.5,

    > 0.5, 0.0, 0.5,
    > 0.5, 0.5, 0.0],

    typat=1,
    xred=[0.0, 0.0, 0.0],
    ntypat=1,
    znucl=13,

)

xred can be replaced with xcart or xangst.


### structure_to_abivars(structure, enforce_znucl=None, enforce_typat=None, \*\*kwargs)
Receives a structure and returns a dictionary with ABINIT variables.


* **Parameters**


    * **enforce_znucl** – List of ntypat entries with the value of Z for each type of atom.
    Used to change the default ordering.


    * **enforce_typat** – List with natom entries with the type index.
    Fortran conventions: start to count from 1.
    Used to change the default ordering.


## pymatgen.io.abinit.abitimer module

This module provides objects for extracting timing data from the ABINIT output files
It also provides tools to analyze and to visualize the parallel efficiency.


### _class_ AbinitTimer(sections, info, cpu_time, wall_time)
Bases: `object`

Container class storing the timing results.


* **Parameters**


    * **sections** – List of sections


    * **info** – Dictionary with extra info.


    * **cpu_time** – Cpu-time in seconds.


    * **wall_time** – Wall-time in seconds.



#### _reduce_sections(keys, operator)

#### cpuwall_histogram(ax: Axes = None, \*\*kwargs)
Plot histogram with cpu- and wall-time on axis ax.


* **Parameters**

    **ax** – matplotlib Axes or None if a new figure should be created.



* **Returns**

    matplotlib figure



* **Return type**

    plt.Figure


Keyword arguments controlling the display of the figure:

| kwargs

 | Meaning

 |
| ------------ | ------------------------------------------------------------------------------------------------- |  |  |  |  |  |  |
| title

        | Title of the plot (Default: None).

                                                                |
| show

         | True to show the figure (default: True).

                                                          |
| savefig

      | “abc.png” or “abc.eps” to save the figure to a file.

                                              |
| size_kwargs

  | Dictionary with options passed to fig.set_size_inches
e.g. size_kwargs=dict(w=3, h=4)

             |
| tight_layout

 | True to call fig.tight_layout (default: False)

                                                    |
| ax_grid

      | True (False) to add (remove) grid from all axes in fig.
Default: None i.e. fig is left unchanged.

 |
| ax_annotate

  | Add labels to  subplots e.g. (a), (b).
Default: False

                                             |
| fig_close

    | Close figure. Default: False.

                                                                     |

#### get_dataframe(sort_key='wall_time', \*\*kwargs)
Return a pandas DataFrame with entries sorted according to sort_key.


#### get_section(section_name)
Return section associated to section_name.


#### get_values(keys)
Return a list of values associated to a particular list of keys.


#### names_and_values(key, minval=None, minfract=None, sorted=True)
Select the entries whose value[key] is >= minval or whose fraction[key] is >= minfract
Return the names of the sections and the corresponding values.


#### _property_ ncpus()
Total number of CPUs employed.


#### order_sections(key, reverse=True)
Sort sections according to the value of key.


#### pie(key='wall_time', minfract=0.05, ax: Axes = None, \*\*kwargs)
Plot pie chart for this timer.


* **Parameters**


    * **key** – Keyword used to extract data from the timer.


    * **minfract** – Don’t show sections whose relative weight is less that minfract.


    * **ax** – matplotlib Axes or None if a new figure should be created.



* **Returns**

    matplotlib figure



* **Return type**

    plt.Figure


Keyword arguments controlling the display of the figure:

| kwargs

       | Meaning

                                                                                           |
| ------------ | ------------------------------------------------------------------------------------------------- |
| title

        | Title of the plot (Default: None).

                                                                |
| show

         | True to show the figure (default: True).

                                                          |
| savefig

      | “abc.png” or “abc.eps” to save the figure to a file.

                                              |
| size_kwargs

  | Dictionary with options passed to fig.set_size_inches
e.g. size_kwargs=dict(w=3, h=4)

             |
| tight_layout

 | True to call fig.tight_layout (default: False)

                                                    |
| ax_grid

      | True (False) to add (remove) grid from all axes in fig.
Default: None i.e. fig is left unchanged.

 |
| ax_annotate

  | Add labels to  subplots e.g. (a), (b).
Default: False

                                             |
| fig_close

    | Close figure. Default: False.

                                                                     |

#### scatter_hist(ax: Axes = None, \*\*kwargs)
Scatter plot + histogram.


* **Parameters**

    **ax** – matplotlib Axes or None if a new figure should be created.



* **Returns**

    matplotlib figure



* **Return type**

    plt.Figure


Keyword arguments controlling the display of the figure:

| kwargs

       | Meaning

                                                                                           |
| ------------ | ------------------------------------------------------------------------------------------------- |
| title

        | Title of the plot (Default: None).

                                                                |
| show

         | True to show the figure (default: True).

                                                          |
| savefig

      | “abc.png” or “abc.eps” to save the figure to a file.

                                              |
| size_kwargs

  | Dictionary with options passed to fig.set_size_inches
e.g. size_kwargs=dict(w=3, h=4)

             |
| tight_layout

 | True to call fig.tight_layout (default: False)

                                                    |
| ax_grid

      | True (False) to add (remove) grid from all axes in fig.
Default: None i.e. fig is left unchanged.

 |
| ax_annotate

  | Add labels to  subplots e.g. (a), (b).
Default: False

                                             |
| fig_close

    | Close figure. Default: False.

                                                                     |

#### sum_sections(keys)
Sum value of keys.


#### to_csv(fileobj=<_io.TextIOWrapper name='<stdout>' mode='w' encoding='utf-8'>)
Write data on file fileobj using CSV format.


#### to_table(sort_key='wall_time', stop=None)
Return a table (list of lists) with timer data.


#### totable(sort_key='wall_time', stop=None)
Return a table (list of lists) with timer data.


### _exception_ AbinitTimerParseError()
Bases: [`ParseError`](pymatgen.io.md#pymatgen.io.core.ParseError)

Errors raised by AbinitTimerParser.


### _class_ AbinitTimerParser()
Bases: `Iterable`

Responsible for parsing a list of output files, extracting the timing results
and analyzing the results.
Assume the Abinit output files have been produced with timopt -1.

### Example

parser = AbinitTimerParser()
parser.parse(list_of_files)

To analyze all

```
*
```

.abo files within top, use:

> parser, paths, okfiles = AbinitTimerParser.walk(top=”.”, ext=”.abo”)

Initialize object.


#### BEGIN_TAG(_ = '-<BEGIN_TIMER_ )

#### END_TAG(_ = '-<END_TIMER>_ )

#### Error()
alias of `AbinitTimerParseError`


#### _abc_impl(_ = <_abc._abc_data object_ )

#### _read(fh, fname)
Parse the TIMER section.


#### _property_ filenames()
List of files that have been parsed successfully.


#### get_sections(section_name)
Return the list of sections stored in self.timers() given section_name
A fake section is returned if the timer does not have section_name.


#### parse(filenames)
Read and parse a filename or a list of filenames.
Files that cannot be opened are ignored. A single filename may also be given.

Return: list of successfully read files.


#### pefficiency()
Analyze the parallel efficiency.

Return: ParallelEfficiency object.


#### plot_all(show=True, \*\*kwargs)
Call all plot methods provided by the parser.


#### plot_efficiency(key='wall_time', what='good+bad', nmax=5, ax: Axes = None, \*\*kwargs)
Plot the parallel efficiency.


* **Parameters**


    * **key** – Parallel efficiency is computed using the wall_time.


    * **what** – Specifies what to plot: good for sections with good parallel efficiency.
    bad for sections with bad efficiency. Options can be concatenated with +.


    * **nmax** – Maximum number of entries in plot


    * **ax** – matplotlib Axes or None if a new figure should be created.


| kwargs

       | Meaning

                                                                                           |
| ------------ | ------------------------------------------------------------------------------------------------- |
| linewidth

    | matplotlib linewidth. Default: 2.0

                                                                |
| markersize

   | matplotlib markersize. Default: 10

                                                                |

* **Returns**

    matplotlib figure


Keyword arguments controlling the display of the figure:

| kwargs

       | Meaning

                                                                                           |
| ------------ | ------------------------------------------------------------------------------------------------- |
| title

        | Title of the plot (Default: None).

                                                                |
| show

         | True to show the figure (default: True).

                                                          |
| savefig

      | “abc.png” or “abc.eps” to save the figure to a file.

                                              |
| size_kwargs

  | Dictionary with options passed to fig.set_size_inches
e.g. size_kwargs=dict(w=3, h=4)

             |
| tight_layout

 | True to call fig.tight_layout (default: False)

                                                    |
| ax_grid

      | True (False) to add (remove) grid from all axes in fig.
Default: None i.e. fig is left unchanged.

 |
| ax_annotate

  | Add labels to  subplots e.g. (a), (b).
Default: False

                                             |
| fig_close

    | Close figure. Default: False.

                                                                     |

#### plot_pie(key='wall_time', minfract=0.05, \*\*kwargs)
Plot pie charts of the different timers.


* **Parameters**


    * **key** – Keyword used to extract data from timers.


    * **minfract** – Don’t show sections whose relative weight is less that minfract.



* **Returns**

    matplotlib figure


Keyword arguments controlling the display of the figure:

| kwargs

       | Meaning

                                                                                           |
| ------------ | ------------------------------------------------------------------------------------------------- |
| title

        | Title of the plot (Default: None).

                                                                |
| show

         | True to show the figure (default: True).

                                                          |
| savefig

      | “abc.png” or “abc.eps” to save the figure to a file.

                                              |
| size_kwargs

  | Dictionary with options passed to fig.set_size_inches
e.g. size_kwargs=dict(w=3, h=4)

             |
| tight_layout

 | True to call fig.tight_layout (default: False)

                                                    |
| ax_grid

      | True (False) to add (remove) grid from all axes in fig.
Default: None i.e. fig is left unchanged.

 |
| ax_annotate

  | Add labels to  subplots e.g. (a), (b).
Default: False

                                             |
| fig_close

    | Close figure. Default: False.

                                                                     |

#### plot_stacked_hist(key='wall_time', nmax=5, ax: Axes = None, \*\*kwargs)
Plot stacked histogram of the different timers.


* **Parameters**


    * **key** – Keyword used to extract data from the timers. Only the first nmax
    sections with largest value are show.


    * **nmax** – Maximum number of sections to show. Other entries are grouped together
    in the others section.


    * **ax** – matplotlib Axes or None if a new figure should be created.



* **Returns**

    matplotlib figure


Keyword arguments controlling the display of the figure:

| kwargs

       | Meaning

                                                                                           |
| ------------ | ------------------------------------------------------------------------------------------------- |
| title

        | Title of the plot (Default: None).

                                                                |
| show

         | True to show the figure (default: True).

                                                          |
| savefig

      | “abc.png” or “abc.eps” to save the figure to a file.

                                              |
| size_kwargs

  | Dictionary with options passed to fig.set_size_inches
e.g. size_kwargs=dict(w=3, h=4)

             |
| tight_layout

 | True to call fig.tight_layout (default: False)

                                                    |
| ax_grid

      | True (False) to add (remove) grid from all axes in fig.
Default: None i.e. fig is left unchanged.

 |
| ax_annotate

  | Add labels to  subplots e.g. (a), (b).
Default: False

                                             |
| fig_close

    | Close figure. Default: False.

                                                                     |

#### section_names(ordkey='wall_time')
Return the names of sections ordered by ordkey.
For the time being, the values are taken from the first timer.


#### summarize(\*\*kwargs)
Return pandas DataFrame with the most important results stored in the timers.


#### timers(filename=None, mpi_rank='0')
Return the list of timers associated to the given filename and MPI rank mpi_rank.


#### _classmethod_ walk(top='.', ext='.abo')
Scan directory tree starting from top, look for files with extension ext and
parse timing data.

Return: (parser, paths, okfiles)

    where parser is the new object, paths is the list of files found and okfiles
    is the list of files that have been parsed successfully.
    (okfiles == paths) if all files have been parsed.


### _class_ AbinitTimerSection(name, cpu_time, cpu_fract, wall_time, wall_fract, ncalls, gflops)
Bases: `object`

Record with the timing results associated to a section of code.


* **Parameters**


    * **name** – Name of the sections.


    * **cpu_time** – CPU time in seconds.


    * **cpu_fract** – Percentage of CPU time.


    * **wall_time** – Wall-time in seconds.


    * **wall_fract** – Percentage of wall-time.


    * **ncalls** – Number of calls


    * **gflops** – Gigaflops.



#### FIELDS(_ = ('name', 'wall_time', 'wall_fract', 'cpu_time', 'cpu_fract', 'ncalls', 'gflops'_ )

#### NUMERIC_FIELDS(_ = ('wall_time', 'wall_fract', 'cpu_time', 'cpu_fract', 'ncalls', 'gflops'_ )

#### STR_FIELDS(_ = ('name',_ )

#### _classmethod_ fake()
Return a fake section. Mainly used to fill missing entries if needed.


#### to_csvline(with_header=False)
Return a string with data in CSV format. Add header if with_header.


#### to_dict()
Convert object to dictionary.


#### to_tuple()
Convert object to tuple.


### _class_ ParallelEfficiency(filenames, ref_idx, \*args, \*\*kwargs)
Bases: `dict`

Store results concerning the parallel efficiency of the job.


* **Parameters**


    * **filennames** – List of filenames


    * **ref_idx** – Index of the Reference time (calculation done with the smallest number of cpus).



#### _order_by_peff(key, criterion, reverse=True)

#### bad_sections(key='wall_time', criterion='mean', nmax=5)
Return first nmax sections with worst value of key key using criterion criterion.


#### good_sections(key='wall_time', criterion='mean', nmax=5)
Return first nmax sections with best value of key key using criterion criterion.


#### totable(stop=None, reverse=True)
Return table (list of lists) with timing results.


* **Parameters**


    * **stop** – Include results up to stop. None for all


    * **reverse** – Put items with highest wall_time in first positions if True.



### alternate(\*iterables)
[a[0], b[0], … , a[1], b[1], …, a[n], b[n] …]
>>> alternate([1,4], [2,5], [3,6])
[1, 2, 3, 4, 5, 6].

## pymatgen.io.abinit.inputs module

This module defines a simplified interface for generating ABINIT input files.
Note that not all the features of Abinit are supported by BasicAbinitInput.
For a more comprehensive implementation, use the AbinitInput object provided by AbiPy.


### _class_ AbstractInput()
Bases: `MutableMapping`

Abstract class defining the methods that must be implemented by Input objects.


#### _abc_impl(_ = <_abc._abc_data object_ )

#### _abstract_ _check_varname(key)
Check if key is a valid name. Raise self.Error if not valid.


#### deepcopy()
Deep copy of the input.


#### pop_vars(keys)
Remove the variables listed in keys.
Return dictionary with the variables that have been removed.
Unlike remove_vars, no exception is raised if the variables are not in the input.


* **Parameters**

    **keys** – string or list of strings with variable names.


### Example

inp.pop_vars([“ionmov”, “optcell”, “ntime”, “dilatmx”])


#### remove_vars(keys, strict=True)
Remove the variables listed in keys.
Return dictionary with the variables that have been removed.


* **Parameters**


    * **keys** – string or list of strings with variable names.


    * **strict** – If True, KeyError is raised if at least one variable is not present.



#### set_vars(\*args, \*\*kwargs)
Set the value of the variables.
Return dict with the variables added to the input.

### Example

input.set_vars(ecut=10, ionmov=3)


#### set_vars_ifnotin(\*args, \*\*kwargs)
Set the value of the variables but only if the variable is not already present.
Return dict with the variables added to the input.

### Example

input.set_vars(ecut=10, ionmov=3)


#### _abstract_ to_str()
Returns a string with the input.


#### _abstract property_ vars()
Dictionary with the input variables. Used to implement dict-like interface.


#### write(filepath='run.abi')
Write the input file to file to filepath.


### _class_ BasicAbinitInput(structure, pseudos, pseudo_dir=None, comment=None, abi_args=None, abi_kwargs=None)
Bases: `AbstractInput`, `MSONable`

This object stores the ABINIT variables for a single dataset.


* **Parameters**


    * **structure** (*file with*) – Parameters defining the crystalline structure. Accepts

    ```
    |Structure|
    ```

     object


    * **structure** –


    * **pseudos** – Pseudopotentials to be used for the calculation. Accepts: string or list of strings
    with the name of the pseudopotential files, list of

    ```
    |Pseudo|
    ```

     objects
    or

    ```
    |PseudoTable|
    ```

     object.


    * **pseudo_dir** – Name of the directory where the pseudopotential files are located.


    * **ndtset** – Number of datasets.


    * **comment** – Optional string with a comment that will be placed at the beginning of the file.


    * **abi_args** – list of tuples (key, value) with the initial set of variables. Default: Empty


    * **abi_kwargs** – Dictionary with the initial set of variables. Default: Empty.



#### Error()
alias of `BasicAbinitInputError`


#### _abc_impl(_ = <_abc._abc_data object_ )

#### _check_varname(key)
Check if key is a valid name. Raise self.Error if not valid.


#### add_abiobjects(\*abi_objects)
This function receive a list of AbiVarable objects and add
the corresponding variables to the input.


#### as_dict()
JSON interface used in pymatgen for easier serialization.


#### _property_ comment()
Optional string with comment. None if comment is not set.


#### _classmethod_ from_dict(d)
JSON interface used in pymatgen for easier serialization.


#### _property_ isnc()
True if norm-conserving calculation.


#### _property_ ispaw()
True if PAW calculation.


#### new_with_vars(\*args, \*\*kwargs)
Return a new input with the given variables.

### Example

new = input.new_with_vars(ecut=20)


#### pop_irdvars()
Remove all the ird\* variables present in self.
Return dictionary with the variables that have been removed.


#### pop_tolerances()
Remove all the tolerance variables present in self.
Return dictionary with the variables that have been removed.


#### _property_ pseudos()
List of

```
|Pseudo|
```

 objects.


#### set_comment(comment)
Set a comment to be included at the top of the file.


#### set_gamma_sampling()
Gamma-only sampling of the BZ.


#### set_kmesh(ngkpt, shiftk, kptopt=1)
Set the variables for the sampling of the BZ.


* **Parameters**


    * **ngkpt** – Monkhorst-Pack divisions


    * **shiftk** – List of shifts.


    * **kptopt** – Option for the generation of the mesh.



#### set_kpath(ndivsm, kptbounds=None, iscf=-2)
Set the variables for the computation of the electronic band structure.


* **Parameters**


    * **ndivsm** – Number of divisions for the smallest segment.


    * **kptbounds** – k-points defining the path in k-space.
    If None, we use the default high-symmetry k-path defined in the pymatgen database.



#### set_spin_mode(spin_mode)
Set the variables used to the treat the spin degree of freedom.
Return dictionary with the variables that have been removed.


* **Parameters**


    * **spin_mode** – SpinMode object or string. Possible values for string are:


    * **polarized** (*-*) –


    * **unpolarized** (*-*) –


    * **afm** (*-*) –


    * **spinor** (*-*) –


    * **spinor_nomag** (*-*) –



#### set_structure(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))
Set structure.


#### _property_ structure()
The

```
|Structure|
```

 object associated to this input.


#### to_str(post=None, with_structure=True, with_pseudos=True, exclude=None)
String representation.


* **Parameters**


    * **post** – String that will be appended to the name of the variables
    Note that post is usually autodetected when we have multiple datatasets
    It is mainly used when we have an input file with a single dataset
    so that we can prevent the code from adding “1” to the name of the variables
    (In this case, indeed, Abinit complains if ndtset=1 is not specified
    and we don’t want ndtset=1 simply because the code will start to add
    _DS1_ to all the input and output files.


    * **with_structure** – False if section with structure variables should not be printed.


    * **with_pseudos** – False if JSON section with pseudo data should not be added.


    * **exclude** – List of variable names that should be ignored.



#### to_string(\*\*kwds)
to_string is deprecated!
Use to_str instead


#### _property_ vars()
Dictionary with variables.


### _exception_ BasicAbinitInputError()
Bases: `Exception`

Base error class for exceptions raised by BasicAbinitInput.


### _class_ BasicMultiDataset(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), pseudos, pseudo_dir='', ndtset=1)
Bases: `object`

This object is essentially a list of BasicAbinitInput objects.
that provides an easy-to-use interface to apply global changes to the
the inputs stored in the objects.

Let’s assume for example that multi contains two BasicAbinitInput objects and we
want to set ecut to 1 in both dictionaries. The direct approach would be:

> for inp in multi:

>     inp.set_vars(ecut=1)

or alternatively:

> for i in range(multi.ndtset):

>     multi[i].set_vars(ecut=1)

BasicMultiDataset provides its own implementation of __getattr__ so that one can simply use:

> multi.set_vars(ecut=1)

> multi.get(“ecut”) returns a list of values. It’s equivalent to:

> > [inp[“ecut”] for inp in multi]

> Note that if “ecut” is not present in one of the input of multi, the corresponding entry is set to None.
> A default value can be specified with:

> > multi.get(“paral_kgb”, 0)

**WARNING**: BasicMultiDataset does not support calculations done with different sets of pseudopotentials.
The inputs can have different crystalline structures (as long as the atom types are equal)
but each input in BasicMultiDataset must have the same set of pseudopotentials.


* **Parameters**


    * **structure** – file with the structure,

    ```
    |Structure|
    ```

     object or dictionary with ABINIT geo variable
    Accepts also list of objects that can be converted to Structure object.
    In this case, however, ndtset must be equal to the length of the list.


    * **pseudos** – String or list of string with the name of the pseudopotential files.


    * **pseudo_dir** – Name of the directory where the pseudopotential files are located.


    * **ndtset** – Number of datasets.



#### Error()
alias of `BasicAbinitInputError`


#### addnew_from(dtindex)
Add a new entry in the multidataset by copying the input with index dtindex.


#### append(abinit_input)
Add a

```
|BasicAbinitInput|
```

 to the list.


#### deepcopy()
Deep copy of the BasicMultiDataset.


#### extend(abinit_inputs)
Extends self with a list of

```
|BasicAbinitInput|
```

 objects.


#### _classmethod_ from_inputs(inputs)
Build object from a list of BasicAbinitInput objects.


#### _property_ has_same_structures()
True if all inputs in BasicMultiDataset are equal.


#### _property_ isnc()
True if norm-conserving calculation.


#### _property_ ispaw()
True if PAW calculation.


#### _property_ ndtset()
Number of inputs in self.


#### _property_ pseudos()
Pseudopotential objects.


#### _classmethod_ replicate_input(input, ndtset)
Construct a multidataset with ndtset from the BasicAbinitInput input.


#### split_datasets()
Return list of

```
|BasicAbinitInput|
```

 objects..


#### to_str(with_pseudos=True)
String representation i.e. the input file read by Abinit.


* **Parameters**

    **with_pseudos** – False if JSON section with pseudo data should not be added.



#### to_string(\*\*kwds)
to_string is deprecated!
Use to_str instead


#### write(filepath='run.abi')
Write ndset input files to disk. The name of the file
is constructed from the dataset index e.g. run0.abi.


### _class_ ShiftMode(value)
Bases: `Enum`

Class defining the mode to be used for the shifts.
G: Gamma centered
M: Monkhorst-Pack ((0.5, 0.5, 0.5))
S: Symmetric. Respects the chksymbreak with multiple shifts
O: OneSymmetric. Respects the chksymbreak with a single shift (as in ‘S’ if a single shift is given, gamma

> centered otherwise.


#### GammaCentered(_ = 'G_ )

#### MonkhorstPack(_ = 'M_ )

#### OneSymmetric(_ = 'O_ )

#### Symmetric(_ = 'S_ )

#### _classmethod_ from_object(obj)
Returns an instance of ShiftMode based on the type of object passed. Converts strings to ShiftMode depending
on the initial letter of the string. G for GammaCentered, M for MonkhorstPack,
S for Symmetric, O for OneSymmetric.
Case insensitive.


### _find_ecut_pawecutdg(ecut, pawecutdg, pseudos, accuracy)
Return a

```
|AttrDict|
```

 with the value of ecut and pawecutdg.


### _find_scf_nband(structure, pseudos, electrons, spinat=None)
Find the value of nband.


### _get_shifts(shift_mode, structure)
Gives the shifts based on the selected shift mode and on the symmetry of the structure.
G: Gamma centered
M: Monkhorst-Pack ((0.5, 0.5, 0.5))
S: Symmetric. Respects the chksymbreak with multiple shifts
O: OneSymmetric. Respects the chksymbreak with a single shift (as in ‘S’ if a single shift is given, gamma

> centered otherwise.

Note: for some cases (e.g. body centered tetragonal), both the Symmetric and OneSymmetric may fail to satisfy the

    chksymbreak condition (Abinit input variable).


### _stopping_criterion(run_level, accuracy)
Return the stopping criterion for this run_level with the given accuracy.


### as_structure(obj)
Convert obj into a Structure. Accepts:

>
> * Structure object.


> * Filename


> * Dictionaries (MSONable format or dictionaries with abinit variables).


### calc_shiftk(structure, symprec: float = 0.01, angle_tolerance=5)
Find the values of shiftk and nshiftk appropriated for the sampling of the Brillouin zone.

When the primitive vectors of the lattice do NOT form a FCC or a BCC lattice,
the usual (shifted) Monkhorst-Pack grids are formed by using nshiftk=1 and shiftk 0.5 0.5 0.5 .
This is often the preferred k point sampling. For a non-shifted Monkhorst-Pack grid,
use nshiftk=1 and shiftk 0.0 0.0 0.0, but there is little reason to do that.

When the primitive vectors of the lattice form a FCC lattice, with rprim:

```default
0.0 0.5 0.5
0.5 0.0 0.5
0.5 0.5 0.0
```

the (very efficient) usual Monkhorst-Pack sampling will be generated by using nshiftk= 4 and shiftk:

```default
0.5 0.5 0.5
0.5 0.0 0.0
0.0 0.5 0.0
0.0 0.0 0.5
```

When the primitive vectors of the lattice form a BCC lattice, with rprim:

```default
-0.5  0.5  0.5
 0.5 -0.5  0.5
 0.5  0.5 -0.5
```

the usual Monkhorst-Pack sampling will be generated by using nshiftk= 2 and shiftk:

```default
 0.25  0.25  0.25
-0.25 -0.25 -0.25
```

However, the simple sampling nshiftk=1 and shiftk 0.5 0.5 0.5 is excellent.

For hexagonal lattices with hexagonal axes, e.g. rprim:

```default
 1.0  0.0       0.0
-0.5  sqrt(3)/2 0.0
 0.0  0.0       1.0
```

one can use nshiftk= 1 and shiftk 0.0 0.0 0.5
In rhombohedral axes, e.g. using angdeg 3\*60., this corresponds to shiftk 0.5 0.5 0.5,
to keep the shift along the symmetry axis.


* **Returns**

    Suggested value of shiftk.



### ebands_input(structure, pseudos, kppa=None, nscf_nband=None, ndivsm=15, ecut=None, pawecutdg=None, scf_nband=None, accuracy='normal', spin_mode='polarized', smearing='fermi_dirac:0.1 eV', charge=0.0, scf_algorithm=None, dos_kppa=None)
Returns a

```
|BasicMultiDataset|
```

 object for band structure calculations.


* **Parameters**


    * **structure** –

    ```
    |Structure|
    ```

     object.


    * **pseudos** – List of filenames or list of

    ```
    |Pseudo|
    ```

     objects or

    ```
    |PseudoTable|
    ```

     object.


    * **kppa** – Defines the sampling used for the SCF run. Defaults to 1000 if not given.


    * **nscf_nband** – Number of bands included in the NSCF run. Set to scf_nband + 10 if None.


    * **ndivsm** – Number of divisions used to sample the smallest segment of the k-path.
    if 0, only the GS input is returned in multi[0].


    * **ecut** – cutoff energy in Ha (if None, ecut is initialized from the pseudos according to accuracy)


    * **pawecutdg** – cutoff energy in Ha for PAW double-grid (if None, pawecutdg is initialized from the pseudos
    according to accuracy)


    * **scf_nband** – Number of bands for SCF run. If scf_nband is None, nband is automatically initialized
    from the list of pseudos, the structure and the smearing option.


    * **accuracy** – Accuracy of the calculation.


    * **spin_mode** – Spin polarization.


    * **smearing** – Smearing technique.


    * **charge** – Electronic charge added to the unit cell.


    * **scf_algorithm** – Algorithm used for solving of the SCF cycle.


    * **dos_kppa** – Scalar or List of integers with the number of k-points per atom
    to be used for the computation of the DOS (None if DOS is not wanted).



### gs_input(structure, pseudos, kppa=None, ecut=None, pawecutdg=None, scf_nband=None, accuracy='normal', spin_mode='polarized', smearing='fermi_dirac:0.1 eV', charge=0.0, scf_algorithm=None)
Returns a

```
|BasicAbinitInput|
```

 for ground-state calculation.


* **Parameters**


    * **structure** –

    ```
    |Structure|
    ```

     object.


    * **pseudos** – List of filenames or list of

    ```
    |Pseudo|
    ```

     objects or

    ```
    |PseudoTable|
    ```

     object.


    * **kppa** – Defines the sampling used for the SCF run. Defaults to 1000 if not given.


    * **ecut** – cutoff energy in Ha (if None, ecut is initialized from the pseudos according to accuracy)


    * **pawecutdg** – cutoff energy in Ha for PAW double-grid (if None, pawecutdg is initialized from the pseudos
    according to accuracy)


    * **scf_nband** – Number of bands for SCF run. If scf_nband is None, nband is automatically initialized
    from the list of pseudos, the structure and the smearing option.


    * **accuracy** – Accuracy of the calculation.


    * **spin_mode** – Spin polarization.


    * **smearing** – Smearing technique.


    * **charge** – Electronic charge added to the unit cell.


    * **scf_algorithm** – Algorithm used for solving of the SCF cycle.



### ion_ioncell_relax_input(structure, pseudos, kppa=None, nband=None, ecut=None, pawecutdg=None, accuracy='normal', spin_mode='polarized', smearing='fermi_dirac:0.1 eV', charge=0.0, scf_algorithm=None, shift_mode='Monkhorst-pack')
Returns a

```
|BasicMultiDataset|
```

 for a structural relaxation. The first dataset optmizes the
atomic positions at fixed unit cell. The second datasets optimizes both ions and unit cell parameters.


* **Parameters**


    * **structure** –

    ```
    |Structure|
    ```

     object.


    * **pseudos** – List of filenames or list of

    ```
    |Pseudo|
    ```

     objects or

    ```
    |PseudoTable|
    ```

     object.


    * **kppa** – Defines the sampling used for the Brillouin zone.


    * **nband** – Number of bands included in the SCF run.


    * **accuracy** – Accuracy of the calculation.


    * **spin_mode** – Spin polarization.


    * **smearing** – Smearing technique.


    * **charge** – Electronic charge added to the unit cell.


    * **scf_algorithm** – Algorithm used for the solution of the SCF cycle.



### num_valence_electrons(structure, pseudos)
Returns the number of valence electrons.


* **Parameters**

    **pseudos** – List of

    ```
    |Pseudo|
    ```

     objects or list of filenames.


## pymatgen.io.abinit.netcdf module

Wrapper for netCDF readers.


### _class_ AbinitHeader(\*args, \*\*kwargs)
Bases: `AttrDict`

Stores the values reported in the Abinit header.


* **Parameters**


    * **args** – Passthrough arguments for standard dict.


    * **kwargs** – Passthrough keyword arguments for standard dict.



#### to_str(verbose=0, title=None, \*\*kwargs)
String representation. kwargs are passed to pprint.pformat.


* **Parameters**


    * **verbose** – Verbosity level


    * **title** – Title string.



#### to_string(\*\*kwds)
to_string is deprecated!
Use to_str instead


### _class_ ETSF_Reader(path)
Bases: `NetcdfReader`

This object reads data from a file written according to the ETSF-IO specifications.

We assume that the netcdf file contains at least the crystallographic section.

Open the Netcdf file specified by path (read mode).


#### chemical_symbols()
Chemical symbols char [number of atom species][symbol length].


#### read_abinit_hdr()
Read the variables associated to the Abinit header.

Return AbinitHeader


#### read_abinit_xcfunc()
Read ixc from an Abinit file. Return XcFunc object.


#### read_structure(cls=<class 'pymatgen.core.structure.Structure'>)
Returns the crystalline structure stored in the rootgrp.


#### typeidx_from_symbol(symbol)
Returns the type index from the chemical symbol. Note python convention.


### _class_ NO_DEFAULT()
Bases: `object`

Signal that read_value should raise an Error.


### _class_ NetcdfReader(path)
Bases: `object`

Wraps and extends netCDF4.Dataset. Read only mode. Supports with statements.

Additional documentation available at:

    [http://netcdf4-python.googlecode.com/svn/trunk/docs/netCDF4-module.html](http://netcdf4-python.googlecode.com/svn/trunk/docs/netCDF4-module.html)

Open the Netcdf file specified by path (read mode).


#### Error()
alias of `NetcdfReaderError`


#### _read_dimensions(\*dim_names, \*\*kwargs)

#### _read_variables(\*var_names, \*\*kwargs)

#### close()
Close the file.


#### print_tree()
Print all the groups in the file.


#### read_dimvalue(dimname, path='/', default=<class 'pymatgen.io.abinit.netcdf.NO_DEFAULT'>)
Returns the value of a dimension.


* **Parameters**


    * **dimname** – Name of the variable


    * **path** – path to the group.


    * **default** – return default if dimname is not present and
    default is not NO_DEFAULT else raise self.Error.



#### read_keys(keys, dict_cls=<class 'monty.collections.AttrDict'>, path='/')
Read a list of variables/dimensions from file. If a key is not present the corresponding
entry in the output dictionary is set to None.


#### read_value(varname, path='/', cmode=None, default=<class 'pymatgen.io.abinit.netcdf.NO_DEFAULT'>)
Returns the values of variable with name varname in the group specified by path.


* **Parameters**


    * **varname** – Name of the variable


    * **path** – path to the group.


    * **cmode** – if cmode==”c”, a complex ndarrays is constructed and returned
    (netcdf does not provide native support from complex datatype).


    * **default** – returns default if varname is not present.
    self.Error is raised if default is set to NO_DEFAULT



* **Returns**

    numpy array if varname represents an array, scalar otherwise.



#### read_variable(varname, path='/')
Returns the variable with name varname in the group specified by path.


#### read_varnames(path='/')
List of variable names stored in the group specified by path.


#### walk_tree(top=None)
Navigate all the groups in the file starting from top.
If top is None, the root group is used.


### _exception_ NetcdfReaderError()
Bases: `Exception`

Base error class for NetcdfReader.


### _class_ _H(name, doc, etsf_name=None)
Bases: `object`


#### doc()

#### etsf_name()

#### name()

### _asreader(file, cls)

### as_etsfreader(file)
Return an ETSF_Reader. Accepts filename or ETSF_Reader.


### as_ncreader(file)
Convert file into a NetcdfReader instance.
Returns reader, closeit where closeit is set to True
if we have to close the file before leaving the procedure.


### structure_from_ncdata(ncdata, site_properties=None, cls=<class 'pymatgen.core.structure.Structure'>)
Reads and returns a pymatgen structure from a NetCDF file
containing crystallographic data in the ETSF-IO format.


* **Parameters**


    * **ncdata** – filename or NetcdfReader instance.


    * **site_properties** – Dictionary with site properties.


    * **cls** – The Structure class to instantiate.


## pymatgen.io.abinit.pseudos module

This module provides objects describing the basic parameters of the
pseudopotentials used in Abinit, and a parser to instantiate pseudopotential objects..


### _class_ AbinitHeader()
Bases: `dict`

Dictionary whose keys can be also accessed as attributes.


### _class_ AbinitPseudo(path, header)
Bases: `Pseudo`

An AbinitPseudo is a pseudopotential whose file contains an abinit header.


* **Parameters**


    * **path** – Filename.


    * **header** – AbinitHeader instance.



#### _property_ Z()
The atomic number of the atom.


#### _property_ Z_val()
Valence charge.


#### _abc_impl(_ = <_abc._abc_data object_ )

#### _property_ l_local()
Angular momentum used for the local part.


#### _property_ l_max()
Maximum angular momentum.


#### _property_ summary()
Summary line reported in the ABINIT header.


#### _property_ supports_soc()
True if the pseudo can be used in a calculation with spin-orbit coupling.
Base classes should provide a concrete implementation that computes this value.


### _class_ Hint(ecut, pawecutdg=None)
Bases: `object`

Suggested value for the cutoff energy [Hartree units]
and the cutoff energy for the dense grid (only for PAW pseudos).


#### as_dict()
Return dictionary for MSONable protocol.


#### _classmethod_ from_dict(d)
Build instance from dictionary (MSONable protocol).


### _class_ NcAbinitHeader(summary, \*\*kwargs)
Bases: `AbinitHeader`

The abinit header found in the NC pseudopotential files.


#### _VARS(_ = {'fchrg': (0.0, <class 'float'>), 'lloc': (None, <class 'int'>), 'lmax': (None, <class 'int'>), 'mmax': (None, <class 'float'>), 'pspcod': (None, <class 'int'>), 'pspdat': (None, <class 'float'>), 'pspxc': (None, <class 'int'>), 'qchrg': (0.0, <class 'float'>), 'r2well': (None, <class 'float'>), 'rchrg': (0.0, <class 'float'>), 'zatom': (None, <function _int_from_str>), 'zion': (None, <class 'float'>)_ )

#### _static_ fhi_header(filename, ppdesc)
Parse the FHI abinit header. Example:

Troullier-Martins psp for element  Sc        Thu Oct 27 17:33:22 EDT 1994

    21.00000   3.00000    940714                zatom, zion, pspdat
    1    1    2    0      2001    .00000      pspcod,pspxc,lmax,lloc,mmax,r2well
    1.80626423934776     .22824404341771    1.17378968127746   rchrg,fchrg,qchrg


#### _static_ gth_header(filename, ppdesc)
Parse the GTH abinit header. Example:

Goedecker-Teter-Hutter  Wed May  8 14:27:44 EDT 1996
1   1   960508                     zatom,zion,pspdat
2   1   0    0    2001    0.       pspcod,pspxc,lmax,lloc,mmax,r2well
0.2000000 -4.0663326  0.6778322 0 0     rloc, c1, c2, c3, c4
0 0 0                              rs, h1s, h2s
0 0                                rp, h1p

> 1.36 .2   0.6                    rcutoff, rloc


#### _static_ hgh_header(filename, ppdesc)
Parse the HGH abinit header. Example:

Hartwigsen-Goedecker-Hutter psp for Ne,  from PRB58, 3641 (1998)

    10   8  010605 zatom,zion,pspdat
    3 1   1 0 2001 0  pspcod,pspxc,lmax,lloc,mmax,r2well


#### _static_ oncvpsp_header(filename, ppdesc)
Parse the ONCVPSP abinit header. Example:

Li    ONCVPSP  r_core=  2.01  3.02

    > > 3.0000      3.0000      140504    zatom,zion,pspd

    > 8     2     1     4   600     0    pspcod,pspxc,lmax,lloc,mmax,r2well

    5.99000000  0.00000000  0.00000000    rchrg fchrg qchrg

        > 2     2     0     0     0    nproj
        > 0                 extension_switch

        0                        -2.5000025868368D+00 -1.2006906995331D+00

            1  0.0000000000000D+00  0.0000000000000D+00  0.0000000000000D+00
            2  1.0000000000000D-02  4.4140499497377D-02  1.9909081701712D-02


#### _static_ tm_header(filename, ppdesc)
Parse the TM abinit header. Example:

Troullier-Martins psp for element Fm         Thu Oct 27 17:28:39 EDT 1994
100.00000  14.00000    940714                zatom, zion, pspdat

> 1    1    3    0      2001    .00000      pspcod,pspxc,lmax,lloc,mmax,r2well
> 0   4.085   6.246    0   2.8786493        l,e99.0,e99.9,nproj,rcpsp
> .00000000    .0000000000    .0000000000    .00000000   rms,ekb1,ekb2,epsatm
> 1   3.116   4.632    1   3.4291849        l,e99.0,e99.9,nproj,rcpsp
> .00000000    .0000000000    .0000000000    .00000000   rms,ekb1,ekb2,epsatm
> 2   4.557   6.308    1   2.1865358        l,e99.0,e99.9,nproj,rcpsp
> .00000000    .0000000000    .0000000000    .00000000   rms,ekb1,ekb2,epsatm
> 3  23.251  29.387    1   2.4776730        l,e99.0,e99.9,nproj,rcpsp
> .00000000    .0000000000    .0000000000    .00000000   rms,ekb1,ekb2,epsatm
> 3.62474762267880     .07409391739104    3.07937699839200   rchrg,fchrg,qchrg


### _class_ NcAbinitPseudo(path, header)
Bases: `NcPseudo`, `AbinitPseudo`

Norm-conserving pseudopotential in the Abinit format.


* **Parameters**


    * **path** – Filename.


    * **header** – AbinitHeader instance.



#### _property_ Z()
The atomic number of the atom.


#### _property_ Z_val()
Number of valence electrons.


#### _abc_impl(_ = <_abc._abc_data object_ )

#### _property_ l_local()
Angular momentum used for the local part.


#### _property_ l_max()
Maximum angular momentum.


#### _property_ nlcc_radius()
Radius at which the core charge vanish (i.e. cut-off in a.u.).
Returns 0.0 if nlcc is not used.


#### _property_ summary()
Summary line reported in the ABINIT header.


### _class_ NcPseudo()
Bases: `object`

Abstract class defining the methods that must be implemented
by the concrete classes representing norm-conserving pseudopotentials.


#### _abc_impl(_ = <_abc._abc_data object_ )

#### _property_ has_nlcc()
True if the pseudo is generated with non-linear core correction.


#### _abstract property_ nlcc_radius()
Radius at which the core charge vanish (i.e. cut-off in a.u.).
Returns 0.0 if nlcc is not used.


#### _property_ rcore()
Radius of the pseudization sphere in a.u.


### _class_ PawAbinitHeader(summary, \*\*kwargs)
Bases: `AbinitHeader`

The abinit header found in the PAW pseudopotential files.


#### _VARS(_ = {'basis_size': (None, <class 'int'>), 'creatorID': (None, <class 'int'>), 'lloc': (None, <class 'int'>), 'lmax': (None, <class 'int'>), 'lmn_size': (None, <class 'int'>), 'mmax': (None, <class 'int'>), 'number_of_meshes': (None, <class 'int'>), 'orbitals': (None, <class 'list'>), 'pspcod': (None, <class 'int'>), 'pspdat': (None, <class 'float'>), 'pspfmt': (None, <class 'str'>), 'pspxc': (None, <class 'int'>), 'r2well': (None, <class 'float'>), 'r_cut': (None, <class 'float'>), 'rshape': (None, <class 'float'>), 'shape_type': (None, <class 'int'>), 'zatom': (None, <function _int_from_str>), 'zion': (None, <class 'float'>)_ )

#### _static_ paw_header(filename, ppdesc)
Parse the PAW abinit header. Examples:

Paw atomic data for element Ni - Generated by AtomPAW (N. Holzwarth) + AtomPAW2Abinit v3.0.5

    > 28.000  18.000 20061204               : zatom,zion,pspdat
    > 7  7  2 0   350 0.                    : pspcod,pspxc,lmax,lloc,mmax,r2well

    paw3 1305

        5 13                                  : basis_size,lmn_size

    0 0 1 1 2                              : orbitals
    3                                      : number_of_meshes
    1 3  350 1.1803778368E-05 3.5000000000E-02 : mesh 1, type,size,rad_step[,log_step]
    2 1  921 2.500000000000E-03                : mesh 2, type,size,rad_step[,log_step]
    3 3  391 1.1803778368E-05 3.5000000000E-02 : mesh 3, type,size,rad_step[,log_step]

    > 2.3000000000                          : r_cut(SPH)

    2 0.

Another format:

C  (US d-loc) - PAW data extracted from US-psp (D.Vanderbilt) - generated by USpp2Abinit v2.3.0

    > > 6.000   4.000 20090106               : zatom,zion,pspdat

    > 7 11  1 0   560 0.                    : pspcod,pspxc,lmax,lloc,mmax,r2well

    paw4 2230

        4  8                                  : basis_size,lmn_size

    0 0 1 1                                : orbitals
    5                                      : number_of_meshes
    1 2  560 1.5198032759E-04 1.6666666667E-02 : mesh 1, type,size,rad_step[,log_step]
    2 2  556 1.5198032759E-04 1.6666666667E-02 : mesh 2, type,size,rad_step[,log_step]
    3 2  576 1.5198032759E-04 1.6666666667E-02 : mesh 3, type,size,rad_step[,log_step]
    4 2  666 1.5198032759E-04 1.6666666667E-02 : mesh 4, type,size,rad_step[,log_step]
    5 2  673 1.5198032759E-04 1.6666666667E-02 : mesh 5, type,size,rad_step[,log_step]

    > 1.5550009124                          : r_cut(PAW)

    3 0.                                   : shape_type,rshape

Yet nnother one:

Paw atomic data for element Si - Generated by atompaw v3.0.1.3 & AtomPAW2Abinit v3.3.1

    > 14.000   4.000 20120814               : zatom,zion,pspdat
    > 7      11  1 0   663 0.               : pspcod,pspxc,lmax,lloc,mmax,r2well

    paw5 1331

        4  8                                  : basis_size,lmn_size

    0 0 1 1                                : orbitals
    5                                      : number_of_meshes
    1 2  663 8.2129718540404674E-04 1.1498160595656655E-02 : mesh 1, type,size,rad_step[,log_step]
    2 2  658 8.2129718540404674E-04 1.1498160595656655E-02 : mesh 2, type,size,rad_step[,log_step]
    3 2  740 8.2129718540404674E-04 1.1498160595656655E-02 : mesh 3, type,size,rad_step[,log_step]
    4 2  819 8.2129718540404674E-04 1.1498160595656655E-02 : mesh 4, type,size,rad_step[,log_step]
    5 2  870 8.2129718540404674E-04 1.1498160595656655E-02 : mesh 5, type,size,rad_step[,log_step]

    > 1.5669671236                          : r_cut(PAW)

    2 0.                                   : shape_type,rshape


### _class_ PawAbinitPseudo(path, header)
Bases: `PawPseudo`, `AbinitPseudo`

Paw pseudopotential in the Abinit format.


* **Parameters**


    * **path** – Filename.


    * **header** – AbinitHeader instance.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### _property_ paw_radius()
Radius of the PAW sphere in a.u.


#### _property_ supports_soc()
True if the pseudo can be used in a calculation with spin-orbit coupling.
Base classes should provide a concrete implementation that computes this value.


### _class_ PawPseudo()
Bases: `object`

Abstract class that defines the methods that must be implemented
by the concrete classes representing PAW pseudopotentials.


#### _abc_impl(_ = <_abc._abc_data object_ )

#### _abstract property_ paw_radius()
Radius of the PAW sphere in a.u.


#### _property_ rcore()
Alias of paw_radius.


### _class_ PawXmlSetup(filepath)
Bases: `Pseudo`, `PawPseudo`

Setup class for PawXml.


* **Parameters**

    **filepath** (*str*) – Path to the XML file.



#### _property_ Z()
The atomic number of the atom.


#### _property_ Z_val()
Number of valence electrons.


#### _abc_impl(_ = <_abc._abc_data object_ )

#### _static_ _eval_grid(grid_params)
This function receives a dictionary with the parameters defining the
radial mesh and returns a ndarray with the mesh.


#### _parse_all_radfuncs(func_name)
Parse all the nodes with tag func_name in the XML file.


#### _parse_radfunc(func_name)
Parse the first occurrence of func_name in the XML file.


#### ae_core_density()
The all-electron radial density.


#### ae_partial_waves()
Dictionary with the AE partial waves indexed by state.


#### _property_ l_local()
Angular momentum used for the local part.


#### _property_ l_max()
Maximum angular momentum.


#### _property_ paw_radius()
Radius of the PAW sphere in a.u.


#### plot_densities(ax: plt.Axes = None, \*\*kwargs)
Plot the PAW densities.


* **Parameters**

    **ax** – matplotlib Axes or None if a new figure should be created.



* **Returns**

    matplotlib figure


Keyword arguments controlling the display of the figure:

| kwargs

       | Meaning

                                                                                           |
| ------------ | ------------------------------------------------------------------------------------------------- |
| title

        | Title of the plot (Default: None).

                                                                |
| show

         | True to show the figure (default: True).

                                                          |
| savefig

      | “abc.png” or “abc.eps” to save the figure to a file.

                                              |
| size_kwargs

  | Dictionary with options passed to fig.set_size_inches
e.g. size_kwargs=dict(w=3, h=4)

             |
| tight_layout

 | True to call fig.tight_layout (default: False)

                                                    |
| ax_grid

      | True (False) to add (remove) grid from all axes in fig.
Default: None i.e. fig is left unchanged.

 |
| ax_annotate

  | Add labels to  subplots e.g. (a), (b).
Default: False

                                             |
| fig_close

    | Close figure. Default: False.

                                                                     |

#### plot_projectors(ax: plt.Axes = None, fontsize=12, \*\*kwargs)
Plot the PAW projectors.


* **Parameters**

    **ax** – matplotlib Axes or None if a new figure should be created.



* **Returns**

    matplotlib figure



* **Return type**

    plt.Figure


Keyword arguments controlling the display of the figure:

| kwargs

       | Meaning

                                                                                           |
| ------------ | ------------------------------------------------------------------------------------------------- |
| title

        | Title of the plot (Default: None).

                                                                |
| show

         | True to show the figure (default: True).

                                                          |
| savefig

      | “abc.png” or “abc.eps” to save the figure to a file.

                                              |
| size_kwargs

  | Dictionary with options passed to fig.set_size_inches
e.g. size_kwargs=dict(w=3, h=4)

             |
| tight_layout

 | True to call fig.tight_layout (default: False)

                                                    |
| ax_grid

      | True (False) to add (remove) grid from all axes in fig.
Default: None i.e. fig is left unchanged.

 |
| ax_annotate

  | Add labels to  subplots e.g. (a), (b).
Default: False

                                             |
| fig_close

    | Close figure. Default: False.

                                                                     |

#### plot_waves(ax: plt.Axes = None, fontsize=12, \*\*kwargs)
Plot the AE and the pseudo partial waves.


* **Parameters**


    * **ax** – matplotlib Axes or None if a new figure should be created.


    * **fontsize** – fontsize for legends and titles



* **Returns**

    matplotlib figure



* **Return type**

    plt.Figure


Keyword arguments controlling the display of the figure:

| kwargs

       | Meaning

                                                                                           |
| ------------ | ------------------------------------------------------------------------------------------------- |
| title

        | Title of the plot (Default: None).

                                                                |
| show

         | True to show the figure (default: True).

                                                          |
| savefig

      | “abc.png” or “abc.eps” to save the figure to a file.

                                              |
| size_kwargs

  | Dictionary with options passed to fig.set_size_inches
e.g. size_kwargs=dict(w=3, h=4)

             |
| tight_layout

 | True to call fig.tight_layout (default: False)

                                                    |
| ax_grid

      | True (False) to add (remove) grid from all axes in fig.
Default: None i.e. fig is left unchanged.

 |
| ax_annotate

  | Add labels to  subplots e.g. (a), (b).
Default: False

                                             |
| fig_close

    | Close figure. Default: False.

                                                                     |

#### projector_functions()
Dictionary with the PAW projectors indexed by state.


#### pseudo_core_density()
The pseudized radial density.


#### _property_ pseudo_partial_waves()
Dictionary with the pseudo partial waves indexed by state.


#### root()
Root tree of XML.


#### _property_ summary()
String summarizing the most important properties.


#### _property_ supports_soc()
Here I assume that the ab-initio code can treat the SOC within the on-site approximation.


#### yield_figs(\*\*kwargs)
This function *generates* a predefined list of matplotlib figures with minimal input from the user.


### _class_ Pseudo()
Bases: `MSONable`

Abstract base class defining the methods that must be
implemented by the concrete pseudo-potential sub-classes.


#### _abstract property_ Z(_: in_ )
The atomic number of the atom.


#### _abstract property_ Z_val(_: in_ )
Valence charge.


#### _abc_impl(_ = <_abc._abc_data object_ )

#### as_dict(\*\*kwargs)
Return dictionary for MSONable protocol.


#### _classmethod_ as_pseudo(obj)
Convert obj into a pseudo. Accepts:

>
> * Pseudo object.


> * string defining a valid path.


#### as_tmpfile(tmpdir=None)
Copy the pseudopotential to a temporary a file and returns a new pseudopotential object.
Useful for unit tests in which we have to change the content of the file.


* **Parameters**

    **tmpdir** – If None, a new temporary directory is created and files are copied here
    else tmpdir is used.



#### _property_ basename(_: st_ )
File basename.


#### compute_md5()
Compute and return MD5 hash value.


#### _property_ djrepo_path()
The path of the djrepo file. None if file does not exist.


#### _property_ element(_: [Element](pymatgen.core.md#pymatgen.core.periodic_table.Element_ )
Pymatgen Element.


#### _property_ filepath(_: st_ )
Absolute path to pseudopotential file.


#### _classmethod_ from_dict(dct)
Build instance from dictionary (MSONable protocol).


#### _static_ from_file(filename)
Build an instance of a concrete Pseudo subclass from filename.
Note: the parser knows the concrete class that should be instantiated
Client code should rely on the abstract interface provided by Pseudo.


#### _property_ has_dojo_report()
True if the pseudo has an associated DOJO_REPORT section.


#### _property_ has_hints()
True if self provides hints on the cutoff energy.


#### hint_for_accuracy(accuracy='normal')
Returns a Hint object with the suggested value of ecut [Ha] and
pawecutdg [Ha] for the given accuracy.
ecut and pawecutdg are set to zero if no hint is available.


* **Parameters**

    **accuracy** – [“low”, “normal”, “high”]



#### _property_ isnc(_: boo_ )
True if norm-conserving pseudopotential.


#### _property_ ispaw(_: boo_ )
True if PAW pseudopotential.


#### _abstract property_ l_local(_: in_ )
Angular momentum used for the local part.


#### _abstract property_ l_max(_: in_ )
Maximum angular momentum.


#### md5()
MD5 hash value.


#### open_pspsfile(ecut=20, pawecutdg=None)
Calls Abinit to compute the internal tables for the application of the
pseudopotential part. Returns PspsFile object providing methods
to plot and analyze the data or None if file is not found or it’s not readable.


* **Parameters**


    * **ecut** – Cutoff energy in Hartree.


    * **pawecutdg** – Cutoff energy for the PAW double grid.



#### _abstract property_ summary(_: st_ )
String summarizing the most important properties.


#### _abstract property_ supports_soc()
True if the pseudo can be used in a calculation with spin-orbit coupling.
Base classes should provide a concrete implementation that computes this value.


#### _property_ symbol(_: st_ )
Element symbol.


#### to_str(verbose=0)
String representation.


#### to_string(\*\*kwds)
to_string is deprecated!
Use to_str instead


#### _property_ type(_: st_ )
Type of pseudo.


### _exception_ PseudoParseError()
Bases: [`ParseError`](pymatgen.io.md#pymatgen.io.core.ParseError)

Base Error class for the exceptions raised by PseudoParser.


### _class_ PseudoParser()
Bases: `object`

Responsible for parsing pseudopotential files and returning pseudopotential objects.

Usage:

```default
pseudo = PseudoParser().parse("filename")
```


#### Error()
alias of `PseudoParseError`


#### _PSPCODES(_ = {1: (1, 'TM', 'NC', None), 2: (2, 'GTH', 'NC', None), 3: (3, 'HGH', 'NC', None), 4: (4, 'Teter', 'NC', None), 6: (6, 'FHI', 'NC', None), 7: (6, 'PAW_abinit_text', 'PAW', None), 8: (8, 'ONCVPSP', 'NC', None), 10: (10, 'HGHK', 'NC', None)_ )

#### parse(filename)
Read and parse a pseudopotential file. Main entry point for client code.


* **Returns**

    pseudopotential object or None if filename is not a valid pseudopotential file.



#### read_ppdesc(filename)
Read the pseudopotential descriptor from filename.


* **Returns**

    Pseudopotential descriptor. None if filename is not a valid pseudopotential file.



* **Raises**

    **PseudoParseError** –



#### scan_directory(dirname, exclude_exts=(), exclude_fnames=())
Analyze the files contained in directory dirname.


* **Parameters**


    * **dirname** – directory path


    * **exclude_exts** – list of file extensions that should be skipped.


    * **exclude_fnames** – list of file names that should be skipped.



* **Returns**

    List of pseudopotential objects.



### _class_ PseudoTable(pseudos: Sequence[Pseudo])
Bases: `Sequence`, `MSONable`

Define the pseudopotentials from the element table.
Individidual elements are accessed by name, symbol or atomic number.

For example, the following all retrieve iron:

print elements[26]
Fe
print elements.Fe
Fe
print elements.symbol(‘Fe’)
Fe
print elements.name(‘iron’)
Fe
print elements.isotope(‘Fe’)
Fe


* **Parameters**

    **pseudos** – List of pseudopotentials or filepaths.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### all_combinations_for_elements(element_symbols)
Return a list with all the possible combination of pseudos
for the given list of element_symbols.
Each item is a list of pseudopotential objects.

Example:

```default
table.all_combinations_for_elements(["Li", "F"])
```


#### _property_ allnc(_: boo_ )
True if all pseudos are norm-conserving.


#### _property_ allpaw()
True if all pseudos are PAW.


#### as_dict(\*\*kwargs)
Return dictionary for MSONable protocol.


#### _classmethod_ as_table(items)
Return an instance of PseudoTable from the iterable items.


#### _classmethod_ from_dict(d)
Build instance from dictionary (MSONable protocol).


#### _classmethod_ from_dir(top, exts=None, exclude_dirs='_\*')
Find all pseudos in the directory tree starting from top.


* **Parameters**


    * **top** – Top of the directory tree


    * **exts** – List of files extensions. if exts == “all_files”
    we try to open all files in top


    * **exclude_dirs** – Wildcard used to exclude directories.


return: PseudoTable sorted by atomic number Z.


#### get_pseudos_for_structure(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))
Return the list of Pseudo objects to be used for this Structure.


* **Parameters**

    **structure** – pymatgen Structure.



* **Raises**


    * **ValueError** –


    * **multiple occurrences are present in the table.** –



#### is_complete(zmax=118)
True if table is complete i.e. all elements with Z < zmax have at least on pseudopotential.


#### print_table(stream=<_io.TextIOWrapper name='<stdout>' mode='w' encoding='utf-8'>, filter_function=None)
A pretty ASCII printer for the periodic table, based on some filter_function.


* **Parameters**


    * **stream** – file-like object


    * **filter_function** – A filtering function that take a Pseudo as input and returns a boolean.
    For example, setting filter_function = lambda p: p.Z_val > 2 will print
    a periodic table containing only pseudos with Z_val > 2.



#### pseudo_with_symbol(symbol, allow_multi=False)
Return the pseudo with the given chemical symbol.


* **Parameters**


    * **symbols** – String with the chemical symbol of the element


    * **allow_multi** – By default, the method raises ValueError
    if multiple occurrences are found. Use allow_multi to prevent this.



* **Raises**

    **ValueError if symbol is not found**** or ****multiple occurrences are present and not allow_multi** –



#### pseudos_with_symbols(symbols)
Return the pseudos with the given chemical symbols.


* **Raises**

    **ValueError if one**** of ****the symbols is not found**** or ****multiple occurrences are present.** –



#### select(condition)
Select only those pseudopotentials for which condition is True.


* **Parameters**

    **condition** – Function that accepts a Pseudo object and returns True or False.



* **Returns**

    New PseudoTable instance with pseudos for which condition is True.



* **Return type**

    PseudoTable



#### select_family(family)
Return PseudoTable with element belonging to the specified family, e.g. family=”alkaline”.


#### select_rows(rows)
Return new class:PseudoTable object with pseudos in the given rows of the periodic table.
rows can be either a int or a list of integers.


#### select_symbols(symbols, ret_list=False)
Return a PseudoTable with the pseudopotentials with the given list of chemical symbols.


* **Parameters**


    * **symbols** – str or list of symbols
    Prepend the symbol string with “-”, to exclude pseudos.


    * **ret_list** – if True a list of pseudos is returned instead of a PseudoTable



#### sort_by_z()
Return a new PseudoTable with pseudos sorted by Z.


#### sorted(attrname, reverse=False)
Sort the table according to the value of attribute attrname.


* **Returns**

    PseudoTable object



* **Return type**

    New class



#### to_table(filter_function=None)
Return string with data in tabular form.


#### with_dojo_report()
Select pseudos containing the DOJO_REPORT section. Return new class:PseudoTable object.


#### _property_ zlist()
Ordered list with the atomic numbers available in the table.


### _class_ RadialFunction(mesh, values)
Bases: `RadialFunction`

Radial Function class.

Create new instance of RadialFunction(mesh, values)


### _dict_from_lines(lines, key_nums, sep=None)
Helper function to parse formatted text structured like:

value1 value2 … sep key1, key2 …

key_nums is a list giving the number of keys for each line. 0 if line should be skipped.
sep is a string denoting the character that separates the keys from the value (None if
no separator is present).


* **Returns**

    value1, key2 : value2, …}



* **Return type**

    dict{key1



* **Raises**

    **ValueError if parsing fails.** –



### _int_from_str(string)
Convert string into integer.


* **Raises**

    **TypeError if string is not a valid integer** –



### _read_nlines(filename: str, n_lines: int)
Read at most nlines lines from file filename.
If nlines is < 0, the entire file is read.


### l2str(l_ang_mom)
Convert the angular momentum l (int) to string.


### str2l(s)
Convert a string to the angular momentum l (int).


### straceback()
Returns a string with the traceback.

## pymatgen.io.abinit.variable module

Support for Abinit input variables.


### _class_ InputVariable(name, value, units='', valperline=3)
Bases: `object`

An Abinit input variable.


* **Parameters**


    * **name** – Name of the variable.


    * **value** – Value of the variable.


    * **units** – String specifying one of the units supported by Abinit. Default: atomic units.


    * **valperline** – Number of items printed per line.



#### _property_ basename()
Return the name trimmed of any dataset index.


#### _property_ dataset()
Return the dataset index in string form.


#### format_list(values, float_decimal=0)
Format a list of values into a string.
The result might be spread among several lines.


#### _static_ format_list2d(values, float_decimal=0)
Format a list of lists.


#### _static_ format_scalar(val, float_decimal=0)
Format a single numerical value into a string
with the appropriate number of decimal.


#### get_value()
Return the value.


#### _property_ name()
Name of the variable.


#### _property_ units()
Return the units.


### flatten(iterable)
Make an iterable flat, i.e. a 1d iterable object.