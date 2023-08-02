---
layout: default
title: pymatgen.io.abinit.abiobjects.md
nav_exclude: true
---

# pymatgen.io.abinit.abiobjects module

Low-level objects providing an abstraction for the objects involved in the calculation.


### _class_ pymatgen.io.abinit.abiobjects.AbivarAble()
Bases: `object`

An AbivarAble object provides a method to_abivars
that returns a dictionary with the abinit variables.


#### _abstract_ to_abivars()
Returns a dictionary with the abinit variables.


### _class_ pymatgen.io.abinit.abiobjects.Constraints()
Bases: `AbivarAble`

This object defines the constraints for structural relaxation.


#### to_abivars()
Dictionary with Abinit variables.


### _class_ pymatgen.io.abinit.abiobjects.Electrons(spin_mode='polarized', smearing='fermi_dirac:0.1 eV', algorithm=None, nband=None, fband=None, charge=0.0, comment=None)
Bases: `AbivarAble`, `MSONable`

The electronic degrees of freedom.

Constructor for Electrons object.


* **Parameters**


    * **comment** – String comment for Electrons


    * **charge** – Total charge of the system. Default is 0.



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


### _class_ pymatgen.io.abinit.abiobjects.ElectronsAlgorithm(\*args, \*\*kwargs)
Bases: `dict`, `AbivarAble`, `MSONable`

Variables controlling the SCF/NSCF algorithm.

Initialize object.


#### as_dict()
Convert object to dict.


#### _classmethod_ from_dict(d)
Build object from dict.


#### to_abivars()
Dictionary with Abinit input variables.


### _class_ pymatgen.io.abinit.abiobjects.ExcHamiltonian(bs_loband, nband, mbpt_sciss, coulomb_mode, ecuteps, spin_mode='polarized', mdf_epsinf=None, exc_type='TDA', algo='haydock', with_lf=True, bs_freq_mesh=None, zcut=None, \*\*kwargs)
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


### _class_ pymatgen.io.abinit.abiobjects.HilbertTransform(nomegasf, domegasf=None, spmeth=1, nfreqre=None, freqremax=None, nfreqim=None, freqremin=None)
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



#### to_abivars()
Returns a dictionary with the abinit variables.


### _class_ pymatgen.io.abinit.abiobjects.KSampling(mode=KSamplingModes.monkhorst, num_kpts=0, kpts=((1, 1, 1),), kpt_shifts=(0.5, 0.5, 0.5), kpts_weights=None, use_symmetries=True, use_time_reversal=True, chksymbreak=None, comment=None)
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

    `KSampling` object.



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

    `KSampling` object.



#### _classmethod_ monkhorst_automatic(structure, ngkpt, use_symmetries=True, use_time_reversal=True, chksymbreak=None, comment=None)
Convenient static constructor for an automatic Monkhorst-Pack mesh.


* **Parameters**


    * **structure** – `Structure` object.


    * **ngkpt** – Subdivisions N_1, N_2 and N_3 along reciprocal lattice vectors.


    * **use_symmetries** – Use spatial symmetries to reduce the number of k-points.


    * **use_time_reversal** – Use time-reversal symmetry to reduce the number of k-points.



* **Returns**

    `KSampling` object.



#### _classmethod_ path_from_structure(ndivsm, structure)
See _path for the meaning of the variables.


#### to_abivars()
Dictionary with Abinit variables.


### _class_ pymatgen.io.abinit.abiobjects.KSamplingModes(value)
Bases: `Enum`

Enum if the different samplings of the BZ.


#### automatic(_ = _ )

#### monkhorst(_ = _ )

#### path(_ = _ )

### _class_ pymatgen.io.abinit.abiobjects.ModelDielectricFunction(mdf_epsinf)
Bases: `AbivarAble`

Model dielectric function used for BSE calculation.


* **Parameters**

    **mdf_epsinf** – Value of epsilon_infinity.



#### to_abivars()
Return dictionary with abinit variables.


### _class_ pymatgen.io.abinit.abiobjects.PPModel(mode='godby', plasmon_freq=None)
Bases: `AbivarAble`, `MSONable`

Parameters defining the plasmon-pole technique.
The common way to instantiate a PPModel object is via the class method PPModel.as_ppmodel(string).


* **Parameters**


    * **mode** – ppmodel type


    * **plasmon_freq** – Plasmon frequency in Ha.



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


### _class_ pymatgen.io.abinit.abiobjects.PPModelModes(value)
Bases: `Enum`

Different kind of plasmon-pole models.


#### farid(_ = _ )

#### godby(_ = _ )

#### hybersten(_ = _ )

#### linden(_ = _ )

#### noppmodel(_ = _ )

### _class_ pymatgen.io.abinit.abiobjects.RelaxationMethod(\*args, \*\*kwargs)
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


### _class_ pymatgen.io.abinit.abiobjects.Screening(ecuteps, nband, w_type='RPA', sc_mode='one_shot', hilbert=None, ecutwfn=None, inclvkb=2)
Bases: `AbivarAble`

This object defines the parameters used for the
computation of the screening function.


* **Parameters**


    * **ecuteps** – Cutoff energy for the screening (Ha units).


    * **function** (*nband Number** of **bands for the Green's*) –


    * **w_type** – Screening type


    * **sc_mode** – Self-consistency mode.


    * **hilbert** – Instance of `HilbertTransform` defining the parameters for the Hilber transform method.


    * **ecutwfn** – Cutoff energy for the wavefunctions (Default: ecutwfn == ecut).


    * **inclvkb** – Option for the treatment of the dipole matrix elements (NC pseudos).



#### to_abivars()
Returns a dictionary with the abinit variables.


#### _property_ use_hilbert()
True if we are using the Hilbert transform method.


### _class_ pymatgen.io.abinit.abiobjects.SelfEnergy(se_type, sc_mode, nband, ecutsigx, screening, gw_qprange=1, ppmodel=None, ecuteps=None, ecutwfn=None, gwpara=2)
Bases: `AbivarAble`

This object defines the parameters used for the computation of the self-energy.


* **Parameters**


    * **se_type** – Type of self-energy (str)


    * **sc_mode** – Self-consistency mode.


    * **nband** – Number of bands for the Green’s function


    * **ecutsigx** – Cutoff energy for the exchange part of the self-energy (Ha units).


    * **screening** – `Screening` instance.


    * **gw_qprange** – Option for the automatic selection of k-points and bands for GW corrections.
    See Abinit docs for more detail. The default value makes the code computie the
    QP energies for all the point in the IBZ and one band above and one band below the Fermi level.


    * **ppmodel** – `PPModel` instance with the parameters used for the plasmon-pole technique.


    * **ecuteps** – Cutoff energy for the screening (Ha units).


    * **ecutwfn** – Cutoff energy for the wavefunctions (Default: ecutwfn == ecut).



#### _property_ gwcalctyp()
Returns the value of the gwcalctyp input variable.


#### _property_ symsigma()
1 if symmetries can be used to reduce the number of q-points.


#### to_abivars()
Returns a dictionary with the abinit variables.


#### _property_ use_ppmodel()
True if we are using the plasmon-pole approximation.


### _class_ pymatgen.io.abinit.abiobjects.Smearing(occopt, tsmear)
Bases: `AbivarAble`, `MSONable`

Variables defining the smearing technique. The preferred way to instantiate
a Smearing object is via the class method Smearing.as_smearing(string).

Build object with occopt and tsmear.


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


### _class_ pymatgen.io.abinit.abiobjects.SpinMode(mode, nsppol, nspinor, nspden)
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


#### as_dict()
Convert object to dict.


#### _classmethod_ as_spinmode(obj)
Converts obj into a SpinMode instance.


#### _classmethod_ from_dict(d)
Build object from dict.


#### to_abivars()
Dictionary with Abinit input variables.


### pymatgen.io.abinit.abiobjects.contract(s)
```python
>>> assert contract("1 1 1 2 2 3") == "3*1 2*2 1*3"
>>> assert contract("1 1 3 2 3") == "2*1 1*3 1*2 1*3".
```


### pymatgen.io.abinit.abiobjects.lattice_from_abivars(cls=None, \*args, \*\*kwargs)
Returns a Lattice object from a dictionary
with the Abinit variables acell and either rprim in Bohr or angdeg
If acell is not given, the Abinit default is used i.e. [1,1,1] Bohr.


* **Parameters**

    **cls** – Lattice class to be instantiated. pymatgen.core.lattice.Lattice if cls is None


### Example

lattice_from_abivars(acell=3\*[10], rprim=np.eye(3))


### pymatgen.io.abinit.abiobjects.species_by_znucl(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))
Return list of unique specie found in structure **ordered according to sites**.

### Example

Site0: 0.5 0 0 O
Site1: 0   0 0 Si

produces [Specie_O, Specie_Si] and not set([Specie_O, Specie_Si]) as in types_of_specie


### pymatgen.io.abinit.abiobjects.structure_from_abivars(cls=None, \*args, \*\*kwargs)
Build a `Structure` object from a dictionary with ABINIT variables.


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


### pymatgen.io.abinit.abiobjects.structure_to_abivars(structure, enforce_znucl=None, enforce_typat=None, \*\*kwargs)
Receives a structure and returns a dictionary with ABINIT variables.


* **Parameters**


    * **enforce_znucl** – List of ntypat entries with the value of Z for each type of atom.
    Used to change the default ordering.


    * **enforce_typat** – List with natom entries with the type index.
    Fortran conventions: start to count from 1.
    Used to change the default ordering.