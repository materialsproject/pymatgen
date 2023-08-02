---
layout: default
title: pymatgen.io.abinit.inputs.md
nav_exclude: true
---

# pymatgen.io.abinit.inputs module

This module defines a simplified interface for generating ABINIT input files.
Note that not all the features of Abinit are supported by BasicAbinitInput.
For a more comprehensive implementation, use the AbinitInput object provided by AbiPy.


### _class_ pymatgen.io.abinit.inputs.AbstractInput()
Bases: `MutableMapping`

Abstract class defining the methods that must be implemented by Input objects.


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
Write the input file to file to `filepath`.


### _class_ pymatgen.io.abinit.inputs.BasicAbinitInput(structure, pseudos, pseudo_dir=None, comment=None, abi_args=None, abi_kwargs=None)
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


#### add_abiobjects(\*abi_objects)
This function receive a list of `AbiVarable` objects and add
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


    * **spin_mode** – `SpinMode` object or string. Possible values for string are:


    * **polarized** (*-*) –


    * **unpolarized** (*-*) –


    * **afm** (*-*) –


    * **spinor** (*-*) –


    * **spinor_nomag** (*-*) –



#### set_structure(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))
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


### _exception_ pymatgen.io.abinit.inputs.BasicAbinitInputError()
Bases: `Exception`

Base error class for exceptions raised by `BasicAbinitInput`.


### _class_ pymatgen.io.abinit.inputs.BasicMultiDataset(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), pseudos, pseudo_dir='', ndtset=1)
Bases: `object`

This object is essentially a list of BasicAbinitInput objects.
that provides an easy-to-use interface to apply global changes to the
the inputs stored in the objects.

Let’s assume for example that multi contains two `BasicAbinitInput` objects and we
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
Add a new entry in the multidataset by copying the input with index `dtindex`.


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
Write `ndset` input files to disk. The name of the file
is constructed from the dataset index e.g. run0.abi.


### _class_ pymatgen.io.abinit.inputs.ShiftMode(value)
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
on the iniital letter of the string. G for GammaCenterd, M for MonkhorstPack,
S for Symmetric, O for OneSymmetric.
Case insensitive.


### pymatgen.io.abinit.inputs.as_structure(obj)
Convert obj into a Structure. Accepts:

>
> * Structure object.


> * Filename


> * Dictionaries (MSONable format or dictionaries with abinit variables).


### pymatgen.io.abinit.inputs.calc_shiftk(structure, symprec: float = 0.01, angle_tolerance=5)
Find the values of `shiftk` and `nshiftk` appropriated for the sampling of the Brillouin zone.

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



### pymatgen.io.abinit.inputs.ebands_input(structure, pseudos, kppa=None, nscf_nband=None, ndivsm=15, ecut=None, pawecutdg=None, scf_nband=None, accuracy='normal', spin_mode='polarized', smearing='fermi_dirac:0.1 eV', charge=0.0, scf_algorithm=None, dos_kppa=None)
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



### pymatgen.io.abinit.inputs.gs_input(structure, pseudos, kppa=None, ecut=None, pawecutdg=None, scf_nband=None, accuracy='normal', spin_mode='polarized', smearing='fermi_dirac:0.1 eV', charge=0.0, scf_algorithm=None)
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



### pymatgen.io.abinit.inputs.ion_ioncell_relax_input(structure, pseudos, kppa=None, nband=None, ecut=None, pawecutdg=None, accuracy='normal', spin_mode='polarized', smearing='fermi_dirac:0.1 eV', charge=0.0, scf_algorithm=None, shift_mode='Monkhorst-pack')
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



### pymatgen.io.abinit.inputs.num_valence_electrons(structure, pseudos)
Returns the number of valence electrons.


* **Parameters**

    **pseudos** – List of

    ```
    |Pseudo|
    ```

     objects or list of filenames.