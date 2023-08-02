---
layout: default
title: pymatgen.io.nwchem.md
nav_exclude: true
---

# pymatgen.io.nwchem module

This module implements input and output processing from Nwchem.

2015/09/21 - Xin Chen ([chenxin13@mails.tsinghua.edu.cn](mailto:chenxin13@mails.tsinghua.edu.cn)):

> NwOutput will read new kinds of data:

> >
> > 1. normal hessian matrix.       [“hessian”]


> > 2. projected hessian matrix.    [“projected_hessian”]


> > 3. normal frequencies.          [“normal_frequencies”]

> For backward compatibility, the key for accessing the projected frequencies
> is still ‘frequencies’.

2015/10/12 - Xin Chen

    NwOutput will read new kinds of data:

    >
    > 1. forces.                      [“forces”]


### _class_ pymatgen.io.nwchem.NwInput(mol, tasks, directives=None, geometry_options=('units', 'angstroms'), symmetry_options=None, memory_options=None)
Bases: `MSONable`

An object representing a Nwchem input file, which is essentially a list
of tasks on a particular molecule.


* **Parameters**


    * **mol** – Input molecule. If molecule is a single string, it is used as a
    direct input to the geometry section of the Gaussian input
    file.


    * **tasks** – List of NwTasks.


    * **directives** – List of root level directives as tuple. E.g.,
    [(“start”, “water”), (“print”, “high”)]


    * **geometry_options** – Additional list of options to be supplied to the
    geometry. E.g., [“units”, “angstroms”, “noautoz”]. Defaults to
    (“units”, “angstroms”).


    * **symmetry_options** – Addition list of option to be supplied to the
    symmetry. E.g. [“c1”] to turn off the symmetry


    * **memory_options** – Memory controlling options. str.
    E.g “total 1000 mb stack 400 mb”.



#### as_dict()
Returns: MSONable dict.


#### _classmethod_ from_dict(d)

* **Parameters**

    **d** (*dict*) – Dict representation.



* **Returns**

    NwInput



#### _classmethod_ from_file(filename)
Read an NwInput from a file. Currently tested to work with
files generated from this class itself.


* **Parameters**

    **filename** – Filename to parse.



* **Returns**

    NwInput object



#### _classmethod_ from_str(string_input)
Read an NwInput from a string. Currently tested to work with
files generated from this class itself.


* **Parameters**

    **string_input** – string_input to parse.



* **Returns**

    NwInput object



#### _classmethod_ from_string(\*args, \*\*kwds)
from_string is deprecated!
Use from_str instead


#### _property_ molecule()
Returns molecule associated with this GaussianInput.


#### write_file(filename)

* **Parameters**

    **filename** (*str*) – Filename.



### _exception_ pymatgen.io.nwchem.NwInputError()
Bases: `Exception`

Error class for NwInput.


### _class_ pymatgen.io.nwchem.NwOutput(filename)
Bases: `object`

A Nwchem output file parser. Very basic for now - supports only dft and
only parses energies and geometries. Please note that Nwchem typically
outputs energies in either au or kJ/mol. All energies are converted to
eV in the parser.


* **Parameters**

    **filename** – Filename to read.



#### get_excitation_spectrum(width=0.1, npoints=2000)
Generate an excitation spectra from the singlet roots of TDDFT
calculations.


* **Parameters**


    * **width** (*float*) – Width for Gaussian smearing.


    * **npoints** (*int*) – Number of energy points. More points => smoother
    curve.



* **Returns**

    (ExcitationSpectrum) which can be plotted using

        pymatgen.vis.plotters.SpectrumPlotter.




#### parse_tddft()
Parses TDDFT roots. Adapted from nw_spectrum.py script.


* **Returns**

    {

        “singlet”: [

            {

                “energy”: float,
                “osc_strength: float

            }

        ],
        “triplet”: [

        > {

        >     “energy”: float

        > }

        ]

    }




### _class_ pymatgen.io.nwchem.NwTask(charge, spin_multiplicity, basis_set, basis_set_option='cartesian', title=None, theory='dft', operation='optimize', theory_directives=None, alternate_directives=None)
Bases: `MSONable`

Base task for Nwchem.

Very flexible arguments to support many types of potential setups.
Users should use more friendly static methods unless they need the
flexibility.


* **Parameters**


    * **charge** – Charge of the molecule. If None, charge on molecule is
    used. Defaults to None. This allows the input file to be set a
    charge independently from the molecule itself.


    * **spin_multiplicity** – Spin multiplicity of molecule. Defaults to None,
    which means that the spin multiplicity is set to 1 if the
    molecule has no unpaired electrons and to 2 if there are
    unpaired electrons.


    * **basis_set** – The basis set used for the task as a dict. E.g.,
    {“C”: “6-311++G\*\*”, “H”: “6-31++G\*\*”}.


    * **basis_set_option** – cartesian (default) | spherical,


    * **title** – Title for the task. Defaults to None, which means a title
    based on the theory and operation of the task is
    autogenerated.


    * **theory** – The theory used for the task. Defaults to “dft”.


    * **operation** – The operation for the task. Defaults to “optimize”.


    * **theory_directives** – A dict of theory directives. For example,
    if you are running dft calculations, you may specify the
    exchange correlation functional using {“xc”: “b3lyp”}.


    * **alternate_directives** – A dict of alternate directives. For
    example, to perform cosmo calculations and dielectric
    constant of 78, you’d supply {‘cosmo’: {“dielectric”: 78}}.



#### as_dict()
Returns: MSONable dict.


#### _classmethod_ dft_task(mol, xc='b3lyp', \*\*kwargs)
A class method for quickly creating DFT tasks with optional
cosmo parameter .


* **Parameters**


    * **mol** – Input molecule


    * **xc** – Exchange correlation to use.


    * **kwargs** – Any of the other kwargs supported by NwTask. Note the
    theory is always “dft” for a dft task.



#### _classmethod_ esp_task(mol, \*\*kwargs)
A class method for quickly creating ESP tasks with RESP
charge fitting.


* **Parameters**


    * **mol** – Input molecule


    * **kwargs** – Any of the other kwargs supported by NwTask. Note the
    theory is always “dft” for a dft task.



#### _classmethod_ from_dict(d)

* **Parameters**

    **d** (*dict*) – Dict representation.



* **Returns**

    NwTask



#### _classmethod_ from_molecule(mol, theory, charge=None, spin_multiplicity=None, basis_set='6-31g', basis_set_option='cartesian', title=None, operation='optimize', theory_directives=None, alternate_directives=None)
Very flexible arguments to support many types of potential setups.
Users should use more friendly static methods unless they need the
flexibility.


* **Parameters**


    * **mol** – Input molecule


    * **charge** – Charge of the molecule. If None, charge on molecule is
    used. Defaults to None. This allows the input file to be set a
    charge independently from the molecule itself.


    * **spin_multiplicity** – Spin multiplicity of molecule. Defaults to None,
    which means that the spin multiplicity is set to 1 if the
    molecule has no unpaired electrons and to 2 if there are
    unpaired electrons.


    * **basis_set** – The basis set to be used as string or a dict. E.g.,
    {“C”: “6-311++G\*\*”, “H”: “6-31++G\*\*”} or “6-31G”. If string,
    same basis set is used for all elements.


    * **basis_set_option** – cartesian (default) | spherical,


    * **title** – Title for the task. Defaults to None, which means a title
    based on the theory and operation of the task is
    autogenerated.


    * **theory** – The theory used for the task. Defaults to “dft”.


    * **operation** – The operation for the task. Defaults to “optimize”.


    * **theory_directives** – A dict of theory directives. For example,
    if you are running dft calculations, you may specify the
    exchange correlation functional using {“xc”: “b3lyp”}.


    * **alternate_directives** – A dict of alternate directives. For
    example, to perform cosmo calculations with DFT, you’d supply
    {‘cosmo’: “cosmo”}.



#### operations(_ = {'': 'dummy', 'dynamics': 'Perform classical molecular dynamics.', 'energy': 'Evaluate the single point energy.', 'freq': 'Same as frequencies.', 'frequencies': 'Compute second derivatives and print out an analysis of molecular vibrations.', 'gradient': 'Evaluate the derivative of the energy with respect to nuclear coordinates.', 'hessian': 'Compute second derivatives.', 'optimize': 'Minimize the energy by varying the molecular structure.', 'property': 'Calculate the properties for the wave function.', 'saddle': 'Conduct a search for a transition state (or saddle point).', 'thermodynamics': 'Perform multi-configuration thermodynamic integration using classical MD.', 'vscf': 'Compute anharmonic contributions to the vibrational modes.'_ )

#### theories(_ = {'band': 'Pseudopotential plane-wave DFT for solids using NWPW', 'ccsd': 'Coupled-cluster single and double excitations', 'ccsd(t)': 'Coupled-cluster linearized triples approximation', 'ccsd+t(ccsd)': 'Fourth order triples contribution', 'dft': 'DFT', 'direct_mp2': 'MP2 using a full-direct algorithm', 'esp': 'ESP', 'g3gn': 'some description', 'mcscf': 'Multiconfiguration SCF', 'md': 'Classical molecular dynamics simulation', 'mp2': 'MP2 using a semi-direct algorithm', 'pspw': 'Pseudopotential plane-wave DFT for molecules and insulating solids using NWPW', 'rimp2': 'MP2 using the RI approximation', 'scf': 'Hartree-Fock', 'selci': 'Selected CI with perturbation correction', 'sodft': 'Spin-Orbit DFT', 'tce': 'Tensor Contraction Engine', 'tddft': 'Time Dependent DFT'_ )