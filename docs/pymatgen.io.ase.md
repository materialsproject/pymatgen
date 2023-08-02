---
layout: default
title: pymatgen.io.ase.md
nav_exclude: true
---

# pymatgen.io.ase module

This module provides conversion between the Atomic Simulation Environment
Atoms object and pymatgen Structure objects.


### _class_ pymatgen.io.ase.AseAtomsAdaptor()
Bases: `object`

Adaptor serves as a bridge between ASE Atoms and pymatgen objects.


#### _static_ get_atoms(structure: [SiteCollection](pymatgen.core.structure.md#pymatgen.core.structure.SiteCollection), \*\*kwargs)
Returns ASE Atoms object from pymatgen structure or molecule.


* **Parameters**


    * **structure** ([*SiteCollection*](pymatgen.core.structure.md#pymatgen.core.structure.SiteCollection)) – pymatgen Structure or Molecule


    * **\*\*kwargs** – passed to the ASE Atoms constructor



* **Returns**

    ASE Atoms object



* **Return type**

    [Atoms](pymatgen.io.feff.inputs.md#pymatgen.io.feff.inputs.Atoms)



#### _static_ get_molecule(atoms: ~ase.atoms.Atoms, cls: type[pymatgen.core.structure.Molecule] = <class 'pymatgen.core.structure.Molecule'>, \*\*cls_kwargs)
Returns pymatgen molecule from ASE Atoms.


* **Parameters**


    * **atoms** – ASE Atoms object


    * **cls** – The Molecule class to instantiate (defaults to pymatgen molecule)


    * **\*\*cls_kwargs** – Any additional kwargs to pass to the cls



* **Returns**

    Equivalent pymatgen.core.structure.Molecule



* **Return type**

    [Molecule](pymatgen.core.structure.md#pymatgen.core.structure.Molecule)



#### _static_ get_structure(atoms: ~ase.atoms.Atoms, cls: type[pymatgen.core.structure.Structure] = <class 'pymatgen.core.structure.Structure'>, \*\*cls_kwargs)
Returns pymatgen structure from ASE Atoms.


* **Parameters**


    * **atoms** – ASE Atoms object


    * **cls** – The Structure class to instantiate (defaults to pymatgen Structure)


    * **\*\*cls_kwargs** – Any additional kwargs to pass to the cls



* **Returns**

    Equivalent pymatgen.core.structure.Structure