---
layout: default
title: pymatgen.io.babel.md
nav_exclude: true
---

# pymatgen.io.babel module

OpenBabel interface module, which opens up access to the hundreds of file
formats supported by OpenBabel. Requires openbabel with python bindings to be
installed. Please consult the
[openbabel documentation](http://openbabel.org/wiki/Main_Page).


### _class_ pymatgen.io.babel.BabelMolAdaptor(mol)
Bases: `object`

Adaptor serves as a bridge between OpenBabel’s Molecule and pymatgen’s
Molecule.

Initializes with pymatgen Molecule or OpenBabel’s OBMol.


* **Parameters**

    **mol** – pymatgen’s Molecule/IMolecule or OpenBabel OBMol



#### add_hydrogen()
Add hydrogens (make all hydrogen explicit).


#### confab_conformers(forcefield='mmff94', freeze_atoms=None, rmsd_cutoff=0.5, energy_cutoff=50.0, conf_cutoff=100000, verbose=False)
Conformer generation based on Confab to generate all diverse low-energy
conformers for molecules. This is different from rotor_conformer or
gen3d_conformer as it aims to not simply to find a low energy
conformation but to generate several different conformations.


* **Parameters**


    * **forcefield** (*str*) – Default is mmff94. Options are ‘gaff’, ‘ghemical’,
    ‘mmff94’, ‘mmff94s’, and ‘uff’.


    * **freeze_atoms** (*[**int**]*) – index of atoms to be freezed when performing
    conformer search, default is None.


    * **rmsd_cutoff** (*float*) – rmsd_cufoff, default is 0.5 Angstrom.


    * **energy_cutoff** (*float*) – energy_cutoff, default is 50.0 kcal/mol.


    * **conf_cutoff** (*float*) – max number of conformers to test,
    default is 1 million.


    * **verbose** (*bool*) – whether to display information on torsions found,
    default is False.



* **Returns**

    list of pymatgen Molecule objects for generated conformers.



* **Return type**

    (list)



#### _static_ from_file(filename, file_format='xyz', return_all_molecules=False)
Uses OpenBabel to read a molecule from a file in all supported formats.


* **Parameters**


    * **filename** – Filename of input file


    * **file_format** – String specifying any OpenBabel supported formats.


    * **return_all_molecules** – If `True`, will return a list of
    `BabelMolAdaptor` instances, one for each molecule found in
    the file. If `False`, will return only the first molecule.



* **Returns**

    BabelMolAdaptor object or list thereof



#### _static_ from_molecule_graph(mol)
Read a molecule from a pymatgen MoleculeGraph object.


* **Parameters**

    **mol** – pymatgen MoleculeGraph object.



* **Returns**

    BabelMolAdaptor object



#### _static_ from_str(string_data, file_format='xyz')
Uses OpenBabel to read a molecule from a string in all supported
formats.


* **Parameters**


    * **string_data** – String containing molecule data.


    * **file_format** – String specifying any OpenBabel supported formats.



* **Returns**

    BabelMolAdaptor object



#### _classmethod_ from_string(\*args, \*\*kwds)
from_string is deprecated!
Use from_str instead


#### gen3d_conformer()
A combined method to first generate 3D structures from 0D or 2D
structures and then find the minimum energy conformer:


1. Use OBBuilder to create a 3D structure using rules and ring templates


2. Do 250 steps of a steepest descent geometry optimization with the
MMFF94 forcefield


3. Do 200 iterations of a Weighted Rotor conformational search
(optimizing each conformer with 25 steps of a steepest descent)


4. Do 250 steps of a conjugate gradient geometry optimization.

Warning from openbabel docs:
For many applications where 100s if not 1000s of molecules need to be
processed, gen3d is rather SLOW. Sometimes this function can cause a
segmentation fault.
A future version of Open Babel will provide options for slow/medium/fast
3D structure generation which will involve different compromises
between speed and finding the global energy minimum.


#### localopt(forcefield='mmff94', steps=500)
A wrapper to pybel’s localopt method to optimize a Molecule.


* **Parameters**


    * **forcefield** – Default is mmff94. Options are ‘gaff’, ‘ghemical’,
    ‘mmff94’, ‘mmff94s’, and ‘uff’.


    * **steps** – Default is 500.



#### make3d(forcefield='mmff94', steps=50)
A wrapper to pybel’s make3D method generate a 3D structure from a
2D or 0D structure.
The 3D structure is made very quickly using a combination of rules
(e.g. sp3 atoms should have four bonds arranged in a tetrahedron) and
ring templates (e.g. cyclohexane is shaped like a chair). Once 3D
coordinates are generated, hydrogens are added and a quick local
optimization is carried out as default.

The generated 3D structure can have clashes or have high energy
structures due to some strain. Please consider to use the conformer
search or geometry optimization to further optimize the structure.


* **Parameters**


    * **forcefield** – Default is mmff94. Options are ‘gaff’, ‘ghemical’,
    ‘mmff94’, ‘mmff94s’, and ‘uff’.


    * **steps** – Default is 50.



#### _property_ openbabel_mol()
Returns OpenBabel’s OBMol.


#### _property_ pybel_mol()
Returns Pybel’s Molecule object.


#### _property_ pymatgen_mol()
Returns pymatgen Molecule object.


#### remove_bond(idx1, idx2)
Remove a bond from an openbabel molecule.


* **Parameters**


    * **idx1** – The atom index of one of the atoms participating the in bond


    * **idx2** – The atom index of the other atom participating in the bond



#### rotor_conformer(\*rotor_args, algo='WeightedRotorSearch', forcefield='mmff94')
Conformer search based on several Rotor Search algorithms of openbabel.
If the input molecule is not 3D, make3d will be called (generate 3D
structure, add hydrogen, a quick localopt). All hydrogen atoms need
to be made explicit.


* **Parameters**


    * **rotor_args** – pass args to Rotor Search in openbabel.
    for “WeightedRotorSearch”: (conformers, geomSteps,
    sampleRingBonds-default False)
    for “SystematicRotorSearch”: (geomSteps-default 2500,
    sampleRingBonds-default False)
    for “RandomRotorSearch”: (conformers, geomSteps-default 2500,
    sampleRingBonds-default False)


    * **algo** (*str*) – Default is “WeightedRotorSearch”. Options are
    “SystematicRotorSearch”, “RandomRotorSearch”, and
    “WeightedRotorSearch”.


    * **forcefield** (*str*) – Default is mmff94. Options are ‘gaff’, ‘ghemical’,
    ‘mmff94’, ‘mmff94s’, and ‘uff’.



#### write_file(filename, file_format='xyz')
Uses OpenBabel to output all supported formats.


* **Parameters**


    * **filename** – Filename of file to output


    * **file_format** – String specifying any OpenBabel supported formats.