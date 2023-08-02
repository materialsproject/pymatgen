---
layout: default
title: pymatgen.analysis.molecule_structure_comparator.md
nav_exclude: true
---

# pymatgen.analysis.molecule_structure_comparator module

This module provides classes to comparison the structures of the two
molecule. As long as the two molecule have the same bond connection tables,
the molecules are deemed to be same. The atom in the two molecule must be
paired accordingly.
This module is supposed to perform rough comparisons with the atom order
correspondence prerequisite, while molecule_matcher is supposed to do exact
comparisons without the atom order correspondence prerequisite.


### _class_ pymatgen.analysis.molecule_structure_comparator.CovalentRadius()
Bases: `object`

Covalent radius of the elements.

Beatriz C. et al. Dalton Trans. 2008, 2832-2838. [https://doi.org/10.1039/b801115j](https://doi.org/10.1039/b801115j)


#### radius(_ = {'Ac': 2.15, 'Ag': 1.45, 'Al': 1.21, 'Am': 1.8, 'Ar': 1.06, 'As': 1.19, 'At': 1.5, 'Au': 1.36, 'B': 0.84, 'Ba': 2.15, 'Be': 0.96, 'Bi': 1.48, 'Br': 1.2, 'C': 0.73, 'Ca': 1.76, 'Cd': 1.44, 'Ce': 2.04, 'Cl': 1.02, 'Cm': 1.69, 'Co': 1.38, 'Cr': 1.39, 'Cs': 2.44, 'Cu': 1.32, 'Dy': 1.92, 'Er': 1.89, 'Eu': 1.98, 'F': 0.57, 'Fe': 1.42, 'Fr': 2.6, 'Ga': 1.22, 'Gd': 1.96, 'Ge': 1.2, 'H': 0.31, 'He': 0.28, 'Hf': 1.75, 'Hg': 1.32, 'Ho': 1.92, 'I': 1.39, 'In': 1.42, 'Ir': 1.41, 'K': 2.03, 'Kr': 1.16, 'La': 2.07, 'Li': 1.28, 'Lu': 1.87, 'Mg': 1.41, 'Mn': 1.5, 'Mo': 1.54, 'N': 0.71, 'Na': 1.66, 'Nb': 1.64, 'Nd': 2.01, 'Ne': 0.58, 'Ni': 1.24, 'Np': 1.9, 'O': 0.66, 'Os': 1.44, 'P': 1.07, 'Pa': 2, 'Pb': 1.46, 'Pd': 1.39, 'Pm': 1.99, 'Po': 1.4, 'Pr': 2.03, 'Pt': 1.36, 'Pu': 1.87, 'Ra': 2.21, 'Rb': 2.2, 'Re': 1.51, 'Rh': 1.42, 'Rn': 1.5, 'Ru': 1.46, 'S': 1.05, 'Sb': 1.39, 'Sc': 1.7, 'Se': 1.2, 'Si': 1.11, 'Sm': 1.98, 'Sn': 1.39, 'Sr': 1.95, 'Ta': 1.7, 'Tb': 1.94, 'Tc': 1.47, 'Te': 1.38, 'Th': 2.06, 'Ti': 1.6, 'Tl': 1.45, 'Tm': 1.9, 'U': 1.96, 'V': 1.53, 'W': 1.62, 'Xe': 1.4, 'Y': 1.9, 'Yb': 1.87, 'Zn': 1.22, 'Zr': 1.75_ )

### _class_ pymatgen.analysis.molecule_structure_comparator.MoleculeStructureComparator(bond_length_cap=0.3, covalent_radius={'Ac': 2.15, 'Ag': 1.45, 'Al': 1.21, 'Am': 1.8, 'Ar': 1.06, 'As': 1.19, 'At': 1.5, 'Au': 1.36, 'B': 0.84, 'Ba': 2.15, 'Be': 0.96, 'Bi': 1.48, 'Br': 1.2, 'C': 0.73, 'Ca': 1.76, 'Cd': 1.44, 'Ce': 2.04, 'Cl': 1.02, 'Cm': 1.69, 'Co': 1.38, 'Cr': 1.39, 'Cs': 2.44, 'Cu': 1.32, 'Dy': 1.92, 'Er': 1.89, 'Eu': 1.98, 'F': 0.57, 'Fe': 1.42, 'Fr': 2.6, 'Ga': 1.22, 'Gd': 1.96, 'Ge': 1.2, 'H': 0.31, 'He': 0.28, 'Hf': 1.75, 'Hg': 1.32, 'Ho': 1.92, 'I': 1.39, 'In': 1.42, 'Ir': 1.41, 'K': 2.03, 'Kr': 1.16, 'La': 2.07, 'Li': 1.28, 'Lu': 1.87, 'Mg': 1.41, 'Mn': 1.5, 'Mo': 1.54, 'N': 0.71, 'Na': 1.66, 'Nb': 1.64, 'Nd': 2.01, 'Ne': 0.58, 'Ni': 1.24, 'Np': 1.9, 'O': 0.66, 'Os': 1.44, 'P': 1.07, 'Pa': 2, 'Pb': 1.46, 'Pd': 1.39, 'Pm': 1.99, 'Po': 1.4, 'Pr': 2.03, 'Pt': 1.36, 'Pu': 1.87, 'Ra': 2.21, 'Rb': 2.2, 'Re': 1.51, 'Rh': 1.42, 'Rn': 1.5, 'Ru': 1.46, 'S': 1.05, 'Sb': 1.39, 'Sc': 1.7, 'Se': 1.2, 'Si': 1.11, 'Sm': 1.98, 'Sn': 1.39, 'Sr': 1.95, 'Ta': 1.7, 'Tb': 1.94, 'Tc': 1.47, 'Te': 1.38, 'Th': 2.06, 'Ti': 1.6, 'Tl': 1.45, 'Tm': 1.9, 'U': 1.96, 'V': 1.53, 'W': 1.62, 'Xe': 1.4, 'Y': 1.9, 'Yb': 1.87, 'Zn': 1.22, 'Zr': 1.75}, priority_bonds=(), priority_cap=0.8, ignore_ionic_bond=True, bond_13_cap=0.05)
Bases: `MSONable`

Class to check whether the connection tables of the two molecules are the
same. The atom in the two molecule must be paired accordingly.


* **Parameters**


    * **bond_length_cap** – The ratio of the elongation of the bond to be
    acknowledged. If the distance between two atoms is less than (
    empirical covalent bond length) X (1 + bond_length_cap), the bond
    between the two atoms will be acknowledged.


    * **covalent_radius** – The covalent radius of the atoms.
    dict (element symbol -> radius)


    * **priority_bonds** – The bonds that are known to be existed in the initial
    molecule. Such bonds will be acknowledged in a loose criteria.
    The index should start from 0.


    * **priority_cap** – The ratio of the elongation of the bond to be
    acknowledged for the priority bonds.



#### are_equal(mol1, mol2)
Compare the bond table of the two molecules.


* **Parameters**


    * **mol1** – first molecule. pymatgen Molecule object.


    * **mol2** – second molecules. pymatgen Molecule object.



#### as_dict()
Returns: MSONable dict.


#### _classmethod_ from_dict(dct)

* **Parameters**

    **d** (*dict*) – Dict representation.



* **Returns**

    MoleculeStructureComparator



#### _static_ get_13_bonds(priority_bonds)

* **Parameters**

    **(****)** (*priority_bonds*) –


Returns:


#### halogen_list(_ = ('F', 'Cl', 'Br', 'I'_ )

#### ionic_element_list(_ = ('Na', 'Mg', 'Al', 'Sc', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Rb', 'Sr'_ )