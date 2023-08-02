---
layout: default
title: pymatgen.io.phonopy.md
nav_exclude: true
---

# pymatgen.io.phonopy module

Module for interfacing with phonopy, see [https://atztogo.github.io/phonopy/](https://atztogo.github.io/phonopy/).


### pymatgen.io.phonopy.eigvec_to_eigdispl(v, q, frac_coords, mass)
Converts a single eigenvector to an eigendisplacement in the primitive cell
according to the formula:

```default
exp(2*pi*i*(frac_coords \\dot q) / sqrt(mass) * v
```

Compared to the modulation option in phonopy, here all the additional
multiplicative and phase factors are set to 1.


* **Parameters**


    * **v** – the vector that should be converted. A 3D complex numpy array.


    * **q** – the q point in fractional coordinates


    * **frac_coords** – the fractional coordinates of the atom


    * **mass** – the mass of the atom



### pymatgen.io.phonopy.get_complete_ph_dos(partial_dos_path, phonopy_yaml_path)
Creates a pymatgen CompletePhononDos from a partial_dos.dat and
phonopy.yaml files.
The second is produced when generating a Dos and is needed to extract
the structure.


* **Parameters**


    * **partial_dos_path** – path to the partial_dos.dat file.


    * **phonopy_yaml_path** – path to the phonopy.yaml file.



### pymatgen.io.phonopy.get_displaced_structures(pmg_structure, atom_disp=0.01, supercell_matrix=None, yaml_fname=None, \*\*kwargs)
Generate a set of symmetrically inequivalent displaced structures for
phonon calculations.


* **Parameters**


    * **pmg_structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – A pymatgen structure object.


    * **atom_disp** (*float*) – Atomic displacement. Default is 0.01 $\\AA$.


    * **supercell_matrix** (*3x3 array*) – Scaling matrix for supercell.


    * **yaml_fname** (*str*) – If not None, it represents the full path to
    the outputting displacement yaml file, e.g. disp.yaml.


    * **\*\*kwargs** – Parameters used in Phonopy.generate_displacement method.



* **Returns**

    A list of symmetrically inequivalent structures with displacements, in
    which the first element is the perfect supercell structure.



### pymatgen.io.phonopy.get_gruneisen_ph_bs_symm_line(gruneisen_path, structure=None, structure_path=None, labels_dict=None, fit=False)
Creates a pymatgen GruneisenPhononBandStructure from a band.yaml file.
The labels will be extracted from the dictionary, if present.
If the ‘eigenvector’ key is found the eigendisplacements will be
calculated according to the formula:
\\exp(2\*pi\*i\*(frac_coords \\dot q) / sqrt(mass) \* v

> and added to the object.


* **Parameters**


    * **gruneisen_path** – path to the band.yaml file


    * **structure** – pymaten Structure object


    * **structure_path** – path to a structure file (e.g., POSCAR)


    * **labels_dict** – dict that links a qpoint in frac coords to a label.


    * **fit** – Substitute Grueneisen parameters close to the gamma point
    with points obtained from a fit to a spline if the derivate from
    a smooth curve (i.e. if the slope changes by more than 200% in the
    range of 10% around the gamma point).
    These derivations occur because of very small frequencies
    (and therefore numerical inaccuracies) close to gamma.



### pymatgen.io.phonopy.get_gruneisenparameter(gruneisen_path, structure=None, structure_path=None)
Get Gruneisen object from gruneisen.yaml file, as obtained from phonopy (Frequencies in THz!).
The order is structure > structure path > structure from gruneisen dict.
Newer versions of phonopy include the structure in the yaml file,
the structure/structure_path is kept for compatibility.


* **Parameters**


    * **gruneisen_path** – Path to gruneisen.yaml file (frequencies have to be in THz!)


    * **structure** – pymatgen Structure object


    * **structure_path** – path to structure in a file (e.g., POSCAR)


Returns: GruneisenParameter object


### pymatgen.io.phonopy.get_gs_ph_bs_symm_line_from_dict(gruneisen_dict, structure=None, structure_path=None, labels_dict=None, fit=False)
Creates a pymatgen GruneisenPhononBandStructure object from the dictionary
extracted by the gruneisen.yaml file produced by phonopy. The labels
will be extracted from the dictionary, if present. If the ‘eigenvector’
key is found the eigendisplacements will be calculated according to the
formula:

```default
exp(2*pi*i*(frac_coords \\dot q) / sqrt(mass) * v
```

and added to the object. A fit algorithm can be used to replace diverging
Gruneisen values close to gamma.


* **Parameters**


    * **gruneisen_dict** (*dict*) – the dictionary extracted from the gruneisen.yaml file


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – pymatgen structure object


    * **structure_path** – path to structure file


    * **labels_dict** (*dict*) – dict that links a qpoint in frac coords to a label.
    Its value will replace the data contained in the band.yaml.


    * **fit** (*bool*) – Substitute Grueneisen parameters close to the gamma point
    with points obtained from a fit to a spline if the derivate from
    a smooth curve (i.e. if the slope changes by more than 200% in the
    range of 10% around the gamma point).
    These derivations occur because of very small frequencies
    (and therefore numerical inaccuracies) close to gamma.



### pymatgen.io.phonopy.get_ph_bs_symm_line(bands_path, has_nac=False, labels_dict=None)
Creates a pymatgen PhononBandStructure from a band.yaml file.
The labels will be extracted from the dictionary, if present.
If the ‘eigenvector’  key is found the eigendisplacements will be
calculated according to the formula:
\\exp(2\*pi\*i\*(frac_coords \\dot q) / sqrt(mass) \* v

> and added to the object.


* **Parameters**


    * **bands_path** – path to the band.yaml file


    * **has_nac** – True if the data have been obtained with the option
    –nac option. Default False.


    * **labels_dict** – dict that links a qpoint in frac coords to a label.



### pymatgen.io.phonopy.get_ph_bs_symm_line_from_dict(bands_dict, has_nac=False, labels_dict=None)
Creates a pymatgen PhononBandStructure object from the dictionary
extracted by the band.yaml file produced by phonopy. The labels
will be extracted from the dictionary, if present. If the ‘eigenvector’
key is found the eigendisplacements will be calculated according to the
formula:

```default
exp(2*pi*i*(frac_coords \\dot q) / sqrt(mass) * v
```

and added to the object.


* **Parameters**


    * **bands_dict** – the dictionary extracted from the band.yaml file


    * **has_nac** – True if the data have been obtained with the option
    –nac option. Default False.


    * **labels_dict** – dict that links a qpoint in frac coords to a label.
    Its value will replace the data contained in the band.yaml.



### pymatgen.io.phonopy.get_ph_dos(total_dos_path)
Creates a pymatgen PhononDos from a total_dos.dat file.


* **Parameters**

    **total_dos_path** – path to the total_dos.dat file.



### pymatgen.io.phonopy.get_phonon_band_structure_from_fc(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), supercell_matrix: ndarray, force_constants: ndarray, mesh_density: float = 100.0, \*\*kwargs)
Get a uniform phonon band structure from phonopy force constants.


* **Parameters**


    * **structure** – A structure.


    * **supercell_matrix** – The supercell matrix used to generate the force
    constants.


    * **force_constants** – The force constants in phonopy format.


    * **mesh_density** – The density of the q-point mesh. See the docstring
    for the `mesh` argument in Phonopy.init_mesh() for more details.


    * **\*\*kwargs** – Additional kwargs passed to the Phonopy constructor.



* **Returns**

    The uniform phonon band structure.



### pymatgen.io.phonopy.get_phonon_band_structure_symm_line_from_fc(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), supercell_matrix: ndarray, force_constants: ndarray, line_density: float = 20.0, symprec: float = 0.01, \*\*kwargs)
Get a phonon band structure along a high symmetry path from phonopy force
constants.


* **Parameters**


    * **structure** – A structure.


    * **supercell_matrix** – The supercell matrix used to generate the force
    constants.


    * **force_constants** – The force constants in phonopy format.


    * **line_density** – The density along the high symmetry path.


    * **symprec** – Symmetry precision passed to phonopy and used for determining
    the band structure path.


    * **\*\*kwargs** – Additional kwargs passed to the Phonopy constructor.



* **Returns**

    The line mode band structure.



### pymatgen.io.phonopy.get_phonon_dos_from_fc(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), supercell_matrix: ndarray, force_constants: ndarray, mesh_density: float = 100.0, num_dos_steps: int = 200, \*\*kwargs)
Get a projected phonon density of states from phonopy force constants.


* **Parameters**


    * **structure** – A structure.


    * **supercell_matrix** – The supercell matrix used to generate the force
    constants.


    * **force_constants** – The force constants in phonopy format.


    * **mesh_density** – The density of the q-point mesh. See the docstring
    for the `mesh` argument in Phonopy.init_mesh() for more details.


    * **num_dos_steps** – Number of frequency steps in the energy grid.


    * **\*\*kwargs** – Additional kwargs passed to the Phonopy constructor.



* **Returns**

    The density of states.



### pymatgen.io.phonopy.get_phonopy_structure(pmg_structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))
Convert a pymatgen Structure object to a PhonopyAtoms object.


* **Parameters**

    **pmg_structure** (*pymatgen Structure*) – A Pymatgen structure object.



### pymatgen.io.phonopy.get_pmg_structure(phonopy_structure: PhonopyAtoms)
Convert a PhonopyAtoms object to pymatgen Structure object.


* **Parameters**

    **phonopy_structure** (*PhonopyAtoms*) – A phonopy structure object.



### pymatgen.io.phonopy.get_structure_from_dict(d)
Extracts a structure from the dictionary extracted from the output
files of phonopy like phonopy.yaml or band.yaml.
Adds “phonopy_masses” in the site_properties of the structures.
Compatible with older phonopy versions.


### pymatgen.io.phonopy.get_thermal_displacement_matrices(thermal_displacements_yaml='thermal_displacement_matrices.yaml', structure_path='POSCAR')
Function to read “thermal_displacement_matrices.yaml” from phonopy and return a list of
ThermalDisplacementMatrices objects
:param thermal_displacements_yaml: path to thermal_displacement_matrices.yaml
:param structure_path: path to POSCAR.

Returns: