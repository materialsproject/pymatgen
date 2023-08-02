---
layout: default
title: pymatgen.phonon.bandstructure.md
nav_exclude: true
---

# pymatgen.phonon.bandstructure module

This module provides classes to define a phonon band structure.


### _class_ pymatgen.phonon.bandstructure.PhononBandStructure(qpoints: list[[Kpoint](pymatgen.electronic_structure.bandstructure.md#pymatgen.electronic_structure.bandstructure.Kpoint)], frequencies: np.ndarray, lattice: [Lattice](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice), nac_frequencies=None, eigendisplacements=None, nac_eigendisplacements=None, labels_dict=None, coords_are_cartesian=False, structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure) | None = None)
Bases: `MSONable`

This is the most generic phonon band structure data possible
it’s defined by a list of qpoints + frequencies for each of them.
Additional information may be given for frequencies at Gamma, where
non-analytical contribution may be taken into account.


* **Parameters**


    * **qpoints** – list of qpoint as numpy arrays, in frac_coords of the
    given lattice by default


    * **frequencies** – list of phonon frequencies in THz as a numpy array with shape
    (3\*len(structure), len(qpoints)). The First index of the array
    refers to the band and the second to the index of the qpoint.


    * **lattice** – The reciprocal lattice as a pymatgen Lattice object.
    Pymatgen uses the physics convention of reciprocal lattice vectors
    WITH a 2\*pi coefficient.


    * **nac_frequencies** – Frequencies with non-analytical contributions at Gamma in THz.
    A list of tuples. The first element of each tuple should be a list
    defining the direction (not necessarily a versor, will be normalized
    internally). The second element containing the 3\*len(structure)
    phonon frequencies with non-analytical correction for that direction.


    * **eigendisplacements** – the phonon eigendisplacements associated to the
    frequencies in Cartesian coordinates. A numpy array of complex
    numbers with shape (3\*len(structure), len(qpoints), len(structure), 3).
    he First index of the array refers to the band, the second to the index
    of the qpoint, the third to the atom in the structure and the fourth
    to the Cartesian coordinates.


    * **nac_eigendisplacements** – the phonon eigendisplacements associated to the
    non-analytical frequencies in nac_frequencies in Cartesian coordinates.
    A list of tuples. The first element of each tuple should be a list
    defining the direction. The second element containing a numpy array of
    complex numbers with shape (3\*len(structure), len(structure), 3).


    * **labels_dict** – (dict) of {} this links a qpoint (in frac coords or
    Cartesian coordinates depending on the coords) to a label.


    * **coords_are_cartesian** – Whether the qpoint coordinates are Cartesian.


    * **structure** – The crystal structure (as a pymatgen Structure object)
    associated with the band structure. This is needed if we
    provide projections to the band structure.



#### as_dict()

* **Returns**

    MSONable dict



#### asr_breaking(tol_eigendisplacements=1e-05)
Returns the breaking of the acoustic sum rule for the three acoustic modes,
if Gamma is present. None otherwise.
If eigendisplacements are available they are used to determine the acoustic
modes: selects the bands corresponding  to the eigendisplacements that
represent to a translation within tol_eigendisplacements. If these are not
identified or eigendisplacements are missing the first 3 modes will be used
(indices [0:3]).


#### _classmethod_ from_dict(dct)

* **Parameters**

    **dct** – Dict representation



* **Returns**

    PhononBandStructure



#### get_nac_eigendisplacements_along_dir(direction)
Returns the nac_eigendisplacements for the given direction (not necessarily a versor).
None if the direction is not present or nac_eigendisplacements has not been calculated.


* **Parameters**

    **direction** – the direction as a list of 3 elements



* **Returns**

    the eigendisplacements as a numpy array of complex numbers with shape
    (3\*len(structure), len(structure), 3). None if not found.



#### get_nac_frequencies_along_dir(direction)
Returns the nac_frequencies for the given direction (not necessarily a versor).
None if the direction is not present or nac_frequencies has not been calculated.


* **Parameters**

    **direction** – the direction as a list of 3 elements



* **Returns**

    the frequencies as a numpy array o(3\*len(structure), len(qpoints)).
    None if not found.



#### _property_ has_eigendisplacements(_: boo_ )
True if eigendisplacements are present.


#### has_imaginary_freq(tol: float = 1e-05)
True if imaginary frequencies are present in the BS.


#### _property_ has_nac(_: boo_ )
True if nac_frequencies are present.


#### min_freq()
Returns the point where the minimum frequency is reached and its value.


### _class_ pymatgen.phonon.bandstructure.PhononBandStructureSymmLine(qpoints, frequencies, lattice, has_nac=False, eigendisplacements=None, labels_dict=None, coords_are_cartesian=False, structure=None)
Bases: `PhononBandStructure`

This object stores phonon band structures along selected (symmetry) lines in the
Brillouin zone. We call the different symmetry lines (ex: \\Gamma to Z)
“branches”.


* **Parameters**


    * **qpoints** – list of qpoints as numpy arrays, in frac_coords of the
    given lattice by default


    * **frequencies** – list of phonon frequencies in eV as a numpy array with shape
    (3\*len(structure), len(qpoints))


    * **lattice** – The reciprocal lattice as a pymatgen Lattice object.
    Pymatgen uses the physics convention of reciprocal lattice vectors
    WITH a 2\*pi coefficient


    * **has_nac** – specify if the band structure has been produced taking into account
    non-analytical corrections at Gamma. If True frequencies at Gamma from
    different directions will be stored in naf. Default False.


    * **eigendisplacements** – the phonon eigendisplacements associated to the
    frequencies in Cartesian coordinates. A numpy array of complex
    numbers with shape (3\*len(structure), len(qpoints), len(structure), 3).
    he First index of the array refers to the band, the second to the index
    of the qpoint, the third to the atom in the structure and the fourth
    to the Cartesian coordinates.


    * **labels_dict** – (dict) of {} this links a qpoint (in frac coords or
    Cartesian coordinates depending on the coords) to a label.


    * **coords_are_cartesian** – Whether the qpoint coordinates are cartesian.


    * **structure** – The crystal structure (as a pymatgen Structure object)
    associated with the band structure. This is needed if we
    provide projections to the band structure.



#### as_dict()
Returns: MSONable dict.


#### as_phononwebsite()
Return a dictionary with the phononwebsite format:
[http://henriquemiranda.github.io/phononwebsite](http://henriquemiranda.github.io/phononwebsite).


#### band_reorder()
Re-order the eigenvalues according to the similarity of the eigenvectors.


#### _classmethod_ from_dict(dct)

* **Parameters**

    **dct** – Dict representation.


Returns: PhononBandStructureSymmLine


#### get_branch(index)
Returns in what branch(es) is the qpoint. There can be several
branches.


* **Parameters**

    **index** – the qpoint index



* **Returns**

    A list of dictionaries [{“name”,”start_index”,”end_index”,”index”}]
    indicating all branches in which the qpoint is. It takes into
    account the fact that one qpoint (e.g., \\Gamma) can be in several
    branches



#### get_equivalent_qpoints(index)
Returns the list of qpoint indices equivalent (meaning they are the
same frac coords) to the given one.


* **Parameters**

    **index** – the qpoint index



* **Returns**

    a list of equivalent indices


TODO: now it uses the label we might want to use coordinates instead
(in case there was a mislabel)


#### write_phononwebsite(filename)
Write a json file for the phononwebsite:
[http://henriquemiranda.github.io/phononwebsite](http://henriquemiranda.github.io/phononwebsite).


### pymatgen.phonon.bandstructure.eigenvectors_from_displacements(disp, masses)
Calculate the eigenvectors from the atomic displacements.


### pymatgen.phonon.bandstructure.estimate_band_connection(prev_eigvecs, eigvecs, prev_band_order)
A function to order the phonon eigenvectors taken from phonopy.


### pymatgen.phonon.bandstructure.get_reasonable_repetitions(n_atoms: int)
Choose the number of repetitions in a supercell
according to the number of atoms in the system.