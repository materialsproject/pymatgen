---
layout: default
title: pymatgen.electronic_structure.bandstructure.md
nav_exclude: true
---

# pymatgen.electronic_structure.bandstructure module

This module provides classes to define everything related to band structures.


### _class_ pymatgen.electronic_structure.bandstructure.BandStructure(kpoints: np.ndarray, eigenvals: dict[[Spin](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Spin), np.ndarray], lattice: [Lattice](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice), efermi: float, labels_dict=None, coords_are_cartesian: bool = False, structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure) | None = None, projections: dict[[Spin](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Spin), np.ndarray] | None = None)
Bases: `object`

This is the most generic band structure data possible
it’s defined by a list of kpoints + energies for each of them.


### kpoints:()

### the list of kpoints (as Kpoint objects) in the band structure()

#### lattice_rec()
the reciprocal lattice of the band structure.


#### efermi()
the fermi energy


#### is_spin_polarized()
True if the band structure is spin-polarized, False otherwise


#### bands()
The energy eigenvalues as a {spin: ndarray}. Note that the use of an
ndarray is necessary for computational as well as memory efficiency
due to the large amount of numerical data. The indices of the ndarray
are [band_index, kpoint_index].


#### nb_bands()
returns the number of bands in the band structure


#### structure()
returns the structure


#### projections()
The projections as a {spin: ndarray}. Note that the use of an
ndarray is necessary for computational as well as memory efficiency
due to the large amount of numerical data. The indices of the ndarray
are [band_index, kpoint_index, orbital_index, ion_index].


* **Parameters**


    * **kpoints** – list of kpoint as numpy arrays, in frac_coords of the
    given lattice by default


    * **eigenvals** – dict of energies for spin up and spin down
    {Spin.up:[][],Spin.down:[][]}, the first index of the array
    [][] refers to the band and the second to the index of the
    kpoint. The kpoints are ordered according to the order of the
    kpoints array. If the band structure is not spin polarized, we
    only store one data set under Spin.up


    * **lattice** – The reciprocal lattice as a pymatgen Lattice object.
    Pymatgen uses the physics convention of reciprocal lattice vectors
    WITH a 2\*pi coefficient


    * **efermi** (*float*) – fermi energy


    * **labels_dict** – (dict) of {} this links a kpoint (in frac coords or
    Cartesian coordinates depending on the coords) to a label.


    * **coords_are_cartesian** – Whether coordinates are cartesian.


    * **structure** – The crystal structure (as a pymatgen Structure object)
    associated with the band structure. This is needed if we
    provide projections to the band structure


    * **projections** – dict of orbital projections as {spin: ndarray}. The
    indices of the ndarrayare [band_index, kpoint_index, orbital_index,
    ion_index].If the band structure is not spin polarized, we only
    store one data set under Spin.up.



#### as_dict()
JSON-serializable dict representation of BandStructure.


#### _classmethod_ from_dict(dct)
Create from dict.


* **Parameters**

    **dct** – A dict with all data for a band structure object.



* **Returns**

    A BandStructure object



#### _classmethod_ from_old_dict(dct)

* **Parameters**

    **dct** (*dict*) – A dict with all data for a band structure symmetry line object.



* **Returns**

    A BandStructureSymmLine object



#### get_band_gap()
Returns band gap data.


* **Returns**

    “energy”: band gap energy
    “direct”: A boolean telling if the gap is direct or not
    “transition”: kpoint labels of the transition (e.g., “\\Gamma-X”)



* **Return type**

    A dict {“energy”,”direct”,”transition”}



#### get_cbm()
Returns data about the CBM.


* **Returns**

    {“band_index”,”kpoint_index”,”kpoint”,”energy”}
    - “band_index”: A dict with spin keys pointing to a list of the
    indices of the band containing the CBM (please note that you
    can have several bands sharing the CBM) {Spin.up:[],
    Spin.down:[]}
    - “kpoint_index”: The list of indices in self.kpoints for the
    kpoint CBM. Please note that there can be several
    kpoint_indices relating to the same kpoint (e.g., Gamma can
    occur at different spots in the band structure line plot)
    - “kpoint”: The kpoint (as a kpoint object)
    - “energy”: The energy of the CBM
    - “projections”: The projections along sites and orbitals of the
    CBM if any projection data is available (else it is an empty
    dictionary). The format is similar to the projections field in
    BandStructure: {spin:{‘Orbital’: [proj]}} where the array
    [proj] is ordered according to the sites in structure



#### get_direct_band_gap()
Returns the direct band gap.


* **Returns**

    the value of the direct band gap



#### get_direct_band_gap_dict()
Returns a dictionary of information about the direct
band gap.


* **Returns**

    a dictionary of the band gaps indexed by spin
    along with their band indices and k-point index



#### get_kpoint_degeneracy(kpoint, cartesian=False, tol: float = 0.01)
Returns degeneracy of a given k-point based on structure symmetry
:param kpoint: coordinate of the k-point
:type kpoint: 1x3 array
:param cartesian: kpoint is in Cartesian or fractional coordinates
:type cartesian: bool
:param tol: tolerance below which coordinates are considered equal.
:type tol: float


* **Returns**

    degeneracy or None if structure is not available



* **Return type**

    (int or None)



#### get_projection_on_elements()
Method returning a dictionary of projections on elements.


* **Returns**

    [][{Element:values}],
    Spin.down:[][{Element:values}]} format
    if there is no projections in the band structure
    returns an empty dict



* **Return type**

    a dictionary in the {Spin.up



#### get_projections_on_elements_and_orbitals(el_orb_spec)
Method returning a dictionary of projections on elements and specific
orbitals.


* **Parameters**

    **el_orb_spec** – A dictionary of Elements and Orbitals for which we want
    to have projections on. It is given as: {Element:[orbitals]},
    e.g., {‘Cu’:[‘d’,’s’]}



* **Returns**

    A dictionary of projections on elements in the
    {Spin.up:[][{Element:{orb:values}}],
    Spin.down:[][{Element:{orb:values}}]} format
    if there is no projections in the band structure returns an empty
    dict.



#### get_sym_eq_kpoints(kpoint, cartesian=False, tol: float = 0.01)
Returns a list of unique symmetrically equivalent k-points.


* **Parameters**


    * **kpoint** (*1x3 array*) – coordinate of the k-point


    * **cartesian** (*bool*) – kpoint is in Cartesian or fractional coordinates


    * **tol** (*float*) – tolerance below which coordinates are considered equal



* **Returns**

    if structure is not available returns None



* **Return type**

    ([1x3 array] or None)



#### get_vbm()
Returns data about the VBM.


* **Returns**

    dict as {“band_index”,”kpoint_index”,”kpoint”,”energy”}
    - “band_index”: A dict with spin keys pointing to a list of the
    indices of the band containing the VBM (please note that you
    can have several bands sharing the VBM) {Spin.up:[],
    Spin.down:[]}
    - “kpoint_index”: The list of indices in self.kpoints for the
    kpoint VBM. Please note that there can be several
    kpoint_indices relating to the same kpoint (e.g., Gamma can
    occur at different spots in the band structure line plot)
    - “kpoint”: The kpoint (as a kpoint object)
    - “energy”: The energy of the VBM
    - “projections”: The projections along sites and orbitals of the
    VBM if any projection data is available (else it is an empty
    dictionary). The format is similar to the projections field in
    BandStructure: {spin:{‘Orbital’: [proj]}} where the array
    [proj] is ordered according to the sites in structure



#### is_metal(efermi_tol=0.0001)
Check if the band structure indicates a metal by looking if the fermi
level crosses a band.


* **Returns**

    True if a metal, False if not



### _class_ pymatgen.electronic_structure.bandstructure.BandStructureSymmLine(kpoints, eigenvals, lattice, efermi, labels_dict, coords_are_cartesian=False, structure=None, projections=None)
Bases: `BandStructure`, `MSONable`

This object stores band structures along selected (symmetry) lines in the
Brillouin zone. We call the different symmetry lines (ex: \\Gamma to Z)
“branches”.


* **Parameters**


    * **kpoints** – list of kpoint as numpy arrays, in frac_coords of the
    given lattice by default


    * **eigenvals** – dict of energies for spin up and spin down
    {Spin.up:[][],Spin.down:[][]}, the first index of the array
    [][] refers to the band and the second to the index of the
    kpoint. The kpoints are ordered according to the order of the
    kpoints array. If the band structure is not spin polarized, we
    only store one data set under Spin.up.


    * **lattice** – The reciprocal lattice.
    Pymatgen uses the physics convention of reciprocal lattice vectors
    WITH a 2\*pi coefficient


    * **efermi** – fermi energy


    * **labels_dict** – (dict) of {} this link a kpoint (in frac coords or
    Cartesian coordinates depending on the coords).


    * **coords_are_cartesian** – Whether coordinates are cartesian.


    * **structure** – The crystal structure (as a pymatgen Structure object)
    associated with the band structure. This is needed if we
    provide projections to the band structure.


    * **projections** – dict of orbital projections as {spin: ndarray}. The
    indices of the ndarray are [band_index, kpoint_index, orbital_index,
    ion_index].If the band structure is not spin polarized, we only
    store one data set under Spin.up.



#### apply_scissor(new_band_gap)
Apply a scissor operator (shift of the CBM) to fit the given band gap.
If it’s a metal, we look for the band crossing the Fermi level
and shift this one up. This will not work all the time for metals!


* **Parameters**

    **new_band_gap** – the band gap the scissor band structure need to have.



* **Returns**

    with the applied scissor shift



* **Return type**

    BandStructureSymmLine



#### as_dict()
JSON-serializable dict representation of BandStructureSymmLine.


#### get_branch(index)
Returns in what branch(es) is the kpoint. There can be several
branches.


* **Parameters**

    **index** – the kpoint index



* **Returns**

    A list of dictionaries [{“name”,”start_index”,”end_index”,”index”}]
    indicating all branches in which the k_point is. It takes into
    account the fact that one kpoint (e.g., \\Gamma) can be in several
    branches



#### get_equivalent_kpoints(index)
Returns the list of kpoint indices equivalent (meaning they are the
same frac coords) to the given one.


* **Parameters**

    **index** – the kpoint index



* **Returns**

    a list of equivalent indices


TODO: now it uses the label we might want to use coordinates instead
(in case there was a mislabel)


### _class_ pymatgen.electronic_structure.bandstructure.Kpoint(coords, lattice, to_unit_cell=False, coords_are_cartesian=False, label=None)
Bases: `MSONable`

Class to store kpoint objects. A kpoint is defined with a lattice and frac
or Cartesian coordinates syntax similar than the site object in
pymatgen.core.structure.


* **Parameters**


    * **coords** – coordinate of the kpoint as a numpy array


    * **lattice** – A pymatgen.core.lattice.Lattice lattice object representing
    the reciprocal lattice of the kpoint


    * **to_unit_cell** – Translates fractional coordinate to the basic unit
    cell, i.e., all fractional coordinates satisfy 0 <= a < 1.
    Defaults to False.


    * **coords_are_cartesian** – Boolean indicating if the coordinates given are
    in Cartesian or fractional coordinates (by default fractional)


    * **label** – the label of the kpoint if any (None by default).



#### _property_ a()
Fractional a coordinate of the kpoint.


#### as_dict()
JSON-serializable dict representation of a kpoint.


#### _property_ b()
Fractional b coordinate of the kpoint.


#### _property_ c()
Fractional c coordinate of the kpoint.


#### _property_ cart_coords()
The Cartesian coordinates of the kpoint as a numpy array.


#### _property_ frac_coords()
The fractional coordinates of the kpoint as a numpy array.


#### _classmethod_ from_dict(dct)
Create from dict.


* **Parameters**

    **dct** (*dict*) – A dict with all data for a kpoint object.



* **Returns**

    A Kpoint object



#### _property_ label()
The label associated with the kpoint.


#### _property_ lattice()
The lattice associated with the kpoint. It’s a
pymatgen.core.lattice.Lattice object.


### _class_ pymatgen.electronic_structure.bandstructure.LobsterBandStructureSymmLine(kpoints, eigenvals, lattice, efermi, labels_dict, coords_are_cartesian=False, structure=None, projections=None)
Bases: `BandStructureSymmLine`

Lobster subclass of BandStructure with customized functions.


* **Parameters**


    * **kpoints** – list of kpoint as numpy arrays, in frac_coords of the
    given lattice by default


    * **eigenvals** – dict of energies for spin up and spin down
    {Spin.up:[][],Spin.down:[][]}, the first index of the array
    [][] refers to the band and the second to the index of the
    kpoint. The kpoints are ordered according to the order of the
    kpoints array. If the band structure is not spin polarized, we
    only store one data set under Spin.up.


    * **lattice** – The reciprocal lattice.
    Pymatgen uses the physics convention of reciprocal lattice vectors
    WITH a 2\*pi coefficient


    * **efermi** – fermi energy


    * **labels_dict** – (dict) of {} this link a kpoint (in frac coords or
    Cartesian coordinates depending on the coords).


    * **coords_are_cartesian** – Whether coordinates are cartesian.


    * **structure** – The crystal structure (as a pymatgen Structure object)
    associated with the band structure. This is needed if we
    provide projections to the band structure.


    * **projections** – dict of orbital projections as {spin: ndarray}. The
    indices of the ndarray are [band_index, kpoint_index, orbital_index,
    ion_index].If the band structure is not spin polarized, we only
    store one data set under Spin.up.



#### as_dict()
JSON-serializable dict representation of BandStructureSymmLine.


#### _classmethod_ from_dict(dct)

* **Parameters**

    **dct** (*dict*) – A dict with all data for a band structure symmetry line
    object.



* **Returns**

    A BandStructureSymmLine object



#### _classmethod_ from_old_dict(dct)

* **Parameters**

    **dct** (*dict*) – A dict with all data for a band structure symmetry line
    object.



* **Returns**

    A BandStructureSymmLine object



#### get_projection_on_elements()
Method returning a dictionary of projections on elements.
It sums over all available orbitals for each element.


* **Returns**

    [][{Element:values}],
    Spin.down:[][{Element:values}]} format
    if there is no projections in the band structure
    returns an empty dict



* **Return type**

    a dictionary in the {Spin.up



#### get_projections_on_elements_and_orbitals(el_orb_spec)
Method returning a dictionary of projections on elements and specific
orbitals.


* **Parameters**

    **el_orb_spec** – A dictionary of Elements and Orbitals for which we want
    to have projections on. It is given as: {Element:[orbitals]},
    e.g., {‘Si’:[‘3s’,’3p’]} or {‘Si’:[‘3s’,’3p_x’, ‘3p_y’, ‘3p_z’]} depending on input files



* **Returns**

    A dictionary of projections on elements in the
    {Spin.up:[][{Element:{orb:values}}],
    Spin.down:[][{Element:{orb:values}}]} format
    if there is no projections in the band structure returns an empty
    dict.



### pymatgen.electronic_structure.bandstructure.get_reconstructed_band_structure(list_bs, efermi=None)
This method takes a list of band structures and reconstructs
one band structure object from all of them.

This is typically very useful when you split non self consistent
band structure runs in several independent jobs and want to merge back
the results


* **Parameters**


    * **list_bs** – A list of BandStructure or BandStructureSymmLine objects.


    * **efermi** – The Fermi energy of the reconstructed band structure. If
    None is assigned an average of all the Fermi energy in each
    object in the list_bs is used.



* **Returns**

    A BandStructure or BandStructureSymmLine object (depending on
    the type of the list_bs objects)