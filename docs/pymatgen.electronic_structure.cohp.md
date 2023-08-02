---
layout: default
title: pymatgen.electronic_structure.cohp.md
nav_exclude: true
---

# pymatgen.electronic_structure.cohp module

This module defines classes to represent crystal orbital Hamilton
populations (COHP) and integrated COHP (ICOHP), but can also be used
for crystal orbital overlap populations (COOP) or crystal orbital bond indices (COBIs).
If you use this module, please cite:
J. George, G. Petretto, A. Naik, M. Esters, A. J. Jackson, R. Nelson, R. Dronskowski, G.-M. Rignanese, G. Hautier,
“Automated Bonding Analysis with Crystal Orbital Hamilton Populations”,
ChemPlusChem 2022, e202200123,
DOI: 10.1002/cplu.202200123.


### _class_ pymatgen.electronic_structure.cohp.Cohp(efermi, energies, cohp, are_coops=False, are_cobis=False, icohp=None)
Bases: `MSONable`

Basic COHP object.


* **Parameters**


    * **are_coops** – Indicates whether this object describes COOPs.


    * **are_cobis** – Indicates whether this object describes COBIs.


    * **efermi** – Fermi energy.


    * **energies** – A sequence of energies.


    * **(****{Spin** (*icohp*) – np.array}): representing the COHP for each spin.


    * **(****{Spin** – np.array}): representing the ICOHP for each spin.



#### as_dict()
JSON-serializable dict representation of COHP.


#### _classmethod_ from_dict(dct)
Returns a COHP object from a dict representation of the COHP.


#### get_cohp(spin=None, integrated=False)
Returns the COHP or ICOHP for a particular spin.


* **Parameters**


    * **spin** – Spin. Can be parsed as spin object, integer (-1/1)
    or str (“up”/”down”)


    * **integrated** – Return COHP (False) or ICOHP (True)



* **Returns**

    Returns the CHOP or ICOHP for the input spin. If Spin is
    None and both spins are present, both spins will be returned
    as a dictionary.



#### get_icohp(spin=None)
Convenient alternative to get the ICOHP for a particular spin.


#### get_interpolated_value(energy, integrated=False)
Returns the COHP for a particular energy.


* **Parameters**


    * **energy** – Energy to return the COHP value for.


    * **integrated** – Return COHP (False) or ICOHP (True)



#### has_antibnd_states_below_efermi(spin=None, limit=0.01)
Returns dict indicating if there are antibonding states below the Fermi level depending on the spin

    spin: Spin
    limit: -COHP smaller -limit will be considered.


### _class_ pymatgen.electronic_structure.cohp.CompleteCohp(structure, avg_cohp, cohp_dict, bonds=None, are_coops=False, are_cobis=False, orb_res_cohp=None)
Bases: `Cohp`

A wrapper class that defines an average COHP, and individual COHPs.

<!-- attribute: are_coops

Indicates whether the object is consisting of COOPs. -->
<!-- attribute: are_cobis

Indicates whether the object is consisting of COBIs. -->
<!-- attribute: efermi

Fermi energy -->
<!-- attribute: energies

Sequence of energies -->
<!-- attribute: structure

Structure associated with the COHPs. -->
<!-- attribute: cohp, icohp

The average COHP/ICOHP. -->
<!-- attribute: all_cohps

A dict of COHPs for individual bonds of the form {label: COHP} -->
<!-- attribute: orb_res_cohp

Orbital-resolved COHPs. -->

* **Parameters**


    * **structure** – Structure associated with this COHP.


    * **avg_cohp** – The average cohp as a COHP object.


    * **cohp_dict** – A dict of COHP objects for individual bonds of the form
    {label: COHP}


    * **bonds** – A dict containing information on the bonds of the form
    {label: {key: val}}. The key-val pair can be any information
    the user wants to put in, but typically contains the sites,
    the bond length, and the number of bonds. If nothing is
    supplied, it will default to an empty dict.


    * **are_coops** – indicates whether the Cohp objects are COOPs.
    Defaults to False for COHPs.


    * **are_cobis** – indicates whether the Cohp objects are COBIs.
    Defaults to False for COHPs.


    * **orb_res_cohp** – Orbital-resolved COHPs.



#### as_dict()
JSON-serializable dict representation of CompleteCohp.


#### _classmethod_ from_dict(d)
Returns CompleteCohp object from dict representation.


#### _classmethod_ from_file(fmt, filename=None, structure_file=None, are_coops=False, are_cobis=False)
Creates a CompleteCohp object from an output file of a COHP
calculation. Valid formats are either LMTO (for the Stuttgart
LMTO-ASA code) or LOBSTER (for the LOBSTER code).


* **Parameters**


    * **fmt** – A string for the code that was used to calculate
    the COHPs so that the output file can be handled
    correctly. Can take the values “LMTO” or “LOBSTER”.


    * **filename** – Name of the COHP output file. Defaults to COPL
    for LMTO and COHPCAR.lobster/COOPCAR.lobster for LOBSTER.


    * **structure_file** – Name of the file containing the structure.
    If no file name is given, use CTRL for LMTO and POSCAR
    for LOBSTER.


    * **are_coops** – Indicates whether the populations are COOPs or
    COHPs. Defaults to False for COHPs.


    * **are_cobis** – Indicates whether the populations are COBIs or
    COHPs. Defaults to False for COHPs.



* **Returns**

    A CompleteCohp object.



#### get_cohp_by_label(label, summed_spin_channels=False)
Get specific COHP object.


* **Parameters**


    * **label** – string (for newer Lobster versions: a number)


    * **summed_spin_channels** – bool, will sum the spin channels and return the sum in Spin.up if true



* **Returns**

    Returns the COHP object to simplify plotting



#### get_orbital_resolved_cohp(label, orbitals, summed_spin_channels=False)
Get orbital-resolved COHP.


* **Parameters**


    * **label** – bond label (Lobster: labels as in ICOHPLIST/ICOOPLIST.lobster).


    * **orbitals** – The orbitals as a label, or list or tuple of the form
    [(n1, orbital1), (n2, orbital2)]. Orbitals can either be str,
    int, or Orbital.


    * **summed_spin_channels** – bool, will sum the spin channels and return the sum in Spin.up if true



* **Returns**

    A Cohp object if CompleteCohp contains orbital-resolved cohp,
    or None if it doesn’t.


Note: It currently assumes that orbitals are str if they aren’t the

    other valid types. This is not ideal, but the easiest way to
    avoid unicode issues between python 2 and python 3.


#### get_summed_cohp_by_label_and_orbital_list(label_list, orbital_list, divisor=1, summed_spin_channels=False)
Returns a COHP object that includes a summed COHP divided by divisor.


* **Parameters**


    * **label_list** – list of labels for the COHP that should be included in the summed cohp


    * **orbital_list** – list of orbitals for the COHPs that should be included in the summed cohp (same order as
    label_list)


    * **divisor** – float/int, the summed cohp will be divided by this divisor


    * **summed_spin_channels** – bool, will sum the spin channels and return the sum in Spin.up if true



* **Returns**

    Returns a COHP object including a summed COHP



#### get_summed_cohp_by_label_list(label_list, divisor=1, summed_spin_channels=False)
Returns a COHP object that includes a summed COHP divided by divisor.


* **Parameters**


    * **label_list** – list of labels for the COHP that should be included in the summed cohp


    * **divisor** – float/int, the summed cohp will be divided by this divisor


    * **summed_spin_channels** – bool, will sum the spin channels and return the sum in Spin.up if true



* **Returns**

    Returns a COHP object including a summed COHP



### _class_ pymatgen.electronic_structure.cohp.IcohpCollection(list_labels, list_atom1, list_atom2, list_length, list_translation, list_num, list_icohp, is_spin_polarized, list_orb_icohp=None, are_coops=False, are_cobis=False)
Bases: `MSONable`

Class to store IcohpValues.


#### are_coops()

### Boolean to indicate if these are ICOOPs()

#### are_cobis()

### Boolean to indicate if these are ICOOPs()

#### is_spin_polarized()

### Boolean to indicate if the Lobster calculation was done spin polarized or not()

* **Parameters**


    * **list_labels** – list of labels for ICOHP/ICOOP values


    * **list_atom1** – list of str of atomnames e.g. “O1”


    * **list_atom2** – list of str of atomnames e.g. “O1”


    * **list_length** – list of lengths of corresponding bonds in Angstrom


    * **list_translation** – list of translation list, e.g. [0,0,0]


    * **list_num** – list of equivalent bonds, usually 1 starting from Lobster 3.0.0


    * **list_icohp** – list of dict={Spin.up: icohpvalue for spin.up, Spin.down: icohpvalue for spin.down}


    * **is_spin_polarized** – Boolean to indicate if the Lobster calculation was done spin polarized or not Boolean to
    indicate if the Lobster calculation was done spin polarized or not


    * **list_orb_icohp** – list of dict={[str(Orbital1)-str(Orbital2)]: {“icohp”:{Spin.up: icohpvalue for spin.up,
    Spin.down: icohpvalue for spin.down}, “orbitals”:[Orbital1, Orbital2]}}


    * **are_coops** – Boolean to indicate whether ICOOPs are stored


    * **are_cobis** – Boolean to indicate whether ICOBIs are stored.



#### _property_ are_cobis(_: boo_ )
Whether this a cobi.


* **Type**

    return



#### _property_ are_coops(_: boo_ )
Whether this is a coop.


* **Type**

    return



#### extremum_icohpvalue(summed_spin_channels=True, spin=Spin.up)
get ICOHP/ICOOP of strongest bond
:param summed_spin_channels: Boolean to indicate whether the ICOHPs/ICOOPs of both spin channels should be summed.
:param spin: if summed_spin_channels is equal to False, this spin indicates which spin channel should be returned


* **Returns**

    lowest ICOHP/largest ICOOP value (i.e. ICOHP/ICOOP value of strongest bond)



#### get_icohp_by_label(label, summed_spin_channels=True, spin=Spin.up, orbitals=None)
get an icohp value for a certain bond as indicated by the label (bond labels starting by “1” as in
ICOHPLIST/ICOOPLIST).


* **Parameters**


    * **label** – label in str format (usually the bond number in Icohplist.lobster/Icooplist.lobster


    * **summed_spin_channels** – Boolean to indicate whether the ICOHPs/ICOOPs of both spin channels should be summed


    * **spin** – if summed_spin_channels is equal to False, this spin indicates which spin channel should be returned


    * **orbitals** – List of Orbital or “str(Orbital1)-str(Orbital2)”



* **Returns**

    float describing ICOHP/ICOOP value



#### get_icohp_dict_by_bondlengths(minbondlength=0.0, maxbondlength=8.0)
get a dict of IcohpValues corresponding to certain bond lengths
:param minbondlength: defines the minimum of the bond lengths of the bonds
:param maxbondlength: defines the maximum of the bond lengths of the bonds


* **Returns**

    dict of IcohpValues, the keys correspond to the values from the initial list_labels.



#### get_icohp_dict_of_site(site, minsummedicohp=None, maxsummedicohp=None, minbondlength=0.0, maxbondlength=8.0, only_bonds_to=None)
get a dict of IcohpValue for a certain site (indicated by integer).


* **Parameters**


    * **site** – integer describing the site of interest, order as in Icohplist.lobster/Icooplist.lobster, starts at 0


    * **minsummedicohp** – float, minimal icohp/icoop of the bonds that are considered. It is the summed ICOHP value
    from both spin channels for spin polarized cases


    * **maxsummedicohp** – float, maximal icohp/icoop of the bonds that are considered. It is the summed ICOHP value
    from both spin channels for spin polarized cases


    * **minbondlength** – float, defines the minimum of the bond lengths of the bonds


    * **maxbondlength** – float, defines the maximum of the bond lengths of the bonds


    * **only_bonds_to** – list of strings describing the bonding partners that are allowed, e.g. [‘O’]



* **Returns**

    dict of IcohpValues, the keys correspond to the values from the initial list_labels



#### get_summed_icohp_by_label_list(label_list, divisor=1.0, summed_spin_channels=True, spin=Spin.up)
get the sum of several ICOHP values that are indicated by a list of labels (labels of the bonds are the same as
in ICOHPLIST/ICOOPLIST).


* **Parameters**


    * **label_list** – list of labels of the ICOHPs/ICOOPs that should be summed


    * **divisor** – is used to divide the sum


    * **summed_spin_channels** – Boolean to indicate whether the ICOHPs/ICOOPs of both spin channels should be summed


    * **spin** – if summed_spin_channels is equal to False, this spin indicates which spin channel should be returned



* **Returns**

    float that is a sum of all ICOHPs/ICOOPs as indicated with label_list



#### _property_ is_spin_polarized(_: boo_ )
Whether it is spin polarized.


* **Type**

    return



### _class_ pymatgen.electronic_structure.cohp.IcohpValue(label, atom1, atom2, length, translation, num, icohp, are_coops=False, are_cobis=False, orbitals=None)
Bases: `MSONable`

Class to store information on an ICOHP or ICOOP value.


#### num_bonds()

### number of bonds used for the average cohp (relevant for Lobster versions <3.0) (int)()

#### are_coops()

### Boolean to indicates ICOOPs()

#### are_cobis()

### Boolean to indicates ICOBIs()

#### icohp()

### dict={Spin.up: icohpvalue for spin.up, Spin.down: icohpvalue for spin.down}()

### summed_icohp:()

### sum of icohp/icoop of both spin channels()

* **Parameters**


    * **label** – label for the icohp


    * **atom1** – str of atom that is contributing to the bond


    * **atom2** – str of second atom that is contributing to the bond


    * **length** – float of bond lengths


    * **translation** – translation list, e.g. [0,0,0]


    * **num** – integer describing how often the bond exists


    * **icohp** – dict={Spin.up: icohpvalue for spin.up, Spin.down: icohpvalue for spin.down}


    * **are_coops** – if True, this are COOPs


    * **are_cobis** – if True, this are COBIs


    * **orbitals** – {[str(Orbital1)-str(Orbital2)]: {“icohp”:{Spin.up: icohpvalue for spin.up, Spin.down:
    icohpvalue for spin.down}, “orbitals”:[Orbital1, Orbital2]}}.



#### _property_ are_cobis(_: boo_ )
tells if ICOBIs or not
:returns: Boolean.


#### _property_ are_coops(_: boo_ )
tells if ICOOPs or not
:returns: Boolean.


#### _property_ icohp()
dict with icohps for spinup and spindown
:returns: icohpvalue for spin.up, Spin.down: icohpvalue for spin.down}.
:rtype: dict={Spin.up


#### icohpvalue(spin=Spin.up)

* **Parameters**

    **spin** – Spin.up or Spin.down



* **Returns**

    icohpvalue (float) corresponding to chosen spin.



#### icohpvalue_orbital(orbitals, spin=Spin.up)

* **Parameters**


    * **orbitals** – List of Orbitals or “str(Orbital1)-str(Orbital2)”


    * **spin** – Spin.up or Spin.down



* **Returns**

    icohpvalue (float) corresponding to chosen spin.



#### _property_ is_spin_polarized(_: boo_ )
tells if spin polarized calculation or not
:returns: Boolean.


#### _property_ num_bonds()
tells the number of bonds for which the ICOHP value is an average
:returns: Int.


#### _property_ summed_icohp()
Sums ICOHPs of both spin channels for spin polarized compounds
:returns: icohp value in eV.


#### _property_ summed_orbital_icohp()
Sums orbitals-resolved ICOHPs of both spin channels for spin-plarized compounds
:returns: icohp value in eV}.
:rtype: {“str(Orbital1)-str(Ortibal2)”


### pymatgen.electronic_structure.cohp.get_integrated_cohp_in_energy_range(cohp, label, orbital=None, energy_range=None, relative_E_Fermi=True, summed_spin_channels=True)
Method that can integrate completecohp objects which include data on integrated COHPs
:param cohp: CompleteCOHP object
:param label: label of the COHP data
:param orbital: If not None, a orbital resolved integrated COHP will be returned
:param energy_range: if None, returns icohp value at Fermi level;

> if float, integrates from this float up to the Fermi level;
> if [float,float], will integrate in between


* **Parameters**


    * **relative_E_Fermi** – if True, energy scale with E_Fermi at 0 eV is chosen


    * **summed_spin_channels** – if True, Spin channels will be summed.



* **Returns**

    float indicating the integrated COHP if summed_spin_channels==True, otherwise dict of the following form {
    Spin.up:float, Spin.down:float}