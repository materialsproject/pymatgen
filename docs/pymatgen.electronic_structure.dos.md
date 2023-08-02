---
layout: default
title: pymatgen.electronic_structure.dos.md
nav_exclude: true
---

# pymatgen.electronic_structure.dos module

This module defines classes to represent the density of states, etc.


### _class_ pymatgen.electronic_structure.dos.CompleteDos(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), total_dos: Dos, pdoss: Mapping[[PeriodicSite](pymatgen.core.sites.md#pymatgen.core.sites.PeriodicSite), Mapping[[Orbital](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Orbital), Mapping[[Spin](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Spin), ArrayLike]]], normalize: bool = False)
Bases: `Dos`

This wrapper class defines a total dos, and also provides a list of PDos.
Mainly used by pymatgen.io.vasp.Vasprun to create a complete Dos from
a vasprun.xml file. You are unlikely to try to generate this object
manually.


#### structure()
Structure associated with the CompleteDos.


#### pdos()
Dict of partial densities of the form {Site:{Orbital:{Spin:Densities}}}


* **Parameters**


    * **structure** – Structure associated with this particular DOS.


    * **total_dos** – total Dos for structure


    * **pdoss** – The pdoss are supplied as an {Site: {Orbital: {Spin:Densities}}}


    * **normalize** – Whether to normalize the densities by the volume of the structure.
    If True, the units of the densities are states/eV/Angstrom^3. Otherwise,
    the units are states/eV.



#### as_dict()
JSON-serializable dict representation of CompleteDos.


#### _static_ fp_to_dict(fp: NamedTuple)
Converts a fingerprint into a dictionary.


* **Parameters**

    **fp** – The DOS fingerprint to be converted into a dictionary



* **Returns**

    A dict of the fingerprint Keys=type, Values=np.ndarray(energies, densities)



* **Return type**

    dict



#### _classmethod_ from_dict(d)
Returns CompleteDos object from dict representation.


#### get_band_center(band: [OrbitalType](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.OrbitalType) = OrbitalType.d, elements: list[SpeciesLike] | None = None, sites: list[[PeriodicSite](pymatgen.core.sites.md#pymatgen.core.sites.PeriodicSite)] | None = None, spin: [Spin](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Spin) | None = None, erange: list[float] | None = None)
Compute the orbital-projected band center, defined as the first moment
relative to the Fermi level

> int_{-inf}^{+inf} rho(E)\*E dE/int_{-inf}^{+inf} rho(E) dE

based on the work of Hammer and Norskov, Surf. Sci., 343 (1995) where the
limits of the integration can be modified by erange and E is the set
of energies taken with respect to the Fermi level. Note that the band center
is often highly sensitive to the selected erange.


* **Parameters**


    * **band** – Orbital type to get the band center of (default is d-band)


    * **elements** – Elements to get the band center of (cannot be used in conjunction with site)


    * **sites** – Sites to get the band center of (cannot be used in conjunction with el)


    * **spin** – Spin channel to use. By default, the spin channels will be combined.


    * **erange** – [min, max] energy range to consider, with respect to the Fermi level.
    Default is None, which means all energies are considered.



* **Returns**

    band center in eV, often denoted epsilon_d for the d-band center



#### get_band_filling(band: [OrbitalType](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.OrbitalType) = OrbitalType.d, elements: list[SpeciesLike] | None = None, sites: list[[PeriodicSite](pymatgen.core.sites.md#pymatgen.core.sites.PeriodicSite)] | None = None, spin: [Spin](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Spin) | None = None)
Compute the orbital-projected band filling, defined as the zeroth moment
up to the Fermi level.


* **Parameters**


    * **band** – Orbital type to get the band center of (default is d-band)


    * **elements** – Elements to get the band center of (cannot be used in conjunction with site)


    * **sites** – Sites to get the band center of (cannot be used in conjunction with el)


    * **spin** – Spin channel to use. By default, the spin channels will be combined.



* **Returns**

    band filling in eV, often denoted f_d for the d-band



#### get_band_kurtosis(band: [OrbitalType](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.OrbitalType) = OrbitalType.d, elements: list[SpeciesLike] | None = None, sites: list[[PeriodicSite](pymatgen.core.sites.md#pymatgen.core.sites.PeriodicSite)] | None = None, spin: [Spin](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Spin) | None = None, erange: list[float] | None = None)
Get the orbital-projected kurtosis, defined as the fourth standardized moment

    int_{-inf}^{+inf} rho(E)\*(E-E_center)^4 dE/int_{-inf}^{+inf} rho(E) dE)
    /
    (int_{-inf}^{+inf} rho(E)\*(E-E_center)^2 dE/int_{-inf}^{+inf} rho(E) dE))^2

where E_center is the orbital-projected band center, the limits of the integration can be
modified by erange, and E is the set of energies taken with respect to the Fermi level.
Note that the skewness is often highly sensitive to the selected erange.


* **Parameters**


    * **band** – Orbital type to get the band center of (default is d-band)


    * **elements** – Elements to get the band center of (cannot be used in conjunction with site)


    * **sites** – Sites to get the band center of (cannot be used in conjunction with el)


    * **spin** – Spin channel to use. By default, the spin channels will be combined.


    * **erange** – [min, max] energy range to consider, with respect to the Fermi level.
    Default is None, which means all energies are considered.



* **Returns**

    Orbital-projected kurtosis in eV



#### get_band_skewness(band: [OrbitalType](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.OrbitalType) = OrbitalType.d, elements: list[SpeciesLike] | None = None, sites: list[[PeriodicSite](pymatgen.core.sites.md#pymatgen.core.sites.PeriodicSite)] | None = None, spin: [Spin](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Spin) | None = None, erange: list[float] | None = None)
Get the orbital-projected skewness, defined as the third standardized moment

    int_{-inf}^{+inf} rho(E)\*(E-E_center)^3 dE/int_{-inf}^{+inf} rho(E) dE)
    /
    (int_{-inf}^{+inf} rho(E)\*(E-E_center)^2 dE/int_{-inf}^{+inf} rho(E) dE))^(3/2)

where E_center is the orbital-projected band center, the limits of the integration can be
modified by erange, and E is the set of energies taken with respect to the Fermi level.
Note that the skewness is often highly sensitive to the selected erange.


* **Parameters**


    * **band** – Orbitals to get the band center of (default is d-band)


    * **elements** – Elements to get the band center of (cannot be used in conjunction with site)


    * **sites** – Sites to get the band center of (cannot be used in conjunction with el)


    * **spin** – Spin channel to use. By default, the spin channels will be combined.


    * **erange** – [min, max] energy range to consider, with respect to the Fermi level.
    Default is None, which means all energies are considered.



* **Returns**

    Orbital-projected skewness in eV



#### get_band_width(band: [OrbitalType](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.OrbitalType) = OrbitalType.d, elements: list[SpeciesLike] | None = None, sites: list[[PeriodicSite](pymatgen.core.sites.md#pymatgen.core.sites.PeriodicSite)] | None = None, spin: [Spin](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Spin) | None = None, erange: list[float] | None = None)
Get the orbital-projected band width, defined as the square root of the second moment

    sqrt(int_{-inf}^{+inf} rho(E)\*(E-E_center)^2 dE/int_{-inf}^{+inf} rho(E) dE)

where E_center is the orbital-projected band center, the limits of the integration can be
modified by erange, and E is the set of energies taken with respect to the Fermi level.
Note that the band width is often highly sensitive to the selected erange.


* **Parameters**


    * **band** – Orbital type to get the band center of (default is d-band)


    * **elements** – Elements to get the band center of (cannot be used in conjunction with site)


    * **sites** – Sites to get the band center of (cannot be used in conjunction with el)


    * **spin** – Spin channel to use. By default, the spin channels will be combined.


    * **erange** – [min, max] energy range to consider, with respect to the Fermi level.
    Default is None, which means all energies are considered.



* **Returns**

    Orbital-projected band width in eV



#### get_dos_fp(type: str = 'summed_pdos', binning: bool = True, min_e: float | None = None, max_e: float | None = None, n_bins: int = 256, normalize: bool = True)
Generates the DOS fingerprint based on work of
F. Knoop, T. A. r Purcell, M. Scheffler, C. Carbogno, J. Open Source Softw. 2020, 5, 2671.
Source - [https://gitlab.com/vibes-developers/vibes/-/tree/master/vibes/materials_fp](https://gitlab.com/vibes-developers/vibes/-/tree/master/vibes/materials_fp)
Copyright (c) 2020 Florian Knoop, Thomas A.R.Purcell, Matthias Scheffler, Christian Carbogno.


* **Parameters**


    * **type** (*str*) – Specify fingerprint type needed can accept ‘{s/p/d/f/}summed_{pdos/tdos}’


    * **summed_pdos****)** (*(**default is*) –


    * **binning** (*bool*) – If true, the DOS fingerprint is binned using np.linspace and n_bins.
    Default is True.


    * **min_e** (*float*) – The minimum mode energy to include in the fingerprint (default is None)


    * **max_e** (*float*) – The maximum mode energy to include in the fingerprint (default is None)


    * **n_bins** (*int*) – Number of bins to be used in the fingerprint (default is 256)


    * **normalize** (*bool*) – If true, normalizes the area under fp to equal to 1. Default is True.



* **Raises**

    **ValueError** – If type is not one of the accepted values {s/p/d/f/}summed_{pdos/tdos}.



* **Returns**

    The electronic density of states fingerprint
    of format (energies, densities, type, n_bins)



* **Return type**

    Fingerprint(namedtuple)



#### _static_ get_dos_fp_similarity(fp1: NamedTuple, fp2: NamedTuple, col: int = 1, pt: int | str = 'All', normalize: bool = False, tanimoto: bool = False)
Calculates the similarity index (dot product) of two fingerprints.


* **Parameters**


    * **fp1** (*NamedTuple*) – The 1st dos fingerprint object


    * **fp2** (*NamedTuple*) – The 2nd dos fingerprint object


    * **col** (*int*) – The item in the fingerprints (0:energies,1: densities) to take the dot product of (default is 1)


    * **pt** (*int** or **str*) – The index of the point that the dot product is to be taken (default is All)


    * **normalize** (*bool*) – If True normalize the scalar product to 1 (default is False)


    * **tanimoto** (*bool*) – If True will compute Tanimoto index (default is False)



* **Raises**

    **ValueError** – If both tanimoto and normalize are set to True.


Returns:
Similarity index (float): The value of dot product


#### get_element_dos()
Get element projected Dos.


* **Returns**

    Dos}



* **Return type**

    dict of {Element



#### get_element_spd_dos(el: SpeciesLike)
Get element and spd projected Dos.


* **Parameters**

    **el** – Element in Structure.composition associated with CompleteDos



* **Returns**

    Dos}, e.g. {OrbitalType.s: Dos object, …}



* **Return type**

    dict of {OrbitalType



#### get_hilbert_transform(band: [OrbitalType](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.OrbitalType) = OrbitalType.d, elements: list[SpeciesLike] | None = None, sites: list[[PeriodicSite](pymatgen.core.sites.md#pymatgen.core.sites.PeriodicSite)] | None = None)
Return the Hilbert transform of the orbital-projected density of states,
often plotted for a Newns-Anderson analysis.


* **Parameters**


    * **elements** – Elements to get the band center of (cannot be used in conjunction with site)


    * **sites** – Sites to get the band center of (cannot be used in conjunction with el)


    * **band** – Orbitals to get the band center of (default is d-band)



* **Returns**

    Hilbert transformation of the projected DOS.



#### get_n_moment(n: int, band: [OrbitalType](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.OrbitalType) = OrbitalType.d, elements: list[SpeciesLike] | None = None, sites: list[[PeriodicSite](pymatgen.core.sites.md#pymatgen.core.sites.PeriodicSite)] | None = None, spin: [Spin](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Spin) | None = None, erange: list[float] | None = None, center: bool = True)
Get the nth moment of the DOS centered around the orbital-projected band center, defined as

    int_{-inf}^{+inf} rho(E)\*(E-E_center)^n dE/int_{-inf}^{+inf} rho(E) dE

where n is the order, E_center is the orbital-projected band center, the limits of the integration can be
modified by erange, and E is the set of energies taken with respect to the Fermi level. If center is False,
then the E_center reference is not used.


* **Parameters**


    * **n** – The order for the moment


    * **band** – Orbital type to get the band center of (default is d-band)


    * **elements** – Elements to get the band center of (cannot be used in conjunction with site)


    * **sites** – Sites to get the band center of (cannot be used in conjunction with el)


    * **spin** – Spin channel to use. By default, the spin channels will be combined.


    * **erange** – [min, max] energy range to consider, with respect to the Fermi level.
    Default is None, which means all energies are considered.


    * **center** – Take moments with respect to the band center



* **Returns**

    Orbital-projected nth moment in eV



#### get_normalized()
Returns a normalized version of the CompleteDos.


#### get_site_dos(site: [PeriodicSite](pymatgen.core.sites.md#pymatgen.core.sites.PeriodicSite))
Get the total Dos for a site (all orbitals).


* **Parameters**

    **site** – Site in Structure associated with CompleteDos.



* **Returns**

    Dos containing summed orbital densities for site.



#### get_site_orbital_dos(site: [PeriodicSite](pymatgen.core.sites.md#pymatgen.core.sites.PeriodicSite), orbital: [Orbital](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Orbital))
Get the Dos for a particular orbital of a particular site.


* **Parameters**


    * **site** – Site in Structure associated with CompleteDos.


    * **orbital** – Orbital in the site.



* **Returns**

    Dos containing densities for orbital of site.



#### get_site_spd_dos(site: [PeriodicSite](pymatgen.core.sites.md#pymatgen.core.sites.PeriodicSite))
Get orbital projected Dos of a particular site.


* **Parameters**

    **site** – Site in Structure associated with CompleteDos.



* **Returns**

    Dos}, e.g. {OrbitalType.s: Dos object, …}



* **Return type**

    dict of {OrbitalType



#### get_site_t2g_eg_resolved_dos(site: [PeriodicSite](pymatgen.core.sites.md#pymatgen.core.sites.PeriodicSite))
Get the t2g, eg projected DOS for a particular site.


* **Parameters**

    **site** – Site in Structure associated with CompleteDos.



* **Returns**

    A dict {“e_g”: Dos, “t2g”: Dos} containing summed e_g and t2g DOS for the site.



* **Return type**

    dict[str, Dos]



#### get_spd_dos()
Get orbital projected Dos.


* **Returns**

    Dos}, e.g. {OrbitalType.s: Dos object, …}



* **Return type**

    dict of {OrbitalType



#### get_upper_band_edge(band: [OrbitalType](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.OrbitalType) = OrbitalType.d, elements: list[SpeciesLike] | None = None, sites: list[[PeriodicSite](pymatgen.core.sites.md#pymatgen.core.sites.PeriodicSite)] | None = None, spin: [Spin](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Spin) | None = None, erange: list[float] | None = None)
Get the orbital-projected upper band edge. The definition by Xin et al.
Phys. Rev. B, 89, 115114 (2014) is used, which is the highest peak position of the
Hilbert transform of the orbital-projected DOS.


* **Parameters**


    * **band** – Orbital type to get the band center of (default is d-band)


    * **elements** – Elements to get the band center of (cannot be used in conjunction with site)


    * **sites** – Sites to get the band center of (cannot be used in conjunction with el)


    * **spin** – Spin channel to use. By default, the spin channels will be combined.


    * **erange** – [min, max] energy range to consider, with respect to the Fermi level.
    Default is None, which means all energies are considered.



* **Returns**

    Upper band edge in eV, often denoted epsilon_u



#### _property_ spin_polarization(_: float | Non_ )
Calculates spin polarization at Fermi level. If the
calculation is not spin-polarized, None will be
returned.

See Sanvito et al., doi: 10.1126/sciadv.1602241 for
an example usage.


* **Return (float)**

    spin polarization in range [0, 1],


will also return NaN if spin polarization ill-defined
(e.g. for insulator)


### _class_ pymatgen.electronic_structure.dos.DOS(energies: ArrayLike, densities: ArrayLike, efermi: float)
Bases: [`Spectrum`](pymatgen.core.spectrum.md#pymatgen.core.spectrum.Spectrum)

Replacement basic DOS object. All other DOS objects are extended versions
of this object. Work in progress.

<!-- attribute: energies

The sequence of energies -->
<!-- attribute: densities

A dict of spin densities, e.g., {Spin.up: [...], Spin.down: [...]} -->
<!-- attribute: efermi

Fermi level -->

* **Parameters**


    * **energies** – A sequence of energies


    * **densities** (*ndarray*) – Either a Nx1 or a Nx2 array. If former, it is
    interpreted as a Spin.up only density. Otherwise, the first column
    is interpreted as Spin.up and the other is Spin.down.


    * **efermi** – Fermi level energy.



#### XLABEL(_ = 'Energy_ )

#### YLABEL(_ = 'Density_ )

#### get_cbm_vbm(tol: float = 0.001, abs_tol: bool = False, spin=None)
Expects a DOS object and finds the cbm and vbm.


* **Parameters**


    * **tol** – tolerance in occupations for determining the gap


    * **abs_tol** – An absolute tolerance (True) and a relative one (False)


    * **spin** – Possible values are None - finds the gap in the summed
    densities, Up - finds the gap in the up spin channel,
    Down - finds the gap in the down spin channel.



* **Returns**

    float in eV corresponding to the gap



* **Return type**

    (cbm, vbm)



#### get_gap(tol: float = 0.001, abs_tol: bool = False, spin: [Spin](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Spin) | None = None)
Expects a DOS object and finds the gap.


* **Parameters**


    * **tol** – tolerance in occupations for determining the gap


    * **abs_tol** – An absolute tolerance (True) and a relative one (False)


    * **spin** – Possible values are None - finds the gap in the summed
    densities, Up - finds the gap in the up spin channel,
    Down - finds the gap in the down spin channel.



* **Returns**

    gap in eV



#### get_interpolated_gap(tol: float = 0.001, abs_tol: bool = False, spin: [Spin](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Spin) | None = None)
Expects a DOS object and finds the gap.


* **Parameters**


    * **tol** – tolerance in occupations for determining the gap


    * **abs_tol** – Set to True for an absolute tolerance and False for a
    relative one.


    * **spin** – Possible values are None - finds the gap in the summed
    densities, Up - finds the gap in the up spin channel,
    Down - finds the gap in the down spin channel.



* **Returns**

    Tuple of floats in eV corresponding to the gap, cbm and vbm.



* **Return type**

    (gap, cbm, vbm)



### _class_ pymatgen.electronic_structure.dos.Dos(efermi: float, energies: ArrayLike, densities: Mapping[[Spin](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Spin), ArrayLike], norm_vol: float | None = None)
Bases: `MSONable`

Basic DOS object. All other DOS objects are extended versions of this
object.

<!-- attribute: energies

The sequence of energies -->
<!-- attribute: densities

A dict of spin densities, e.g., {Spin.up: [...], Spin.down: [...]} -->
<!-- attribute: efermi

Fermi level -->

* **Parameters**


    * **efermi** – Fermi level energy


    * **energies** – A sequences of energies


    * **(****dict****[****Spin** (*densities*) – np.array]): representing the density of states for each Spin.


    * **norm_vol** – The volume used to normalize the densities. Defaults to 1 if None which will not perform any
    normalization. If not None, the resulting density will have units of states/eV/Angstrom^3, otherwise
    the density will be in states/eV.



#### as_dict()
JSON-serializable dict representation of Dos.


#### _classmethod_ from_dict(d)
Returns Dos object from dict representation of Dos.


#### get_cbm_vbm(tol: float = 0.001, abs_tol: bool = False, spin: [Spin](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Spin) | None = None)
Expects a DOS object and finds the cbm and vbm.


* **Parameters**


    * **tol** – tolerance in occupations for determining the gap


    * **abs_tol** – An absolute tolerance (True) and a relative one (False)


    * **spin** – Possible values are None - finds the gap in the summed
    densities, Up - finds the gap in the up spin channel,
    Down - finds the gap in the down spin channel.



* **Returns**

    float in eV corresponding to the gap



* **Return type**

    (cbm, vbm)



#### get_densities(spin: [Spin](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Spin) | None = None)
Returns the density of states for a particular spin.


* **Parameters**

    **spin** – Spin



* **Returns**

    Returns the density of states for a particular spin. If Spin is
    None, the sum of all spins is returned.



#### get_gap(tol: float = 0.001, abs_tol: bool = False, spin: [Spin](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Spin) | None = None)
Expects a DOS object and finds the gap.


* **Parameters**


    * **tol** – tolerance in occupations for determining the gap


    * **abs_tol** – An absolute tolerance (True) and a relative one (False)


    * **spin** – Possible values are None - finds the gap in the summed
    densities, Up - finds the gap in the up spin channel,
    Down - finds the gap in the down spin channel.



* **Returns**

    gap in eV



#### get_interpolated_gap(tol: float = 0.001, abs_tol: bool = False, spin: [Spin](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Spin) | None = None)
Expects a DOS object and finds the gap.


* **Parameters**


    * **tol** – tolerance in occupations for determining the gap


    * **abs_tol** – Set to True for an absolute tolerance and False for a
    relative one.


    * **spin** – Possible values are None - finds the gap in the summed
    densities, Up - finds the gap in the up spin channel,
    Down - finds the gap in the down spin channel.



* **Returns**

    Tuple of floats in eV corresponding to the gap, cbm and vbm.



* **Return type**

    (gap, cbm, vbm)



#### get_interpolated_value(energy: float)
Returns interpolated density for a particular energy.


* **Parameters**

    **energy** – Energy to return the density for.



#### get_smeared_densities(sigma: float)
Returns the Dict representation of the densities, {Spin: densities},
but with a Gaussian smearing of std dev sigma.


* **Parameters**

    **sigma** – Std dev of Gaussian smearing function.



* **Returns**

    Dict of Gaussian-smeared densities.



### _class_ pymatgen.electronic_structure.dos.FermiDos(dos: Dos, structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure) | None = None, nelecs: float | None = None, bandgap: float | None = None)
Bases: `Dos`, `MSONable`

This wrapper class helps relate the density of states, doping levels
(i.e. carrier concentrations) and corresponding fermi levels. A negative
doping concentration indicates the majority carriers are electrons
(n-type doping); a positive doping concentration indicates holes are the
majority carriers (p-type doping).


* **Parameters**


    * **dos** – Pymatgen Dos object.


    * **structure** – A structure. If not provided, the structure
    of the dos object will be used. If the dos does not have an
    associated structure object, an error will be thrown.


    * **nelecs** – The number of electrons included in the energy range of
    dos. It is used for normalizing the densities. Default is the total
    number of electrons in the structure.


    * **bandgap** – If set, the energy values are scissored so that the electronic
    band gap matches this value.



#### as_dict()
JSON-serializable dict representation of Dos.


#### _classmethod_ from_dict(d)
Returns Dos object from dict representation of Dos.


#### get_doping(fermi_level: float, temperature: float)
Calculate the doping (majority carrier concentration) at a given
Fermi level  and temperature. A simple Left Riemann sum is used for
integrating the density of states over energy & equilibrium Fermi-Dirac
distribution.


* **Parameters**


    * **fermi_level** – The fermi_level level in eV.


    * **temperature** – The temperature in Kelvin.



* **Returns**

    The doping concentration in units of 1/cm^3. Negative values
    indicate that the majority carriers are electrons (n-type doping)
    whereas positive values indicates the majority carriers are holes
    (p-type doping).



#### get_fermi(concentration: float, temperature: float, rtol: float = 0.01, nstep: int = 50, step: float = 0.1, precision: int = 8)
Finds the Fermi level at which the doping concentration at the given
temperature (T) is equal to concentration. A greedy algorithm is used
where the relative error is minimized by calculating the doping at a
grid which continually becomes finer.


* **Parameters**


    * **concentration** – The doping concentration in 1/cm^3. Negative values
    represent n-type doping and positive values represent p-type
    doping.


    * **temperature** – The temperature in Kelvin.


    * **rtol** – The maximum acceptable relative error.


    * **nstep** – The number of steps checked around a given Fermi level.


    * **step** – Initial step in energy when searching for the Fermi level.


    * **precision** – Essentially the decimal places of calculated Fermi level.



* **Raises**

    **ValueError** – If the Fermi level cannot be found.



* **Returns**

    The Fermi level in eV. Note that this is different from the default
    dos.efermi.



#### get_fermi_interextrapolated(concentration: float, temperature: float, warn: bool = True, c_ref: float = 10000000000.0, \*\*kwargs)
Similar to get_fermi except that when get_fermi fails to converge,
an interpolated or extrapolated fermi is returned with the assumption
that the Fermi level changes linearly with log(abs(concentration)).


* **Parameters**


    * **concentration** – The doping concentration in 1/cm^3. Negative values
    represent n-type doping and positive values represent p-type
    doping.


    * **temperature** – The temperature in Kelvin.


    * **warn** – Whether to give a warning the first time the fermi cannot be
    found.


    * **c_ref** – A doping concentration where get_fermi returns a
    value without error for both c_ref and -c_ref.


    * **\*\*kwargs** – Keyword arguments passed to the get_fermi function.



* **Returns**

    The Fermi level. Note, the value is possibly interpolated or
    extrapolated and must be used with caution.



### _class_ pymatgen.electronic_structure.dos.LobsterCompleteDos(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), total_dos: Dos, pdoss: Mapping[[PeriodicSite](pymatgen.core.sites.md#pymatgen.core.sites.PeriodicSite), Mapping[[Orbital](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Orbital), Mapping[[Spin](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Spin), ArrayLike]]], normalize: bool = False)
Bases: `CompleteDos`

Extended CompleteDOS for Lobster.


* **Parameters**


    * **structure** – Structure associated with this particular DOS.


    * **total_dos** – total Dos for structure


    * **pdoss** – The pdoss are supplied as an {Site: {Orbital: {Spin:Densities}}}


    * **normalize** – Whether to normalize the densities by the volume of the structure.
    If True, the units of the densities are states/eV/Angstrom^3. Otherwise,
    the units are states/eV.



#### _classmethod_ from_dict(d)
Hydrate CompleteDos object from dict representation.


#### get_element_spd_dos(el: SpeciesLike)
Get element and spd projected Dos.


* **Parameters**

    **el** – Element in Structure.composition associated with LobsterCompleteDos



* **Returns**

    densities, OrbitalType.p: densities, OrbitalType.d: densities}



* **Return type**

    dict of {OrbitalType.s



#### get_site_orbital_dos(site: [PeriodicSite](pymatgen.core.sites.md#pymatgen.core.sites.PeriodicSite), orbital: str)
Get the Dos for a particular orbital of a particular site.


* **Parameters**


    * **site** – Site in Structure associated with CompleteDos.


    * **orbital** – principal quantum number and orbital in string format, e.g. “4s”.
    possible orbitals are: “s”, “p_y”, “p_z”, “p_x”, “d_xy”, “d_yz”, “d_z^2”,
    “d_xz”, “d_x^2-y^2”, “f_y(3x^2-y^2)”, “f_xyz”,
    “f_yz^2”, “f_z^3”, “f_xz^2”, “f_z(x^2-y^2)”, “f_x(x^2-3y^2)”
    In contrast to the Cohpcar and the Cohplist objects, the strings from the Lobster files are used



* **Returns**

    Dos containing densities of an orbital of a specific site.



#### get_site_t2g_eg_resolved_dos(site: [PeriodicSite](pymatgen.core.sites.md#pymatgen.core.sites.PeriodicSite))
Get the t2g, eg projected DOS for a particular site.


* **Parameters**

    **site** – Site in Structure associated with CompleteDos.



* **Returns**

    Dos, “t2g”: Dos} containing summed e_g and t2g DOS
    for the site.



* **Return type**

    A dict {“e_g”



#### get_spd_dos()
Get orbital projected Dos.
For example, if 3s and 4s are included in the basis of some element, they will be both summed in the orbital
projected DOS.


* **Returns**

    Dos}, e.g. {“s”: Dos object, …}



* **Return type**

    dict of {orbital



### pymatgen.electronic_structure.dos.add_densities(density1: Mapping[[Spin](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Spin), ArrayLike], density2: Mapping[[Spin](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Spin), ArrayLike])
Sum two densities.


* **Parameters**


    * **density1** – First density.


    * **density2** – Second density.



* **Returns**

    dict[Spin, np.ndarray]



### pymatgen.electronic_structure.dos.f0(E, fermi, T)
Return the equilibrium fermi-dirac.


* **Parameters**


    * **E** (*float*) – energy in eV


    * **fermi** (*float*) – the Fermi level in eV


    * **T** (*float*) – the temperature in kelvin



* **Returns**

    float