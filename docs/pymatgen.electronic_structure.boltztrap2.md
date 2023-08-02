---
layout: default
title: pymatgen.electronic_structure.boltztrap2.md
nav_exclude: true
---

# pymatgen.electronic_structure.boltztrap2 module

BoltzTraP2 is a python software interpolating band structures and
computing materials properties from dft band structure using Boltzmann
semi-classical transport theory.
This module provides a pymatgen interface to BoltzTraP2.
Some of the code is written following the examples provided in BoltzTraP2.

BoltzTraP2 has been developed by Georg Madsen, Jesús Carrete, Matthieu J. Verstraete.

[https://gitlab.com/sousaw/BoltzTraP2](https://gitlab.com/sousaw/BoltzTraP2)
[https://www.sciencedirect.com/science/article/pii/S0010465518301632](https://www.sciencedirect.com/science/article/pii/S0010465518301632)

References are:

> Georg K.H.Madsen, Jesús Carrete, Matthieu J.Verstraete
> BoltzTraP2, a program for interpolating band structures and
> calculating semi-classical transport coefficients
> Computer Physics Communications 231, 140-145, 2018

> Madsen, G. K. H., and Singh, D. J. (2006).
> BoltzTraP. A code for calculating band-structure dependent quantities.
> Computer Physics Communications, 175, 67-71

Todo:
- DONE: spin polarized bands
- read first derivative of the eigenvalues from vasprun.xml (mommat)
- handle magnetic moments (magmom)


### _class_ pymatgen.electronic_structure.boltztrap2.BandstructureLoader(bs_obj, structure=None, nelect=None, mommat=None, magmom=None)
Bases: `object`

Loader for Bandstructure object.


* **Parameters**


    * **bs_obj** – BandStructure object.


    * **structure** – Structure object. It is needed if it is not contained in the BandStructure obj.


    * **nelect** – Number of electrons in the calculation.


    * **mommat** – Matrix of derivatives of energy eigenvalues. TODO Not implemented yet.


    * **magmom** – Matrix of magnetic moments in non collinear calculations. Not implemented yet.


### Example

vrun = Vasprun(‘vasprun.xml’)
bs = vrun.get_band_structure()
st = vrun.final_structure
ne = vrun.parameters[‘NELECT’]
data = BandstructureLoader(bs,st,ne)


#### bandana(emin=-inf, emax=inf)
Cut out bands outside the range (emin,emax).


#### get_lattvec()

* **Returns**

    The lattice vectors.



#### get_volume()

* **Returns**

    Volume



#### set_upper_lower_bands(e_lower, e_upper)
Set fake upper/lower bands, useful to set the same energy
range in the spin up/down bands when calculating the DOS.


### _class_ pymatgen.electronic_structure.boltztrap2.BztInterpolator(data, lpfac=10, energy_range=1.5, curvature=True, save_bztInterp=False, load_bztInterp=False, save_bands=False, fname='bztInterp.json.gz')
Bases: `object`

Interpolate the dft band structures.


* **Parameters**


    * **data** – A loader


    * **lpfac** – the number of interpolation points in the real space. By
    default 10 gives 10 time more points in the real space than
    the number of kpoints given in reciprocal space.


    * **energy_range** – usually the interpolation is not needed on the entire energy
    range but on a specific range around the Fermi level.
    This energy in eV fix the range around the Fermi level
    (E_fermi-energy_range,E_fermi+energy_range) of
    bands that will be interpolated
    and taken into account to calculate the transport properties.


    * **curvature** – boolean value to enable/disable the calculation of second
    derivative related transport properties (Hall coefficient).


    * **save_bztInterp** – Default False. If True coefficients and equivalences are
    saved in fname file.


    * **load_bztInterp** – Default False. If True the coefficients and equivalences
    are loaded from fname file, not calculated. It can be faster than
    re-calculate them in some cases.


    * **save_bands** – Default False. If True interpolated bands are also stored.
    It can be slower than interpolate them. Not recommended.


    * **fname** – File path where to store/load from the coefficients and equivalences.


### Example

data = VasprunLoader().from_file(‘vasprun.xml’)
bztInterp = BztInterpolator(data)


#### get_band_structure(kpaths=None, kpoints_lbls_dict=None, density=20)
Return a BandStructureSymmLine object interpolating bands along a
High symmetry path calculated from the structure using HighSymmKpath
function. If kpaths and kpoints_lbls_dict are provided, a custom
path is interpolated.
kpaths: List of lists of following kpoints labels defining

> the segments of the path. E.g. [[‘L’,’M’],[‘L’,’X’]]

kpoints_lbls_dict: Dict where keys are the kpoint labels used in kpaths

    and values are their fractional coordinates.
    E.g. {‘L’:np.array(0.5,0.5,0.5)},

    > ‘M’:np.array(0.5,0.,0.5),
    > ‘X’:np.array(0.5,0.5,0.)}

density: Number of points in each segment.


#### get_dos(partial_dos=False, npts_mu=10000, T=None, progress=False)
Return a Dos object interpolating bands.


* **Parameters**


    * **partial_dos** – if True, projections will be interpolated as well
    and partial doses will be return. Projections must be available
    in the loader.


    * **npts_mu** – number of energy points of the Dos


    * **T** – parameter used to smooth the Dos


    * **progress** – Default False, If True a progress bar is shown when
    partial dos are computed.



#### get_partial_doses(tdos, eband_ud, spins, enr, npts_mu, T, progress)
Return a CompleteDos object interpolating the projections.

tdos: total dos previously calculated
npts_mu: number of energy points of the Dos
T: parameter used to smooth the Dos
progress: Default False, If True a progress bar is shown.


#### load(fname='bztInterp.json.gz')
Load the coefficient, equivalences, bands from fname.


#### save(fname='bztInterp.json.gz', bands=False)
Save the coefficient, equivalences to fname.
If bands is True, also interpolated bands are stored.


### _class_ pymatgen.electronic_structure.boltztrap2.BztPlotter(bzt_transP=None, bzt_interp=None)
Bases: `object`

Plotter to plot transport properties, interpolated bands along some high
symmetry k-path, and DOS.

### Example

bztPlotter = BztPlotter(bztTransp,bztInterp)
fig = self.bztPlotter.plot_props(‘S’, ‘mu’, ‘temp’, temps=[300, 500])
fig.show()


* **Parameters**


    * **bzt_transP** –


    * **bzt_interp** –



#### plot_bands()
Plot a band structure on symmetry line using BSPlotter().


#### plot_dos(T=None, npoints=10000)
Plot the total Dos using DosPlotter().


#### plot_props(prop_y, prop_x, prop_z='temp', output='avg_eigs', dop_type='n', doping=None, temps=None, xlim=(-2, 2), ax=None)
Function to plot the transport properties.


* **Parameters**


    * **prop_y** – property to plot among (“Conductivity”,”Seebeck”,”Kappa”,”Carrier_conc”,
    “Hall_carrier_conc_trace”). Abbreviations are possible, like “S” for “Seebeck”


    * **prop_x** – independent variable in the x-axis among (‘mu’,’doping’,’temp’)


    * **prop_z** – third variable to plot multiple curves (‘doping’,’temp’)


    * **output** – ‘avg_eigs’ to plot the average of the eigenvalues of the properties
    tensors; ‘eigs’ to plot the three eigenvalues of the properties
    tensors.


    * **dop_type** – ‘n’ or ‘p’ to specify the doping type in plots that use doping
    levels as prop_x or prop_z


    * **doping** – list of doping level to plot, useful to reduce the number of curves
    when prop_z=’doping’


    * **temps** – list of temperatures to plot, useful to reduce the number of curves
    when prop_z=’temp’


    * **xlim** – chemical potential range in eV, useful when prop_x=’mu’


    * **ax** – figure.axes where to plot. If None, a new figure is produced.


Example:
bztPlotter.plot_props(‘S’,’mu’,’temp’,temps=[600,900,1200]).show()
more example are provided in the notebook
“How to use Boltztra2 interface.ipynb”.


### _class_ pymatgen.electronic_structure.boltztrap2.BztTransportProperties(BztInterpolator, temp_r=None, doping=None, npts_mu=4000, CRTA=1e-14, margin=None, save_bztTranspProps=False, load_bztTranspProps=False, fname='bztTranspProps.json.gz')
Bases: `object`

Compute Seebeck, Conductivity, Electrical part of thermal conductivity
and Hall coefficient, conductivity effective mass, Power Factor tensors
w.r.t. the chemical potential and temperatures, from dft band structure via
interpolation.


* **Parameters**


    * **BztInterpolator** – a BztInterpolator previously generated


    * **temp_r** – numpy array of temperatures at which to calculate transport properties


    * **doping** – doping levels at which to calculate transport properties. If provided,
    transport properties w.r.t. these doping levels are also computed. See
    compute_properties_doping() method for details.


    * **npts_mu** – number of energy points at which to calculate transport properties


    * **CRTA** – constant value of the relaxation time


    * **margin** – The energy range of the interpolation is extended by this value on both sides.
    Defaults to 9 \* units.BOLTZMANN \* temp_r.max().


    * **save_bztTranspProps** – Default False. If True all computed transport properties
    will be stored in fname file.


    * **load_bztTranspProps** – Default False. If True all computed transport properties
    will be loaded from fname file.


    * **fname** – File path where to save/load transport properties.


Upon creation, it contains properties tensors w.r.t. the chemical potential
of size (len(temp_r),npts_mu,3,3):

> Conductivity_mu (S/m), Seebeck_mu (microV/K), Kappa_mu (W/(m\*K)),
> Power_Factor_mu (milliW/K m);
> cond_Effective_mass_mu (m_e) calculated as Ref.

Also:

    Carrier_conc_mu: carrier concentration of size (len(temp_r),npts_mu)
    Hall_carrier_conc_trace_mu: trace of Hall carrier concentration of size

    > (len(temp_r),npts_mu)

    mu_r_eV: array of energies in eV and with E_fermi at 0.0

        where all the properties are calculated.

### Example

bztTransp = BztTransportProperties(bztInterp,temp_r = np.arange(100,1400,100))


#### compute_properties_doping(doping, temp_r=None)
Calculate all the properties w.r.t. the doping levels in input.


* **Parameters**


    * **doping** – numpy array specifying the doping levels


    * **temp_r** – numpy array specifying the temperatures


When executed, it add the following variable at the BztTransportProperties
object:

> Conductivity_doping, Seebeck_doping, Kappa_doping, Power_Factor_doping,
> cond_Effective_mass_doping are dictionaries with ‘n’ and ‘p’ keys and
> arrays of dim (len(temp_r),len(doping),3,3) as values.
> Carriers_conc_doping: carriers concentration for each doping level and T.
> mu_doping_eV: the chemical potential corrispondent to each doping level.


#### load(fname='bztTranspProps.json.gz')
Load the transport properties from fname file.


#### save(fname='bztTranspProps.json.gz')
Save the transport properties to fname file.


### _class_ pymatgen.electronic_structure.boltztrap2.VasprunBSLoader(obj, structure=None, nelect=None)
Bases: `object`

Loader for Bandstructure and Vasprun pmg objects.


* **Parameters**


    * **obj** – Either a pmg Vasprun or a BandStructure object.


    * **structure** – Structure object in case is not included in the BandStructure object.


    * **nelect** – number of electrons in case a BandStructure obj is provided.


### Example

vrun = Vasprun(‘vasprun.xml’)
data = VasprunBSLoader(vrun)


#### bandana(emin=-inf, emax=inf)
Cut out bands outside the range (emin,emax).


#### _classmethod_ from_file(vasprun_file)
Get a vasprun.xml file and return a VasprunBSLoader.


#### get_lattvec()

* **Returns**

    The lattice vectors.



#### get_volume()

* **Returns**

    Volume



### _class_ pymatgen.electronic_structure.boltztrap2.VasprunLoader(vrun_obj=None)
Bases: `object`

Loader for Vasprun object.

vrun_obj: Vasprun object.


#### bandana(emin=-inf, emax=inf)
Cut out bands outside the range (emin,emax).


#### _classmethod_ from_file(vasprun_file)
Get a vasprun.xml file and return a VasprunLoader.


#### get_lattvec()

* **Returns**

    Lattice vectors



#### get_volume()

* **Returns**

    Volume of cell



### pymatgen.electronic_structure.boltztrap2.merge_up_down_doses(dos_up, dos_dn)
Merge the up and down DOSs.

Args:
dos_up: Up DOS.
dos_dn: Down DOS
Return:
CompleteDos object