---
layout: default
title: pymatgen.electronic_structure.plotter.md
nav_exclude: true
---

# pymatgen.electronic_structure.plotter module

This module implements plotter for DOS and band structure.


### _class_ pymatgen.electronic_structure.plotter.BSDOSPlotter(bs_projection: Literal['elements'] | None = 'elements', dos_projection: str = 'elements', vb_energy_range: float = 4, cb_energy_range: float = 4, fixed_cb_energy: bool = False, egrid_interval: float = 1, font: str = 'Times New Roman', axis_fontsize: float = 20, tick_fontsize: float = 15, legend_fontsize: float = 14, bs_legend: str = 'best', dos_legend: str = 'best', rgb_legend: bool = True, fig_size: tuple[float, float] = (11, 8.5))
Bases: `object`

A joint, aligned band structure and density of states plot. Contributions
from Jan Pohls as well as the online example from Germain Salvato-Vallverdu:
[http://gvallver.perso.univ-pau.fr/?p=587](http://gvallver.perso.univ-pau.fr/?p=587).

Instantiate plotter settings.


* **Parameters**


    * **bs_projection** (*'elements'** | **None*) – Whether to project the bands onto elements.


    * **dos_projection** (*str*) – “elements”, “orbitals”, or None


    * **vb_energy_range** (*float*) – energy in eV to show of valence bands


    * **cb_energy_range** (*float*) – energy in eV to show of conduction bands


    * **fixed_cb_energy** (*bool*) – If true, the cb_energy_range will be interpreted
    as constant (i.e., no gap correction for cb energy)


    * **egrid_interval** (*float*) – interval for grid marks


    * **font** (*str*) – font family


    * **axis_fontsize** (*float*) – font size for axis


    * **tick_fontsize** (*float*) – font size for axis tick labels


    * **legend_fontsize** (*float*) – font size for legends


    * **bs_legend** (*str*) – matplotlib string location for legend or None


    * **dos_legend** (*str*) – matplotlib string location for legend or None


    * **rgb_legend** (*bool*) – (T/F) whether to draw RGB triangle/bar for element proj.


    * **fig_size** (*tuple*) – dimensions of figure size (width, height)



#### get_plot(bs: [BandStructureSymmLine](pymatgen.electronic_structure.bandstructure.md#pymatgen.electronic_structure.bandstructure.BandStructureSymmLine), dos: [Dos](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.Dos) | [CompleteDos](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.CompleteDos) | None = None)
Get a matplotlib plot object.


* **Parameters**


    * **bs** ([*BandStructureSymmLine*](pymatgen.electronic_structure.bandstructure.md#pymatgen.electronic_structure.bandstructure.BandStructureSymmLine)) – the bandstructure to plot. Projection
    data must exist for projected plots.


    * **dos** ([*Dos*](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.Dos)) – the Dos to plot. Projection data must exist (i.e.,
    CompleteDos) for projected plots.



* **Returns**

    matplotlib.pyplot object on which you can call commands like show()
    and savefig()



### _class_ pymatgen.electronic_structure.plotter.BSPlotter(bs: [BandStructureSymmLine](pymatgen.electronic_structure.bandstructure.md#pymatgen.electronic_structure.bandstructure.BandStructureSymmLine))
Bases: `object`

Class to plot or get data to facilitate the plot of band structure objects.


* **Parameters**

    **bs** – A BandStructureSymmLine object.



#### add_bs(bs: [BandStructureSymmLine](pymatgen.electronic_structure.bandstructure.md#pymatgen.electronic_structure.bandstructure.BandStructureSymmLine) | list[[BandStructureSymmLine](pymatgen.electronic_structure.bandstructure.md#pymatgen.electronic_structure.bandstructure.BandStructureSymmLine)])
Method to add bands objects to the BSPlotter.


#### bs_plot_data(zero_to_efermi=True, bs=None, bs_ref=None, split_branches=True)
Get the data nicely formatted for a plot.


* **Parameters**


    * **zero_to_efermi** – Automatically set the Fermi level as the plot’s origin (i.e. subtract E - E_f).
    Defaults to True.


    * **bs** – the bandstructure to get the data from. If not provided, the first
    one in the self._bs list will be used.


    * **bs_ref** – is the bandstructure of reference when a rescale of the distances
    is need to plot multiple bands


    * **split_branches** – if True distances and energies are split according to the
    branches. If False distances and energies are split only where branches
    are discontinuous (reducing the number of lines to plot).



* **Returns**

    A dictionary of the following format:
    ticks: A dict with the ‘distances’ at which there is a kpoint (the
    x axis) and the labels (None if no label).
    energy: A dict storing bands for spin up and spin down data
    {Spin:[np.array(nb_bands,kpoints),…]} as a list of discontinuous kpath
    of energies. The energy of multiple continuous branches are stored together.
    vbm: A list of tuples (distance,energy) marking the vbms. The
    energies are shifted with respect to the Fermi level is the
    option has been selected.
    cbm: A list of tuples (distance,energy) marking the cbms. The
    energies are shifted with respect to the Fermi level is the
    option has been selected.
    lattice: The reciprocal lattice.
    zero_energy: This is the energy used as zero for the plot.
    band_gap:A string indicating the band gap and its nature (empty if
    it’s a metal).
    is_metal: True if the band structure is metallic (i.e., there is at
    least one band crossing the Fermi level).



* **Return type**

    dict



#### get_plot(zero_to_efermi=True, ylim=None, smooth=False, vbm_cbm_marker=False, smooth_tol=0, smooth_k=3, smooth_np=100, bs_labels=None)
Get a matplotlib object for the bandstructures plot.
Multiple bandstructure objs are plotted together if they have the
same high symm path.


* **Parameters**


    * **zero_to_efermi** – Automatically set the Fermi level as the plot’s origin (i.e. subtract E - E_f).
    Defaults to True.


    * **ylim** – Specify the y-axis (energy) limits; by default None let
    the code choose. It is vbm-4 and cbm+4 if insulator
    efermi-10 and efermi+10 if metal


    * **smooth** (*bool** or **list**(**bools**)*) – interpolates the bands by a spline cubic.
    A single bool values means to interpolate all the bandstructure objs.
    A list of bools allows to select the bandstructure obs to interpolate.


    * **vbm_cbm_marker** (*bool*) – if True, a marker is added to the vbm and cbm.


    * **smooth_tol** (*float*) – tolerance for fitting spline to band data.
    Default is None such that no tolerance will be used.


    * **smooth_k** (*int*) – degree of splines 1<k<5


    * **smooth_np** (*int*) – number of interpolated points per each branch.


    * **bs_labels** – labels for each band for the plot legend.



#### get_ticks()
Get all ticks and labels for a band structure plot.


* **Returns**

    A dictionary with ‘distance’: a list of distance at which
    ticks should be set and ‘label’: a list of label for each of those
    ticks.



* **Return type**

    dict



#### get_ticks_old()
Get all ticks and labels for a band structure plot.


* **Returns**

    A dictionary with ‘distance’: a list of distance at which
    ticks should be set and ‘label’: a list of label for each of those
    ticks.



* **Return type**

    dict



#### plot_brillouin()
Plot the Brillouin zone.


#### plot_compare(other_plotter, legend=True)
Plot two band structure for comparison. One is in red the other in blue
(no difference in spins). The two band structures need to be defined
on the same symmetry lines! and the distance between symmetry lines is
the one of the band structure used to build the BSPlotter.


* **Parameters**


    * **other_plotter** – Another band structure object defined along the same symmetry lines


    * **legend** – True to add a legend to the plot



* **Returns**

    a matplotlib object with both band structures



#### save_plot(filename, img_format='eps', ylim=None, zero_to_efermi=True, smooth=False)
Save matplotlib plot to a file.


* **Parameters**


    * **filename** – Filename to write to.


    * **img_format** – Image format to use. Defaults to EPS.


    * **ylim** – Specifies the y-axis limits.


    * **zero_to_efermi** – Automatically set the Fermi level as the plot’s origin (i.e. subtract E - E_f).
    Defaults to True.


    * **smooth** – Cubic spline interpolation of the bands.



#### show(zero_to_efermi=True, ylim=None, smooth=False, smooth_tol=None)
Show the plot using matplotlib.


* **Parameters**


    * **zero_to_efermi** – Automatically set the Fermi level as the plot’s origin (i.e. subtract E - E_f).
    Defaults to True.


    * **ylim** – Specify the y-axis (energy) limits; by default None let
    the code choose. It is vbm-4 and cbm+4 if insulator
    efermi-10 and efermi+10 if metal


    * **smooth** – interpolates the bands by a spline cubic


    * **smooth_tol** (*float*) – tolerance for fitting spline to band data.
    Default is None such that no tolerance will be used.



### _class_ pymatgen.electronic_structure.plotter.BSPlotterProjected(bs)
Bases: `BSPlotter`

Class to plot or get data to facilitate the plot of band structure objects
projected along orbitals, elements or sites.


* **Parameters**

    **bs** – A BandStructureSymmLine object with projections.



#### get_elt_projected_plots(zero_to_efermi=True, ylim=None, vbm_cbm_marker=False)
Method returning a plot composed of subplots along different elements.


* **Returns**

    a pyplot object with different subfigures for each projection
    The blue and red colors are for spin up and spin down
    The bigger the red or blue dot in the band structure the higher
    character for the corresponding element and orbital



#### get_elt_projected_plots_color(zero_to_efermi=True, elt_ordered=None)
Returns a pyplot plot object with one plot where the band structure
line color depends on the character of the band (along different
elements). Each element is associated with red, green or blue
and the corresponding rgb color depending on the character of the band
is used. The method can only deal with binary and ternary compounds.

spin up and spin down are differientiated by a ‘-’ and a ‘–’ line


* **Parameters**


    * **zero_to_efermi** – Automatically set the Fermi level as the plot’s origin (i.e. subtract E - E_f).
    Defaults to True.


    * **elt_ordered** – A list of Element ordered. The first one is red,
    second green, last blue



* **Returns**

    a pyplot object



#### get_projected_plots_dots(dictio, zero_to_efermi=True, ylim=None, vbm_cbm_marker=False)
Method returning a plot composed of subplots along different elements
and orbitals.


* **Parameters**


    * **dictio** – The element and orbitals you want a projection on. The
    format is {Element:[Orbitals]} for instance
    {‘Cu’:[‘d’,’s’],’O’:[‘p’]} will give projections for Cu on
    d and s orbitals and on oxygen p.
    If you use this class to plot LobsterBandStructureSymmLine,
    the orbitals are named as in the FATBAND filename, e.g.
    “2p” or “2p_x”


    * **zero_to_efermi** – Automatically set the Fermi level as the plot’s origin (i.e. subtract E - E_f).
    Defaults to True.


    * **ylim** – Specify the y-axis limits. Defaults to None.


    * **vbm_cbm_marker** – Add markers for the VBM and CBM. Defaults to False.



* **Returns**

    a pyplot object with different subfigures for each projection
    The blue and red colors are for spin up and spin down.
    The bigger the red or blue dot in the band structure the higher
    character for the corresponding element and orbital.



#### get_projected_plots_dots_patom_pmorb(dictio, dictpa, sum_atoms=None, sum_morbs=None, zero_to_efermi=True, ylim=None, vbm_cbm_marker=False, selected_branches=None, w_h_size=(12, 8), num_column=None)
Method returns a plot composed of subplots for different atoms and
orbitals (subshell orbitals such as ‘s’, ‘p’, ‘d’ and ‘f’ defined by
azimuthal quantum numbers l = 0, 1, 2 and 3, respectively or
individual orbitals like ‘px’, ‘py’ and ‘pz’ defined by magnetic
quantum numbers m = -1, 1 and 0, respectively).
This is an extension of “get_projected_plots_dots” method.


* **Parameters**


    * **dictio** – The elements and the orbitals you need to project on. The
    format is {Element:[Orbitals]}, for instance:
    {‘Cu’:[‘dxy’,’s’,’px’],’O’:[‘px’,’py’,’pz’]} will give projections for Cu on
    orbitals dxy, s, px and for O on orbitals px, py, pz. If you want to sum over all
    individual orbitals of subshell orbitals, for example, ‘px’, ‘py’ and ‘pz’ of O,
    just simply set {‘Cu’:[‘dxy’,’s’,’px’],’O’:[‘p’]} and set sum_morbs (see
    explanations below) as {‘O’:[p],…}. Otherwise, you will get an error.


    * **dictpa** – The elements and their sites (defined by site numbers) you
    need to project on. The format is {Element: [Site numbers]}, for instance:
    {‘Cu’:[1,5],’O’:[3,4]} will give projections for Cu on site-1 and on site-5, O on
    site-3 and on site-4 in the cell. The correct site numbers of atoms are consistent
    with themselves in the structure computed. Normally, the structure should be totally
    similar with POSCAR file, however, sometimes VASP can rotate or translate the cell.
    Thus, it would be safe if using Vasprun class to get the final_structure and as a
    result, correct index numbers of atoms.


    * **sum_atoms** – Sum projection of the similar atoms together (e.g.: Cu
    on site-1 and Cu on site-5). The format is {Element: [Site numbers]}, for instance:

    > {‘Cu’: [1,5], ‘O’: [3,4]} means summing projections over Cu on site-1 and Cu on
    > site-5 and O on site-3 and on site-4. If you do not want to use this functional,
    > just turn it off by setting sum_atoms = None.



    * **sum_morbs** – Sum projections of individual orbitals of similar atoms
    together (e.g.: ‘dxy’ and ‘dxz’). The format is {Element: [individual orbitals]},
    for instance: {‘Cu’: [‘dxy’, ‘dxz’], ‘O’: [‘px’, ‘py’]} means summing projections
    over ‘dxy’ and ‘dxz’ of Cu and ‘px’ and ‘py’ of O. If you do not want to use this
    functional, just turn it off by setting sum_morbs = None.


    * **zero_to_efermi** – Automatically set the Fermi level as the plot’s origin (i.e. subtract E - E_f).
    Defaults to True.


    * **ylim** – The y-axis limit. Defaults to None.


    * **vbm_cbm_marker** – Whether to plot points to indicate valence band maxima and conduction
    band minima positions. Defaults to False.


    * **selected_branches** – The index of symmetry lines you chose for
    plotting. This can be useful when the number of symmetry lines (in KPOINTS file) are
    manny while you only want to show for certain ones. The format is [index of line],
    for instance: [1, 3, 4] means you just need to do projection along lines number 1, 3
    and 4 while neglecting lines number 2 and so on. By default, this is None type and
    all symmetry lines will be plotted.


    * **w_h_size** – This variable help you to control the width and height
    of figure. By default, width = 12 and height = 8 (inches). The width/height ratio is
    kept the same for subfigures and the size of each depends on how many number of
    subfigures are plotted.


    * **num_column** – This variable help you to manage how the subfigures are
    arranged in the figure by setting up the number of columns of subfigures. The value
    should be an int number. For example, num_column = 3 means you want to plot
    subfigures in 3 columns. By default, num_column = None and subfigures are aligned in
    2 columns.



* **Returns**

    A pyplot object with different subfigures for different projections.
    The blue and red colors lines are bands
    for spin up and spin down. The green and cyan dots are projections
    for spin up and spin down. The bigger
    the green or cyan dots in the projected band structures, the higher
    character for the corresponding elements
    and orbitals. List of individual orbitals and their numbers (set up
    by VASP and no special meaning):
    s = 0; py = 1 pz = 2 px = 3; dxy = 4 dyz = 5 dz2 = 6 dxz = 7 dx2 = 8;
    f_3 = 9 f_2 = 10 f_1 = 11 f0 = 12 f1 = 13 f2 = 14 f3 = 15



### _class_ pymatgen.electronic_structure.plotter.BoltztrapPlotter(bz)
Bases: `object`

class containing methods to plot the data from Boltztrap.


* **Parameters**

    **bz** – a BoltztrapAnalyzer object.



#### plot_carriers(temp=300)
Plot the carrier concentration in function of Fermi level.


* **Parameters**

    **temp** – the temperature



* **Returns**

    a matplotlib object



#### plot_complexity_factor_mu(temps=(300,), output='average', Lambda=0.5)
Plot respect to the chemical potential of the Fermi surface complexity
factor calculated as explained in Ref.
Gibbs, Z. M. et al., Effective mass and fermi surface complexity factor
from ab initio band structure calculations.
npj Computational Materials 3, 8 (2017).


* **Parameters**


    * **output** – ‘average’ returns the complexity factor calculated using the average
    of the three diagonal components of the seebeck and conductivity tensors.
    ‘tensor’ returns the complexity factor respect to the three
    diagonal components of seebeck and conductivity tensors.


    * **temps** – list of temperatures of calculated seebeck and conductivity.


    * **Lambda** – fitting parameter used to model the scattering (0.5 means constant
    relaxation time).



* **Returns**

    a matplotlib object



#### plot_conductivity_dop(temps='all', output='average', relaxation_time=1e-14)
Plot the conductivity in function of doping levels for different
temperatures.


* **Parameters**


    * **temps** – the default ‘all’ plots all the temperatures in the analyzer.
    Specify a list of temperatures if you want to plot only some.


    * **output** – with ‘average’ you get an average of the three directions
    with ‘eigs’ you get all the three directions.


    * **relaxation_time** – specify a constant relaxation time value



* **Returns**

    a matplotlib object



#### plot_conductivity_mu(temp: float = 600, output: str = 'eig', relaxation_time: float = 1e-14, xlim: Sequence[float] | None = None)
Plot the conductivity in function of Fermi level. Semi-log plot.


* **Parameters**


    * **temp** (*float*) – the temperature


    * **output** (*str*) – “eig” or “average”


    * **relaxation_time** (*float*) – A relaxation time in s. Defaults to 1e-14 and the plot is in
    units of relaxation time


    * **xlim** (*tuple**[**float**, **float**]*) – a 2-tuple of min and max fermi energy. Defaults to (0, band gap)



* **Returns**

    a matplotlib object



#### plot_conductivity_temp(doping='all', output='average', relaxation_time=1e-14)
Plot the conductivity in function of temperature for different doping levels.


* **Parameters**


    * **doping** (*str*) – the default ‘all’ plots all the doping levels in the analyzer.
    Specify a list of doping levels if you want to plot only some.


    * **output** – with ‘average’ you get an average of the three directions
    with ‘eigs’ you get all the three directions.


    * **relaxation_time** – specify a constant relaxation time value



* **Returns**

    a matplotlib object



#### plot_dos(sigma=0.05)
Plot dos.


* **Parameters**

    **sigma** – a smearing



* **Returns**

    a matplotlib object



#### plot_eff_mass_dop(temps='all', output='average')
Plot the average effective mass in function of doping levels
for different temperatures.


* **Parameters**


    * **temps** – the default ‘all’ plots all the temperatures in the analyzer.
    Specify a list of temperatures if you want to plot only some.


    * **output** – with ‘average’ you get an average of the three directions
    with ‘eigs’ you get all the three directions.


    * **relaxation_time** – specify a constant relaxation time value



* **Returns**

    a matplotlib object



#### plot_eff_mass_temp(doping='all', output='average')
Plot the average effective mass in function of temperature
for different doping levels.


* **Parameters**


    * **doping** (*str*) – the default ‘all’ plots all the doping levels in the analyzer.
    Specify a list of doping levels if you want to plot only some.


    * **output** – with ‘average’ you get an average of the three directions
    with ‘eigs’ you get all the three directions.



* **Returns**

    a matplotlib object



#### plot_hall_carriers(temp=300)
Plot the Hall carrier concentration in function of Fermi level.


* **Parameters**

    **temp** – the temperature



* **Returns**

    a matplotlib object



#### plot_power_factor_dop(temps='all', output='average', relaxation_time=1e-14)
Plot the Power Factor in function of doping levels for different temperatures.


* **Parameters**


    * **temps** – the default ‘all’ plots all the temperatures in the analyzer.
    Specify a list of temperatures if you want to plot only some.


    * **output** – with ‘average’ you get an average of the three directions
    with ‘eigs’ you get all the three directions.


    * **relaxation_time** – specify a constant relaxation time value



* **Returns**

    a matplotlib object



#### plot_power_factor_mu(temp: float = 600, output: str = 'eig', relaxation_time: float = 1e-14, xlim: Sequence[float] | None = None)
Plot the power factor in function of Fermi level. Semi-log plot.


* **Parameters**


    * **temp** (*float*) – the temperature


    * **output** (*str*) – “eig” or “average”


    * **relaxation_time** (*float*) – A relaxation time in s. Defaults to 1e-14 and the plot is in
    units of relaxation time


    * **xlim** (*tuple**[**float**, **float**]*) – a 2-tuple of min and max fermi energy. Defaults to (0, band gap)



* **Returns**

    a matplotlib object



#### plot_power_factor_temp(doping='all', output='average', relaxation_time=1e-14)
Plot the Power Factor in function of temperature for different doping levels.


* **Parameters**


    * **doping** (*str*) – the default ‘all’ plots all the doping levels in the analyzer.
    Specify a list of doping levels if you want to plot only some.


    * **output** – with ‘average’ you get an average of the three directions
    with ‘eigs’ you get all the three directions.


    * **relaxation_time** – specify a constant relaxation time value



* **Returns**

    a matplotlib object



#### plot_seebeck_dop(temps='all', output='average')
Plot the Seebeck in function of doping levels for different temperatures.


* **Parameters**


    * **temps** – the default ‘all’ plots all the temperatures in the analyzer.
    Specify a list of temperatures if you want to plot only some.


    * **output** – with ‘average’ you get an average of the three directions
    with ‘eigs’ you get all the three directions.



* **Returns**

    a matplotlib object



#### plot_seebeck_eff_mass_mu(temps=(300,), output='average', Lambda=0.5)
Plot respect to the chemical potential of the Seebeck effective mass
calculated as explained in Ref.
Gibbs, Z. M. et al., Effective mass and fermi surface complexity factor
from ab initio band structure calculations.
npj Computational Materials 3, 8 (2017).


* **Parameters**


    * **output** – ‘average’ returns the seebeck effective mass calculated
    using the average of the three diagonal components of the
    seebeck tensor. ‘tensor’ returns the seebeck effective mass
    respect to the three diagonal components of the seebeck tensor.


    * **temps** – list of temperatures of calculated seebeck.


    * **Lambda** – fitting parameter used to model the scattering (0.5 means
    constant relaxation time).



* **Returns**

    a matplotlib object



#### plot_seebeck_mu(temp: float = 600, output: str = 'eig', xlim: Sequence[float] | None = None)
Plot the seebeck coefficient in function of Fermi level.


* **Parameters**


    * **temp** (*float*) – the temperature


    * **output** (*str*) – “eig” or “average”


    * **xlim** (*tuple**[**float**, **float**]*) – a 2-tuple of min and max fermi energy. Defaults to (0, band gap)



* **Returns**

    a matplotlib object



#### plot_seebeck_temp(doping='all', output='average')
Plot the Seebeck coefficient in function of temperature for different
doping levels.


* **Parameters**


    * **doping** (*str*) – the default ‘all’ plots all the doping levels in the analyzer.
    Specify a list of doping levels if you want to plot only some.


    * **output** – with ‘average’ you get an average of the three directions
    with ‘eigs’ you get all the three directions.



* **Returns**

    a matplotlib object



#### plot_zt_dop(temps='all', output='average', relaxation_time=1e-14)
Plot the figure of merit zT in function of doping levels for different
temperatures.


* **Parameters**


    * **temps** – the default ‘all’ plots all the temperatures in the analyzer.
    Specify a list of temperatures if you want to plot only some.


    * **output** – with ‘average’ you get an average of the three directions
    with ‘eigs’ you get all the three directions.


    * **relaxation_time** – specify a constant relaxation time value



* **Returns**

    a matplotlib object



#### plot_zt_mu(temp: float = 600, output: str = 'eig', relaxation_time: float = 1e-14, xlim: Sequence[float] | None = None)
Plot the ZT in function of Fermi level.


* **Parameters**


    * **temp** (*float*) – the temperature


    * **output** (*str*) – “eig” or “average”


    * **relaxation_time** (*float*) – A relaxation time in s. Defaults to 1e-14 and the plot is in
    units of relaxation time


    * **xlim** (*tuple**[**float**, **float**]*) – a 2-tuple of min and max fermi energy. Defaults to (0, band gap)



* **Returns**

    matplotlib.pyplot module



#### plot_zt_temp(doping='all', output: Literal['average', 'eigs'] = 'average', relaxation_time=1e-14)
Plot the figure of merit zT in function of temperature for different doping levels.


* **Parameters**


    * **doping** (*str*) – the default ‘all’ plots all the doping levels in the analyzer.
    Specify a list of doping levels if you want to plot only some.


    * **output** – with ‘average’ you get an average of the three directions
    with ‘eigs’ you get all the three directions.


    * **relaxation_time** – specify a constant relaxation time value



* **Raises**

    **ValueError** – if output is not ‘average’ or ‘eigs’



* **Returns**

    a matplotlib object



### _class_ pymatgen.electronic_structure.plotter.CohpPlotter(zero_at_efermi=True, are_coops=False, are_cobis=False)
Bases: `object`

Class for plotting crystal orbital Hamilton populations (COHPs) or
crystal orbital overlap populations (COOPs). It is modeled after the
DosPlotter object.


* **Parameters**


    * **zero_at_efermi** – Whether to shift all populations to have zero
    energy at the Fermi level. Defaults to True.


    * **are_coops** – Switch to indicate that these are COOPs, not COHPs.
    Defaults to False for COHPs.


    * **are_cobis** – Switch to indicate that these are COBIs, not COHPs/COOPs.
    Defaults to False for COHPs.



#### add_cohp(label, cohp)
Adds a COHP for plotting.


* **Parameters**


    * **label** – Label for the COHP. Must be unique.


    * **cohp** – COHP object.



#### add_cohp_dict(cohp_dict, key_sort_func=None)
Adds a dictionary of COHPs with an optional sorting function
for the keys.


* **Parameters**


    * **cohp_dict** – dict of the form {label: Cohp}


    * **key_sort_func** – function used to sort the cohp_dict keys.



#### get_cohp_dict()
Returns the added COHPs as a json-serializable dict. Note that if you
have specified smearing for the COHP plot, the populations returned
will be the smeared and not the original populations.


* **Returns**

    Dict of COHP data of the form {label: {“efermi”: efermi,
    “energies”: …, “COHP”: {Spin.up: …}, “ICOHP”: …}}.



* **Return type**

    dict



#### get_plot(xlim=None, ylim=None, plot_negative=None, integrated=False, invert_axes=True)
Get a matplotlib plot showing the COHP.


* **Parameters**


    * **xlim** – Specifies the x-axis limits. Defaults to None for
    automatic determination.


    * **ylim** – Specifies the y-axis limits. Defaults to None for
    automatic determination.


    * **plot_negative** – It is common to plot -COHP(E) so that the
    sign means the same for COOPs and COHPs. Defaults to None
    for automatic determination: If are_coops is True, this
    will be set to False, else it will be set to True.


    * **integrated** – Switch to plot ICOHPs. Defaults to False.


    * **invert_axes** – Put the energies onto the y-axis, which is
    common in chemistry.



* **Returns**

    A matplotlib object.



#### save_plot(filename, img_format='eps', xlim=None, ylim=None)
Save matplotlib plot to a file.


* **Parameters**


    * **filename** – File name to write to.


    * **img_format** – Image format to use. Defaults to EPS.


    * **xlim** – Specifies the x-axis limits. Defaults to None for
    automatic determination.


    * **ylim** – Specifies the y-axis limits. Defaults to None for
    automatic determination.



#### show(xlim=None, ylim=None)
Show the plot using matplotlib.


* **Parameters**


    * **xlim** – Specifies the x-axis limits. Defaults to None for
    automatic determination.


    * **ylim** – Specifies the y-axis limits. Defaults to None for
    automatic determination.



### _class_ pymatgen.electronic_structure.plotter.DosPlotter(zero_at_efermi: bool = True, stack: bool = False, sigma: float | None = None)
Bases: `object`

Class for plotting phonon DOSs. The interface is extremely flexible given there are many
different ways in which people want to view DOS.
Typical usage is:

> # Initializes plotter with some optional args. Defaults are usually fine
> plotter = PhononDosPlotter().

> # Add DOS with a label
> plotter.add_dos(“Total DOS”, dos)

> # Alternatively, you can add a dict of DOSes. This is the typical form
> # returned by CompletePhononDos.get_element_dos().
> plotter.add_dos_dict({“dos1”: dos1, “dos2”: dos2})
> plotter.add_dos_dict(complete_dos.get_spd_dos())


* **Parameters**


    * **zero_at_efermi** (*bool*) – Whether to shift all Dos to have zero energy at the
    fermi energy. Defaults to True.


    * **stack** (*bool*) – Whether to plot the DOS as a stacked area graph


    * **sigma** (*float*) – Specify a standard deviation for Gaussian smearing
    the DOS for nicer looking plots. Defaults to None for no
    smearing.



#### add_dos(label: str, dos: [Dos](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.Dos))
Adds a dos for plotting.


* **Parameters**


    * **label** – label for the DOS. Must be unique.


    * **dos** – Dos object



#### add_dos_dict(dos_dict, key_sort_func=None)
Add a dictionary of doses, with an optional sorting function for the
keys.


* **Parameters**


    * **dos_dict** – dict of {label: Dos}


    * **key_sort_func** – function used to sort the dos_dict keys.



#### get_dos_dict()
Returns the added doses as a json-serializable dict. Note that if you
have specified smearing for the DOS plot, the densities returned will
be the smeared densities, not the original densities.


* **Returns**

    Dict of dos data. Generally of the form
    {label: {‘energies’:…, ‘densities’: {‘up’:…}, ‘efermi’:efermi}}



* **Return type**

    dict



#### get_plot(xlim: tuple[float, float] | None = None, ylim: tuple[float, float] | None = None, invert_axes: bool = False, beta_dashed: bool = False)
Get a matplotlib plot showing the DOS.


* **Parameters**


    * **xlim** (*tuple**[**float**, **float**]*) – The energy axis limits. Defaults to None for automatic
    determination.


    * **ylim** (*tuple**[**float**, **float**]*) – The y-axis limits. Defaults to None for automatic determination.


    * **invert_axes** (*bool*) – Whether to invert the x and y axes. Enables chemist style DOS plotting.
    Defaults to False.


    * **beta_dashed** (*bool*) – Plots the beta spin channel with a dashed line. Defaults to False.



#### save_plot(filename, img_format='eps', xlim=None, ylim=None, invert_axes=False, beta_dashed=False)
Save matplotlib plot to a file.


* **Parameters**


    * **filename** – Filename to write to.


    * **img_format** – Image format to use. Defaults to EPS.


    * **xlim** – Specifies the x-axis limits. Set to None for automatic
    determination.


    * **ylim** – Specifies the y-axis limits.


    * **invert_axes** (*bool*) – Whether to invert the x and y axes. Enables chemist style DOS plotting.
    Defaults to False.


    * **beta_dashed** (*bool*) – Plots the beta spin channel with a dashed line. Defaults to False.



#### show(xlim=None, ylim=None, invert_axes=False, beta_dashed=False)
Show the plot using matplotlib.


* **Parameters**


    * **xlim** – Specifies the x-axis limits. Set to None for automatic
    determination.


    * **ylim** – Specifies the y-axis limits.


    * **invert_axes** (*bool*) – Whether to invert the x and y axes. Enables chemist style DOS plotting.
    Defaults to False.


    * **beta_dashed** (*bool*) – Plots the beta spin channel with a dashed line. Defaults to False.



### pymatgen.electronic_structure.plotter.fold_point(p, lattice, coords_are_cartesian=False)
Folds a point with coordinates p inside the first Brillouin zone of the lattice.


* **Parameters**


    * **p** – coordinates of one point


    * **lattice** – Lattice object used to convert from reciprocal to Cartesian coordinates


    * **coords_are_cartesian** – Set to True if you are providing
    coordinates in Cartesian coordinates. Defaults to False.



* **Returns**

    The Cartesian coordinates folded inside the first Brillouin zone



### pymatgen.electronic_structure.plotter.plot_brillouin_zone(bz_lattice, lines=None, labels=None, kpoints=None, fold=False, coords_are_cartesian=False, ax=None, \*\*kwargs)
Plots a 3D representation of the Brillouin zone of the structure.
Can add to the plot paths, labels and kpoints.


* **Parameters**


    * **bz_lattice** – Lattice object of the Brillouin zone


    * **lines** – list of lists of coordinates. Each list represent a different path


    * **labels** – dict containing the label as a key and the coordinates as value.


    * **kpoints** – list of coordinates


    * **fold** – whether the points should be folded inside the first Brillouin Zone.
    Defaults to False. Requires lattice if True.


    * **coords_are_cartesian** – Set to True if you are providing
    coordinates in Cartesian coordinates. Defaults to False.


    * **ax** – matplotlib `Axes` or None if a new figure should be created.


    * **kwargs** – provided by add_fig_kwargs decorator



* **Returns**

    matplotlib figure

    Keyword arguments controlling the display of the figure:

    | kwargs

     | Meaning

     |
    | ------------ | ------------------------------------------------------------------------------------------------- |  |  |  |  |
    | title

            | Title of the plot (Default: None).

                                                                    |
    | show

             | True to show the figure (default: True).

                                                              |
    | savefig

          | ”abc.png” or “abc.eps” to save the figure to a file.

                                                  |
    | size_kwargs

      | Dictionary with options passed to fig.set_size_inches
    e.g. size_kwargs=dict(w=3, h=4)

                 |
    | tight_layout

     | True to call fig.tight_layout (default: False)

                                                        |
    | ax_grid

          | True (False) to add (remove) grid from all axes in fig.
    Default: None i.e. fig is left unchanged.

     |
    | ax_annotate

      | Add labels to  subplots e.g. (a), (b).
    Default: False

                                                 |
    | fig_close

        | Close figure. Default: False.

                                                                         |



### pymatgen.electronic_structure.plotter.plot_brillouin_zone_from_kpath(kpath, ax=None, \*\*kwargs)
Gives the plot (as a matplotlib object) of the symmetry line path in

    the Brillouin Zone.


* **Parameters**


    * **kpath** ([*HighSymmKpath*](pymatgen.symmetry.bandstructure.md#pymatgen.symmetry.bandstructure.HighSymmKpath)) – a HighSymmKPath object


    * **ax** – matplotlib `Axes` or None if a new figure should be created.


    * **\*\*kwargs** – provided by add_fig_kwargs decorator



* **Returns**

    matplotlib figure

    Keyword arguments controlling the display of the figure:

    | kwargs

           | Meaning

                                                                                               |
    | ------------ | ------------------------------------------------------------------------------------------------- |
    | title

            | Title of the plot (Default: None).

                                                                    |
    | show

             | True to show the figure (default: True).

                                                              |
    | savefig

          | ”abc.png” or “abc.eps” to save the figure to a file.

                                                  |
    | size_kwargs

      | Dictionary with options passed to fig.set_size_inches
    e.g. size_kwargs=dict(w=3, h=4)

                 |
    | tight_layout

     | True to call fig.tight_layout (default: False)

                                                        |
    | ax_grid

          | True (False) to add (remove) grid from all axes in fig.
    Default: None i.e. fig is left unchanged.

     |
    | ax_annotate

      | Add labels to  subplots e.g. (a), (b).
    Default: False

                                                 |
    | fig_close

        | Close figure. Default: False.

                                                                         |



### pymatgen.electronic_structure.plotter.plot_ellipsoid(hessian, center, lattice=None, rescale=1.0, ax=None, coords_are_cartesian=False, arrows=False, \*\*kwargs)
Plots a 3D ellipsoid rappresenting the Hessian matrix in input.
Useful to get a graphical visualization of the effective mass
of a band in a single k-point.


* **Parameters**


    * **hessian** – the Hessian matrix


    * **center** – the center of the ellipsoid in reciprocal coords (Default)


    * **lattice** – Lattice object of the Brillouin zone


    * **rescale** – factor for size scaling of the ellipsoid


    * **ax** – matplotlib `Axes` or None if a new figure should be created.


    * **coords_are_cartesian** – Set to True if you are providing a center in
    Cartesian coordinates. Defaults to False.


    * **arrows** – whether to plot arrows for the principal axes of the ellipsoid. Defaults to False.


    * **\*\*kwargs** – passed to the matplotlib function ‘plot_wireframe’.
    Color defaults to blue, rstride and cstride
    default to 4, alpha defaults to 0.2.



* **Returns**

    matplotlib figure and matplotlib ax


Example of use:

    fig,ax=plot_wigner_seitz(struct.reciprocal_lattice)
    plot_ellipsoid(hessian,[0.0,0.0,0.0], struct.reciprocal_lattice,ax=ax)


### pymatgen.electronic_structure.plotter.plot_fermi_surface(data, structure, cbm, energy_levels=None, multiple_figure=True, mlab_figure=None, kpoints_dict=None, colors=None, transparency_factor=None, labels_scale_factor=0.05, points_scale_factor=0.02, interactive=True)
Plot the Fermi surface at specific energy value using Boltztrap 1 FERMI
mode.

The easiest way to use this plotter is:

>
> 1. Run boltztrap in ‘FERMI’ mode using BoltztrapRunner,


> 2. Load BoltztrapAnalyzer using your method of choice (e.g., from_files)


> 3. Pass in your BoltztrapAnalyzer’s fermi_surface_data as this

>     function’s data argument.


* **Parameters**


    * **data** – energy values in a 3D grid from a CUBE file via read_cube_file
    function, or from a BoltztrapAnalyzer.fermi_surface_data


    * **structure** – structure object of the material


    * **energy_levels** (*[**float**]*) – Energy values for plotting the fermi surface(s)
    By default 0 eV correspond to the VBM, as in the plot of band
    structure along symmetry line.
    Default: One surface, with max energy value + 0.01 eV


    * **cbm** (*bool*) – Boolean value to specify if the considered band is a
    conduction band or not


    * **multiple_figure** (*bool*) – If True a figure for each energy level will be
    shown. If False all the surfaces will be shown in the same figure.
    In this last case, tune the transparency factor.


    * **mlab_figure** (*mayavi.mlab.figure*) – A previous figure to plot a new
    surface on.


    * **kpoints_dict** (*dict*) – dictionary of kpoints to label in the plot.
    Example: {“K”:[0.5,0.0,0.5]}, coords are fractional


    * **colors** (*[**tuple**]*) – Iterable of 3-tuples (r,g,b) of integers to define
    the colors of each surface (one per energy level).
    Should be the same length as the number of surfaces being plotted.
    Example (3 surfaces): colors=[(1,0,0), (0,1,0), (0,0,1)]
    Example (2 surfaces): colors=[(0, 0.5, 0.5)]


    * **transparency_factor** (*float*) – Values in the range [0,1] to tune the
    opacity of each surface. Should be one transparency_factor per
    surface.


    * **labels_scale_factor** (*float*) – factor to tune size of the kpoint labels


    * **points_scale_factor** (*float*) – factor to tune size of the kpoint points


    * **interactive** (*bool*) – if True an interactive figure will be shown.
    If False a non interactive figure will be shown, but it is possible
    to plot other surfaces on the same figure. To make it interactive,
    run mlab.show().



* **Returns**

    The mlab plotter and an interactive

        figure to control the plot.




* **Return type**

    ((mayavi.mlab.figure, mayavi.mlab))


Note: Experimental.

    Please, double check the surface shown by using some
    other software and report issues.


### pymatgen.electronic_structure.plotter.plot_labels(labels, lattice=None, coords_are_cartesian=False, ax=None, \*\*kwargs)
Adds labels to a matplotlib Axes.


* **Parameters**


    * **labels** – dict containing the label as a key and the coordinates as value.


    * **lattice** – Lattice object used to convert from reciprocal to Cartesian coordinates


    * **coords_are_cartesian** – Set to True if you are providing.
    coordinates in Cartesian coordinates. Defaults to False.
    Requires lattice if False.


    * **ax** – matplotlib `Axes` or None if a new figure should be created.


    * **kwargs** – kwargs passed to the matplotlib function ‘text’. Color defaults to blue
    and size to 25.



* **Returns**

    matplotlib figure and matplotlib ax



### pymatgen.electronic_structure.plotter.plot_lattice_vectors(lattice, ax=None, \*\*kwargs)
Adds the basis vectors of the lattice provided to a matplotlib Axes.


* **Parameters**


    * **lattice** – Lattice object


    * **ax** – matplotlib `Axes` or None if a new figure should be created.


    * **kwargs** – kwargs passed to the matplotlib function ‘plot’. Color defaults to green
    and linewidth to 3.



* **Returns**

    matplotlib figure and matplotlib ax



### pymatgen.electronic_structure.plotter.plot_path(line, lattice=None, coords_are_cartesian=False, ax=None, \*\*kwargs)
Adds a line passing through the coordinates listed in ‘line’ to a matplotlib Axes.


* **Parameters**


    * **line** – list of coordinates.


    * **lattice** – Lattice object used to convert from reciprocal to Cartesian coordinates


    * **coords_are_cartesian** – Set to True if you are providing
    coordinates in Cartesian coordinates. Defaults to False.
    Requires lattice if False.


    * **ax** – matplotlib `Axes` or None if a new figure should be created.


    * **kwargs** – kwargs passed to the matplotlib function ‘plot’. Color defaults to red
    and linewidth to 3.



* **Returns**

    matplotlib figure and matplotlib ax



### pymatgen.electronic_structure.plotter.plot_points(points, lattice=None, coords_are_cartesian=False, fold=False, ax=None, \*\*kwargs)
Adds Points to a matplotlib Axes.


* **Parameters**


    * **points** – list of coordinates


    * **lattice** – Lattice object used to convert from reciprocal to Cartesian coordinates


    * **coords_are_cartesian** – Set to True if you are providing
    coordinates in Cartesian coordinates. Defaults to False.
    Requires lattice if False.


    * **fold** – whether the points should be folded inside the first Brillouin Zone.
    Defaults to False. Requires lattice if True.


    * **ax** – matplotlib `Axes` or None if a new figure should be created.


    * **kwargs** – kwargs passed to the matplotlib function ‘scatter’. Color defaults to blue



* **Returns**

    matplotlib figure and matplotlib ax



### pymatgen.electronic_structure.plotter.plot_wigner_seitz(lattice, ax=None, \*\*kwargs)
Adds the skeleton of the Wigner-Seitz cell of the lattice to a matplotlib Axes.


* **Parameters**


    * **lattice** – Lattice object


    * **ax** – matplotlib `Axes` or None if a new figure should be created.


    * **kwargs** – kwargs passed to the matplotlib function ‘plot’. Color defaults to black
    and linewidth to 1.



* **Returns**

    matplotlib figure and matplotlib ax