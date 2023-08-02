---
layout: default
title: pymatgen.phonon.plotter.md
nav_exclude: true
---

# pymatgen.phonon.plotter module

This module implements plotter for DOS and band structure.


### _class_ pymatgen.phonon.plotter.FreqUnits(factor, label)
Bases: `tuple`

Create new instance of FreqUnits(factor, label)


#### factor()
Alias for field number 0


#### label()
Alias for field number 1


### _class_ pymatgen.phonon.plotter.GruneisenPhononBSPlotter(bs)
Bases: `PhononBSPlotter`

Class to plot or get data to facilitate the plot of band structure objects.


* **Parameters**

    **bs** – A GruneisenPhononBandStructureSymmLine object.



#### bs_plot_data()
Get the data nicely formatted for a plot.


* **Returns**

    ticks: A dict with the ‘distances’ at which there is a qpoint (the
    x axis) and the labels (None if no label)
    frequencies: A list (one element for each branch) of frequencies for
    each qpoint: [branch][qpoint][mode]. The data is
    stored by branch to facilitate the plotting
    gruneisen: GruneisenPhononBandStructureSymmLine
    lattice: The reciprocal lattice.



* **Return type**

    A dict of the following format



#### get_plot_gs(ylim=None)
Get a matplotlib object for the gruneisen bandstructure plot.


* **Parameters**

    **ylim** – Specify the y-axis (gruneisen) limits; by default None let
    the code choose.



#### plot_compare_gs(other_plotter: GruneisenPhononBSPlotter)
plot two band structure for comparison. One is in red the other in blue.
The two band structures need to be defined on the same symmetry lines!
and the distance between symmetry lines is
the one of the band structure used to build the PhononBSPlotter.


* **Parameters**

    **other_plotter** (*GruneisenPhononBSPlotter*) – another phonon DOS plotter defined along
    the same symmetry lines.



* **Returns**

    a matplotlib object with both band structures



#### save_plot_gs(filename, img_format='eps', ylim=None)
Save matplotlib plot to a file.


* **Parameters**


    * **filename** – Filename to write to.


    * **img_format** – Image format to use. Defaults to EPS.


    * **ylim** – Specifies the y-axis limits.



#### show_gs(ylim=None)
Show the plot using matplotlib.


* **Parameters**

    **ylim** – Specifies the y-axis limits.



### _class_ pymatgen.phonon.plotter.GruneisenPlotter(gruneisen)
Bases: `object`

Class to plot Gruneisenparameter Object.

Class to plot information from Gruneisenparameter Object
:param gruneisen: GruneisenParameter Object.


#### get_plot(marker='o', markersize=None, units='thz')
will produce a plot
:param marker: marker for the depiction
:param markersize: size of the marker
:param units: unit for the plots, accepted units: thz, ev, mev, ha, cm-1, cm^-1.

Returns: plot


#### save_plot(filename, img_format='pdf', units='thz')
Will save the plot to a file
:param filename: name of the filename
:param img_format: format of the saved plot
:param units: accepted units: thz, ev, mev, ha, cm-1, cm^-1.


#### show(units='thz')
will show the plot
:param units: units for the plot, accepted units: thz, ev, mev, ha, cm-1, cm^-1.

Returns: plot


### _class_ pymatgen.phonon.plotter.PhononBSPlotter(bs)
Bases: `object`

Class to plot or get data to facilitate the plot of band structure objects.


* **Parameters**

    **bs** – A PhononBandStructureSymmLine object.



#### bs_plot_data()
Get the data nicely formatted for a plot.


* **Returns**

    ticks: A dict with the ‘distances’ at which there is a qpoint (the
    x axis) and the labels (None if no label)
    frequencies: A list (one element for each branch) of frequencies for
    each qpoint: [branch][qpoint][mode]. The data is
    stored by branch to facilitate the plotting
    lattice: The reciprocal lattice.



* **Return type**

    A dict of the following format



#### get_plot(ylim=None, units='thz')
Get a matplotlib object for the bandstructure plot.


* **Parameters**


    * **ylim** – Specify the y-axis (frequency) limits; by default None let
    the code choose.


    * **units** – units for the frequencies. Accepted values thz, ev, mev, ha, cm-1, cm^-1.



#### get_proj_plot(site_comb: str | list[list[int]] = 'element', ylim: tuple[None | float, None | float] | None = None, units: str = 'thz', rgb_labels: tuple[None | str] | None = None)
Get a matplotlib object for the bandstructure plot projected along atomic
sites.


* **Parameters**


    * **site_comb** – a list of list, for example, [[0],[1],[2,3,4]];
    the numbers in each sublist represents the indices of atoms;
    the atoms in a same sublist will be plotted in a same color;
    if not specified, unique elements are automatically grouped.


    * **ylim** – Specify the y-axis (frequency) limits; by default None let
    the code choose.


    * **units** – units for the frequencies. Accepted values thz, ev, mev, ha, cm-1, cm^-1.
    Defaults to “thz”.


    * **rgb_labels** – a list of rgb colors for the labels; if not specified,
    the colors will be automatically generated.



#### get_ticks()
Get all ticks and labels for a band structure plot.


* **Returns**

    a list of distance at which ticks should
    be set and ‘label’: a list of label for each of those ticks.



* **Return type**

    A dict with ‘distance’



#### plot_brillouin()
Plot the Brillouin zone.


#### plot_compare(other_plotter, units='thz')
plot two band structure for comparison. One is in red the other in blue.
The two band structures need to be defined on the same symmetry lines!
and the distance between symmetry lines is the one of the band structure
used to build the PhononBSPlotter.


* **Parameters**


    * **other_plotter** – another PhononBSPlotter object defined along the same symmetry lines


    * **units** – units for the frequencies. Accepted values thz, ev, mev, ha, cm-1, cm^-1.
    Defaults to ‘thz’.



* **Returns**

    a matplotlib object with both band structures



#### save_plot(filename, img_format='eps', ylim=None, units='thz')
Save matplotlib plot to a file.


* **Parameters**


    * **filename** – Filename to write to.


    * **img_format** – Image format to use. Defaults to EPS.


    * **ylim** – Specifies the y-axis limits.


    * **units** – units for the frequencies. Accepted values thz, ev, mev, ha, cm-1, cm^-1.



#### show(ylim=None, units='thz')
Show the plot using matplotlib.


* **Parameters**


    * **ylim** – Specify the y-axis (frequency) limits; by default None let
    the code choose.


    * **units** – units for the frequencies. Accepted values thz, ev, mev, ha, cm-1, cm^-1.



#### show_proj(site_comb: str | list[list[int]] = 'element', ylim: tuple[None | float, None | float] | None = None, units: str = 'thz', rgb_labels: tuple[str] | None = None)
Show the projected plot using matplotlib.


* **Parameters**


    * **site_comb** – A list of list of indices of sites to combine. For example,
    [[0, 1], [2, 3]] will combine the projections of sites 0 and 1,
    and sites 2 and 3. Defaults to “element”, which will combine
    sites by element.


    * **ylim** – Specify the y-axis (frequency) limits; by default None let
    the code choose.


    * **units** – units for the frequencies. Accepted values thz, ev, mev, ha, cm-1, cm^-1.
    Defaults to “thz”.


    * **rgb_labels** – A list of labels for the rgb triangle. Defaults to None,
    which will use the element symbols.



### _class_ pymatgen.phonon.plotter.PhononDosPlotter(stack=False, sigma=None)
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


* **Parameters**


    * **stack** – Whether to plot the DOS as a stacked area graph


    * **sigma** – A float specifying a standard deviation for Gaussian smearing
    the DOS for nicer looking plots. Defaults to None for no smearing.



#### add_dos(label, dos)
Adds a dos for plotting.


* **Parameters**


    * **label** – label for the DOS. Must be unique.


    * **dos** – PhononDos object



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

    {‘frequencies’:..,
    ‘densities’: …}}



* **Return type**

    Dict of dos data. Generally of the form, {label



#### get_plot(xlim=None, ylim=None, units='thz')
Get a matplotlib plot showing the DOS.


* **Parameters**


    * **xlim** – Specifies the x-axis limits. Set to None for automatic
    determination.


    * **ylim** – Specifies the y-axis limits.


    * **units** – units for the frequencies. Accepted values thz, ev, mev, ha, cm-1, cm^-1.



#### save_plot(filename, img_format='eps', xlim=None, ylim=None, units='thz')
Save matplotlib plot to a file.


* **Parameters**


    * **filename** – Filename to write to.


    * **img_format** – Image format to use. Defaults to EPS.


    * **xlim** – Specifies the x-axis limits. Set to None for automatic
    determination.


    * **ylim** – Specifies the y-axis limits.


    * **units** – units for the frequencies. Accepted values thz, ev, mev, ha, cm-1, cm^-1



#### show(xlim=None, ylim=None, units='thz')
Show the plot using matplotlib.


* **Parameters**


    * **xlim** – Specifies the x-axis limits. Set to None for automatic
    determination.


    * **ylim** – Specifies the y-axis limits.


    * **units** – units for the frequencies. Accepted values thz, ev, mev, ha, cm-1, cm^-1.



### _class_ pymatgen.phonon.plotter.ThermoPlotter(dos, structure=None)
Bases: `object`

Plotter for thermodynamic properties obtained from phonon DOS.
If the structure corresponding to the DOS, it will be used to extract the formula unit and provide
the plots in units of mol instead of mole-cell.


* **Parameters**


    * **dos** – A PhononDos object.


    * **structure** – A Structure object corresponding to the structure used for the calculation.



#### plot_cv(tmin, tmax, ntemp, ylim=None, \*\*kwargs)
Plots the constant volume specific heat C_v in a temperature range.


* **Parameters**


    * **tmin** – minimum temperature


    * **tmax** – maximum temperature


    * **ntemp** – number of steps


    * **ylim** – tuple specifying the y-axis limits.


    * **kwargs** – kwargs passed to the matplotlib function ‘plot’.



* **Returns**

    matplotlib figure


Keyword arguments controlling the display of the figure:

| kwargs

 | Meaning

 |
| ------------ | ------------------------------------------------------------------------------------------------- |  |  |  |  |  |  |  |  |  |  |  |  |
| title

        | Title of the plot (Default: None).

                                                                |
| show

         | True to show the figure (default: True).

                                                          |
| savefig

      | “abc.png” or “abc.eps” to save the figure to a file.

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

#### plot_entropy(tmin, tmax, ntemp, ylim=None, \*\*kwargs)
Plots the vibrational entrpy in a temperature range.


* **Parameters**


    * **tmin** – minimum temperature


    * **tmax** – maximum temperature


    * **ntemp** – number of steps


    * **ylim** – tuple specifying the y-axis limits.


    * **kwargs** – kwargs passed to the matplotlib function ‘plot’.



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

      | “abc.png” or “abc.eps” to save the figure to a file.

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

#### plot_helmholtz_free_energy(tmin, tmax, ntemp, ylim=None, \*\*kwargs)
Plots the vibrational contribution to the Helmoltz free energy in a temperature range.


* **Parameters**


    * **tmin** – minimum temperature


    * **tmax** – maximum temperature


    * **ntemp** – number of steps


    * **ylim** – tuple specifying the y-axis limits.


    * **kwargs** – kwargs passed to the matplotlib function ‘plot’.



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

      | “abc.png” or “abc.eps” to save the figure to a file.

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

#### plot_internal_energy(tmin, tmax, ntemp, ylim=None, \*\*kwargs)
Plots the vibrational internal energy in a temperature range.


* **Parameters**


    * **tmin** – minimum temperature


    * **tmax** – maximum temperature


    * **ntemp** – number of steps


    * **ylim** – tuple specifying the y-axis limits.


    * **kwargs** – kwargs passed to the matplotlib function ‘plot’.



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

      | “abc.png” or “abc.eps” to save the figure to a file.

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

#### plot_thermodynamic_properties(tmin, tmax, ntemp, ylim=None, \*\*kwargs)
Plots all the thermodynamic properties in a temperature range.


* **Parameters**


    * **tmin** – minimum temperature


    * **tmax** – maximum temperature


    * **ntemp** – number of steps


    * **ylim** – tuple specifying the y-axis limits.


    * **kwargs** – kwargs passed to the matplotlib function ‘plot’.



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

      | “abc.png” or “abc.eps” to save the figure to a file.

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

### pymatgen.phonon.plotter.freq_units(units)

* **Parameters**

    **units** – str, accepted values: thz, ev, mev, ha, cm-1, cm^-1.



* **Returns**

    Returns conversion factor from THz to the required units and the label in the form of a namedtuple