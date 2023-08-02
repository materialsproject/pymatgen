---
layout: default
title: pymatgen.vis.plotters.md
nav_exclude: true
---

# pymatgen.vis.plotters module

This module defines generic plotters.


### _class_ pymatgen.vis.plotters.SpectrumPlotter(xshift=0.0, yshift=0.0, stack=False, color_cycle=('qualitative', 'Set1_9'))
Bases: `object`

Class for plotting Spectrum objects and subclasses. Note that the interface
is extremely flexible given that there are many different ways in which
people want to view spectra. The typical usage is:

```default
# Initializes plotter with some optional args. Defaults are usually
# fine,
plotter = SpectrumPlotter()

# Adds a DOS (A kind of spectra) with a label.
plotter.add_spectrum("Total DOS", dos)

# Alternatively, you can add a dict of DOSs. This is the typical
# form returned by CompleteDos.get_spd/element/others_dos().
plotter.add_spectra({"dos1": dos1, "dos2": dos2})
```


* **Parameters**


    * **xshift** (*float*) – A shift that is applied to the x values. This is
    commonly used to shift to an arbitrary zero. E.g., zeroing at the
    Fermi energy in DOS, or at the absorption edge in XAS spectra. The
    same xshift is applied to all spectra.


    * **yshift** (*float*) – A shift that is applied to the y values. This is
    commonly used to displace spectra for easier visualization.
    Successive spectra are applied successive shifts.


    * **stack** (*bool*) – Whether to stack plots rather than simply plot them.
    For example, DOS plot can usually be stacked to look at the
    contribution of each orbital.


    * **color_cycle** (*str*) – Default color cycle to use. Note that this can be
    overridden.



#### add_spectra(spectra_dict, key_sort_func=None)
Add a dictionary of doses, with an optional sorting function for the
keys.


* **Parameters**


    * **dos_dict** – dict of {label: Dos}


    * **key_sort_func** – function used to sort the dos_dict keys.



#### add_spectrum(label, spectrum, color=None)
Adds a Spectrum for plotting.


* **Parameters**


    * **label** (*str*) – Label for the Spectrum. Must be unique.


    * **spectrum** – Spectrum object


    * **color** (*str*) – This is passed on to matplotlib. E.g., “k–” indicates
    a dashed black line. If None, a color will be chosen based on
    the default color cycle.



#### get_plot(xlim=None, ylim=None)
Get a matplotlib plot showing the DOS.


* **Parameters**


    * **xlim** – Specifies the x-axis limits. Set to None for automatic
    determination.


    * **ylim** – Specifies the y-axis limits.



#### save_plot(filename, img_format='eps', \*\*kwargs)
Save matplotlib plot to a file.


* **Parameters**


    * **filename** – Filename to write to.


    * **img_format** – Image format to use. Defaults to EPS.



#### show(\*\*kwargs)
Show the plot using matplotlib.