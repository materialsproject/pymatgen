---
layout: default
title: pymatgen.apps.battery.plotter.md
nav_exclude: true
---

# pymatgen.apps.battery.plotter module

This module provides plotting capabilities for battery related applications.


### _class_ pymatgen.apps.battery.plotter.VoltageProfilePlotter(xaxis='capacity', hide_negative=False)
Bases: `object`

A plotter to make voltage profile plots for batteries.


* **Parameters**


    * **xaxis** – The quantity to use as the xaxis. Can be either


    * **capacity_grav** (*-*) – the graviometric capcity


    * **capacity_vol** (*-*) – the volumetric capacity


    * **x_form** (*-*) – the number of working ions per formula unit of the host


    * **frac_x** (*-*) – the atomic fraction of the working ion


    * **hide_negative** – If True only plot the voltage steps above zero.



#### add_electrode(electrode, label=None)
Add an electrode to the plot.


* **Parameters**


    * **electrode** – An electrode. All electrodes satisfying the
    AbstractElectrode interface should work.


    * **label** – A label for the electrode. If None, defaults to a counting
    system, i.e. ‘Electrode 1’, ‘Electrode 2’, …



#### get_plot(width=8, height=8, term_zero=True)
Returns a plot object.


* **Parameters**


    * **width** – Width of the plot. Defaults to 8 in.


    * **height** – Height of the plot. Defaults to 6 in.


    * **term_zero** – If True append zero voltage point at the end



* **Returns**

    A matplotlib plot object.



#### get_plot_data(electrode, term_zero=True)

* **Parameters**


    * **electrode** – Electrode object


    * **term_zero** – If True append zero voltage point at the end.



* **Returns**

    Plot data in x, y.



#### get_plotly_figure(width=800, height=600, font_dict=None, term_zero=True, \*\*kwargs)
Return plotly Figure object.


* **Parameters**


    * **width** – Width of the plot. Defaults to 800 px.


    * **height** – Height of the plot. Defaults to 600 px.


    * **font_dict** – define the font. Defaults to {“family”: “Arial”, “size”: 24, “color”: “#000000”}


    * **term_zero** – If True append zero voltage point at the end


    * **\*\*kwargs** – passed to plotly.graph_objects.Layout



#### save(filename, image_format='eps', width=8, height=6)
Save the plot to an image file.


* **Parameters**


    * **filename** – Filename to save to.


    * **image_format** – Format to save to. Defaults to eps.


    * **width** – Width of the plot. Defaults to 8 in.


    * **height** – Height of the plot. Defaults to 6 in.



#### show(width=8, height=6)
Show the voltage profile plot.


* **Parameters**


    * **width** – Width of the plot. Defaults to 8 in.


    * **height** – Height of the plot. Defaults to 6 in.