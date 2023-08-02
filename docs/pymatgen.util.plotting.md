---
layout: default
title: pymatgen.util.plotting.md
nav_exclude: true
---

# pymatgen.util.plotting module

Utilities for generating nicer plots.


### pymatgen.util.plotting.add_fig_kwargs(func)
Decorator that adds keyword arguments for functions returning matplotlib
figures.

The function should return either a matplotlib figure or None to signal
some sort of error/unexpected event.
See doc string below for the list of supported options.


### pymatgen.util.plotting.format_formula(formula)
Converts str of chemical formula into
latex format for labelling purposes.


* **Parameters**

    **formula** (*str*) – Chemical formula



### pymatgen.util.plotting.get_ax3d_fig_plt(ax=None, \*\*kwargs)
Helper function used in plot functions supporting an optional Axes3D
argument. If ax is None, we build the matplotlib figure and create the
Axes3D else we return the current active figure.


* **Parameters**


    * **ax** (*Axes3D**, **optional*) – Axes3D object. Defaults to None.


    * **kwargs** – keyword arguments are passed to plt.figure if ax is not None.



* **Returns**

    `Axes` object
    figure: matplotlib figure
    plt: matplotlib pyplot module.



* **Return type**

    ax



### pymatgen.util.plotting.get_ax_fig_plt(ax=None, \*\*kwargs)
Helper function used in plot functions supporting an optional Axes argument.
If ax is None, we build the matplotlib figure and create the Axes else
we return the current active figure.


* **Parameters**


    * **ax** (*Axes**, **optional*) – Axes object. Defaults to None.


    * **kwargs** – keyword arguments are passed to plt.figure if ax is not None.



* **Returns**

    `Axes` object
    figure: matplotlib figure
    plt: matplotlib pyplot module.



* **Return type**

    ax



### pymatgen.util.plotting.get_axarray_fig_plt(ax_array, nrows=1, ncols=1, sharex=False, sharey=False, squeeze=True, subplot_kw=None, gridspec_kw=None, \*\*fig_kw)
Helper function used in plot functions that accept an optional array of Axes
as argument. If ax_array is None, we build the matplotlib figure and
create the array of Axes by calling plt.subplots else we return the
current active figure.


* **Returns**

    Array of `Axes` objects
    figure: matplotlib figure
    plt: matplotlib pyplot module.



* **Return type**

    ax



### pymatgen.util.plotting.periodic_table_heatmap(elemental_data=None, cbar_label='', cbar_label_size=14, show_plot=False, cmap='YlOrRd', cmap_range=None, blank_color='grey', edge_color='white', value_format=None, value_fontsize=10, symbol_fontsize=14, max_row=9, readable_fontcolor=False, pymatviz: bool = True, \*\*kwargs)
A static method that generates a heat map overlaid on a periodic table.


* **Parameters**


    * **elemental_data** (*dict*) – A dictionary with the element as a key and a
    value assigned to it, e.g. surface energy and frequency, etc.
    Elements missing in the elemental_data will be grey by default
    in the final table elemental_data={“Fe”: 4.2, “O”: 5.0}.


    * **cbar_label** (*str*) – Label of the color bar. Default is “”.


    * **cbar_label_size** (*float*) – Font size for the color bar label. Default is 14.


    * **cmap_range** (*tuple*) – Minimum and maximum value of the color map scale.
    If None, the color map will automatically scale to the range of the
    data.


    * **show_plot** (*bool*) – Whether to show the heatmap. Default is False.


    * **value_format** (*str*) – Formatting string to show values. If None, no value
    is shown. Example: “%.4f” shows float to four decimals.


    * **value_fontsize** (*float*) – Font size for values. Default is 10.


    * **symbol_fontsize** (*float*) – Font size for element symbols. Default is 14.


    * **cmap** (*str*) – Color scheme of the heatmap. Default is ‘YlOrRd’.
    Refer to the matplotlib documentation for other options.


    * **blank_color** (*str*) – Color assigned for the missing elements in
    elemental_data. Default is “grey”.


    * **edge_color** (*str*) – Color assigned for the edge of elements in the
    periodic table. Default is “white”.


    * **max_row** (*int*) – Maximum number of rows of the periodic table to be
    shown. Default is 9, which means the periodic table heat map covers
    the standard 7 rows of the periodic table + 2 rows for the lanthanides
    and actinides. Use a value of max_row = 7 to exclude the lanthanides and
    actinides.


    * **readable_fontcolor** (*bool*) – Whether to use readable font color depending
    on background color. Default is False.


    * **pymatviz** (*bool*) – Whether to use pymatviz to generate the heatmap. Defaults to True.
    See [https://github.com/janosh/pymatviz](https://github.com/janosh/pymatviz).


    * **kwargs** – Passed to pymatviz.ptable_heatmap_plotly



### pymatgen.util.plotting.pretty_plot(width=8, height=None, plt=None, dpi=None, color_cycle=('qualitative', 'Set1_9'))
Provides a publication quality plot, with nice defaults for font sizes etc.


* **Parameters**


    * **width** (*float*) – Width of plot in inches. Defaults to 8in.


    * **height** (*float*) – Height of plot in inches. Defaults to width \* golden
    ratio.


    * **plt** (*matplotlib.pyplot*) – If plt is supplied, changes will be made to an
    existing plot. Otherwise, a new plot will be created.


    * **dpi** (*int*) – Sets dot per inch for figure. Defaults to 300.


    * **color_cycle** (*tuple*) – Set the color cycle for new plots to one of the
    color sets in palettable. Defaults to a qualitative Set1_9.



* **Returns**

    Matplotlib plot object with properly sized fonts.



### pymatgen.util.plotting.pretty_plot_two_axis(x, y1, y2, xlabel=None, y1label=None, y2label=None, width=8, height=None, dpi=300, \*\*plot_kwargs)
Variant of pretty_plot that does a dual axis plot. Adapted from matplotlib
examples. Makes it easier to create plots with different axes.


* **Parameters**


    * **x** (*np.ndarray/list*) – Data for x-axis.


    * **y1** (*dict/np.ndarray/list*) – Data for y1 axis (left). If a dict, it will
    be interpreted as a {label: sequence}.


    * **y2** (*dict/np.ndarray/list*) – Data for y2 axis (right). If a dict, it will
    be interpreted as a {label: sequence}.


    * **xlabel** (*str*) – If not None, this will be the label for the x-axis.


    * **y1label** (*str*) – If not None, this will be the label for the y1-axis.


    * **y2label** (*str*) – If not None, this will be the label for the y2-axis.


    * **width** (*float*) – Width of plot in inches. Defaults to 8in.


    * **height** (*float*) – Height of plot in inches. Defaults to width \* golden
    ratio.


    * **dpi** (*int*) – Sets dot per inch for figure. Defaults to 300.


    * **plot_kwargs** – Passthrough kwargs to matplotlib’s plot method. E.g.,
    linewidth, etc.



* **Returns**

    matplotlib.pyplot



### pymatgen.util.plotting.pretty_polyfit_plot(x, y, deg=1, xlabel=None, ylabel=None, \*\*kwargs)
Convenience method to plot data with trend lines based on polynomial fit.


* **Parameters**


    * **x** – Sequence of x data.


    * **y** – Sequence of y data.


    * **deg** (*int*) – Degree of polynomial. Defaults to 1.


    * **xlabel** (*str*) – Label for x-axis.


    * **ylabel** (*str*) – Label for y-axis.


    * **kwargs** – Keyword args passed to pretty_plot.



* **Returns**

    matplotlib.pyplot object.



### pymatgen.util.plotting.van_arkel_triangle(list_of_materials, annotate=True)
A static method that generates a binary van Arkel-Ketelaar triangle to

    quantify the ionic, metallic and covalent character of a compound
    by plotting the electronegativity difference (y) vs average (x).
    See:

    > A.E. van Arkel, Molecules and Crystals in Inorganic Chemistry,

    >     Interscience, New York (1956)

    and

        J.A.A Ketelaar, Chemical Constitution (2nd edition), An Introduction

            to the Theory of the Chemical Bond, Elsevier, New York (1958).


* **Parameters**


    * **list_of_materials** (*list*) – A list of computed entries of binary
    materials or a list of lists containing two elements (str).


    * **annotate** (*bool*) – Whether or not to label the points on the
    triangle with reduced formula (if list of entries) or pair
    of elements (if list of list of str).