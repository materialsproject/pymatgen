# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
"""
Utilities for generating nicer plots.
"""
import math
import numpy as np

from pymatgen.core.periodic_table import Element


__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Mar 13, 2012"


def pretty_plot(width=8, height=None, plt=None, dpi=None,
                color_cycle=("qualitative", "Set1_9")):
    """
    Provides a publication quality plot, with nice defaults for font sizes etc.

    Args:
        width (float): Width of plot in inches. Defaults to 8in.
        height (float): Height of plot in inches. Defaults to width * golden
            ratio.
        plt (matplotlib.pyplot): If plt is supplied, changes will be made to an
            existing plot. Otherwise, a new plot will be created.
        dpi (int): Sets dot per inch for figure. Defaults to 300.
        color_cycle (tuple): Set the color cycle for new plots to one of the
            color sets in palettable. Defaults to a qualitative Set1_9.

    Returns:
        Matplotlib plot object with properly sized fonts.
    """
    ticksize = int(width * 2.5)

    golden_ratio = (math.sqrt(5) - 1) / 2

    if not height:
        height = int(width * golden_ratio)

    if plt is None:
        import matplotlib.pyplot as plt
        import importlib
        mod = importlib.import_module("palettable.colorbrewer.%s" %
                                      color_cycle[0])
        colors = getattr(mod, color_cycle[1]).mpl_colors
        from cycler import cycler

        plt.figure(figsize=(width, height), facecolor="w", dpi=dpi)
        ax = plt.gca()
        ax.set_prop_cycle(cycler('color', colors))
    else:
        fig = plt.gcf()
        fig.set_size_inches(width, height)
    plt.xticks(fontsize=ticksize)
    plt.yticks(fontsize=ticksize)

    ax = plt.gca()
    ax.set_title(ax.get_title(), size=width * 4)

    labelsize = int(width * 3)

    ax.set_xlabel(ax.get_xlabel(), size=labelsize)
    ax.set_ylabel(ax.get_ylabel(), size=labelsize)

    return plt


def pretty_plot_two_axis(x, y1, y2, xlabel=None, y1label=None, y2label=None,
                         width=8, height=None, dpi=300, **plot_kwargs):
    """
    Variant of pretty_plot that does a dual axis plot. Adapted from matplotlib
    examples. Makes it easier to create plots with different axes.

    Args:
        x (np.ndarray/list): Data for x-axis.
        y1 (dict/np.ndarray/list): Data for y1 axis (left). If a dict, it will
            be interpreted as a {label: sequence}.
        y2 (dict/np.ndarray/list): Data for y2 axis (right). If a dict, it will
            be interpreted as a {label: sequence}.
        xlabel (str): If not None, this will be the label for the x-axis.
        y1label (str): If not None, this will be the label for the y1-axis.
        y2label (str): If not None, this will be the label for the y2-axis.
        width (float): Width of plot in inches. Defaults to 8in.
        height (float): Height of plot in inches. Defaults to width * golden
            ratio.
        dpi (int): Sets dot per inch for figure. Defaults to 300.
        plot_kwargs: Passthrough kwargs to matplotlib's plot method. E.g.,
            linewidth, etc.

    Returns:
        matplotlib.pyplot
    """

    import palettable.colorbrewer.diverging

    colors = palettable.colorbrewer.diverging.RdYlBu_4.mpl_colors
    c1 = colors[0]
    c2 = colors[-1]

    golden_ratio = (math.sqrt(5) - 1) / 2

    if not height:
        height = int(width * golden_ratio)

    import matplotlib.pyplot as plt
    width = 12
    labelsize = int(width * 3)
    ticksize = int(width * 2.5)
    styles = ["-", "--", "-.", "."]

    fig, ax1 = plt.subplots()
    fig.set_size_inches((width, height))
    if dpi:
        fig.set_dpi(dpi)
    if isinstance(y1, dict):
        for i, (k, v) in enumerate(y1.items()):
            ax1.plot(x, v, c=c1, marker='s', ls=styles[i % len(styles)],
                     label=k, **plot_kwargs)
        ax1.legend(fontsize=labelsize)
    else:
        ax1.plot(x, y1, c=c1, marker='s', ls='-', **plot_kwargs)

    if xlabel:
        ax1.set_xlabel(xlabel, fontsize=labelsize)

    if y1label:
        # Make the y-axis label, ticks and tick labels match the line color.
        ax1.set_ylabel(y1label, color=c1, fontsize=labelsize)

    ax1.tick_params('x', labelsize=ticksize)
    ax1.tick_params('y', colors=c1, labelsize=ticksize)

    ax2 = ax1.twinx()
    if isinstance(y2, dict):
        for i, (k, v) in enumerate(y2.items()):
            ax2.plot(x, v, c=c2, marker='o', ls=styles[i % len(styles)],
                     label=k)
        ax2.legend(fontsize=labelsize)
    else:
        ax2.plot(x, y2, c=c2, marker='o', ls='-')

    if y2label:
        # Make the y-axis label, ticks and tick labels match the line color.
        ax2.set_ylabel(y2label, color=c2, fontsize=labelsize)

    ax2.tick_params('y', colors=c2, labelsize=ticksize)
    return plt


def pretty_polyfit_plot(x, y, deg=1, xlabel=None, ylabel=None, **kwargs):
    r"""
    Convenience method to plot data with trend lines based on polynomial fit.

    Args:
        x: Sequence of x data.
        y: Sequence of y data.
        deg (int): Degree of polynomial. Defaults to 1.
        xlabel (str): Label for x-axis.
        ylabel (str): Label for y-axis.
        \\*\\*kwargs: Keyword args passed to pretty_plot.

    Returns:
        matplotlib.pyplot object.
    """
    plt = pretty_plot(**kwargs)
    pp = np.polyfit(x, y, deg)
    xp = np.linspace(min(x), max(x), 200)
    plt.plot(xp, np.polyval(pp, xp), 'k--', x, y, 'o')
    if xlabel:
        plt.xlabel(xlabel)
    if ylabel:
        plt.ylabel(ylabel)
    return plt


def periodic_table_heatmap(elemental_data, cbar_label="", cbar_label_size=14,
                           show_plot=False, cmap="YlOrRd", cmap_range=None, blank_color="grey",
                           value_format=None, max_row=9):
    """
    A static method that generates a heat map overlayed on a periodic table.

    Args:
         elemental_data (dict): A dictionary with the element as a key and a
            value assigned to it, e.g. surface energy and frequency, etc.
            Elements missing in the elemental_data will be grey by default
            in the final table elemental_data={"Fe": 4.2, "O": 5.0}.
         cbar_label (string): Label of the colorbar. Default is "".
         cbar_label_size (float): Font size for the colorbar label. Default is 14.
         cmap_range (tuple): Minimum and maximum value of the colormap scale.
            If None, the colormap will autotmatically scale to the range of the
            data.
         show_plot (bool): Whether to show the heatmap. Default is False.
         value_format (str): Formatting string to show values. If None, no value
            is shown. Example: "%.4f" shows float to four decimals.
         cmap (string): Color scheme of the heatmap. Default is 'YlOrRd'.
            Refer to the matplotlib documentation for other options.
         blank_color (string): Color assigned for the missing elements in
            elemental_data. Default is "grey".
         max_row (integer): Maximum number of rows of the periodic table to be
            shown. Default is 9, which means the periodic table heat map covers
            the first 9 rows of elements.
    """

    # Convert primitive_elemental data in the form of numpy array for plotting.
    if cmap_range is not None:
        max_val = cmap_range[1]
        min_val = cmap_range[0]
    else:
        max_val = max(elemental_data.values())
        min_val = min(elemental_data.values())

    max_row = min(max_row, 9)

    if max_row <= 0:
        raise ValueError("The input argument 'max_row' must be positive!")

    value_table = np.empty((max_row, 18)) * np.nan
    blank_value = min_val - 0.01

    for el in Element:
        if el.row > max_row:
            continue
        value = elemental_data.get(el.symbol, blank_value)
        value_table[el.row - 1, el.group - 1] = value

    # Initialize the plt object
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    plt.gcf().set_size_inches(12, 8)

    # We set nan type values to masked values (ie blank spaces)
    data_mask = np.ma.masked_invalid(value_table.tolist())
    heatmap = ax.pcolor(data_mask, cmap=cmap, edgecolors='w', linewidths=1,
                        vmin=min_val - 0.001, vmax=max_val + 0.001)
    cbar = fig.colorbar(heatmap)

    # Grey out missing elements in input data
    cbar.cmap.set_under(blank_color)

    # Set the colorbar label and tick marks
    cbar.set_label(cbar_label, rotation=270, labelpad=25, size=cbar_label_size)
    cbar.ax.tick_params(labelsize=cbar_label_size)

    # Refine and make the table look nice
    ax.axis('off')
    ax.invert_yaxis()

    # Label each block with corresponding element and value
    for i, row in enumerate(value_table):
        for j, el in enumerate(row):
            if not np.isnan(el):
                symbol = Element.from_row_and_group(i + 1, j + 1).symbol
                plt.text(j + 0.5, i + 0.25, symbol,
                         horizontalalignment='center',
                         verticalalignment='center', fontsize=14)
                if el != blank_value and value_format is not None:
                    plt.text(j + 0.5, i + 0.5, value_format % el,
                             horizontalalignment='center',
                             verticalalignment='center', fontsize=10)

    plt.tight_layout()

    if show_plot:
        plt.show()

    return plt


def format_formula(formula):
    """
    Converts str of chemical formula into
    latex format for labelling purposes

    Args:
        formula (str): Chemical formula
    """

    formatted_formula = ""
    number_format = ""
    for i, s in enumerate(formula):
        if s.isdigit():
            if not number_format:
                number_format = "_{"
            number_format += s
            if i == len(formula) - 1:
                number_format += "}"
                formatted_formula += number_format
        else:
            if number_format:
                number_format += "}"
                formatted_formula += number_format
                number_format = ""
            formatted_formula += s

    return r"$%s$" % (formatted_formula)


def van_arkel_triangle(list_of_materials, annotate=True):
    """
    A static method that generates a binary van Arkel-Ketelaar triangle to
        quantify the ionic, metallic and covalent character of a compound
        by plotting the electronegativity difference (y) vs average (x).
        See:
            A.E. van Arkel, Molecules and Crystals in Inorganic Chemistry,
                Interscience, New York (1956)
        and
            J.A.A Ketelaar, Chemical Constitution (2nd edn.), An Introduction
                to the Theory of the Chemical Bond, Elsevier, New York (1958)

    Args:
         list_of_materials (list): A list of computed entries of binary
            materials or a list of lists containing two elements (str).
         annotate (bool): Whether or not to lable the points on the
            triangle with reduced formula (if list of entries) or pair
            of elements (if list of list of str).
    """

    # F-Fr has the largest X difference. We set this
    # as our top corner of the triangle (most ionic)
    pt1 = np.array([(Element("F").X + Element("Fr").X) / 2,
                    abs(Element("F").X - Element("Fr").X)])
    # Cs-Fr has the lowest average X. We set this as our
    # bottom left corner of the triangle (most metallic)
    pt2 = np.array([(Element("Cs").X + Element("Fr").X) / 2,
                    abs(Element("Cs").X - Element("Fr").X)])
    # O-F has the highest average X. We set this as our
    # bottom right corner of the triangle (most covalent)
    pt3 = np.array([(Element("O").X + Element("F").X) / 2,
                    abs(Element("O").X - Element("F").X)])

    # get the parameters for the lines of the triangle
    d = np.array(pt1) - np.array(pt2)
    slope1 = d[1] / d[0]
    b1 = pt1[1] - slope1 * pt1[0]
    d = pt3 - pt1
    slope2 = d[1] / d[0]
    b2 = pt3[1] - slope2 * pt3[0]

    # Initialize the plt object
    import matplotlib.pyplot as plt

    # set labels and appropriate limits for plot
    plt.xlim(pt2[0] - 0.45, -b2 / slope2 + 0.45)
    plt.ylim(-0.45, pt1[1] + 0.45)
    plt.annotate("Ionic", xy=[pt1[0] - 0.3, pt1[1] + 0.05], fontsize=20)
    plt.annotate("Covalent", xy=[-b2 / slope2 - 0.65, -0.4], fontsize=20)
    plt.annotate("Metallic", xy=[pt2[0] - 0.4, -0.4], fontsize=20)
    plt.xlabel(r"$\frac{\chi_{A}+\chi_{B}}{2}$", fontsize=25)
    plt.ylabel(r"$|\chi_{A}-\chi_{B}|$", fontsize=25)

    # Set the lines of the triangle
    chi_list = [el.X for el in Element]
    plt.plot([min(chi_list), pt1[0]], [slope1 * min(chi_list) + b1, pt1[1]], 'k-', linewidth=3)
    plt.plot([pt1[0], -b2 / slope2], [pt1[1], 0], 'k-', linewidth=3)
    plt.plot([min(chi_list), -b2 / slope2], [0, 0], 'k-', linewidth=3)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)

    # Shade with appropriate colors corresponding to ionic, metallci and covalent
    ax = plt.gca()
    # ionic filling
    ax.fill_between([min(chi_list), pt1[0]],
                    [slope1 * min(chi_list) + b1, pt1[1]], facecolor=[1, 1, 0],
                    zorder=-5, edgecolor=[1, 1, 0])
    ax.fill_between([pt1[0], -b2 / slope2],
                    [pt1[1], slope2 * min(chi_list) - b1], facecolor=[1, 1, 0],
                    zorder=-5, edgecolor=[1, 1, 0])
    # metal filling
    XPt = Element("Pt").X
    ax.fill_between([min(chi_list), (XPt + min(chi_list)) / 2],
                    [0, slope1 * (XPt + min(chi_list)) / 2 + b1],
                    facecolor=[1, 0, 0], zorder=-3, alpha=0.8)
    ax.fill_between([(XPt + min(chi_list)) / 2, XPt],
                    [slope1 * ((XPt + min(chi_list)) / 2) + b1, 0],
                    facecolor=[1, 0, 0], zorder=-3, alpha=0.8)
    # covalent filling
    ax.fill_between([(XPt + min(chi_list)) / 2, ((XPt + min(chi_list)) / 2 + -b2 / slope2) / 2],
                    [0, slope2 * (((XPt + min(chi_list)) / 2 + -b2 / slope2) / 2) + b2],
                    facecolor=[0, 1, 0], zorder=-4, alpha=0.8)
    ax.fill_between([((XPt + min(chi_list)) / 2 + -b2 / slope2) / 2, -b2 / slope2],
                    [slope2 * (((XPt + min(chi_list)) / 2 + -b2 / slope2) / 2) + b2, 0],
                    facecolor=[0, 1, 0], zorder=-4, alpha=0.8)

    # Label the triangle with datapoints
    for entry in list_of_materials:
        if type(entry).__name__ not in ['ComputedEntry', 'ComputedStructureEntry']:
            X_pair = [Element(el).X for el in entry]
            formatted_formula = "%s-%s" % tuple(entry)
        else:
            X_pair = [Element(el).X for el in entry.composition.as_dict().keys()]
            formatted_formula = format_formula(entry.composition.reduced_formula)
        plt.scatter(np.mean(X_pair), abs(X_pair[0] - X_pair[1]), c='b', s=100)
        if annotate:
            plt.annotate(formatted_formula, fontsize=15,
                         xy=[np.mean(X_pair) + 0.005, abs(X_pair[0] - X_pair[1])])

    plt.tight_layout()
    return plt


def get_ax_fig_plt(ax=None, **kwargs):
    """
    Helper function used in plot functions supporting an optional Axes argument.
    If ax is None, we build the `matplotlib` figure and create the Axes else
    we return the current active figure.

    Args:
        kwargs: keyword arguments are passed to plt.figure if ax is not None.

    Returns:
        ax: :class:`Axes` object
        figure: matplotlib figure
        plt: matplotlib pyplot module.
    """
    import matplotlib.pyplot as plt
    if ax is None:
        fig = plt.figure(**kwargs)
        ax = fig.add_subplot(1, 1, 1)
    else:
        fig = plt.gcf()

    return ax, fig, plt


def get_ax3d_fig_plt(ax=None, **kwargs):
    """
    Helper function used in plot functions supporting an optional Axes3D
    argument. If ax is None, we build the `matplotlib` figure and create the
    Axes3D else we return the current active figure.

    Args:
        kwargs: keyword arguments are passed to plt.figure if ax is not None.

    Returns:
        ax: :class:`Axes` object
        figure: matplotlib figure
        plt: matplotlib pyplot module.
    """
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import axes3d
    if ax is None:
        fig = plt.figure(**kwargs)
        ax = axes3d.Axes3D(fig)
    else:
        fig = plt.gcf()

    return ax, fig, plt


def get_axarray_fig_plt(ax_array, nrows=1, ncols=1, sharex=False, sharey=False,
                        squeeze=True, subplot_kw=None, gridspec_kw=None,
                        **fig_kw):
    """
    Helper function used in plot functions that accept an optional array of Axes
    as argument. If ax_array is None, we build the `matplotlib` figure and
    create the array of Axes by calling plt.subplots else we return the
    current active figure.

    Returns:
        ax: Array of :class:`Axes` objects
        figure: matplotlib figure
        plt: matplotlib pyplot module.
    """
    import matplotlib.pyplot as plt

    if ax_array is None:
        fig, ax_array = plt.subplots(nrows=nrows, ncols=ncols, sharex=sharex,
                                     sharey=sharey, squeeze=squeeze,
                                     subplot_kw=subplot_kw,
                                     gridspec_kw=gridspec_kw, **fig_kw)
    else:
        fig = plt.gcf()
        ax_array = np.reshape(np.array(ax_array), (nrows, ncols))
        if squeeze:
            if ax_array.size == 1:
                ax_array = ax_array[0]
            elif any(s == 1 for s in ax_array.shape):
                ax_array = ax_array.ravel()

    return ax_array, fig, plt


def add_fig_kwargs(func):
    """
    Decorator that adds keyword arguments for functions returning matplotlib
    figures.

    The function should return either a matplotlib figure or None to signal
    some sort of error/unexpected event.
    See doc string below for the list of supported options.
    """
    from functools import wraps

    @wraps(func)
    def wrapper(*args, **kwargs):
        # pop the kwds used by the decorator.
        title = kwargs.pop("title", None)
        size_kwargs = kwargs.pop("size_kwargs", None)
        show = kwargs.pop("show", True)
        savefig = kwargs.pop("savefig", None)
        tight_layout = kwargs.pop("tight_layout", False)
        ax_grid = kwargs.pop("ax_grid", None)
        ax_annotate = kwargs.pop("ax_annotate", None)
        fig_close = kwargs.pop("fig_close", False)

        # Call func and return immediately if None is returned.
        fig = func(*args, **kwargs)
        if fig is None:
            return fig

        # Operate on matplotlib figure.
        if title is not None:
            fig.suptitle(title)

        if size_kwargs is not None:
            fig.set_size_inches(size_kwargs.pop("w"), size_kwargs.pop("h"),
                                **size_kwargs)

        if ax_grid is not None:
            for ax in fig.axes:
                ax.grid(bool(ax_grid))

        if ax_annotate:
            from string import ascii_letters
            tags = ascii_letters
            if len(fig.axes) > len(tags):
                tags = (1 + len(ascii_letters) // len(fig.axes)) * ascii_letters
            for ax, tag in zip(fig.axes, tags):
                ax.annotate("(%s)" % tag, xy=(0.05, 0.95), xycoords="axes fraction")

        if tight_layout:
            try:
                fig.tight_layout()
            except Exception as exc:
                # For some unknown reason, this problem shows up only on travis.
                # https://stackoverflow.com/questions/22708888/valueerror-when-using-matplotlib-tight-layout
                print("Ignoring Exception raised by fig.tight_layout\n", str(exc))

        if savefig:
            fig.savefig(savefig)

        import matplotlib.pyplot as plt
        if show:
            plt.show()
        if fig_close:
            plt.close(fig=fig)

        return fig

    # Add docstring to the decorated method.
    s = "\n\n" + """\
        Keyword arguments controlling the display of the figure:

        ================  ====================================================
        kwargs            Meaning
        ================  ====================================================
        title             Title of the plot (Default: None).
        show              True to show the figure (default: True).
        savefig           "abc.png" or "abc.eps" to save the figure to a file.
        size_kwargs       Dictionary with options passed to fig.set_size_inches
                          e.g. size_kwargs=dict(w=3, h=4)
        tight_layout      True to call fig.tight_layout (default: False)
        ax_grid           True (False) to add (remove) grid from all axes in fig.
                          Default: None i.e. fig is left unchanged.
        ax_annotate       Add labels to  subplots e.g. (a), (b).
                          Default: False
        fig_close         Close figure. Default: False.
        ================  ====================================================

"""

    if wrapper.__doc__ is not None:
        # Add s at the end of the docstring.
        wrapper.__doc__ += "\n" + s
    else:
        # Use s
        wrapper.__doc__ = s

    return wrapper
