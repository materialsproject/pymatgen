"""Utilities for generating nicer plots."""

from __future__ import annotations

import importlib
import math
from functools import wraps
from string import ascii_letters
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np
import palettable.colorbrewer.diverging
from matplotlib import cm, colors

from pymatgen.core import Element

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Literal

    from matplotlib.axes import Axes
    from matplotlib.figure import Figure
    from mpl_toolkits.mplot3d.axes3d import Axes3D
    from numpy.typing import NDArray


def pretty_plot(
    width: float = 8,
    height: float | None = None,
    ax: Axes = None,
    dpi: float | None = None,
    color_cycle: tuple[str, str] = ("qualitative", "Set1_9"),
) -> Axes:
    """Get a publication quality plot, with nice defaults for font sizes etc.

    Args:
        width (float): Width of plot in inches. Defaults to 8in.
        height (float): Height of plot in inches. Defaults to width * golden
            ratio.
        ax (Axes): If ax is supplied, changes will be made to an
            existing plot. Otherwise, a new plot will be created.
        dpi (float): Sets dot per inch for figure. Defaults to 300.
        color_cycle (tuple): Set the color cycle for new plots to one of the
            color sets in palettable. Defaults to a qualitative Set1_9.

    Returns:
        Axes: matplotlib axes object with properly sized fonts.
    """
    tick_size = int(width * 2.5)
    golden_ratio = (math.sqrt(5) - 1) / 2

    if not height:
        height = int(width * golden_ratio)

    if ax is None:
        mod = importlib.import_module(f"palettable.colorbrewer.{color_cycle[0]}")
        colors = getattr(mod, color_cycle[1]).mpl_colors
        from cycler import cycler

        plt.figure(figsize=(width, height), facecolor="w", dpi=dpi)
        ax = plt.gca()
        ax.set_prop_cycle(cycler("color", colors))
    else:
        fig = plt.gcf()
        fig.set_size_inches(width, height)

    plt.xticks(fontsize=tick_size)
    plt.yticks(fontsize=tick_size)

    ax.set_title(ax.get_title(), size=width * 4)

    label_size = int(width * 3)
    ax.set_xlabel(ax.get_xlabel(), size=label_size)
    ax.set_ylabel(ax.get_ylabel(), size=label_size)

    return ax


def pretty_plot_two_axis(
    x,
    y1,
    y2,
    xlabel=None,
    y1label=None,
    y2label=None,
    width: float = 8,
    height: float | None = None,
    dpi=300,
    **plot_kwargs,
):
    """Variant of pretty_plot that does a dual axis plot. Adapted from matplotlib
    examples. Makes it easier to create plots with different axes.

    Args:
        x (Sequence[float]): Data for x-axis.
        y1 (Sequence[float] | dict[str, Sequence[float]]): Data for y1 axis (left). If a dict, it will
            be interpreted as a {label: sequence}.
        y2 (Sequence[float] | dict[str, Sequence[float]]): Data for y2 axis (right). If a dict, it will
            be interpreted as a {label: sequence}.
        xlabel (str): If not None, this will be the label for the x-axis.
        y1label (str): If not None, this will be the label for the y1-axis.
        y2label (str): If not None, this will be the label for the y2-axis.
        width (float): Width of plot in inches. Defaults to 8in.
        height (float): Height of plot in inches. Defaults to width * golden
            ratio.
        dpi (int): Sets dot per inch for figure. Defaults to 300.
        plot_kwargs: Passthrough kwargs to matplotlib's plot method. e.g.
            linewidth, etc.

    Returns:
        plt.Axes: matplotlib axes object with properly sized fonts.
    """
    colors = palettable.colorbrewer.diverging.RdYlBu_4.mpl_colors
    c1 = colors[0]
    c2 = colors[-1]

    golden_ratio = (math.sqrt(5) - 1) / 2

    if not height:
        height = int(width * golden_ratio)

    width = 12
    label_size = int(width * 3)
    tick_size = int(width * 2.5)
    styles = ["-", "--", "-.", "."]

    fig, ax1 = plt.subplots()
    fig.set_size_inches((width, height))
    if dpi:
        fig.set_dpi(dpi)
    if isinstance(y1, dict):
        for idx, (key, val) in enumerate(y1.items()):
            ax1.plot(
                x,
                val,
                c=c1,
                marker="s",
                ls=styles[idx % len(styles)],
                label=key,
                **plot_kwargs,
            )
        ax1.legend(fontsize=label_size)
    else:
        ax1.plot(x, y1, c=c1, marker="s", ls="-", **plot_kwargs)

    if xlabel:
        ax1.set_xlabel(xlabel, fontsize=label_size)

    if y1label:
        # Make the y-axis label, ticks and tick labels match the line color.
        ax1.set_ylabel(y1label, color=c1, fontsize=label_size)

    ax1.tick_params("x", labelsize=tick_size)
    ax1.tick_params("y", colors=c1, labelsize=tick_size)

    ax2 = ax1.twinx()
    if isinstance(y2, dict):
        for idx, (key, val) in enumerate(y2.items()):
            ax2.plot(x, val, c=c2, marker="o", ls=styles[idx % len(styles)], label=key)
        ax2.legend(fontsize=label_size)
    else:
        ax2.plot(x, y2, c=c2, marker="o", ls="-")

    if y2label:
        # Make the y-axis label, ticks and tick labels match the line color.
        ax2.set_ylabel(y2label, color=c2, fontsize=label_size)

    ax2.tick_params("y", colors=c2, labelsize=tick_size)
    return ax1


def pretty_polyfit_plot(x: NDArray, y: NDArray, deg: int = 1, xlabel=None, ylabel=None, **kwargs):
    """Convenience method to plot data with trend lines based on polynomial fit.

    Args:
        x: Sequence of x data.
        y: Sequence of y data.
        deg (int): Degree of polynomial. Defaults to 1.
        xlabel (str): Label for x-axis.
        ylabel (str): Label for y-axis.
        kwargs: Keyword args passed to pretty_plot.

    Returns:
        plt.Axes
    """
    ax = pretty_plot(**kwargs)
    pp = np.polyfit(x, y, deg)
    xp = np.linspace(min(x), max(x), 200)
    ax.plot(xp, np.polyval(pp, xp), "k--", x, y, "o")
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    return ax


def _decide_fontcolor(rgba: tuple) -> Literal["black", "white"]:
    red, green, blue, _ = rgba
    if red * 0.299 + green * 0.587 + blue * 0.114 > (186 / 255):
        return "black"

    return "white"


def periodic_table_heatmap(
    elemental_data=None,
    cbar_label="",
    cbar_label_size=14,
    show_plot: bool = False,
    cmap="YlOrRd",
    cmap_range=None,
    blank_color="grey",
    edge_color="white",
    value_format=None,
    value_fontsize=10,
    symbol_fontsize=14,
    max_row: int = 9,
    readable_fontcolor=False,
    pymatviz: bool = True,
    **kwargs,
):
    """A static method that generates a heat map overlaid on a periodic table.

    Args:
        elemental_data (dict): A dictionary with the element as a key and a
            value assigned to it, e.g. surface energy and frequency, etc.
            Elements missing in the elemental_data will be grey by default
            in the final table elemental_data={"Fe": 4.2, "O": 5.0}.
        cbar_label (str): Label of the color bar. Default is "".
        cbar_label_size (float): Font size for the color bar label. Default is 14.
        cmap_range (tuple): Minimum and maximum value of the color map scale.
            If None, the color map will automatically scale to the range of the
            data.
        show_plot (bool): Whether to show the heatmap. Default is False.
        value_format (str): Formatting string to show values. If None, no value
            is shown. Example: "%.4f" shows float to four decimals.
        value_fontsize (float): Font size for values. Default is 10.
        symbol_fontsize (float): Font size for element symbols. Default is 14.
        cmap (str): Color scheme of the heatmap. Default is 'YlOrRd'.
            Refer to the matplotlib documentation for other options.
        blank_color (str): Color assigned for the missing elements in
            elemental_data. Default is "grey".
        edge_color (str): Color assigned for the edge of elements in the
            periodic table. Default is "white".
        max_row (int): Maximum number of rows of the periodic table to be
            shown. Default is 9, which means the periodic table heat map covers
            the standard 7 rows of the periodic table + 2 rows for the lanthanides
            and actinides. Use a value of max_row = 7 to exclude the lanthanides and
            actinides.
        readable_fontcolor (bool): Whether to use readable font color depending
            on background color. Default is False.
        pymatviz (bool): Whether to use pymatviz to generate the heatmap. Defaults to True.
            See https://github.com/janosh/pymatviz.
        kwargs: Passed to pymatviz.ptable_heatmap_plotly

    Returns:
        plt.Axes: matplotlib Axes object
    """
    if pymatviz:
        try:
            from pymatviz import ptable_heatmap_plotly

            if elemental_data:
                kwargs.setdefault("values", elemental_data)
                print('elemental_data is deprecated, use values={"Fe": 4.2, "O": 5.0} instead')
            if cbar_label:
                kwargs.setdefault("color_bar", {}).setdefault("title", cbar_label)
                print('cbar_label is deprecated, use color_bar={"title": cbar_label} instead')
            if cbar_label_size != 14:
                kwargs.setdefault("color_bar", {}).setdefault("title.font", {}).setdefault("size", cbar_label_size)
                print('cbar_label_size is deprecated, use color_bar={"title.font": {"size": cbar_label_size}} instead')
            if cmap:
                kwargs.setdefault("colorscale", cmap)
                print("cmap is deprecated, use colorscale=cmap instead")
            if cmap_range:
                kwargs.setdefault("cscale_range", cmap_range)
                print("cmap_range is deprecated, use cscale_range instead")
            if value_format:
                kwargs.setdefault("fmt", value_format)
                print("value_format is deprecated, use fmt instead")
            if blank_color != "grey":
                print("blank_color is deprecated")
            if edge_color != "white":
                print("edge_color is deprecated")
            if symbol_fontsize != 14:
                print("symbol_fontsize is deprecated, use font_size instead")
                kwargs.setdefault("font_size", symbol_fontsize)
            if value_fontsize != 10:
                print("value_fontsize is deprecated, use font_size instead")
                kwargs.setdefault("font_size", value_fontsize)
            if max_row != 9:
                print("max_row is deprecated, use max_row instead")
            if readable_fontcolor:
                print("readable_fontcolor is deprecated, use font_colors instead, e.g. ('black', 'white')")

            return ptable_heatmap_plotly(**kwargs)
        except ImportError:
            print(
                "You're using a deprecated version of periodic_table_heatmap(). Consider `pip install pymatviz` which "
                "offers an interactive plotly periodic table heatmap. You can keep calling this same function from "
                "pymatgen. Some of the arguments have changed which you'll be warned about. "
                "To disable this warning, pass pymatviz=False."
            )

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
        value = elemental_data.get(el.symbol, blank_value)
        if 57 <= el.Z <= 71:
            plot_row = 8
            plot_group = (el.Z - 54) % 32
        elif 89 <= el.Z <= 103:
            plot_row = 9
            plot_group = (el.Z - 54) % 32
        else:
            plot_row = el.row
            plot_group = el.group
        if plot_row > max_row:
            continue
        value_table[plot_row - 1, plot_group - 1] = value

    fig, ax = plt.subplots()
    plt.gcf().set_size_inches(12, 8)

    # We set nan type values to masked values (ie blank spaces)
    data_mask = np.ma.masked_invalid(value_table.tolist())
    heatmap = ax.pcolor(
        data_mask,
        cmap=cmap,
        edgecolors=edge_color,
        linewidths=1,
        vmin=min_val - 0.001,
        vmax=max_val + 0.001,
    )
    cbar = fig.colorbar(heatmap)

    # Grey out missing elements in input data
    cbar.cmap.set_under(blank_color)

    # Set the color bar label and tick marks
    cbar.set_label(cbar_label, rotation=270, labelpad=25, size=cbar_label_size)
    cbar.ax.tick_params(labelsize=cbar_label_size)

    # Refine and make the table look nice
    ax.axis("off")
    ax.invert_yaxis()

    # Set the scalarmap for fontcolor
    norm = colors.Normalize(vmin=min_val, vmax=max_val)
    scalar_cmap = cm.ScalarMappable(norm=norm, cmap=cmap)

    # Label each block with corresponding element and value
    for ii, row in enumerate(value_table):
        for jj, el in enumerate(row):  # type: ignore[arg-type]
            if not np.isnan(el):
                symbol = Element.from_row_and_group(ii + 1, jj + 1).symbol
                rgba = scalar_cmap.to_rgba(el)
                fontcolor = _decide_fontcolor(tuple(rgba)) if readable_fontcolor else "black"
                plt.text(
                    jj + 0.5,
                    ii + 0.25,
                    symbol,
                    horizontalalignment="center",
                    verticalalignment="center",
                    fontsize=symbol_fontsize,
                    color=fontcolor,
                )
                if el != blank_value and value_format is not None:
                    plt.text(
                        jj + 0.5,
                        ii + 0.5,
                        value_format % el,
                        horizontalalignment="center",
                        verticalalignment="center",
                        fontsize=value_fontsize,
                        color=fontcolor,
                    )

    plt.tight_layout()

    if show_plot:
        plt.show()

    return ax


def format_formula(formula: str) -> str:
    """Convert str of chemical formula into
    latex format for labelling purposes.

    Args:
        formula (str): Chemical formula
    """
    formatted_formula = ""
    number_format = ""
    for idx, char in enumerate(formula, start=1):
        if char.isdigit():
            if not number_format:
                number_format = "_{"
            number_format += char
            if idx == len(formula):
                number_format += "}"
                formatted_formula += number_format
        else:
            if number_format:
                number_format += "}"
                formatted_formula += number_format
                number_format = ""
            formatted_formula += char

    return f"${formatted_formula}$"


def van_arkel_triangle(list_of_materials: Sequence, annotate: bool = True):
    """A static method that generates a binary van Arkel-Ketelaar triangle to
    quantify the ionic, metallic and covalent character of a compound
    by plotting the electronegativity difference (y) vs average (x).
    See:
        A.E. van Arkel, Molecules and Crystals in Inorganic Chemistry,
            Interscience, New York (1956)
    and
        J.A.A Ketelaar, Chemical Constitution (2nd edition), An Introduction
            to the Theory of the Chemical Bond, Elsevier, New York (1958).

    Args:
        list_of_materials (list): A list of computed entries of binary
            materials or a list of lists containing two elements (str).
        annotate (bool): Whether or not to label the points on the
            triangle with reduced formula (if list of entries) or pair
            of elements (if list of list of str).

    Returns:
        plt.Axes: matplotlib Axes object
    """
    # F-Fr has the largest X difference. We set this
    # as our top corner of the triangle (most ionic)
    pt1 = np.array([(Element("F").X + Element("Fr").X) / 2, abs(Element("F").X - Element("Fr").X)])
    # Cs-Fr has the lowest average X. We set this as our
    # bottom left corner of the triangle (most metallic)
    pt2 = np.array(
        [
            (Element("Cs").X + Element("Fr").X) / 2,
            abs(Element("Cs").X - Element("Fr").X),
        ]
    )
    # O-F has the highest average X. We set this as our
    # bottom right corner of the triangle (most covalent)
    pt3 = np.array([(Element("O").X + Element("F").X) / 2, abs(Element("O").X - Element("F").X)])

    # get the parameters for the lines of the triangle
    d = np.array(pt1) - np.array(pt2)
    slope1 = d[1] / d[0]
    b1 = pt1[1] - slope1 * pt1[0]
    d = pt3 - pt1
    slope2 = d[1] / d[0]
    b2 = pt3[1] - slope2 * pt3[0]

    # set labels and appropriate limits for plot
    plt.xlim(pt2[0] - 0.45, -b2 / slope2 + 0.45)
    plt.ylim(-0.45, pt1[1] + 0.45)
    plt.annotate("Ionic", xy=(pt1[0] - 0.3, pt1[1] + 0.05), fontsize=20)
    plt.annotate("Covalent", xy=(-b2 / slope2 - 0.65, -0.4), fontsize=20)
    plt.annotate("Metallic", xy=(pt2[0] - 0.4, -0.4), fontsize=20)
    plt.xlabel(r"$\frac{\chi_{A}+\chi_{B}}{2}$", fontsize=25)
    plt.ylabel(r"$|\chi_{A}-\chi_{B}|$", fontsize=25)

    # Set the lines of the triangle
    chi_list = [el.X for el in Element]
    plt.plot(
        [min(chi_list), pt1[0]],
        [slope1 * min(chi_list) + b1, pt1[1]],
        "k-",
        linewidth=3,
    )
    plt.plot([pt1[0], -b2 / slope2], [pt1[1], 0], "k-", linewidth=3)
    plt.plot([min(chi_list), -b2 / slope2], [0, 0], "k-", linewidth=3)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)

    # Shade with appropriate colors corresponding to ionic, metallic and covalent
    ax = plt.gca()
    # ionic filling
    ax.fill_between(
        [min(chi_list), pt1[0]],
        [slope1 * min(chi_list) + b1, pt1[1]],
        facecolor=[1, 1, 0],
        zorder=-5,
        edgecolor=[1, 1, 0],
    )
    ax.fill_between(
        [pt1[0], -b2 / slope2],
        [pt1[1], slope2 * min(chi_list) - b1],
        facecolor=[1, 1, 0],
        zorder=-5,
        edgecolor=[1, 1, 0],
    )
    # metal filling
    x_pt = Element("Pt").X
    ax.fill_between(
        [min(chi_list), (x_pt + min(chi_list)) / 2],
        [0, slope1 * (x_pt + min(chi_list)) / 2 + b1],
        facecolor=[1, 0, 0],
        zorder=-3,
        alpha=0.8,
    )
    ax.fill_between(
        [(x_pt + min(chi_list)) / 2, x_pt],
        [slope1 * ((x_pt + min(chi_list)) / 2) + b1, 0],
        facecolor=[1, 0, 0],
        zorder=-3,
        alpha=0.8,
    )
    # covalent filling
    ax.fill_between(
        [(x_pt + min(chi_list)) / 2, ((x_pt + min(chi_list)) / 2 + -b2 / slope2) / 2],
        [0, slope2 * (((x_pt + min(chi_list)) / 2 + -b2 / slope2) / 2) + b2],
        facecolor=[0, 1, 0],
        zorder=-4,
        alpha=0.8,
    )
    ax.fill_between(
        [((x_pt + min(chi_list)) / 2 + -b2 / slope2) / 2, -b2 / slope2],
        [slope2 * (((x_pt + min(chi_list)) / 2 + -b2 / slope2) / 2) + b2, 0],
        facecolor=[0, 1, 0],
        zorder=-4,
        alpha=0.8,
    )

    # Label the triangle with data points
    for entry in list_of_materials:
        if type(entry).__name__ not in ["ComputedEntry", "ComputedStructureEntry"]:
            X_pair = [Element(el).X for el in entry]
            el_1, el_2 = entry
            formatted_formula = f"{el_1}-{el_2}"
        else:
            X_pair = [Element(el).X for el in entry.composition.as_dict()]
            formatted_formula = format_formula(entry.reduced_formula)
        plt.scatter(np.mean(X_pair), abs(X_pair[0] - X_pair[1]), c="b", s=100)
        if annotate:
            plt.annotate(
                formatted_formula,
                fontsize=15,
                xy=(np.mean(X_pair) + 0.005, abs(X_pair[0] - X_pair[1])),  # type: ignore[arg-type]
            )

    plt.tight_layout()
    return ax


def get_ax_fig(ax: Axes = None, **kwargs) -> tuple[Axes, Figure]:
    """Helper function used in plot functions supporting an optional Axes argument.
    If ax is None, we build the `matplotlib` figure and create the Axes else
    we return the current active figure.

    Args:
        ax (Axes, optional): Axes object. Defaults to None.
        kwargs: keyword arguments are passed to plt.figure if ax is not None.

    Returns:
        tuple[Axes, Figure]: matplotlib Axes object and Figure objects
    """
    if ax is None:
        fig = plt.figure(**kwargs)
        ax = fig.gca()
    else:
        fig = plt.gcf()

    return ax, fig


def get_ax3d_fig(ax: Axes = None, **kwargs) -> tuple[Axes3D, Figure]:
    """Helper function used in plot functions supporting an optional Axes3D
    argument. If ax is None, we build the `matplotlib` figure and create the
    Axes3D else we return the current active figure.

    Args:
        ax (Axes3D, optional): Axes3D object. Defaults to None.
        kwargs: keyword arguments are passed to plt.figure if ax is not None.

    Returns:
        tuple[Axes3D, Figure]: matplotlib Axes3D and corresponding figure objects
    """
    if ax is None:
        fig = plt.figure(**kwargs)
        ax = fig.add_subplot(projection="3d")
    else:
        fig = plt.gcf()

    return ax, fig


def get_axarray_fig_plt(
    ax_array,
    nrows=1,
    ncols=1,
    sharex: bool = False,
    sharey: bool = False,
    squeeze: bool = True,
    subplot_kw=None,
    gridspec_kw=None,
    **fig_kw,
):
    """Helper function used in plot functions that accept an optional array of Axes
    as argument. If ax_array is None, we build the `matplotlib` figure and
    create the array of Axes by calling plt.subplots else we return the
    current active figure.

    Returns:
        ax: Array of Axes objects
        figure: matplotlib figure
        plt: matplotlib pyplot module.
    """
    if ax_array is None:
        fig, ax_array = plt.subplots(
            nrows=nrows,
            ncols=ncols,
            sharex=sharex,
            sharey=sharey,
            squeeze=squeeze,
            subplot_kw=subplot_kw,
            gridspec_kw=gridspec_kw,
            **fig_kw,
        )
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
    """Decorator that adds keyword arguments for functions returning matplotlib
    figures.

    The function should return either a matplotlib figure or None to signal
    some sort of error/unexpected event.
    See doc string below for the list of supported options.
    """

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
            fig.set_size_inches(size_kwargs.pop("w"), size_kwargs.pop("h"), **size_kwargs)

        if ax_grid is not None:
            for ax in fig.axes:
                ax.grid(bool(ax_grid))

        if ax_annotate:
            tags = ascii_letters
            if len(fig.axes) > len(tags):
                tags = (1 + len(ascii_letters) // len(fig.axes)) * ascii_letters
            for ax, tag in zip(fig.axes, tags, strict=True):
                ax.annotate(f"({tag})", xy=(0.05, 0.95), xycoords="axes fraction")

        if tight_layout:
            try:
                fig.tight_layout()
            except Exception as exc:
                # For some unknown reason, this problem shows up only on travis.
                # https://stackoverflow.com/questions/22708888/valueerror-when-using-matplotlib-tight-layout
                print("Ignoring Exception raised by fig.tight_layout\n", str(exc))

        if savefig:
            fig.savefig(savefig)

        if show:
            plt.show()
        if fig_close:
            plt.close(fig=fig)

        return fig

    # Add docstring to the decorated method.
    doc_str = """\n\n
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
        ax_annotate       Add labels to subplots e.g. (a), (b).
                          Default: False
        fig_close         Close figure. Default: False.
        ================  ====================================================

"""

    if wrapper.__doc__ is not None:
        # Add s at the end of the docstring.
        wrapper.__doc__ += f"\n{doc_str}"
    else:
        # Use s
        wrapper.__doc__ = doc_str

    return wrapper
