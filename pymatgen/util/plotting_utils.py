# coding: utf-8

from __future__ import division, unicode_literals

"""
Utilities for generating nicer plots.
"""


__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Mar 13, 2012"

import math


def get_publication_quality_plot(width=8, height=None, plt=None):
    """
    Provides a publication quality plot, with nice defaults for font sizes etc.

    Args:
        width: Width of plot in inches. Defaults to 8in.
        height. Height of plot in inches. Defaults to width * golden ratio.
        plt: If plt is supplied, changes will be made to an existing plot.
            Otherwise, a new plot will be created.

    Returns:
        Matplotlib plot object with properly sized fonts.
    """
    ticksize = int(width * 2.5)

    golden_ratio = (math.sqrt(5) - 1.0) / 2.0

    if not height:
        height = int(width * golden_ratio)

    if plt is None:
        import matplotlib.pyplot as plt
        plt.figure(figsize=(width, height), facecolor="w")
    else:
        fig = plt.gcf()
        fig.set_size_inches(width, height)
    plt.xticks(fontsize=ticksize)
    plt.yticks(fontsize=ticksize)

    axes = plt.gca()
    axes.set_title(axes.get_title(), size=width * 4)

    labelsize = int(width * 3)

    axes.set_xlabel(axes.get_xlabel(), size=labelsize)
    axes.set_ylabel(axes.get_ylabel(), size=labelsize)

    return plt


def add_fig_kwargs(func):
    """
    Decorator that adds keyword arguments for functions returning matplotlib figure.
    See doc string below for the list of supported options.
    """ 
    from functools import wraps
    @wraps(func)
    def wrapper(*args, **kwargs):
        # pop the kwds used by the decorator.
        title = kwargs.pop("title", None)
        show = kwargs.pop("show", True)
        savefig = kwargs.pop("savefig", None)

        # Call func
        import matplotlib.pyplot as plt
        fig = func(*args, **kwargs)

        # Operate on matplotlib figure.
        if title is not None: fig.suptitle(title)
        if savefig: fig.savefig(savefig)
        if show: plt.show()
        return fig

    wrapper.__doc__ += "\n" + """\
    keyword arguments controlling the display of the figure:

    ================  ====================================================
    kwargs            Meaning
    ================  ====================================================
    title             Title of the plot (Default: None).
    show              True to show the figure (Default True).
    savefig           'abc.png' or 'abc.eps' to save the figure to a file.
    ================  ===================================================="""
    return wrapper
