#!/usr/bin/env python

"""
Utilities for generating nicer plots.
"""

from __future__ import division

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
        width:
            Width of plot in inches. Defaults to 8in.
        height.
            Height of plot in inches. Defaults to width * golden ratio.
        plt:
            If plt is supplied, changes will be made to an existing plot.
            Otherwise, a new plot will be created.

    Returns:
        Matplotlib plot object with properly sized fonts.
    """
    golden_ratio = (math.sqrt(5) - 1.0) / 2.0
    if not height:
        height = int(width * golden_ratio)
    import matplotlib as mpl
    mpl.rcParams["axes.titlesize"] = width * 4
    mpl.rcParams["axes.labelsize"] = width * 3

    if plt is None:
        import matplotlib.pyplot as plt
        plt.figure(figsize=(width, height), facecolor="w")
    else:
        fig = plt.gcf()
        fig.set_size_inches(width, height)
    plt.xticks(fontsize=width * 2)
    plt.yticks(fontsize=width * 2)
    return plt
