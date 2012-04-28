#!/usr/bin/env python

'''
Utilities for generating nicer plots.
'''

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Mar 13, 2012"

import math

import matplotlib.pyplot as plt


def get_publication_quality_plot(width=8, height=None):
    """
    Provides a publication quality plot, with nice defaults for font sizes etc.
    
    Args:
        width:
            Width of plot in inches. Defaults to 8in.
        height.
            Height of plot in inches. Defaults to width * golden ratio.
    """
    golden_ratio = (math.sqrt(5) - 1.0) / 2.0
    if not height:
        height = int(width * golden_ratio)
    plt.figure(figsize=(width, height), facecolor='w')
    plt.ylabel('', fontsize=36)
    plt.xlabel('', fontsize=36)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.title('', fontsize=26)

    return plt
