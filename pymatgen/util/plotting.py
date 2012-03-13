#!/usr/bin/env python

'''
Created on Mar 13, 2012
'''

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Mar 13, 2012"


import matplotlib.pyplot as plt
import math


def get_publication_quality_plot(width = 8, height = None):
    golden_ratio = (math.sqrt(5) - 1.0) / 2.0
    if not height:
        height = int(width * golden_ratio)
    plt.figure(figsize = (width, height), facecolor = 'w')
    plt.ylabel('', fontsize = 26, fontweight = 'bold')
    plt.xlabel('', fontsize = 26, fontweight = 'bold')
    plt.xticks(fontsize = 24)
    plt.yticks(fontsize = 24)
    plt.title('', fontsize = 26)

    return plt
