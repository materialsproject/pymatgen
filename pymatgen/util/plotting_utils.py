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
import numpy as np

def gaussian_smear(x, y, mean_x, sigma):
    """
    Convolute densities with gaussian smearing function.
    
    Args:
        x:
            Sequence of x values (independent variable). Must contain mean_x
            within sequence.
        y:
            Sequence of y values (dependent variable). Must be same length as 
            x.
        mean_x:
            Mean of smearing Gaussian function.
        sigma:
            Standard deviation of smearing Gaussian function.
            
    Returns:
        Smeared y values.
    """
    from scipy.signal import convolve
    gauss = [1 / sigma / math.sqrt(2 * math.pi) * math.exp(-(val - mean_x) ** 2 / (2 * sigma ** 2)) for val in x]
    #need to normalize the gaussian pdf values to avoid scaling output
    gauss = np.array(gauss) / sum(gauss)
    smeared_y = convolve(gauss, y, 'full')
    #Locate index of the mean level, needed to perform shift
    shift = None
    for i in xrange(len(x)):
        if x[i] > mean_x:
            shift = i
            break
    if not shift:
        raise ValueError("mean_x not in range of x array.")
    return smeared_y[shift:shift + len(y)]


def get_publication_quality_plot(width=8, height=None):
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
