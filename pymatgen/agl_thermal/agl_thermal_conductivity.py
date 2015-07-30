#!/usr/bin/env python

"""
This module calculates the thermal conductivity using the Leibfried-Schloemann equation.
"""

from __future__ import division
import warnings
import sys
import subprocess
import unittest
import pymatgen
import numpy as np
import os
from numpy import matrix
from numpy import linalg
import math

__author__ = "Cormac Toher"
__copyright__ = "Copyright 2014, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Cormac Toher"
__email__ = "cormac.toher@duke.edu"
__date__ = "April 17, 2014"
       
 

def thermalconductivity(thetaD, temp, ntemp, kappaT, gammaD, voleqD, avmass):
    kboltz = 1.3807e-23
    hbar = 1.05459e-34
    dpi = 3.141592653589793238462643
    L = 5.72e+7
    tol = 1e-12
    third = 1.0/3.0
    kth = kboltz * thetaD / hbar
    kths = kth * kth
    vDcr = voleqD**third
    kappaD = ((0.849 * 3 * (4**third)) / (20 * (dpi**3) * (1.0 - (0.514 / gammaD) + (0.228 / (gammaD * gammaD)) ))) * kths * kboltz * vDcr * avmass / (hbar * gammaD * gammaD)
    for i in xrange(ntemp):
        if (math.fabs(temp[i]) < tol):
            kappaT.append(0.0)
        else: 
            kappaT.append(kappaD * thetaD / temp[i])
    return kappaD
  

