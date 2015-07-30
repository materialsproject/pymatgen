#!/usr/bin/env python

"""
This module finds the Debye temperature which produces the best fit of the Debye model to the calculated heat capacity.
"""

from __future__ import division
import warnings
import sys
import subprocess
import unittest
import pymatgen
from pymatgen.agl_thermal.agl_thermal import gauleg
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

#
# cvdebfit: determines which value of the Debye temperature produces best fit to heat capacity data
# 
# Runs through all integer values of Debye temperature in range between minimum and maximum values
# Generates heat capacity curve for each Debye temperature using Debye model (see eqn 23.26, Solid State Physics, Ashcroft and Mermin):
#
# c_V = 9 n k_B (T / Theta_D)^3 \int_0^{Theta_D / T} (x^4 e^x) / (e^x - 1)^2 dx
# 
def cvdebfit(Cv, tdmin, tdmax, nkb, Temp, npoints, agl_data):
    itmin = int(tdmin)
    itmax = int(tdmax)
    rmsmin = 1e30
    i = itmin
    while (i <= itmax):
        tdtrial = float(i)
        smsq = 0.0
        j = 1
        while (j < npoints):
            y = tdtrial / Temp[j]
            Deb = debint(y, agl_data)
            cvt = 9.0 * nkb * Deb;
            smsq = smsq + (cvt - Cv[j])**2
            j = j + 1
        rms = math.sqrt(smsq)
        if (rms < rmsmin):
            rmsmin = rms
            tdbest = tdtrial
        i = i + 1
    return tdbest
  


#
# debint: Evaluation of the Debye integral:   
#                         
#                                      |       x^4      |
#     Debye (y) = 3*y^(-3) * INT (0,y) | -------------- | dx
#                                      | (exp(x) - 1)^2 |
#
#    where y=ThetaD/T, being ThetaD Debye's temperature (K) and T the
#     absolute (thermodynamic) temperature. The integral is evaluated
#     using a Gauss-Legendre quadrature.
#
def fdebye(z):
    integz = ((z**4) * math.exp(z)) / ((math.exp(z) - 1)**2)
    return integz


def debint (y,agl_data):
    eps=1e-12
    cero=0.0
    maxnl=100
    pi = 3.14159265358979323846
    #
    #.....error condition controls
    #
    x = [0.0 for i in range(maxnl)]
    w = [0.0 for i in range(maxnl)]
    debye=3.0*pi*pi*pi*pi/y/y/y/15.0
    if (y <= 250):
        #
        #.....Loop with increasing number of Legendre points.
        #
        debye0=1e30
        nl = 5
        xabs=math.fabs(debye-debye0)
        while ((nl <= maxnl) and (xabs >= eps)):
            gauleg (cero, y, x, w, nl, agl_data)
            sum=0.0
            for i in xrange(nl):
                sum=sum + w[i] * (fdebye(x[i]))
            debye=sum/y/y/y
            xabs=math.fabs(debye-debye0)
            debye0=debye
            nl = nl + 5
    return debye



