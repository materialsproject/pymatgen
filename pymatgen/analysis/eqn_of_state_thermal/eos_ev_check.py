#!/usr/bin/env python

"""
This module sorts and checks the (E, V) data for use in the MP-AGL module.
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
__date__ = "April 8, 2014"



# ************************************************************************************************
#  This class contains a set of functions to sort and check the (E, V) data
# ************************************************************************************************

class eos_ev_check:
    def _init_(self):
        pass
        

    def checkminpos(self, vol_inp, energ_inp, ninp):
        minpos = 0
        vol_inp, energ_inp = self.evsort(vol_inp, energ_inp, ninp)
        imin = 0
        energmin = energ_inp[0]
        i = 1
        while i < ninp:
            if energ_inp[i] < energmin:
                energmin = energ_inp[i]
                imin = i
            i = i + 1
        if (imin == 0) or (imin == ninp - 1):
            minpos = 1
        return minpos


    def checkevconcav(self, vol_inp, energ_inp, ninp):
        concav = 0
        vol_inp, energ_inp = self.evsort(vol_inp, energ_inp, ninp)
        energmin = energ_inp[0]
        i = 0
        itmin = 0
        while i < ninp:
            if energ_inp[i] < energmin:
                energmin = energ_inp[i]
                itmin = i
            i = i + 1
        j = 0
        while ((energ_inp[j] - energ_inp[j-1])/(vol_inp[j] - vol_inp[j-1]) >= (energ_inp[j] - energ_inp[j+1])/(vol_inp[j] - vol_inp[j+1])) and (i < ninp - 2):
            j = j + 1
            if (j == ninp - 2) or (j >= itmin):
                concav = 1
                j = j + 1
            jtmax = 2
            i = j + 1
            while i < ninp:
                if (energ_inp[j] - energ_inp[j-1])/(vol_inp[j] - vol_inp[j-1]) < (energ_inp[j] - energ_inp[i])/(vol_inp[j] - vol_inp[i]):
                    j = j + 1
                    jtmax = i
                i = i + 1
            if jtmax <= itmin:
                concav = 1
        return concav



    def evsort(self, vol_inp, energ_inp, ninp):
        # First check if data is already in correct order
        i = 0
        icheck = 0
        while i < (ninp - 1):
            if (vol_inp[i+1] < vol_inp[i]):
                icheck = 1
            i = i + 1
        if (icheck == 0):
            #        eos_thermal_data.logstr = eos_thermal_data.logstr + "MP Eqn of State Thermal: (E, V) data already in correct order \n"
            return
        else:
            arrEV = np.zeros(ninp, dtype = {'names':['volume', 'energy'], 'formats':[float, float]})
        for i in xrange(ninp):
            arrEV[i] = (vol_inp[i], energ_inp[i])
        arrEVsorted = np.sort(arrEV, order='volume')
        for i in xrange(ninp):
            vol_inp[i] = arrEVsorted[i][0]
            energ_inp[i] = arrEVsorted[i][1]
        return vol_inp, energ_inp


    def evtsort(self, vol_inp, energ_inp, tdebye, ninp):
        # First check if data is already in correct order
        i = 0
        icheck = 0
        while i < (ninp - 1):
            if (vol_inp[i+1] < vol_inp[i]):
                icheck = 1
            i = i + 1
        if (icheck == 0):
        #        eos_thermal_data.logstr = eos_thermal_data.logstr + "MP Eqn of State Thermal: (E, V) data already in correct order \n"
            return
        else:
            arrEV = np.zeros(ninp, dtype = {'names':['volume', 'energy', 'debyetemp'], 'formats':[float, float, float]})
            for i in xrange(ninp):
                arrEV[i] = (vol_inp[i], energ_inp[i], tdebye[i])
            arrEVsorted = np.sort(arrEV, order='volume')
            for i in xrange(ninp):
                vol_inp[i] = arrEVsorted[i][0]
                energ_inp[i] = arrEVsorted[i][1]
                tdebye[i] = arrEVsorted[i][2]
        return vol_inp, energ_inp, tdebye

