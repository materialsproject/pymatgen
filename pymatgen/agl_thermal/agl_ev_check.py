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
#  This set of functions sorts and checks the (E, V) data
# ************************************************************************************************


def checkminpos(agl_data):
    minpos = 0
    iord = []
    i = 0
    while i < agl_data.ninp:
        iord.append(i)
        i = i + 1
    ifirst = 0
    ilast = agl_data.ninp - 1
    evsort(agl_data)
    imin = 0
    energmin = agl_data.energ_inp[0]
    i = 1
    while i < agl_data.ninp:
        if agl_data.energ_inp[i] < energmin:
            energmin = agl_data.energ_inp[i]
            imin = i
        i = i + 1
    if (imin == 0) or (imin == agl_data.ninp - 1):
        minpos = 1
    return minpos


def checkevconcav(agl_data):
    concav = 0
    evsort(agl_data)
    energmin = agl_data.energ_inp[0]
    i = 0
    itmin = 0
    while i < agl_data.ninp:
        if agl_data.energ_inp[i] < energmin:
            energmin = agl_data.energ_inp[i]
            itmin = i
        i = i + 1
    j = 0
    while ((agl_data.energ_inp[j] - agl_data.energ_inp[j-1])/(agl_data.vol_inp[j] - agl_data.vol_inp[j-1]) >= (agl_data.energ_inp[j] - agl_data.energ_inp[j+1])/(agl_data.vol_inp[j] - agl_data.vol_inp[j+1])) and (i < agl_data.ninp - 2):
        j = j + 1
    if (j == agl_data.ninp - 2) or (j >= itmin):
        concav = 1
    j = j + 1
    jtmax = 2
    i = j + 1
    while i < agl_data.ninp:
        if (agl_data.energ_inp[j] - agl_data.energ_inp[j-1])/(agl_data.vol_inp[j] - agl_data.vol_inp[j-1]) < (agl_data.energ_inp[j] - agl_data.energ_inp[i])/(agl_data.vol_inp[j] - agl_data.vol_inp[i]):
            j = j + 1
            jtmax = i
        i = i + 1
    if jtmax <= itmin:
        concav = 1
    return concav



def evsort(agl_data):
# First check if data is already in correct order
    i = 0
    icheck = 0
    while i < (agl_data.ninp - 1):
        if (agl_data.vol_inp[i+1] < agl_data.vol_inp[i]):
            icheck = 1
        i = i + 1
    if (icheck == 0):
#        agl_data.logstr = agl_data.logstr + "MP AGL: (E, V) data already in correct order \n"
        return
    else:
        arrEV = np.zeros(agl_data.ninp, dtype = {'names':['volume', 'energy'], 'formats':[float, float]})
        for i in xrange(agl_data.ninp):
            arrEV[i] = (agl_data.vol_inp[i], agl_data.energ_inp[i])
        arrEVsorted = np.sort(arrEV, order='volume')
        for i in xrange(agl_data.ninp):
            agl_data.vol_inp[i] = arrEVsorted[i][0]
            agl_data.energ_inp[i] = arrEVsorted[i][1]
    return


def evtsort(agl_data):
# First check if data is already in correct order
    i = 0
    icheck = 0
    while i < (agl_data.ninp - 1):
        if (agl_data.vol_inp[i+1] < agl_data.vol_inp[i]):
            icheck = 1
        i = i + 1
    if (icheck == 0):
#        agl_data.logstr = agl_data.logstr + "MP AGL: (E, V) data already in correct order \n"
        return
    else:
        arrEV = np.zeros(agl_data.ninp, dtype = {'names':['volume', 'energy', 'debyetemp'], 'formats':[float, float, float]})
        for i in xrange(agl_data.ninp):
            arrEV[i] = (agl_data.vol_inp[i], agl_data.energ_inp[i], agl_data.tdebye[i])
        arrEVsorted = np.sort(arrEV, order='volume')
        for i in xrange(agl_data.ninp):
            agl_data.vol_inp[i] = arrEVsorted[i][0]
            agl_data.energ_inp[i] = arrEVsorted[i][1]
            agl_data.tdebye[i] = arrEVsorted[i][2]
    return

