#!/usr/bin/env python

"""
This module writes the input file for the original Fortran version of GIBBS (useful for testing and debugging).
"""

from __future__ import division
import warnings
import sys
import subprocess
import unittest
import pymatgen
import numpy as np
import os

__author__ = "Cormac Toher"
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Cormac Toher"
__email__ = "cormac.toher@duke.edu"
__date__ = "June 7, 2013"
    
def Gibbs_input_write(agl_data, amu2au):
    fo=open(agl_data.inpfname, "w")
    fo.write(agl_data.title)
    fo.write("\n")
    fo.write(agl_data.outfname)
    fo.write("\n")
    fo.write(str(agl_data.natoms))
    fo.write("\n")
    fo.write(str(agl_data.cellmass/amu2au))
    fo.write("\n")
    fo.write(str(agl_data.energ_inf))
    fo.write("\n")
    fo.write(str(agl_data.ieos))
    fo.write("\n")
    fo.write(str(agl_data.idebye))
    fo.write(" ")
    fo.write(str(agl_data.poisson))
    fo.write("\n")
    fo.write(str(agl_data.npressure))
    fo.write(" ")
    i = 0
    while i < agl_data.npressure:
        fo.write(str(agl_data.pressure[i]))
        fo.write(" ")   
        i = i + 1
    fo.write("\n")
    fo.write(str(agl_data.ntemperature))
    fo.write(" ")
    i = 0
    while i < agl_data.ntemperature:
        fo.write(str(agl_data.temperature[i]))
        fo.write(" ")   
        i = i + 1
    fo.write("\n")
    fo.write(str(agl_data.ninp))
    fo.write("\n")
    i = 0
    while i < agl_data.ninp:
        fo.write("\t")
        fo.write(str(agl_data.vol_inp[i]))
        fo.write("   ") 
        fo.write(str(agl_data.energ_inp[i]))
        fo.write("\n")   
        i = i + 1  
    fo.close()


if __name__ == "__main__":
    Gibbs_input_write()
