#!/usr/bin/env python

"""
This module is a wrapper for the GIBBS module.
"""

from __future__ import division
import warnings
import sys
import subprocess
import unittest
import pymatgen
from pymatgen.io.vaspio import Poscar
from pymatgen.io.vaspio import Vasprun
from pymatgen.io.cifio import CifWriter
from pymatgen.io.cifio import CifParser
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.io.vaspio_set import MITVaspInputSet
from pymatgen.io.vaspio_set import AbstractVaspInputSet
from pymatgen.agl_thermal.agl_gibbs import gibbs
from pymatgen.agl_thermal.agl_vasp_output_read import Vasp_output_read
from pymatgen.agl_thermal.agl_ev_check import checkminpos
from pymatgen.agl_thermal.agl_ev_check import checkevconcav
from pymatgen.agl_thermal.agl_input_write import Gibbs_input_write
from pymatgen.agl_thermal.agl_cv_debye_fit import cvdebfit
from pymatgen.agl_thermal.agl_thermal_conductivity import thermalconductivity
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

def agl_run(inpfilename, nstructs=28, stepsize=0.01, ieos = 0, idebye = 0, poisson = 0.25, npres = 11, spres = 2.0, ntemp = 201, stemp = 10.0, energy_inf=0.0):
    agl_data.vol_inp = []
    agl_data.energ_inp = []
    agl_data.tdebye = []
    agl_data.pressure = []
    agl_data.temperature = []
    agl_data.iord = []
    agl_data.mfit = 500
    agl_data.mpar = 20
    agl_data.maxloops = 250
    agl_data.maiG = 2
    agl_data.opt_g = False
    agl_data.iopt_g = 0
    agl_data.fittype = 0
    agl_data.fterr = 0
    agl_data.third = 1.0 / 3.0
    agl_data.hy2kjmol = 2625.5 
    agl_data.hart2ev = 27.211383
    agl_data.ev2hartree = 0.036749309
    agl_data.au2gpa = 29421.4  
    agl_data.pckbau = 3.166830e-6 
    agl_data.pi = 3.14159265358979323846
    agl_data.kj2unit = 0.1202707236
    agl_data.energ_inf = energy_inf
    pcamu = 1.6605402e-24
    pcme = 9.1093897e-28
    amu2au = pcamu/pcme
    amu2kg = 1.660538921e-27
    mod_struct = []
#    initstruct = Structure()
    initstruct = pymatgen.read_structure(inpfilename)
    Vasp_output_read(agl_data, initstruct, nstructs, stepsize)
    for i in xrange(nstructs):
        print("Volume = ", agl_data.vol_inp[i], "Energy = ", agl_data.energ_inp[i])
    print(agl_data.energ_inf)
    agl_data.natoms = len(initstruct.sites)
    cellmass = 0.0
    for i in xrange(agl_data.natoms):
        cellmass = cellmass + initstruct.species[i].atomic_mass
    agl_data.title = ""
    for i in xrange(agl_data.natoms):
        agl_data.title = agl_data.title + initstruct.species[i].symbol
    agl_data.npressure = npres
    agl_data.ntemperature = ntemp
    for i in xrange(npres):
        idoub = i
        agl_data.pressure.append(idoub * spres)
    for i in xrange(ntemp):
        idoub = i
        agl_data.temperature.append(idoub * stemp)
    agl_data.savpres = True
    agl_data.outfname = agl_data.title + '.out'
    agl_data.inpfname = agl_data.title + '.inp'
    agl_data.cellmass = cellmass * amu2au
    agl_data.ieos = int(ieos)
    agl_data.idebye = int(idebye)
    agl_data.poisson = float(poisson)
    agl_data.ninp = int(nstructs)
    agl_data.ndata = agl_data.ninp
    i = 0
    while i < agl_data.ninp:
        agl_data.iord.append(int(i))
        i = i + 1
    minpos = checkminpos(agl_data)
    concav = checkevconcav(agl_data)
    if (minpos != 0):
        print("Error: minimum is not contained in calculated (E, V) data")
    if (concav != 0):
        print("Error: (E, V) data is not concave")
    # Write input file for original Fortran version of GIBBS algorithm (useful for debugging and testing)
    Gibbs_input_write(agl_data, amu2au)
    # Run GIBBS algorithm to calculate thermal properties
    gibbs(agl_data)
    print(agl_data.logstr)
    print(agl_data.outfname)
    fo=open(agl_data.outfname, "w")  
    fo.write(agl_data.outstr)
    fot=open("agl_thermal_properties.out", "w")
    fot.write(agl_data.outstr_thermal)
    fotd=open("Debye_temperature", "w")
    fotd.write("# Debye temperature as a function of temperature \n")
    for j in xrange(agl_data.ntemperature):
        jstr = str(agl_data.temperature[j]) + "\t" + str(agl_data.tdeb0p[j]) + "\n"
        fotd.write(jstr)
    tdmin = agl_data.tdeb0p[0]
    tdmax = agl_data.tdeb0p[0]
    j = 1
    while (j < agl_data.ntemperature):
        if (agl_data.tdeb0p[j] > tdmax):
            tdmax = agl_data.tdeb0p[j]
            jtdmax = j
        if (agl_data.tdeb0p[j] < tdmin):
            tdmin = agl_data.tdeb0p[j]
            jtdmin = j
        j = j + 1
    nkb = 1.0 * agl_data.natoms
    tdbestfit = cvdebfit(agl_data.cvu0p, tdmin, tdmax, nkb, agl_data.temperature, agl_data.ntemperature, agl_data)
    print("Best fit Debye temperature = ", tdbestfit)
    difftdmin = math.fabs(tdbestfit - agl_data.tdeb0p[0])
    jtdbest = 0
    j = 1
    while (j < agl_data.ntemperature):
        difftd = math.fabs(tdbestfit - agl_data.tdeb0p[j])
        if (difftd < difftdmin):
            jtdbest = j
            difftdmin = difftd
        j = j + 1
    jt300 = 0
    j = 1
    difft300min = math.fabs(agl_data.temperature[j] - 300.0) 
    while (j < agl_data.ntemperature):
        difft300 = math.fabs(agl_data.temperature[j] - 300.0)
        if (difft300 < difft300min):
            jt300 = j
            difft300min = difft300
        j = j + 1    
    kappaT =[]
    print("jtdbest = ", jtdbest)
    print("jt300 = ", jt300)
    ang32bohr3 = 6.74833303710415
    voleqD = (agl_data.volminsav[jtdbest][0] / ang32bohr3) * 1e-30
    avmasskg = (cellmass * amu2kg) / agl_data.natoms
    thetaa = tdbestfit * (agl_data.natoms**(-agl_data.third))
    kappaD = thermalconductivity(thetaa, agl_data.temperature, agl_data.ntemperature, kappaT, agl_data.ga0p[jtdbest], voleqD, avmasskg)
    print("Thermal conductivity at Debye temperature = ", kappaD)
    print("Thermal conductivity at ", agl_data.temperature[jt300], "K = ", kappaT[jt300])
    print("agl_run complete")



class agl_data:
    def _init_(self):
        self.vol_inp = []
        self.energ_inp = []
        self.tdebye = []
        self.pressure = []
        self.temperature = []
        self.xconfigvector = []
        self.datatofit = []
       
 





if __name__ == "__main__":
    import sys
    rungibbs()
