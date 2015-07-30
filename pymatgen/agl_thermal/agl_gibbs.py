#!/usr/bin/env python

"""
This module is a Python implementation of the GIBBS algorithm.
"""

from __future__ import division
import warnings
import sys
import subprocess
import unittest
import pymatgen
from pymatgen.agl_thermal.agl_ev_check import evsort
from pymatgen.agl_thermal.agl_ev_check import evtsort
from pymatgen.agl_thermal.agl_polynomial import polynomial_fit
from pymatgen.agl_thermal.agl_polynomial import polminbrack
from pymatgen.agl_thermal.agl_polynomial import minbrack
from pymatgen.agl_thermal.agl_polynomial import polmin
from pymatgen.agl_thermal.agl_polynomial import polin0
from pymatgen.agl_thermal.agl_polynomial import polin1
from pymatgen.agl_thermal.agl_polynomial import polin2
from pymatgen.agl_thermal.agl_polynomial import polin3
from pymatgen.agl_thermal.agl_polynomial import polin4
from pymatgen.agl_thermal.agl_eqn_state import numer
from pymatgen.agl_thermal.agl_eqn_state import vinet
from pymatgen.agl_thermal.agl_eqn_state import birch
from pymatgen.agl_thermal.agl_eqn_state import bcnt
from pymatgen.agl_thermal.agl_thermal import scdebye
from pymatgen.agl_thermal.agl_thermal import debfitt
from pymatgen.agl_thermal.agl_thermal import thermal
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

#
#.....gibbs - Debye-Grunseisen model for thermal properties
#
# Adapted from original Fortran version written by M. A. Blanco et al.
# See Computer Physics Communications 158, 57-72 (2004) and Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details
#
# This automated version described in Phys. Rev. B 90, 174107 (2014)
#
def gibbs(agl_data):
    #
    #.....Compute the function f(poisson) in the Debye model.
    #(See Miguel A. Blanco, PhD Thesis)
    #
    fx=2*(1 + agl_data.poisson)/3.0/(1 - 2*agl_data.poisson)
    gx=(1 + agl_data.poisson)/3.0/(1 - agl_data.poisson)
    hx=2.0*math.sqrt(fx**3)+math.sqrt(gx**3)
    agl_data.poratio=math.exp(-math.log(hx/3)/3)
    # Open main output file and write appropriate header for GIBBS output data for all p, T values
    if (agl_data.ieos >= 0):
        agl_data.outstr = 'MP AGL - (P,T) thermodynamics of crystals from (E,V) data \n'
        agl_data.outstr = agl_data.outstr + '\n'
        agl_data.outstr = agl_data.outstr + agl_data.title + '\n'
        agl_data.outstr = agl_data.outstr + 'Number of data points: ' + str(agl_data.ndata) + ' \n'
        agl_data.outstr = agl_data.outstr + ' \n'
    else:
        agl_data.outstr = 'MP AGL - (P,T) thermodynamics of crystals from (E,V) data \n'
        agl_data.outstr = agl_data.outstr + ' \n'
        agl_data.outstr = agl_data.outstr + str(agl_data.title) + ' \n'
        agl_data.outstr = agl_data.outstr + ' \n' 
        agl_data.outstr = agl_data.outstr + '# P(GPa) \t  T(K) \t V(bohr^3)  \t  G(kJ/mol) \t Bt(GPa)' + '\n'
    agl_data.logstr = 'MP AGL: Log file: status information, errors and warnings \n'
    for i in xrange(agl_data.ninp):
        agl_data.energ_inp[i] = agl_data.energ_inp[i] - agl_data.energ_inf
    if (agl_data.idebye == 1):
        evtsort(agl_data)
    else:
        evsort(agl_data)
    #
    #.....Find global minimum of energy data points
    #
    etref = agl_data.energ_inp[0]
    itmin = 0
    for i in xrange(agl_data.ndata):
        if (agl_data.energ_inp[i] < etref):
            etref = agl_data.energ_inp[i]
            itmin = i
    #
    #.....Check that all points show concave patterns
    #.....First check where the concave pattern begins
    #
    j = 1
    while ( ( (agl_data.energ_inp[j]-agl_data.energ_inp[j-1])/(agl_data.vol_inp[j]-agl_data.vol_inp[j-1]) >= (agl_data.energ_inp[j]-agl_data.energ_inp[j+1])/(agl_data.vol_inp[j]-agl_data.vol_inp[j+1]) ) and j < agl_data.ndata-2):
        j = j + 1
    # If only the last three points are concave, then MP AGL will not be able to use them to fit a polynomial
    # MP AGL exits giving an error message
    if (j == agl_data.ndata-2):
        agl_data.logstr = agl_data.logstr + "MP AGL: All points show convex patterns \n"
        agl_data.grerr = 2;
        return  
    agl_data.energ_inp[0] = agl_data.energ_inp[j-1]
    agl_data.vol_inp[0] = agl_data.vol_inp[j-1]
    agl_data.energ_inp[1] = agl_data.energ_inp[j]
    agl_data.vol_inp[1] = agl_data.vol_inp[j]
    agl_data.energ_inp[2] = agl_data.energ_inp[j+1]
    agl_data.vol_inp[2] = agl_data.vol_inp[j+1]
    if (agl_data.idebye == 1):
        agl_data.tdebye[0] = agl_data.tdebye[j-1]
        agl_data.tdebye[1] = agl_data.tdebye[j]
        agl_data.tdebye[2] = agl_data.tdebye[j+1]
    j = j + 1
    jtmax = 2
    #
    #.....j marks the last accepted point, i the new trial point
    #
    i = j + 1
    while (i < (agl_data.ndata-1)):
        if (agl_data.fittype == 0 or agl_data.fittype == 3):
            if ( (agl_data.energ_inp[j]-agl_data.energ_inp[j-1])/(agl_data.vol_inp[j]-agl_data.vol_inp[j-1]) < (agl_data.energ_inp[j]-agl_data.energ_inp[i])/(agl_data.vol_inp[j]-agl_data.vol_inp[i]) ):
                j = j + 1
                jtmax = i
                agl_data.energ_inp[j] = agl_data.energ_inp[i]
                agl_data.vol_inp[j] = agl_data.vol_inp[i]
                if (agl_data.idebye == 1): 
                    agl_data.tdebye[j] = agl_data.tdebye[i]
        i = i + 1
    agl_data.ndata = j + 1
    #
    #.....search for the minimum of the accepted data
    #
    volref = agl_data.vol_inp[0]
    eref = agl_data.energ_inp[0]
    imin = 0
    i = 0
    while (i < agl_data.ndata):
        if (agl_data.energ_inp[i] < eref):
            volref = agl_data.vol_inp[i]
            eref = agl_data.energ_inp[i]
            imin = i
        i = i + 1
    agl_data.xconfigvector = []
    for i in xrange(agl_data.ndata):
        agl_data.xconfigvector.append((agl_data.vol_inp[i]/volref)**agl_data.third)
    # If the lowest energy corresponds to the smallest or largest volume, then the minimum lies outside of the accepted input data
    # MP AGL exits giving an error message
    if (imin == 0 or imin == (agl_data.ndata - 1)):
        agl_data.logstr = agl_data.logstr + "MP AGL: static minimum outside input data \n"
        agl_data.logstr = agl_data.logstr + "MP AGL: structure with lowest energy is " + str(imin) + " \n"
        agl_data.logstr = agl_data.logstr + "MP AGL: energy of this structure = " + str(eref) + " \n"
        agl_data.logstr = agl_data.logstr + "MP AGL: total number of structures = " + str(agl_data.ndata) + " \n"
        agl_data.grerr = 2
        return
    #
    #.....Obtain the polynomial fit of E(static) vs. x
    #     x = V^(1/3)
    #     epol contains coefficients of fitted polynomial
    #
    agl_data.datatofit = []
    for i in xrange(agl_data.ndata):
        agl_data.datatofit.append(agl_data.energ_inp[i])
    epol = []
    eerr = []
    nepol = polynomial_fit (imin, epol, eerr, agl_data)
    if (agl_data.fterr != 0):
        if (agl_data.fterr == 2):
            agl_data.logstr = agl_data.logstr + "MP AGL: Problem inverting matrix to fit polynomial \n"
            agl_data.grerr = 4
        else:
            agl_data.logstr = agl_data.logstr + "MP AGL: No polynomial fit found with a minimum within bounds of input data \n"
            agl_data.grerr = 2
        return
    if(nepol == 2):
        nepol = nepol + 1
        epol.append(0.0)
    #
    # Find minimum of polynomial epol to find minimum of (E, V) curve
    # First bracket minimum of (E, V) data
    itry = imin
    itry = minbrack (itry, agl_data)
    if (agl_data.ierr != 0):
        agl_data.logstr = agl_data.logstr + "MP AGL: Cannot find minimum of (E, V) data \n"
        agl_data.grerr = 2
        return
    else: 
        agl_data.logstr = agl_data.logstr + "MP AGL: Minimum of (E, V) data is at point imin = " + str(itry) + " \n"
    #
    # Find minimum of epol, using itry as initial guess
    agl_data.pmerr = 0
    xmin = polmin (agl_data.xconfigvector[itry], agl_data.xconfigvector[max(itry-2,0)], agl_data.xconfigvector[min(itry+2,agl_data.ndata-1)], nepol, epol, agl_data)
    # If polmin returns an error, then MP AGL tries shifting the trial point to try to correct the error
    # If this does not correct the error, then MP AGL exits giving an error message
    if (agl_data.pmerr != 0):
        if (agl_data.pmerr == 1): 
            while ((agl_data.pmerr == 1) and ((itry + 2) < agl_data.ndata)):
                itry = itry + 1
                agl_data.logstr = agl_data.logstr + "MP AGL: Resetting itry to itry = " + str(itry) + " \n"
                xmin = polmin(agl_data.xconfigvector[itry], agl_data.xconfigvector[max(itry-2,0)], agl_data.xconfigvector[min(itry+2,agl_data.ndata-1)], nepol, epol, agl_data)
            # If error indicator has changed from 1 to 2, then the bracket has shifted from one side of the minimum to the other
            # Need to expand the size of the bracket to incorporate the minimum
            if (agl_data.pmerr == 2):
                xmin = polmin (agl_data.xconfigvector[itry], agl_data.xconfigvector[max(itry-4,0)], agl_data.xconfigvector[min(itry+4,agl_data.ndata-1)], nepol, epol, xmin, agl_data)
            if (agl_data.pmerr == 0):
                agl_data.logstr = agl_data.logstr + "MP AGL: Minimum of (E, V) data is at point xmin = " + str(xmin) + " Bohr \n"
            else:
                agl_data.logstr = agl_data.logstr + "MP AGL: Cannot find minimum of (E, V) data \n"
                agl_data.grerr = 2
                return
        elif (agl_data.pmerr == 2):
            while ((agl_data.pmerr == 2) and ((itry - 2) >= 0)):
                itry = itry - 1
                agl_data.logstr = agl_data.logstr + "MP AGL: Resetting itry to itry = " + str(itry) + " \n"
                xmin = polmin (agl_data.xconfigvector[itry], agl_data.xconfigvector[max(itry-2,0)], agl_data[min(itry+2,agl_data.ndata-1)], nepol, epol, agl_data)
            # If error indicator has changed from 2 to 1, then the bracket has shifted from one side of the minimum to the other
            # Need to expand the size of the bracket to incorporate the minimum
            if (agl_data.pmerr == 1) :
                polmin (agl_data.xconfigvector[itry], agl_data.xconfigvector[max(itry-4,0)], agl_data.xconfigvector[min(itry+4,agl_data.ndata-1)], nepol, epol, agl_data)     
            if (pmerr == 0):
                agl_data.logstr = agl_data.logstr + "MP AGL: Minimum of (E, V) data is at point xmin = " + str(xmin) + " Bohr \n"
            else:
                agl_data.logstr = agl_data.logstr + "MP AGL: Cannot find minimum of (E, V) data \n"
                agl_data.grerr = 2
                return
        else: 
            agl_data.logstr = agl_data.logstr + "MP AGL: Cannot find minimum of (E, V) data \n"
            agl_data.grerr = 2
            return 
    Vmin = xmin**3
    # Evaluate polynomial at minimum point xmin to get Emin 
    Emin = polin0 (xmin, nepol, epol)
    # Write out Emin in different units
    agl_data.logstr = agl_data.logstr + "MP AGL: Minimum of (E, V) data is Emin = " + str(Emin) + " Hartree/cell \n"
    agl_data.logstr = agl_data.logstr + "MP AGL: Minimum of (E, V) data is Emin = " + str(Emin * agl_data.hy2kjmol) + " kJ/mol \n"
    agl_data.logstr = agl_data.logstr + "MP AGL: Minimum of (E, V) data is Emin = " + str(Emin * agl_data.hart2ev) + " eV/cell \n"
    #
    #.....Minimize G(static) - V; G = Gibbs free energy
    #
    gpol = []
    for i in xrange(nepol+1):
        gpol.append(epol[i]) 
    ngpol = nepol
    itry = imin
    agl_data.voleqmin = []
    g = []
    agl_data.bulkmod = []
    rerr = []
    for k in xrange(agl_data.npressure):
        #
        #.....bracket the minimum with the numerical function
        #
        for i in xrange(agl_data.ndata):
            agl_data.datatofit[i] = agl_data.energ_inp[i] + agl_data.vol_inp[i]*agl_data.pressure[k]/agl_data.au2gpa
        itry = minbrack (itry, agl_data)
        #
        #.....obtain the minimum of the fitted function
        #
        # Adds pV term onto coefficient of x^3; x^3 ~ V
        gpol[3] = epol[3] + agl_data.pressure[k]/agl_data.au2gpa * volref
        xmin = polmin (agl_data.xconfigvector[itry], agl_data.xconfigvector[max(itry-2,0)],agl_data.xconfigvector[min(itry+2,agl_data.ndata-1)], ngpol, gpol, agl_data)
        # Evaluates polynomial for (G, V) at minimum to get equilibrium values for V, G, and B
        agl_data.voleqmin.append((xmin**3.0) * volref)
        plnv = polin0 (xmin, ngpol, gpol)
        g.append(plnv * agl_data.hy2kjmol)
        plnv = polin2 (xmin, ngpol, gpol)
        agl_data.bulkmod.append(plnv * agl_data.au2gpa * xmin*xmin/(9.0*agl_data.voleqmin[k]))
        eact = polin0 (xmin, nepol, epol)
        twonepol = 2*nepol
        plnv = polin0 (xmin,twonepol,eerr)
        aerr = math.sqrt (math.fabs( plnv - eact*eact ))
        rerr.append(aerr / max (math.fabs(eact), math.fabs(eact)+aerr/2.0))
    #
    #.....Write static EOS results to stringstream
    #
    vol0pres=agl_data.voleqmin[0]
    gfe0pres=g[0]
    binp_bcnt = agl_data.bulkmod[0]
    if (agl_data.ieos >= 0):
        agl_data.outstr = agl_data.outstr + 'Static EOS calculation - Numerical results \n'
        agl_data.outstr = agl_data.outstr + 'Vmin(static;  P=0)    = ' + str(vol0pres) + ' bohr^3 \n'
        agl_data.outstr = agl_data.outstr + 'Gmin(static;  P=0)    = ' + str(gfe0pres) + ' kJ/mol \n'
        agl_data.outstr = agl_data.outstr + '\n'
        agl_data.outstr = agl_data.outstr + 'NUMERICAL EQUILIBRIUM PROPERTIES \n'
        agl_data.outstr = agl_data.outstr + '================================'
        agl_data.outstr = agl_data.outstr + '\n'
        agl_data.outstr = agl_data.outstr + '  P(GPa) \t G(kJ/mol) \t V(bohr^3)  \t   V/V0 \t B(GPa) \t rel.err. \n'
        agl_data.outstr = agl_data.outstr + ' ------------------------------------------------------------------------------------------- \n'
        for k in xrange(agl_data.npressure):
            agl_data.outstr = agl_data.outstr + '  ' + str(agl_data.pressure[k]).rjust(6) + '\t' + str(g[k]).rjust(10)[:10] + '\t' + str(agl_data.voleqmin[k]).rjust(10)[:10] + '\t' + str(agl_data.voleqmin[k]/vol0pres).rjust(7)[:7] + '\t     ' + str(agl_data.bulkmod[k]).rjust(10)[:10] + '\t     ' + str(round(rerr[k], 10)).rjust(12)[:12] + '\n'
    elif (agl_data.idebye == -1):
        for k in xrange(agl_data.npressure):
            agl_data.outstr = agl_data.outstr + '  ' + str(agl_data.pressure[k]).rjust(6) + '\t' + str(0.0).rjust(10)[:10] + '\t' + str(agl_data.voleqmin[k]).rjust(10)[:10] + '\t' + str(g[k]).rjust(10)[:10] + '\t     ' + str(agl_data.bulkmod[k]).rjust(10)[:10] + '\n' 
    else:
        for k in xrange(agl_data.npressure):
            agl_data.outstr = agl_data.outstr + '  ' + str(agl_data.pressure[k]).rjust(6) + '\t' + 'static'.rjust(10)[:10] + '\t' + str(agl_data.voleqmin[k]).rjust(10)[:10] + '\t' + str(g[k]).rjust(10)[:10] + '\t     ' + str(agl_data.bulkmod[k]).rjust(10)[:10] + '\n'
    #
    #.....Arrays to save EOS variables
    #
    agl_data.uder = []
    agl_data.ust = []
    F = [0.0 for i in range(agl_data.ndata)]
    agl_data.udyn = [0.0 for k in range(agl_data.npressure)]
    agl_data.pstatic = [0.0 for k in range(agl_data.npressure)]
    agl_data.gamma_G = [0.0 for k in range(agl_data.npressure)]
    agl_data.pfit = [0.0 for k in range(agl_data.npressure)]
    agl_data.alpha = [0.0 for k in range(agl_data.npressure)]
    theta = [0.0 for k in range(agl_data.npressure)]
    Cv = [0.0 for k in range(agl_data.npressure)]
    agl_data.astatic = []
    pst = []
    fpol = []
    ferr = []
    agl_data.gfepev = [[0.0 for k in range(agl_data.npressure)] for j in range(agl_data.ntemperature)]
    agl_data.xminsav = [[0.0 for k in range(agl_data.npressure)] for j in range(agl_data.ntemperature)]
    agl_data.volminsav = [[0.0 for k in range(agl_data.npressure)] for j in range(agl_data.ntemperature)]
    vol0p = []
    g0p = []
    b0p = []
    bp0p = []
    bpp0p = []
    agl_data.gfe0p = []
    agl_data.gfe0pev = []
    tpol = []
    bcnt_xk_sp = []
    bcnt_xp_sp = []
    bcnt_xv_sp = [] 
    bcnt_xg_sp = []
    agl_data.gamma_poly = []
    agl_data.cv0p = []
    agl_data.cvu0p = []
    agl_data.uvib0p = []
    agl_data.uvib0pmev = []
    agl_data.ent0p = []
    agl_data.ent0pmev = []
    agl_data.ent0pu = []
    agl_data.tdeb0p = []
    agl_data.ga0p = []
    agl_data.he0p = []
    agl_data.he0pmev = []
    agl_data.gvfe0p = []
    agl_data.gvfe0pev = []
    agl_data.cvp = [[0.0 for k in range(agl_data.npressure)] for j in range(agl_data.ntemperature)]
    agl_data.cvup = [[0.0 for k in range(agl_data.npressure)] for j in range(agl_data.ntemperature)]
    agl_data.uvibp = [[0.0 for k in range(agl_data.npressure)] for j in range(agl_data.ntemperature)]
    agl_data.uvibpmev = [[0.0 for k in range(agl_data.npressure)] for j in range(agl_data.ntemperature)]
    agl_data.entp = [[0.0 for k in range(agl_data.npressure)] for j in range(agl_data.ntemperature)]
    agl_data.entpmev = [[0.0 for k in range(agl_data.npressure)] for j in range(agl_data.ntemperature)]
    agl_data.entup = [[0.0 for k in range(agl_data.npressure)] for j in range(agl_data.ntemperature)]
    agl_data.tdebp = [[0.0 for k in range(agl_data.npressure)] for j in range(agl_data.ntemperature)]
    agl_data.gap = [[0.0 for k in range(agl_data.npressure)] for j in range(agl_data.ntemperature)]
    agl_data.hep = [[0.0 for k in range(agl_data.npressure)] for j in range(agl_data.ntemperature)]
    agl_data.hepmev = [[0.0 for k in range(agl_data.npressure)] for j in range(agl_data.ntemperature)]
    agl_data.a0p = []
    agl_data.bs0p = []
    agl_data.cp0p = []
    #
    #.....Debye temperature: numerical evaluation of thermal properties
    #
    if (agl_data.ieos == 0):
        numer (volref, nepol, epol, nepol, epol, True, agl_data)
        agl_data.outstr = agl_data.outstr + '\n'
        agl_data.outstr = agl_data.outstr + 'Debye temperature - numerical derivatives \n'
    #
    #.....Vinet equation of state.
    #
    elif (agl_data.ieos == 1):
        vinet (vol0pres, gfe0pres/agl_data.hy2kjmol, True, agl_data)
        agl_data.outstr = agl_data.outstr + '\n'
        agl_data.outstr = agl_data.outstr + 'Debye temperature - Vinet EOS derivatives \n'
    #
    #.....Birch-Murnaghan equation of state.
    #
    elif (agl_data.ieos == 2):
        iG = 2
        birch (vol0pres, gfe0pres/agl_data.hy2kjmol, iG, True, agl_data)
        agl_data.outstr = agl_data.outstr + '\n'
        agl_data.outstr = agl_data.outstr + 'Debye temperature - Birch-Murnaghan EOS derivatives \n'
    #
    #.....Vinet equation of state with numerical calculation to obtain all properties.
    #
    elif (agl_data.ieos == 3):
        vinet (vol0pres, gfe0pres/agl_data.hy2kjmol, True, agl_data)
        numer (volref, nepol, epol, nepol, epol, True, agl_data)
        agl_data.outstr = agl_data.outstr + '\n'
        agl_data.outstr = agl_data.outstr + 'Debye temperature - numerical derivatives \n'
    #
    #.....Birch-Murnaghan equation of state with numerical calculation to obtain all properties.
    #
    elif (agl_data.ieos == 4):
        iG = 2
        birch (vol0pres, gfe0pres/agl_data.hy2kjmol, iG, True, agl_data)
        numer (volref, nepol, epol, nepol, epol, True, agl_data)
        agl_data.outstr = agl_data.outstr + '\n'
        agl_data.outstr = agl_data.outstr + 'Debye temperature - numerical derivatives \n'
    #
    #.....Equation of state of Baonza-Caceres-Nunez-Taravillo.
    #
    elif (agl_data.ieos == 5):
        bcnt (vol0pres, gfe0pres/agl_data.hy2kjmol, binp_bcnt, True, agl_data)
        agl_data.outstr = agl_data.outstr + '\n'
        agl_data.outstr = agl_data.outstr + 'Debye temperature - BCNT EOS derivatives \n'
    #
    #.....Equation of state of Baonza-Caceres-Nunez-Taravillo, 
    #     but numerical calculation to obtain all properties.
    #
    elif (agl_data.ieos == 6):
        bcnt (vol0pres, gfe0pres/agl_data.hy2kjmol, binp_bcnt, True, agl_data)
        numer (volref, nepol, epol, nepol, epol, True, agl_data)
        agl_data.outstr = agl_data.outstr + '\n'
        agl_data.outstr = agl_data.outstr + 'Debye temperature - numerical EOS derivatives \n'
    #
    #.....Write Poisson coefficient and poisson ratio function
    #
    if (agl_data.ieos >= 0 and (agl_data.idebye == 0 or agl_data.idebye >= 2)):
        agl_data.outstr = agl_data.outstr + 'Poisson coefficient: ' + str(agl_data.poisson) + ', Poisson ratio function: ' + str(agl_data.poratio) + '\n'
        agl_data.outstr = agl_data.outstr + '\n'
    #
    #.....Calculate Debye temperatures at each volume
    #
    if (agl_data.bu1 < 0):
        agl_data.outstr = agl_data.outstr + 'gibbs: Warning! B''<0, will use idebye=3 \n'
        agl_data.idebye = 3
    if (agl_data.idebye == 1): 
        if (agl_data.ieos >= 0):
            agl_data.outstr = agl_data.outstr + '  V(bohr^3) \t   TDebye(K) \t   Computed(K) \n'
            agl_data.outstr = agl_data.outstr + '----------- \t ----------- \t ------------- \n'
        debfitt (ntpol, tpol, agl_data)
        if (agl_data.fterr != 0):
            if (agl_data.fterr == 2):
                agl_data.logstr = agl_data.logstr + "MP AGL: Problem inverting matrix to fit polynomial \n"
                agl_data.grerr = 4;
            else:
                agl_data.logstr = agl_data.logstr + "MP AGL: No polynomial fit found with a minimum within bounds of input data \n"
                agl_data.grerr = 2
            return
    elif (agl_data.idebye == 3):
        tdebyemin = ((6*agl_data.pi*agl_data.pi*agl_data.natoms*agl_data.vol_inp[imin]*agl_data.vol_inp[imin])**agl_data.third) / agl_data.pckbau * agl_data.poratio * math.sqrt(math.fabs(agl_data.uder[imin])/agl_data.cellmass)
    else:
        if (agl_data.ieos >= 0):
            agl_data.outstr = agl_data.outstr + "   V(bohr^3) \t  TDebye(K) \n"
            agl_data.outstr = agl_data.outstr + " ----------- \t ----------\n"
    ij = 0
    # Checks second derivative is positive at all points
    # If not, the function is convex at that point and that point is skipped
    for i in xrange(agl_data.ndata):
        if (agl_data.idebye != 1 and agl_data.uder[i] <= 0.0):
            if (i > itry):
                agl_data.ndata = i-1
                agl_data.logstr = agl_data.logstr + "MP AGL: Warning! convex function at i = " + str(i) + " \n"
                agl_data.logstr = agl_data.logstr + "MP AGL: Warning! uder = " + str(agl_data.uder[i]) + " \n"
                agl_data.logstr = agl_data.logstr + "MP AGL: The points following this one will be discarded \n"
            else:
                agl_data.logstr = agl_data.logstr + "MP AGL: Warning! convex function at i = " + str(i) + " \n"
                agl_data.logstr = agl_data.logstr + "MP AGL: Warning! uder = " + str(agl_data.uder[i]) + " \n"
                agl_data.logstr = agl_data.logstr + "MP AGL: This point will be skipped and the next concave one taken as the next point \n"
                agl_data.logstr = agl_data.logstr + "MP AGL: Recommend increasing the number of k-points and rerunning MP AGL \n"
        else: 
            # tmp is the Debye temperature for the structure with volume agl_data.vol_inp[ij]
            tmp = ((6*agl_data.pi*agl_data.pi*agl_data.natoms*agl_data.vol_inp[ij]*agl_data.vol_inp[ij])**agl_data.third) / agl_data.pckbau * agl_data.poratio * math.sqrt(math.fabs(agl_data.uder[ij])/agl_data.cellmass)
            if (agl_data.idebye == 3):
                agl_data.tdebye.append(tdebyemin)
                if (agl_data.ieos >= 0):
                    agl_data.outstr = agl_data.outstr + '  ' + str(agl_data.vol_inp[ij]).rjust(10)[:10] + '\t ' + str(agl_data.tdebye[ij]).rjust(10)[:10] + '\t' + str(tmp).rjust(10)[:10] + '\n'
            elif (agl_data.idebye == 1):
                if (agl_data.ieos >= 0):
                    agl_data.outstr = agl_data.outstr + '  ' + str(agl_data.vol_inp[ij]).rjust(10)[:10] + '\t  ' + str(agl_data.tdebye[ij]).rjust(10)[:10] + '\n'	
            else:
                agl_data.tdebye.append(tmp)
            if (agl_data.ieos >= 0):  
                agl_data.outstr = agl_data.outstr + '  ' + str(agl_data.vol_inp[ij]).rjust(10)[:10] + '\t ' + str(agl_data.tdebye[ij]).rjust(10)[:10] + '\n'
            ij = ij + 1  
    if (agl_data.ndata != ij):
        agl_data.ndata = ij
        agl_data.logstr = agl_data.logstr + "MP AGL: Number of (E, V) points set to = " + str(ij) + " \n"
    #
    #.....End of static calculation: exit if IDEBYE = -1
    #
    if (agl_data.idebye == -1):
        agl_data.logstr = agl_data.logstr + "MP AGL: End of static run ok! \n"
        return
    #
    #.....Initialize self-consistent Debye variables
    #     Debye temperature is calculated self-consistently if IDEBYE = 2
    #
    elif (agl_data.idebye == 2):
        scerr = scdebye (dzero, epol, eerr, nepol, imin, true, pst, agl_data)
        if(scerr):
            agl_data.logstr = agl_data.logstr + "MP AGL: Problem in self-consistent Debye! \n"
            scerr = False
            agl_data.grerr = 3
            return
    #
    #.....If iopt_g = 2 the variable opt_g is changed to false
    #    Optimization of beta is only performed for the static calculation
    #    Optimized value of beta from the static calculation used for the finite temperature calculations
    #
    if (agl_data.ieos == 5 or agl_data.ieos == 6):
        if (agl_data.iopt_g == 2):
            agl_data.opt_g = False
    # String to save thermal properties in a plottable format
    if (agl_data.idebye >= 0):
        agl_data.outstr_thermal = "#   T(K)    U(meV/cell)     F(meV/cell)      S(kB/cell)     Cv(kB/cell)      Theta_D(K)     Gruneisen parameter \n"
    #
    #.....Loop over temperatures
    #
    j = 0
    while j < agl_data.ntemperature:
        agl_data.logstr = agl_data.logstr + "MP AGL: Temperature = " + str(agl_data.temperature[j]) + "K \n"
        #
        #.....Self consistent Debye temperatures
        #
        if (agl_data.idebye == 2):
            scerr = scdebye(agl_data.temperature[j], fpol, ferr, nfpol, imin, false, pst, agl_data)
            # If scdebye returns an error for zero temperature, MP AGL exits giving an error
            # If scdebye returns an error for T > zero, MP AGL resets the maximum temperature to the previous value and skips the rest of the temperature loop
            # It then continues to complete the remainder of the MP AGL algorithm
            if(scerr):
                if (j == 0):      
                    agl_data.logstr = agl_data.logstr + "MP AGL: ffitt,  T = " + str(agl_data.temperature[j]) + ", itry = " + str(itry) + ", imin = " + str(imin) + ", agl_data.ndata = " + str(agl_data.ndata) + " \n"
                    agl_data.logstr = agl_data.logstr + "MP AGL: agl_data.datatofit = " + str(agl_data.datatofit) + " \n"
                    agl_data.logstr = agl_data.logstr + "MP AGL: Cannot find minimum for temperature = " + str(agl_data.temperature[j]) + " \n"
                    agl_data.grerr = 2
                    return
                else: 
                    agl_data.logstr = agl_data.logstr + "MP AGL: ffitt,  T = " + str(agl_data.temperature[j]) + ", itry = " + str(itry) + ", imin = " + str(imin) + ", agl_data.ndata = " + str(agl_data.ndata) + " \n"
                    agl_data.logstr = agl_data.logstr + "MP AGL: agl_data.datatofit = " + str(agl_data.datatofit) + " \n"
                    agl_data.logstr = agl_data.logstr + "MP AGL: Cannot find minimum of E + A at T = " + str(agl_data.temperature[j]) + " \n"
                    agl_data.logstr = agl_data.logstr + "MP AGL: Resetting maximum value of T to T = " + str(agl_data.temperature[j]) + " \n"
                    agl_data.logstr = agl_data.logstr + "MP AGL: Resetting number of temperature points to " + str(j) + " \n"
                    agl_data.ntemperature = j;
                    break
            ntpol = debfitt (tpol, agl_data)
            # If debfitt returns an error for zero temperature, MP AGL exits giving an error
            # If debfitt returns an error for T > zero, MP AGL resets the maximum temperature to the previous value and skips the rest of the temperature loop 
            # It then continues to complete the remainder of the MP AGL algorithm
            if (agl_data.fterr != 0):
                if (j == 0):     
                    agl_data.logstr = agl_data.logstr + "MP AGL: No polynomial fit found with a minimum within bounds of input data \n"
                    agl_data.logstr = agl_data.logstr + "MP AGL: ffitt,  T = " + str(agl_data.temperature[j]) + ", itry = " + str(itry) + ", imin = " + str(imin) + ", agl_data.ndata = " + str(agl_data.ndata)  + " \n"
                    agl_data.logstr = agl_data.logstr + "MP AGL: Cannot find minimum for temperature = " + str(agl_data.temperature[j]) + " \n"
                    if (agl_data.fterr == 2):
                        agl_data.logstr = agl_data.logstr + "MP AGL: Problem inverting matrix to fit polynomial \n"
                        agl_data.grerr = 4
                    else: 
                        agl_data.grerr = 2
                    return
                else: 
                    agl_data.logstr = agl_data.logstr + "MP AGL: No polynomial fit found with a minimum within bounds of input data \n"
                    agl_data.logstr = agl_data.logstr + "MP AGL: ffitt,  T = " + str(agl_data.temperature[j]) + ", itry = " + str(itry) + ", imin = " + str(imin) + ", agl_data.ndata = "  + str(agl_data.ndata) + " \n"
                    agl_data.logstr = agl_data.logstr + "MP AGL: Resetting maximum value of T to T = " + str(agl_data.temperature[j]) + " \n"
                    agl_data.logstr = agl_data.logstr + "MP AGL: Resetting number of temperature points to " + str(j) + " \n"
                    agl_data.ntemperature = j
                    break
        #
        #.....Obtain vibrational Helmholtz function fit (Debye model) for this value of the temperature, T
        #     At each volume, use the Debye temperature to get the vibrational Helmholtz energy (A = U - TS) 
        #     Add the value of A at each volume to the DFT final energy at that volume to agl_data.datatofit
        #     Fit (A + E, V) data to get the polynomial fpol (zero pressure, constant volume)
        #
        D = 0.0
        Derr = 0.0
        U = 0.0
        Uvib = 0.0
        Cvt = 0.0
        S = 0.0
        Fhelm = 0.0
        helm = 0.0
        ent = 0.0
        agl_data.thermal_helmholtz = 0.0
        agl_data.thermal_entropy = 0.0
        agl_data.thermal_energ = 0.0
        agl_data.thermal_cv = 0.0
        for i in xrange(agl_data.ndata):
            thermal (agl_data.tdebye[i], agl_data.temperature[j], D, Derr, agl_data)
            agl_data.datatofit[i] = agl_data.energ_inp[i] + agl_data.thermal_helmholtz
            F[i] = agl_data.thermal_helmholtz
        imin = minbrack (imin, agl_data)
        # If minbrack returns an error for zero temperature, MP AGL exits giving an error
        # If minbrack returns an error for T > zero, MP AGL resets the maximum temperature to the previous value and skips the rest of the temperature loop
        # It then continues to complete the remainder of the MP AGL algorithm
        if (agl_data.ierr != 0):
            if (j == 0):    
                agl_data.logstr = agl_data.logstr + "MP AGL: ffitt,  T = " + str(agl_data.temperature[j]) + ", itry = " + str(itry) + ", imin = " + str(imin) + ", agl_data.ndata = " + str(agl_data.ndata) + " \n"
                agl_data.logstr = agl_data.logstr + "MP AGL: agl_data.datatofit = " + str(agl_data.datatofit) + " \n"
                agl_data.logstr = agl_data.logstr + "MP AGL: Cannot find minimum for temperature = " + str(agl_data.temperature[j]) + " \n"
                agl_data.grerr = 2
                return 
            else: 
                agl_data.logstr = agl_data.logstr + "MP AGL: ffitt,  T = " + str(agl_data.temperature[j]) + ", itry = " + str(itry) + ", imin = " + str(imin) + ", agl_data.ndata = " + str(agl_data.ndata) + " \n"
                agl_data.logstr = agl_data.logstr + "MP AGL: agl_data.datatofit = " + str(agl_data.datatofit) + " \n"
                agl_data.logstr = agl_data.logstr + "MP AGL: Cannot find minimum of E + A at T = " + str(agl_data.temperature[j]) + " \n"
                agl_data.logstr = agl_data.logstr + "MP AGL: Resetting maximum value of T to T = " + str(agl_data.temperature[j]) + " \n"
                agl_data.logstr = agl_data.logstr + "MP AGL: Resetting number of temperature points to " + str(j) + " \n"
                agl_data.ntemperature = j
                break
        itry = imin
        # Fits polynomial to E + A
        nfpol = polynomial_fit (imin, fpol, ferr, agl_data)
        # If fitt returns an error for zero temperature, MP AGL exits giving an error
        # If fitt returns an error for T > zero, MP AGL resets the maximum temperature to the previous value and skips the rest of the temperature loop
        # It then continues to complete the remainder of the MP AGL algorithm
        if (agl_data.fterr != 0):
            if (j == 0):      
                agl_data.logstr = agl_data.logstr + "MP AGL: No polynomial fit found with a minimum within bounds of input data \n"
                agl_data.logstr = agl_data.logstr + "MP AGL: ffitt,  T = " + str(agl_data.temperature[j]) + ", itry = " + str(itry) + ", imin = " + str(imin) + ", agl_data.ndata = " + str(agl_data.ndata) + " \n"
                agl_data.logstr = agl_data.logstr + "MP AGL: Cannot find minimum for temperature = " + str(agl_data.temperature[j]) + " \n"
                if (agl_data.fterr == 2):
                    agl_data.logstr = agl_data.logstr + "MP AGL: Problem inverting matrix to fit polynomial \n"
                    agl_data.grerr = 4
                else: 
                    agl_data.grerr = 2
                return
            else: 
                agl_data.logstr = agl_data.logstr + "MP AGL: No polynomial fit found with a minimum within bounds of input data \n"
                agl_data.logstr = agl_data.logstr + "MP AGL: ffitt,  T = " + str(agl_data.temperature[j]) + ", itry = " + str(itry) + ", imin = " + str(imin) + ", agl_data.ndata = " + str(agl_data.ndata) + " \n"
                agl_data.logstr = agl_data.logstr + "MP AGL: Resetting maximum value of T to T = " + str(agl_data.temperature[j-1]) + " \n"
                agl_data.logstr = agl_data.logstr + "MP AGL: Resetting number of temperature points to " + str(j) + " \n"
                agl_data.ntemperature = j
                break
        if(nfpol == 2):
            nfpol = nfpol + 1
            fpol.append(0.0)
            ferr.append(0.0)
            ferr.append(0.0)
        #
        #.....Loop over pressures
        #
        for i in xrange(nfpol+1):
            gpol[i] = fpol[i]
        ngpol = nfpol
        k = 0
#        for k in xrange(agl_data.npressure):
        while k < agl_data.npressure:
            agl_data.logstr = agl_data.logstr + "MP AGL: Pressure = " + str(agl_data.pressure[k]) + " GPa \n"
            #
            #.....Calculate the value of the Gibbs free energy E + A + pV for each volume agl_data.vol_inp.at(i) and store in agl_data.datatofit
            #     Bracket the minimum of E + A + pV  with the numerical function
            #     Gibbs free energy = A + pV; constant pressure, variable volume
            #
            for i in xrange(agl_data.ndata):
                agl_data.datatofit[i] = agl_data.energ_inp[i] + F[i] + agl_data.vol_inp[i]*agl_data.pressure[k]/agl_data.au2gpa
            imin = minbrack (itry, agl_data)
            # If minbrack returns an error for zero pressure and zero temperature, MP AGL exits giving an error
            # If minbrack returns an error for p = 0, T > 0, MP AGL resets the maximum temperature to the previous value and skips the rest of the loop
            # If minbrack returns an error for p > zero, MP AGL resets the maximum pressure to the previous value and skips the rest of the pressure loop
            # It then continues to complete the remainder of the MP AGL algorithm
            if (agl_data.ierr != 0):
                if (k == 0):
                    if (j == 0):
                        agl_data.logstr = agl_data.logstr + "MP AGL: gmin, T = " + str(agl_data.temperature[j]) + ", P = " + str(agl_data.pressure[k]) + ", trial point = " + str(itry) + ", minimum point = " + str(imin) + ", total points = " + str(agl_data.ndata) + " \n"
                        agl_data.logstr = agl_data.logstr + "MP AGL: func = " + str(agl_data.datatofit[i]) + " \n"
                        agl_data.logstr = agl_data.logstr + "MP AGL: Cannot find minimum for pressure = " + str(agl_data.pressure[k]) + " and temperature = " + str(agl_data.temperature[j]) + " \n"
                        agl_data.grerr = 2
                        return
                    else:
                        agl_data.logstr = agl_data.logstr + "MP AGL: gmin, T = " + str(agl_data.temperature[j]) + ", P = " + str(agl_data.pressure[k]) + ", trial point = " + str(itry) + ", minimum point = " + str(imin) + ", total points = " + str(agl_data.ndata) + " \n"
                        agl_data.logstr = agl_data.logstr + "MP AGL: agl_data.datatofit = " + str(agl_data.datatofit[i]) + " \n"
                        agl_data.logstr = agl_data.logstr + "MP AGL: Cannot find minimum of E + pV + A at T = " + str(agl_data.temperature[j]) + ", p = " + str(agl_data.pressure[k]) + " \n"
                        agl_data.logstr = agl_data.logstr + "MP AGL: Resetting maximum value of T to T = " + str(agl_data.temperature[j-1]) + " \n"
                        agl_data.logstr = agl_data.logstr + "MP AGL: Resetting number of pressure points to " + str(j) + " \n"
                        agl_data.ntemperature = j
                        break	 
                else:
                    agl_data.logstr = agl_data.logstr + "MP AGL: gmin, T = " + str(agl_data.temperature[j]) + ", P = " + str(agl_data.pressure[k]) + ", trial point = " + str(itry) + ", minimum point = " + str(imin) + ", total points = " + str(agl_data.ndata) + " \n"
                    agl_data.logstr = agl_data.logstr + "MP AGL: agl_data.datatofit = " + str(agl_data.datatofit[i]) + " \n"
                    agl_data.logstr = agl_data.logstr + "MP AGL: Cannot find minimum of E + pV + A at T = " + str(agl_data.temperature[j]) + ", p = " + str(agl_data.pressure[k]) + " \n"
                    agl_data.logstr = agl_data.logstr + "MP AGL: Resetting maximum value of p to p = " + str(agl_data.pressure[k-1]) + " \n"
                    agl_data.logstr = agl_data.logstr + "MP AGL: Resetting number of pressure points to " + str(k) + " \n"
                    agl_data.npressure = k
                    break	 
            #
            #.....Find the volume which minimizes the fitted G function at (p, T) (polynomial fit is in gpol)
            #     G function is Gibbs free energy: G = E + pV + A; E = DFT energy
            #     For a given temperature and pressure, the equilibrium system is the one which minimizes the Gibbs free energy
            #
            gpol[3] = fpol[3] + agl_data.pressure[k]/agl_data.au2gpa * volref
            xmin = polmin (agl_data.xconfigvector[itry], agl_data.xconfigvector[max(itry-2,0)], agl_data.xconfigvector[min(itry+2,agl_data.ndata-1)], ngpol, gpol, agl_data)
            # If polmin returns an error, then MP AGL tries shifting the trial point to try to correct the error
            if (agl_data.pmerr != 0):
                agl_data.logstr = agl_data.logstr + "MP AGL: itry = " + str(itry) + " \n"
                agl_data.logstr = agl_data.logstr + "MP AGL: pmerr = " + str(agl_data.pmerr) + " \n"
                if (agl_data.pmerr == 1):
                    while ((agl_data.pmerr == 1) and ((itry + 2) < agl_data.ndata)):
                        itry = itry + 1
                        agl_data.logstr = agl_data.logstr + "MP AGL: Resetting itry to itry = " + str(itry) + " \n"
                        xmin = polmin (agl_data.xconfigvector[itry], agl_data.xconfigvector[max(itry-2,0)],agl_data.xconfigvector[min(itry+2,agl_data.ndata-1)], ngpol, gpol, agl_data)
                        agl_data.logstr = agl_data.logstr + "MP AGL: pmerr = " + str(agl_data.pmerr) + " \n"
                    # If error indicator has changed from 1 to 2, then the bracket has shifted from one side of the minimum to the other
                    # Need to expand the size of the bracket to incorporate the minimum
                    if (agl_data.pmerr == 2):
                        xmin = polmin (agl_data.xconfigvector[itry], agl_data.xconfigvector[max(itry-4,0)],agl_data.xconfigvector[min(itry+4,agl_data.ndata-1)], ngpol, gpol, agl_data)
                    # If polynomial minimum has still not been found successfully, writes polynomial and its first derivative to a file to aid debugging 
                    if (agl_data.pmerr != 0): 
                        agl_data.logstr = agl_data.logstr + "MP AGL: List of polynomial values \n"
                        for i in xrange(agl_data.ndata):
                            plnv = polin0(agl_data.xconfigvector[i], ngpol, gpol)
                            agl_data.logstr = agl_data.logstr + str(agl_data.xconfigvector[i]) + '\t' + str(plnv) + '\n'
                        agl_data.logstr = agl_data.logstr + "MP AGL: List of first derivative values of polynomial \n"
                        for i in xrange(agl_data.ndata):
                            plnv = polin1(agl_data.xconfigvector[i], ngpol, gpol)
                            agl_data.logstr = agl_data.logstr + str(agl_data.xconfigvector[i]) + '\t' + str(plnv) + '\n'
                elif (agl_data.pmerr == 2):
                    while ((agl_data.pmerr == 2) and ((itry - 2) >= 0)):
                        itry = itry - 1
                        agl_data.logstr = agl_data.logstr + "MP AGL: Resetting itry to itry = " + str(itry) + " \n"
                        xmin = polmin (agl_data.xconfigvector[itry], agl_data.xconfigvector[max(itry-2,0)],agl_data.xconfigvector[min(itry+2,agl_data.ndata-1)], ngpol, gpol, agl_data)
                        agl_data.logstr = agl_data.logstr + "MP AGL: pmerr = " + str(agl_data.pmerr) + " \n"
                    # If error indicator has changed from 2 to 1, then the bracket has shifted from one side of the minimum to the other
                    # Need to expand the size of the bracket to incorporate the minimum
                    if (agl_data.pmerr == 1):
                        xmin = polmin (agl_data.xconfigvector[itry], agl_data.xconfigvector[max(itry-4,0)],agl_data.xconfigvector[min(itry+4,agl_data.ndata-1)], ngpol, gpol, agl_data)
                    # If polynomial minimum has still not been found successfully, writes polynomial and its first derivative to a file to aid debugging 
                    if (agl_data.pmerr != 0): 
                        agl_data.logstr = agl_data.logstr + "MP AGL: List of polynomial values \n"
                        for i in xrange(agl_data.ndata):
                            plnv = polin0(agl_data.xconfigvector[i],ngpol,gpol)
                            agl_data.logstr = agl_data.logstr + str(agl_data.xconfigvector[i]) + '\t' + str(plnv) + '\n'
                        agl_data.logstr = agl_data.logstr + "MP AGL: List of first derivative values of polynomial \n"
                        for i in xrange(agl_data.ndata):
                            plnv = polin1(agl_data.xconfigvector[i], ngpol, gpol)
                            agl_data.logstr = agl_data.logstr + str(agl_data.xconfigvector[i]) + '\t' + str(plnv) + '\n'
            # Check that polmin has actually found a minimum and not a maximum
            plnv = polin2 (xmin, ngpol, gpol)
            if (plnv < 0.0):
                agl_data.logstr = agl_data.logstr + "MP AGL: Minimum polynomial fitted to E + pV is actually a maximum \n"
                agl_data.logstr = agl_data.logstr + "MP AGL: Rebracketing minimum of polynomial points \n"
                agl_data.ierr, itry = polminbrack (itry, ngpol, agl_data.ndata, agl_data.xconfigvector, gpol)
                if (agl_data.ierr != 0):
                    if (k == 0):
                        if (j == 0):
                            agl_data.logstr = agl_data.logstr + "MP AGL: Cannot find minimum of E + A + pV for p = " + str(agl_data.pressure[k]) + ", T = " + str(agl_data.temperature[j]) + " \n"
                            agl_data.grerr = 2
                            return
                        else :
                            agl_data.logstr = agl_data.logstr + "MP AGL: Cannot find minimum of E + A + pV for p = " + str(agl_data.pressure[k]) + ", T = " + str( agl_data.temperature[j]) + " \n"
                            agl_data.logstr = agl_data.logstr + "MP AGL: Resetting maximum value of T to T = " + str(agl_data.temperature[j-1]) + " \n"
                            agl_data.logstr = agl_data.logstr + "MP AGL: Resetting maximum number of temperature points to " + str(j) + " \n"
                            agl_data.ntemperature = j
                            break
                    else:
                        agl_data.logstr = agl_data.logstr + "MP AGL: Cannot find minimum of E + A + pV for p = " + str(agl_data.pressure[k]) + ", T = " + (agl_data.temperature[j]) + " \n"
                        agl_data.logstr = agl_data.logstr + "MP AGL: Resetting maximum value of p to p = " + str(agl_data.pressure[k-1]) + " \n"
                        agl_data.logstr = agl_data.logstr + "MP AGL: Resetting maximum number of pressure points to " + str(k) + " \n"
                        agl_data.npressure = k
                        break
                else: 
                    agl_data.logstr = agl_data.logstr + "MP AGL: Minimum of (E, V) data is at point imin = " + str(itry) + " \n"
                    xmin = polmin (agl_data.xconfigvector[itry], agl_data.xconfigvector[max(itry-2,0)],agl_data.xconfigvector[min(itry+2,agl_data.ndata-1)], ngpol, gpol, agl_data)
            # If polmin still returns an error for zero temperature and pressure after shifting the trial point, MP AGL exits giving an error
            # If polmin returns an error for T > zero, MP AGL resets the maximum temperature or pressure to the previous value and skips the rest of that loop
            # It then continues to complete the remainder of the MP AGL algorithm
            if (agl_data.pmerr != 0):
                if (k == 0):
                    if (j == 0):
                       agl_data.logstr = agl_data.logstr + "MP AGL: Cannot find minimum of polynomial fit for E + pV + A for p = " + str(agl_data.pressure[k]) + " and T = " + str(agl_data.temperature[j]) + " \n"
                       agl_data.grerr = 2
                       return
                    else: 
                        agl_data.logstr = agl_data.logstr + "MP AGL: Cannot find minimum of polynomial fit for E + pV + A at T = " + str(agl_data.temperature[j]) + ", p = " + str(agl_data.pressure[k]) + " \n"
                        agl_data.logstr = agl_data.logstr + "MP AGL: Resetting maximum temperature value to T = " + str(agl_data.temperature[j-1]) + " \n"
                        agl_data.logstr = agl_data.logstr + "MP AGL: Resetting number of temperature points to " + str(j) + " \n"
                        agl_data.ntemperature = j   
                        break
                else:
                    agl_data.logstr = agl_data.logstr + "MP AGL: Cannot find minimum of polynomial fit for E + pV + A at T = " + str(agl_data.temperature[j]) + ", p = " + str(agl_data.pressure[k]) + " \n"
                    agl_data.logstr = agl_data.logstr + "MP AGL: Resetting maximum pressure value to p = " + str(agl_data.pressure[k-1]) + " \n"
                    agl_data.logstr = agl_data.logstr + "MP AGL: Resetting number of pressure points to " + str(k) + " \n"
                    agl_data.npressure = k
                    break
            # Evaluates polynomial at minimum to get equilibrium V, G, and B
            agl_data.voleqmin[k] = (xmin**3) * volref
            plnv = polin0(xmin,ngpol,gpol)
            g[k] = plnv * agl_data.hy2kjmol
            agl_data.gfepev[j][k] = plnv*agl_data.hart2ev
            plnv = polin2(xmin,ngpol,gpol)
            agl_data.bulkmod[k] = plnv * agl_data.au2gpa * xmin*xmin/(9.0*agl_data.voleqmin[k])
            eact = polin0 (xmin,nfpol,fpol)
            twonfpol = 2 * nfpol
            plnv = polin0 (xmin, twonfpol, ferr)
            aerr = math.sqrt (math.fabs( plnv - eact*eact ))
            rerr[k] = aerr / max (math.fabs(eact), math.fabs(eact)+aerr/2.0)
            # Writes information about equilibrium volume to string
            agl_data.xminsav[j][k] = xmin
            agl_data.volminsav[j][k] = agl_data.voleqmin[k]
            k = k + 1
        #
        #.....ieos < 0 means minimum output
        #
        if (agl_data.ieos < 0):
            for k in xrange (agl_data.npressure):
                agl_data.outstr = agl_data.outstr + str(agl_data.pressure[k]) + "\t" + str(agl_data.temperature[j]) + "\t" + str(agl_data.voleqmin[k]) + "\t" + str(g[k]) + "\t" + str(agl_data.bulkmod[k]) + "\n"
        else: 
            #
            #.....Numerical energetic, geometric and elastic properties
            #
            vol0pres=agl_data.voleqmin[0]
            gfe0pres=g[0]
            binp_bcnt = agl_data.bulkmod[0]
            agl_data.outstr = agl_data.outstr + "\n"
            agl_data.outstr = agl_data.outstr + "Temperature:  T = " + str(agl_data.temperature[j]) + "K \n" 
            agl_data.outstr = agl_data.outstr + "Vmin(T; P=0) = " + str(vol0pres) + "bohr^3 \n" 
            agl_data.outstr = agl_data.outstr + "Gmin(T; P=0) = " + str(gfe0pres) + "kJ/mol \n"
            agl_data.outstr = agl_data.outstr + "\n"
            agl_data.outstr = agl_data.outstr + "NUMERICAL EQUILIBRIUM PROPERTIES \n" 
            agl_data.outstr = agl_data.outstr + "================================ \n"
            agl_data.outstr = agl_data.outstr + "\n"
            agl_data.outstr = agl_data.outstr + "  P(GPa) \t G(kJ/mol) \t V(bohr^3) \t    V/V0 \t    B(GPa) \t   rel. err.  \n"
            agl_data.outstr = agl_data.outstr + " -------------------------------------------------------------------------------------------- \n"
            for k in xrange (agl_data.npressure):
                agl_data.outstr = agl_data.outstr + '  ' + str(agl_data.pressure[k]).rjust(6) + "\t" + str(g[k]).rjust(10)[:10] + "\t" + str(agl_data.voleqmin[k]).rjust(10)[:10] + "\t" + str(agl_data.voleqmin[k]/vol0pres).rjust(8)[:8] + "\t" + str(agl_data.bulkmod[k]).rjust(10)[:10] + "\t" + str(round(rerr[k],10)).rjust(12)[:12] + "\n"
      
        #
        #.....Finite temperature numerical dynamic results.
        #
        if (agl_data.ieos == 0):
            numer (volref, nepol, epol, nfpol, fpol, False, agl_data)
        #
        #.....Vinet equation of state.
        #
        elif (agl_data.ieos == 1):
            vinet (vol0pres, gfe0pres, False, agl_data)
        #
        #.....Birch-Murnaghan equation of state.
        #
        elif (agl_data.ieos == 2):
            iG = 2
            birch (vol0pres, gfe0pres, iG, False, agl_data)
        #
        #.....Vinet equation of state with numerical calculation of properties.
        #
        elif (agl_data.ieos == 3):
            vinet (vol0pres, gfe0pres/agl_data.hy2kjmol, False, agl_data)
            numer (volref, nepol, epol, nepol, epol, False, agl_data)
        #
        #.....Birch-Murnaghan equation of state with numerical calculation of properties.
        #
        elif (agl_data.ieos == 4):
            iG = 2
            birch (vol0pres, gfe0pres/agl_data.hy2kjmol, iG, False, agl_data)
            numer (volref, nepol, epol, nepol, epol, False, agl_data)
        #
        #.....Equation of state of Baonza-Caceres-Nunez-Taravillo.
        #
        elif (agl_data.ieos == 5):
            bcnt (vol0pres, gfe0pres/agl_data.hy2kjmol, binp_bcnt, False, agl_data)
            bcnt_xk_sp.append(agl_data.xsupa);
            bcnt_xp_sp.append(agl_data.pspin);
            bcnt_xv_sp.append(agl_data.vspin);
            bcnt_xg_sp.append(agl_data.beta);
        #
        #.....Equation of state of Baonza-Caceres-Nunez-Taravillo, 
        #     but numerical calculation to obtain all properties.
        #
        elif (agl_data.ieos == 6):
            bcnt (vol0pres, gfe0pres/agl_data.hy2kjmol, binp_bcnt, False, agl_data)
            numer (volref, nepol, epol, nepol, epol, False, agl_data)
            bcnt_xk_sp.append(agl_data.xsupa);
            bcnt_xp_sp.append(agl_data.pspin);
            bcnt_xv_sp.append(agl_data.vspin);
            bcnt_xg_sp.append(agl_data.beta);
        #
        #.....Save properties at P=0 for all of the temperatures.
        #
        vol0p.append(agl_data.voleqmin[0])
        g0p.append(g[0])
        b0p.append(agl_data.bu0)
        bp0p.append(agl_data.bu1)
        bpp0p.append(agl_data.bu2)
        agl_data.gfe0p.append(g[0])
        g0pev = (g[0] / agl_data.hy2kjmol) * agl_data.hart2ev
        agl_data.gfe0pev.append(g0pev)
        #
        #.....Compute vibrational properties at the equilibrium volumes.
        #
        agl_data.outstr = agl_data.outstr + "\n"
        agl_data.outstr = agl_data.outstr + "VIBRATIONAL PROPERTIES \n"
        agl_data.outstr = agl_data.outstr + "====================== \n"
        agl_data.outstr = agl_data.outstr + "\n"
        agl_data.outstr = agl_data.outstr + "  P(GPa) \t  U(kJ/mol) \t Cv(J/mol*K) \t A(kJ/mol) \t S(J/mol*K) \t  Theta(K) \t gamma \n" 
        agl_data.outstr = agl_data.outstr + " ------------------------------------------------------------------------------------------------------- \n"
        k = 0
#        for k in xrange(agl_data.npressure):
        while k < agl_data.npressure:
            if (agl_data.idebye == 1 or agl_data.idebye == 2):
                tmp = math.log (agl_data.voleqmin[k])
                plnv = polin0 (tmp, ntpol, tpol)
                theta[k] = math.exp(plnv)
                plnv = polin1 (tmp, ntpol, tpol)
                agl_data.gamma_G[k] = - plnv
                thermal (theta[k], agl_data.temperature[j], D, Derr, agl_data)
                Uvib = agl_data.thermal_energ
                Cv[k] = agl_data.thermal_cv
                helm = agl_data.thermal_helmholtz
                ent = agl_data.thermal_entropy
                vibu=Uvib*hy2kjmol
                vibumev=Uvib*hart2ev*1000
                vibcv=Cv.at(k)*hy2kjmol*1000
                vibcvu=vibcv*kj2unit
                vibhe=helm*hy2kjmol
                vibhemev = helm*hart2ev*1000
                vibent=ent*hy2kjmol*1000
                vibentmev=ent*hart2ev*1000
                vibentu=vibent*kj2unit
                plnv = polin0(xmin,nepol,epol)
                vibg = plnv + helm + ((agl_data.pressure[k]/agl_data.au2gpa) * (xmin**3))
                vibgev = vibg * hart2ev
                vibg = vibg * hy2kjmol
            else:
                # Calculate and save Gruneisen parameter using derivative as a check
                if (k == 0):
                    ntpol = debfitt (tpol, agl_data)
                    tmp = math.log (agl_data.voleqmin[k]);
                    plnv = polin1 (tmp, ntpol, tpol)
                    agl_data.gamma_poly.append(-plnv)
                #
                #.....if udyn is not consistent, make a Gruneisen interpolation
                #
                if (agl_data.udyn[k] <= 0.0):
                    agl_data.logstr = agl_data.logstr + "MP AGL: inconsistent derivative, v[k] = " + str(agl_data.voleqmin[k]) + ", p[k] = " + str(agl_data.pressure[k]) + " \n"
                    agl_data.logstr = agl_data.logstr + "MP AGL: k = ", k, ", T = " + str(agl_data.temperature[j]) + ", udyn[k] = " + str(agl_data.udyn[k]) + " \n"
                    agl_data.logstr = agl_data.logstr + "MP AGL: making a Gruneisen interpolation \n"
                    ilow = 0
                    while (agl_data.voleqmin[k] >= agl_data.vol_inp[ilow] and ilow < (agl_data.ndata-1)):
                        ilow = ilow + 1	    
                    if ((ilow == (agl_data.ndata-1) and agl_data.voleqmin[k] >= agl_data.vol_inp[ilow]) or ilow == 0):
                        if (k == 0):
                            agl_data.logstr = agl_data.logstr + "MP AGL: inconsistent volume, v[k] = " + str(agl_data.voleqmin[k]) + ", p[k] = " + str(agl_data.pressure[k]) + " \n"
                            agl_data.grerr = 3
                            return
                        else: 
                            agl_data.logstr = agl_data.logstr + "MP AGL: Inconsistent volume, v[k] = " + str(agl_data.voleqmin[k]) + ", p[k] = " + str(agl_data.pressure[k]) + ", k = " + str(k) + " \n"
                            agl_data.logstr = agl_data.logstr + "MP AGL: Resetting maximum pressure value to p = " + str(agl_data.pressure[k-1]) + " \n"
                            agl_data.logstr = agl_data.logstr + "MP AGL: Resetting number of pressure points to " + str(k) + " \n"
                            agl_data.npressure = k   
                            break     
                    ilow = ilow - 1
                    agl_data.udyn[k] = agl_data.uder[ilow] * ((agl_data.voleqmin[k]/agl_data.vol_inp[ilow])**((math.log(agl_data.uder[ilow+1]/agl_data.uder[ilow])) / (math.log(agl_data.vol_inp[ilow+1]/agl_data.vol_inp[ilow]))))
                #
                #.....isotropic Debye model properties
                #
                theta[k] = ((6*agl_data.pi*agl_data.pi*agl_data.natoms*agl_data.voleqmin[k]*agl_data.voleqmin[k])**agl_data.third) / agl_data.pckbau * agl_data.poratio * math.sqrt(agl_data.udyn[k]/agl_data.cellmass)
                thermal(theta[k], agl_data.temperature[j], D, Derr, agl_data)
                Uvib = agl_data.thermal_energ
                Cv[k] = agl_data.thermal_cv
                helm = agl_data.thermal_helmholtz
                ent = agl_data.thermal_entropy
                vibu=Uvib*agl_data.hy2kjmol
                vibumev=Uvib*agl_data.hart2ev*1000
                vibcv=Cv[k]*agl_data.hy2kjmol*1000
                vibcvu=vibcv*agl_data.kj2unit
                vibhe=helm*agl_data.hy2kjmol
                vibhemev=helm*agl_data.hart2ev*1000
                vibent=ent*agl_data.hy2kjmol*1000
                vibentmev=ent*agl_data.hart2ev*1000
                vibentu=vibent*agl_data.kj2unit
                polin0(xmin,nepol,epol)
                vibg = plnv + helm + ((agl_data.pressure[k]/agl_data.au2gpa) * (xmin**3))
                vibgev = vibg * agl_data.hart2ev
                vibg = vibg * agl_data.hy2kjmol
            agl_data.outstr = agl_data.outstr + '  ' + str(agl_data.pressure[k]).rjust(6) + "\t " + str(vibu).rjust(10)[:10] + "\t  " + str(vibcv).rjust(10)[:10] + "\t" + str(vibhe).rjust(10)[:10] + "\t " + str(vibent).rjust(10)[:10] + "\t" + str(theta[k]).rjust(10)[:10] + "\t" + str(agl_data.gamma_G[k]).rjust(6)[:6] + "\n"
            if (k == 0):
                # Save heat capacity, internal energy, entropy, Debye temperature, Gruneisen parameter, Helmholtz free energy and Gibbs free energy at zero pressure for analysis and plotting later on
                agl_data.cv0p.append(vibcv)
                agl_data.cvu0p.append(vibcvu)
                agl_data.uvib0p.append(vibu)
                agl_data.uvib0pmev.append(vibumev)
                agl_data.ent0p.append(vibent)
                agl_data.ent0pmev.append(vibentmev)
                agl_data.ent0pu.append(vibentu)
                agl_data.tdeb0p.append(theta[k])
                agl_data.ga0p.append(agl_data.gamma_G[k])
                agl_data.he0p.append(vibhe)
                agl_data.he0pmev.append(vibhemev)
                agl_data.gvfe0p.append(vibg)
                agl_data.gvfe0pev.append(vibgev)
            if (agl_data.savpres):
                # Save heat capacity, internal energy, entropy, Debye temperature, Gruneisen parameter, and Helmholtz free energy at all pressures for analysis and plotting later on
                # Saving these values is optional and can be activated by setting 
                agl_data.cvp[j][k] = vibcv
                agl_data.cvup[j][k] = vibcvu
                agl_data.uvibp[j][k] = vibu
                agl_data.uvibpmev[j][k] = vibumev
                agl_data.entp[j][k] = vibent
                agl_data.entpmev[j][k] = vibentmev
                agl_data.entup[j][k] = vibentu
                agl_data.tdebp[j][k] = theta[k]
                agl_data.gap[j][k] = agl_data.gamma_G[k]
                agl_data.hep[j][k] = vibhe
                agl_data.hepmev[j][k] = vibhemev
            k = k + 1
        #
        #.....Compute EOS derivatives
        #
        agl_data.outstr = agl_data.outstr + "\n"
        agl_data.outstr = agl_data.outstr + "THERMAL EOS DERIVATIVES \n"
        agl_data.outstr = agl_data.outstr + "======================= \n"
        agl_data.outstr = agl_data.outstr + "\n"
        agl_data.outstr = agl_data.outstr + "  P(GPa) \t alpha(10^5/K) \t dp/dt(GPa/K) \t Bs(GPa) \t Cp(J/mol*K) \n"
        agl_data.outstr = agl_data.outstr + " ---------------------------------------------------------------------------- \n" 
        for k in xrange(agl_data.npressure):
            Pbeta=Cv[k]*agl_data.gamma_G[k]/agl_data.voleqmin[k] * agl_data.au2gpa
            agl_data.alpha[k]=Pbeta/agl_data.bulkmod[k]
            tmp = 1.0 + agl_data.gamma_G[k]*agl_data.alpha[k]*agl_data.temperature[j]
            Cp = Cv[k] * tmp * agl_data.hy2kjmol * 1000
            Bs = agl_data.bulkmod[k] * tmp
            agl_data.outstr = agl_data.outstr + '  ' + str(agl_data.pressure[k]).rjust(6) + "\t" + str(round(agl_data.alpha[k]*1e5,12)).rjust(14)[:14] + "\t " + str(round(Pbeta,10)).rjust(12)[:12] + "\t" + str(Bs).rjust(8)[:8] + "\t" + str(Cp).rjust(12)[:12] + "\n"
            if (k == 0):
                agl_data.a0p.append(agl_data.alpha[k])
                agl_data.bs0p.append(Bs)
                agl_data.cp0p.append(Cp)
        agl_data.outstr_thermal = agl_data.outstr_thermal + str(agl_data.temperature[j]).rjust(6) + "\t" + str(agl_data.uvib0pmev[j]).rjust(10)[:10] + "\t" + str(agl_data.he0pmev[j]).rjust(10)[:10] + "\t" + str(agl_data.ent0pu[j]).rjust(10)[:10] + "\t" + str(agl_data.cvu0p[j]).rjust(10)[:10] + "\t" + str(agl_data.tdeb0p[j]).rjust(10)[:10] + "\t" + str(agl_data.ga0p[j]).rjust(8)[:8] + "\n"
        j = j + 1
    #
    #.....Close temperature loop
    #
    #.....Write all of the results at P=0 for all of the temperatures to the stringstream for the main MP AGL output file.
    #
    if (agl_data.ieos >= 0):
        agl_data.outstr = agl_data.outstr + "\n"
        agl_data.outstr = agl_data.outstr + "RESULTS AT P=0 FOR ALL TEMPERATURES \n"
        agl_data.outstr = agl_data.outstr +  "=================================== \n" 
        agl_data.outstr = agl_data.outstr + "\n"
        agl_data.outstr = agl_data.outstr +  "    T(K) \t V(bohr^3) \t G(kJ/mol) \t U(kJ/mol) \t S(J/mol K) \t Cv(J/mol K) \n"
        agl_data.outstr = agl_data.outstr +  " -------------------------------------------------------------------------------------------- \n"
        for j in xrange(agl_data.ntemperature):
            agl_data.outstr = agl_data.outstr + '  ' + str(agl_data.temperature[j]).rjust(6) + "\t" + str(vol0p[j]).rjust(10)[:10] + "\t" + str(g0p[j]).rjust(10)[:10] + "\t" + str(agl_data.uvib0p[j]).rjust(10)[:10] + "\t " + str(agl_data.ent0p[j]).rjust(10)[:10] + "\t  " + str(agl_data.cv0p[j]).rjust(10)[:10] + "\n"
        agl_data.outstr = agl_data.outstr + "\n"
        agl_data.outstr = agl_data.outstr + "OTHER THERMODYNAMIC PROPERTIES AT P=0 \n"
        agl_data.outstr = agl_data.outstr + "===================================== \n"
        agl_data.outstr = agl_data.outstr + "\n"
        agl_data.outstr = agl_data.outstr + "    T(K) \t B0(GPa) \t     B0' \t B0''(GPa-1) \t   Bs(GPa) \t alpha(10^5/K) \n" 
        agl_data.outstr = agl_data.outstr + " -------------------------------------------------------------------------------------------------- \n"
        for j in xrange(agl_data.ntemperature):
            agl_data.outstr = agl_data.outstr + '  ' + str(agl_data.temperature[j]).rjust(6) + "\t" + str(b0p[j]).rjust(8)[:8] + "\t" + str(bp0p[j]).rjust(8)[:8] + "\t  " + str(bpp0p[j]).rjust(10)[:10] + "\t" + str(agl_data.bs0p[j]).rjust(10)[:10] + "\t  " + str(agl_data.a0p[j]*100000.0).rjust(12)[:12] + "\n"
        if (agl_data.ieos == 5 or agl_data.ieos == 6):
            agl_data.outstr = agl_data.outstr + "\n"
            agl_data.outstr = agl_data.outstr + "SPINODAL PROPERTIES \n" 
            agl_data.outstr = agl_data.outstr + "=================== \n" 
            agl_data.outstr = agl_data.outstr + "\n"
            agl_data.outstr = agl_data.outstr + "    T(K) \t     K* \t  Psp(GPa) \t Vsp(bohr^3) \t    beta \n" 
            agl_data.outstr = agl_data.outstr + " ------------------------------------------------------------------------ \n" 
            for j in xrange(agl_data.ntemperature):
                agl_data.outstr = agl_data.outstr + '  ' + str(agl_data.temperature[j]).rjust(6) + "\t" + str(bcnt_xk_sp[j]).rjust(10)[:10] + "\t" + str(-bcnt_xp_sp[j]).rjust(10)[:10] + "\t  " + str(bcnt_xv_sp[j]).rjust(10)[:10] + "\t  " + str(bcnt_xg_sp[j]).rjust(6)[:6] + "\n"
        agl_data.outstr = agl_data.outstr + "\n"
        agl_data.outstr = agl_data.outstr + "DEBYE MODEL RELATED PROPERTIES AT P=0 \n"
        agl_data.outstr = agl_data.outstr + "===================================== \n"
        agl_data.outstr = agl_data.outstr + "\n"
        agl_data.outstr = agl_data.outstr + "    T(K) \t Theta(K) \t gamma \n"
        agl_data.outstr = agl_data.outstr + " -------------------------------------- \n" 
        for j in xrange(agl_data.ntemperature):
            agl_data.outstr = agl_data.outstr + '  ' + str(agl_data.temperature[j]).rjust(6) + "\t" + str(agl_data.tdeb0p[j]).rjust(9)[:9] + "\t" + str(agl_data.ga0p[j]).rjust(6)[:6] + "\n"
    #
    #.....end of the program
    #
    agl_data.logstr = agl_data.logstr + "MP AGL: End of run ok! \n"
    return

