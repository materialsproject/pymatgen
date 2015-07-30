#!/usr/bin/env python

"""
This module calculates the Debye temperature self-consistently (optional), fits the Debye temperature, and calculates thermal properties from the Debye temperature.
"""

from __future__ import division
import warnings
import sys
import subprocess
import unittest
import pymatgen
from pymatgen.agl_thermal.agl_polynomial import polfit
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





# **************************************************************************************
#  This set of functions determine the Debye temperature self-consistently
# **************************************************************************************

#
#.....scdebye - calculates the table of self consistent Debye
#     temperatures at T for all the input volumes.
#
#     If firsttime is true, only calculates the static pressures table
#
# Adapted from original Fortran version written by M. A. Blanco et al.
# See Computer Physics Communications 158, 57-72 (2004) and Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details
#
def scdebye (T, pol, err, npol, imin, firsttime, pst, agl_data): 
    eps=5.0
    mloops=25
    scerr = False
    #
    #.....static pressures?
    #
    if (firsttime):
        for i in xrange(agl_data.ndata):
            plnv = polin1 (agl_data.xconfigvector[i], npol, pol)
            pst[i] = -agl_data.xconfigvector[i] * plnv / (3.0*agl_data.vol_inp[i])
        return scerr
    #
    #.....loop until convergence is achieved
    #
    converged = False
    iloops = 0
    itry = imin
    while (((not converged) or (iloops < mloops)) and (iloops < agl_data.maxloops)):
        iloops = iloops + 1
        for i in xrange(agl_data.ndata):
            thermal (agl_data.tdebye[i], T, D, Derr, agl_data)
            agl_data.datatofit[i] = agl_data.energ_inp[i] + agl_data.thermal_helmholtz
        imin = minbrack (imin, agl_data)
        if (agl_data.ierr != 0):
            agl_data.logstr = agl_data.logstr + "MP AGL scdebye: T = "  + str(T) + ", minimum point = " + str(imin) + ", trial point = " + str(itry) + ", total points = " + str(agl_data.ndata) + " \n"
            agl_data.logstr = agl_data.logstr + "MP AGL scdebye: func = " + str(agl_data.datatofit) + " \n"
            scerr = True
            return scerr
        npol = polynomial_fit (imin, pol, err, agl_data)
        if (agl_data.fterr != 0):
            agl_data.logstr = agl_data.logstr + "MP AGL scdebye: No polynomial fit found with a minimum within bounds of input data \n"
            agl_data.logstr = agl_data.logstr + "MP AGL scdebye: T = " + str(T) + ", minimum point = " + str(imin) + ", trial point = " + str(itry) + ", total points = " + str(agl_data.ndata) + " \n"
            scerr = True
            return scerr
        #
        #.....calculate the new set of debye temperatures and its convergence
        #
        converged = True
        theta0 = agl_data.tdebye[0]
        for i in xrange(agl_data.ndata):
            thermal (agl_data.tdebye[i], T, D, Derr, agl_data)
            U = agl_data.thermal_energ
            Cvt = agl_data.thermal_cv
            f1 = polin1 (agl_data.xconfigvector[i], npol, pol)
            f2 = polin2 (agl_data.xconfigvector[i], npol, pol)
            p = -agl_data.xconfigvector[i] * f1 / (3.0*agl_data.vol_inp[i])
            bt = -agl_data.xconfigvector[i] / (9.0*agl_data.vol_inp[i]) * (2.0*f1 - agl_data.xconfigvector[i]*f2)
            gm = (p-pst[i]) * agl_data.vol_inp[i] / U
            bsv = agl_data.vol_inp[i] * bt + T * gm*gm * Cvt
            if (bsv < 0.0):
                if (i <= imin):
                    agl_data.logstr = agl_data.logstr + "MP AGL scdebye: Bs < 0 for equilibrium V! \n"
                    return scerr
                if (i > 0):
                    theta = agl_data.tdebye[i]/theta0 * agl_data.tdebye[i-1]
                dt = 0.0
            else: 
                theta = ((6.0*pi*pi*agl_data.natoms/agl_data.vol_inp[i])**agl_data.third) / agl_data.pckbau * agl_data.poratio * math.sqrt(bsv/agl_data.cellmass)
                dt = theta - agl_data.tdebye[i]
                if (i > 0):
                    if (theta > agl_data.tdebye[i-1]):
                        theta = agl_data.tdebye[i]/theta0 * agl_data.tdebye[i-1]
                        agl_data.logstr = agl_data.logstr + "MP AGL scdebye: Warning! gamma < 0, i = " + str(i) + ", T = " + str(T) + " \n"
                        dt = 0.0
            theta0 = agl_data.tdebye[i]
            agl_data.tdebye[i] = theta;
            if (agl_data.vol_inp[i]-agl_data.vol_inp[imin] < agl_data.vol_inp[imin]*0.1 and i >= max(3,agl_data.ndata/10)):
                converged = converged and (math.fabs(dt) < eps)

    #
    #....warn of convergence failure
    #
    if (iloops >= agl_data.maxloops and (not converged)):
        agl_data.logstr = agl_data.logstr + "MP AGL scdebye: Maximum number of convergence iterations exceeded \n"
        agl_data.logstr = agl_data.logstr + "MP AGL scdebye: Number of convergence loops = " + str(iloops) + " \n"
    if (not converged):
        agl_data.logstr = agl_data.logstr + "MP AGL scdebye: Warning! convergence not achieved \n"
    #
    #.....end
    #
    return scerr


#
#.....debfitt - fitts Ln theta - Ln V data to polynomials, and
#     averages them.
#
#     It returns the averaged polynomial coefficients in polynomialcoeffs.
#
#
# Adapted from original Fortran version written by V. Luana and M. A. Blanco, Departamento de Quimica Fisica y Analitica, Universidad de Oviedo, 33006-Oviedo, Spain.  
# See Computer Physics Communications 158, 57-72 (2004) and Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details
#
def debfitt (polynomialcoeffs, agl_data):
    polycoeffswork = []
    weight = []
    x = []
    y = []
    #
    #.....promedium mechanism;
    #
    npolycoeffsfit = []
    ndatafit = []
    #
    #.....service variables:
    #
    #
    #.....limit in the dimensions of polynomials (# of points - limit)
    #
    limit = 4
    #
    #.....number of pairs of data points to delete in the fitting procedure
    #
    ndel = 3
    agl_data.fterr = 0
    #
    #.....initialize the weights and the variables
    #
    for i in xrange (agl_data.ndata):
        x.append(math.log(agl_data.vol_inp[i]))
        y.append(math.log(agl_data.tdebye[i]))
        weight.append(1.0)
    for i in xrange (agl_data.mpar+1):
        polycoeffswork.append(0.0)
    for i in xrange(agl_data.mfit+1):
        npolycoeffsfit.append(0.0)
        ndatafit.append(0.0)
    xdata_tofit = []
    ydata_tofit = []
    weight_tofit = []
    #
    #.....fitting loops through different polynomials:
    #
    nfit = 0
    npolynomialcoeffs = 0
    npolycoeffsmin = 2
    iinf = 0
    isup = agl_data.ndata - 1
    ndatamin = max (agl_data.ndata-2*ndel, npolycoeffsmin+1)
    ndataact = agl_data.ndata
    rmsmin = 1e33
    #
    #.....loop eliminating pairs of outermost elements ndel times:
    #
    while (ndataact >= ndatamin):
        #
        #.....loop increasing the number of parameters until limit is reached:
        #
        npolycoeffsmax = min (ndataact-limit-1, agl_data.mpar)
        npolycoeffs = npolycoeffsmin
        while npolycoeffs <= npolycoeffsmax:
            nfit = nfit + 1
            del xdata_tofit[:]
            del ydata_tofit[:]
            del weight_tofit[:]
            i = iinf
            while i < (isup + 1):
                xdata_tofit.append(x[i])
                ydata_tofit.append(y[i])
                weight_tofit.append(weight[i])
                i = i + 1
            rms, polycoeffswork = polfit (iinf, isup, x, y, weight, npolycoeffs)
            #
            #.....discard fits over the maximum declared dimensions
            #
            if (nfit > agl_data.mfit):
                nfit = agl_data.mfit
                agl_data.logstr = agl_data.logstr + "MP AGL debfitt: Warning! maximum number of fits exceeded \n"
            #
            #.....save fit parameters
            #
            else:
                npolynomialcoeffs = max (npolynomialcoeffs,npolycoeffs)
                npolycoeffsfit[nfit] = npolycoeffs
                ndatafit[nfit] = ndataact
                rmsmin = min(rmsmin,rms*npolycoeffs/ndataact)
            npolycoeffs = npolycoeffs + 1
        iinf = iinf + 1
        isup = isup - 1
        ndataact = ndataact - 2
    #
    #.....number of fits control
    #
    if (nfit == 0):
        agl_data.logstr = agl_data.logstr + "MP AGL fitt: no fits to average! \n"
        agl_data.fterr = 1
        return npolynomialcoeffs
    #
    #.....average the polynomial coefficients (repeating the fits)
    #
    wnorm = 0.0
    for i in xrange(npolynomialcoeffs+1):
        polynomialcoeffs.append(0.0)
    ifit = 1
    while ifit <= nfit:
        ndataact = ndatafit[ifit]
        npolycoeffs = npolycoeffsfit[ifit]
        iinf = int(1 + ((agl_data.ndata-ndataact) / 2) - 1)
        isup = agl_data.ndata - iinf - 1
        del xdata_tofit[:]
        del ydata_tofit[:]
        del weight_tofit[:]
        i = iinf
        while i < (isup + 1):
            xdata_tofit.append(x[i])
            ydata_tofit.append(y[i])
            weight_tofit.append(weight[i])
            i = i + 1
        rms, polycoeffswork = polfit (iinf, isup, x, y, weight, npolycoeffs)
        wtmp = rms*npolycoeffs/(rmsmin*ndataact)
        wtmp = math.exp(-wtmp*wtmp)
        for i in xrange(npolycoeffs + 1):
            polynomialcoeffs[i] = polynomialcoeffs[i] + wtmp * polycoeffswork[i];
        wnorm = wnorm + wtmp
        ifit = ifit + 1
    #
    #.....put the proper weight into the polynomials
    #
    for i in xrange(npolynomialcoeffs + 1):
        polynomialcoeffs[i] = polynomialcoeffs[i] / wnorm
    #
    #.....end of routine
    #
    return npolynomialcoeffs







# **************************************************************************************
#  This set of functions compute the Debye model vibrational properties
# **************************************************************************************

#-----------------------------------------------------------------------
#
#.....thermal - compute Debye model vibrational properties.
#
#     This routine obtains the molar vibrational properties of a given
#     crystal by means of the Debye model: internal energy (U), heat
#     capacity at constant volume (Cv), Helmholtz's free energy (F),
#     and vibrational entropy (S).
#
#     To evaluate this properties, the following integral is needed:
#                         
#                                      |    x^3     |
#     Debye (y) = 3*y^(-3) * INT (0,y) | ---------- | dx
#                                      | exp(x) - 1 |
#
#    where y=ThetaD/T, being ThetaD Debye's temperature (K) and T the
#     absolute (thermodynamic) temperature. The integral is evaluated
#     using a Gauss-Legendre quadrature.
#
#-----INPUT-------------------------------------------------------------
#       ThetaD : Debye's temperature (K).     
#            T : Absolute temperature (K).
#       natoms : Number of atoms in the unit cell.
#-----OUTPUT-------------------------------------------------------------
#           en : Vibrational internal energy, U (hartree/molecule).
#           cv : Constant V heat capacity, Cv (hartree/K molecule).
#           he : Helmholtz's free energy (hartree/molecule).
#          ent : Entropy (hartree/K molecule).
#        Debye : Debye's integral.
#         xabs : Maximum error in Debye's integral evaluation.
#----------------------------------------------------------------------------
#
# Adapted from original Fortran version written by M. A. Blanco et al.
# See Computer Physics Communications 158, 57-72 (2004) and Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details
#
def thermal (ThetaD, T, debye, xabs, agl_data):
    eps=1e-12
    cero=0.0;
    maxnl=100;
    tpi=3.14159265358979323846
    x = [0.0 for i in range(maxnl)]
    w = [0.0 for i in range(maxnl)]
    tol2=1e-8
    #
    #.....error condition controls
    #
    debye=cero
    xabs=cero
    if (ThetaD < 0.0 or T < 0.0):
        agl_data.logstr = agl_data.logstr + "MP AGL thermal: Bad ThetaD or T value \n"
        return
    if (T < tol2):
        agl_data.thermal_energ = agl_data.natoms*9.0*agl_data.pckbau*ThetaD/8.0
        agl_data.thermal_cv = cero
        agl_data.thermal_helmholtz = agl_data.natoms*9.0*agl_data.pckbau*ThetaD/8.0
        agl_data.thermal_entropy = cero
        return
    if (ThetaD < tol2):
        agl_data.thermal_energ = agl_data.natoms*3.0*agl_data.pckbau*T
        agl_data.thermal_cv = agl_data.natoms*3.0*agl_data.pckbau
        agl_data.thermal_entropy = 1e100
        agl_data.thermal_helmholtz = agl_data.thermal_energ - T * agl_data.thermal_entropy
        return
    y=ThetaD/T
    debye=3.0*tpi*tpi*tpi*tpi/y/y/y/15.0
    if (y <= 250):
        #
        #.....Loop with increasing number of Legendre points to evaluate the Debye integral.
        #
        debye0=1e30
        nl = 5
        xabs=math.fabs(debye-debye0)
        while ((nl <= maxnl) and (xabs >= eps)):
            gauleg (cero, y, x, w, nl, agl_data)
            sum=0.0
            i = 0
            while (i < nl):
                sum=sum+w[i]*fdebye(x[i])
                i = i + 1
            debye=sum*3.0/y/y/y
            xabs=math.fabs(debye-debye0)
            debye0=debye
            nl = nl + 5
    #
    #.....thermodynamic vibrational properties
    #     energ = internal energy, cv = heat capacity, entropy = entropy, helmholtz = Helmholtz vibrational free energy
    #
    agl_data.thermal_energ  = agl_data.natoms * 3.0 * agl_data.pckbau * (ThetaD*3.0/8.0 + T*debye)
    agl_data.thermal_cv  = agl_data.natoms * 3.0 * agl_data.pckbau * (4.0*debye - 3.0*y/(math.exp(y)-1.0))
    agl_data.thermal_entropy = agl_data.natoms * 3.0 * agl_data.pckbau * (debye*4.0/3.0 - math.log(1.0-math.exp(-y)))
    agl_data.thermal_helmholtz  = agl_data.thermal_energ - T * agl_data.thermal_entropy
    # End of function
    return

def fdebye(z):
    return z*z*z / (math.exp(z) - 1.0)



#
#.....gauleg - Gauss-Legendre quadrature coefficients.
#
#.....Given the lower and upper limits of integration x1 and x2, 
#     and given n, this routine returns arrays x and w of length n, 
#     containing the abscissas and weights of the Gauss-Legendre 
#     n-point quadrature formula.
#-----------------------------------------------------------------------
#.....High precision is a good idea for this routine
#
# Adapted from original Fortran version written by M. A. Blanco et al., which was based on the Numerical Recipes routine
# See Computer Physics Communications 158, 57-72 (2004) and Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details
#    
def gauleg(x1, x2, x, w, n, agl_data):
    dpi = 3.141592653589793238462643
#.....Increase eps if you don't have this floating precision.
    eps=3.0e-16
    #
    #.....The roots are symmetric in the interval, so we only have to find
    #     half of them.
    #
    m=int((n+1)/2)
    xm=0.5*(x2+x1)
    xl=0.5*(x2-x1)
    #
    #.....Loop over the desired roots.
    for i in xrange(m):
        idoub = i + 1
        z=math.cos(dpi*(idoub-0.25)/(n+0.5))
        # Initialize convergence test
        zdiff = 1.0
        iloops = 0
        while ((zdiff > eps) and (iloops < agl_data.maxloops)):
            iloops = iloops + 1
            p1=1.0
            p2=0.0
            #
            #...Loop up the recurrence relation to get the Legendre polynomial 
            #   evaluated at z.
            j = 1
            while(j <= n): 
                p3=p2
                p2=p1
                p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j
                j = j + 1
            #
            #.........p1 is now the desired Legendre polynomial. We next compute pp,
            #         derivative , by a standard relation involving p2, the polyn-
            #         omial of one lower order.
            #
            pp=n*(z*p1-p2)/(z*z-1.0)
            z1=z
            #
            #.........Newton's method.
            #
            z=z1-p1/pp
            # Convergence test
            zdiff = math.fabs(z-z1)
        if ((iloops >= agl_data.maxloops) and (zdiff > eps)):
            agl_data.logstr = agl_data.logstr + "MP AGL gauleg: Maximum number of convergence loops exceeded \n"
            agl_data.logstr = agl_data.logstr + "MP AGL gauleg: Number of loops for convergence = " + str(iloops) + " \n"
        #.......Scale the root to the desired interval.
        x[i] = xm - xl * z
        #     
        #.......and put in its symmetric counterpart.
        x[n-idoub] = xm + xl * z
        #
        #......compute the weight.
        w[i]=2.0*xl/((1.0-z*z)*pp*pp)
        #.......and its symmetric counterpart.
        w[n-idoub]=w[i]
    return



