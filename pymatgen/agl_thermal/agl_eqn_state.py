#!/usr/bin/env python

"""
This module calculates thermal properties using different equations of state.
"""

from __future__ import division
import warnings
import sys
import subprocess
import unittest
import pymatgen
from pymatgen.agl_thermal.agl_polynomial import polfit
from pymatgen.agl_thermal.agl_polynomial import polin0
from pymatgen.agl_thermal.agl_polynomial import polin1
from pymatgen.agl_thermal.agl_polynomial import polin2
from pymatgen.agl_thermal.agl_polynomial import polin3
from pymatgen.agl_thermal.agl_polynomial import polin4
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
__date__ = "April 8, 2014"


# **************************************************************************************
#  These functions calculate the thermal properties using different equations of state
# **************************************************************************************

#
#.....numer - numerical EOS calculation.
#
#.....Numer computes the derivatives of the Helmholtz function and the
#     static energy needed to obtain Debye's temperature, the static
#     pressure, and succesive derivatives of the bulk modulus.
#
# Adapted from original Fortran version written by M. A. Blanco et al.
# See Computer Physics Communications 158, 57-72 (2004) and Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details
#
def numer (volref, nepol, epol, nfpol, fpol, statcalc, agl_data):
    #
    #.....Compute Pfit(P), B(P), B'(P), and B''(P)
    #
    if (agl_data.ieos >= 0):
        agl_data.outstr = agl_data.outstr + '\n'
        agl_data.outstr = agl_data.outstr + 'NUMERICAL EOS PRESSURE DERIVATIVES \n'
        agl_data.outstr = agl_data.outstr + '================================== \n' 
        agl_data.outstr = agl_data.outstr + "  P(GPa) \t V(bohr^3) \t    V/V0 \t Pfit(GPa) \t    B(GPa) \t    B' \t  B''(GPa-1) \n"
        agl_data.outstr = agl_data.outstr + ' ------------------------------------------------------------------------------------------------------ \n'
    for k in xrange(agl_data.npressure):
        xeqmin = (agl_data.voleqmin[k]/volref)**agl_data.third
        f1 = polin1 (xeqmin, nfpol, fpol)
        f2 = polin2 (xeqmin, nfpol, fpol)
        f3 = polin3 (xeqmin, nfpol, fpol)
        f4 = polin4 (xeqmin, nfpol, fpol)
        pt = -xeqmin * f1 / (3.0*agl_data.voleqmin[k]) * agl_data.au2gpa
        tmp = 2.0 * f1 - xeqmin * f2
        agl_data.bulkmod[k] = -xeqmin / (9.0*agl_data.voleqmin[k]) * tmp * agl_data.au2gpa
        tmp2 = (f2 - xeqmin * f3) / tmp
        b1 = agl_data.third * (2.0 - xeqmin * tmp2)
        b2 = -agl_data.voleqmin[k] * (tmp2*(1.0-xeqmin*tmp2) - xeqmin*xeqmin*f4/tmp) / (agl_data.au2gpa*tmp)
        if (k == 0):
            agl_data.bu0 = agl_data.bulkmod[k]
            agl_data.bu1 = b1
            agl_data.bu2 = b2
        if (agl_data.ieos >= 0):
            agl_data.outstr = agl_data.outstr + '  ' + str(agl_data.pressure[k]).rjust(6) + '\t' + str(agl_data.voleqmin[k]).rjust(10)[:10] + '\t' + str(agl_data.voleqmin[k]/agl_data.voleqmin[0]).rjust(8)[:8] + '\t' + str(pt).rjust(10)[:10] + '\t' + str(agl_data.bulkmod[k]).rjust(10)[:10] + '\t' + str(b1).rjust(6)[:6] + '\t     ' + str(b2).rjust(7)[:7] + '\n'
    #
    #.....Static calculation: get second derivative of static energy
    #
    if (statcalc):
        if (agl_data.ieos >= 0):
            agl_data.outstr = agl_data.outstr +'\n'
            agl_data.outstr = agl_data.outstr + 'INPUT AND FITTED VALUES OF THE LATTICE ENERGY \n'
            agl_data.outstr = agl_data.outstr + '============================================= \n'
            agl_data.outstr = agl_data.outstr + '\n'
            agl_data.outstr = agl_data.outstr + '   V(bohr^3)     E_inp(hartree)     E_fit(hartree) \n'
            agl_data.outstr = agl_data.outstr + ' --------------------------------------------------\n'
        for i in xrange(agl_data.ndata):
            f0 = polin0 (agl_data.xconfigvector[i], nepol, epol)
            f1 = polin1 (agl_data.xconfigvector[i], nepol, epol)
            f2 = polin2 (agl_data.xconfigvector[i], nepol, epol)
            tmp = agl_data.xconfigvector[i] * f2 - 2.0 * f1
            v3 = 3.0 * agl_data.vol_inp[i]
            agl_data.uder.append(tmp * agl_data.xconfigvector[i] / (v3*v3))
            if (agl_data.ieos >= 0):
                agl_data.outstr = agl_data.outstr + '  ' + str(agl_data.vol_inp[i]).rjust(10)[:10] + '\t ' + str(agl_data.energ_inp[i]).rjust(14)[:14] + '\t    ' + str(f0).rjust(14)[:14] + '\n'
        #
        #.....Dynamic calculation: get static pressure and second derivative of the energy
        #
    else:
        for k in xrange(agl_data.npressure):
            xeqmin = (agl_data.voleqmin[k]/volref)**agl_data.third
            f1 = polin1 (xeqmin, nepol, epol)
            f2 = polin2 (xeqmin, nepol, epol)
            f3 = polin3 (xeqmin, nepol, epol)
            f4 = polin4 (xeqmin, nepol, epol)
            v3 = 3.0 * agl_data.voleqmin[k]
            agl_data.pstatic[k] = -f1*agl_data.au2gpa * xeqmin / v3
            tmp = xeqmin * f2 - 2.0 * f1;
            agl_data.udyn[k] = tmp * xeqmin / (v3*v3)
            tmp2 = f2 - xeqmin * f3;
            agl_data.gamma_G[k] = (1.0 + xeqmin * tmp2 / tmp)/6.0
    return




#...................................................................
#.....vinet - computes Vinet EOS from (P,V) data.
#
#.....VINET computes the EOS from the (P,V) data. The EOS has the
#     following expresion:
#     log H = A + B(1-x)
#     being H = Px**2/(3(1-x))
#           A = log Bo
#           B = 3/2((Bo)'-1)
#           X = (V/Vo)**(1/3)
#
# Adapted from original Fortran version written by M. A. Blanco et al.
# See Computer Physics Communications 158, 57-72 (2004) and Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details
#...................................................................
def vinet (vol0pres, gfe0pres, statcalc, agl_data):
    #
    #.....fit Log H vs. (1-x)
    #
    logh = [0.0 for k in range(agl_data.npressure)]
    x = [0.0 for k in range(agl_data.npressure)]
    db = []
    d2b = []
    sumz=0.0
    sumy=0.0
    sumzy=0.0
    sumz2=0.0
    sumy2=0.0
    n=0
    x[0] = 1.0
    i = 1
    while (i < agl_data.npressure):
        x[i] = (agl_data.voleqmin[i]/vol0pres)**agl_data.third
        h = agl_data.pressure[i]*x[i]*x[i]/(3.0*(1-x[i]))
        logh[i] = math.log(h)
        z = 1-x[i]
        n=n+1
        sumz = sumz + z
        sumy = sumy + logh[i]
        sumzy = sumzy + z*logh[i]
        sumz2 = sumz2 + z*z
        sumy2 = sumy2 + logh[i]*logh[i]
        i = i + 1
    lnb0=(sumy*sumz2 - sumzy*sumz)/(n*sumz2 - sumz*sumz)
    A=(n*sumzy - sumz*sumy)/(n*sumz2 - sumz*sumz)
    raiz=math.sqrt((sumz2 - sumz*sumz/n)*(sumy2 - sumy*sumy/n))
    rfit=(sumzy-sumz*sumy/n)/raiz
    logh[0] = lnb0
    #
    #.....obtain B0, B0', B0''
    #
    agl_data.bu0=math.exp(lnb0)
    agl_data.bu1=2.0*A*agl_data.third+1.0
    agl_data.bu2 = -(2.0+A*(A+6.0))/(9.0*agl_data.bu0)
    #
    #.....save static values
    #
    if (statcalc):
        agl_data.g00k = gfe0pres
        agl_data.v00k = vol0pres
        agl_data.b00k = agl_data.bu0/agl_data.au2gpa
        agl_data.A00k = A
    #
    #.....Compute Pfit(P), B(P), B'(P), and B''(P)
    #
    agl_data.bulkmod[0]=agl_data.bu0
    db.append(agl_data.bu1)
    d2b.append(agl_data.bu2)
    agl_data.pfit[0]=0.0
    i = 1
    while (i < agl_data.npressure):
        a1x = A * (1.0 - x[i])
        ax1 = A * x[i] + 1.0
        f0x = x[i] * (1.0-a1x) - 2.0
        f1x = ax1 - a1x
        f2x = 2.0 * A
        f1f0 = f1x / f0x
        f2f0 = f2x / f0x
        fnw = 1.0 - x[i] * f1f0
        x2inv = 1.0 / (x[i]*x[i])
        b0exp = agl_data.bu0 * math.exp(a1x)
        agl_data.bulkmod[i] = -b0exp * f0x * x2inv
        db.append(agl_data.third * (ax1+fnw))
        d2b.append(x[i]/(9.0*agl_data.bulkmod[i]) * (x[i]*f2f0 - A + f1f0*fnw))
        agl_data.pfit[i] = 3.0 * (1.0-x[i]) * x2inv * b0exp
        i = i + 1
    #
    #.....output
    #
    agl_data.outstr = agl_data.outstr + "\n"
    agl_data.outstr = agl_data.outstr + "VINET EOS PRESSURE DERIVATIVES \n"
    agl_data.outstr = agl_data.outstr + "============================== \n"
    agl_data.outstr = agl_data.outstr + "\n"
    agl_data.outstr = agl_data.outstr + "  1-V/V0 \t Vinet-Func \t P(GPa) \t Pfit(GPa) \t    B(GPa) \t        B' \t  B''(GPa-1) \n"
    agl_data.outstr = agl_data.outstr + " ------------------------------------------------------------------------------------------------------------ \n" 
    for i in xrange(agl_data.npressure):
        agl_data.outstr = agl_data.outstr + '  ' + str(1.0-x[i]).rjust(6)[:6] + "\t " + str(logh[i]).rjust(10)[:10] + "\t " + str(agl_data.pressure[i]).rjust(6)[:6] + "\t        " + str(agl_data.pfit[i]).rjust(10)[:10] + "\t" + str(agl_data.bulkmod[i]).rjust(10)[:10] + "\t" + str(db[i]).rjust(10)[:10] + "\t  " + str(d2b[i]).rjust(10)[:10] + "\n"
    agl_data.outstr = agl_data.outstr + "\n"
    agl_data.outstr = agl_data.outstr + "B0 = " + str(agl_data.bu0) + ", B0' = " + str(agl_data.bu1) + ", B0'' = " + str(agl_data.bu2)  + " reg.coef = " + str(rfit) + "\n"
    agl_data.outstr = agl_data.outstr + "\n"

    #
    #.....Static calculation: get static energy and its second derivative
    #
    if (statcalc):
        for i in xrange(agl_data.ndata):
            x00k = (agl_data.vol_inp[i]/agl_data.v00k)**agl_data.third
            a1x = agl_data.A00k * (1.0 - x00k)
            f0x = x00k * (1.0-a1x) - 2.0
            b0exp = agl_data.b00k * math.exp(a1x)
            agl_data.ust.append(agl_data.g00k + 9.0*agl_data.v00k/(agl_data.A00k*agl_data.A00k) * (b0exp*(a1x-1.0)+agl_data.b00k))
            agl_data.uder.append(-f0x / (x00k*x00k*agl_data.vol_inp[i]) * b0exp)
        #
        #.......Print input and fitted values of the lattice energy.
        #
        agl_data.outstr = agl_data.outstr + "\n"
        agl_data.outstr = agl_data.outstr + "INPUT AND FITTED VALUES OF THE LATTICE ENERGY \n"
        agl_data.outstr = agl_data.outstr + "============================================= \n"
        agl_data.outstr = agl_data.outstr + "\n"
        agl_data.outstr = agl_data.outstr + "   V(bohr^3)     E_inp(hartree)     E_fit(hartree) \n" 
        agl_data.outstr = agl_data.outstr + " -------------------------------------------------- \n"
        for i in xrange(agl_data.ndata):
            agl_data.outstr = agl_data.outstr + '  ' + str(agl_data.vol_inp[i]).rjust(10)[:10] + "\t " + str(agl_data.energ_inp[i]).rjust(14)[:14] + "\t    " + str(agl_data.ust[i]).rjust(14)[:14] + "\n"
        #
        #.....Dynamic calculation: get static pressure and second derivative
        #     of the energy
        #
    else:
        for i in xrange(agl_data.npressure):
            x00k = (agl_data.voleqmin[i]/agl_data.v00k)**agl_data.third
            a1x = agl_data.A00k * (1.0 - x00k)
            ax1 = agl_data.A00k * x00k + 1.0
            f0x = x00k * (1.0-a1x) - 2.0
            f1x = ax1 - a1x
            f2x = 2.0 * agl_data.A00k
            f1f0 = f1x / f0x
            fnw = 1.0 - x00k * f1f0
            x2inv = 1.0 / (x00k*x00k)
            b0exp = agl_data.b00k * math.exp(a1x)
            agl_data.pstatic[i] = 3.0 * (1.0-x00k) * x2inv * b0exp * agl_data.au2gpa
            agl_data.udyn[i] = -f0x * x2inv * x2inv / (x00k*agl_data.v00k) * b0exp
            agl_data.gamma_G[i] = (ax1 + fnw - 1.0)/6.0
    return





#-----------------------------------------------------------------------------
#.....birch - computes the Birch-Murnaghan EOS of order iG from the
#     (P,V) data.
#
#     The EOS has the following expression:
#
#     F = Sum (i=0,iG) a(i)*f^i
# 
#     being : F = P/[3f(1+2f)^(5/2)]
#             f = [x^(-2)-1]/2
#             x = [V(i)/V(1)]^(1/3)
#
#-----INPUT
#  vol0pres      : Molecular volume (bohr^3/mol) at P=0.
#  gfe0pres      : Gibbs energy (or 0k static energy) at v0 (hartree).
#  iG      : order of the fitting.
#  press() : Pressure values (GPa). common /eos/.
#  vinp()  : Initial values of the volume (bohr^3/mol). common /input/.
#  statcalc: Logical variable that determines if the calculation is
#            static or dynamic. In the first case the second derivative
#            of the static energy (uder) is computed for all the input
#            values of the volume. In the second case the second
#            derivative of the static energy (udyn) is computed for
#            the equilibrium volumes at the different pressures.
#
#-----OUTPUT
#  pstatic() : Static pressures in GPa (only on dynamic calculations).
#  uder()  : Second derivative of ust(k) for each vinp(). Hy/bohr^6
#  udyn()  : Second derivative of ust(k) for each V(). Hy/bohr^6
#  rms     : Root mean square deviation.
#  bu0,bu1,bu2 : Bulk modulus and their derivatives at P=0.
#
#.....The output is stored in common /eos/
#
# Adapted from original Fortran version written by M. A. Blanco et al.
# See Computer Physics Communications 158, 57-72 (2004) and Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details
#-----------------------------------------------------------------------
def birch (vol0pres, gfe0pres, iG, statcalc, agl_data):
    tol = 1e-12
    npresm2 = agl_data.npressure - 2
    izero = 0
    if (iG > agl_data.maiG):
        agl_data.logstr = agl_data.logstr + "MP AGL birch : Too high fitting order \n"
        agl_data.brerr = 1
        return
    if (math.fabs(agl_data.pressure[0]) > tol):
        agl_data.logstr = agl_data.logstr + "MP AGL birch : P(0) must be 0.0 \n"
        agl_data.brerr = 1
        return
    acoef = [0.0 for i in range(agl_data.maiG+1)]
    fstr = []
    ybir = []
    weight = []
    db = []
    d2b = []
    #
    #.....Compute the Birch function F and strain variable f.
    #
    weight.append(0.0)
    i = 1
    while (i < agl_data.npressure):
        rr0 = (agl_data.voleqmin[i]/vol0pres)**agl_data.third
        fstr.append((rr0**(-2)-1)/2.0)
        ybir.append(agl_data.pressure[i]/agl_data.au2gpa/(3*fstr[i-1]*((1+2*fstr[i-1]**2.5))))
        weight.append(1.0)
        i = i + 1
    #
    #.....Fitting to a polynomial of order iG.
    #
    rms, acoef = polfit (izero, npresm2, fstr, ybir, weight, iG)
    #
    #.....Compute B0,B0',B0''.
    #
    agl_data.bu0=acoef[0]*agl_data.au2gpa
    if (iG == 0):
        agl_data.bu1=4.0
        agl_data.bu2=-35.0/(9.0*agl_data.bu0)
    elif (iG == 1): 
        agl_data.bu1=4.0+2.0*acoef[1]*agl_data.au2gpa/(3.0*agl_data.bu0)
        agl_data.bu2=(-agl_data.bu1*(agl_data.bu1-7.0)-143.0/9.0)/agl_data.bu0
    elif (iG >= 2):
        agl_data.bu1=4.0+2.0*acoef[1]*agl_data.au2gpa/(3.0*agl_data.bu0);
        agl_data.bu2=(2.0*acoef[2]/(agl_data.bu0*3.0)-agl_data.bu1*(agl_data.bu1-7.0)-143.0/9.0)/agl_data.bu0
    #
    #.....Compute B(P), B'(P), and B''(P). (b(), db(), and d2b().
    #
    for i in xrange(agl_data.npressure):
        if (i == 0):
            agl_data.pfit[i]=0.0
            agl_data.bulkmod[i]=agl_data.bu0;
            db.append(agl_data.bu1)
            d2b.append(agl_data.bu2)
        else: 
            st=fstr[i-1]
            stsq=math.sqrt(1.0+2.0*st)
            s2=stsq*stsq
            st32=stsq*s2
            st52=st32*s2
            pol0 = polin0(st, iG, acoef)
            pol1 = polin1(st, iG, acoef)
            pol2 = polin2(st, iG, acoef)
            pol3=0.0
            if (iG > 2):
                pol3 = polin3(st, iG, acoef)
            #
            #.........Fitted pressure and B(P).
            #
            agl_data.pfit[i]=3.0*st*st52*pol0*agl_data.au2gpa
            sum1=st32*(st*s2*pol1+(1.0+7.0*st)*pol0)
            agl_data.bulkmod[i]=s2*sum1
            sum2=st52*(s2*pol1+2*st*pol1+st*s2*pol2+7*pol0+(1+7*st)*pol1)
            den=3*st*st52*pol1+(3.0*st52+15.0*st*st32)*pol0
            #
            #.........B'(P).
            #
            db.append((5*sum1+sum2)/den)
            d2bdf2=25*stsq*(st*s2*pol1+(1.0+7.0*st)*pol0)
            d2bdf2=d2bdf2+10.0*st32*((2.0+11.0*st)*pol1+7*pol0+st*s2*pol2)
            d2bdf2=d2bdf2+st52*((3.0+15.0*st)*pol2+18*pol1+st*s2*pol3)
            d2pdf2=3*st52*pol1+15*st*st32*pol1+3*st*st52*pol2
            d2pdf2=d2pdf2+(30*st32+45*st*stsq)*pol0
            d2pdf2=d2pdf2+(3*st52+15*st*st32)*pol1
            #
            #.........B''(P).
            #
            d2b.append((den*d2bdf2-(5*sum1+sum2)*d2pdf2)/(den**3))
            agl_data.bulkmod[i]=agl_data.bulkmod[i]*agl_data.au2gpa
            d2b[i]=d2b[i]/agl_data.au2gpa
    #
    #.....Output.
    #
    agl_data.outstr = agl_data.outstr + "\n"
    agl_data.outstr = agl_data.outstr + "BIRCH-MURNAGHAN EOS PRESSURE DERIVATIVES \n"
    agl_data.outstr = agl_data.outstr + "======================================== \n"
    agl_data.outstr = agl_data.outstr + "\n"
    agl_data.outstr = agl_data.outstr + "  Strain \t Birch-Func \t P(GPa) \t Pfit(GPa) \t    B(GPa) \t        B' \t B''(GPa-1) \n"
    agl_data.outstr = agl_data.outstr + " ----------------------------------------------------------------------------------------------------------- \n"
    for i in xrange(agl_data.npressure):
        if (i == 0):
            agl_data.outstr = agl_data.outstr + '  ' + str(0.0).rjust(6)[:6] + "\t " + str(agl_data.bu0/agl_data.au2gpa).rjust(10)[:10] + "\t " + str(agl_data.pressure[i]).rjust(6)[:6] + "\t    " + str(agl_data.pfit[i]).rjust(14)[:14] + "\t" + str(agl_data.bulkmod[i]).rjust(10)[:10] + "\t" + str(db[i]).rjust(10)[:10] + "\t " + str(d2b[i]).rjust(10)[:10] + "\n"
        else: 
            agl_data.outstr = agl_data.outstr + '  ' + str(fstr[i-1]).rjust(6)[:6] + "\t " + str(ybir[i-1]).rjust(10)[:10] + "\t " + str(agl_data.pressure[i]).rjust(6)[:6] + "\t    " + str(agl_data.pfit[i]).rjust(14)[:14] + "\t" + str(agl_data.bulkmod[i]).rjust(10)[:10] + "\t" +  str(db[i]).rjust(10)[:10] + "\t " +  str(d2b[i]).rjust(10)[:10] + "\n"
    agl_data.outstr = agl_data.outstr + "\n"
    agl_data.outstr = agl_data.outstr + "B0 = " + str(agl_data.bu0) + ", B0' = " + str(agl_data.bu1) + ", B0'' = " + str(agl_data.bu2) + ", reg.coef = " + str(rms) + "\n"
    agl_data.outstr = agl_data.outstr + "\n"
    if (statcalc):
        #
        #.......Compute the static potential energy U(V) and its second
        #       derivative U''(V) with respect to V for all the input
        #       values of the volume.
        #
        agl_data.v00k=vol0pres
        agl_data.g00k=gfe0pres
        for k in xrange(iG + 1): 
            agl_data.astatic.append(acoef[k])
        for k in xrange(agl_data.ndata):
            agl_data.ust.append(agl_data.g00k)
            agl_data.uder.append(0.0)
            st=(agl_data.vol_inp[k]/agl_data.v00k)**agl_data.third
            st=((st**(-2))-1)/2.0
            s2=(1.0+2.0*st)
            pol0 = polin0(st, iG, agl_data.astatic)
            pol1 = polin1(st, iG, agl_data.astatic)
            v9=9.0*agl_data.v00k
            for j in xrange(iG + 1): 
                agl_data.ust[k]=agl_data.ust[k]+v9*agl_data.astatic[j]/(j+2)*(st**(j+2)) 
            agl_data.uder[k]=s2*s2*s2*s2/agl_data.v00k*(st*s2*pol1+(1.0+7.0*st)*pol0)
        #
        #.......Print input and fitted values of the lattice energy.
        #
        agl_data.outstr = agl_data.outstr + "\n"
        agl_data.outstr = agl_data.outstr + "INPUT AND FITTED VALUES OF THE LATTICE ENERGY \n"
        agl_data.outstr = agl_data.outstr + "============================================= \n"
        agl_data.outstr = agl_data.outstr + "\n"
        agl_data.outstr = agl_data.outstr + "   V(bohr^3)     E_inp(hartree)     E_fit(hartree) \n"
        agl_data.outstr = agl_data.outstr + " -------------------------------------------------- \n"
        for i in xrange(agl_data.ndata):
            agl_data.outstr = agl_data.outstr + '  ' +  str(agl_data.vol_inp[i]).rjust(10)[:10] + "\t " + str(agl_data.energ_inp[i]).rjust(14)[:14] + "\t    " + str(agl_data.ust[i]).rjust(14)[:14] + "\n"
        return
    else: 
        #
        #.......Compute the second derivative U''(V) with respect to V 
        #       for all the equilibrium values of the volume at the 
        #       different pressures.
        #
        for k in xrange(agl_data.npressure):
            st=(agl_data.voleqmin[k]/agl_data.v00k)**agl_data.third
            st=((st)**(-2)-1)/2.0
            s2=(1.0+2.0*st)
            s22=s2*s2
            pol0 = polin0(st, iG, agl_data.astatic)
            pol1 = polin1(st, iG, agl_data.astatic)
            pol2 = polin2(st, iG, agl_data.astatic)
            pol3=0.0
            if (iG > 2):
                pol3 = polin3(st, iG, agl_data.astatic) 
            agl_data.pstatic[k]=agl_data.au2gpa*3.0*st*(s2**2.5)*pol0
            tmp = (1.0+7.0*st)*pol0 + st*s2*pol1
            agl_data.udyn[k] = s22*s22 / agl_data.v00k * tmp
            tmp = 1.0 / tmp
            tmp2 = s2*tmp * (7.0*pol0 + (2.0+11.0*st)*pol1 + st*s2*pol2)
            v3 = agl_data.voleqmin[k] / (3.0*agl_data.v00k)
            agl_data.gamma_G[k] = -2.0*agl_data.third + 0.5*s2*math.sqrt(s2)*v3*(8.0+tmp2)
    #
    #.....end
    #
    return




#...................................................................
#
#.....bcnt - compute the Spinodal (BCNT) EOS from (B,p) data.
#
#     The EOS has the following expresion:
#
#                       g
#     B(p) = ( p - Psp ) / K
#
#     where//c
#     g = 0.85  (If opt_g = .true. ===> g is optimized)
#     (-Psp) and K are the parameter to optimize.
#     
#     These parameters bear the following relation with Bo and Bo'.
#                 g  -1
#     Bo  = (-Psp)  K
#
#                      -1
#     Bo' = g Bo (-Psp)
#
#-----Input parameters:
#     lg     : Logical unit for results output.
#     vol0pres     : Zero pressure volume, either static or dynamic.
#     gfe0pres     : Zero pressure Gibbs function.
#     B0     : Bulk modulus used to compute the initial value of
#              -Psp (GPa).
#     opt_g  : if .true.  ==> g is optimized.
#     static : if .true.  ==> static calculation.
#
# Adapted from original Fortran version written by M. A. Blanco et al.
# See Computer Physics Communications 158, 57-72 (2004) and Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details
# ...................................................................
def bcnt (vol0pres, gfe0pres, b0, statcalc, agl_data):
    eps=1e-10 
    maxnl = 100
    tol = 1e-12
    agl_data.volp0 = vol0pres
    db = []
    d2b = []
    xg = [0.0 for i in range(maxnl)]
    wg = [0.0 for i in range(maxnl)]
    #
    #.....Initial values of properties to optimize.
    #
    agl_data.gbao = 0.85
    if (not statcalc):
        if (agl_data.iopt_g == 2):
            agl_data.gbao = agl_data.gbao0
    x_Psp = agl_data.gbao*b0/4.0
    #
    #.....Optimize g and x_Psp.
    #
    ax = x_Psp*0.5
    bx = x_Psp
    cx = x_Psp*2.0
    ax, bx, cx, fa, fb, fc = mnbrak (ax, bx, cx, agl_data)
    x_Psp, desv = brent (ax, bx, cx, tol, x_Psp, agl_data)
    #
    #.....Final properties.
    #
    xx = math.exp((agl_data.gbao-1)*math.log(x_Psp))
    agl_data.xmopt = agl_data.xkopt/xx/(1.0-agl_data.gbao)
    agl_data.bu0 = xx * x_Psp / agl_data.xkopt
    agl_data.bu1 = agl_data.gbao * agl_data.bu0 / x_Psp
    agl_data.bu2 = agl_data.gbao * (agl_data.gbao - 1.0) * agl_data.bu0 / x_Psp / x_Psp
    vsp = vol0pres * math.exp(agl_data.gbao/(1-agl_data.gbao)/agl_data.bu1)
    agl_data.pspin = x_Psp
    agl_data.xsupa = agl_data.xkopt
    agl_data.vspin = vsp
    agl_data.beta = agl_data.gbao
    #
    #.....save static values
    #
    if (statcalc):
        agl_data.g00k = gfe0pres
        agl_data.b00k = agl_data.bu0/agl_data.au2gpa
        agl_data.v00k = vol0pres
        agl_data.vsp0k = vsp    
        agl_data.xkopt0 = agl_data.xkopt     
        agl_data.xmopt0 = agl_data.xmopt     
        agl_data.x_Psp0 = x_Psp     
        agl_data.gbao0 = agl_data.gbao     
    #
    #.....Compute Pfit(P), B(P), B'(P), and B''(P)
    #
    for i in xrange(agl_data.npressure):
        xxx = (agl_data.xmopt + math.log (vol0pres/agl_data.voleqmin[i]))/agl_data.xmopt
        ug = 1.0/(1-agl_data.gbao)
        agl_data.pfit[i] = x_Psp * (math.exp(ug*math.log(xxx)) - 1.0);
        xdu = math.exp((agl_data.gbao-1)*math.log(agl_data.pressure[i]+x_Psp));
        agl_data.bulkmod[i] = xdu * (agl_data.pressure[i]+x_Psp) / agl_data.xkopt;
        db.append(agl_data.gbao * xdu / agl_data.xkopt);
        d2b.append(agl_data.gbao * (agl_data.gbao - 1) * xdu / agl_data.xkopt / (agl_data.pressure[i]+x_Psp));
    #
    #.....output
    #
    agl_data.outstr = agl_data.outstr + "\n"
    agl_data.outstr = agl_data.outstr + "SPINODAL EOS PRESSURE DERIVATIVES \n"
    agl_data.outstr = agl_data.outstr + "================================= \n" 
    agl_data.outstr = agl_data.outstr + "\n"
    agl_data.outstr = agl_data.outstr + "  P(GPa) \t Pfit(GPa) \t    B(GPa) \t      B' \t B''(GPa-1) \n" 
    agl_data.outstr = agl_data.outstr + " --------------------------------------------------------------------------- \n" 
    for i in xrange(agl_data.npressure): 
        agl_data.outstr = agl_data.outstr + '  ' + str(agl_data.pressure[i]).rjust(6)[:6] + "\t" + str(agl_data.pfit[i]).rjust(10)[:10] + "\t" + str(agl_data.bulkmod[i]).rjust(10)[:10] +  "\t" + str(db[i]).rjust(8)[:8] + "\t " + str(d2b[i]).rjust(10)[:10] + "\n"
    agl_data.outstr = agl_data.outstr + "\n"
    if (agl_data.opt_g):
        agl_data.outstr = agl_data.outstr + "B0 = " + str(agl_data.bu0) + ", B0' = " + str(db[0]) + ", B0'' = " + str(d2b[0]) + ", reg.coef = " + str(desv) + ", Vsp = " + str(vsp) + ", -Psp = " + str(x_Psp) + ", K* = " + str(agl_data.xkopt) + ", gamma = " + str(agl_data.gbao) + "\n"
    else: 
        agl_data.outstr = agl_data.outstr + "B0 = " + str(agl_data.bulkmod[0]) + ", B0' = " + str(db[0]) + ", B0'' = " + str(d2b[0]) + ", reg.coef = " + str(desv) + ", Vsp = " + str(vsp) + ", -Psp = " + str(x_Psp) + ", K* = " + str(agl_data.xkopt) + ", gamma = " + str(agl_data.gbao) + "\n"
    #
    #.....Static calculation: get static energy and its second derivative.
    #
    if (statcalc):
        for i in xrange(agl_data.ndata):
            x = (agl_data.xmopt0+math.log(agl_data.v00k/agl_data.vol_inp[i]))/agl_data.xmopt0
            auxg = 1.0/(1.0-agl_data.gbao0)
            auxg2 = agl_data.gbao0/(1.0-agl_data.gbao0)
            #
            #.........Compute numerically the integrated Helmholtz function by means
            #         of a loop with increasing number of Legendre points.
            #
            xinf = 1.0;
            xsup = x;
            factor = 1.0;
            if (xsup < xinf):
                aux = xinf
                xinf = xsup
                xsup = aux
                factor = -1.0
            #
            #.........Iterative loop.
            #
            sum0=1e30
            nl = 5 
            xabs=1.0
            while ((nl <= maxnl) and (xabs >= eps)):
                gauleg (xinf, xsup, xg, wg, nl, agl_data)
                sum=0.0
                for ii in xrange(nl):
                    term = math.exp(agl_data.xmopt0*(1.0-xg[ii])) * (math.exp(auxg*math.log(xg[ii])) - 1.0)
                    sum = sum + wg[ii] * factor * term * agl_data.xmopt0 * agl_data.v00k * agl_data.x_Psp0	
                xabs = math.fabs(sum-sum0)
                sum0 = sum
                nl = nl + 5
            agl_data.ust.append(agl_data.g00k + sum / agl_data.au2gpa)
            agl_data.uder.append(agl_data.x_Psp0 / (agl_data.xmopt0 * agl_data.vol_inp[i] * (1.0 - agl_data.gbao0)) * math.exp(auxg2*math.log(x)) / agl_data.au2gpa)  
        #
        #.......Print input and fitted values of the lattice energy.
        #
        agl_data.outstr = agl_data.outstr + "\n"
        agl_data.outstr = agl_data.outstr + "INPUT AND FITTED VALUES OF THE LATTICE ENERGY \n"
        agl_data.outstr = agl_data.outstr + "============================================= \n"
        agl_data.outstr = agl_data.outstr + "\n"
        agl_data.outstr = agl_data.outstr + "   V(bohr^3)     E_inp(hartree)     E_fit(hartree) \n"
        agl_data.outstr = agl_data.outstr + " -------------------------------------------------- \n" 
        for i in xrange(agl_data.ndata): 
            agl_data.outstr = agl_data.outstr + '  ' + str(agl_data.vol_inp[i]).rjust(10)[:10] + "\t " + str(agl_data.energ_inp[i]).rjust(14)[:14] + "\t    " + str(agl_data.ust[i]).rjust(14)[:14] + "\n"
    #
    #.....Dynamic calculation: get static pressure and second derivative
    #     of the energy
    #
    else:
        for i in xrange(agl_data.npressure): 
            xxx = (agl_data.xmopt0 + math.log (agl_data.v00k/agl_data.voleqmin[i]))/agl_data.xmopt0
            ug = 1.0/(1.0-agl_data.gbao0);
            auxg2 = agl_data.gbao0/(1.0-agl_data.gbao0);
            agl_data.pstatic[i] = agl_data.x_Psp0 * (math.exp(ug*math.log(xxx)) - 1.0)
            agl_data.udyn[i] = agl_data.x_Psp0 / (agl_data.xmopt0 * agl_data.v00k * (1.0 - agl_data.gbao0)) * math.exp(-agl_data.xmopt0*(1.0-xxx)) * math.exp(auxg2*math.log(xxx)) / agl_data.au2gpa
            agl_data.gamma_G[i] = -1.0/6.0 + agl_data.gbao0 * ug / 2.0/ agl_data.xmopt0 / xxx      
    #
    #.....end
    #
    return





# **************************************************************************************
#  This set of functions implement routines required for the BCNT EOS
# **************************************************************************************

#
#.....mnbrak - brackets a minimum of the function f.
#
#     Given a function, and two distinct initial points ax and bx,
#     this routine searches in the downhill direction (defined by the
#     function as evaluated at the initial points) and returns new
#     points ax, bx, and cx which bracket a minimum of the function.
#     Also returned are the function values at the three points: fa, fb,
#     and fc.
#
# Adapted from original Fortran version written by M. A. Blanco et al.
# See Computer Physics Communications 158, 57-72 (2004) and Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details
#
def mnbrak (ax, bx, cx, agl_data):
    gold = 1.618034
    glimit = 100.0
    tiny = 1e-20
    fa = optm(ax, agl_data)
    fb = optm(bx, agl_data)
    if (fb > fa):
        dum = ax
        ax = bx
        bx = dum
        dum = fb
        fb = fa
        fa = dum
    cx = bx+gold*(bx-ax)
    fc = optm(cx, agl_data)
    while ( fb >= fc ):
        endpart = True
        r = (bx-ax)*(fb-fc)
        q = (bx-cx)*(fb-fa)
        u = bx-((bx-cx)*q-(bx-ax)*r)/(2.0*math.copysign(max(math.fabs(q-r),tiny),q-r))
        ulim = bx+glimit*(cx-bx)
        if ((bx-u)*(u-cx) > 0.0):
            fu = optm(u, agl_data)
            if (fu < fc):
                ax = bx
                fa = fb
                bx = u
                fb = fu
                endpart = False
            elif (fu > fb):
                cx = u
                fc = fu
                endpart = False
            else:
                u = cx+gold*(cx-bx)
                fu = optm(u, agl_data)
        elif ((cx-u)*(u-ulim) > 0.0):
            fu = optm(u, agl_data)	  
            if (fu < fc):
                bx = cx
                cx = u
                u = cx+gold*(cx-bx)
                fb = fc
                fc = fu
                fu = optm(u, agl_data)
        elif ((u-ulim)*(ulim-cx) >= 0.0):
            u = ulim
            fu = optm(u, agl_data)
        else:
            u = cx+gold*(cx-bx)
            fu = optm(u, agl_data)
        if (endpart):
            ax = bx
            bx = cx
            cx = u
            fa = fb
            fb = fc
            fc = fu
    return ax, bx, cx, fa, fb, fc


#
#.....brent - unidimensional minimization of f in the range [ax,cx].
#
#    Given a function, and a bracketing triplet of abscissas this
#    routine isolates the minimum to a fractional precission of tol
#    using Brent's method. The bracketing triplet must be such that bx
#    is between ax and cx, and that f(bx) is less than both f(ax) and
#    f(cx). The abscissa of the minimum is returned as xmin, and the
#    minimum function value as brent.
#
# Adapted from original Fortran version written by M. A. Blanco et al. 
# See Computer Physics Communications 158, 57-72 (2004) and Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details
#
def brent (ax, bx, cx, tol, xmin, agl_data):
    itmax = 100
    cgold = 0.3819660;
    zeps = 1.0e-10;
    tol3 = 1e-12;
    notskip1 = True
    a = min(ax,cx)
    b = max(ax,cx)
    v = bx
    w = v
    x = v
    e = 0.0
    d = 0.0
    fx = optm(x, agl_data)
    fv = fx
    fw = fx
    iter = 1
    while (iter <= itmax): 
        xm = 0.5*(a+b)
        tol1 = tol*math.fabs(x)+zeps
        tol2 = 2.0*tol1
        if (math.fabs(x-xm) <= (tol2-0.5*(b-a))):
            xmin = x
            brentx = fx
            return xmin, brentx
        else:
            if (math.fabs(e) > tol1):
                r = (x-w)*(fx-fv)
                q = (x-v)*(fx-fw)
                p = (x-v)*q-(x-w)*r
                q = 2.0*(q-r)
                if (q > 0.0): 
                    p = -p
                q = math.fabs(q)
                etemp = e
                e = d
                if (math.fabs(p) >= math.fabs(0.5*q*etemp) or p <= q*(a-x) or p >= q*(b-x)):
                    notskip1 = True
                else:
                    d = p/q
                    u = x+d
                    if (u-a < tol2 or b-u < tol2):
                        d = math.copysign(tol1,xm-x)
                        if (math.fabs(d) >= tol1):
                            u = x+d
                        else: 
                            u = x + math.copysign(tol1, d)
                        notskip1 = False
            if (notskip1):
                if (x >= xm):
                    e = a-x
                else: 
                    e = b-x
                d = cgold*e;
            if (math.fabs(d) >= tol1): 
                u = x+d
            else:
                u = x + math.copysign(tol1, d)
            fu = optm(u, agl_data)
            if (fu <= fx): 
                if (u >= x): 
                    a = x
                else: 
                    b = x
                v = w
                fv = fw
                w = x
                fw = fx
                x = u
                fx = fu
            else:
                if (u < x): 
                    a = u
                else: 
                    b = u
            if (fu <= fw or math.fabs(w - x) < tol3 ):
                v = w
                fv = fw
                w = u
                fw = fu
            elif (fu <= fv or math.fabs(v - x) < tol3 or math.fabs(v - w) < tol3 ):
                v = u;
                fv = fu;
        iter = iter + 1
    agl_data.logstr = agl_data.logstr + "MP AGL brent: exceeded maximum iterations. \n"
    xmin = x
    brentx = fx
    return xmin, brentx

	  

#
#.....optm - optimization of exponent parameter "g" (aka "beta") required for BCNT EOS.
#
# Adapted from original Fortran version written by M. A. Blanco et al.  
# See Computer Physics Communications 158, 57-72 (2004) and Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details
#
def optm (x_Psp, agl_data):
    xfunc = [0.0 for k in range(agl_data.npressure)]
    yfunc = [0.0 for k in range(agl_data.npressure)]
    if (x_Psp < 0.0):
        agl_data.logstr = agl_data.logstr + "MP AGL optm: Warning: Spinodal pressure is negative \n"
        desv = 1e30
        return desv
    if (agl_data.opt_g):
        a11 = 0.0
        a12 = 0.0
        a21 = 0.0
        a22 = 0.0
        z1 = 0.0
        z2 = 0.0
        for i in xrange(agl_data.npressure):
            xfunc[i] = math.log(agl_data.pressure[i]+x_Psp)
            yfunc[i] = math.log(agl_data.bulkmod[i])
            a11 = a11 + 1.0
            a12 = a12 + xfunc[i]
            a21 = a12
            a22 = a22 + xfunc[i]*xfunc[i]
            z1 = z1 + yfunc[i]
            z2 = z2 + xfunc[i]*yfunc[i]
        det = a11 * a22 - a12 * a21
        x_in = (z1 * a22 - z2 * a12)/det
        x_de = (a11 * z2 - z1 * a21)/det
        desv = 0.0
        for i in xrange(agl_data.npressure):
            desv = desv + (yfunc[i] - x_in - x_de * xfunc[i])**2
        agl_data.xkopt = math.exp(-x_in)
        agl_data.gbao = x_de;      
        return desv
    else: 
        a12 = 0.0
        z1 = 0.0
        for i in xrange(agl_data.npressure):
            xfunc[i] = math.log (agl_data.pressure[i]+x_Psp);
            yfunc[i] = math.log(agl_data.bulkmod[i]);
            a12 = a12 + xfunc[i]
            z1 = z1 + yfunc[i]
        x_in = (z1-agl_data.gbao*a12)/agl_data.npressure
        agl_data.xkopt = math.exp(-x_in)
        desv = 0.0
        for i in xrange(agl_data.npressure):
            desv = desv + (yfunc[i] - x_in - agl_data.gbao * xfunc[i])**2
    return desv


