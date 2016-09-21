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


class eos_thermal_functions:
    def _init_():
        pass


    # **************************************************************************************
    #  This set of functions compute the Debye model vibrational properties
    # **************************************************************************************


    def thermal(self, ThetaD, T, natoms, pckbau, maxloops, logstr):
        """.....thermal - compute Debye model vibrational properties.
    
         This routine obtains the molar vibrational properties of a given
         crystal by means of the Debye model: internal energy (U), heat
         capacity at constant volume (Cv), Helmholtz's free energy (F),
         and vibrational entropy (S).
    
         To evaluate this properties, the following integral is needed:
    
                                          |    x^3     |
         Debye (y) = 3*y^(-3) * INT (0,y) | ---------- | dx
                                          | exp(x) - 1 |
    
        where y=ThetaD/T, being ThetaD Debye's temperature (K) and T the
         absolute (thermodynamic) temperature. The integral is evaluated
         using a Gauss-Legendre quadrature.
    
         -----INPUT-------------------------------------------------------------
           ThetaD : Debye's temperature (K).
                T : Absolute temperature (K).
           natoms : Number of atoms in the unit cell.
        -----OUTPUT-------------------------------------------------------------
            energ : Vibrational internal energy, U (hartree/molecule).
               cv : Constant V heat capacity, Cv (hartree/K molecule).
        helmholtz : Helmholtz's free energy (hartree/molecule).
          entropy : Entropy (hartree/K molecule).
            Debye : Debye's integral.
             xabs : Maximum error in Debye's integral evaluation.
        ----------------------------------------------------------------------------
    
        Adapted from original Fortran version written by M. A. Blanco et al.
        See Computer Physics Communications 158, 57-72 (2004) and Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details
        """
        eps = 1e-12
        cero = 0.0
        maxnl = 100
        tpi = 3.14159265358979323846
        xabscissa = [0.0 for i in range(maxnl)]
        weight = [0.0 for i in range(maxnl)]
        tol2 = 1e-8
        #
        # .....error condition controls
        #
        debye = cero
        xabs = cero
        if ThetaD < 0.0 or T < 0.0:
            logstr = logstr + "MP Eqn of State Thermal thermal: Bad ThetaD or T value \n"
            thermal_energ = 0.0
            thermal_cv = 0.0
            thermal_helmholtz = 0.0
            thermal_entropy = 0.0
            return thermal_energ, thermal_entropy, thermal_helmholtz, thermal_cv, debye, xabs
        if T < tol2:
            thermal_energ = natoms*9.0*pckbau*ThetaD/8.0
            thermal_cv = cero
            thermal_helmholtz = natoms*9.0*pckbau*ThetaD/8.0
            thermal_entropy = cero
            return thermal_energ, thermal_entropy, thermal_helmholtz, thermal_cv, debye, xabs
        if ThetaD < tol2:
            thermal_energ = natoms*3.0*pckbau*T
            thermal_cv = natoms*3.0*pckbau
            thermal_entropy = 1e100
            thermal_helmholtz = thermal_energ - T * thermal_entropy
            return thermal_energ, thermal_entropy, thermal_helmholtz, thermal_cv, debye, xabs
        y = ThetaD/T
        debye = 3.0*tpi*tpi*tpi*tpi/y/y/y/15.0
        if y <= 250:
            #
            # .....Loop with increasing number of Legendre points to evaluate the Debye integral.
            #
            debye0 = 1e30
            nl = 5
            xabs = math.fabs(debye-debye0)
            while (nl <= maxnl) and (xabs >= eps):
                xabscissa, weight, logstr = self.gauss_legendre(cero, y, xabscissa, weight, nl, maxloops, logstr)
                sum = 0.0
                i = 0
                while i < nl:
                    sum = sum+weight[i]*self.fdebye(xabscissa[i])
                    i = i + 1
                debye = sum*3.0/y/y/y
                xabs = math.fabs(debye-debye0)
                debye0 = debye
                nl = nl + 5
        #
        # .....thermodynamic vibrational properties
        #     energ = internal energy, cv = heat capacity, entropy = entropy, helmholtz = Helmholtz vibrational free energy
        #
        thermal_energ = natoms * 3.0 * pckbau * (ThetaD*3.0/8.0 + T*debye)
        thermal_cv = natoms * 3.0 * pckbau * (4.0*debye - 3.0*y/(math.exp(y)-1.0))
        thermal_entropy = natoms * 3.0 * pckbau * (debye*4.0/3.0 - math.log(1.0-math.exp(-y)))
        thermal_helmholtz = thermal_energ - T * thermal_entropy
        # End of function
        return thermal_energ, thermal_entropy, thermal_helmholtz, thermal_cv, debye, xabs

    def fdebye(self, z):
        return z*z*z / (math.exp(z) - 1.0)



    def gauss_legendre(self, xlower, xupper, xabscissa, weight, n, maxloops, logstr):
        """.....gauss_legendre - Gauss-Legendre quadrature coefficients.
        
        .....Given the lower and upper limits of integration xupper and xlower,
        and given n, this routine returns arrays xabscissa and weight of length n,
        containing the abscissas and weights of the Gauss-Legendre
        n-point quadrature formula.
        -----------------------------------------------------------------------
        .....High precision is a good idea for this routine
        
        Adapted from original Fortran version written by M. A. Blanco et al., which was based on the Numerical Recipes routine
        See Computer Physics Communications 158, 57-72 (2004) and Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details
        """
        dpi = 3.141592653589793238462643
    # .....Increase eps if you don't have this floating precision.
        eps = 3.0e-16
        #
        # .....The roots are symmetric in the interval, so we only have to find
        #     half of them.
        #
        midpoint = int((n+1)/2)
        xm = 0.5 * (xupper + xlower)
        xl = 0.5 * (xupper - xlower)
        #
        # .....Loop over the desired roots.
        for i in xrange(midpoint):
            idoub = i + 1
            z = math.cos(dpi*(idoub-0.25)/(n+0.5))
        # Initialize convergence test
            zdiff = 1.0
            iloops = 0
            while (zdiff > eps) and (iloops < maxloops):
                iloops = iloops + 1
                polynom1 = 1.0
                polynom2 = 0.0
                #
                # ...Loop up the recurrence relation to get the Legendre polynomial
                #   evaluated at z.
                j = 1
                while j <= n:
                    polynom3 = polynom2
                    polynom2 = polynom1
                    polynom1 = ((2.0 * j - 1.0) * z * polynom2 - (j - 1.0) * polynom3)/j
                    j = j + 1
                #
                # .........p1 is now the desired Legendre polynomial. We next compute pp,
                #         derivative , by a standard relation involving p2, the polyn-
                #         omial of one lower order.
                #
                polynomderiv = n * (z * polynom1 - polynom2)/(z * z - 1.0)
                z1 = z
                #
                # .........Newton's method.
                #
                z = z1 - polynom1 / polynomderiv
                # Convergence test
                zdiff = math.fabs(z - z1)
            if (iloops >= maxloops) and (zdiff > eps):
                logstr = logstr + "MP Eqn of State Thermal gauleg: Maximum number of convergence loops exceeded \n"
                logstr = logstr + "MP Eqn of State Thermal gauleg: Number of loops for convergence = " + str(iloops) + " \n"
            # .......Scale the root to the desired interval.
            xabscissa[i] = xm - xl * z
            #
            # .......and put in its symmetric counterpart.
            xabscissa[n-idoub] = xm + xl * z
            #
            # ......compute the weight.
            weight[i] = 2.0 * xl / ((1.0 - z * z)* polynomderiv * polynomderiv)
            # .......and its symmetric counterpart.
            weight[n-idoub] = weight[i]
        return xabscissa, weight, logstr




    def cvdebfit(self, Cv, tdmin, tdmax, nkb, temp, npoints, maxloops, logstr):
        """
        cvdebfit: determines which value of the Debye temperature produces
        best fit to heat capacity data

        Runs through all integer values of Debye temperature in range
        between minimum and maximum values

        Generates heat capacity curve for each Debye temperature using
        Debye model (see eqn 23.26, Solid State Physics, Ashcroft and Mermin):

        c_V = 9 n k_B (T / Theta_D)^3 \int_0^{Theta_D / T} (x^4 e^x) / (e^x - 1)^2 dx
        """
        itmin = int(tdmin)
        itmax = int(tdmax)
        rmsmin = 1e30
        i = itmin
        while i <= itmax:
            tdtrial = float(i)
            smsq = 0.0
            j = 1
            while j < npoints:
                y = tdtrial / temp[j]
                deb, logstr = self.debint(y, maxloops, logstr)
                cvt = 9.0 * nkb * deb
                smsq = smsq + (cvt - Cv[j])**2
                j = j + 1
            rms = math.sqrt(smsq)
            if rms < rmsmin:
                rmsmin = rms
                tdbest = tdtrial
            i = i + 1
        return tdbest, logstr


    def fdebyeint(self, z):
        """
        Evaluates Debye integral
        """
        integz = ((z**4) * math.exp(z)) / ((math.exp(z) - 1)**2)
        return integz


    def debint(self, y, maxloops, logstr):
        """
        debint: Evaluation of the Debye integral:

                                          |       x^4      |
         Debye (y) = 3*y^(-3) * INT (0,y) | -------------- | dx
                                          | (exp(x) - 1)^2 |

         where y=ThetaD/T, being ThetaD Debye's temperature (K) and T the
         absolute (thermodynamic) temperature. The integral is evaluated
         using a Gauss-Legendre quadrature.
         """
        eps = 1e-12
        cero = 0.0
        maxnl = 100
        pi = 3.14159265358979323846
        #
        # .....error condition controls
        #
        xabscissa = [0.0 for i in range(maxnl)]
        weight = [0.0 for i in range(maxnl)]
        debye = 3.0*pi*pi*pi*pi/y/y/y/15.0
        if y <= 250:
            #
            # .....Loop with increasing number of Legendre points.
            #
            debye0 = 1e30
            nl = 5
            xabs = math.fabs(debye-debye0)
            while (nl <= maxnl) and (xabs >= eps):
                xabscissa, weight, logstr = self.gauss_legendre(cero, y, xabscissa, weight, nl, maxloops, logstr)
                total = 0.0 
                for i in xrange(nl):
                    total = total + weight[i] * (self.fdebyeint(xabscissa[i]))
                debye = total/y/y/y
                xabs = math.fabs(debye-debye0)
                debye0 = debye
                nl = nl + 5
        return debye, logstr




    def thermalconductivity(self, thetaD, temp, ntemp, gammaD, voleqD, avmass):
        """
        Calculates thermal conductivity using Slack formulation of Liebfried-Schloemann eqn.
        """
        kboltz = 1.3807e-23
        hbar = 1.05459e-34
        dpi = 3.141592653589793238462643
        L = 5.72e+7
        tol = 1e-12
        third = 1.0/3.0
        kappaT = []
        kth = kboltz * thetaD / hbar
        kths = kth * kth
        vDcr = voleqD**third
        kappaD = ((0.849 * 3 * (4**third)) / (20 * (dpi**3) * (1.0 - (0.514 / gammaD) + (0.228 / (gammaD * gammaD))))) * kths * kboltz * vDcr * avmass / (hbar * gammaD * gammaD)
        for i in xrange(ntemp):
            if math.fabs(temp[i]) < tol:
                kappaT.append(0.0)
            else:
                kappaT.append(kappaD * thetaD / temp[i])
        return kappaD, kappaT




