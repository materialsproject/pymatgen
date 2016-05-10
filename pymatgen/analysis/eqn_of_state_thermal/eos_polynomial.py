#!/usr/bin/env python

"""
This module fits (E, V) by a polynomial, evaluates polynomials and their derivatives at given points, and finds the minimum.
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
#  This set of functions fits (E, V) by a polynomial, and finds the derivatives and minimum
# ************************************************************************************************

class eos_polynomial:
    def _init_():
        pass

    #
    #.....polynomial_fit - fits polynomials to (variable, function) data and averages
    #     them, weighted by its chi-square test probabilities. It returns
    #     the averaged polynomial coefficients in polynomialcoeffs, and the
    #     coefficients of its square in polynomialerrors.
    #
    # Adapted from original Fortran version written by M. A. Blanco et al.
    # See Computer Physics Communications 158, 57-72 (2004) and Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details
    def polynomial_fit(self, imin, xconfigvector, datatofit, ndata, mpar, mfit):
        # Data containers
        polycoeffswork = []
        weight = []
        # Promedium mechanism
        mfit = 500
        npolycoeffsfit = []
        ndatafit = []
        xdata_tofit = []
        ydata_tofit = []
        weight_tofit = []
        polynomialcoeffs = []
        polynomialerrors = []
        # Error check
        fterr = 0
        # Service variables
        limit = 4
        ndel = 3
        # Initialize the weights
        for i in xrange(ndata):
            weight.append(1.0)
        for i in xrange(mpar+1):
            polycoeffswork.append(0.0)
        for i in xrange(mfit):
            npolycoeffsfit.append(0.0)
            ndatafit.append(0.0)
        #.....Fitting loops through different polynomials
        nfit = 0
        npolynomialcoeffs = 0
        npolycoeffsmin = 2
        iinf = 0
        isup = ndata - 1
        ndatamin = max(ndata-2*ndel, npolycoeffsmin+1)
        ndataact = ndata
        rmsmin = 1.0e33
        #.....Loop eliminating pairs of outermost elements ndel times:
        while (ndataact >= ndatamin) and ((iinf < imin) and (isup > imin)):
        #.....Loop increasing the number of parameters until limit is reached:
            npolycoeffsmax = min(ndataact-limit-1, mpar)
            npolycoeffs = npolycoeffsmin
            while npolycoeffs <= npolycoeffsmax:
                nfit = nfit + 1
                rms, polycoeffswork = self.polfit(iinf, isup, xconfigvector, datatofit, weight, npolycoeffs)
                #.....Discard fits that don't give a minimum in the input bracket
                p1m1 = self.polin1(xconfigvector[imin-1], npolycoeffs, polycoeffswork)
                p1p1 = self.polin1(xconfigvector[imin+1], npolycoeffs, polycoeffswork)
                tmp = p1m1 * p1p1
                if tmp >= 0.0:
                    nfit = nfit - 1
                #.....Discard fits over the maximum declared dimensions
                elif nfit > (mfit-1):
                    nfit = mfit-1
                #.....Save fit parameters
                npolynomialcoeffs = max(npolynomialcoeffs, npolycoeffs)
                npolycoeffsfit[nfit] = npolycoeffs
                ndatafit[nfit] = ndataact
                rmsmin = min(rmsmin, rms*npolycoeffs/ndataact)
                npolycoeffs = npolycoeffs + 1
            iinf = iinf + 1
            isup = isup - 1
            ndataact = ndataact - 2
        #.....Number of fits control
        if nfit == 0:
            fterr = 1
            return fterr, npolycoeffs, polynomialcoeffs, polynomialerrors
        #.....average the polynomial coefficients (repeating the fits)
        wnorm = 0.0
        del polynomialcoeffs[:]
        del polynomialerrors[:]
        i = 0
        while i <= npolynomialcoeffs:
            polynomialcoeffs.append(0.0)
            i = i + 1
            i = 0
        while i <= 2*npolynomialcoeffs:
            polynomialerrors.append(0.0)
            i = i + 1
        ifit = 1
        while ifit <= nfit:
            ndataact = ndatafit[ifit]
            npolycoeffs = npolycoeffsfit[ifit]
            iinf = int(1 + ((eos_thermal_data.ndata-ndataact) / 2) - 1)
            isup = eos_thermal_data.ndata - iinf - 1
            rms, polycoeffswork = self.polfit(iinf, isup, eos_thermal_data.xconfigvector, eos_thermal_data.datatofit, weight, npolycoeffs)
            wtmp = rms*(npolycoeffs+1)/(rmsmin*(ndataact+1))
            wtmp = math.exp(-wtmp*wtmp)
            i = 0
            while i <= npolycoeffs:
                polynomialcoeffs[i] = polynomialcoeffs[i] + wtmp * polycoeffswork[i]
                polynomialerrors[2*i] = polynomialerrors[2*i] + wtmp * polycoeffswork[i] * polycoeffswork[i]
                j = 0
                while j < i:
                    polynomialerrors[i+j] = polynomialerrors[i+j] + 2.0 * wtmp * polycoeffswork[j] * polycoeffswork[i]
                    j = j + 1
                i = i + 1
            wnorm = wnorm + wtmp
            ifit = ifit + 1
        #.....Put the proper weight into the polynomials
        i = 0
        while i <= npolynomialcoeffs:
            polynomialcoeffs[i] = polynomialcoeffs[i] / wnorm
            polynomialerrors[i] = polynomialerrors[i] / wnorm
            i = i + 1
        while i <= 2*npolynomialcoeffs:
            polynomialerrors[i] = polynomialerrors[i] / wnorm
            i = i + 1
        #.....end of routine
        return fterr, npolynomialcoeffs, polynomialcoeffs, polynomialerrors


    #
    # .....debfitt - fitts Ln theta - Ln V data to polynomials, and
    #     averages them.
    #
    #     It returns the averaged polynomial coefficients in polynomialcoeffs.
    #
    #
    # Adapted from original Fortran version written by V. Luana and M. A. Blanco, Departamento de Quimica Fisica y Analitica, Universidad de Oviedo, 33006-Oviedo, Spain.
    # See Computer Physics Communications 158, 57-72 (2004) and Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details
    #
    def debfitt(self, volinp, tdebye, ndata, mpar, mfit, logstr):
        polycoeffswork = []
        weight = []
        x = []
        y = []
        #
        # .....promedium mechanism;
        #
        npolycoeffsfit = []
        ndatafit = []
        polynomialcoeffs = []
        polynomialerrors = []
        #
        # .....service variables:
        #
        #
        # .....limit in the dimensions of polynomials (# of points - limit)
        #
        limit = 4
        #
        # .....number of pairs of data points to delete in the fitting procedure
        #
        ndel = 3
        fterr = 0
        #
        # .....initialize the weights and the variables
        #
        for i in xrange(ndata):
            x.append(math.log(volinp[i]))
            y.append(math.log(tdebye[i]))
            weight.append(1.0)
        for i in xrange(mpar+1):
            polycoeffswork.append(0.0)
        for i in xrange(mfit+1):
            npolycoeffsfit.append(0.0)
            ndatafit.append(0.0)
            xdata_tofit = []
        ydata_tofit = []
        weight_tofit = []
        #
        # .....fitting loops through different polynomials:
        #
        nfit = 0
        npolynomialcoeffs = 0
        npolycoeffsmin = 2
        iinf = 0
        isup = ndata - 1
        ndatamin = max(ndata-2*ndel, npolycoeffsmin+1)
        ndataact = ndata
        rmsmin = 1e33
        #
        # .....loop eliminating pairs of outermost elements ndel times:
        #
        while ndataact >= ndatamin:
            #
            # .....loop increasing the number of parameters until limit is reached:
            #
            npolycoeffsmax = min(ndataact-limit-1, mpar)
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
                rms, polycoeffswork = self.polfit(iinf, isup, x, y, weight, npolycoeffs)
                #
                # .....discard fits over the maximum declared dimensions
                #
                if nfit > mfit:
                    nfit = mfit
                    logstr = logstr + "MP Eqn of State Thermal debfitt: Warning! maximum number of fits exceeded \n"
            #
            #.....save fit parameters
            #
            else:
                npolynomialcoeffs = max(npolynomialcoeffs, npolycoeffs)
                npolycoeffsfit[nfit] = npolycoeffs
                ndatafit[nfit] = ndataact
                rmsmin = min(rmsmin, rms*npolycoeffs/ndataact)
            npolycoeffs = npolycoeffs + 1
        iinf = iinf + 1
        isup = isup - 1
        ndataact = ndataact - 2
        #
        # .....number of fits control
        #
        if nfit == 0:
            logstr = logstr + "MP Eqn of State Thermal fitt: no fits to average! \n"
            fterr = 1
            return npolynomialcoeffs, polynomialcoeffs, logstr, fterr
        #
        # .....average the polynomial coefficients (repeating the fits)
        #
        wnorm = 0.0
        for i in xrange(npolynomialcoeffs+1):
            polynomialcoeffs.append(0.0)
        ifit = 1
        while ifit <= nfit:
            ndataact = ndatafit[ifit]
            npolycoeffs = npolycoeffsfit[ifit]
            iinf = int(1 + ((ndata-ndataact) / 2) - 1)
            isup = ndata - iinf - 1
            del xdata_tofit[:]
            del ydata_tofit[:]
            del weight_tofit[:]
            i = iinf
            while i < (isup + 1):
                xdata_tofit.append(x[i])
                ydata_tofit.append(y[i])
                weight_tofit.append(weight[i])
                i = i + 1
            rms, polycoeffswork = self.polfit(iinf, isup, x, y, weight, npolycoeffs)
            wtmp = rms*npolycoeffs/(rmsmin*ndataact)
            wtmp = math.exp(-wtmp*wtmp)
            for i in xrange(npolycoeffs + 1):
                polynomialcoeffs[i] = polynomialcoeffs[i] + wtmp * polycoeffswork[i]
            wnorm = wnorm + wtmp
            ifit = ifit + 1
        #
        # .....put the proper weight into the polynomials
        #
        for i in xrange(npolynomialcoeffs + 1):
            polynomialcoeffs[i] = polynomialcoeffs[i] / wnorm
        #
        # .....end of routine
        #
        return npolynomialcoeffs, polynomialcoeffs, logstr, fterr




    #
    #.....minbrack - given a table of values, and an initial index imin,
    #     the routine searches for the minimum on the list by a downhill
    #     algorithm, assuming the index as an abscissa.
    #
    # Adapted from original Fortran version written by M. A. Blanco et al.
    # See Computer Physics Communications 158, 57-72 (2004) and Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details
    #
    def minbrack(self, imin, data_to_fit, ndata, logstr):
        ierr = 0
        tol = 1e-12
        #
        #.....imin outside limits
        #
        if imin <= 1 or imin >= ndata:
            logstr = logstr + "MP Eqn of State Thermal minbrack: trial point outside limits \n"
            logstr = logstr + "MP Eqn of State Thermal minbrack: imin = " + str(imin) + " \n"
            ierr = 1
            return ierr, imin, logstr
        #
        #.....decide the direction of the search
        #
        istep = 0
        if datatofit[imin] < datatofit[imin+1] and datatofit[imin] < datatofit[imin-1]:
            return ierr, imin, logstr
        elif math.fabs(datatofit[imin] - datatofit[imin+1]) < tol and math.fabs(datatofit[imin] - datatofit[imin-1]) < tol:
            logstr = logstr + "MP Eqn of State Thermal minbrack: flat function, impossible to search \n"
            ierr = 2
            return ierr, imin, logstr
        elif datatofit[imin] >= datatofit[imin+1] and datatofit[imin] <= datatofit[imin-1]:
            istep = 1
        elif datatofit[imin] <= datatofit[imin+1] and datatofit[imin] >= datatofit[imin-1]:
            istep = -1
        else:
            logstr = logstr + "MP Eqn of State Thermal minbrack: trial point is a maximum \n"
            ierr = 3
            logstr = logstr + "MP Eqn of State Thermal minbrack: imin = " + str(imin)
            logstr = logstr + "MP Eqn of State Thermal minbrack: func[imin - 1] = " + str(datatofit[imin-1]) + " \n"
            logstr = logstr + "MP Eqn of State Thermal minbrack: func[imin] = " + str(datatofit[imin]) + " \n"
            logstr = logstr + "MP Eqn of State Thermal minbrack: func[imin + 1] = " + str(datatofit[imin+1]) + " \n"
            return ierr, imin, logstr
        imin = imin + istep
        #
        #.....search for a minimum pattern
        #
        while imin > 1 and imin < ndata:
            if datatofit[imin] > datatofit[imin+istep]:
                imin = imin + istep
            #
            #.....minimum pattern found
            #
            else:
                return ierr, imin, logstr
        #
        #.....no minimum pattern
        #
        logstr = logstr + "MP Eqn of State Thermal minbrack: monotonic function, there's no minimum \n"
        ierr = 4
        return ierr, imin, logstr




    #
    #.....polminbrack - given a polynomial and a set of x values, evaluates the polynomial at each x-point,
    #     then searches for the global minimum of the polynomial by checking all of the calculated values
    #     all points on the list
    #
    def polminbrack(self, imin, npol, ndata, x, polynomialcoeffs):
        y = []
        for istep in xrange(ndata):
            plnv = polin0(x[istep], npol, polynomialcoeffs)
            y.append(plnv)
        imin = 0
        funcmin = y[imin]
        istep = 1
        while istep < ndata:
            if y[istep] < funcmin:
                imin = istep
            istep = istep + 1
        #
        #.....imin outside limits
        #
        if imin < 1 or imin > (ndata-1):
            return 1, imin
        else:
            return 0, imin




    #
    #.....polmin - gets the minimum of a polymomial, using Newton-Raphson
    #     method to zero its first derivative.
    #
    #     The minimum must be in the interval [xmin,xmax], and whenever
    #     Newton's method doesn't converge, a bisection step will be given
    #     for this interval
    #
    # Adapted from original Fortran version written by M. A. Blanco et al.
    # See Computer Physics Communications 158, 57-72 (2004) and Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details
    #
    def polmin(self, xinitial, xmin, xmax, npolycoeffs, polycoeffswork, pmerr, logstr):
        x = xinitial
        a = xmin
        b = xmax
        f_x = polin1(x, npolycoeffs, polycoeffswork)
        f_a = polin1(a, npolycoeffs, polycoeffswork)
        f_b = polin1(b, npolycoeffs, polycoeffswork)
        dx = b - a
        tol = dx*1e-6
        pmerr = 0
        if a > b or (not self.isbr(x - a, x - b)) or (not self.isbrack(f_a, f_b)):
            logstr = logstr + "MP Eqn of State Thermal polmin: wrong input data \n"
            logstr = logstr + "MP Eqn of State Thermal polmin: x = " + str(x) + ', a = ' + str(a) + ', b = ' + str(b) + ' \n'
            logstr = logstr + "MP Eqn of State Thermal polmin: f_x= " + str(f_x) + ' f_a= ' + str(f_a) + ', f_b= ' + str(f_b) + ' \n'
            logstr = logstr + "MP Eqn of State Thermal: npolycoeffs = " + str(npolycoeffs) + " \n"
            logstr = logstr + "MP Eqn of State Thermal: polycoeffswork = " + str(polycoeffswork) + " \n"
            # Checks error type and sets pmerr appropriately to tell gibbsrun function which correction steps to attempt
            # Errors are checked in reverse order of severity (i.e. most severe error is checked last)
            # This sets pmerr signal to most severe error warning so this can be corrected first (if possible)
            if not self.isbrack(f_a, f_b):
                if f_a > 0.0 and f_b > 0.0:
                    logstr = logstr + "MP Eqn of State Thermal polmin: f_a and f_b are both greater than 0 \n"
                    if f_a > f_b:
                        pmerr = 1
                    else:
                        pmerr = 2
                elif f_a < 0.0 and f_b < 0.0:
                    logstr = logstr + "MP Eqn of State Thermal polmin: f_a and f_b are both less than 0 \n"
                    if f_a > f_b:
                        pmerr = 2
                    else:
                        pmerr = 1
            if not self.isbr(x - a, x - b):
                if a > x:
                    logstr = logstr + "MP EOS_THERMAL polmin: a > x \n"
                    pmerr = 3
                elif b < x:
                    logstr = logstr + "MP Eqn of State Thermal polmin: b < x \n"
                    pmerr = 4
            if a > b:
                logstr = logstr + "MP Eqn of State Thermal polmin: a > b \n"
                pmerr = 5
                return pmerr, x, logstr
        while (math.fabs(b - a) > tol) and (math.fabs(dx) > tol):
            if self.isbrack(f_a, f_x):
                b = x
                f_b = fx
            else:
                a = x
                f_a = f_x
                dx = polin2(x, npolycoeffs, polycoeffswork)
                f_x = polin1(x, npolycoeffs, polycoeffswork)
            x_0 = x
            if not self.isbr(dx * (x_0 - a) - f_x, dx * (x_0 - b) - f_x):
                x = 0.5 * (b + a)
                dx = x - x_0
            else:
                dx = - f_x / dx
                x = x_0 + dx
        return pmerr, x, logstr


    def isbr(self, fa, fb):
        if fa >= 0.0 and fb <= 0.0:
            return True
        else:
            return False


    def isbrack(self, fa, fb):
        if self.isbr(fa, fb) or self.isbr(fb, fa):
            return True
        else:
            return False




    #
    # ..... polfit - fit data by a polynomial of order npolycoeffs
    #
    #      Calculates the coefficients of a polynomial of order npolycoeffs to fit a x-y data set
    #      Generates linear matrix system c from data set to be fitted
    #      Calls Gaussian elimination function to solve this system and obtain coefficients
    #      Calculates RMS deviation for the fit between the polynomial and the input data
    #
    def polfit(self, iinf, isup, x, y, w, npolycoeffs):
        c = []
        # .....Calculate the norm of the weights
        wnorm = 0.0
        k = iinf
        while k <= isup:
            wnorm = wnorm + w[k]
            k = k + 1
            wnorm = 1.0 / wnorm
        for i in xrange(npolycoeffs+1):
            c.append([])
            for j in xrange(npolycoeffs+2):
                c[i].append(0.0)
        # .....Construct the linear system matrix c:
        j = 0
        while j <= npolycoeffs:
            c[j][npolycoeffs+1] = 0.0
            if j > 0:
                k = iinf
                while k <= isup:
                    c[j][npolycoeffs+1] = c[j][npolycoeffs+1] + w[k] * y[k] * (x[k]**j)
                    k = k + 1
            else:
                k = iinf
                while k <= isup:
                    c[j][npolycoeffs+1] = c[j][npolycoeffs+1] + w[k] * y[k]
                    k = k + 1
            c[j][npolycoeffs+1] = wnorm * c[j][npolycoeffs+1]
            i = j
            while i <= npolycoeffs:
                c[i][j] = 0.0
                ij = i + j
                if ij > 0:
                    k = iinf
                    while k <= isup:
                        c[i][j] = c[i][j] + w[k] * (x[k]**ij)
                        k = k + 1
                        c[i][j] = wnorm * c[i][j]
                else:
                    c[i][j] = 1.0
                c[j][i] = c[i][j]
                i = i + 1
            j = j + 1
        # .....Solve the linear system for the best A()'s:
        amatrix = []
        for i in xrange(npolycoeffs+1):
            amatrix.append([])
            for j in xrange(npolycoeffs+1):
                amatrix[i].append(0.0)
        for i in xrange(npolycoeffs+1):
            for j in xrange(npolycoeffs+1):
                amatrix[i][j] = c[i][j]
        bmatrix = []
        for i in xrange(npolycoeffs+1):
            bmatrix.append(0.0)
        for i in xrange(npolycoeffs+1):
            bmatrix[i] = c[i][npolycoeffs+1]
        polycoeffswork = np.linalg.solve(amatrix, bmatrix)
        list(polycoeffswork)
        # .....Compute the rms deviation:
        s2 = 0.0
        k = iinf
        while k <= isup:
            s2 = s2 + w[k] * ((y[k] - polin0(x[k], npolycoeffs, polycoeffswork))**2)
            k = k + 1
        rms = math.sqrt(s2 * wnorm)
        return rms, polycoeffswork



    # **************************************************************************************
    #  This set of functions evaluate polynomials and their derivatives
    # **************************************************************************************

    #
    # .....polin0 - Horner's evaluation of a polynomial.
    #
    #     The polynomial is given by:
    #
    #       y(x) = SUM(i=0,npolycoeffs) polycoeffswork(i) * x**i
    #
    #     It is assumed that npolycoeffs>=0.
    #
    def polin0(self, x, npolycoeffs, polycoeffswork):
        y = polycoeffswork[npolycoeffs]
        i = npolycoeffs
        while i > 0:
            y = y * x + polycoeffswork[i-1]
            i = i - 1
        return y

    #
    # .....polin1 - Horner's evaluation of the first derivative of a
    #     polynomial.
    #
    #     It is assumed that npolycoeffs>=1.
    #
    def polin1(self, x, npolycoeffs, polycoeffswork):
        y = polycoeffswork[npolycoeffs] * (npolycoeffs)
        i = npolycoeffs
        while i > 1:
            y = y * x + polycoeffswork[i-1] * (i-1)
            i = i - 1
        return y


    #
    # .....polin2 - Horner's evaluation of the second derivative of a
    #     polynomial.
    #
    #     It is assumed that npolycoeffs>=2.
    #
    def polin2(self, x, npolycoeffs, polycoeffswork):
        y = polycoeffswork[npolycoeffs] * (npolycoeffs) * (npolycoeffs-1)
        i = npolycoeffs
        while i > 2:
            y = y * x + polycoeffswork[i-1] * (i-1) * (i-2)
            i = i - 1
        return y


    #
    # .....polin3 - Horner's evaluation of the third derivative of a
    #    polynomial.
    #
    #     It is assumed that npolycoeffs>=3.
    #
    def polin3(self, x, npolycoeffs, polycoeffswork):
        if npolycoeffs < 3:
            return 0.0
        y = polycoeffswork[npolycoeffs] * (npolycoeffs) * (npolycoeffs-1) * (npolycoeffs-2)
        i = npolycoeffs
        while i > 3:
            y = y * x + polycoeffswork[i-1] * (i-1) * (i-2) * (i-3)
            i = i - 1
        return y


    #
    # .....polin4 - Horner's evaluation of the fourth derivative of a
    #     polynomial.
    #
    #     It is assumed that npolycoeffs>=4.
    #
    def polin4(self, x, npolycoeffs, polycoeffswork):
        if npolycoeffs < 4:
            return 0.0
        y = polycoeffswork[npolycoeffs] * (npolycoeffs) * (npolycoeffs-1) * (npolycoeffs-2) * (npolycoeffs-3)
        i = npolycoeffs
        while i > 4:
            y = y * x + polycoeffswork[i-1] * (i-1) * (i-2) * (i-3) * (i-4)
            i = i - 1
        return y


