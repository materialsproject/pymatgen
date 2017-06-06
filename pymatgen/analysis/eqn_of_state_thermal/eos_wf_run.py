#!/usr/bin/env python

"""
This module is a wrapper for the MP-Eqn of State thermal properties module.
"""

from __future__ import division
import warnings
import sys
import subprocess
import unittest
import pymatgen
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.analysis.eqn_of_state_thermal.eos_ev_check import eos_ev_check
from pymatgen.analysis.eqn_of_state_thermal.eos_polynomial import eos_polynomial
from pymatgen.analysis.eqn_of_state_thermal.eos_thermal_functions import eos_thermal_functions
import numpy as np
import os
import warnings
from numpy import matrix
from numpy import linalg
import math

__author__ = "Cormac Toher"
__copyright__ = "Copyright 2014, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Cormac Toher"
__email__ = "cormac.toher@duke.edu"
__date__ = "January 14, 2015"


class eos_thermal_properties:
    def _init_(self):
        self.vol_inp = []
        self.energ_inp = []
        self.tdebye = []
        self.pressure = []
        self.temperature = []
        self.xconfigvector = []
        self.datatofit = []
        self.bmcoeffs_static = []
        self.mfit = 500
        self.mpar = 20
        self.maxloops = 250
        self.birchfitorder_iG = 2
        self.opt_g = False
        self.iopt_g = 0
        self.fittype = 0
        self.fterr = 0
        self.third = 1.0 / 3.0
        self.hy2kjmol = 2625.5 
        self.hart2ev = 27.211383
        self.ev2hartree = 0.036749309
        self.au2gpa = 29421.4  
        self.pckbau = 3.166830e-6 
        self.pi = 3.14159265358979323846
        self.kj2unit = 0.1202707236

    def eos_thermal_run(self, initstruct, volume_values, energy_values, poissonratio = 0.25, ieos = 0, idebye = 0, npres = 11, spres = 2.0, ntemp = 201, stemp = 10.0, energy_inf = 0.0, fit_type = 0):
        nstructs = len(volume_values)
        self.vol_inp = volume_values
        self.energ_inp = energy_values
        if len(self.vol_inp) != nstructs:
            raise ValueError("Wrong number of volumes")
            return
        if len(self.energ_inp) != nstructs:
            raise ValueError("Wrong number of energies")
            return
        for i in xrange(nstructs):
            print("Volume (in Ang^3) = ", self.vol_inp[i], "Energy (in eV) = ", self.energ_inp[i])
        self.energ_inf = energy_inf
        self.tdebye = []
        self.natoms = len(initstruct.sites)
        self.title = ""
        for i in xrange(self.natoms):
            self.title = self.title + initstruct.species[i].symbol
        self.mfit = 500
        self.mpar = 20
        self.maxloops = 250
        self.birchfitorder_iG = 2
        self.opt_g = False
        self.iopt_g = 0
        self.fittype = 0
        self.fterr = 0
        self.third = 1.0 / 3.0
        self.hy2kjmol = 2625.5 
        self.hart2ev = 27.211383
        self.ev2hartree = 0.036749309
        self.au2gpa = 29421.4  
        self.pckbau = 3.166830e-6 
        self.pi = 3.14159265358979323846
        self.kj2unit = 0.1202707236
        self.pcamu = 1.6605402e-24
        self.pcme = 9.1093897e-28
        self.amu2au = self.pcamu/self.pcme
        self.amu2kg = 1.660538921e-27
        self.ang32bohr3 = 6.74833303710415
        self.npressure = npres
        self.ntemperature = ntemp
        self.fittype = fit_type
        self.pressure = []
        self.temperature = []
        for i in xrange(npres):
            idoub = i
            self.pressure.append(idoub * spres)
        for i in xrange(ntemp):
            idoub = i
            self.temperature.append(idoub * stemp)
        self.savpres = True
        self.outfname = self.title + '.out'
        self.inpfname = self.title + '.inp'
        self.ieos = int(ieos)
        self.idebye = int(idebye)
        self.poissonratio = float(poissonratio)
        self.ninp = int(nstructs)
        self.ndata = self.ninp
        self.outstr = ""
        self.outstr_thermal = ""
        self.bmcoeffs_static = []
        for i in xrange(nstructs):
            self.vol_inp[i] = self.vol_inp[i] * self.ang32bohr3
            self.energ_inp[i] = self.energ_inp[i] * self.ev2hartree
            print("Volume (in Bohr^3) = ", self.vol_inp[i], "Energy (in Hartree) = ", self.energ_inp[i])
        self.cellmassamu, self.cellmassval = self.cellmass(initstruct)
        eos_ev_check_inst = eos_ev_check()
        minpos = eos_ev_check_inst.checkminpos(self.vol_inp, self.energ_inp, self.ninp)
        concav = eos_ev_check_inst.checkevconcav(self.vol_inp, self.energ_inp, self.ninp)
        if (minpos != 0):
            raise ValueError("Minimum is not contained in calculated (E, V) data")
            return
        if (concav != 0):
            raise ValueError("E(V) data is not concave")
            return
# Run GIBBS algorithm to calculate thermal properties
        self.eos_thermal_gibbs()
        print(self.logstr)
        print(self.outfname)
        fo=open(self.outfname, "w")  
        fo.write(self.outstr)
        fot=open("eos_thermal_properties.out", "w")
        fot.write(self.outstr_thermal)
        fotd=open("Debye_temperature", "w")
        fotd.write("# Debye temperature as a function of temperature \n")
        for j in xrange(self.ntemperature):
            jstr = str(self.temperature[j]) + "\t" + str(self.tdeb0p[j]) + "\n"
            fotd.write(jstr)
        tdmin = self.tdeb0p[0]
        tdmax = self.tdeb0p[0]
        j = 1
        while (j < self.ntemperature):
            if (self.tdeb0p[j] > tdmax):
                tdmax = self.tdeb0p[j]
                jtdmax = j
            if (self.tdeb0p[j] < tdmin):
                tdmin = self.tdeb0p[j]
                jtdmin = j
            j = j + 1
        print("tdmin = ", tdmin)
        print("tdmax = ", tdmax)
        nkb = 1.0 * self.natoms
        eos_thermal_functions_inst = eos_thermal_functions()
        tdbestfit, self.logstr = eos_thermal_functions_inst.cvdebfit(self.cvu0p, tdmin, tdmax, nkb, self.temperature, self.ntemperature, self.maxloops, self.logstr)
        print("Best fit Debye temperature = ", tdbestfit)
        difftdmin = math.fabs(tdbestfit - self.tdeb0p[0])
        jtdbest = 0
        j = 1
        while (j < self.ntemperature):
            difftd = math.fabs(tdbestfit - self.tdeb0p[j])
            if (difftd < difftdmin):
                jtdbest = j
                difftdmin = difftd
            j = j + 1
        jt300 = 0
        j = 1
        difft300min = math.fabs(self.temperature[j] - 300.0) 
        while (j < self.ntemperature):
            difft300 = math.fabs(self.temperature[j] - 300.0)
            if (difft300 < difft300min):
                jt300 = j
                difft300min = difft300
            j = j + 1    
        print("jtdbest = ", jtdbest)
        print("jt300 = ", jt300)
        voleqD = (self.volminsav[jtdbest][0] / self.ang32bohr3) * 1e-30
        avmasskg = (self.cellmassamu * self.amu2kg) / self.natoms
        thetaa = tdbestfit * (self.natoms**(-self.third))
        print("thetaa = ", thetaa)
        print("gamma = ", self.ga0p[jtdbest])
        kappaD, kappaT = eos_thermal_functions_inst.thermalconductivity(thetaa, self.temperature, self.ntemperature, self.ga0p[jtdbest], voleqD, avmasskg)
        print("Thermal conductivity at Debye temperature = ", kappaD)
        print("Thermal conductivity at ", self.temperature[jt300], "K = ", kappaT[jt300])
        print("eos_thermal_run complete")
        eos_results_dict = {}
        eos_results_dict["temperature"] = self.temperature
        eos_results_dict["pressure"] = self.pressure
        eos_results_dict["thermal_conductivity"] = kappaT
        eos_results_dict["Best_fit_temperature"] = jtdbest
        eos_results_dict["300K_point"] = jt300
        eos_results_dict["Debye_temperature"] = self.tdeb0p
        eos_results_dict["Gruneisen_parameter"] = self.ga0p
        eos_results_dict["Heat_capacity_Cv"] = self.cvu0p
        eos_results_dict["Heat_capacity_Cp"] = self.cp0p
        eos_results_dict["Volume"] = self.vol0p
        eos_results_dict["Bulk_modulus"] = self.b0p
        eos_results_dict["BM_coeffs"] = self.bmcoeffs_static
        return eos_results_dict

    def cellmass(self, initstruct):
        cellmassamu = 0.0
        for i in xrange(self.natoms):
            cellmassamu = cellmassamu + initstruct.species[i].atomic_mass
        return cellmassamu, cellmassamu * self.amu2au


    def eos_thermal_gibbs(self):
        """
        .....gibbs - Debye-Grunseisen model for thermal properties

        Adapted from original Fortran version written by M. A. Blanco et al.
        See Computer Physics Communications 158, 57-72 (2004) and
        Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details

        This automated version described in Phys. Rev. B 90, 174107 (2014)
        """
        #
        # .....Compute the scaling factor function f(poisson) in the Debye model.
        # (See for example CMS 98, 34 (2015), Poirier - Intro. to the
        #  Physics of the Earth's Interior, and Miguel A. Blanco, PhD Thesis)
        #
        #        fx = 2*(1 + self.poisson)/3.0/(1 - 2*self.poisson)
        #        gx = (1 + self.poisson)/3.0/(1 - self.poisson)
        #        hx = 2.0*math.sqrt(fx**3)+math.sqrt(gx**3)
        #        self.poratio = math.exp(-math.log(hx/3)/3)
        poissonratiofunc1 = (2.0/3.0) * ((1.0 + self.poissonratio) / (1.0 - 2.0 * self.poissonratio))
        poissonratiofunc2 = (1.0/3.0) * ((1.0 + self.poissonratio) / (1.0 - self.poissonratio))
        poissonratiofunc3 = 2.0 * math.sqrt(poissonratiofunc1**3) + math.sqrt(poissonratiofunc2**3)
        self.poissonratiofunction = (3.0 / poissonratiofunc3)**(1.0 / 3.0)
        # String for holding data to be written to main output file
        # Write appropriate header for GIBBS output data for all p, T values
        if self.ieos >= 0:
            self.outstr = ('MP Eqn of State Thermal - (P,T) '+
                                   'thermodynamics of crystals from (E,V) data \n')
            self.outstr = self.outstr + '\n'
            self.outstr = self.outstr + self.title + '\n'
            self.outstr = (self.outstr + 'Number of data points: ' +
                                   str(self.ndata) + ' \n')
            self.outstr = self.outstr + ' \n'
        else:
            self.outstr = ('MP Eqn of State Thermal - (P,T) ' +
                                   'thermodynamics of crystals from (E,V) data \n')
            self.outstr = self.outstr + ' \n'
            self.outstr = self.outstr + str(self.title) + ' \n'
            self.outstr = self.outstr + ' \n'
            self.outstr = (self.outstr + '# P(GPa) \t  T(K) \t V(bohr^3)' +
                                   '\t  G(kJ/mol) \t Bt(GPa)' + '\n')
        # String for holding log information on calculation status, error messages,
        # warnings and information about automatic parameter adjustments
        self.logstr = ('MP Eqn of State Thermal: Log file: ' +
                               'status information, errors and warnings \n')
        for i in xrange(self.ninp):
            self.energ_inp[i] = self.energ_inp[i] - self.energ_inf
        eos_ev_check_inst = eos_ev_check()
        if self.idebye == 1:
            self.vol_inp, self.energ_inp, self.tdebye = eos_ev_check_inst.evtsort(self.vol_inp, self.energ_inp, self.tdebye, self.ninp)
        else:
            self.vol_inp, self.energ_inp = eos_ev_check_inst.evsort(self.vol_inp, self.energ_inp, self.ninp)
        #
        # .....Find global minimum of energy data points
        #
        etref = self.energ_inp[0]
        itmin = 0
        for i in xrange(self.ndata):
            if self.energ_inp[i] < etref:
                etref = self.energ_inp[i]
                itmin = i
        #
        # .....Check that all points show concave patterns
        # .....First check where the concave pattern begins
        #
        j = 1
        while ((self.energ_inp[j]-self.energ_inp[j-1])/(self.vol_inp[j]-self.vol_inp[j-1]) >= (self.energ_inp[j]-self.energ_inp[j+1])/(self.vol_inp[j]-self.vol_inp[j+1])) and j < self.ndata-2:
            j = j + 1
        # If only the last three points are concave, then MP Eqn of State Thermal will not be able to use them to fit a polynomial
        # MP Eqn of State Thermal exits giving an error message
        if j == self.ndata-2:
            raise ValueError("MP Eqn of State Thermal: All points show convex patterns")
#            self.logstr = self.logstr + "MP Eqn of State Thermal: All points show convex patterns \n"
            self.grerr = 2
            return
        self.energ_inp[0] = self.energ_inp[j-1]
        self.vol_inp[0] = self.vol_inp[j-1]
        self.energ_inp[1] = self.energ_inp[j]
        self.vol_inp[1] = self.vol_inp[j]
        self.energ_inp[2] = self.energ_inp[j+1]
        self.vol_inp[2] = self.vol_inp[j+1]
        if self.idebye == 1:
            self.tdebye[0] = self.tdebye[j-1]
            self.tdebye[1] = self.tdebye[j]
            self.tdebye[2] = self.tdebye[j+1]
        j = j + 1
        jtmax = 2
        #
        # .....j marks the last accepted point, i the new trial point
        #
        i = j + 1
        while i < (self.ndata-1):
            if self.fittype == 0 or self.fittype == 3:
                if (self.energ_inp[j]-self.energ_inp[j-1])/(self.vol_inp[j]-self.vol_inp[j-1]) < (self.energ_inp[j]-self.energ_inp[i])/(self.vol_inp[j]-self.vol_inp[i]):
                    j = j + 1
                    jtmax = i
                    self.energ_inp[j] = self.energ_inp[i]
                    self.vol_inp[j] = self.vol_inp[i]
                    if self.idebye == 1:
                        self.tdebye[j] = self.tdebye[i]
            i = i + 1
        self.ndata = j + 1
        #
        # .....search for the minimum of the accepted data
        #
        volref = self.vol_inp[0]
        eref = self.energ_inp[0]
        imin = 0
        i = 0
        while i < self.ndata:
            if self.energ_inp[i] < eref:
                volref = self.vol_inp[i]
                eref = self.energ_inp[i]
                imin = i
            i = i + 1
        self.xconfigvector = []
        for i in xrange(self.ndata):
            self.xconfigvector.append((self.vol_inp[i]/volref)**self.third)
        # If the lowest energy corresponds to the smallest or largest volume, then the minimum lies outside of the accepted input data
        # MP EOS_THERMAL exits giving an error message
        if imin == 0 or imin == (self.ndata - 1):
            self.logstr = self.logstr + "MP Eqn of State Thermal: static minimum outside input data \n"
            self.logstr = self.logstr + "MP Eqn of State Thermal: structure with lowest energy is " + str(imin) + " \n"
            self.logstr = self.logstr + "MP Eqn of State Thermal: energy of this structure = " + str(eref) + " \n"
            self.logstr = self.logstr + "MP Eqn of State Thermal: total number of structures = " + str(self.ndata) + " \n"
            self.grerr = 2
            raise ValueError("MP Eqn of State Thermal: static minimum outside input data")
            return
        # 
        # .....Obtain the polynomial fit of E(static) vs. x
        #     x = V^(1/3)
        #     epol contains coefficients of fitted polynomial
        #
        eos_polynomial_inst = eos_polynomial()
        self.datatofit = []
        for i in xrange(self.ndata):
            self.datatofit.append(self.energ_inp[i])
        epol = []
        eerr = []
        self.fterr, nepol, epol, eerr = eos_polynomial_inst.polynomial_fit(imin, self.xconfigvector, self.datatofit, self.ndata, self.mpar, self.mfit)
        if self.fterr != 0:
            if self.fterr == 2:
#                self.logstr = self.logstr + "MP Eqn of State Thermal: Problem inverting matrix to fit polynomial \n"
                self.grerr = 4
                raise ValueError("MP Eqn of State Thermal: Problem inverting matrix to fit polynomial")
            else:
#                self.logstr = self.logstr + "MP Eqn of State Thermal: No polynomial fit found with a minimum within bounds of input data \n"
                self.grerr = 2
                raise ValueError("MP Eqn of State Thermal: No polynomial fit found with a minimum within bounds of input data")
            return
        if nepol == 2:
            nepol = nepol + 1
            epol.append(0.0)
        #
        # Find minimum of polynomial epol to find minimum of (E, V) curve
        # First bracket minimum of (E, V) data
        itry = imin
        self.ierr, itry, self.logstr = eos_polynomial_inst.minbrack(itry, self.datatofit, self.ndata, self.logstr)
        if self.ierr != 0:
#            self.logstr = self.logstr + "MP Eqn of State Thermal: Cannot find minimum of (E, V) data \n"
            self.grerr = 2
            raise ValueError("MP Eqn of State Thermal: Cannot find minimum of (E, V) data")
            return
        else:
            self.logstr = self.logstr + "MP Eqn of State Thermal: Minimum of (E, V) data is at point imin = " + str(itry) + " \n"
        #
        # Find minimum of epol, using itry as initial guess
        self.pmerr = 0
        self.pmerr, xmin, self.logstr = eos_polynomial_inst.polmin(self.xconfigvector[itry], self.xconfigvector[max(itry-2, 0)], self.xconfigvector[min(itry+2, self.ndata-1)], nepol, epol, self.pmerr, self.logstr)
        # If polmin returns an error, then MP Eqn of State Thermal tries shifting the trial point to try to correct the error
        # If this does not correct the error, then MP Eqn of State Thermal exits giving an error message
        if self.pmerr != 0:
            if self.pmerr == 1:
                while (self.pmerr == 1) and ((itry + 2) < (self.ndata-1)):
                    itry = itry + 1
                    self.logstr = self.logstr + "MP Eqn of State Thermal: Resetting itry to itry = " + str(itry) + " \n"
                    self.pmerr, xmin, self.logstr = eos_polynomial_inst.polmin(self.xconfigvector[itry], self.xconfigvector[max(itry-2, 0)], self.xconfigvector[min(itry+2, self.ndata-1)], nepol, epol, self.pmerr, self.logstr)
                # If error indicator has changed from 1 to 2, then the bracket has shifted from one side of the minimum to the other
                # Need to expand the size of the bracket to incorporate the minimum
                if self.pmerr == 2:
                    self.pmerr, xmin, self.logstr = eos_polynomial_inst.polmin(self.xconfigvector[itry], self.xconfigvector[max(itry-4, 0)], self.xconfigvector[min(itry+4, self.ndata-1)], nepol, epol, xmin, self.pmerr, self.logstr)
                if self.pmerr == 0:
                    self.logstr = self.logstr + "MP Eqn of State Thermal: Minimum of (E, V) data is at point xmin = " + str(xmin) + " Bohr \n"
                else:
                    self.logstr = self.logstr + "MP Eqn of State Thermal: Cannot find minimum of (E, V) data \n"
                    self.grerr = 2
                    return
            elif self.pmerr == 2:
                while (self.pmerr == 2) and ((itry - 2) >= 0):
                    itry = itry - 1
                    self.logstr = self.logstr + "MP Eqn of State Thermal: Resetting itry to itry = " + str(itry) + " \n"
                    self.pmerr, xmin, self.logstr = eos_polynomial_inst.polmin(self.xconfigvector[itry], self.xconfigvector[max(itry-2, 0)], self[min(itry+2, self.ndata-1)], nepol, epol, self.pmerr, self.logstr)
                # If error indicator has changed from 2 to 1, then the bracket has shifted from one side of the minimum to the other
                # Need to expand the size of the bracket to incorporate the minimum
                if self.pmerr == 1:
                    self.pmerr, xmin, self.logstr = eos_polynomial_inst.polmin(self.xconfigvector[itry], self.xconfigvector[max(itry-4, 0)], self.xconfigvector[min(itry+4, self.ndata-1)], nepol, epol, self.pmerr, self.logstr)
                if pmerr == 0:
                    self.logstr = self.logstr + "MP Eqn of State Thermal: Minimum of (E, V) data is at point xmin = " + str(xmin) + " Bohr \n"
                else:
#                    self.logstr = self.logstr + "MP Eqn of State Thermal: Cannot find minimum of (E, V) data \n"
                    self.grerr = 2
                    raise ValueError("MP Eqn of State Thermal: Cannot find minimum of (E, V) data")
                    return
            else:
#                self.logstr = self.logstr + "MP Eqn of State Thermal: Cannot find minimum of (E, V) data \n"
                self.grerr = 2
                raise ValueError("MP Eqn of State Thermal: Cannot find minimum of (E, V) data")
                return
        Vmin = xmin**3
        # Evaluate polynomial at minimum point xmin to get Emin
        Emin = eos_polynomial_inst.polin0(xmin, nepol, epol)
        # Write out Emin in different units
        self.logstr = self.logstr + "MP Eqn of State Thermal: Minimum of (E, V) data is Emin = " + str(Emin) + " Hartree/cell \n"
        self.logstr = self.logstr + "MP Eqn of State Thermal: Minimum of (E, V) data is Emin = " + str(Emin * self.hy2kjmol) + " kJ/mol \n"
        self.logstr = self.logstr + "MP Eqn of State Thermal: Minimum of (E, V) data is Emin = " + str(Emin * self.hart2ev) + " eV/cell \n"
        #
        # .....Minimize G(static) - V; G = Gibbs free energy
        #
        gpol = []
        for i in xrange(nepol+1):
            gpol.append(epol[i])
        ngpol = nepol
        itry = imin
        self.voleqmin = []
        g = []
        self.bulkmod = []
        rerr = []
        for k in xrange(self.npressure):
            #
            # .....bracket the minimum with the numerical function
            #
            for i in xrange(self.ndata):
                self.datatofit[i] = self.energ_inp[i] + self.vol_inp[i]*self.pressure[k]/self.au2gpa
            self.ierr, itry, self.logstr = eos_polynomial_inst.minbrack(itry, self.datatofit, self.ndata, self.logstr)
            #
            # .....obtain the minimum of the fitted function
            #
            # Adds pV term onto coefficient of x^3; x^3 ~ V
            gpol[3] = epol[3] + self.pressure[k]/self.au2gpa * volref
            self.pmerr, xmin, self.logstr = eos_polynomial_inst.polmin(self.xconfigvector[itry], self.xconfigvector[max(itry-2, 0)], self.xconfigvector[min(itry+2, self.ndata-1)], ngpol, gpol, self.pmerr, self.logstr)
            # Evaluates polynomial for (G, V) at minimum to get equilibrium values for V, G, and B
            self.voleqmin.append((xmin**3.0) * volref)
            plnv = eos_polynomial_inst.polin0(xmin, ngpol, gpol)
            g.append(plnv * self.hy2kjmol)
            plnv = eos_polynomial_inst.polin2(xmin, ngpol, gpol)
            self.bulkmod.append(plnv * self.au2gpa * xmin*xmin/(9.0*self.voleqmin[k]))
            eact = eos_polynomial_inst.polin0(xmin, nepol, epol)
            twonepol = 2*nepol
            plnv = eos_polynomial_inst.polin0(xmin, twonepol, eerr)
            aerr = math.sqrt(math.fabs(plnv - eact*eact))
            rerr.append(aerr / max(math.fabs(eact), math.fabs(eact)+aerr/2.0))
        #
        # .....Write static EOS results to stringstream
        #
        vol0pres = self.voleqmin[0]
        gfe0pres = g[0]
        binp_bcnt = self.bulkmod[0]
        if self.ieos >= 0:
            self.outstr = self.outstr + 'Static EOS calculation - Numerical results \n'
            self.outstr = self.outstr + 'Vmin(static;  P=0)    = ' + str(vol0pres) + ' bohr^3 \n'
            self.outstr = self.outstr + 'Gmin(static;  P=0)    = ' + str(gfe0pres) + ' kJ/mol \n'
            self.outstr = self.outstr + '\n'
            self.outstr = self.outstr + 'NUMERICAL EQUILIBRIUM PROPERTIES \n'
            self.outstr = self.outstr + '================================'
            self.outstr = self.outstr + '\n'
            self.outstr = self.outstr + '  P(GPa) \t G(kJ/mol) \t V(bohr^3)  \t   V/V0 \t B(GPa) \t rel.err. \n'
            self.outstr = self.outstr + ' ------------------------------------------------------------------------------------------- \n'
            for k in xrange(self.npressure):
                self.outstr = self.outstr + '  ' + str(self.pressure[k]).rjust(6) + '\t' + str(g[k]).rjust(10)[:10] + '\t' + str(self.voleqmin[k]).rjust(10)[:10] + '\t' + str(self.voleqmin[k]/vol0pres).rjust(7)[:7] + '\t     ' + str(self.bulkmod[k]).rjust(10)[:10] + '\t     ' + str(round(rerr[k], 10)).rjust(12)[:12] + '\n'
        elif self.idebye == -1:
            for k in xrange(self.npressure):
                self.outstr = self.outstr + '  ' + str(self.pressure[k]).rjust(6) + '\t' + str(0.0).rjust(10)[:10] + '\t' + str(self.voleqmin[k]).rjust(10)[:10] + '\t' + str(g[k]).rjust(10)[:10] + '\t     ' + str(self.bulkmod[k]).rjust(10)[:10] + '\n'
        else:
            for k in xrange(self.npressure):
                self.outstr = self.outstr + '  ' + str(self.pressure[k]).rjust(6) + '\t' + 'static'.rjust(10)[:10] + '\t' + str(self.voleqmin[k]).rjust(10)[:10] + '\t' + str(g[k]).rjust(10)[:10] + '\t     ' + str(self.bulkmod[k]).rjust(10)[:10] + '\n'
        #
        # .....Arrays to save EOS variables
        #
        self.d2EnergydVolume2_static = []
        self.ust = []
        F = [0.0 for i in range(self.ndata)]
        self.d2EnergydVolume2_dynamic = [0.0 for k in range(self.npressure)]
        self.pstatic = [0.0 for k in range(self.npressure)]
        self.gamma_G = [0.0 for k in range(self.npressure)]
        self.pfit = [0.0 for k in range(self.npressure)]
        self.alpha = [0.0 for k in range(self.npressure)]
        theta = [0.0 for k in range(self.npressure)]
        Cv = [0.0 for k in range(self.npressure)]
        self.astatic = [0.0 for i in range(self.birchfitorder_iG+1)]
        pst = []
        fpol = []
        ferr = []
        self.gfepev = [[0.0 for k in range(self.npressure)] for j in range(self.ntemperature)]
        self.xminsav = [[0.0 for k in range(self.npressure)] for j in range(self.ntemperature)]
        self.volminsav = [[0.0 for k in range(self.npressure)] for j in range(self.ntemperature)]
        self.vol0p = []
        g0p = []
        self.b0p = []
        bp0p = []
        bpp0p = []
        self.gfe0p = []
        self.gfe0pev = []
        tpol = []
        bcnt_xk_sp = []
        bcnt_xp_sp = []
        bcnt_xv_sp = []
        bcnt_xg_sp = []
        self.gamma_poly = []
        self.cv0p = []
        self.cvu0p = []
        self.uvib0p = []
        self.uvib0pmev = []
        self.ent0p = []
        self.ent0pmev = []
        self.ent0pu = []
        self.tdeb0p = []
        self.ga0p = []
        self.he0p = []
        self.he0pmev = []
        self.gvfe0p = []
        self.gvfe0pev = []
        self.cvp = [[0.0 for k in range(self.npressure)] for j in range(self.ntemperature)]
        self.cvup = [[0.0 for k in range(self.npressure)] for j in range(self.ntemperature)]
        self.uvibp = [[0.0 for k in range(self.npressure)] for j in range(self.ntemperature)]
        self.uvibpmev = [[0.0 for k in range(self.npressure)] for j in range(self.ntemperature)]
        self.entp = [[0.0 for k in range(self.npressure)] for j in range(self.ntemperature)]
        self.entpmev = [[0.0 for k in range(self.npressure)] for j in range(self.ntemperature)]
        self.entup = [[0.0 for k in range(self.npressure)] for j in range(self.ntemperature)]
        self.tdebp = [[0.0 for k in range(self.npressure)] for j in range(self.ntemperature)]
        self.gap = [[0.0 for k in range(self.npressure)] for j in range(self.ntemperature)]
        self.hep = [[0.0 for k in range(self.npressure)] for j in range(self.ntemperature)]
        self.hepmev = [[0.0 for k in range(self.npressure)] for j in range(self.ntemperature)]
        self.a0p = []
        self.bs0p = []
        self.cp0p = []
        #
        #.....Debye temperature: numerical evaluation of thermal properties
        #
        if self.ieos == 0:
            iG = 2
            # Get static Birch-Murnaghan coefficients
            self.birch_murnaghan_eos(vol0pres, gfe0pres/self.hy2kjmol, iG, True)
            # Reinitialize d2E/dV2 
            self.d2EnergydVolume2_static = []
            # Call numerical equation of state fitting
            self.numerical_eos(volref, nepol, epol, nepol, epol, True)
            self.outstr = self.outstr + '\n'
            self.outstr = self.outstr + 'Debye temperature - numerical derivatives \n'
        #
        #.....Vinet equation of state.
        #
        elif self.ieos == 1:
            self.vinet_eos(vol0pres, gfe0pres/self.hy2kjmol, True)
            self.outstr = self.outstr + '\n'
            self.outstr = self.outstr + 'Debye temperature - Vinet EOS derivatives \n'
        #
        #.....Birch-Murnaghan equation of state.
        #
        elif self.ieos == 2:
            iG = 2
            self.birch_murnaghan_eos(vol0pres, gfe0pres/self.hy2kjmol, iG, True)
            self.outstr = self.outstr + '\n'
            self.outstr = self.outstr + 'Debye temperature - Birch-Murnaghan EOS derivatives \n'
        #
        #.....Vinet equation of state with numerical calculation to obtain all properties.
        #
        elif self.ieos == 3:
            self.vinet_eos(vol0pres, gfe0pres/self.hy2kjmol, True)
            self.numerical_eos(volref, nepol, epol, nepol, epol, True)
            self.outstr = self.outstr + '\n'
            self.outstr = self.outstr + 'Debye temperature - numerical derivatives \n'
        #
        #.....Birch-Murnaghan equation of state with numerical calculation to obtain all properties.
        #
        elif self.ieos == 4:
            iG = 2
            self.birch_murnaghan_eos(vol0pres, gfe0pres/self.hy2kjmol, iG, True)
            self.numerical_eos(volref, nepol, epol, nepol, epol, True)
            self.outstr = self.outstr + '\n'
            self.outstr = self.outstr + 'Debye temperature - numerical derivatives \n'
        #
        #.....Equation of state of Baonza-Caceres-Nunez-Taravillo.
        #
        elif self.ieos == 5:
            self.bcn_eos(vol0pres, gfe0pres/self.hy2kjmol, binp_bcnt, True)
            self.outstr = self.outstr + '\n'
            self.outstr = self.outstr + 'Debye temperature - BCNT EOS derivatives \n'
        #
        #.....Equation of state of Baonza-Caceres-Nunez-Taravillo,
        #     but numerical calculation to obtain all properties.
        #
        elif self.ieos == 6:
            self.bcn_eos(vol0pres, gfe0pres/self.hy2kjmol, binp_bcnt, True)
            self.numerical_eos(volref, nepol, epol, nepol, epol, True)
            self.outstr = self.outstr + '\n'
            self.outstr = self.outstr + 'Debye temperature - numerical EOS derivatives \n'
        #
        #.....Write Poisson coefficient and poisson ratio function
        #
        if self.ieos >= 0 and (self.idebye == 0 or self.idebye >= 2):
            self.outstr = self.outstr + 'Poisson coefficient: ' + str(self.poissonratio) + ', Poisson ratio function: ' + str(self.poissonratiofunction) + '\n'
            self.outstr = self.outstr + '\n'
        #
        #.....Calculate Debye temperatures at each volume
        #
        if self.dbulkmoddp0pres < 0:
            self.outstr = self.outstr + 'gibbs: Warning! B''<0, will use idebye=3 \n'
            self.idebye = 3
        if self.idebye == 1:
            if self.ieos >= 0:
                self.outstr = self.outstr + '  V(bohr^3) \t   TDebye(K) \t   Computed(K) \n'
                self.outstr = self.outstr + '----------- \t ----------- \t ------------- \n'
            ntpol, tpol, logstr, fterr = eos_polynomial_inst.debfitt(self.vol_inp, self.tdebye, self.ndata, self.mpar, self.mfit, self.logstr)
            if self.fterr != 0:
                if self.fterr == 2:
#                    self.logstr = self.logstr + "MP Eqn of State Thermal: Problem inverting matrix to fit polynomial \n"
                    self.grerr = 4
                    raise ValueError("MP Eqn of State Thermal: Problem inverting matrix to fit polynomial")
                else:
#                    self.logstr = self.logstr + "MP Eqn of State Thermal: No polynomial fit found with a minimum within bounds of input data \n"
                    self.grerr = 2
                    raise ValueError("MP Eqn of State Thermal: No polynomial fit found with a minimum within bounds of input data")
                return
        elif self.idebye == 3:
            tdebyemin = ((6*self.pi*self.pi*self.natoms*self.vol_inp[imin]*self.vol_inp[imin])**self.third) / self.pckbau * self.poissonratiofunction * math.sqrt(math.fabs(self.d2EnergydVolume2_static[imin])/self.cellmassval)
        else:
            if self.ieos >= 0:
                self.outstr = self.outstr + "   V(bohr^3) \t  TDebye(K) \n"
                self.outstr = self.outstr + " ----------- \t ----------\n"
        ij = 0
        # Checks second derivative is positive at all points
        # If not, the function is convex at that point and that point is skipped
        for i in xrange(self.ndata):
            if self.idebye != 1 and self.d2EnergydVolume2_static[i] <= 0.0:
                if i > itry:
                    self.ndata = i-1
                    self.logstr = self.logstr + "MP Eqn of State Thermal: Warning! convex function at i = " + str(i) + " \n"
                    self.logstr = self.logstr + "MP Eqn of State Thermal: Warning! d2EnergydVolume2_static = " + str(self.d2EnergydVolume2_static[i]) + " \n"
                    self.logstr = self.logstr + "MP Eqn of State Thermal: The points following this one will be discarded \n"
                else:
                    self.logstr = self.logstr + "MP Eqn of State Thermal: Warning! convex function at i = " + str(i) + " \n"
                    self.logstr = self.logstr + "MP Eqn of State Thermal: Warning! d2EnergydVolume2_static = " + str(self.d2EnergydVolume2_static[i]) + " \n"
                    self.logstr = self.logstr + "MP Eqn of State Thermal: This point will be skipped and the next concave one taken as the next point \n"
                    self.logstr = self.logstr + "MP Eqn of State Thermal: Recommend increasing the number of k-points and rerunning MP Eqn of State Thermal \n"
            else:
                # tmp is the Debye temperature for the structure with volume self.vol_inp[ij]
                tmp = ((6*self.pi*self.pi*self.natoms*self.vol_inp[ij]*self.vol_inp[ij])**self.third) / self.pckbau * self.poissonratiofunction * math.sqrt(math.fabs(self.d2EnergydVolume2_static[ij])/self.cellmassval)
                if self.idebye == 3:
                    self.tdebye.append(tdebyemin)
                    if self.ieos >= 0:
                        self.outstr = self.outstr + '  ' + str(self.vol_inp[ij]).rjust(10)[:10] + '\t ' + str(self.tdebye[ij]).rjust(10)[:10] + '\t' + str(tmp).rjust(10)[:10] + '\n'
                elif self.idebye == 1:
                    if self.ieos >= 0:
                        self.outstr = self.outstr + '  ' + str(self.vol_inp[ij]).rjust(10)[:10] + '\t  ' + str(self.tdebye[ij]).rjust(10)[:10] + '\n'
                else:
                    self.tdebye.append(tmp)
                if self.ieos >= 0:
                    self.outstr = self.outstr + '  ' + str(self.vol_inp[ij]).rjust(10)[:10] + '\t ' + str(self.tdebye[ij]).rjust(10)[:10] + '\n'
                ij = ij + 1
        if self.ndata != ij:
            self.ndata = ij
            self.logstr = self.logstr + "MP Eqn of State Thermal: Number of (E, V) points set to = " + str(ij) + " \n"
        #
        #.....End of static calculation: exit if IDEBYE = -1
        #
        if self.idebye == -1:
            self.logstr = self.logstr + "MP Eqn of State Thermal: End of static run ok! \n"
            return
        #
        #.....Initialize self-consistent Debye variables
        #     Debye temperature is calculated self-consistently if IDEBYE = 2
        #
        elif self.idebye == 2:
            dzero = 0.0
            scerr = self.scdebye(dzero, epol, eerr, nepol, imin, true, pst)
            if scerr:
#                self.logstr = self.logstr + "MP Eqn of State Thermal: Problem in self-consistent Debye! \n"
                scerr = False
                self.grerr = 3
                raise ValueError("MP Eqn of State Thermal: Problem in self-consistent Debye!")
                return
        #
        #.....If iopt_g = 2 the variable opt_g is changed to false
        #    Optimization of beta is only performed for the static calculation
        #    Optimized value of beta from the static calculation used for the finite temperature calculations
        #
        if self.ieos == 5 or self.ieos == 6:
            if self.iopt_g == 2:
                self.opt_g = False
        # String to save thermal properties in a plottable format
        if self.idebye >= 0:
            self.outstr_thermal = "#   T(K)    U(meV/cell)     F(meV/cell)      S(kB/cell)     Cv(kB/cell)      Theta_D(K)     Gruneisen parameter \n"
        #
        #.....Loop over temperatures
        #
        j = 0
        eos_thermal_functions_inst = eos_thermal_functions()
        while j < self.ntemperature:
            self.logstr = self.logstr + "MP Eqn of State Thermal: Temperature = " + str(self.temperature[j]) + "K \n"
            #
            #.....Self consistent Debye temperatures
            #
            if self.idebye == 2:
                scerr = self.scdebye(self.temperature[j], fpol, ferr, nfpol, imin, false, pst)
                # If scdebye returns an error for zero temperature, MP Eqn of State Thermal exits giving an error
                # If scdebye returns an error for T > zero, MP Eqn of State Thermal resets the maximum temperature to the previous value and skips the rest of the temperature loop
                # It then continues to complete the remainder of the MP Eqn of State Thermal algorithm
                if scerr:
                    if j == 0:
                        self.logstr = self.logstr + "MP Eqn of State Thermal: ffitt,  T = " + str(self.temperature[j]) + ", itry = " + str(itry) + ", imin = " + str(imin) + ", self.ndata = " + str(self.ndata) + " \n"
                        self.logstr = self.logstr + "MP Eqn of State Thermal: self.datatofit = " + str(self.datatofit) + " \n"
                        self.logstr = self.logstr + "MP Eqn of State Thermal: Cannot find minimum for temperature = " + str(self.temperature[j]) + " \n"
                        self.grerr = 2
                        return
                    else:
                        self.logstr = self.logstr + "MP Eqn of State Thermal: ffitt,  T = " + str(self.temperature[j]) + ", itry = " + str(itry) + ", imin = " + str(imin) + ", self.ndata = " + str(self.ndata) + " \n"
                        self.logstr = self.logstr + "MP Eqn of State Thermal: self.datatofit = " + str(self.datatofit) + " \n"
                        self.logstr = self.logstr + "MP Eqn of State Thermal: Cannot find minimum of E + A at T = " + str(self.temperature[j]) + " \n"
                        self.logstr = self.logstr + "MP Eqn of State Thermal: Resetting maximum value of T to T = " + str(self.temperature[j]) + " \n"
                        self.logstr = self.logstr + "MP Eqn of State Thermal: Resetting number of temperature points to " + str(j) + " \n"
                        self.ntemperature = j
                        break
                ntpol, tpol, self.logstr, self.fterr = eos_polynomial_inst.debfitt(self.vol_inp, self.tdebye, self.ndata, self.mpar, self.mfit, self.logstr)
                # If debfitt returns an error for zero temperature, MP Eqn of State Thermal exits giving an error
                # If debfitt returns an error for T > zero, MP Eqn of State Thermal resets the maximum temperature to the previous value and skips the rest of the temperature loop
                # It then continues to complete the remainder of the MP Eqn of State Thermal algorithm
                if self.fterr != 0:
                    if j == 0:
                        self.logstr = self.logstr + "MP Eqn of State Thermal: No polynomial fit found with a minimum within bounds of input data \n"
                        self.logstr = self.logstr + "MP Eqn of State Thermal: ffitt,  T = " + str(self.temperature[j]) + ", itry = " + str(itry) + ", imin = " + str(imin) + ", self.ndata = " + str(self.ndata)  + " \n"
                        self.logstr = self.logstr + "MP Eqn of State Thermal: Cannot find minimum for temperature = " + str(self.temperature[j]) + " \n"
                        if self.fterr == 2:
                            self.logstr = self.logstr + "MP Eqn of State Thermal: Problem inverting matrix to fit polynomial \n"
                            self.grerr = 4
                        else:
                            self.grerr = 2
                        return
                    else:
                        self.logstr = self.logstr + "MP Eqn of State Thermal: No polynomial fit found with a minimum within bounds of input data \n"
                        self.logstr = self.logstr + "MP Eqn of State Thermal: ffitt,  T = " + str(self.temperature[j]) + ", itry = " + str(itry) + ", imin = " + str(imin) + ", self.ndata = "  + str(self.ndata) + " \n"
                        self.logstr = self.logstr + "MP Eqn of State Thermal: Resetting maximum value of T to T = " + str(self.temperature[j]) + " \n"
                        self.logstr = self.logstr + "MP Eqn of State Thermal: Resetting number of temperature points to " + str(j) + " \n"
                        self.ntemperature = j
                        break
            #
            #.....Obtain vibrational Helmholtz function fit (Debye model) for this value of the temperature, T
            #     At each volume, use the Debye temperature to get the vibrational Helmholtz energy (A = U - TS)
            #     Add the value of A at each volume to the DFT final energy at that volume to self.datatofit
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
            self.thermal_helmholtz = 0.0
            self.thermal_entropy = 0.0
            self.thermal_energ = 0.0
            self.thermal_cv = 0.0
            for i in xrange(self.ndata):
                self.thermal_energ, self.thermal_entropy, self.thermal_helmholtz, self.thermal_cv, D, Derr = eos_thermal_functions_inst.thermal(self.tdebye[i], self.temperature[j], self.natoms, self.pckbau, self.maxloops, self.logstr)
                self.datatofit[i] = self.energ_inp[i] + self.thermal_helmholtz
                F[i] = self.thermal_helmholtz
            self.ierr, imin, self.logstr = eos_polynomial_inst.minbrack(imin, self.datatofit, self.ndata, self.logstr)
            # If minbrack returns an error for zero temperature, MP Eqn of State Thermal exits giving an error
            # If minbrack returns an error for T > zero, MP Eqn of State Thermal resets the maximum temperature to the previous value and skips the rest of the temperature loop
            # It then continues to complete the remainder of the MP Eqn of State Thermal algorithm
            if self.ierr != 0:
                if j == 0:
                    self.logstr = self.logstr + "MP Eqn of State Thermal: ffitt,  T = " + str(self.temperature[j]) + ", itry = " + str(itry) + ", imin = " + str(imin) + ", self.ndata = " + str(self.ndata) + " \n"
                    self.logstr = self.logstr + "MP Eqn of State Thermal: self.datatofit = " + str(self.datatofit) + " \n"
                    self.logstr = self.logstr + "MP Eqn of State Thermal: Cannot find minimum for temperature = " + str(self.temperature[j]) + " \n"
                    self.grerr = 2
                    raise ValueError("MP Eqn of State Thermal: Cannot find minimum of Gibbs free energy")
                    return
                else:
                    self.logstr = self.logstr + "MP Eqn of State Thermal: ffitt,  T = " + str(self.temperature[j]) + ", itry = " + str(itry) + ", imin = " + str(imin) + ", self.ndata = " + str(self.ndata) + " \n"
                    self.logstr = self.logstr + "MP Eqn of State Thermal: self.datatofit = " + str(self.datatofit) + " \n"
                    self.logstr = self.logstr + "MP Eqn of State Thermal: Cannot find minimum of E + A at T = " + str(self.temperature[j]) + " \n"
                    self.logstr = self.logstr + "MP Eqn of State Thermal: Resetting maximum value of T to T = " + str(self.temperature[j]) + " \n"
                    self.logstr = self.logstr + "MP Eqn of State Thermal: Resetting number of temperature points to " + str(j) + " \n"
                    self.ntemperature = j
                    break
            itry = imin
            # Fits polynomial to E + A
            self.fterr, nfpol, fpol, ferr = eos_polynomial_inst.polynomial_fit(imin, self.xconfigvector, self.datatofit, self.ndata, self.mpar, self.mfit)
            # If fitt returns an error for zero temperature, MP Eqn of State Thermal exits giving an error
            # If fitt returns an error for T > zero, MP Eqn of State Thermal resets the maximum temperature to the previous value and skips the rest of the temperature loop
            # It then continues to complete the remainder of the MP Eqn of State Thermal algorithm
            if self.fterr != 0:
                if j == 0:
                    self.logstr = self.logstr + "MP Eqn of State Thermal: No polynomial fit found with a minimum within bounds of input data \n"
                    self.logstr = self.logstr + "MP Eqn of State Thermal: ffitt,  T = " + str(self.temperature[j]) + ", itry = " + str(itry) + ", imin = " + str(imin) + ", self.ndata = " + str(self.ndata) + " \n"
                    self.logstr = self.logstr + "MP Eqn of State Thermal: Cannot find minimum for temperature = " + str(self.temperature[j]) + " \n"
                    if self.fterr == 2:
#                        self.logstr = self.logstr + "MP Eqn of State Thermal: Problem inverting matrix to fit polynomial \n"
                        self.grerr = 4
                        raise ValueError("MP Eqn of State Thermal: Problem inverting matrix to fit polynomial")
                    else:
                        self.grerr = 2
                        raise ValueError("MP Eqn of State Thermal: Cannot find minimum of Gibbs free energy")
                    return
                else:
                    self.logstr = self.logstr + "MP Eqn of State Thermal: No polynomial fit found with a minimum within bounds of input data \n"
                    self.logstr = self.logstr + "MP Eqn of State Thermal: ffitt,  T = " + str(self.temperature[j]) + ", itry = " + str(itry) + ", imin = " + str(imin) + ", self.ndata = " + str(self.ndata) + " \n"
                    self.logstr = self.logstr + "MP Eqn of State Thermal: Resetting maximum value of T to T = " + str(self.temperature[j-1]) + " \n"
                    self.logstr = self.logstr + "MP Eqn of State Thermal: Resetting number of temperature points to " + str(j) + " \n"
                    self.ntemperature = j
                    break
            if nfpol == 2:
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
#        for k in xrange(self.npressure):
            while k < self.npressure:
                self.logstr = self.logstr + "MP Eqn of State Thermal: Pressure = " + str(self.pressure[k]) + " GPa \n"
                #
                #.....Calculate the value of the Gibbs free energy E + A + pV for each volume self.vol_inp.at(i) and store in self.datatofit
                #     Bracket the minimum of E + A + pV  with the numerical function
                #     Gibbs free energy = A + pV; constant pressure, variable volume
                #
                for i in xrange(self.ndata):
                    self.datatofit[i] = self.energ_inp[i] + F[i] + self.vol_inp[i]*self.pressure[k]/self.au2gpa
                self.ierr, imin, self.logstr = eos_polynomial_inst.minbrack(itry, self.datatofit, self.ndata, self.logstr)
                # If minbrack returns an error for zero pressure and zero temperature, MP Eqn of State Thermal exits giving an error
                # If minbrack returns an error for p = 0, T > 0, MP Eqn of State Thermal resets the maximum temperature to the previous value and skips the rest of the loop
                # If minbrack returns an error for p > zero, MP Eqn of State Thermal resets the maximum pressure to the previous value and skips the rest of the pressure loop
                # It then continues to complete the remainder of the MP Eqn of State Thermal algorithm
                if self.ierr != 0:
                    if k == 0:
                        if j == 0:
                            self.logstr = self.logstr + "MP Eqn of State Thermal: gmin, T = " + str(self.temperature[j]) + ", P = " + str(self.pressure[k]) + ", trial point = " + str(itry) + ", minimum point = " + str(imin) + ", total points = " + str(self.ndata) + " \n"
                            self.logstr = self.logstr + "MP Eqn of State Thermal: func = " + str(self.datatofit[i]) + " \n"
                            self.logstr = self.logstr + "MP Eqn of State Thermal: Cannot find minimum for pressure = " + str(self.pressure[k]) + " and temperature = " + str(self.temperature[j]) + " \n"
                            self.grerr = 2
                            return
                        else:
                            self.logstr = self.logstr + "MP Eqn of State Thermal: gmin, T = " + str(self.temperature[j]) + ", P = " + str(self.pressure[k]) + ", trial point = " + str(itry) + ", minimum point = " + str(imin) + ", total points = " + str(self.ndata) + " \n"
                            self.logstr = self.logstr + "MP Eqn of State Thermal: self.datatofit = " + str(self.datatofit[i]) + " \n"
                            self.logstr = self.logstr + "MP Eqn of State Thermal: Cannot find minimum of E + pV + A at T = " + str(self.temperature[j]) + ", p = " + str(self.pressure[k]) + " \n"
                            self.logstr = self.logstr + "MP Eqn of State Thermal: Resetting maximum value of T to T = " + str(self.temperature[j-1]) + " \n"
                            self.logstr = self.logstr + "MP Eqn of State Thermal: Resetting number of pressure points to " + str(j) + " \n"
                            self.ntemperature = j
                            break
                    else:
                        self.logstr = self.logstr + "MP Eqn of State Thermal: gmin, T = " + str(self.temperature[j]) + ", P = " + str(self.pressure[k]) + ", trial point = " + str(itry) + ", minimum point = " + str(imin) + ", total points = " + str(self.ndata) + " \n"
                        self.logstr = self.logstr + "MP Eqn of State Thermal: self.datatofit = " + str(self.datatofit[i]) + " \n"
                        self.logstr = self.logstr + "MP Eqn of State Thermal: Cannot find minimum of E + pV + A at T = " + str(self.temperature[j]) + ", p = " + str(self.pressure[k]) + " \n"
                        self.logstr = self.logstr + "MP Eqn of State Thermal: Resetting maximum value of p to p = " + str(self.pressure[k-1]) + " \n"
                        self.logstr = self.logstr + "MP Eqn of State Thermal: Resetting number of pressure points to " + str(k) + " \n"
                        self.npressure = k
                        break
                #
                #.....Find the volume which minimizes the fitted G function at (p, T) (polynomial fit is in gpol)
                #     G function is Gibbs free energy: G = E + pV + A; E = DFT energy
                #     For a given temperature and pressure, the equilibrium system is the one which minimizes the Gibbs free energy
                #
                gpol[3] = fpol[3] + self.pressure[k]/self.au2gpa * volref
                self.pmerr, xmin, self.logstr = eos_polynomial_inst.polmin(self.xconfigvector[itry], self.xconfigvector[max(itry-2, 0)], self.xconfigvector[min(itry+2, self.ndata-1)], ngpol, gpol, self.pmerr, self.logstr)
                # If polmin returns an error, then MP Eqn of State Thermal tries shifting the trial point to try to correct the error
                if self.pmerr != 0:
                    self.logstr = self.logstr + "MP Eqn of State Thermal: itry = " + str(itry) + " \n"
                    self.logstr = self.logstr + "MP Eqn of State Thermal: pmerr = " + str(self.pmerr) + " \n"
                    if self.pmerr == 1:
                        while (self.pmerr == 1) and ((itry + 2) < (self.ndata-1)):
                            itry = itry + 1
                            self.logstr = self.logstr + "MP Eqn of State Thermal: Resetting itry to itry = " + str(itry) + " \n"
                            self.pmerr, xmin, self.logstr = eos_polynomial_inst.polmin(self.xconfigvector[itry], self.xconfigvector[max(itry-2, 0)], self.xconfigvector[min(itry+2, self.ndata-1)], ngpol, gpol, self.pmerr, self.logstr)
                            self.logstr = self.logstr + "MP Eqn of State Thermal: pmerr = " + str(self.pmerr) + " \n"
                        # If error indicator has changed from 1 to 2, then the bracket has shifted from one side of the minimum to the other
                        # Need to expand the size of the bracket to incorporate the minimum
                        if self.pmerr == 2:
                            self.pmerr, xmin, self.logstr = eos_polynomial_inst.polmin(self.xconfigvector[itry], self.xconfigvector[max(itry-4, 0)], self.xconfigvector[min(itry+4, self.ndata-1)], ngpol, gpol, self.pmerr, self.logstr)
                        # If polynomial minimum has still not been found successfully, writes polynomial and its first derivative to a file to aid debugging
                        if self.pmerr != 0:
                            self.logstr = self.logstr + "MP Eqn of State Thermal: List of polynomial values \n"
                            for i in xrange(self.ndata):
                                plnv = eos_polynomial_inst.polin0(self.xconfigvector[i], ngpol, gpol)
                                self.logstr = self.logstr + str(self.xconfigvector[i]) + '\t' + str(plnv) + '\n'
                            self.logstr = self.logstr + "MP Eqn of State Thermal: List of first derivative values of polynomial \n"
                            for i in xrange(self.ndata):
                                plnv = eos_polynomial_inst.polin1(self.xconfigvector[i], ngpol, gpol)
                                self.logstr = self.logstr + str(self.xconfigvector[i]) + '\t' + str(plnv) + '\n'
                    elif self.pmerr == 2:
                        while (self.pmerr == 2) and ((itry - 2) >= 0):
                            itry = itry - 1
                            self.logstr = self.logstr + "MP Eqn of State Thermal: Resetting itry to itry = " + str(itry) + " \n"
                            self.pmerr, xmin, self.logstr = eos_polynomial_inst.polmin(self.xconfigvector[itry], self.xconfigvector[max(itry-2, 0)], self.xconfigvector[min(itry+2, self.ndata-1)], ngpol, gpol, self.pmerr, self.logstr)
                            self.logstr = self.logstr + "MP Eqn of State Thermal: pmerr = " + str(self.pmerr) + " \n"
                        # If error indicator has changed from 2 to 1, then the bracket has shifted from one side of the minimum to the other
                        # Need to expand the size of the bracket to incorporate the minimum
                        if self.pmerr == 1:
                            self.pmerr, xmin, self.logstr = eos_polynomial_inst.polmin(self.xconfigvector[itry], self.xconfigvector[max(itry-4, 0)], self.xconfigvector[min(itry+4, self.ndata-1)], ngpol, gpol, self.pmerr, self.logstr)
                        # If polynomial minimum has still not been found successfully, writes polynomial and its first derivative to a file to aid debugging
                        if self.pmerr != 0:
                            self.logstr = self.logstr + "MP Eqn of State Thermal: List of polynomial values \n"
                            for i in xrange(self.ndata):
                                plnv = eos_polynomial_inst.polin0(self.xconfigvector[i], ngpol, gpol)
                                self.logstr = self.logstr + str(self.xconfigvector[i]) + '\t' + str(plnv) + '\n'
                                self.logstr = self.logstr + "MP Eqn of State Thermal: List of first derivative values of polynomial \n"
                        for i in xrange(self.ndata):
                            plnv = eos_polynomial_inst.polin1(self.xconfigvector[i], ngpol, gpol)
                            self.logstr = self.logstr + str(self.xconfigvector[i]) + '\t' + str(plnv) + '\n'
                # Check that polmin has actually found a minimum and not a maximum
                plnv = eos_polynomial_inst.polin2(xmin, ngpol, gpol)
                if plnv < 0.0:
                    self.logstr = self.logstr + "MP Eqn of State Thermal: Minimum polynomial fitted to E + pV is actually a maximum \n"
                    self.logstr = self.logstr + "MP Eqn of State Thermal: Rebracketing minimum of polynomial points \n"
                    self.ierr, itry = eos_polynomial_inst.polminbrack(itry, ngpol, self.ndata, self.xconfigvector, gpol)
                    if self.ierr != 0:
                        if k == 0:
                            if j == 0:
                                self.logstr = self.logstr + "MP Eqn of State Thermal: Cannot find minimum of E + A + pV for p = " + str(self.pressure[k]) + ", T = " + str(self.temperature[j]) + " \n"
                                self.grerr = 2
                                raise ValueError("MP Eqn of State Thermal: Cannot find minimum of Gibbs free energy")
                                return
                            else:
                                self.logstr = self.logstr + "MP Eqn of State Thermal: Cannot find minimum of E + A + pV for p = " + str(self.pressure[k]) + ", T = " + str(self.temperature[j]) + " \n"
                                self.logstr = self.logstr + "MP Eqn of State Thermal: Resetting maximum value of T to T = " + str(self.temperature[j-1]) + " \n"
                                self.logstr = self.logstr + "MP Eqn of State Thermal: Resetting maximum number of temperature points to " + str(j) + " \n"
                                self.ntemperature = j
                                break
                        else:
                            self.logstr = self.logstr + "MP Eqn of State Thermal: Cannot find minimum of E + A + pV for p = " + str(self.pressure[k]) + ", T = " + (self.temperature[j]) + " \n"
                            self.logstr = self.logstr + "MP Eqn of State Thermal: Resetting maximum value of p to p = " + str(self.pressure[k-1]) + " \n"
                            self.logstr = self.logstr + "MP Eqn of State Thermal: Resetting maximum number of pressure points to " + str(k) + " \n"
                            self.npressure = k
                            break
                    else:
                        self.logstr = self.logstr + "MP Eqn of State Thermal: Minimum of (E, V) data is at point imin = " + str(itry) + " \n"
                        self.pmerr, xmin, self.logstr = eos_polynomial_inst.polmin(self.xconfigvector[itry], self.xconfigvector[max(itry-2, 0)], self.xconfigvector[min(itry+2, self.ndata-1)], ngpol, gpol, self.pmerr, self.logstr)
                # If polmin still returns an error for zero temperature and pressure after shifting the trial point, MP Eqn of State Thermal exits giving an error
                # If polmin returns an error for T > zero, MP Eqn of State Thermal resets the maximum temperature or pressure to the previous value and skips the rest of that loop
                # It then continues to complete the remainder of the MP Eqn of State Thermal algorithm
                if self.pmerr != 0:
                    if k == 0:
                        if j == 0:
                            self.logstr = self.logstr + "MP Eqn of State Thermal: Cannot find minimum of polynomial fit for E + pV + A for p = " + str(self.pressure[k]) + " and T = " + str(self.temperature[j]) + " \n"
                            self.grerr = 2
                            raise ValueError("MP Eqn of State Thermal: Cannot find minimum of Gibbs free energy")
                            return
                        else:
                            self.logstr = self.logstr + "MP Eqn of State Thermal: Cannot find minimum of polynomial fit for E + pV + A at T = " + str(self.temperature[j]) + ", p = " + str(self.pressure[k]) + " \n"
                            self.logstr = self.logstr + "MP Eqn of State Thermal: Resetting maximum temperature value to T = " + str(self.temperature[j-1]) + " \n"
                            self.logstr = self.logstr + "MP Eqn of State Thermal: Resetting number of temperature points to " + str(j) + " \n"
                            self.ntemperature = j
                            break
                    else:
                        self.logstr = self.logstr + "MP Eqn of State Thermal: Cannot find minimum of polynomial fit for E + pV + A at T = " + str(self.temperature[j]) + ", p = " + str(self.pressure[k]) + " \n"
                        self.logstr = self.logstr + "MP Eqn of State Thermal: Resetting maximum pressure value to p = " + str(self.pressure[k-1]) + " \n"
                        self.logstr = self.logstr + "MP Eqn of State Thermal: Resetting number of pressure points to " + str(k) + " \n"
                        self.npressure = k
                        break
                # Evaluates polynomial at minimum to get equilibrium V, G, and B
                self.voleqmin[k] = (xmin**3) * volref
                plnv = eos_polynomial_inst.polin0(xmin, ngpol, gpol)
                g[k] = plnv * self.hy2kjmol
                self.gfepev[j][k] = plnv*self.hart2ev
                plnv = eos_polynomial_inst.polin2(xmin, ngpol, gpol)
                self.bulkmod[k] = plnv * self.au2gpa * xmin*xmin/(9.0*self.voleqmin[k])
                eact = eos_polynomial_inst.polin0(xmin, nfpol, fpol)
                twonfpol = 2 * nfpol
                plnv = eos_polynomial_inst.polin0(xmin, twonfpol, ferr)
                aerr = math.sqrt(math.fabs(plnv - eact*eact))
                rerr[k] = aerr / max(math.fabs(eact), math.fabs(eact)+aerr/2.0)
                # Writes information about equilibrium volume to string
                self.xminsav[j][k] = xmin
                self.volminsav[j][k] = self.voleqmin[k]
                k = k + 1
            #
            #.....ieos < 0 means minimum output
            #
            if self.ieos < 0:
                for k in xrange(self.npressure):
                    self.outstr = self.outstr + str(self.pressure[k]) + "\t" + str(self.temperature[j]) + "\t" + str(self.voleqmin[k]) + "\t" + str(g[k]) + "\t" + str(self.bulkmod[k]) + "\n"
            else:
                #
                #.....Numerical energetic, geometric and elastic properties
                #
                vol0pres = self.voleqmin[0]
                gfe0pres = g[0]
                binp_bcnt = self.bulkmod[0]
                self.outstr = self.outstr + "\n"
                self.outstr = self.outstr + "Temperature:  T = " + str(self.temperature[j]) + "K \n"
                self.outstr = self.outstr + "Vmin(T; P=0) = " + str(vol0pres) + "bohr^3 \n"
                self.outstr = self.outstr + "Gmin(T; P=0) = " + str(gfe0pres) + "kJ/mol \n"
                self.outstr = self.outstr + "\n"
                self.outstr = self.outstr + "NUMERICAL EQUILIBRIUM PROPERTIES \n"
                self.outstr = self.outstr + "================================ \n"
                self.outstr = self.outstr + "\n"
                self.outstr = self.outstr + "  P(GPa) \t G(kJ/mol) \t V(bohr^3) \t    V/V0 \t    B(GPa) \t   rel. err.  \n"
                self.outstr = self.outstr + " -------------------------------------------------------------------------------------------- \n"
                for k in xrange(self.npressure):
                    self.outstr = self.outstr + '  ' + str(self.pressure[k]).rjust(6) + "\t" + str(g[k]).rjust(10)[:10] + "\t" + str(self.voleqmin[k]).rjust(10)[:10] + "\t" + str(self.voleqmin[k]/vol0pres).rjust(8)[:8] + "\t" + str(self.bulkmod[k]).rjust(10)[:10] + "\t" + str(round(rerr[k], 10)).rjust(12)[:12] + "\n"
            #
            #.....Finite temperature numerical dynamic results.
            #
            if self.ieos == 0:
                self.numerical_eos(volref, nepol, epol, nfpol, fpol, False)
            #
            #.....Vinet equation of state.
            #
            elif self.ieos == 1:
                self.vinet_eos(vol0pres, gfe0pres, False)
            #
            #.....Birch-Murnaghan equation of state.
            #
            elif self.ieos == 2:
                iG = 2
                self.birch_murnaghan_eos(vol0pres, gfe0pres, iG, False)
            #
            #.....Vinet equation of state with numerical calculation of properties.
            #
            elif self.ieos == 3:
                self.vinet_eos(vol0pres, gfe0pres/self.hy2kjmol, False)
                self.numerical_eos(volref, nepol, epol, nepol, epol, False)
            #
            #.....Birch-Murnaghan equation of state with numerical calculation of properties.
            #
            elif self.ieos == 4:
                iG = 2
                self.birch_murnaghan_eos(vol0pres, gfe0pres/self.hy2kjmol, iG, False)
                self.numerical_eos(volref, nepol, epol, nepol, epol, False)
            #
            #.....Equation of state of Baonza-Caceres-Nunez-Taravillo.
            #
            elif self.ieos == 5:
                self.bcn_eos(vol0pres, gfe0pres/self.hy2kjmol, binp_bcnt, False)
                bcnt_xk_sp.append(self.xsupa)
                bcnt_xp_sp.append(self.pspin)
                bcnt_xv_sp.append(self.vspin)
                bcnt_xg_sp.append(self.beta)
            #
            #.....Equation of state of Baonza-Caceres-Nunez-Taravillo,
            #     but numerical calculation to obtain all properties.
            #
            elif self.ieos == 6:
                self.bcn_eos(vol0pres, gfe0pres/self.hy2kjmol, binp_bcnt, False)
                self.numerical_eos(volref, nepol, epol, nepol, epol, False)
                bcnt_xk_sp.append(self.xsupa)
                bcnt_xp_sp.append(self.pspin)
                bcnt_xv_sp.append(self.vspin)
                bcnt_xg_sp.append(self.beta)
            #
            #.....Save properties at P=0 for all of the temperatures.
            #
            self.vol0p.append(self.voleqmin[0])
            g0p.append(g[0])
            self.b0p.append(self.bulkmod0pres)
            bp0p.append(self.dbulkmoddp0pres)
            bpp0p.append(self.d2bulkmoddp20pres)
            self.gfe0p.append(g[0])
            g0pev = (g[0] / self.hy2kjmol) * self.hart2ev
            self.gfe0pev.append(g0pev)
            #
            #.....Compute vibrational properties at the equilibrium volumes.
            #
            self.outstr = self.outstr + "\n"
            self.outstr = self.outstr + "VIBRATIONAL PROPERTIES \n"
            self.outstr = self.outstr + "====================== \n"
            self.outstr = self.outstr + "\n"
            self.outstr = self.outstr + "  P(GPa) \t  U(kJ/mol) \t Cv(J/mol*K) \t A(kJ/mol) \t S(J/mol*K) \t  Theta(K) \t gamma \n"
            self.outstr = self.outstr + " ------------------------------------------------------------------------------------------------------- \n"
            k = 0
#        for k in xrange(self.npressure):
            while k < self.npressure:
                if self.idebye == 1 or self.idebye == 2:
                    tmp = math.log(self.voleqmin[k])
                    plnv = eos_polynomial_inst.polin0(tmp, ntpol, tpol)
                    theta[k] = math.exp(plnv)
                    plnv = eos_polynomial_inst.polin1(tmp, ntpol, tpol)
                    self.gamma_G[k] = - plnv
                    self.thermal_energ, self.thermal_entropy, self.thermal_helmholtz, self.thermal_cv, D, Derr = eos_thermal_functions_inst.thermal(theta[k], self.temperature[j], self.natoms, self.pckbau, self.maxloops, self.logstr)
                    Uvib = self.thermal_energ
                    Cv[k] = self.thermal_cv
                    helm = self.thermal_helmholtz
                    ent = self.thermal_entropy
                    vibu = Uvib*hy2kjmol
                    vibumev = Uvib*hart2ev*1000
                    vibcv = Cv.at(k)*hy2kjmol*1000
                    vibcvu = vibcv*kj2unit
                    vibhe = helm*hy2kjmol
                    vibhemev = helm*hart2ev*1000
                    vibent = ent*hy2kjmol*1000
                    vibentmev = ent*hart2ev*1000
                    vibentu = vibent*kj2unit
                    plnv = eos_polynomial_inst.polin0(xmin, nepol, epol)
                    vibg = plnv + helm + ((self.pressure[k]/self.au2gpa) * (xmin**3))
                    vibgev = vibg * hart2ev
                    vibg = vibg * hy2kjmol
                else:
                    # Calculate and save Gruneisen parameter using derivative as a check
                    if k == 0:
                        ntpol, tpol, self.logstr, self.fterr = eos_polynomial_inst.debfitt(self.vol_inp, self.tdebye, self.ndata, self.mpar, self.mfit, self.logstr)
                        tmp = math.log(self.voleqmin[k])
                        plnv = eos_polynomial_inst.polin1(tmp, ntpol, tpol)
                        self.gamma_poly.append(-plnv)
                    #
                    #.....if d2EnergydVolume2_dynamic is not consistent, make a Gruneisen interpolation
                    #
                    if self.d2EnergydVolume2_dynamic[k] <= 0.0:
                        self.logstr = self.logstr + "MP Eqn of State Thermal: inconsistent derivative, v[k] = " + str(self.voleqmin[k]) + ", p[k] = " + str(self.pressure[k]) + " \n"
                        self.logstr = self.logstr + "MP Eqn of State Thermal: k = ", k, ", T = " + str(self.temperature[j]) + ", d2EnergydVolume2_dynamic[k] = " + str(self.d2EnergydVolume2_dynamic[k]) + " \n"
                        self.logstr = self.logstr + "MP Eqn of State Thermal: making a Gruneisen interpolation \n"
                        ilow = 0
                        while self.voleqmin[k] >= self.vol_inp[ilow] and ilow < (self.ndata-1):
                            ilow = ilow + 1
                        if (ilow == (self.ndata-1) and self.voleqmin[k] >= self.vol_inp[ilow]) or ilow == 0:
                            if k == 0:
                                self.logstr = self.logstr + "MP Eqn of State Thermal: inconsistent volume, v[k] = " + str(self.voleqmin[k]) + ", p[k] = " + str(self.pressure[k]) + " \n"
                                self.grerr = 3
                                raise ValueError("MP Eqn of State Thermal: Inconsistent volume")
                                return
                            else:
                                self.logstr = self.logstr + "MP Eqn of State Thermal: Inconsistent volume, v[k] = " + str(self.voleqmin[k]) + ", p[k] = " + str(self.pressure[k]) + ", k = " + str(k) + " \n"
                                self.logstr = self.logstr + "MP EOS_THERMAL: Resetting maximum pressure value to p = " + str(self.pressure[k-1]) + " \n"
                                self.logstr = self.logstr + "MP Eqn of State Thermal: Resetting number of pressure points to " + str(k) + " \n"
                                self.npressure = k
                                break
                        ilow = ilow - 1
                        self.d2EnergydVolume2_dynamic[k] = self.d2EnergydVolume2_static[ilow] * ((self.voleqmin[k]/self.vol_inp[ilow])**((math.log(self.d2EnergydVolume2_static[ilow+1]/self.d2EnergydVolume2_static[ilow])) / (math.log(self.vol_inp[ilow+1]/self.vol_inp[ilow]))))
                    #
                    #.....isotropic Debye model properties
                    #
                    theta[k] = ((6*self.pi*self.pi*self.natoms*self.voleqmin[k]*self.voleqmin[k])**self.third) / self.pckbau * self.poissonratiofunction * math.sqrt(self.d2EnergydVolume2_dynamic[k]/self.cellmassval)
                    self.thermal_energ, self.thermal_entropy, self.thermal_helmholtz, self.thermal_cv, D, Derr = eos_thermal_functions_inst.thermal(theta[k], self.temperature[j], self.natoms, self.pckbau, self.maxloops, self.logstr)
                    Uvib = self.thermal_energ
                    Cv[k] = self.thermal_cv
                    helm = self.thermal_helmholtz
                    ent = self.thermal_entropy
                    vibu = Uvib*self.hy2kjmol
                    vibumev = Uvib*self.hart2ev*1000
                    vibcv = Cv[k]*self.hy2kjmol*1000
                    vibcvu = vibcv*self.kj2unit
                    vibhe = helm*self.hy2kjmol
                    vibhemev = helm*self.hart2ev*1000
                    vibent = ent*self.hy2kjmol*1000
                    vibentmev = ent*self.hart2ev*1000
                    vibentu = vibent*self.kj2unit
                    eos_polynomial_inst.polin0(xmin, nepol, epol)
                    vibg = plnv + helm + ((self.pressure[k]/self.au2gpa) * (xmin**3))
                    vibgev = vibg * self.hart2ev
                    vibg = vibg * self.hy2kjmol
                self.outstr = self.outstr + '  ' + str(self.pressure[k]).rjust(6) + "\t " + str(vibu).rjust(10)[:10] + "\t  " + str(vibcv).rjust(10)[:10] + "\t" + str(vibhe).rjust(10)[:10] + "\t " + str(vibent).rjust(10)[:10] + "\t" + str(theta[k]).rjust(10)[:10] + "\t" + str(self.gamma_G[k]).rjust(6)[:6] + "\n"
                if k == 0:
                    # Save heat capacity, internal energy, entropy, Debye temperature, Gruneisen parameter, Helmholtz free energy and Gibbs free energy at zero pressure for analysis and plotting later on
                    self.cv0p.append(vibcv)
                    self.cvu0p.append(vibcvu)
                    self.uvib0p.append(vibu)
                    self.uvib0pmev.append(vibumev)
                    self.ent0p.append(vibent)
                    self.ent0pmev.append(vibentmev)
                    self.ent0pu.append(vibentu)
                    self.tdeb0p.append(theta[k])
                    self.ga0p.append(self.gamma_G[k])
                    self.he0p.append(vibhe)
                    self.he0pmev.append(vibhemev)
                    self.gvfe0p.append(vibg)
                    self.gvfe0pev.append(vibgev)
                if self.savpres:
                    # Save heat capacity, internal energy, entropy, Debye temperature, Gruneisen parameter, and Helmholtz free energy at all pressures for analysis and plotting later on
                    # Saving these values is optional and can be activated by setting
                    self.cvp[j][k] = vibcv
                    self.cvup[j][k] = vibcvu
                    self.uvibp[j][k] = vibu
                    self.uvibpmev[j][k] = vibumev
                    self.entp[j][k] = vibent
                    self.entpmev[j][k] = vibentmev
                    self.entup[j][k] = vibentu
                    self.tdebp[j][k] = theta[k]
                    self.gap[j][k] = self.gamma_G[k]
                    self.hep[j][k] = vibhe
                    self.hepmev[j][k] = vibhemev
                k = k + 1
            #
            #.....Compute EOS derivatives
            #
            self.outstr = self.outstr + "\n"
            self.outstr = self.outstr + "THERMAL EOS DERIVATIVES \n"
            self.outstr = self.outstr + "======================= \n"
            self.outstr = self.outstr + "\n"
            self.outstr = self.outstr + "  P(GPa) \t alpha(10^5/K) \t dp/dt(GPa/K) \t Bs(GPa) \t Cp(J/mol*K) \n"
            self.outstr = self.outstr + " ---------------------------------------------------------------------------- \n"
            for k in xrange(self.npressure):
                Pbeta = Cv[k]*self.gamma_G[k]/self.voleqmin[k] * self.au2gpa
                self.alpha[k] = Pbeta/self.bulkmod[k]
                tmp = 1.0 + self.gamma_G[k]*self.alpha[k]*self.temperature[j]
                Cp = Cv[k] * tmp * self.hy2kjmol * 1000
                Bs = self.bulkmod[k] * tmp
                self.outstr = self.outstr + '  ' + str(self.pressure[k]).rjust(6) + "\t" + str(round(self.alpha[k]*1e5, 12)).rjust(14)[:14] + "\t " + str(round(Pbeta, 10)).rjust(12)[:12] + "\t" + str(Bs).rjust(8)[:8] + "\t" + str(Cp).rjust(12)[:12] + "\n"
                if k == 0:
                    self.a0p.append(self.alpha[k])
                    self.bs0p.append(Bs)
                    self.cp0p.append(Cp)
            self.outstr_thermal = self.outstr_thermal + str(self.temperature[j]).rjust(6) + "\t" + str(self.uvib0pmev[j]).rjust(10)[:10] + "\t" + str(self.he0pmev[j]).rjust(10)[:10] + "\t" + str(self.ent0pu[j]).rjust(10)[:10] + "\t" + str(self.cvu0p[j]).rjust(10)[:10] + "\t" + str(self.tdeb0p[j]).rjust(10)[:10] + "\t" + str(self.ga0p[j]).rjust(8)[:8] + "\n"
            j = j + 1
        #
        #.....Close temperature loop
        #
        #.....Write all of the results at P=0 for all of the temperatures to the stringstream for the main MP Eqn of State Thermal output file.
        #
        if self.ieos >= 0:
            self.outstr = self.outstr + "\n"
            self.outstr = self.outstr + "RESULTS AT P=0 FOR ALL TEMPERATURES \n"
            self.outstr = self.outstr +  "=================================== \n"
            self.outstr = self.outstr + "\n"
            self.outstr = self.outstr +  "    T(K) \t V(bohr^3) \t G(kJ/mol) \t U(kJ/mol) \t S(J/mol K) \t Cv(J/mol K) \n"
            self.outstr = self.outstr +  " -------------------------------------------------------------------------------------------- \n"
            for j in xrange(self.ntemperature):
                self.outstr = self.outstr + '  ' + str(self.temperature[j]).rjust(6) + "\t" + str(self.vol0p[j]).rjust(10)[:10] + "\t" + str(g0p[j]).rjust(10)[:10] + "\t" + str(self.uvib0p[j]).rjust(10)[:10] + "\t " + str(self.ent0p[j]).rjust(10)[:10] + "\t  " + str(self.cv0p[j]).rjust(10)[:10] + "\n"
            self.outstr = self.outstr + "\n"
            self.outstr = self.outstr + "OTHER THERMODYNAMIC PROPERTIES AT P=0 \n"
            self.outstr = self.outstr + "===================================== \n"
            self.outstr = self.outstr + "\n"
            self.outstr = self.outstr + "    T(K) \t B0(GPa) \t     B0' \t B0''(GPa-1) \t   Bs(GPa) \t alpha(10^5/K) \n"
            self.outstr = self.outstr + " -------------------------------------------------------------------------------------------------- \n"
            for j in xrange(self.ntemperature):
                self.outstr = self.outstr + '  ' + str(self.temperature[j]).rjust(6) + "\t" + str(self.b0p[j]).rjust(8)[:8] + "\t" + str(bp0p[j]).rjust(8)[:8] + "\t  " + str(bpp0p[j]).rjust(10)[:10] + "\t" + str(self.bs0p[j]).rjust(10)[:10] + "\t  " + str(self.a0p[j]*100000.0).rjust(12)[:12] + "\n"
            if self.ieos == 5 or self.ieos == 6:
                self.outstr = self.outstr + "\n"
                self.outstr = self.outstr + "SPINODAL PROPERTIES \n"
                self.outstr = self.outstr + "=================== \n"
                self.outstr = self.outstr + "\n"
                self.outstr = self.outstr + "    T(K) \t     K* \t  Psp(GPa) \t Vsp(bohr^3) \t    beta \n"
                self.outstr = self.outstr + " ------------------------------------------------------------------------ \n"
                for j in xrange(self.ntemperature):
                    self.outstr = self.outstr + '  ' + str(self.temperature[j]).rjust(6) + "\t" + str(bcnt_xk_sp[j]).rjust(10)[:10] + "\t" + str(-bcnt_xp_sp[j]).rjust(10)[:10] + "\t  " + str(bcnt_xv_sp[j]).rjust(10)[:10] + "\t  " + str(bcnt_xg_sp[j]).rjust(6)[:6] + "\n"
            self.outstr = self.outstr + "\n"
            self.outstr = self.outstr + "DEBYE MODEL RELATED PROPERTIES AT P=0 \n"
            self.outstr = self.outstr + "===================================== \n"
            self.outstr = self.outstr + "\n"
            self.outstr = self.outstr + "    T(K) \t Theta(K) \t gamma \n"
            self.outstr = self.outstr + " -------------------------------------- \n"
            for j in xrange(self.ntemperature):
                self.outstr = self.outstr + '  ' + str(self.temperature[j]).rjust(6) + "\t" + str(self.tdeb0p[j]).rjust(9)[:9] + "\t" + str(self.ga0p[j]).rjust(6)[:6] + "\n"
        #
        #.....end of the program
        #
        self.logstr = self.logstr + "MP Eqn of State Thermal: End of run ok! \n"
        return

 



    # **************************************************************************************
    #  These functions calculate the thermal properties using different equations of state
    # **************************************************************************************


    def numerical_eos(self, volref, nepol, epol, nfpol, fpol, statcalc):
        """
        .....numerical_eos - numerical EOS calculation.
        
        .....Numerical_Eos computes the derivatives of the Helmholtz function and the
        static energy needed to obtain Debye's temperature, the static
        pressure, and succesive derivatives of the bulk modulus.
        
        Adapted from original Fortran version written by M. A. Blanco et al.
        See Computer Physics Communications 158, 57-72 (2004) and Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details
        """
        #
        #.....Compute Pfit(P), B(P), B'(P), and B''(P)
        #
        if self.ieos >= 0:
            self.outstr = self.outstr + '\n'
            self.outstr = self.outstr + 'NUMERICAL EOS PRESSURE DERIVATIVES \n'
            self.outstr = self.outstr + '================================== \n'
            self.outstr = self.outstr + "  P(GPa) \t V(bohr^3) \t    V/V0 \t Pfit(GPa) \t    B(GPa) \t    B' \t  B''(GPa-1) \n"
            self.outstr = self.outstr + ' ------------------------------------------------------------------------------------------------------ \n'
        eos_polynomial_inst = eos_polynomial()
        for k in xrange(self.npressure):
            xeqmin = (self.voleqmin[k]/volref)**self.third
            dFdx = eos_polynomial_inst.polin1(xeqmin, nfpol, fpol)
            d2Fdx2 = eos_polynomial_inst.polin2(xeqmin, nfpol, fpol)
            d3Fdx3 = eos_polynomial_inst.polin3(xeqmin, nfpol, fpol)
            d4Fdx4 = eos_polynomial_inst.polin4(xeqmin, nfpol, fpol)
            pt = -xeqmin * dFdx / (3.0*self.voleqmin[k]) * self.au2gpa
            tmp = 2.0 * dFdx - xeqmin * d2Fdx2
            self.bulkmod[k] = -xeqmin / (9.0*self.voleqmin[k]) * tmp * self.au2gpa
            tmp2 = (d2Fdx2 - xeqmin * d3Fdx3) / tmp
            b1 = self.third * (2.0 - xeqmin * tmp2)
            b2 = -self.voleqmin[k] * (tmp2*(1.0-xeqmin*tmp2) - xeqmin*xeqmin*d4Fdx4/tmp) / (self.au2gpa*tmp)
            if k == 0:
                self.bulkmod0pres = self.bulkmod[k]
                self.dbulkmoddp0pres = b1
                self.d2bulkmoddp20pres = b2
            if self.ieos >= 0:
                self.outstr = self.outstr + '  ' + str(self.pressure[k]).rjust(6) + '\t' + str(self.voleqmin[k]).rjust(10)[:10] + '\t' + str(self.voleqmin[k]/self.voleqmin[0]).rjust(8)[:8] + '\t' + str(pt).rjust(10)[:10] + '\t' + str(self.bulkmod[k]).rjust(10)[:10] + '\t' + str(b1).rjust(6)[:6] + '\t     ' + str(b2).rjust(7)[:7] + '\n'
        #
        #.....Static calculation: get second derivative of static energy
        #
        if statcalc:
            if self.ieos >= 0:
                self.outstr = self.outstr +'\n'
                self.outstr = self.outstr + 'INPUT AND FITTED VALUES OF THE LATTICE ENERGY \n'
                self.outstr = self.outstr + '============================================= \n'
                self.outstr = self.outstr + '\n'
                self.outstr = self.outstr + '   V(bohr^3)     E_inp(hartree)     E_fit(hartree) \n'
                self.outstr = self.outstr + ' --------------------------------------------------\n'
            for i in xrange(self.ndata):
                f0 = eos_polynomial_inst.polin0(self.xconfigvector[i], nepol, epol)
                dEdx = eos_polynomial_inst.polin1(self.xconfigvector[i], nepol, epol)
                d2Edx2 = eos_polynomial_inst.polin2(self.xconfigvector[i], nepol, epol)
                tmp = self.xconfigvector[i] * d2Edx2 - 2.0 * dEdx
                volumex3 = 3.0 * self.vol_inp[i]
                self.d2EnergydVolume2_static.append(tmp * self.xconfigvector[i] / (volumex3*volumex3))
                if self.ieos >= 0:
                    self.outstr = self.outstr + '  ' + str(self.vol_inp[i]).rjust(10)[:10] + '\t ' + str(self.energ_inp[i]).rjust(14)[:14] + '\t    ' + str(f0).rjust(14)[:14] + '\n'
            #
            #.....Dynamic calculation: get static pressure and second derivative of the energy
            #
        else:
            for k in xrange(self.npressure):
                xeqmin = (self.voleqmin[k]/volref)**self.third
                dEdx = eos_polynomial_inst.polin1(xeqmin, nepol, epol)
                d2Edx2 = eos_polynomial_inst.polin2(xeqmin, nepol, epol)
                d3Edx3 = eos_polynomial_inst.polin3(xeqmin, nepol, epol)
                d4Edx4 = eos_polynomial_inst.polin4(xeqmin, nepol, epol)
                volumex3 = 3.0 * self.voleqmin[k]
                self.pstatic[k] = -dEdx*self.au2gpa * xeqmin / volumex3
                tmp = xeqmin * d2Edx2 - 2.0 * dEdx
                self.d2EnergydVolume2_dynamic[k] = tmp * xeqmin / (volumex3*volumex3)
                tmp2 = d2Edx2 - xeqmin * d3Edx3
                self.gamma_G[k] = (1.0 + xeqmin * tmp2 / tmp)/6.0
        return



    

    def vinet_eos(self, vol0pres, gfe0pres, statcalc):
        """
        .....vinet_eos - computes Vinet EOS from (P,V) data.
        
        .....VINET_EOS computes the EOS from the (P,V) data. The EOS has the
        following expresion:
        log H = A + B(1-x)
        being H = Px**2/(3(1-x))
        A = log Bo
        B = 3/2((Bo)'-1)
        X = (V/Vo)**(1/3)
        Adapted from original Fortran version written by M. A. Blanco et al.
        See Computer Physics Communications 158, 57-72 (2004) and Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details
        """
        #
        #.....fit Log H vs. (1-x)
        #
        logh = [0.0 for k in range(self.npressure)]
        x = [0.0 for k in range(self.npressure)]
        dbdp = []
        d2bdp2 = []
        sumz = 0.0
        sumy = 0.0
        sumzy = 0.0
        sumz2 = 0.0
        sumy2 = 0.0
        n = 0
        x[0] = 1.0
        i = 1
        while i < self.npressure:
            x[i] = (self.voleqmin[i]/vol0pres)**self.third
            h = self.pressure[i]*x[i]*x[i]/(3.0*(1-x[i]))
            logh[i] = math.log(h)
            z = 1-x[i]
            n = n+1
            sumz = sumz + z
            sumy = sumy + logh[i]
            sumzy = sumzy + z*logh[i]
            sumz2 = sumz2 + z*z
            sumy2 = sumy2 + logh[i]*logh[i]
            i = i + 1
        lnb0 = (sumy*sumz2 - sumzy*sumz)/(n*sumz2 - sumz*sumz)
        A = (n*sumzy - sumz*sumy)/(n*sumz2 - sumz*sumz)
        raiz = math.sqrt((sumz2 - sumz*sumz/n)*(sumy2 - sumy*sumy/n))
        rfit = (sumzy-sumz*sumy/n)/raiz
        logh[0] = lnb0
        #
        #.....obtain B0, B0', B0''
        #
        self.bulkmod0pres = math.exp(lnb0)
        self.dbulkmoddp0pres = 2.0*A*self.third+1.0
        self.d2bulkmoddp20pres = -(2.0+A*(A+6.0))/(9.0*self.bulkmod0pres)
        #
        #.....save static values
        #
        if statcalc:
            self.g00k = gfe0pres
            self.volstatcalc0p = vol0pres
            self.b00k = self.bulkmod0pres/self.au2gpa
            self.A00k = A
        #
        #.....Compute Pfit(P), B(P), B'(P), and B''(P)
        #
        self.bulkmod[0] = self.bulkmod0pres
        dbdp.append(self.dbulkmoddp0pres)
        d2bdp2.append(self.d2bulkmoddp20pres)
        self.pfit[0] = 0.0
        i = 1
        while i < self.npressure:
            a1x = A * (1.0 - x[i])
            ax1 = A * x[i] + 1.0
            f0x = x[i] * (1.0-a1x) - 2.0
            d1fdx1 = ax1 - a1x
            d2fdx2 = 2.0 * A
            f1f0 = d1fdx1 / f0x
            f2f0 = d2fdx2 / f0x
            fnw = 1.0 - x[i] * f1f0
            x2inv = 1.0 / (x[i]*x[i])
            b0exp = self.bulkmod0pres * math.exp(a1x)
            self.bulkmod[i] = -b0exp * f0x * x2inv
            dbdp.append(self.third * (ax1+fnw))
            d2bdp2.append(x[i]/(9.0*self.bulkmod[i]) * (x[i]*f2f0 - A + f1f0*fnw))
            self.pfit[i] = 3.0 * (1.0-x[i]) * x2inv * b0exp
            i = i + 1
        #
        #.....output
        #
        self.outstr = self.outstr + "\n"
        self.outstr = self.outstr + "VINET EOS PRESSURE DERIVATIVES \n"
        self.outstr = self.outstr + "============================== \n"
        self.outstr = self.outstr + "\n"
        self.outstr = self.outstr + "  1-V/V0 \t Vinet-Func \t P(GPa) \t Pfit(GPa) \t    B(GPa) \t        B' \t  B''(GPa-1) \n"
        self.outstr = self.outstr + " ------------------------------------------------------------------------------------------------------------ \n"
        for i in xrange(self.npressure):
            self.outstr = self.outstr + '  ' + str(1.0-x[i]).rjust(6)[:6] + "\t " + str(logh[i]).rjust(10)[:10] + "\t " + str(self.pressure[i]).rjust(6)[:6] + "\t        " + str(self.pfit[i]).rjust(10)[:10] + "\t" + str(self.bulkmod[i]).rjust(10)[:10] + "\t" + str(dbdp[i]).rjust(10)[:10] + "\t  " + str(d2bdp2[i]).rjust(10)[:10] + "\n"
            self.outstr = self.outstr + "\n"
            self.outstr = self.outstr + "B0 = " + str(self.bulkmod0pres) + ", B0' = " + str(self.dbulkmoddp0pres) + ", B0'' = " + str(self.d2bulkmoddp20pres)  + " reg.coef = " + str(rfit) + "\n"
            self.outstr = self.outstr + "\n"

        #
        #.....Static calculation: get static energy and its second derivative
        #
        if statcalc:
            for i in xrange(self.ndata):
                xstatcalc = (self.vol_inp[i]/self.volstatcalc0p)**self.third
                a1x = self.A00k * (1.0 - xstatcalc)
                f0x = xstatcalc * (1.0-a1x) - 2.0
                b0exp = self.b00k * math.exp(a1x)
                self.ust.append(self.g00k + 9.0*self.volstatcalc0p/(self.A00k*self.A00k) * (b0exp*(a1x-1.0)+self.b00k))
                self.d2EnergydVolume2_static.append(-f0x / (xstatcalc*xstatcalc*self.vol_inp[i]) * b0exp)
            #
            #.......Print input and fitted values of the lattice energy.
            #
            self.outstr = self.outstr + "\n"
            self.outstr = self.outstr + "INPUT AND FITTED VALUES OF THE LATTICE ENERGY \n"
            self.outstr = self.outstr + "============================================= \n"
            self.outstr = self.outstr + "\n"
            self.outstr = self.outstr + "   V(bohr^3)     E_inp(hartree)     E_fit(hartree) \n"
            self.outstr = self.outstr + " -------------------------------------------------- \n"
            for i in xrange(self.ndata):
                self.outstr = self.outstr + '  ' + str(self.vol_inp[i]).rjust(10)[:10] + "\t " + str(self.energ_inp[i]).rjust(14)[:14] + "\t    " + str(self.ust[i]).rjust(14)[:14] + "\n"
        #
        #.....Dynamic calculation: get static pressure and second derivative
        #     of the energy
        #
        else:
            for i in xrange(self.npressure):
                xstatcalc = (self.voleqmin[i]/self.volstatcalc0p)**self.third
                a1x = self.A00k * (1.0 - xstatcalc)
                ax1 = self.A00k * xstatcalc + 1.0
                f0x = xstatcalc * (1.0-a1x) - 2.0
                f1x = ax1 - a1x
                f2x = 2.0 * self.A00k
                f1f0 = f1x / f0x
                fnw = 1.0 - xstatcalc * f1f0
                x2inv = 1.0 / (xstatcalc*xstatcalc)
                b0exp = self.b00k * math.exp(a1x)
                self.pstatic[i] = 3.0 * (1.0-xstatcalc) * x2inv * b0exp * self.au2gpa
                self.d2EnergydVolume2_dynamic[i] = -f0x * x2inv * x2inv / (xstatcalc*self.volstatcalc0p) * b0exp
                self.gamma_G[i] = (ax1 + fnw - 1.0)/6.0
        return


    def birch_murnaghan_eos(self, vol0pres, gfe0pres, iG, statcalc):
        """
        .....birch_murnaghan_eos - computes the Birch-Murnaghan EOS of order iG from the
        (P,V) data.
        
        The EOS has the following expression:
        
        F = Sum (i=0,iG) a(i)*f^i
         
        being : F = P/[3f(1+2f)^(5/2)]
        f = [x^(-2)-1]/2
        x = [V(i)/V(1)]^(1/3)
    
        -----INPUT
        vol0pres      : Molecular volume (bohr^3/mol) at P=0.
        gfe0pres      : Gibbs energy (or 0k static energy) at v0 (hartree).
        iG      : order of the fitting.
        press() : Pressure values (GPa). common /eos/.
        vinp()  : Initial values of the volume (bohr^3/mol). common /input/.
        statcalc: Logical variable that determines if the calculation is
                static or dynamic. In the first case the second derivative
                of the static energy (d2EnergydVolume2_static) is computed for all the input
                values of the volume. In the second case the second
                derivative of the static energy (d2EnergydVolume2_dynamic) is computed for
                the equilibrium volumes at the different pressures.
    
        -----OUTPUT
        pstatic() : Static pressures in GPa (only on dynamic calculations).
        d2EnergydVolume2_static()  : Second derivative of ust(k) for each vinp(). Hy/bohr^6
        d2EnergydVolume2_dynamic()  : Second derivative of ust(k) for each V(). Hy/bohr^6
        rms     : Root mean square deviation.
        bulkmod0pres,dbulkmoddp0pres,d2bulkmoddp20pres : Bulk modulus and their derivatives at P=0.
        
        
        Adapted from original Fortran version written by M. A. Blanco et al.
        See Computer Physics Communications 158, 57-72 (2004) and Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details
        """
        tol = 1e-12
        npresm2 = self.npressure - 2
        izero = 0
        if iG > self.birchfitorder_iG:
#            self.logstr = self.logstr + "MP Eqn of State Thermal birch_murnaghan_eos : Too high fitting order \n"
            self.brerr = 1
            raise ValueError("MP Eqn of State Thermal: Fitting order for Birch-Murnaghan EoS is too high")
            return
        if math.fabs(self.pressure[0]) > tol:
#            self.logstr = self.logstr + "MP Eqn of State Thermal birch_murnaghan_eos : P(0) must be 0.0 \n"
            self.brerr = 1
            raise ValueError("MP Eqn of State Thermal: First pressure value must be 0.0")
            return
        birchcoef = [0.0 for i in range(self.birchfitorder_iG+1)]
        fbirch = []
        ybirch = []
        weight = []
        dbdp = [0.0 for i in range(self.npressure)]
        d2bdp2 = [0.0 for i in range(self.npressure)]
        #
        # .....Compute the Birch function F and strain variable f.
        #
        weight.append(0.0)
        i = 1
        while i < self.npressure:
            x0birch = (self.voleqmin[i]/vol0pres)**self.third
            fbirch.append((x0birch**(-2)-1)/2.0)
            ybirch.append(self.pressure[i]/self.au2gpa/(3*fbirch[i-1]*((1+2*fbirch[i-1]**2.5))))
            weight.append(1.0)
            i = i + 1
        #
        # .....Fitting to a polynomial of order iG.
        #
        eos_polynomial_inst = eos_polynomial()
        rms, birchcoef = eos_polynomial_inst.polfit(izero, npresm2, fbirch, ybirch, weight, iG)
        #
        # .....Compute B0,B0',B0''.
        #
        self.bulkmod0pres = birchcoef[0]*self.au2gpa
        if iG == 0:
            self.dbulkmoddp0pres = 4.0
            self.d2bulkmoddp20pres = -35.0/(9.0*self.bulkmod0pres)
        elif iG == 1:
            self.dbulkmoddp0pres = 4.0+2.0*birchcoef[1]*self.au2gpa/(3.0*self.bulkmod0pres)
            self.d2bulkmoddp20pres = (-self.dbulkmoddp0pres*(self.dbulkmoddp0pres-7.0)-143.0/9.0)/self.bulkmod0pres
        elif iG >= 2:
            self.dbulkmoddp0pres = 4.0+2.0*birchcoef[1]*self.au2gpa/(3.0*self.bulkmod0pres)
            self.d2bulkmoddp20pres = (2.0*birchcoef[2]/(self.bulkmod0pres*3.0)-self.dbulkmoddp0pres*(self.dbulkmoddp0pres-7.0)-143.0/9.0)/self.bulkmod0pres
        #
        # .....Compute B(P), B'(P), and B''(P). (b(), dbdp(), and d2bdp2().
        #
        for i in xrange(self.npressure):
            if i == 0:
                self.pfit[i] = 0.0
                self.bulkmod[i] = self.bulkmod0pres
                dbdp[i] = self.dbulkmoddp0pres
                d2bdp2[i] = self.d2bulkmoddp20pres
            else:
                st = fbirch[i-1]
                stsq = math.sqrt(1.0+2.0*st)
                s2 = stsq*stsq
                st32 = stsq*s2
                st52 = st32*s2
                pol0 = eos_polynomial_inst.polin0(st, iG, birchcoef)
                pol1 = eos_polynomial_inst.polin1(st, iG, birchcoef)
                pol2 = eos_polynomial_inst.polin2(st, iG, birchcoef)
                pol3 = 0.0
                if iG > 2:
                    pol3 = eos_polynomial_inst.polin3(st, iG, birchcoef)
                #
                # .........Fitted pressure and B(P).
                #
                self.pfit[i] = 3.0*st*st52*pol0*self.au2gpa
                sum1 = st32*(st*s2*pol1+(1.0+7.0*st)*pol0)
                self.bulkmod[i] = s2*sum1
                sum2 = st52*(s2*pol1+2*st*pol1+st*s2*pol2+7*pol0+(1+7*st)*pol1)
                den = 3*st*st52*pol1+(3.0*st52+15.0*st*st32)*pol0
                #
                # .........B'(P).
                #
                dbdp[i] = (5*sum1+sum2)/den
                d2bdf2 = 25*stsq*(st*s2*pol1+(1.0+7.0*st)*pol0)
                d2bdf2 = d2bdf2+10.0*st32*((2.0+11.0*st)*pol1+7*pol0+st*s2*pol2)
                d2bdf2 = d2bdf2+st52*((3.0+15.0*st)*pol2+18*pol1+st*s2*pol3)
                d2pdf2 = 3*st52*pol1+15*st*st32*pol1+3*st*st52*pol2
                d2pdf2 = d2pdf2+(30*st32+45*st*stsq)*pol0
                d2pdf2 = d2pdf2+(3*st52+15*st*st32)*pol1
                #
                # .........B''(P).
                #
                d2bdp2[i] = (den*d2bdf2-(5*sum1+sum2)*d2pdf2)/(den**3)
                self.bulkmod[i] = self.bulkmod[i]*self.au2gpa
                d2bdp2[i] = d2bdp2[i]/self.au2gpa
        #
        # .....Output.
        #
        self.outstr = self.outstr + "\n"
        self.outstr = self.outstr + "BIRCH-MURNAGHAN EOS PRESSURE DERIVATIVES \n"
        self.outstr = self.outstr + "======================================== \n"
        self.outstr = self.outstr + "\n"
        self.outstr = self.outstr + "  Strain \t Birch-Func \t P(GPa) \t Pfit(GPa) \t    B(GPa) \t        B' \t B''(GPa-1) \n"
        self.outstr = self.outstr + " ----------------------------------------------------------------------------------------------------------- \n"
        for i in xrange(self.npressure):
            if i == 0:
                self.outstr = self.outstr + '  ' + str(0.0).rjust(6)[:6] + "\t " + str(self.bulkmod0pres/self.au2gpa).rjust(10)[:10] + "\t " + str(self.pressure[i]).rjust(6)[:6] + "\t    " + str(self.pfit[i]).rjust(14)[:14] + "\t" + str(self.bulkmod[i]).rjust(10)[:10] + "\t" + str(dbdp[i]).rjust(10)[:10] + "\t " + str(d2bdp2[i]).rjust(10)[:10] + "\n"
            else:
                self.outstr = self.outstr + '  ' + str(fbirch[i-1]).rjust(6)[:6] + "\t " + str(ybirch[i-1]).rjust(10)[:10] + "\t " + str(self.pressure[i]).rjust(6)[:6] + "\t    " + str(self.pfit[i]).rjust(14)[:14] + "\t" + str(self.bulkmod[i]).rjust(10)[:10] + "\t" +  str(dbdp[i]).rjust(10)[:10] + "\t " +  str(d2bdp2[i]).rjust(10)[:10] + "\n"
        self.outstr = self.outstr + "\n"
        self.outstr = self.outstr + "B0 = " + str(self.bulkmod0pres) + ", B0' = " + str(self.dbulkmoddp0pres) + ", B0'' = " + str(self.d2bulkmoddp20pres) + ", reg.coef = " + str(rms) + "\n"
        self.outstr = self.outstr + "\n"
        if statcalc:
            #
            # .......Compute the static potential energy U(V) and its second
            #       derivative U''(V) with respect to V for all the input
            #       values of the volume.
            #
            self.volstatcalc0p = vol0pres
            self.g00k = gfe0pres
            for k in xrange(iG + 1):
                self.astatic[k] = birchcoef[k]
            for k in xrange(self.ndata):
                self.ust.append(self.g00k)
                self.d2EnergydVolume2_static.append(0.0)
                st = (self.vol_inp[k]/self.volstatcalc0p)**self.third
                st = ((st**(-2))-1)/2.0
                s2 = (1.0+2.0*st)
                pol0 = eos_polynomial_inst.polin0(st, iG, self.astatic)
                pol1 = eos_polynomial_inst.polin1(st, iG, self.astatic)
                v9 = 9.0*self.volstatcalc0p
                for j in xrange(iG + 1):
                    self.ust[k] = self.ust[k]+v9*self.astatic[j]/(j+2)*(st**(j+2))
                self.d2EnergydVolume2_static[k] = s2*s2*s2*s2/self.volstatcalc0p*(st*s2*pol1+(1.0+7.0*st)*pol0)
            #
            # .......Print input and fitted values of the lattice energy.
            #
            self.outstr = self.outstr + "\n"
            self.outstr = self.outstr + "INPUT AND FITTED VALUES OF THE LATTICE ENERGY \n"
            self.outstr = self.outstr + "============================================= \n"
            self.outstr = self.outstr + "\n"
            self.outstr = self.outstr + "   V(bohr^3)     E_inp(hartree)     E_fit(hartree) \n"
            self.outstr = self.outstr + " -------------------------------------------------- \n"
            for i in xrange(self.ndata):
                self.outstr = self.outstr + '  ' +  str(self.vol_inp[i]).rjust(10)[:10] + "\t " + str(self.energ_inp[i]).rjust(14)[:14] + "\t    " + str(self.ust[i]).rjust(14)[:14] + "\n"
            self.bmcoeffs_static = birchcoef
            return
        else:
            #
            # .......Compute the second derivative U''(V) with respect to V
            #       for all the equilibrium values of the volume at the
            #       different pressures.
            #
            for k in xrange(self.npressure):
                st = (self.voleqmin[k]/self.volstatcalc0p)**self.third
                st = ((st)**(-2)-1)/2.0
                s2 = (1.0+2.0*st)
                s22 = s2*s2
                pol0 = eos_polynomial_inst.polin0(st, iG, self.astatic)
                pol1 = eos_polynomial_inst.polin1(st, iG, self.astatic)
                pol2 = eos_polynomial_inst.polin2(st, iG, self.astatic)
                pol3 = 0.0
                if iG > 2:
                    pol3 = polin3(st, iG, self.astatic)
                self.pstatic[k] = self.au2gpa*3.0*st*(s2**2.5)*pol0
                tmp = (1.0+7.0*st)*pol0 + st*s2*pol1
                self.d2EnergydVolume2_dynamic[k] = s22*s22 / self.volstatcalc0p * tmp
                tmp = 1.0 / tmp
                tmp2 = s2*tmp * (7.0*pol0 + (2.0+11.0*st)*pol1 + st*s2*pol2)
                volumex3 = self.voleqmin[k] / (3.0*self.volstatcalc0p)
                self.gamma_G[k] = -2.0*self.third + 0.5*s2*math.sqrt(s2)*volumex3*(8.0+tmp2)
#            self.bmcoeffs_dynamic.append(birchcoef)
        #
        # .....end
        #
        return



    def bcn_eos(vol0pres, gfe0pres, b0, statcalc):
        """
        .....bcn_eos - compute the Spinodal (BCN) EOS from (B,p) data.
        
        The EOS has the following expresion:
        
                           g
         B(p) = ( p - Psp ) / K
         
         where//c
         g = 0.85  (If opt_g = .true. ===> g is optimized)
         (-Psp) and K are the parameter to optimize.
         
         These parameters bear the following relation with Bo and Bo'.
                     g  -1
         Bo  = (-Psp)  K
     
                          -1
         Bo' = g Bo (-Psp)
     
         -----Input parameters:
         lg     : Logical unit for results output.
         vol0pres     : Zero pressure volume, either static or dynamic.
         gfe0pres     : Zero pressure Gibbs function.
         B0     : Bulk modulus used to compute the initial value of
                   -Psp (GPa).
         opt_g  : if .true.  ==> g is optimized.
         static : if .true.  ==> static calculation.
     
         Adapted from original Fortran version written by M. A. Blanco et al.
         See Computer Physics Communications 158, 57-72 (2004) and Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details
         """
        eps = 1e-10
        maxnl = 100
        tol = 1e-12
        self.volp0 = vol0pres
        dbdp = []
        d2bdp2 = []
        xg = [0.0 for i in range(maxnl)]
        wg = [0.0 for i in range(maxnl)]
        #
        # .....Initial values of properties to optimize.
        #
        self.gbao = 0.85
        if not statcalc:
            if self.iopt_g == 2:
                self.gbao = self.gbao0
        x_Psp = self.gbao*b0/4.0
        #
        # .....Optimize g and x_Psp.
        #
        ax = x_Psp*0.5
        bx = x_Psp
        cx = x_Psp*2.0
        ax, bx, cx, fa, fb, fc = self.mnbrak(ax, bx, cx)
        x_Psp, desv = self.brent(ax, bx, cx, tol, x_Psp)
        #
        # .....Final properties.
        #
        xx = math.exp((self.gbao-1)*math.log(x_Psp))
        self.xmopt = self.xkopt/xx/(1.0-self.gbao)
        self.bulkmod0pres = xx * x_Psp / self.xkopt
        self.dbulkmoddp0pres = self.gbao * self.bulkmod0pres / x_Psp
        self.d2bulkmoddp20pres = self.gbao * (self.gbao - 1.0) * self.bulkmod0pres / x_Psp / x_Psp
        vsp = vol0pres * math.exp(self.gbao/(1-self.gbao)/self.dbulkmoddp0pres)
        self.pspin = x_Psp
        self.xsupa = self.xkopt
        self.vspin = vsp
        self.beta = self.gbao
        #
        # .....save static values
        #
        if statcalc:
            self.g00k = gfe0pres
            self.b00k = self.bulkmod0pres/self.au2gpa
            self.volstatcalc0p = vol0pres
            self.vsp0k = vsp
            self.xkopt0 = self.xkopt
            self.xmopt0 = self.xmopt
            self.x_Psp0 = x_Psp
            self.gbao0 = self.gbao
        #
        # .....Compute Pfit(P), B(P), B'(P), and B''(P)
        #
        for i in xrange(self.npressure):
            xxx = (self.xmopt + math.log(vol0pres/self.voleqmin[i]))/self.xmopt
            ug = 1.0/(1-self.gbao)
            self.pfit[i] = x_Psp * (math.exp(ug*math.log(xxx)) - 1.0)
            xdu = math.exp((self.gbao-1)*math.log(self.pressure[i]+x_Psp))
            self.bulkmod[i] = xdu * (self.pressure[i]+x_Psp) / self.xkopt
            dbdp.append(self.gbao * xdu / self.xkopt)
            d2bdp2.append(self.gbao * (self.gbao - 1) * xdu / self.xkopt / (self.pressure[i]+x_Psp))
        #
        # .....output
        #
        self.outstr = self.outstr + "\n"
        self.outstr = self.outstr + "SPINODAL EOS PRESSURE DERIVATIVES \n"
        self.outstr = self.outstr + "================================= \n"
        self.outstr = self.outstr + "\n"
        self.outstr = self.outstr + "  P(GPa) \t Pfit(GPa) \t    B(GPa) \t      B' \t B''(GPa-1) \n"
        self.outstr = self.outstr + " --------------------------------------------------------------------------- \n"
        for i in xrange(self.npressure):
            self.outstr = self.outstr + '  ' + str(self.pressure[i]).rjust(6)[:6] + "\t" + str(self.pfit[i]).rjust(10)[:10] + "\t" + str(self.bulkmod[i]).rjust(10)[:10] +  "\t" + str(dbdp[i]).rjust(8)[:8] + "\t " + str(d2bdp2[i]).rjust(10)[:10] + "\n"
        self.outstr = self.outstr + "\n"
        if self.opt_g:
            self.outstr = self.outstr + "B0 = " + str(self.bulkmod0pres) + ", B0' = " + str(dbdp[0]) + ", B0'' = " + str(d2bdp2[0]) + ", reg.coef = " + str(desv) + ", Vsp = " + str(vsp) + ", -Psp = " + str(x_Psp) + ", K* = " + str(self.xkopt) + ", gamma = " + str(self.gbao) + "\n"
        else:
            self.outstr = self.outstr + "B0 = " + str(self.bulkmod[0]) + ", B0' = " + str(dbdp[0]) + ", B0'' = " + str(d2bdp2[0]) + ", reg.coef = " + str(desv) + ", Vsp = " + str(vsp) + ", -Psp = " + str(x_Psp) + ", K* = " + str(self.xkopt) + ", gamma = " + str(self.gbao) + "\n"
        #
        # .....Static calculation: get static energy and its second derivative.
        #
        if statcalc:
            for i in xrange(self.ndata):
                x = (self.xmopt0+math.log(self.volstatcalc0p/self.vol_inp[i]))/self.xmopt0
                auxg = 1.0/(1.0-self.gbao0)
                auxg2 = self.gbao0/(1.0-self.gbao0)
                #
                # .........Compute numerically the integrated Helmholtz function by means
                #         of a loop with increasing number of Legendre points.
                #
                xinf = 1.0
                xsup = x
                factor = 1.0
                if xsup < xinf:
                    aux = xinf
                    xinf = xsup
                    xsup = aux
                    factor = -1.0
                #
                # .........Iterative loop.
                #
                sum0 = 1e30
                nl = 5
                xabs = 1.0
                while (nl <= maxnl) and (xabs >= eps):
                    gauleg(xinf, xsup, xg, wg, nl, self)
                    sum = 0.0
                    for ii in xrange(nl):
                        term = math.exp(self.xmopt0*(1.0-xg[ii])) * (math.exp(auxg*math.log(xg[ii])) - 1.0)
                        sum = sum + wg[ii] * factor * term * self.xmopt0 * self.volstatcalc0p * self.x_Psp0
                    xabs = math.fabs(sum-sum0)
                    sum0 = sum
                    nl = nl + 5
                self.ust.append(self.g00k + sum / self.au2gpa)
                self.d2EnergydVolume2_static.append(self.x_Psp0 / (self.xmopt0 * self.vol_inp[i] * (1.0 - self.gbao0)) * math.exp(auxg2*math.log(x)) / self.au2gpa)
            #
            # .......Print input and fitted values of the lattice energy.
            #
            self.outstr = self.outstr + "\n"
            self.outstr = self.outstr + "INPUT AND FITTED VALUES OF THE LATTICE ENERGY \n"
            self.outstr = self.outstr + "============================================= \n"
            self.outstr = self.outstr + "\n"
            self.outstr = self.outstr + "   V(bohr^3)     E_inp(hartree)     E_fit(hartree) \n"
            self.outstr = self.outstr + " -------------------------------------------------- \n"
            for i in xrange(self.ndata):
                self.outstr = self.outstr + '  ' + str(self.vol_inp[i]).rjust(10)[:10] + "\t " + str(self.energ_inp[i]).rjust(14)[:14] + "\t    " + str(self.ust[i]).rjust(14)[:14] + "\n"
        #
        # .....Dynamic calculation: get static pressure and second derivative
        #     of the energy
        #
        else:
            for i in xrange(self.npressure):
                xxx = (self.xmopt0 + math.log(self.volstatcalc0p/self.voleqmin[i]))/self.xmopt0
                ug = 1.0/(1.0-self.gbao0)
                auxg2 = self.gbao0/(1.0-self.gbao0)
                self.pstatic[i] = self.x_Psp0 * (math.exp(ug*math.log(xxx)) - 1.0)
                self.d2EnergydVolume2_dynamic[i] = self.x_Psp0 / (self.xmopt0 * self.volstatcalc0p * (1.0 - self.gbao0)) * math.exp(-self.xmopt0*(1.0-xxx)) * math.exp(auxg2*math.log(xxx)) / self.au2gpa
                self.gamma_G[i] = -1.0/6.0 + self.gbao0 * ug / 2.0/ self.xmopt0 / xxx
        #
        # .....end
        #
        return





    # **************************************************************************************
    #  This set of functions implement routines required for the BCNT EOS
    # **************************************************************************************

    def mnbrak(self, ax, bx, cx):
        """
        .....mnbrak - brackets a minimum of the function f.
        
        Given a function, and two distinct initial points ax and bx,
        this routine searches in the downhill direction (defined by the
        function as evaluated at the initial points) and returns new
        points ax, bx, and cx which bracket a minimum of the function.
        Also returned are the function values at the three points: fa, fb,
        and fc.
        
        Adapted from original Fortran version written by M. A. Blanco et al.
        See Computer Physics Communications 158, 57-72 (2004) and Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details
        """
        gold = 1.618034
        glimit = 100.0
        tiny = 1e-20
        fa = optm(ax, self)
        fb = optm(bx, self)
        if fb > fa:
            dum = ax
            ax = bx
            bx = dum
            dum = fb
            fb = fa
            fa = dum
        cx = bx+gold*(bx-ax)
        fc = optm(cx, self)
        while fb >= fc:
            endpart = True
            r = (bx-ax)*(fb-fc)
            q = (bx-cx)*(fb-fa)
            u = bx-((bx-cx)*q-(bx-ax)*r)/(2.0*math.copysign(max(math.fabs(q-r), tiny), q-r))
            ulim = bx+glimit*(cx-bx)
            if (bx-u)*(u-cx) > 0.0:
                fu = optm(u, self)
                if fu < fc:
                    ax = bx
                    fa = fb
                    bx = u
                    fb = fu
                    endpart = False
                elif fu > fb:
                    cx = u
                    fc = fu
                    endpart = False
                else:
                    u = cx+gold*(cx-bx)
                    fu = optm(u, self)
            elif (cx-u)*(u-ulim) > 0.0:
                fu = optm(u, self)
                if fu < fc:
                    bx = cx
                    cx = u
                    u = cx+gold*(cx-bx)
                    fb = fc
                    fc = fu
                    fu = optm(u, self)
            elif (u-ulim)*(ulim-cx) >= 0.0:
                u = ulim
                fu = optm(u, self)
            else:
                u = cx+gold*(cx-bx)
                fu = optm(u, self)
            if endpart:
                ax = bx
                bx = cx
                cx = u
                fa = fb
                fb = fc
                fc = fu
        return ax, bx, cx, fa, fb, fc


    def brent(self, ax, bx, cx, tol, xmin):
        """
        .....brent - unidimensional minimization of f in the range [ax,cx].
        
        Given a function, and a bracketing triplet of abscissas this
        routine isolates the minimum to a fractional precission of tol
        using Brent's method. The bracketing triplet must be such that bx
        is between ax and cx, and that f(bx) is less than both f(ax) and
        f(cx). The abscissa of the minimum is returned as xmin, and the
        minimum function value as BRENT.
        
        Adapted from original Fortran version written by M. A. Blanco et al.
        See Computer Physics Communications 158, 57-72 (2004) and Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details
        """
        itmax = 100
        cgold = 0.3819660
        zeps = 1.0e-10
        tol3 = 1e-12
        notskip1 = True
        a = min(ax, cx)
        b = max(ax, cx)
        v = bx
        w = v
        x = v
        e = 0.0
        d = 0.0
        fx = self.optm(x)
        fv = fx
        fw = fx
        iter = 1
        while iter <= itmax:
            xm = 0.5*(a+b)
            tol1 = tol*math.fabs(x)+zeps
            tol2 = 2.0*tol1
            if math.fabs(x-xm) <= (tol2-0.5*(b-a)):
                xmin = x
                brentx = fx
                return xmin, brentx
            else:
                if math.fabs(e) > tol1:
                    r = (x-w)*(fx-fv)
                    q = (x-v)*(fx-fw)
                    p = (x-v)*q-(x-w)*r
                    q = 2.0*(q-r)
                    if q > 0.0:
                        p = -p
                        q = math.fabs(q)
                        etemp = e
                        e = d
                    if math.fabs(p) >= math.fabs(0.5*q*etemp) or p <= q*(a-x) or p >= q*(b-x):
                        notskip1 = True
                    else:
                        d = p/q
                        u = x+d
                        if u-a < tol2 or b-u < tol2:
                            d = math.copysign(tol1, xm-x)
                            if math.fabs(d) >= tol1:
                                u = x+d
                            else:
                                u = x + math.copysign(tol1, d)
                            notskip1 = False
                if notskip1:
                    if x >= xm:
                        e = a-x
                    else:
                        e = b-x
                    d = cgold*e
                if math.fabs(d) >= tol1:
                    u = x+d
                else:
                    u = x + math.copysign(tol1, d)
                fu = optm(u, self)
                if fu <= fx:
                    if u >= x:
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
                    if u < x:
                        a = u
                    else:
                        b = u
                if fu <= fw or math.fabs(w - x) < tol3:
                    v = w
                    fv = fw
                    w = u
                    fw = fu
                elif fu <= fv or math.fabs(v - x) < tol3 or math.fabs(v - w) < tol3:
                    v = u
                    fv = fu
            iter = iter + 1
        self.logstr = self.logstr + "MP Eqn of State Thermal brent: exceeded maximum iterations. \n"
        xmin = x
        brentx = fx
        return xmin, brentx



    def optm(self, x_Psp):
        """
        .....optm - optimization of exponent parameter "g" (aka "beta") required for BCNT EOS.
        
        Adapted from original Fortran version written by M. A. Blanco et al.
        See Computer Physics Communications 158, 57-72 (2004) and Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details
        """
        xfunc = [0.0 for k in range(self.npressure)]
        yfunc = [0.0 for k in range(self.npressure)]
        if x_Psp < 0.0:
            self.logstr = self.logstr + "MP Eqn of State Thermal optm: Warning: Spinodal pressure is negative \n"
            desv = 1e30
            return desv
        if self.opt_g:
            a11 = 0.0
            a12 = 0.0
            a21 = 0.0
            a22 = 0.0
            z1 = 0.0
            z2 = 0.0
            for i in xrange(self.npressure):
                xfunc[i] = math.log(self.pressure[i]+x_Psp)
                yfunc[i] = math.log(self.bulkmod[i])
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
            for i in xrange(self.npressure):
                desv = desv + (yfunc[i] - x_in - x_de * xfunc[i])**2
            self.xkopt = math.exp(-x_in)
            self.gbao = x_de
            return desv
        else:
            a12 = 0.0
            z1 = 0.0
            for i in xrange(self.npressure):
                xfunc[i] = math.log(self.pressure[i]+x_Psp)
                yfunc[i] = math.log(self.bulkmod[i])
                a12 = a12 + xfunc[i]
                z1 = z1 + yfunc[i]
            x_in = (z1-self.gbao*a12)/self.npressure
            self.xkopt = math.exp(-x_in)
            desv = 0.0
            for i in xrange(self.npressure):
                desv = desv + (yfunc[i] - x_in - self.gbao * xfunc[i])**2
        return desv





    # **************************************************************************************
    #  This set of functions determine the Debye temperature self-consistently
    # **************************************************************************************


    def scdebye(self, T, pol, err, npol, imin, firsttime, pst):
        """
        .....scdebye - calculates the table of self consistent Debye
        temperatures at T for all the input volumes.
        
        If firsttime is true, only calculates the static pressures table
        
        Adapted from original Fortran version written by M. A. Blanco et al.
        See Computer Physics Communications 158, 57-72 (2004) and Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details
        """
        eps = 5.0
        mloops = 25
        scerr = False
        #
        # .....static pressures?
        #
        eos_polynomial_inst = eos_polynomial()
        eos_thermal_functions_inst = eos_thermal_functions()
        if firsttime:
            for i in xrange(self.ndata):
                plnv = eos_polynomial_inst.polin1(self.xconfigvector[i], npol, pol)
                pst[i] = -self.xconfigvector[i] * plnv / (3.0*self.vol_inp[i])
            return scerr
        #
        # .....loop until convergence is achieved
        #
        converged = False
        iloops = 0
        itry = imin
        while ((not converged) or (iloops < mloops)) and (iloops < self.maxloops):
            iloops = iloops + 1
            for i in xrange(self.ndata):
                self.thermal_energ, self.thermal_entropy, self.thermal_helmholtz, self.thermal_cv, D, Derr = eos_thermal_functions_inst.thermal(self.tdebye[i], T, self.natoms, self.pckbau, self.maxloops, self.logstr)
                self.datatofit[i] = self.energ_inp[i] + self.thermal_helmholtz
            self.ierr, imin, self.logstr = eos_polynomial_inst.minbrack(imin, self.datatofit, self.ndata, self.logstr)
            if self.ierr != 0:
                self.logstr = self.logstr + "MP Eqn of State Thermal scdebye: T = "  + str(T) + ", minimum point = " + str(imin) + ", trial point = " + str(itry) + ", total points = " + str(self.ndata) + " \n"
                self.logstr = self.logstr + "MP Eqn of State Thermal scdebye: func = " + str(self.datatofit) + " \n"
                scerr = True
                return scerr
            self.fterr, npol, pol, err = eos_polynomial_inst.polynomial_fit(imin, self.xconfigvector, self.datatofit, self.ndata, self.mpar, self.mfit)
            if self.fterr != 0:
                self.logstr = self.logstr + "MP Eqn of State Thermal scdebye: No polynomial fit found with a minimum within bounds of input data \n"
                self.logstr = self.logstr + "MP Eqn of State Thermal scdebye: T = " + str(T) + ", minimum point = " + str(imin) + ", trial point = " + str(itry) + ", total points = " + str(self.ndata) + " \n"
                scerr = True
                return scerr
            #
            # .....calculate the new set of debye temperatures and its convergence
            #
            converged = True
            theta0 = self.tdebye[0]
            for i in xrange(self.ndata):
                self.thermal_energ, self.thermal_entropy, self.thermal_helmholtz, self.thermal_cv, D, Derr = eos_thermal_functions_inst.thermal(self.tdebye[i], T, self.natoms, self.pckbau, self.maxloops, self.logstr)
                U = self.thermal_energ
                Cvt = self.thermal_cv
                f1 = polin1(self.xconfigvector[i], npol, pol)
                f2 = polin2(self.xconfigvector[i], npol, pol)
                p = -self.xconfigvector[i] * f1 / (3.0*self.vol_inp[i])
                bt = -self.xconfigvector[i] / (9.0*self.vol_inp[i]) * (2.0*f1 - self.xconfigvector[i]*f2)
                gm = (p-pst[i]) * self.vol_inp[i] / U
                bsv = self.vol_inp[i] * bt + T * gm*gm * Cvt
                if bsv < 0.0:
                    if i <= imin:
                        self.logstr = self.logstr + "MP Eqn of State Thermal scdebye: Bs < 0 for equilibrium V! \n"
                        return scerr
                    if i > 0:
                        theta = self.tdebye[i]/theta0 * self.tdebye[i-1]
                    dt = 0.0
                else:
                    theta = ((6.0*pi*pi*self.natoms/self.vol_inp[i])**self.third) / self.pckbau * self.poissonratiofunction * math.sqrt(bsv/self.cellmassval)
                    dt = theta - self.tdebye[i]
                    if i > 0:
                        if theta > self.tdebye[i-1]:
                            theta = self.tdebye[i]/theta0 * self.tdebye[i-1]
                            self.logstr = self.logstr + "MP Eqn of State Thermal scdebye: Warning! gamma < 0, i = " + str(i) + ", T = " + str(T) + " \n"
                            dt = 0.0
                theta0 = self.tdebye[i]
                self.tdebye[i] = theta
                if self.vol_inp[i]-self.vol_inp[imin] < self.vol_inp[imin]*0.1 and i >= max(3, self.ndata/10):
                    converged = converged and (math.fabs(dt) < eps)
        #
        # ....warn of convergence failure
        #
        if iloops >= self.maxloops and not converged:
            self.logstr = self.logstr + "MP Eqn of State Thermal scdebye: Maximum number of convergence iterations exceeded \n"
            self.logstr = self.logstr + "MP Eqn of State Thermal scdebye: Number of convergence loops = " + str(iloops) + " \n"
        if not converged:
            self.logstr = self.logstr + "MP Eqn of State Thermal scdebye: Warning! convergence not achieved \n"
        #
        # .....end
        #
        return scerr


