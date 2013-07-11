#!/usr/bin/env python

"""
This module defines useful physical constants and conversion factors.
All units are in SI units except for conversion factors.

.. attribute:: ELECTRON_CHARGE or e

    Charge of an electron in coulombs.

.. attribute:: EPSILON_0

    Permittivity of vacuum

.. attribute:: BOLTZMANN_CONST or k_b

    Boltzmann's constant

.. attribute:: R

    Gas constant in J K-1 mol-1

.. attribute:: F

    Faraday's constant in C / mol

.. attribute:: ELECTRON_VOLT

    eV in Joules.

.. attribute:: AVOGADROS_CONST or N_a

    Avogardo's constant

The following are conversion factors.

.. attribute:: EV_PER_ATOM_TO_J_PER_MOL

    Conversion from ev/atom to J/mol

.. attribute:: EV_PER_ATOM_TO_KJ_PER_MOL

    Conversion from ev/atom to kJ/mol

.. attribute:: ELECTRON_TO_AMPERE_HOURS

    Conversion from electron charge to Amphere-hours

.. attribute:: AMU_TO_KG

    Conversion from atomic mass unit to kg
"""

import numpy as _np

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__status__ = "Production"
__date__ = "Sep 23, 2011"


#Constants. Note that some of these may replicate functionality in
#scipy.constants. However, given the difficulty in installing scipy on many
#systems, the replication of these constants minimizes scipy dependency.

ELECTRON_CHARGE = 1.602176565e-19
EPSILON_0 = 8.85418781762e-12
BOLTZMANN_CONST = 1.3806488e-23
ELECTRON_VOLT = 1.602176565e-19
AVOGADROS_CONST = 6.02214129e23

#Some useful aliases
N_a = AVOGADROS_CONST
k_b = BOLTZMANN_CONST
e = ELECTRON_CHARGE
R = AVOGADROS_CONST * BOLTZMANN_CONST
F = AVOGADROS_CONST * ELECTRON_CHARGE

#Conversion factors

EV_PER_ATOM_TO_J_PER_MOL = ELECTRON_VOLT * AVOGADROS_CONST
EV_PER_ATOM_TO_KJ_PER_MOL = EV_PER_ATOM_TO_J_PER_MOL / 1000
ELECTRON_TO_AMPERE_HOURS = EV_PER_ATOM_TO_J_PER_MOL / 3600
AMU_TO_KG = 1.660538921e-27


#: 1 Hartree, in eV
Ha_eV = 27.21138386

eV_Ha = 1./Ha_eV

#: 1 Bohr, in Angstrom
Bohr_Ang = 0.52917720859

Ang_Bohr = 1./Bohr_Ang

# Conversion factor eV/A**3 --> GPa
eVA3_GPa = 160.21773

###############################################################################
# Conversion tools.
###############################################################################


def Ha2eV(Ha):
    """
    Convert Hartree to eV

    Args:
        Ha:
            array_like with Hartree energies(s) to be converted.

    Returns:
        Array of floats with equivalent eV energies.

    >>> Ha2eV([1, 2])
    array([ 27.21138386,  54.42276772])
    """
    return _np.asanyarray(Ha) * Ha_eV


def Ha2meV(Ha):
    """
    Convert Hartree to meV
    """
    return Ha2eV(Ha) * 1000


def eV2Ha(eV):
    """
    Convert eV to Hartree

    Args:
        eV:
            array_like with eV energies to be converted.

    Returns:
        Array of floats with equivalent Hartree energies.

    >>> eV2Ha([ 27.21138386, 1])
    array([ 1.        ,  0.03674933])
    """
    return _np.asanyarray(eV) / Ha_eV


def Bohr2Ang(Bohr):
    """
    Convert Bohr to Angstrom

    Args:
        Bohr:
            array_like with Bohr lengths to be converted.

    Returns:
        array of floats with equivalent Angstrom lengths.

    >>> Bohr2Ang([1, 2])
    array([ 0.52917721,  1.05835442])
    """
    return _np.asanyarray(Bohr) * Bohr_Ang


def Ang2Bohr(Ang):
    """
    Convert Angstrom to Bohr.

    Args
        Ang:
            Array_like with Angstrom lengths to be converted.

    Returns:
        Array of floats with equivalent Bohr lengths.

    >>> Ang2Bohr(Bohr2Ang([1, 2]))
    array([ 1.,  2.])
    """
    return _np.asanyarray(Ang) / Bohr_Ang
