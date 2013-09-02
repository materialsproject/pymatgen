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
"""

import numpy as _np

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__status__ = "Production"
__date__ = "Sep 23, 2011"


#Constants. Note that some of these may replicate functionality in
#scipy.constants. However, given the difficulty in installing scipy on many
#systems, the replication of these constants minimizes scipy dependency.

ELECTRON_CHARGE = 1.602176565e-19
ELECTRON_MASS = 9.10938291e-31
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
RY_TO_EV = 13.605698066
BOHR_TO_ANGS = 0.5291772083

# 1 Hartree in eV
HA_TO_EV = 27.21138386
EV_TO_HA = 1./HA_TO_EV

# Conversion factor eV/A**3 --> GPa
EV_ANGS3_TO_GPA = 160.21773

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
    return _np.array(Ha) * HA_TO_EV


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
    return _np.array(eV) / HA_TO_EV


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
    return _np.array(Bohr) * BOHR_TO_ANGS


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
    return _np.array(Ang) / BOHR_TO_ANGS
