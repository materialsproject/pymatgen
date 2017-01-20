# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import warnings

import scipy.constants as constants

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
"""

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__status__ = "Production"
__date__ = "Sep 23, 2011"


warnings.warn("The pymatgen.core.physical_constants module is deprecated and "
              "will be removed in pymatgen 4.0. Pls use scipy.constants.",
              DeprecationWarning)

# Constants. Note that some of these may replicate functionality in
# scipy.constants. However, given the difficulty in installing scipy on many
# systems, the replication of these constants minimizes scipy dependency.

ELECTRON_CHARGE = constants.e
ELECTRON_MASS = constants.m_e
EPSILON_0 = constants.epsilon_0
BOLTZMANN_CONST = constants.k
ELECTRON_VOLT = constants.e
AVOGADROS_CONST = constants.N_A
HARTREE_TO_ELECTRON_VOLT = 1/constants.physical_constants[
    "electron volt-hartree relationship"][0]
SPEED_OF_LIGHT = constants.c
PLANCK_CONSTANT = constants.h

# Some useful aliases
N_a = AVOGADROS_CONST
k_b = BOLTZMANN_CONST
e = ELECTRON_CHARGE
R = AVOGADROS_CONST * BOLTZMANN_CONST
F = AVOGADROS_CONST * ELECTRON_CHARGE
c = SPEED_OF_LIGHT
h = PLANCK_CONSTANT
