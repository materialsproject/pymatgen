#!/usr/bin/env python

"""
This module defines useful physical constants and conversion factors.
All units are in SI units except for conversion factors.

.. attribute:: ELECTRON_CHARGE

    Charge of an electron in coulombs.

.. attribute:: EPSILON_0

    Permittivity of vacuum

.. attribute:: BOLTZMANN_CONST

    Boltzmann's constant

.. attribute:: ELECTRON_VOLT

    eV in Joules.

.. attribute:: AVOGADROS_CONST

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

#Conversion factors

EV_PER_ATOM_TO_J_PER_MOL = ELECTRON_VOLT * AVOGADROS_CONST
EV_PER_ATOM_TO_KJ_PER_MOL = EV_PER_ATOM_TO_J_PER_MOL / 1000
ELECTRON_TO_AMPERE_HOURS = EV_PER_ATOM_TO_J_PER_MOL / 3600
AMU_TO_KG = 1.660538921e-27
