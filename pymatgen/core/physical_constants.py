#!/usr/bin/env python

"""
This module defines useful physical constants and conversion factors 
that may not be part of the scipy constants package. All units are in SI units
except for conversion factors.
"""

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__status__ = "Production"
__date__ = "Sep 23, 2011"


"""
Constants
"""
BOLTZMANN_CONST = 1.3806488e-23
ELECTRON_VOLT = 1.602176565e-19
AVOGADROS_CONST = 6.02214129e23

"""
Conversion factors
"""
EV_PER_ATOM_TO_J_PER_MOL = ELECTRON_VOLT * AVOGADROS_CONST
EV_PER_ATOM_TO_KJ_PER_MOL = EV_PER_ATOM_TO_J_PER_MOL / 1000
ELECTRON_TO_AMPERE_HOURS = EV_PER_ATOM_TO_J_PER_MOL / 3600
AMU_TO_KG = 1.660538921e-27



