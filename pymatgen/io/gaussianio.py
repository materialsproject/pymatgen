# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
#!/usr/bin/env python

from __future__ import division, unicode_literals

"""
This module has been moved pymatgen.io.gaussian. This sub module will
be removed in pymatgen 4.0.
"""

__author__ = 'Shyue Ping Ong'
__copyright__ = 'Copyright 2013, The Materials Virtual Lab'
__version__ = '0.1'
__maintainer__ = 'Shyue Ping Ong'
__email__ = 'ongsp@ucsd.edu'
__date__ = '8/1/15'


import warnings
warnings.warn("pymatgen.io.gaussianio has been moved pymatgen.io.gaussian. "
              "This stub will be removed in pymatgen 4.0.")
from .gaussian import *
