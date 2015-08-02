# coding: utf-8
#!/usr/bin/env python

from __future__ import division, unicode_literals

"""
#TODO: Write module doc.
"""

__author__ = 'Shyue Ping Ong'
__copyright__ = 'Copyright 2013, The Materials Virtual Lab'
__version__ = '0.1'
__maintainer__ = 'Shyue Ping Ong'
__email__ = 'ongsp@ucsd.edu'
__date__ = '8/1/15'


import warnings
warnings.warn("pymatgen.io.vaspio_set has been moved pymatgen.io.vasp.sets. "
              "This stub will be removed in pymatgen 4.0.")
from .vasp.sets import *
