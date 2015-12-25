# coding: utf-8
from __future__ import division, unicode_literals

"""
This module implements an interface to the Freysoldt et al.'s excellent
code for calculating a correction for formation energy of a charged 
defect.

This module depends on a compiled sxdefectalign executable available in 
the path. Please download the executable and manual at 
http://sxlib.mpie.de/wiki/AddOns.

If you use this module, please cite the following:

Christoph Freysoldt, Jörg Neugebauer, and Chris G. Van de Walle, 
Phys. Rev. Lett. 102, 016402, 2009.
"""

from six.moves import map
from six.moves import zip
__author__ = "Bharat Medasani"
__version__ = "0.1"
__maintainer__ = "Bharat Medasani"
__email__ = "mbkumar@gmail.com"
__status__ = "Beta"
__date__ = "8/5/15"

import os
import subprocess
import shutil
from pymatgen.io.vaspio.vasp_output import Chgcar
from pymatgen.io.vaspio.vasp_input import Potcar
from monty.os.path import which
from monty.dev import requires
from monty.tempfile import ScratchDir


@requires(which("sxdefecatalign"),
          "Charge defect correction requires the executable sxdefectalign to "
          "be in the path. Please download the executable at "
          "http://sxlib.mpie.de/wiki/AddOns.")
class FreysoldtCorrection(object):
    def __init__():
        pass

