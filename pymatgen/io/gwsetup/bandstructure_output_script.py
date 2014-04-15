"""
Script to test writing GW Input for VASP.
Reads the POSCAR_name in the the current folder and outputs GW input to
subfolders name
"""

from __future__ import division

__author__ = "Michiel van Setten"
__copyright__ = " "
__version__ = "0.9"
__maintainer__ = "Michiel van Setten"
__email__ = "mjvansetten@gmail.com"
__date__ = "Oct 23, 2013"

import os
import os.path
import pymatgen as pmg
import pymatgen.io.vaspio.vasp_output
from pymatgen.core import structure

from pymatgen.io.vaspio.vasp_input import Kpoints
from pymatgen.io.vaspio_set import DictVaspInputSet
from pymatgen.matproj.rest import MPRester
from pymatgen.io.vaspio.vasp_output import Vasprun

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))

"""
MPGWVaspInputSet.joson contains the standards for GW calculations. This set contains all
parameters for the first sc dft calculation. The modifications for the subsequent
sub calculations are made below.
For many settings the number of cores on which the calculations will be run is needed, this
number is assumed to be on the environment variable NPARGWCALC.
"""

data = Vasprun('vasprun.xml', ionic_step_skip=1)
print data.converged
bandstructure = data.get_band_structure('../IBZKPT')
print 'gap: ', bandstructure.get_band_gap()['energy'], ' direct : ', bandstructure.get_band_gap()['direct']
print 'cbm: ', bandstructure.get_cbm()['energy'], data.actual_kpoints[bandstructure.get_cbm()['kpoint_index'][0]]
print 'vbm: ', bandstructure.get_vbm()['energy'], data.actual_kpoints[bandstructure.get_vbm()['kpoint_index'][0]]
print 'gap: ', bandstructure.get_cbm()['energy'] - bandstructure.get_vbm()['energy']
print 'direct gap: ', bandstructure.get_direct_band_gap()
