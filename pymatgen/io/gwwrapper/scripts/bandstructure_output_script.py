# coding: utf-8

from __future__ import unicode_literals, division, print_function

"""
Script to parse bandstructure output
"""

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
from abipy.electrons.scissors import ScissorsBuilder

import abipy.data as abidata
from abipy.abilab import abiopen, ElectronBandsPlotter

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))


try:
    data = Vasprun('vasprun.xml', ionic_step_skip=1)
    print(data.converged)
    bandstructure = data.get_band_structure('../IBZKPT')
    print('gap: ', bandstructure.get_band_gap()['energy'], ' direct : ', bandstructure.get_band_gap()['direct'])
    print('cbm: ', bandstructure.get_cbm()['energy'], data.actual_kpoints[bandstructure.get_cbm()['kpoint_index'][0]])
    print('vbm: ', bandstructure.get_vbm()['energy'], data.actual_kpoints[bandstructure.get_vbm()['kpoint_index'][0]])
    print('gap: ', bandstructure.get_cbm()['energy'] - bandstructure.get_vbm()['energy'])
    print('direct gap: ', bandstructure.get_direct_band_gap())
except (IOError, OSError):
    print('no vasp output found')

try:
    # Get the quasiparticle results from the SIGRES.nc database.
    sigma_file = abiopen(abidata.ref_file("tgw1_9o_DS4_SIGRES.nc"))
    #sigma_file.plot_qps_vs_e0()
    qplist_spin = sigma_file.qplist_spin

    print(qplist_spin)

    # Construct the scissors operator
    domains = [[-10, 6.1], [6.1, 18]]
    scissors = qplist_spin[0].build_scissors(domains, bounds=None)

    # Read the KS band energies computed on the k-path
    ks_bands = abiopen(abidata.ref_file("si_nscf_GSR.nc")).ebands

    # Read the KS band energies computed on the Monkhorst-Pack (MP) mesh
    # and compute the DOS with the Gaussian method
    ks_mpbands = abiopen(abidata.ref_file("si_scf_GSR.nc")).ebands
    ks_dos = ks_mpbands.get_edos()

    # Apply the scissors operator first on the KS band structure
    # along the k-path then on the energies computed with the MP mesh.
    qp_bands = ks_bands.apply_scissors(scissors)

    qp_mpbands = ks_mpbands.apply_scissors(scissors)

    # Compute the DOS with the modified QPState energies.
    qp_dos = qp_mpbands.get_edos()

    # Plot the LDA and the QPState band structure with matplotlib.
    plotter = ElectronBandsPlotter()

    plotter.add_ebands("LDA", ks_bands, dos=ks_dos)

    plotter.add_ebands("LDA+scissors(e)", qp_bands, dos=qp_dos)

    # Define the mapping reduced_coordinates -> name of the k-point.
    klabels = {
        (0.5,  0.0,  0.0): "L",
        (0.0,  0.0,  0.0): "$\Gamma$",
        (0.0,  0.5,  0.5): "X",
    }

    plotter.plot(title="Silicon band structure", klabels=klabels)

except (IOError, OSError):
    print('no abinit output found')
