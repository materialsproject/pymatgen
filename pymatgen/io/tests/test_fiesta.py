# coding: utf-8

from __future__ import division, unicode_literals

"""
Created on July 2015
"""


__author__ = "ndardenne"


import unittest2 as unittest
import os

from pymatgen import Molecule
from pymatgen.io.fiesta import FiestaInput, FiestaOutput

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", ".." , 'test_files')


class FiestaInputTest(unittest.TestCase):

    def setUp(self):

        coords = [[0.000000, 0.000000, 0.000000],
                  [0.000000, 0.000000, 1.089000],
                  [1.026719, 0.000000, -0.363000],
                  [-0.513360, -0.889165, -0.363000],
                  [-0.513360, 0.889165, -0.363000]]
        self.coords = coords
        mol = Molecule(["C", "H", "H", "H", "H"], coords)
        self.cellin = FiestaInput(mol, correlation_grid={'dE_grid': u'0.500', 'n_grid': u'14'},
                        Exc_DFT_option={'rdVxcpsi': u'1'},
                        COHSEX_options={'eigMethod': u'C', 'mix_cohsex': u'0.500', 'nc_cohsex': u'0', 'nit_cohsex': u'0',
                                        'nv_cohsex': u'0', 'resMethod': u'V', 'scf_cohsex_wf': u'0'},
                        GW_options={'nc_corr': u'10', 'nit_gw': u'3', 'nv_corr': u'10'},
                        BSE_TDDFT_options={'do_bse': u'1', 'do_tddft': u'0', 'nc_bse': u'382', 'nit_bse': u'50',
                                             'npsi_bse': u'1', 'nv_bse': u'21'})

    def test_init(self):
        mol = Molecule(["C", "H", "H", "H", "H"], self.coords)
        cellin = FiestaInput(mol)
        self.assertEqual(cellin.molecule.spin_multiplicity, 1)

    def test_str_and_from_string(self):
        ans = '# number of atoms and species\n   5    2\n# number of valence bands\n    5\n# number of points and spacing in eV for correlation grid\n    14    0.500\n# relire=1 ou recalculer=0 Exc DFT\n    1\n# number of COHSEX corrected occp and unoccp bands: C=COHSEX  H=HF\n    0    0   C\n# number of COHSEX iter, scf on wfns, mixing coeff; V=RI-V  I=RI-D\n    0   V       0       0.500\n# number of GW corrected occp and unoccp bands\n   10   10\n# number of GW iterations\n    3\n# dumping for BSE and TDDFT\n    1    0\n# number of occp. and virtual bands fo BSE: nocore and up to 40 eVs\n    21   382\n# number of excitations needed and number of iterations\n    1   50\n# list of symbols in order\nC\nH\n# scaling factor\n    1.000\n# atoms x,y,z cartesian .. will be multiplied by scale\n 0.0 0.0 0.0 1\n 0.0 0.0 1.089 2\n 1.026719 0.0 -0.363 2\n -0.51336 -0.889165 -0.363 2\n -0.51336 0.889165 -0.363 2\n            '
        self.assertEqual(str(self.cellin), ans)
        cellin = FiestaInput.from_string(ans)
        self.assertEqual(cellin.GW_options['nc_corr'], '10')
        self.assertEqual(cellin.COHSEX_options['eigMethod'], 'C')


class FiestaOutputTest(unittest.TestCase):


    def setUp(self):
        self.logfiesta = FiestaOutput(os.path.join(test_dir, "log_fiesta"))

    def test_props(self):
        out = self.logfiesta
        self.assertEqual(out.data[0]["Gaps"]["Egap_QP_Linear"], u'10.4135')
        self.assertEqual(out.data[0]["HOMO"], {'band': u'HOMO',
 'eKS': u'-7.3029',
 'eQP_Linear': u'-9.5142',
 'eQP_SCF': u'-8.9264',
 'eQP_old': u'-7.7188',
 'eXX': u'-15.9483',
 'sigma_c_Linear': u'-0.4587',
 'sigma_c_SCF': u'0.3900',
 'z': u'0.87'})



if __name__ == "__main__":
    unittest.main()
