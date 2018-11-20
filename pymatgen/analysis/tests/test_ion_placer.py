# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import copy
import os
import unittest

from pymatgen.core.structure import Molecule
from pymatgen.io.qchem.outputs import QCOutput
from pymatgen.util.testing import PymatgenTest
from pymatgen.analysis.ion_placer import IonPlacer


__author__ = "Samuel Blau"
__email__ = "samblau1@gmail.com"

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files', 'ion_placer_files')


# Still need to test placing different ions, placing a 2nd ion, and all five internal functions
class TestIonPlacer(PymatgenTest):

    def test_ec(self):
        ec_out = QCOutput(os.path.join(test_dir, "EC-12.qout")).data
        mol = ec_out['initial_molecule']
        mol.add_site_property("charge",ec_out['Mulliken'][0][::,0])
        mol.add_site_property("spin",ec_out['Mulliken'][0][::,1])
        tmp = IonPlacer(mol, "Li", 100000)
        save_points = copy.deepcopy(tmp.accepted_points)
        tmp.accepted_points = []
        self.assertEqual(len(save_points)!=0,True)
        for point in save_points:
            self.assertEqual(tmp._check_acceptance(point),True)
            tmp.accepted_points.append(point)

    def test_pc(self):
        ec_out = QCOutput(os.path.join(test_dir, "PC01.qout")).data
        mol = ec_out['initial_molecule']
        mol.add_site_property("charge",ec_out['Mulliken'][0])
        tmp = IonPlacer(mol, "Li", 100000)
        save_points = copy.deepcopy(tmp.accepted_points)
        tmp.accepted_points = []
        self.assertEqual(len(save_points)!=0,True)
        for point in save_points:
            self.assertEqual(tmp._check_acceptance(point),True)
            tmp.accepted_points.append(point)

if __name__ == "__main__":
    unittest.main()
