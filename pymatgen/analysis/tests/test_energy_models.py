# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

"""
TODO: Modify module doc.
"""


__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "11/19/13"

from pymatgen.analysis.energy_models import EwaldElectrostaticModel, \
    SymmetryModel, IsingModel
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure

import os
import unittest


test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')


class EwaldElectrostaticModelTest(unittest.TestCase):

    def test_get_energy(self):
        coords = [[0, 0, 0], [0.75, 0.75, 0.75], [0.5, 0.5, 0.5],
                  [0.25, 0.25, 0.25]]
        lattice = Lattice([[3.0, 0.0, 0.0],
                           [1.0, 3.0, 0.00],
                           [0.00, -2.0, 3.0]])
        s = Structure(lattice,
                      [{"Si4+": 0.5, "O2-": 0.25, "P5+": 0.25},
                       {"Si4+": 0.5, "O2-": 0.25, "P5+": 0.25},
                       {"Si4+": 0.5, "O2-": 0.25, "P5+": 0.25},
                       {"Si4+": 0.5, "O2-": 0.25, "P5+": 0.25}], coords)

        m = EwaldElectrostaticModel()
        self.assertAlmostEqual(m.get_energy(s), 44.1070954178)
        s2 = Structure.from_file(os.path.join(test_dir, "Li2O.cif"))
        self.assertAlmostEqual(m.get_energy(s2), -36.3476248117)

    def test_to_from_dict(self):
        m = EwaldElectrostaticModel()
        d = m.as_dict()
        self.assertIsInstance(EwaldElectrostaticModel.from_dict(d),
                              EwaldElectrostaticModel)

class SymmetryModelTest(unittest.TestCase):

    def test_get_energy(self):
        m = SymmetryModel()
        s2 = Structure.from_file(os.path.join(test_dir, "Li2O.cif"))
        self.assertAlmostEqual(m.get_energy(s2), -225)

    def test_to_from_dict(self):
        m = SymmetryModel(symprec=0.2)
        d = m.as_dict()
        o = SymmetryModel.from_dict(d)
        self.assertIsInstance(o, SymmetryModel)
        self.assertAlmostEqual(o.symprec, 0.2)


class IsingModelTest(unittest.TestCase):

    def test_get_energy(self):
        m = IsingModel(5, 6)
        from pymatgen.core.periodic_table import Specie
        s = Structure.from_file(os.path.join(test_dir, "LiFePO4.cif"))
        s.replace_species({"Fe": Specie("Fe", 2, {"spin": 4})})
        self.assertAlmostEqual(m.get_energy(s), 172.81260515787977)
        s[4] = Specie("Fe", 2, {"spin": -4})
        s[5] = Specie("Fe", 2, {"spin": -4})
        self.assertAlmostEqual(m.get_energy(s), 51.97424405382921)

    def test_to_from_dict(self):
        m = IsingModel(5, 4)
        d = m.as_dict()
        o = IsingModel.from_dict(d)
        self.assertIsInstance(o, IsingModel)
        self.assertAlmostEqual(o.j, 5)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
