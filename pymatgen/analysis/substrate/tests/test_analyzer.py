# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

"""
TODO: Modify unittest doc.
"""

__author__ = "Shyam Dwaraknath"
__copyright__ = "Copyright 2016, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyam Dwaraknath"
__email__ = "shyamd@lbl.gov"
__date__ = "2/5/16"

import unittest
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.analysis.substrate.analyzer import ZurSuperLatticeGenerator as ZSLGen
from pymatgen.util.testing import PymatgenTest
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer



class ZSLGenTest(PymatgenTest):

    def test_ZSLGen(self):
        substrate = SpacegroupAnalyzer(self.get_structure("TiO2"),symprec=0.1).get_conventional_standard_structure()
        film = SpacegroupAnalyzer(self.get_structure("VO2"),symprec=0.1).get_conventional_standard_structure()

        z = ZSLGen(substrate,film)
        return z.generate()

    def runTest(self):

        self.test_ZSLGen()
        return True

if __name__ == '__main__':
    unittest.main()
