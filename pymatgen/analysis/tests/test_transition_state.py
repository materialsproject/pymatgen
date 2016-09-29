# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

import os
import unittest2 as unittest
from pymatgen.util.testing import PymatgenTest
from pymatgen.analysis.transition_state import NEBAnalysis
import json

"""
TODO: Modify unittest doc.
"""

__author__ = "David Waroquiers"
__copyright__ = "Copyright 2016, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyam Dwaraknath"
__email__ = "shyamd@lbl.gov"
__date__ = "2/5/16"


test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files', 'neb_analysis')


class NEBAnalysisTest(PymatgenTest):
    def runTest(self):
        neb_analysis = NEBAnalysis.from_dir(os.path.join(test_dir, 'neb'),
                                            relaxation_dirs=(os.path.join(test_dir, 'start'),
                                                             os.path.join(test_dir, 'end')))
        neb_analysis2 = NEBAnalysis.from_dict(neb_analysis.as_dict())
        json_data = json.dumps(neb_analysis.as_dict())

        neb_dict = json.loads(json_data)
        neb_analysis3 = NEBAnalysis.from_dict(neb_dict)

        self.assertArrayAlmostEqual(neb_analysis.r, neb_analysis2.r)
        self.assertArrayAlmostEqual(neb_analysis.energies, neb_analysis2.energies)
        self.assertArrayAlmostEqual(neb_analysis.forces, neb_analysis2.forces)
        self.assertEqual(neb_analysis.structures, neb_analysis2.structures)

        self.assertArrayAlmostEqual(neb_analysis.r, neb_analysis3.r)
        self.assertArrayAlmostEqual(neb_analysis.energies, neb_analysis3.energies)
        self.assertArrayAlmostEqual(neb_analysis.forces, neb_analysis3.forces)
        self.assertEqual(neb_analysis.structures, neb_analysis3.structures)

        self.assertAlmostEqual(neb_analysis.get_extrema()[1][0], (0.50023335723480078, 325.20043063935128))

if __name__ == '__main__':
    unittest.main()
