#!/usr/bin/env python

"""
TODO: Change the module doc.
"""

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Jul 22, 2012"

import unittest
import os

from nose.exc import SkipTest

from pymatgen.command_line.bader_caller import BaderAnalysis
from pymatgen.util.io_utils import which


bader_present = which('bader')


class BaderAnalysisTest(unittest.TestCase):

    def test_init(self):
        if not bader_present:
            raise SkipTest("bader executable not present. Skipping...")
        test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                                'test_files')
        analysis = BaderAnalysis(os.path.join(test_dir, "CHGCAR.noncubic"))
        self.assertEqual(len(analysis.data), 8)
        self.assertAlmostEqual(analysis.data[0]["charge"], 7.4168)
        self.assertAlmostEqual(analysis.nelectrons, 52)
        self.assertAlmostEqual(analysis.vacuum_charge, 0)


if __name__ == '__main__':
    unittest.main()
