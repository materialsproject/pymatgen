# coding: utf-8

from __future__ import division, unicode_literals

__author__ = "Danny Broberg"
__copyright__ = "Copyright 2015, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Danny Broberg"
__email__ = "dbroberg@berkeley.edu"
__date__ = "Aug 11, 2015"

import unittest
import os

from pymatgen.command_line.sxdefectalign_caller import FreysoldtCorrection
from monty.os.path import which


@unittest.skipIf(not which('sxdefectalign'), "sxdefectalign executable not present.")
class sxdefectalignTest(unittest.TestCase):

    def test_init(self):
        test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                                'test_files')
        analysis = FreysoldtCorrection(os.path.join(test_dir, "LOCPOT_vref"),
                                 os.path.join(test_dir, "LOCPOT_vdef"),-1,5.765,520,[0.0,0.0,0.0],
                                 align=[0,0,0],lengths=[3.536141526,3.536141526,3.536141526])
        self.assert_Equal(analysis.errors['code'],None)
        full=analysis.get_full()
        for i in range(3):
            self.assertAlmostEqual(full[0][i],1.00207)
            self.assertAlmostEqual(full[1][i],1.4471626)
        s1=analysis.get_avgpotalign()
        s2=analysis.get_avgPCen()
        s3=analysis.get_singlecorrection()
        self.assertAlmostEqual(s1,1.4471626)
        self.assertAlmostEqual(s2,1.00207)
        self.assertAlmostEqual(s3,2.4492326)


if __name__ == '__main__':
    unittest.main()
