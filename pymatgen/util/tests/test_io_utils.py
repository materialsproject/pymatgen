#!/usr/bin/env python

'''
Created on Nov 14, 2012
'''
from pymatgen.util.io_utils import read_backwards

__author__ = "Anubhav Jain"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Anubhav Jain"
__email__ = "ajain@lbl.gov"
__date__ = "Nov 14, 2012"

import unittest


class BackwardsReaderTest(unittest.TestCase):
    NUMLINES = 3000
    
    def test_read_backwards(self):
        """
        We are making sure a file containing line numbers is read in reverse order, i.e. the first line
        that is read corresponds to the last line number
        """
        with open("three_thousand_lines.txt") as f:
            for idx, line in enumerate(read_backwards(f)):
                self.assertEqual(int(line), self.NUMLINES - idx, "read_backwards read {} whereas it should have read {}".format(int(line), self.NUMLINES - idx))
        
            
if __name__ == "__main__":
    unittest.main()
