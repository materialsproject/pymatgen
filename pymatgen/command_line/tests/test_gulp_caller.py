'''
Created on Apr 5, 2013

@author: wenhao
'''

import unittest

from nose.exc import SkipTest

from pymatgen.command_line.gulp_caller import get_binoxi_gulp_energy,Tersoff_pot,gulpduplicatecheck
from pymatgen.core.structure import Lattice, Structure
from pymatgen.core.periodic_table import Element
from pymatgen.util.io_utils import which
from pymatgen.io.vaspio.vasp_input import Poscar
import os

test_dir = os.path.join(os.path.dirname(__file__))
gulp_present = which('gulp')

class GulpCallerTest(unittest.TestCase):

    def setUp(self):
        
        if not gulp_present:
            raise SkipTest("aconvasp not present. Skipping...")
        self.struct=Poscar.from_file("POSCAR").structure

        
    def test_gulpduplicatestructure(self):
        print gulpduplicatecheck(self.struct)

if __name__ == '__main__':
    unittest.main()
