import unittest
import os

import pymatgen.io.vaspio 
from pymatgen.core.structure_modifier import OxidationStateDecorator
from pymatgen.analysis.ewald import EwaldSummation
from pymatgen.io.vaspio import Poscar
import numpy as np

class EwaldSummationTest(unittest.TestCase):

    def test_init(self):
        module_path = os.path.dirname(pymatgen.io.vaspio.__file__)
        filepath = os.path.join(module_path, 'tests','vasp_testfiles', 'POSCAR')
        p = Poscar.from_file(filepath)
        s = p.struct

        modifier = OxidationStateDecorator(s,{"Li":1, "Fe":2, "P":5, "O":-2})
        s = modifier.modified_structure
        ham = EwaldSummation(s)
        self.assertAlmostEqual(ham.real_space_energy, -354.91294268, 4, "Real space energy incorrect!")
        self.assertAlmostEqual(ham.reciprocal_space_energy, 25.475754801, 4, "Reciprocal space energy incorrect!")
        self.assertAlmostEqual(ham.point_energy, -790.463835033, 4, "Point space energy incorrect!")
        self.assertAlmostEqual(ham.total_energy, -1119.90102291, 2, "Total space energy incorrect!")
        self.assertAlmostEqual(sum(sum(abs(ham.forces))), 915.925354346, 4, "Forces incorrect")
        self.assertAlmostEqual(sum(sum(ham.real_space_energy_matrix)), -354.91294268, 4, "Real space energy matrix incorrect!")
        self.assertAlmostEqual(sum(sum(ham.reciprocal_space_energy_matrix)), 25.475754801, 4, "Reciprocal space energy matrix incorrect!")
        self.assertAlmostEqual(sum(ham.point_energy_matrix), -790.463835033, 4, "Point space energy matrix incorrect!")
        self.assertAlmostEqual(sum(sum(ham.total_energy_matrix)), -1119.90102291, 2, "Total space energy matrix incorrect!")
        #note that forces are not individually tested, but should work fine.
        

if __name__ == "__main__":
    unittest.main()
    
