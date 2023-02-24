import os
import unittest

from pymatgen.util.testing import PymatgenTest
from pymatgen.io.cif import CifParser

test_dir = os.path.join(os.path.dirname(__file__))

class Test_Skip_Checks(PymatgenTest):

    def setup(self):

        self.cif_list = [i for i in os.listdir(test_dir) if i.endswith('.cif')]
        self.structures = [CifParser(file, occupancy_tolerance=10000000).get_structures(primitive=False, symmetrized=True,skip_checks=True)[0] for file in self.cif_list]
    
    def test_skip_checks(self):

        self.assertEqual(self.structures[0][0].species.as_dict()['O'],1.36)


if __name__ == '__main__':
    unittest.main()
