import os
import unittest
import glob

from pymatgen.io.vaspio import Poscar
import pymatgen.command_line.platon_caller as platon

class TestPlatonCallerCase(unittest.TestCase):
    def test_get_space_group(self):
        dir_name = os.path.dirname(__file__)
        name = 'data'
        data_dir_name = os.path.join(dir_name, name, "POSCAR_*")
        for file_name in glob.glob(data_dir_name):
            number = int(file_name.split('_')[-1])
            poscar = Poscar.from_file(file_name)
            s = poscar.struct
            # SG_NB, SG_HM
            sg_dict =  platon.get_space_group(s)
            self.assertEqual(sg_dict['SG_NB'], number)

if __name__ == '__main__':
    unittest.main()
