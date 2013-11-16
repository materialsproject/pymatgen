import os
from unittest import TestCase
from pymatgen import Molecule
from pymatgen.io.qchemio import QcInput

__author__ = 'xiaohuiqu'


test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files', "molecules")


coords = [[0.000000, 0.000000, 0.000000],
          [0.000000, 0.000000, 1.089000],
          [1.026719, 0.000000, -0.363000],
          [-0.513360, -0.889165, -0.363000],
          [-0.513360, 0.889165, -0.363000]]
mol = Molecule(["C", "H", "H", "H", "Cl"], coords)


class TestQcInput(TestCase):
    def test_str(self):
        qcinp = QcInput(mol, title="Test Methane", exchange="xygjos",
                        job_type="Freq",
                        basis_set={"C": "6-31G*", "h": "6-31g*",
                                   "CL":"6-31+g*"},
                        aux_basis_set={"c": "rimp2-cc-pvdz",
                                       "H": "rimp2-cc-pvdz",
                                       "Cl": "rimp2-aug-cc-pvdz"})
        print qcinp
        print "#"*20