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
mol = Molecule(["C", "H", "H", "H", "H"], coords)


class TestQcInput(TestCase):
    def test_str(self):
        qcinp = QcInput(mol, title="Test Methane", exchange="XYGJOS",
                        job_type="freq")
        print qcinp
        print "#"*20