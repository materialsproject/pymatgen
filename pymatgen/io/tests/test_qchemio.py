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

    def test_simple_basis_str(self):
        ans = '''$comments
 Test Methane
$end


$molecule
 0  1
 C           0.00000000        0.00000000        0.00000000
 H           0.00000000        0.00000000        1.08900000
 H           1.02671900        0.00000000       -0.36300000
 H          -0.51336000       -0.88916500       -0.36300000
 Cl         -0.51336000        0.88916500       -0.36300000
$end


$rem
  job_type = sp
  exchange = b3lyp
     basis = 6-31+g*
$end

'''
        qcinp = QcInput(mol, title="Test Methane", exchange="B3LYP",
                        job_type="SP",
                        basis_set="6-31+G*")
        self.assertEqual(str(qcinp), ans)


    def test_aux_basis_str(self):
        ans = '''$comments
 Test Methane
$end


$molecule
 0  1
 C           0.00000000        0.00000000        0.00000000
 H           0.00000000        0.00000000        1.08900000
 H           1.02671900        0.00000000       -0.36300000
 H          -0.51336000       -0.88916500       -0.36300000
 Cl         -0.51336000        0.88916500       -0.36300000
$end


$rem
   job_type = freq
   exchange = xygjos
      basis = gen
  aux_basis = gen
$end


$aux_basis
 C
 rimp2-cc-pvdz
 ****
 Cl
 rimp2-aug-cc-pvdz
 ****
 H
 rimp2-cc-pvdz
 ****
$end


$basis
 C
 6-31g*
 ****
 Cl
 6-31+g*
 ****
 H
 6-31g*
 ****
$end

'''
        qcinp = QcInput(mol, title="Test Methane", exchange="xygjos",
                        job_type="Freq",
                        basis_set={"C": "6-31G*", "h": "6-31g*",
                                   "CL":"6-31+g*"},
                        aux_basis_set={"c": "rimp2-cc-pvdz",
                                       "H": "rimp2-cc-pvdz",
                                       "Cl": "rimp2-aug-cc-pvdz"})
        self.assertEqual(str(qcinp), ans)

