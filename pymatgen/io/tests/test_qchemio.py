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

coords2 = [[0.0, 0.0, -2.4],
          [0.0, 0.0, 0.0],
          [0.0, 0.0, 2.4]]
heavy_mol = Molecule(["Br", "Cd", "Br"], coords2)


class TestQcInput(TestCase):

    def test_to_and_from_dict(self):
        qcinp = QcInput(mol, title="Test Methane", exchange="xygjos",
                        job_type="Freq",
                        basis_set={"C": "6-31G*", "h": "6-31g*",
                                   "CL":"6-31+g*"},
                        aux_basis_set={"c": "rimp2-cc-pvdz",
                                       "H": "rimp2-cc-pvdz",
                                       "Cl": "rimp2-aug-cc-pvdz"},
                        ecp={"C":"srlc", "Cl":"srsc"})
        d1 = qcinp.to_dict
        qc2 = QcInput.from_dict(d1)
        d2 = qc2.to_dict
        self.assertEqual(d1, d2)


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


    def test_ecp_str(self):
        ans = '''$comments
 Test ECP
$end


$molecule
 0  1
 Br          0.00000000        0.00000000       -2.40000000
 Cd          0.00000000        0.00000000        0.00000000
 Br          0.00000000        0.00000000        2.40000000
$end


$rem
  job_type = opt
  exchange = b3lyp
     basis = gen
       ecp = gen
$end


$basis
 Br
 srlc
 ****
 Cd
 srsc
 ****
$end


$ecp
 Br
 srlc
 ****
 Cd
 srsc
 ****
$end

'''
        qcinp = QcInput(heavy_mol, title="Test ECP", exchange="B3LYP",
                        job_type="Opt",
                        basis_set={"Br":"srlc", "Cd": "srsc"},
                        ecp={"Br": "SrlC", "Cd": "srsc"})
        self.assertEqual(str(qcinp), ans)

    def test_set_memory(self):
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
  mem_static = 500
   mem_total = 18000
$end

'''
        qcinp = QcInput(mol, title="Test Methane", exchange="B3LYP",
                        job_type="SP",
                        basis_set="6-31+G*")
        qcinp.set_memory(total=18000, static=500)
        self.assertEqual(str(qcinp), ans)