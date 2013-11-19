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

    def to_and_from_dict_test(self, qcinp):
        """
        Helper function. This function should be called in each specific test.
        """
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
   jobtype = sp
  exchange = b3lyp
     basis = 6-31+g*
$end

'''
        qcinp = QcInput(mol, title="Test Methane", exchange="B3LYP",
                        jobtype="SP",
                        basis_set="6-31+G*")
        self.assertEqual(str(qcinp), ans)
        self.to_and_from_dict_test(qcinp)

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
    jobtype = freq
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
                        jobtype="Freq",
                        basis_set={"C": "6-31G*", "h": "6-31g*",
                                   "CL": "6-31+g*"},
                        aux_basis_set={"c": "rimp2-cc-pvdz",
                                       "H": "rimp2-cc-pvdz",
                                       "Cl": "rimp2-aug-cc-pvdz"})
        self.assertEqual(str(qcinp), ans)
        self.to_and_from_dict_test(qcinp)

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
   jobtype = opt
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
                        jobtype="Opt",
                        basis_set={"Br": "srlc", "Cd": "srsc"},
                        ecp={"Br": "SrlC", "Cd": "srsc"})
        self.assertEqual(str(qcinp), ans)
        self.to_and_from_dict_test(qcinp)

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
     jobtype = sp
    exchange = b3lyp
       basis = 6-31+g*
  mem_static = 500
   mem_total = 18000
$end

'''
        qcinp = QcInput(mol, title="Test Methane", exchange="B3LYP",
                        jobtype="SP",
                        basis_set="6-31+G*")
        qcinp.set_memory(total=18000, static=500)
        self.assertEqual(str(qcinp), ans)
        self.to_and_from_dict_test(qcinp)

    def test_set_max_num_of_scratch_files(self):
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
           jobtype = sp
          exchange = b3lyp
             basis = 6-31+g*
  max_sub_file_num = 500
$end

'''
        qcinp = QcInput(mol, title="Test Methane", exchange="B3LYP",
                        jobtype="SP",
                        basis_set="6-31+G*")
        qcinp.set_max_num_of_scratch_files(500)
        self.assertEqual(str(qcinp), ans)
        self.to_and_from_dict_test(qcinp)

    def test_set_max_scf_iterations(self):
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
         jobtype = sp
        exchange = b3lyp
           basis = 6-31+g*
  max_scf_cycles = 100
   scf_algorithm = diis_gdm
$end

'''
        qcinp = QcInput(mol, title="Test Methane", exchange="B3LYP",
                        jobtype="SP",
                        basis_set="6-31+G*")
        qcinp.set_scf_algorithm_and_iterations(algorithm="diis_gdm",
                                               iterations=100)
        self.assertEqual(str(qcinp), ans)
        self.to_and_from_dict_test(qcinp)

    def test_set_scf_convergence_threshold(self):
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
          jobtype = sp
         exchange = b3lyp
            basis = 6-31+g*
  scf_convergence = 8
$end

'''
        qcinp = QcInput(mol, title="Test Methane", exchange="B3LYP",
                        jobtype="SP",
                        basis_set="6-31+G*")
        qcinp.set_scf_convergence_threshold(exponent=8)
        self.assertEqual(str(qcinp), ans)
        self.to_and_from_dict_test(qcinp)

    def test_set_integral_threshold(self):
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
   jobtype = sp
  exchange = b3lyp
     basis = 6-31+g*
    thresh = 14
$end

'''
        qcinp = QcInput(mol, title="Test Methane", exchange="B3LYP",
                        jobtype="SP",
                        basis_set="6-31+G*")
        qcinp.set_integral_threshold(thresh=14)
        self.assertEqual(str(qcinp), ans)
        self.to_and_from_dict_test(qcinp)

    def test_set_dft_grid(self):
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
   jobtype = sp
  exchange = b3lyp
     basis = 6-31+g*
   xc_grid = 000110000590
$end

'''
        qcinp = QcInput(mol, title="Test Methane", exchange="B3LYP",
                        jobtype="SP",
                        basis_set="6-31+G*")
        qcinp.set_dft_grid(radical_points=110, angular_points=590)
        self.assertEqual(str(qcinp), ans)
        self.to_and_from_dict_test(qcinp)

    def test_set_scf_initial_guess(self):
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
    jobtype = sp
   exchange = b3lyp
      basis = 6-31+g*
  scf_guess = gwh
$end

'''
        qcinp = QcInput(mol, title="Test Methane", exchange="B3LYP",
                        jobtype="SP",
                        basis_set="6-31+G*")
        qcinp.set_scf_initial_guess("GWH")
        self.assertEqual(str(qcinp), ans)
        self.to_and_from_dict_test(qcinp)

    def test_geom_opt_max_cycles(self):
        ans = '''$comments
 Test Methane
$end


$molecule
 1  2
 C           0.00000000        0.00000000        0.00000000
 H           0.00000000        0.00000000        1.08900000
 H           1.02671900        0.00000000       -0.36300000
 H          -0.51336000       -0.88916500       -0.36300000
 Cl         -0.51336000        0.88916500       -0.36300000
$end


$rem
              jobtype = sp
             exchange = b3lyp
                basis = 6-31+g*
  geom_opt_max_cycles = 100
$end

'''
        qcinp = QcInput(mol, title="Test Methane", exchange="B3LYP",
                        jobtype="SP", charge=1, spin_multiplicity=2,
                        basis_set="6-31+G*")
        qcinp.set_geom_max_iterations(100)
        self.assertEqual(str(qcinp), ans)
        self.to_and_from_dict_test(qcinp)

    def test_set_geom_opt_coords_type(self):
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
          jobtype = sp
         exchange = b3lyp
            basis = 6-31+g*
  geom_opt_coords = 0
$end

'''
        qcinp = QcInput(mol, title="Test Methane", exchange="B3LYP",
                        jobtype="SP",
                        basis_set="6-31+G*")
        qcinp.set_geom_opt_coords_type("cartesian")
        self.assertEqual(str(qcinp), ans)
        self.to_and_from_dict_test(qcinp)

    def test_scale_geom_opt_threshold(self):
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
                    jobtype = sp
                   exchange = b3lyp
                      basis = 6-31+g*
  geom_opt_tol_displacement = 120
        geom_opt_tol_energy = 10
      geom_opt_tol_gradient = 30
$end

'''
        qcinp = QcInput(mol, title="Test Methane", exchange="B3LYP",
                        jobtype="SP",
                        basis_set="6-31+G*")
        qcinp.scale_geom_opt_threshold(gradient=0.1, displacement=0.1,
                                       energy=0.1)
        self.assertEqual(str(qcinp), ans)
        self.to_and_from_dict_test(qcinp)

    def test_set_geom_opt_use_gdiis(self):
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
            jobtype = sp
           exchange = b3lyp
              basis = 6-31+g*
  geom_opt_max_diis = -1
$end

'''
        qcinp = QcInput(mol, title="Test Methane", exchange="B3LYP",
                        jobtype="SP",
                        basis_set="6-31+G*")
        qcinp.set_geom_opt_use_gdiis()
        self.assertEqual(str(qcinp), ans)
        self.to_and_from_dict_test(qcinp)

    def test_disable_symmetry(self):
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
     jobtype = sp
    exchange = b3lyp
       basis = 6-31+g*
  sym_ignore = True
    symmetry = False
$end

'''
        qcinp = QcInput(mol, title="Test Methane", exchange="B3LYP",
                        jobtype="SP",
                        basis_set="6-31+G*")
        qcinp.disable_symmetry()
        self.assertEqual(str(qcinp), ans)
        self.to_and_from_dict_test(qcinp)