import os
from unittest import TestCase
from pymatgen import Molecule
from pymatgen.io.qchemio import QcInput, QcBatchInput

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

    def elementary_io_verify(self, text, qcinp):
        self.to_and_from_dict_verify(qcinp)
        self.from_string_verify(contents=text, ref_dict=qcinp.to_dict)

    def to_and_from_dict_verify(self, qcinp):
        """
        Helper function. This function should be called in each specific test.
        """
        d1 = qcinp.to_dict
        qc2 = QcInput.from_dict(d1)
        d2 = qc2.to_dict
        self.assertEqual(d1, d2)

    def from_string_verify(self, contents, ref_dict):
        qcinp = QcInput.from_string(contents)
        d2 = qcinp.to_dict
        self.assertEqual(ref_dict, d2)

    def test_read_zmatrix(self):
        contents = '''$moLEcule

 1 2

 S
 C  1 1.726563
 H  2 1.085845 1 119.580615
 C  2 1.423404 1 114.230851 3 -180.000000 0
 H  4 1.084884 2 122.286346 1 -180.000000 0
 C  4 1.381259 2 112.717365 1 0.000000 0
 H  6 1.084731 4 127.143779 2 -180.000000 0
 C  6 1.415867 4 110.076147 2 0.000000 0
 F  8 1.292591 6 124.884374 4 -180.000000 0
$end

$reM
   BASIS  =  6-31+G*
   EXCHANGE  =  B3LYP
   jobtype  =  freq
$end

'''
        qcinp = QcInput.from_string(contents)
        ans = '''$molecule
 1  2
 S           0.00000000        0.00000000        0.00000000
 C           0.00000000        0.00000000        1.72656300
 H          -0.94431813        0.00000000        2.26258784
 C           1.29800105       -0.00000002        2.31074808
 H           1.45002821       -0.00000002        3.38492732
 C           2.30733813       -0.00000003        1.36781908
 H           3.37622632       -0.00000005        1.55253338
 C           1.75466906       -0.00000003        0.06427152
 F           2.44231414       -0.00000004       -1.03023099
$end


$rem
   jobtype = freq
  exchange = b3lyp
     basis = 6-31+g*
$end

'''
        ans_tokens = ans.split('\n')
        ans_text_part = ans_tokens[:2] + ans_tokens[11:]
        ans_coords_part = ans_tokens[2:11]
        converted_tokens = str(qcinp).split('\n')
        converted_text_part = converted_tokens[:2] + converted_tokens[11:]
        converted_coords_part = converted_tokens[2:11]
        self.assertEqual(ans_text_part, converted_text_part)
        for ans_coords, converted_coords in zip(ans_coords_part,
                                                converted_coords_part):
            ans_coords_tokens = ans_coords.split()
            converted_coords_tokens = converted_coords.split()
            self.assertEqual(ans_coords_tokens[0], converted_coords_tokens[0])
            xyz1 = ans_coords_tokens[1:]
            xyz2 = converted_coords_tokens[1:]
            for t1, t2 in zip(xyz1, xyz2):
                self.assertTrue(abs(float(t1)-float(t2)) < 0.0001)

    def test_no_mol(self):
        ans = '''$comments
 Test Methane
$end


$molecule
 -1  2
 read
$end


$rem
   jobtype = sp
  exchange = b3lyp
     basis = 6-31+g*
$end

'''
        qcinp = QcInput(molecule="READ", title="Test Methane",
                        exchange="B3LYP", jobtype="SP", charge=-1,
                        spin_multiplicity=2,
                        basis_set="6-31+G*")
        self.assertEqual(str(qcinp), ans)
        self.elementary_io_verify(ans, qcinp)

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
        self.elementary_io_verify(ans, qcinp)

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
        self.elementary_io_verify(ans, qcinp)

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
        self.elementary_io_verify(ans, qcinp)

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
        self.elementary_io_verify(ans, qcinp)

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
        self.elementary_io_verify(ans, qcinp)

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
        self.elementary_io_verify(ans, qcinp)

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
        self.elementary_io_verify(ans, qcinp)

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
        self.elementary_io_verify(ans, qcinp)

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
        self.elementary_io_verify(ans, qcinp)

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
        self.elementary_io_verify(ans, qcinp)

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
        self.elementary_io_verify(ans, qcinp)

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
        self.elementary_io_verify(ans, qcinp)

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
        self.elementary_io_verify(ans, qcinp)

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
        self.elementary_io_verify(ans, qcinp)

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
        self.elementary_io_verify(ans, qcinp)

    def test_use_cosmo(self):
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
  solvent_dielectric = 35.0
      solvent_method = cosmo
$end

'''
        qcinp = QcInput(mol, title="Test Methane", exchange="B3LYP",
                        jobtype="SP",
                        basis_set="6-31+G*")
        qcinp.use_cosmo(dielectric_constant=35.0)
        self.assertEqual(str(qcinp), ans)
        self.elementary_io_verify(ans, qcinp)


class TestQcBatchInput(TestCase):
    def test_str(self):
        ans = '''$comments
 Test Methane Opt
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
   jobtype = opt
  exchange = b3lyp
     basis = 6-31+g*
$end


@@@


$comments
 Test Methane Frequency
$end


$molecule
 read
$end


$rem
   jobtype = freq
  exchange = b3lyp
     basis = 6-31+g*
$end


@@@


$comments
 Test Methane Single Point Energy
$end


$molecule
 read
$end


$rem
   jobtype = sp
  exchange = b3lyp
     basis = 6-311+g(3df,2p)
$end

'''
        qcinp1 = QcInput(mol, title="Test Methane Opt", exchange="B3LYP",
                         jobtype="Opt", basis_set="6-31+G*")
        qcinp2 = QcInput(molecule="read", title="Test Methane Frequency",
                         exchange="B3LYP", jobtype="Freq", basis_set="6-31+G*")
        qcinp3 = QcInput(title="Test Methane Single Point Energy",
                         exchange="B3LYP", jobtype="SP",
                         basis_set="6-311+G(3df,2p)")
        qcbat = QcBatchInput(jobs=[qcinp1, qcinp2, qcinp3])
        self.assertEqual(str(qcbat), ans)

    def test_to_and_from_dict(self):
        qcinp1 = QcInput(mol, title="Test Methane Opt", exchange="B3LYP",
                         jobtype="Opt", basis_set="6-31+G*")
        qcinp2 = QcInput(molecule="read", title="Test Methane Frequency",
                         exchange="B3LYP", jobtype="Freq",
                         basis_set="6-31+G*")
        qcinp3 = QcInput(title="Test Methane Single Point Energy",
                         exchange="B3LYP", jobtype="SP",
                         basis_set="6-311+G(3df,2p)")
        qcbat1 = QcBatchInput(jobs=[qcinp1, qcinp2, qcinp3])
        d1 = qcbat1.to_dict
        qcbat2 = QcBatchInput.from_dict(d1)
        d2 = qcbat2.to_dict
        self.assertEqual(d1, d2)