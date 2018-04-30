# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

import logging
import os
import unittest
from pymatgen import Molecule
from pymatgen.util.testing import PymatgenTest
from pymatgen.io.qchem_io.inputs import QCInput

__author__ = "Brandon Wood, Samuel Blau, Shyam Dwaraknath, Julian Self"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "0.1"
__email__ = "b.wood@berkeley.edu"
__credits__ = "Xiaohui Qu"

logger = logging.getLogger(__name__)

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", "test_files", "qchem")


class TestQCInput(PymatgenTest):

    # ef setUpClass(cls):
    # add things that show up over and over again


    def test_molecule_template(self):
        species = ["C", "O"]
        coords = [[-9.5782000000, 0.6241500000, 0.0000000000], [-7.5827400000, 0.5127000000, -0.0000000000]]
        mol = Molecule(species=species, coords=coords)
        molecule_test = QCInput.molecule_template(mol)
        molecule_actual = """$molecule
 0 1
 C     -9.5782000000      0.6241500000      0.0000000000
 O     -7.5827400000      0.5127000000     -0.0000000000
$end"""

        self.assertEqual(molecule_actual, molecule_test)

    def test_rem_template(self):
        rem_params = {
            "job_type": "opt",
            "method": "wB97M-V",
            "basis": "def2-QZVPPD",
            "max_scf_cycles": 300,
            "gen_scfman": "true"
        }
        rem_test = QCInput.rem_template(rem_params)
        rem_actual = """$rem
   job_type = opt
   method = wB97M-V
   basis = def2-QZVPPD
   max_scf_cycles = 300
   gen_scfman = true
$end"""

        self.assertEqual(rem_actual, rem_test)

    def test_opt_template(self):
        opt_params = {
            "CONSTRAINT": ["tors 2 3 4 5 25.0", "bend 2 1 4 110.0"],
            "FIXED": ["x y 2 4 5"],
            "DUMMY": ["M 2 3 4 5"],
            "CONNECT": ["4 3 2 3 5 6"]
        }
        opt_test = QCInput.opt_template(opt_params)
        opt_actual = """$opt
CONSTRAINT
   tors 2 3 4 5 25.0
   bend 2 1 4 110.0
ENDCONSTRAINT

FIXED
   x y 2 4 5
ENDFIXED

DUMMY
   M 2 3 4 5
ENDDUMMY

CONNECT
   4 3 2 3 5 6
ENDCONNECT
$end"""
        self.assertEqual(opt_actual, opt_test)

    def test_pcm_template(self):
        pcm_params = {"theory": "cpcm"}
        pcm_test = QCInput.pcm_template(pcm_params)
        pcm_actual = """$pcm
   theory cpcm
$end"""
        self.assertEqual(pcm_actual, pcm_test)

    def test_solvent_template(self):
        solvent_params = {"dielectric": "5.0"}
        solvent_test = QCInput.solvent_template(solvent_params)
        solvent_actual = """$solvent
   dielectric 5.0
$end"""
        self.assertEqual(solvent_actual, solvent_test)

    def test_find_sections(self):
        str_single_job_input = """$molecule
 0  1
 S          -0.00250959       -0.05817469       -0.02921636
 C           1.70755408       -0.03033788       -0.01382912
 H           2.24317221       -0.05215019        0.92026728
 C           2.21976393        0.01718014       -1.27293235
 H           3.27786220        0.04082146       -1.48539646
 C           1.20867399        0.04478540       -2.27007793
 H           1.40292257        0.10591684       -3.33110912
 C          -0.05341046        0.01577217       -1.74839343
 C          -1.32843436        0.03545064       -2.45531187
 C          -1.55195156        0.08743920       -3.80184635
 H          -0.75245172        0.10267657       -4.52817967
 C          -2.93293778        0.08408786       -4.13352169
 H          -3.31125108        0.11340328       -5.14405819
 C          -3.73173288        0.02741365       -3.03412864
 H          -4.80776535        0.00535688       -2.99564645
 S          -2.81590978       -0.00516172       -1.58990580
$end


$rem
             job_type = opt
               method = wb97m-v
                basis = def2-tzvppd
           gen_scfman = true
  geom_opt_max_cycles = 75
       max_scf_cycles = 300
        scf_algorithm = diis
            scf_guess = sad
           sym_ignore = true
             symmetry = false
               thresh = 14
$end


$opt
CONSTRAINT
tors 6 8 9 10 0.0
ENDCONSTRAINT
$end
"""
        sections_test = QCInput.find_sections(str_single_job_input)
        section_actual = ["molecule", "rem", "opt"]
        self.assertEqual(section_actual, sections_test)

    def test_read_molecule(self):
        str_molecule = """$molecule
 0 1
 C     -9.5782000000      0.6241500000      0.0000000000
 O     -7.5827400000      0.5127000000     -0.0000000000
$end"""
        molecule_test = QCInput.read_molecule(str_molecule)
        species = ["C", "O"]
        coords = [[-9.5782000000, 0.6241500000, 0.0000000000],
                  [-7.5827400000, 0.5127000000, -0.0000000000]]
        molecule_actual = Molecule(species, coords)
        self.assertEqual(molecule_actual, molecule_test)

    def test_read_rem(self):
        str_rem = """Trying to break you!

$rem
   job_type  opt
  method  wB97M-V
   basis  def2-QZVPPD
   max_scf_cycles  300
  gen_scfman = true
$end"""
        rem_test = QCInput.read_rem(str_rem)
        rem_actual = {
            "job_type": "opt",
            "method": "wB97M-V",
            "basis": "def2-QZVPPD",
            "max_scf_cycles": "300",
            "gen_scfman": "true"
        }
        self.assertDictEqual(rem_actual, rem_test)


    def test_read_only_rem(self):
        str_rem = """Trying to break you!

$rem
   job_type  opt
  method  wB97M-V
   basis  def2-QZVPPD
   max_scf_cycles  300
  gen_scfman = true
$end

$pcm
heavypoints   194
hpoints   194
radii   uff
theory   cpcm
vdwscale   1.1
$end


$solvent
dielectric   10.0
$end


"""
        rem_test = QCInput.read_rem(str_rem)
        rem_actual = {
            "job_type": "opt",
            "method": "wB97M-V",
            "basis": "def2-QZVPPD",
            "max_scf_cycles": "300",
            "gen_scfman": "true"
        }
        self.assertDictEqual(rem_actual, rem_test)

    def test_read_opt(self):
        str_opt = """$opt
CONSTRAINT
  tors 2 3 4 5 25.0
   bend 2 1 4 110.0
ENDCONSTRAINT

FIXED
x y 2 4 5
ENDFIXED

DUMMY
   M 2 3 4 5
ENDDUMMY

CONNECT
4 3 2 3 5 6
ENDCONNECT
$end"""
        opt_test = QCInput.read_opt(str_opt)
        opt_actual = {
            "CONSTRAINT": ["tors 2 3 4 5 25.0", "bend 2 1 4 110.0"],
            "FIXED": ["x y 2 4 5"],
            "DUMMY": ["M 2 3 4 5"],
            "CONNECT": ["4 3 2 3 5 6"]
        }
        self.assertDictEqual(opt_actual, opt_test)

    def test__str__(self):
        species = ["C", "O"]
        coords = [[-9.5782000000, 0.6241500000, 0.0000000000],
                  [-7.5827400000, 0.5127000000, -0.0000000000]]
        molecule = Molecule(species=species, coords=coords)
        rem = {
            "jobtype": "opt",
            "method": "wB97M-V",
            "basis": "def2-QZVPPD",
            "max_scf_cycles": "300",
            "gen_scfman": "true"
        }
        str_test = QCInput(molecule=molecule, rem=rem).__str__()
        str_actual = """$molecule
 0 1
 C     -9.5782000000      0.6241500000      0.0000000000
 O     -7.5827400000      0.5127000000     -0.0000000000
$end

$rem
   job_type = opt
   method = wB97M-V
   basis = def2-QZVPPD
   max_scf_cycles = 300
   gen_scfman = true
$end
"""
        self.assertEqual(str_actual, str_test)

    def test_from_string(self):
        string = """$molecule
 0  1
 S          -0.00250959       -0.05817469       -0.02921636
 C           1.70755408       -0.03033788       -0.01382912
 H           2.24317221       -0.05215019        0.92026728
 C           2.21976393        0.01718014       -1.27293235
 H           3.27786220        0.04082146       -1.48539646
 C           1.20867399        0.04478540       -2.27007793
 H           1.40292257        0.10591684       -3.33110912
 C          -0.05341046        0.01577217       -1.74839343
 C          -1.32843436        0.03545064       -2.45531187
 C          -1.55195156        0.08743920       -3.80184635
 H          -0.75245172        0.10267657       -4.52817967
 C          -2.93293778        0.08408786       -4.13352169
 H          -3.31125108        0.11340328       -5.14405819
 C          -3.73173288        0.02741365       -3.03412864
 H          -4.80776535        0.00535688       -2.99564645
 S          -2.81590978       -0.00516172       -1.58990580
$end


$rem
              jobtype = opt
               method = wb97m-v
                basis = def2-tzvppd
           gen_scfman = true
  geom_opt_max_cycles = 75
       max_scf_cycles = 300
        scf_algorithm = diis
            scf_guess = sad
           sym_ignore = true
             symmetry = false
               thresh = 14
$end


$opt
CONSTRAINT
tors 6 8 9 10 0.0
ENDCONSTRAINT
$end
"""
        qcinput_test = QCInput.from_string(string)
        species = ["S", "C", "H", "C", "H", "C", "H", "C", "C", "C", "H", "C", "H", "C", "H", "S"]
        coords = [[-0.00250959, -0.05817469, -0.02921636],
                  [1.70755408, -0.03033788, -0.01382912],
                  [2.24317221, -0.05215019, 0.92026728],
                  [2.21976393, 0.01718014, -1.27293235],
                  [3.27786220, 0.04082146, -1.48539646],
                  [1.20867399, 0.04478540, -2.27007793],
                  [1.40292257, 0.10591684, -3.33110912],
                  [-0.05341046, 0.01577217, -1.74839343],
                  [-1.32843436, 0.03545064, -2.45531187],
                  [-1.55195156, 0.08743920, -3.80184635],
                  [-0.75245172, 0.10267657, -4.52817967],
                  [-2.93293778, 0.08408786, -4.13352169],
                  [-3.31125108, 0.11340328, -5.14405819],
                  [-3.73173288, 0.02741365, -3.03412864],
                  [-4.80776535, 0.00535688, -2.99564645],
                  [-2.81590978, -0.00516172, -1.58990580]]
        molecule_actual = Molecule(species, coords)
        self.assertEqual(molecule_actual, qcinput_test.molecule)
        rem_actual = {
            "job_type": "opt",
            "method": "wb97m-v",
            "basis": "def2-tzvppd",
            "gen_scfman": "true",
            "geom_opt_max_cycles": "75",
            "max_scf_cycles": "300",
            "scf_algorithm": "diis",
            "scf_guess": "sad",
            "sym_ignore": "true",
            "symmetry": "false",
            "thresh": "14"
        }
        self.assertDictEqual(rem_actual, qcinput_test.rem)
        opt_actual = {"CONSTRAINT": ["tors 6 8 9 10 0.0"]}
        self.assertDictEqual(opt_actual, qcinput_test.opt)

    def test_multi_job_string(self):
        species = ["S", "C", "H", "C", "H", "C", "H", "C", "C", "C", "H", "C", "H", "C", "H", "S"]
        coords = [[-0.00250959, -0.05817469, -0.02921636],
                  [1.70755408, -0.03033788, -0.01382912],
                  [2.24317221, -0.05215019, 0.92026728],
                  [2.21976393, 0.01718014, -1.27293235],
                  [3.27786220, 0.04082146, -1.48539646],
                  [1.20867399, 0.04478540, -2.27007793],
                  [1.40292257, 0.10591684, -3.33110912],
                  [-0.05341046, 0.01577217, -1.74839343],
                  [-1.32843436, 0.03545064, -2.45531187],
                  [-1.55195156, 0.08743920, -3.80184635],
                  [-0.75245172, 0.10267657, -4.52817967],
                  [-2.93293778, 0.08408786, -4.13352169],
                  [-3.31125108, 0.11340328, -5.14405819],
                  [-3.73173288, 0.02741365, -3.03412864],
                  [-4.80776535, 0.00535688, -2.99564645],
                  [-2.81590978, -0.00516172, -1.58990580]]
        molecule_1 = Molecule(species, coords)
        rem_1 = {
            "jobtype": "opt",
            "method": "wb97m-v",
            "basis": "def2-tzvppd",
            "gen_scfman": "true",
            "geom_opt_max_cycles": "75",
            "max_scf_cycles": "300",
            "scf_algorithm": "diis",
            "scf_guess": "sad",
            "sym_ignore": "true",
            "symmetry": "false",
            "thresh": "14"
        }
        opt_1 = {"CONSTRAINT": ["tors 6 8 9 10 0.0"]}
        job_1 = QCInput(molecule=molecule_1, rem=rem_1, opt=opt_1)
        molecule_2 = "read"
        rem_2 = {
            "jobtype": "sp",
            "method": "wb97m-v",
            "basis": "def2-tzvppd",
            "gen_scfman": "true",
            "geom_opt_max_cycles": "75",
            "max_scf_cycles": "300",
            "scf_algorithm": "diis",
            "scf_guess": "read",
            "sym_ignore": "true",
            "symmetry": "false",
            "thresh": "14"
        }
        job_2 = QCInput(molecule=molecule_2, rem=rem_2)
        job_list = [job_1, job_2]
        multi_job_str_test = QCInput.multi_job_string(job_list=job_list)
        multi_job_str_actual = """$molecule
 0 1
 S     -0.0025095900     -0.0581746900     -0.0292163600
 C      1.7075540800     -0.0303378800     -0.0138291200
 H      2.2431722100     -0.0521501900      0.9202672800
 C      2.2197639300      0.0171801400     -1.2729323500
 H      3.2778622000      0.0408214600     -1.4853964600
 C      1.2086739900      0.0447854000     -2.2700779300
 H      1.4029225700      0.1059168400     -3.3311091200
 C     -0.0534104600      0.0157721700     -1.7483934300
 C     -1.3284343600      0.0354506400     -2.4553118700
 C     -1.5519515600      0.0874392000     -3.8018463500
 H     -0.7524517200      0.1026765700     -4.5281796700
 C     -2.9329377800      0.0840878600     -4.1335216900
 H     -3.3112510800      0.1134032800     -5.1440581900
 C     -3.7317328800      0.0274136500     -3.0341286400
 H     -4.8077653500      0.0053568800     -2.9956464500
 S     -2.8159097800     -0.0051617200     -1.5899058000
$end

$rem
   job_type = opt
   method = wb97m-v
   basis = def2-tzvppd
   gen_scfman = true
   geom_opt_max_cycles = 75
   max_scf_cycles = 300
   scf_algorithm = diis
   scf_guess = sad
   sym_ignore = true
   symmetry = false
   thresh = 14
$end

$opt
CONSTRAINT
   tors 6 8 9 10 0.0
ENDCONSTRAINT
$end

@@@

$molecule
 read
$end

$rem
   job_type = sp
   method = wb97m-v
   basis = def2-tzvppd
   gen_scfman = true
   geom_opt_max_cycles = 75
   max_scf_cycles = 300
   scf_algorithm = diis
   scf_guess = read
   sym_ignore = true
   symmetry = false
   thresh = 14
$end
"""
        self.assertEqual(multi_job_str_actual, multi_job_str_test)

    def test_from_multi_jobs_file(self):
        job_list_test = QCInput.from_multi_jobs_file(os.path.join(test_dir, "pt_n2_wb97mv_0.0.in"))
        species = ["S", "C", "H", "C", "H", "C", "H", "C", "C", "C", "H", "C", "H", "C", "H", "S"]
        coords = [[-0.00250959, -0.05817469, -0.02921636],
                  [1.70755408, -0.03033788, -0.01382912],
                  [2.24317221, -0.05215019, 0.92026728],
                  [2.21976393, 0.01718014, -1.27293235],
                  [3.27786220, 0.04082146, -1.48539646],
                  [1.20867399, 0.04478540, -2.27007793],
                  [1.40292257, 0.10591684, -3.33110912],
                  [-0.05341046, 0.01577217, -1.74839343],
                  [-1.32843436, 0.03545064, -2.45531187],
                  [-1.55195156, 0.08743920, -3.80184635],
                  [-0.75245172, 0.10267657, -4.52817967],
                  [-2.93293778, 0.08408786, -4.13352169],
                  [-3.31125108, 0.11340328, -5.14405819],
                  [-3.73173288, 0.02741365, -3.03412864],
                  [-4.80776535, 0.00535688, -2.99564645],
                  [-2.81590978, -0.00516172, -1.58990580]]
        molecule_1_actual = Molecule(species, coords)
        rem_1_actual = {
            "job_type": "opt",
            "method": "wb97m-v",
            "basis": "def2-tzvppd",
            "gen_scfman": "true",
            "geom_opt_max_cycles": "75",
            "max_scf_cycles": "300",
            "scf_algorithm": "diis",
            "scf_guess": "sad",
            "sym_ignore": "true",
            "symmetry": "false",
            "thresh": "14"
        }
        opt_1_actual = {"CONSTRAINT": ["tors 6 8 9 10 0.0"]}
        self.assertEqual(molecule_1_actual, job_list_test[0].molecule)
        self.assertEqual(rem_1_actual, job_list_test[0].rem)
        self.assertEqual(opt_1_actual, job_list_test[0].opt)

        molecule_2_actual = "read"
        rem_2_actual = {
            "job_type": "sp",
            "method": "wb97m-v",
            "basis": "def2-tzvppd",
            "gen_scfman": "true",
            "geom_opt_max_cycles": "75",
            "max_scf_cycles": "300",
            "scf_algorithm": "diis",
            "scf_guess": "read",
            "sym_ignore": "true",
            "symmetry": "false",
            "thresh": "14"
        }
        self.assertEqual(molecule_2_actual, job_list_test[1].molecule)
        self.assertEqual(rem_2_actual, job_list_test[1].rem)

    def test_read_pcm(self):
        str_pcm = """I'm once again trying to break you!

$pcm
   theory cpcm
   radii uff
   vdwscale 1.1
$end"""
        pcm_test = QCInput.read_pcm(str_pcm)
        pcm_actual = {
             "theory": "cpcm",
             "radii": "uff",
             "vdwscale": "1.1"
        }
        self.assertDictEqual(pcm_actual, pcm_test)

    def test_read_bad_pcm(self):
        str_pcm = """I'm once again trying to break you!

$pcm
   theory = cpcm
   radii = uff
   vdwscale = 1.1
$end"""
        pcm_test = QCInput.read_pcm(str_pcm)
        pcm_actual = {}
        self.assertDictEqual(pcm_actual, pcm_test)

    def test_read_solvent(self):
        str_solvent = """Once again, I'm trying to break you!

$solvent
   dielectric 5.0
$end"""
        solvent_test = QCInput.read_solvent(str_solvent)
        solvent_actual = {
             "dielectric": "5.0",
        }
        self.assertDictEqual(solvent_actual, solvent_test)

    def test_read_bad_solvent(self):
        str_solvent = """Once again, I'm trying to break you!

$solvent
   dielectric = 5.0
$end"""
        solvent_test = QCInput.read_solvent(str_solvent)
        solvent_actual = {}
        self.assertDictEqual(solvent_actual, solvent_test)

if __name__ == "__main__":
    unittest.main()
