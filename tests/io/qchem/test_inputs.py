from __future__ import annotations

import logging
import os

import pytest
from monty.serialization import loadfn

from pymatgen.core.structure import Molecule
from pymatgen.io.qchem.inputs import QCInput
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest

__author__ = "Brandon Wood, Samuel Blau, Shyam Dwaraknath, Julian Self, Evan Spotte-Smith, Ryan Kingsbury"
__copyright__ = "Copyright 2018-2022, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Samuel Blau"
__email__ = "samblau1@gmail.com"
__credits__ = "Xiaohui Qu"

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
logger = logging.getLogger(__name__)


class TestQCInput(PymatgenTest):
    # ef setUpClass(cls):
    # add things that show up over and over again

    def test_molecule_template(self):
        species = ["C", "O"]
        coords = [
            [-9.5782000000, 0.6241500000, 0.0000000000],
            [-7.5827400000, 0.5127000000, -0.0000000000],
        ]
        mol = Molecule(species=species, coords=coords)
        molecule_test = QCInput.molecule_template(mol)
        molecule_actual = """$molecule
 0 1
 C     -9.5782000000      0.6241500000      0.0000000000
 O     -7.5827400000      0.5127000000     -0.0000000000
$end"""

        assert molecule_actual == molecule_test

    def test_multi_molecule_template(self):
        self.maxDiff = None
        species = ["C", "C", "H", "H", "H", "H"]
        coords_1 = [
            [0.000000, 0.000000, 0.000000],
            [1.332000, 0.000000, 0.000000],
            [-0.574301, 0.000000, -0.928785],
            [-0.574301, 0.000000, 0.928785],
            [1.906301, 0.000000, 0.928785],
            [1.906301, 0.000000, -0.928785],
        ]
        coords_2 = [
            [0.000000, 4.000000, 0.000000],
            [1.332000, 4.000000, 0.000000],
            [-0.574301, 4.000000, -0.928785],
            [-0.574301, 4.000000, 0.928785],
            [1.906301, 4.000000, 0.928785],
            [1.906301, 4.000000, -0.928785],
        ]

        mol_1 = Molecule(species, coords_1, charge=1)
        mol_2 = Molecule(species, coords_2)
        molecule_test = QCInput.molecule_template([mol_1, mol_2])
        molecule_actual = """$molecule
 1 2
--
 1 2
 C      0.0000000000      0.0000000000      0.0000000000
 C      1.3320000000      0.0000000000      0.0000000000
 H     -0.5743010000      0.0000000000     -0.9287850000
 H     -0.5743010000      0.0000000000      0.9287850000
 H      1.9063010000      0.0000000000      0.9287850000
 H      1.9063010000      0.0000000000     -0.9287850000
--
 0 1
 C      0.0000000000      4.0000000000      0.0000000000
 C      1.3320000000      4.0000000000      0.0000000000
 H     -0.5743010000      4.0000000000     -0.9287850000
 H     -0.5743010000      4.0000000000      0.9287850000
 H      1.9063010000      4.0000000000      0.9287850000
 H      1.9063010000      4.0000000000     -0.9287850000
$end"""

        assert molecule_test == molecule_actual

    # TODO improve this test maybe add ordered dicts
    def test_rem_template(self):
        rem_params = {
            "job_type": "opt",
            "method": "wb97m-v",
            "basis": "def2-qzvppd",
            "max_scf_cycles": 300,
            "gen_scfman": "true",
        }
        rem_test = QCInput.rem_template(rem_params).split("\n")
        rem_actual_list = [
            "$rem",
            "   job_type = opt",
            "   method = wb97m-v",
            "   basis = def2-qzvppd",
            "   max_scf_cycles = 300",
            "   gen_scfman = true",
            "$end",
        ]

        for i_rem in rem_actual_list:
            assert i_rem in rem_test

    def test_opt_template(self):
        opt_params = {
            "CONSTRAINT": ["tors 2 3 4 5 25.0", "bend 2 1 4 110.0"],
            "FIXED": ["x y 2 4 5"],
            "DUMMY": ["M 2 3 4 5"],
            "CONNECT": ["4 3 2 3 5 6"],
        }
        opt_test = QCInput.opt_template(opt_params).split("\n")
        opt_actual_list = [
            "$opt",
            "CONSTRAINT",
            "   tors 2 3 4 5 25.0",
            "   bend 2 1 4 110.0",
            "ENDCONSTRAINT",
            "FIXED",
            "   x y 2 4 5",
            "ENDFIXED",
            "DUMMY",
            "   M 2 3 4 5",
            "ENDDUMMY",
            "CONNECT",
            "   4 3 2 3 5 6",
            "ENDCONNECT",
            "$end",
        ]

        for i_opt in opt_actual_list:
            assert i_opt in opt_test

    def test_pcm_template(self):
        pcm_params = {"theory": "cpcm"}
        pcm_test = QCInput.pcm_template(pcm_params)
        pcm_actual = """$pcm
   theory cpcm
$end"""
        assert pcm_actual == pcm_test

    def test_pcm_nonels_template(self):
        # make sure values that are None get skipped in the output
        pcm_nonels = {
            "A": "-0.006736",
            "B": "0.032698",
            "C": "-1249.6",
            "D": None,
            "Delta": "7.0",
            "Gamma": "3.7",
            "SolvRho": "0.05",
            "GauLag_N": "40",
        }
        pcm_nonels_test = QCInput.pcm_nonels_template(pcm_nonels)
        pcm_nonels_actual = """$pcm_nonels
   A -0.006736
   B 0.032698
   C -1249.6
   Delta 7.0
   Gamma 3.7
   SolvRho 0.05
   GauLag_N 40
$end"""
        assert pcm_nonels_actual == pcm_nonels_test

    def test_solvent_template(self):
        solvent_params = {"dielectric": "5.0"}
        solvent_test = QCInput.solvent_template(solvent_params)
        solvent_actual = """$solvent
   dielectric 5.0
$end"""
        assert solvent_actual == solvent_test

    def test_smx_template(self):
        smx_params = {"solvent": "water"}
        smx_test = QCInput.smx_template(smx_params)
        smx_actual = """$smx
   solvent water
$end"""
        assert smx_actual == smx_test

        smx_params = {"solvent": "dimethyl sulfoxide"}
        smx_test = QCInput.smx_template(smx_params)
        smx_actual = """$smx
   solvent dmso
$end"""
        assert smx_actual == smx_test

    def test_svp_template(self):
        svp_params = {
            "RHOISO": 0.001,
            "DIELST": 78.36,
            "NPTLEB": 1202,
            "ITRNGR": 2,
            "IROTGR": 2,
            "IPNRF": 1,
            "IDEFESR": 1,
        }
        # svp_params = lower_and_check_unique(svp_params)
        svp_test = QCInput.svp_template(svp_params)
        svp_actual = """$svp
RHOISO=0.001, DIELST=78.36, NPTLEB=1202, ITRNGR=2, IROTGR=2, IPNRF=1, IDEFESR=1
$end"""
        assert svp_actual == svp_test

    def test_scan_template(self):
        scan_params = {"stre": ["3 6 1.5 1.9 0.01"], "tors": ["1 2 3 4 -180 180 30"]}
        scan_test = QCInput.scan_template(scan_params)
        scan_actual = """$scan
   stre 3 6 1.5 1.9 0.01
   tors 1 2 3 4 -180 180 30
$end"""
        assert scan_test == scan_actual

        bad_scan = {"stre": ["1 2 1.0 2.0 0.05", "3 4 1.5 2.0 0.05"], "bend": ["7 8 9 90 120 10"]}
        with pytest.raises(ValueError, match="Q-Chem only supports PES_SCAN with two or less variables"):
            QCInput.scan_template(bad_scan)

    def test_van_der_waals_template(self):
        vdw_params = {1: 1.20, 12: 1.72}
        vdw_test_atomic = QCInput.van_der_waals_template(vdw_params, mode="atomic")
        vdw_actual_atomic = """$van_der_waals
1
   1 1.2
   12 1.72
$end"""
        assert vdw_test_atomic == vdw_actual_atomic

        vdw_test_sequential = QCInput.van_der_waals_template(vdw_params, mode="sequential")
        vdw_actual_sequential = """$van_der_waals
2
   1 1.2
   12 1.72
$end"""
        assert vdw_test_sequential == vdw_actual_sequential
        mode = "mymode"
        with pytest.raises(ValueError, match=f"Invalid {mode=}, must be 'atomic' or 'sequential'"):
            QCInput.van_der_waals_template(vdw_params, mode=mode)

    def test_cdft_template(self):
        cdft = [
            [
                {"value": 1.0, "coefficients": [1.0], "first_atoms": [1], "last_atoms": [27], "types": ["c"]},
                {"value": 0.0, "coefficients": [1.0], "first_atoms": [1], "last_atoms": [27], "types": ["s"]},
            ],
            [
                {"value": 0.0, "coefficients": [1.0], "first_atoms": [1], "last_atoms": [27], "types": ["c"]},
                {"value": -1.0, "coefficients": [1.0], "first_atoms": [1], "last_atoms": [27], "types": ["s"]},
            ],
        ]

        cdft_test = QCInput.cdft_template(cdft)
        cdft_actual = """$cdft
   1.0
   1.0 1 27
   0.0
   1.0 1 27 s
--------------
   0.0
   1.0 1 27
   -1.0
   1.0 1 27 s
$end"""

        assert cdft_test == cdft_actual

    def test_almo_template(self):
        almo = [[(1, 2), (0, 1)], [(0, 1), (1, 2)]]
        almo_test = QCInput.almo_template(almo)
        almo_actual = """$almo_coupling
   1 2
   0 1
   --
   0 1
   1 2
$end"""

        assert almo_test == almo_actual

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
        assert section_actual == sections_test

    def test_read_molecule(self):
        str_molecule = """$molecule
 0 1
 C     -9.5782000000      0.6241500000      0.0000000000
 O     -7.5827400000      0.5127000000     -0.0000000000
$end"""
        molecule_test = QCInput.read_molecule(str_molecule)
        species = ["C", "O"]
        coords = [
            [-9.5782000000, 0.6241500000, 0.0000000000],
            [-7.5827400000, 0.5127000000, -0.0000000000],
        ]
        molecule_actual = Molecule(species, coords)
        assert molecule_actual == molecule_test

    def test_read_multi_molecule(self):
        str_molecule = """$molecule
 1 2
--
 1 2
 C      0.0000000000      0.0000000000      0.0000000000
 C      1.3320000000      0.0000000000      0.0000000000
 H     -0.5743010000      0.0000000000     -0.9287850000
 H     -0.5743010000      0.0000000000      0.9287850000
 H      1.9063010000      0.0000000000      0.9287850000
 H      1.9063010000      0.0000000000     -0.9287850000
--
 0 1
 C      0.0000000000      4.0000000000      0.0000000000
 C      1.3320000000      4.0000000000      0.0000000000
 H     -0.5743010000      4.0000000000     -0.9287850000
 H     -0.5743010000      4.0000000000      0.9287850000
 H      1.9063010000      4.0000000000      0.9287850000
 H      1.9063010000      4.0000000000     -0.9287850000
$end"""

        species = ["C", "C", "H", "H", "H", "H"]
        coords_1 = [
            [0.000000, 0.000000, 0.000000],
            [1.332000, 0.000000, 0.000000],
            [-0.574301, 0.000000, -0.928785],
            [-0.574301, 0.000000, 0.928785],
            [1.906301, 0.000000, 0.928785],
            [1.906301, 0.000000, -0.928785],
        ]
        coords_2 = [
            [0.000000, 4.000000, 0.000000],
            [1.332000, 4.000000, 0.000000],
            [-0.574301, 4.000000, -0.928785],
            [-0.574301, 4.000000, 0.928785],
            [1.906301, 4.000000, 0.928785],
            [1.906301, 4.000000, -0.928785],
        ]

        mol_1 = Molecule(species, coords_1, charge=1)
        mol_2 = Molecule(species, coords_2)

        parsed = QCInput.read_molecule(str_molecule)
        assert parsed[0], mol_1
        assert parsed[1], mol_2

    def test_read_rem(self):
        str_rem = """Trying to break you!

$rem
   job_type  opt
  method  wb97m-v
   basis  def2-qzvppd
   max_scf_cycles  300
  gen_scfman = true
$end"""
        rem_test = QCInput.read_rem(str_rem)
        rem_actual = {
            "job_type": "opt",
            "method": "wb97m-v",
            "basis": "def2-qzvppd",
            "max_scf_cycles": "300",
            "gen_scfman": "true",
        }
        assert rem_actual == rem_test

    def test_read_only_rem(self):
        str_rem = """Trying to break you!

$rem
   job_type  opt
  method  wb97m-v
   basis  def2-qzvppd
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
            "method": "wb97m-v",
            "basis": "def2-qzvppd",
            "max_scf_cycles": "300",
            "gen_scfman": "true",
        }
        assert rem_actual == rem_test

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
            "CONNECT": ["4 3 2 3 5 6"],
        }
        assert opt_actual == opt_test

    def test__str__(self):
        species = ["C", "O"]
        coords = [
            [-9.5782000000, 0.6241500000, 0.0000000000],
            [-7.5827400000, 0.5127000000, -0.0000000000],
        ]
        molecule = Molecule(species=species, coords=coords)
        rem = {
            "jobtype": "opt",
            "method": "wb97m-v",
            "basis": "def2-qzvppd",
            "max_scf_cycles": "300",
            "gen_scfman": "true",
        }
        str_test = str(QCInput(molecule=molecule, rem=rem)).split("\n")
        str_actual_list = [
            "$molecule",
            " 0 1",
            " C     -9.5782000000      0.6241500000      0.0000000000",
            " O     -7.5827400000      0.5127000000     -0.0000000000",
            "$end",
            "$rem",
            "   job_type = opt",
            "   method = wb97m-v",
            "   basis = def2-qzvppd",
            "   max_scf_cycles = 300",
            "   gen_scfman = true",
            "$end",
        ]

        for i_str in str_actual_list:
            assert i_str in str_test

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
        qcinput_test = QCInput.from_str(string)
        species = [
            "S",
            "C",
            "H",
            "C",
            "H",
            "C",
            "H",
            "C",
            "C",
            "C",
            "H",
            "C",
            "H",
            "C",
            "H",
            "S",
        ]
        coords = [
            [-0.00250959, -0.05817469, -0.02921636],
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
            [-2.81590978, -0.00516172, -1.58990580],
        ]
        molecule_actual = Molecule(species, coords)
        assert molecule_actual == qcinput_test.molecule
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
            "thresh": "14",
        }
        assert rem_actual == qcinput_test.rem
        opt_actual = {"CONSTRAINT": ["tors 6 8 9 10 0.0"]}
        assert opt_actual == qcinput_test.opt

    # TODO this test needs an update, the assertion doesn't differentiate between the different rem sections
    def test_multi_job_string(self):
        species = [
            "S",
            "C",
            "H",
            "C",
            "H",
            "C",
            "H",
            "C",
            "C",
            "C",
            "H",
            "C",
            "H",
            "C",
            "H",
            "S",
        ]
        coords = [
            [-0.00250959, -0.05817469, -0.02921636],
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
            [-2.81590978, -0.00516172, -1.58990580],
        ]
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
            "thresh": "14",
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
            "thresh": "14",
        }
        job_2 = QCInput(molecule=molecule_2, rem=rem_2)
        job_list = [job_1, job_2]
        multi_job_str_test = QCInput.multi_job_string(job_list=job_list).split("\n")
        multi_job_str_actual_list = [
            "$molecule",
            " 0 1",
            " S     -0.0025095900     -0.0581746900     -0.0292163600",
            " C      1.7075540800     -0.0303378800     -0.0138291200",
            " H      2.2431722100     -0.0521501900      0.9202672800",
            " C      2.2197639300      0.0171801400     -1.2729323500",
            " H      3.2778622000      0.0408214600     -1.4853964600",
            " C      1.2086739900      0.0447854000     -2.2700779300",
            " H      1.4029225700      0.1059168400     -3.3311091200",
            " C     -0.0534104600      0.0157721700     -1.7483934300",
            " C     -1.3284343600      0.0354506400     -2.4553118700",
            " C     -1.5519515600      0.0874392000     -3.8018463500",
            " H     -0.7524517200      0.1026765700     -4.5281796700",
            " C     -2.9329377800      0.0840878600     -4.1335216900",
            " H     -3.3112510800      0.1134032800     -5.1440581900",
            " C     -3.7317328800      0.0274136500     -3.0341286400",
            " H     -4.8077653500      0.0053568800     -2.9956464500",
            " S     -2.8159097800     -0.0051617200     -1.5899058000",
            "$end",
            "$rem",
            "   job_type = opt",
            "   method = wb97m-v",
            "   basis = def2-tzvppd",
            "   gen_scfman = true",
            "   geom_opt_max_cycles = 75",
            "   max_scf_cycles = 300",
            "   scf_algorithm = diis",
            "   scf_guess = sad",
            "   sym_ignore = true",
            "   symmetry = false",
            "   thresh = 14",
            "$end",
            "$opt",
            "CONSTRAINT",
            "   tors 6 8 9 10 0.0",
            "ENDCONSTRAINT",
            "$end",
            "@@@",
            "$molecule",
            " read",
            "$end",
            "$rem",
            "   job_type = opt",
            "   method = wb97m-v",
            "   basis = def2-tzvppd",
            "   gen_scfman = true",
            "   geom_opt_max_cycles = 75",
            "   max_scf_cycles = 300",
            "   scf_algorithm = diis",
            "   scf_guess = sad",
            "   sym_ignore = true",
            "   symmetry = false",
            "   thresh = 14",
            "$end",
        ]

        for i_str in multi_job_str_actual_list:
            assert i_str in multi_job_str_test

    def test_from_multi_jobs_file(self):
        job_list_test = QCInput.from_multi_jobs_file(f"{TEST_FILES_DIR}/qchem/pt_n2_wb97mv_0.0.in")
        species = [
            "S",
            "C",
            "H",
            "C",
            "H",
            "C",
            "H",
            "C",
            "C",
            "C",
            "H",
            "C",
            "H",
            "C",
            "H",
            "S",
        ]
        coords = [
            [-0.00250959, -0.05817469, -0.02921636],
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
            [-2.81590978, -0.00516172, -1.58990580],
        ]
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
            "thresh": "14",
        }
        opt_1_actual = {"CONSTRAINT": ["tors 6 8 9 10 0.0"]}
        assert molecule_1_actual == job_list_test[0].molecule
        assert rem_1_actual == job_list_test[0].rem
        assert opt_1_actual == job_list_test[0].opt

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
            "thresh": "14",
        }
        assert molecule_2_actual == job_list_test[1].molecule
        assert rem_2_actual == job_list_test[1].rem

    def test_read_pcm(self):
        str_pcm = """I'm once again trying to break you!

$pcm
   theory cpcm
   radii uff
   vdwscale 1.1
$end"""
        pcm_test = QCInput.read_pcm(str_pcm)
        pcm_actual = {"theory": "cpcm", "radii": "uff", "vdwscale": "1.1"}
        assert pcm_actual == pcm_test

    def test_read_pcm_nonels(self):
        str_pcm_nonels = """$pcm_nonels
   A         -0.006736
   B          0.032698
   C      -1249.6
   D        -21.405
   Delta    7.0
   Gamma    3.7
   SolvRho  0.05
   GauLag_N 40
$end"""
        pcm_nonels_test = QCInput.read_pcm_nonels(str_pcm_nonels)
        pcm_nonels_actual = {
            "A": "-0.006736",
            "B": "0.032698",
            "C": "-1249.6",
            "D": "-21.405",
            "Delta": "7.0",
            "Gamma": "3.7",
            "SolvRho": "0.05",
            "GauLag_N": "40",
        }
        assert pcm_nonels_actual == pcm_nonels_test

    def test_read_bad_pcm(self):
        str_pcm = """I'm once again trying to break you!

$pcm
   theory = cpcm
   radii = uff
   vdwscale = 1.1
$end"""
        pcm_test = QCInput.read_pcm(str_pcm)
        pcm_actual = {}
        assert pcm_actual == pcm_test

    def test_read_solvent(self):
        str_solvent = """Once again, I'm trying to break you!

$solvent
   dielectric 5.0
$end"""
        solvent_test = QCInput.read_solvent(str_solvent)
        solvent_actual = {
            "dielectric": "5.0",
        }
        assert solvent_actual == solvent_test

    def test_read_bad_solvent(self):
        str_solvent = """Once again, I'm trying to break you!

$solvent
   dielectric = 5.0
$end"""
        solvent_test = QCInput.read_solvent(str_solvent)
        solvent_actual = {}
        assert solvent_actual == solvent_test

    def test_read_smx(self):
        str_smx = """Once again, I'm trying to break you!

$smx
   solvent water
$end"""
        smx_test = QCInput.read_smx(str_smx)
        smx_actual = {
            "solvent": "water",
        }
        assert smx_actual == smx_test

    def test_read_svp(self):
        str_svp = """$svp
RHOISO=0.001, DIELST=78.36, NPTLEB=1202, ITRNGR=2, IROTGR=2, IPNRF=1, IDEFESR=1
$end"""
        svp_test = QCInput.read_svp(str_svp)
        svp_actual = {
            "RHOISO": "0.001",
            "DIELST": "78.36",
            "NPTLEB": "1202",
            "ITRNGR": "2",
            "IROTGR": "2",
            "IPNRF": "1",
            "IDEFESR": "1",
        }
        assert svp_actual == svp_test

    def test_read_bad_smx(self):
        str_smx = """Once again, I'm trying to break you!

$solvent
   solvent = water
$end"""
        smx_test = QCInput.read_smx(str_smx)
        smx_actual = {}
        assert smx_actual == smx_test

    def test_read_scan(self):
        str_scan = """Once more, I'm trying to break you!

$scan
   stre 1 2 1.1 1.4 0.03
   bend 3 4 5 60 90 5
$end"""
        scan_test = QCInput.read_scan(str_scan)
        scan_actual = {"stre": ["1 2 1.1 1.4 0.03"], "bend": ["3 4 5 60 90 5"], "tors": []}

        assert scan_test == scan_actual

    def test_read_bad_scan(self):
        str_scan_1 = """Once more, I"m trying to break you!
$scan
   boo 1 4 1.2 1.5 0.02
   tors = 3 6 1.5 1.9 0.01
$end
"""
        scan_test_1 = QCInput.read_scan(str_scan_1)
        scan_actual_1 = {}
        assert scan_test_1 == scan_actual_1

        str_scan_2 = """Once more, I'm trying to break you!

$scan
   stre 1 2 1.1 1.4 0.03
   bend 3 4 5 60 90 5
   tors 6 7 8 9 -180 180 30
$end"""

        with pytest.raises(ValueError, match="No more than two variables are allows in the scan section"):
            QCInput.read_scan(str_scan_2)

    def test_read_negative(self):
        str_molecule = """$molecule
 -1 1
 S     -1.1516880000      0.8568110000     -0.0787470000
 S      1.1527500000     -0.8580450000     -0.0786430000
 O     -1.6523520000      1.8607750000     -1.0252100000
 O     -0.9052880000      1.2448490000      1.3156410000
 O      0.9072410000     -1.2461780000      1.3158760000
 O      1.6543670000     -1.8616640000     -1.0249090000
 C     -2.5841130000     -0.3746500000      0.0297340000
 C      2.5833220000      0.3755850000      0.0296900000
 F     -3.6480730000      0.2204040000      0.6112110000
 F     -2.2609850000     -1.4531020000      0.7616580000
 F     -2.9656640000     -0.7966010000     -1.1900330000
 F      3.6467050000     -0.2152590000      0.6163310000
 F      2.2560700000      1.4560310000      0.7568190000
 F      2.9672080000      0.7933560000     -1.1908790000
 N     -0.0001900000     -0.0016540000     -0.8250640000
$end

$rem
   job_type = opt
   basis = 6-311++g*
   max_scf_cycles = 200
   gen_scfman = true
   scf_algorithm = diis
   method = wb97xd
   geom_opt_max_cycles = 200
$end
"""
        qcinp = QCInput.from_str(str_molecule)
        assert str_molecule == str(qcinp)

    def test_read_plots(self):
        str_molecule = """$molecule
 0 2
 O      1.6159947668      0.3522275191      0.3343192028
 O     -0.5921658045      1.4368355787      1.2632324885
 C      0.4160355545     -0.4617433561      0.2180766834
 C     -0.7655230468      0.4776728409      0.1826587618
 C      2.8437090411     -0.3853724291      0.0935770045
 C     -1.7918488579      2.2003569978      1.5593659974
 H      0.4649228147     -1.0347597878     -0.7097270414
 H      3.6714833661      0.3051154983      0.2509025369
 H      2.8395611019     -0.7401009356     -0.9372741555
 H     -2.1017802975      2.7482577804      0.6678359687
 H     -1.5445030956      2.8894960726      2.3658396091
 Mg      1.2856817013      1.9249743897      1.4285694502
$end

$rem
   job_type = sp
   basis = def2-tzvppd
   max_scf_cycles = 200
   gen_scfman = true
   xc_grid = 3
   scf_algorithm = gdm
   resp_charges = true
   symmetry = false
   sym_ignore = true
   method = wb97xv
   solvent_method = smd
   ideriv = 1
   thresh = 14
   scf_guess_always = true
   plots = true
   make_cube_files = true
$end

$smx
   solvent thf
$end

$plots
   grid_spacing 0.05
   total_density 0
$end
"""
        qcinp = QCInput.from_str(str_molecule)
        assert str_molecule == str(qcinp)

    def test_read_nbo(self):
        str_molecule = """$molecule
 0 2
 C     -2.0338520000      0.0865500000     -1.4158570000
 C     -1.2819580000      0.3850830000     -0.1564990000
 C     -2.0067300000      1.1271820000      0.9225950000
 C      0.1219120000     -0.0366190000      0.0148810000
 C      0.6767790000     -1.0507090000     -0.7802400000
 C      2.0072450000     -1.4517610000     -0.6185380000
 C      2.8079970000     -0.8434840000      0.3427930000
 C      2.2778880000      0.1645690000      1.1416530000
 C      0.9468200000      0.5630060000      0.9784410000
 H     -1.3919850000      0.1591240000     -2.2995570000
 H     -2.4671570000     -0.9174600000     -1.3722490000
 H     -2.8505080000      0.8017250000     -1.5613060000
 H     -3.0889210000      0.9823990000      0.8362370000
 H     -1.7216740000      0.7761670000      1.9194500000
 H     -1.8021560000      2.1999010000      0.8510710000
 H      0.0793240000     -1.5592640000     -1.5324310000
 H      2.4136820000     -2.2421190000     -1.2440900000
 H      3.8415290000     -1.1539430000      0.4689660000
 H      2.8984450000      0.6464300000      1.8925800000
 H      0.5733200000      1.3632210000      1.6120990000
$end

$rem
   job_type = sp
   max_scf_cycles = 200
   gen_scfman = true
   xc_grid = 3
   scf_algorithm = diis
   method = wb97xv
   basis = def2-tzvp
   symmetry = false
   sym_ignore = true
   nbo = true
$end

$nbo
$end
"""
        qcinp = QCInput.from_str(str_molecule)
        assert str_molecule == str(qcinp)

        str_molecule = """$molecule
 0 2
 C     -2.0338520000      0.0865500000     -1.4158570000
 C     -1.2819580000      0.3850830000     -0.1564990000
 C     -2.0067300000      1.1271820000      0.9225950000
 C      0.1219120000     -0.0366190000      0.0148810000
 C      0.6767790000     -1.0507090000     -0.7802400000
 C      2.0072450000     -1.4517610000     -0.6185380000
 C      2.8079970000     -0.8434840000      0.3427930000
 C      2.2778880000      0.1645690000      1.1416530000
 C      0.9468200000      0.5630060000      0.9784410000
 H     -1.3919850000      0.1591240000     -2.2995570000
 H     -2.4671570000     -0.9174600000     -1.3722490000
 H     -2.8505080000      0.8017250000     -1.5613060000
 H     -3.0889210000      0.9823990000      0.8362370000
 H     -1.7216740000      0.7761670000      1.9194500000
 H     -1.8021560000      2.1999010000      0.8510710000
 H      0.0793240000     -1.5592640000     -1.5324310000
 H      2.4136820000     -2.2421190000     -1.2440900000
 H      3.8415290000     -1.1539430000      0.4689660000
 H      2.8984450000      0.6464300000      1.8925800000
 H      0.5733200000      1.3632210000      1.6120990000
$end

$rem
   job_type = sp
   max_scf_cycles = 200
   gen_scfman = true
   xc_grid = 3
   scf_algorithm = diis
   method = wb97xv
   basis = def2-tzvp
   symmetry = false
   sym_ignore = true
   nbo = true
$end

$nbo
   print = 1
$end
"""
        qcinp = QCInput.from_str(str_molecule)
        assert str_molecule == str(qcinp)

    def test_read_cdft(self):
        str_cdft = """Once again, I'm trying to break you!

$cdft
   1.0
   1.0 1 27
   0.0
   1.0 1 27 s
--------------
   0.0
   1.0 1 27
   -1.0
   1.0 1 27 s
$end
"""

        result = [
            [
                {"value": 1.0, "coefficients": [1.0], "first_atoms": [1], "last_atoms": [27], "types": [None]},
                {"value": 0.0, "coefficients": [1.0], "first_atoms": [1], "last_atoms": [27], "types": ["s"]},
            ],
            [
                {"value": 0.0, "coefficients": [1.0], "first_atoms": [1], "last_atoms": [27], "types": [None]},
                {"value": -1.0, "coefficients": [1.0], "first_atoms": [1], "last_atoms": [27], "types": ["s"]},
            ],
        ]

        parsed = QCInput.read_cdft(str_cdft)
        assert parsed[0][0] == result[0][0]
        assert parsed[0][1] == result[0][1]
        assert parsed[1][0] == result[1][0]
        assert parsed[1][1] == result[1][1]

    def test_read_almo(self):
        str_almo = """I remain resolute in trying to break you!

$almo_coupling
   1 2
   0 1
   --
   0 1
   1 2
$end"""

        result = [[(1, 2), (0, 1)], [(0, 1), (1, 2)]]

        assert QCInput.read_almo(str_almo) == result

    def test_write_file_from_OptSet(self):
        from pymatgen.io.qchem.sets import OptSet

        odd_dict = loadfn(f"{module_dir}/odd.json")
        odd_mol = odd_dict["spec"]["_tasks"][0]["molecule"]
        qcinp = OptSet(odd_mol)
        qcinp.write_file(f"{module_dir}/test.qin")
        test_path = f"{module_dir}/test.qin"
        ref_path = f"{module_dir}/test_ref.qin"

        with open(ref_path) as ref_file, open(test_path) as test_file:
            for l_test, l_ref in zip(test_file, ref_file):
                # By default, if this statement fails the offending line will be printed
                assert l_test == l_ref

        os.remove(f"{module_dir}/test.qin")

    def test_write_file_from_OptSet_with_vdw(self):
        from pymatgen.io.qchem.sets import OptSet

        odd_dict = loadfn(f"{module_dir}/odd.json")
        odd_mol = odd_dict["spec"]["_tasks"][0]["molecule"]
        qcinp = OptSet(odd_mol, overwrite_inputs={"van_der_waals": {"16": 3.14159}})
        qcinp.write_file(f"{module_dir}/test_vdw.qin")
        test_path = f"{module_dir}/test_vdw.qin"
        ref_path = f"{module_dir}/test_ref_vdw.qin"

        with open(ref_path) as ref_file, open(test_path) as test_file:
            for l_test, l_ref in zip(test_file, ref_file):
                # By default, if this statement fails the offending line will be printed
                assert l_test == l_ref

        os.remove(f"{module_dir}/test_vdw.qin")

    def test_read_write_nbo7(self):
        test_path = f"{module_dir}/test_nbo7.qin"
        ref_path = f"{TEST_FILES_DIR}/molecules/new_qchem_files/nbo7.qin"
        qcinp = QCInput.from_file(f"{TEST_FILES_DIR}/molecules/new_qchem_files/nbo7.qin")
        qcinp.write_file(test_path)

        with open(test_path) as ref_file, open(ref_path) as test_file:
            for l_test, l_ref in zip(test_file, ref_file):
                # By default, if this statement fails the offending line will be printed
                assert l_test == l_ref

        os.remove(test_path)

    def test_read_write_nbo_e2pert(self):
        qcinp = QCInput.from_file(f"{TEST_FILES_DIR}/molecules/new_qchem_files/e2pert.qin")
        qcinp.write_file(f"{module_dir}/test_e2pert.qin")
        test_path = f"{TEST_FILES_DIR}/molecules/new_qchem_files/e2pert.qin"
        ref_path = f"{module_dir}/test_e2pert.qin"

        with open(ref_path) as ref_file, open(test_path) as test_file:
            for l_test, l_ref in zip(test_file, ref_file):
                assert l_test == l_ref

        os.remove(f"{module_dir}/test_e2pert.qin")

    def test_read_write_custom_smd(self):
        qcinp = QCInput.from_file(f"{TEST_FILES_DIR}/molecules/new_qchem_files/custom_smd.qin")
        qcinp.write_file(f"{module_dir}/test_custom_smd.qin")
        test_path = f"{TEST_FILES_DIR}/molecules/new_qchem_files/custom_smd.qin"
        ref_path = f"{module_dir}/test_custom_smd.qin"

        with open(ref_path) as ref_file, open(test_path) as test_file:
            for l_test, l_ref in zip(test_file, ref_file):
                assert l_test == l_ref

        os.remove(f"{module_dir}/test_custom_smd.qin")
