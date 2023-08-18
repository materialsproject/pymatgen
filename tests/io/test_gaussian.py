from __future__ import annotations

import unittest

import pytest
from pytest import approx

from pymatgen.core.structure import Molecule
from pymatgen.electronic_structure.core import Spin
from pymatgen.io.gaussian import GaussianInput, GaussianOutput
from pymatgen.util.testing import TEST_FILES_DIR

test_dir = f"{TEST_FILES_DIR}/molecules"


class TestGaussianInput(unittest.TestCase):
    def setUp(self):
        coords = [
            [0.000000, 0.000000, 0.000000],
            [0.000000, 0.000000, 1.089000],
            [1.026719, 0.000000, -0.363000],
            [-0.513360, -0.889165, -0.363000],
            [-0.513360, 0.889165, -0.363000],
        ]
        self.coords = coords
        mol = Molecule(["C", "H", "H", "H", "H"], coords)
        self.gau = GaussianInput(
            mol,
            route_parameters={"SP": "", "SCF": "Tight"},
            input_parameters={"EPS": 12},
        )

    def test_init(self):
        mol = Molecule(["C", "H", "H", "H", "H"], self.coords)
        gau = GaussianInput(mol, charge=1, route_parameters={"SP": "", "SCF": "Tight"})
        assert gau.spin_multiplicity == 2
        mol = Molecule(["C", "H", "H", "H", "H"], self.coords, charge=-1)
        gau = GaussianInput(mol, route_parameters={"SP": "", "SCF": "Tight"})
        assert gau.spin_multiplicity == 2
        with pytest.raises(
            ValueError, match="Charge of -1 and spin multiplicity of 1 is not possible for this molecule"
        ):
            GaussianInput(mol, spin_multiplicity=1)

    def test_str_and_from_string(self):
        answer = """#P HF/6-31G(d) SCF=Tight SP

H4 C1

0 1
C
H 1 B1
H 1 B2 2 A2
H 1 B3 2 A3 3 D3
H 1 B4 2 A4 4 D4

B1=1.089000
B2=1.089000
A2=109.471221
B3=1.089000
A3=109.471213
D3=120.000017
B4=1.089000
A4=109.471213
D4=119.999966

EPS=12

"""
        assert str(self.gau) == answer
        gau = GaussianInput.from_str(answer)
        assert gau.functional == "HF"
        assert gau.input_parameters["EPS"] == "12"

    def test_from_cart_coords(self):
        answer = """#P HF/6-31G(d) SCF=Tight SP

H4 C1

0 1
C 0.000000 0.000000 0.000000
H 0.000000 0.000000 1.089000
H 1.026719 0.000000 -0.363000
H -0.513360 -0.889165 -0.363000
H -0.513360 0.889165 -0.363000

EPS=12

"""
        assert self.gau.to_str(cart_coords=True) == answer
        gau = GaussianInput.from_str(answer)
        assert gau.functional == "HF"
        assert gau.charge == 0
        assert gau.spin_multiplicity == 1
        assert gau.input_parameters["EPS"] == "12"

    def test_from_file(self):
        filepath = f"{test_dir}/MethylPyrrolidine_drawn.gjf"
        gau = GaussianInput.from_file(filepath)
        assert gau.molecule.composition.formula == "H11 C5 N1"
        assert "opt" in gau.route_parameters
        assert gau.route_parameters["geom"] == "connectivity"
        assert gau.functional == "b3lyp"
        assert gau.basis_set == "6-311+g(d,p)"
        filepath = f"{test_dir}/g305_hb.txt"
        with open(filepath) as f:
            txt = f.read()
        toks = txt.split("--link1--")
        for idx, tok in enumerate(toks):
            lines = [line.strip() for line in tok.strip().split("\n")]
            gau = GaussianInput.from_str("\n".join(lines))
            assert gau.molecule is not None
            if idx == 0:
                mol = gau.molecule
        answer = """Full Formula (H4 O2)
Reduced Formula: H2O
Charge = 0, Spin Mult = 1
Sites (6)
0 O     0.000000     0.000000     0.000000
1 O     0.000000     0.000000     2.912902
2 H     0.892596     0.000000    -0.373266
3 H     0.143970     0.000219     0.964351
4 H    -0.582554     0.765401     3.042783
5 H    -0.580711    -0.766761     3.043012"""
        assert str(mol) == answer

    def test_from_string(self):
        gau_str = """%mem=5000000
        %chk=filename
        # mp2/6-31g* scf=direct
        opt freq

        SIH4+ H2---SIH2+ CS //MP2(full)/6-31G* MP2=-290.9225259

        1,2
        Si
        X,1,1.
        H,1,R1,2,HALF1
        H,1,R1,2,HALF1,3,180.,0
        X,1,1.,2,90.,3,90.,0
        X,1,1.,5,THETA,2,180.,0
        H,1,R3,6,HALF3,5,0.,0
        H,1,R4,6,HALF3,7,180.,0

        R1=1.47014
        R3=1.890457
        R4=1.83514
        HALF1=60.633314
        THETA=10.35464
        HALF3=11.861807"""

        gau = GaussianInput.from_str(gau_str)
        assert gau.molecule.composition.reduced_formula == "X3SiH4"
        assert set(gau.route_parameters) == {"opt", "freq", "scf"}

    def test_gen_basis(self):
        gau_str = """#N B3LYP/Gen Pseudo=Read

Test

0 1
C
H 1 B1
H 1 B2 2 A2
H 1 B3 2 A3 3 D3
H 1 B4 2 A4 4 D4

B1=1.089000
B2=1.089000
A2=109.471221
B3=1.089000
A3=109.471213
D3=120.000017
B4=1.089000
A4=109.471213
D4=119.999966

C 0
6-31G(d,p)
****
H 0
6-31G
****



"""
        mol = Molecule(["C", "H", "H", "H", "H"], self.coords)
        gen_basis = "C 0\n6-31G(d,p)\n****\nH 0\n6-31G\n****"
        gau = GaussianInput(
            mol,
            functional="B3LYP",
            gen_basis=gen_basis,
            dieze_tag="#N",
            route_parameters={"Pseudo": "Read"},
            title="Test",
        )
        assert gau.to_str(cart_coords=False) == gau_str

    def test_multiple_paramaters(self):
        """
        This test makes sure that input files with multi-parameter keywords
        and route cards with multiple lines can be parsed accurately.
        """
        filepath = f"{test_dir}/l-cysteine.inp"
        route = {
            "test": None,
            "integral": {"grid": "UltraFine"},
            "opt": {"Z-Matrix": None, "maxcycles": "80", "tight": None},
        }
        gin = GaussianInput.from_file(filepath)
        assert gin.dieze_tag == "#n"
        assert gin.functional == "B3LYP"
        assert gin.basis_set == "6-31+G**"
        assert gin.route_parameters == route
        assert gin.title == "L-cysteine neutral"
        assert gin.charge == 0
        assert gin.spin_multiplicity == 1

    def test_no_molecule(self):
        """Test that we can write input files without a geometry."""
        # Makes a file without geometry
        input_file = GaussianInput(None, charge=0, spin_multiplicity=2)
        input_str = input_file.to_str().strip()
        assert input_str.endswith("0 2")

    def test_no_molecule_func_bset_charge_mult(self):
        """
        Test that we can write inputs files without a geometry,
        functional, basis set, charge or multiplicity
        (mainly for postprocessing jobs where this info is read from .chk).
        """
        gau_str = "#P chkbasis geom=allcheck guess=(only,read) pop=naturalorbital\n"
        gau_str += "\n"
        gau_str += "Restart"

        input_file = GaussianInput(
            None,
            charge=None,
            spin_multiplicity=None,
            functional=None,
            basis_set=None,
            route_parameters={
                "chkbasis": None,
                "geom": "allcheck",
                "guess": {"only": None, "read": None},
                "pop": "naturalorbital",
            },
        )
        input_str = input_file.to_str().strip()
        assert input_str == gau_str


class TestGaussianOutput(unittest.TestCase):
    # TODO: Add unittest for PCM type output.

    def setUp(self):
        self.gauout = GaussianOutput(f"{test_dir}/methane.log")

    def test_resume(self):
        resume = self.gauout.resumes[0]
        methane_resume = r"""1\1\GINC-SHYUE-LAPTOP\FOpt\RHF\3-21G\C1H4\SHYUE\27-Feb-2008\0\\#p hf/3
        -21G opt\\Title Card Required\\0,1\C,0.,0.,0.\H,0.,0.,1.0829014152\H,1
        .0209692454,0.,-0.3609671384\H,-0.5104846227,-0.884185303,-0.360967138
        4\H,-0.5104846227,0.884185303,-0.3609671384\\Version=IA32L-G03RevD.01\
        State=1-A1\HF=-39.9768776\RMSD=3.210e-09\RMSF=5.014e-08\Thermal=0.\Dip
        ole=0.,0.,0.\PG=TD [O(C1),4C3(H1)]\\@"""
        methane_resume = "".join(r.strip() for r in methane_resume.split("\n"))

        assert resume == methane_resume

    def test_props(self):
        gau = self.gauout
        assert len(gau.energies) == 3
        assert gau.energies[-1] == approx(-39.9768775602)
        assert len(gau.structures) == 4
        for mol in gau.structures:
            assert mol.formula == "H4 C1"
        assert "opt" in gau.route_parameters
        assert gau.stationary_type == "Minimum"
        assert gau.functional == "hf"
        assert gau.basis_set == "3-21G"
        assert gau.num_basis_func == 17
        d = gau.as_dict()
        assert d["input"]["functional"] == "hf"
        assert d["output"]["final_energy"] == approx(-39.9768775602)
        assert len(gau.cart_forces) == 3
        assert gau.cart_forces[0][5] == 0.009791094
        assert gau.cart_forces[0][-1] == -0.003263698
        assert gau.cart_forces[2][-1] == -0.000000032
        assert gau.eigenvalues[Spin.up][-1] == 1.95586
        assert gau.num_basis_func == 17
        assert gau.is_spin is False

        ch2o_co2 = GaussianOutput(f"{test_dir}/CH2O_CO2.log")
        assert len(ch2o_co2.frequencies) == 2
        assert len(ch2o_co2.frequencies[0]) == 6
        assert len(ch2o_co2.frequencies[1]) == 4
        assert ch2o_co2.frequencies[0][0]["frequency"] == 1203.1940
        assert ch2o_co2.frequencies[0][0]["symmetry"] == 'A"'
        assert ch2o_co2.frequencies[0][3]["IR_intensity"] == 60.9575
        assert ch2o_co2.frequencies[0][3]["r_mass"] == 3.7543
        assert ch2o_co2.frequencies[0][4]["f_constant"] == 5.4175
        assert ch2o_co2.frequencies[0][1]["mode"] == [
            0.15,
            0.00,
            0.00,
            -0.26,
            0.65,
            0.00,
            -0.26,
            -0.65,
            0.00,
            -0.08,
            0.00,
            0.00,
        ]
        assert ch2o_co2.frequencies[1][3]["mode"] == [0.00, 0.00, 0.88, 0.00, 0.00, -0.33, 0.00, 0.00, -0.33]
        assert ch2o_co2.frequencies[1][3]["symmetry"] == "SGU"
        assert ch2o_co2.eigenvalues[Spin.up][3] == -1.18394

        h2o = GaussianOutput(f"{test_dir}/H2O_gau_vib.out")
        assert len(h2o.frequencies[0]) == 3
        assert h2o.frequencies[0][0]["frequency"] == 1662.8033
        assert h2o.frequencies[0][1]["symmetry"] == "A'"
        assert h2o.hessian[0, 0] == 0.356872
        assert h2o.hessian.shape == (9, 9)
        assert h2o.hessian[8, :].tolist() == [
            -0.143692e-01,
            0.780136e-01,
            -0.362637e-01,
            -0.176193e-01,
            0.277304e-01,
            -0.583237e-02,
            0.319885e-01,
            -0.105744e00,
            0.420960e-01,
        ]

    def test_pop(self):
        gau = GaussianOutput(f"{test_dir}/H2O_gau.out")
        assert gau.num_basis_func == 13
        assert gau.electrons == (5, 5)
        assert gau.is_spin
        assert gau.eigenvalues[Spin.down] == [
            -20.55343,
            -1.35264,
            -0.72655,
            -0.54824,
            -0.49831,
            0.20705,
            0.30297,
            1.10569,
            1.16144,
            1.16717,
            1.20460,
            1.38903,
            1.67608,
        ]
        mo = gau.molecular_orbital
        assert len(mo) == 2  # la 6
        assert len(mo[Spin.down]) == 13
        assert len(mo[Spin.down][0]) == 3
        assert mo[Spin.down][5][0]["1S"] == -0.08771
        assert mo[Spin.down][5][0]["2PZ"] == -0.21625
        assert gau.eigenvectors[Spin.up][:, 5].tolist() == [
            -0.08771,
            0.10840,
            0.00000,
            0.00000,
            -0.21625,
            1.21165,
            0.00000,
            0.00000,
            -0.44481,
            -0.06348,
            -1.00532,
            -0.06348,
            -1.00532,
        ]

        assert gau.atom_basis_labels[0] == ["1S", "2S", "2PX", "2PY", "2PZ", "3S", "3PX", "3PY", "3PZ"]
        assert gau.atom_basis_labels[2] == ["1S", "2S"]

        gau = GaussianOutput(f"{test_dir}/H2O_gau_vib.out")

        assert gau.bond_orders[(0, 1)] == 0.7582
        assert gau.bond_orders[(1, 2)] == 0.0002

    def test_scan(self):
        gau = GaussianOutput(f"{test_dir}/so2_scan.log")
        d = gau.read_scan()
        assert approx(d["energies"][-1]) == -548.02102
        assert len(d["coords"]) == 1
        assert len(d["energies"]) == len(gau.energies)
        assert len(d["energies"]) == 21
        gau = GaussianOutput(f"{test_dir}/so2_scan_opt.log")
        assert len(gau.opt_structures) == 21
        d = gau.read_scan()
        assert approx(d["energies"][-1]) == -548.02336
        assert len(d["coords"]) == 2
        assert len(d["energies"]) == 21
        assert approx(d["coords"]["DSO"][6]) == 1.60000
        assert approx(d["coords"]["ASO"][2]) == 124.01095
        gau = GaussianOutput(f"{test_dir}/H2O_scan_G16.out")
        assert len(gau.opt_structures) == 21
        coords = [
            [0.000000, 0.000000, 0.094168],
            [0.000000, 0.815522, -0.376673],
            [0.000000, -0.815522, -0.376673],
        ]
        assert gau.opt_structures[-1].cart_coords.tolist() == coords
        d = gau.read_scan()
        assert approx(d["energies"][-1]) == -0.00523
        assert len(d["coords"]) == 3
        assert len(d["energies"]) == 21
        assert approx(d["coords"]["R1"][6]) == 0.94710
        assert approx(d["coords"]["R2"][17]) == 0.94277

    def test_geo_opt(self):
        """Test an optimization where no "input orientation" is outputted."""
        gau = GaussianOutput(f"{test_dir}/acene-n_gaussian09_opt.out")
        assert approx(gau.energies[-1]) == -1812.58399675
        assert len(gau.structures) == 6
        # Test the first 3 atom coordinates
        coords = [
            [-13.642932, 0.715060, 0.000444],
            [-13.642932, -0.715060, 0.000444],
            [-12.444202, 1.416837, 0.000325],
        ]
        assert gau.opt_structures[-1].cart_coords[:3].tolist() == coords

    def test_td(self):
        gau = GaussianOutput(f"{test_dir}/so2_td.log")
        transitions = gau.read_excitation_energies()
        assert len(transitions) == 4
        assert transitions[0] == approx((3.9281, 315.64, 0.0054))

    def test_multiple_parameters(self):
        """
        This test makes sure that input files with multi-parameter keywords
        and route cards with multiple lines can be parsed accurately.
        """
        filepath = f"{test_dir}/l-cysteine.out"
        route = {
            "test": None,
            "integral": {"grid": "UltraFine"},
            "opt": {"Z-Matrix": None, "maxcycles": "80", "tight": None},
        }
        gout = GaussianOutput(filepath)
        assert gout.dieze_tag == "#n"
        assert gout.functional == "B3LYP"
        assert gout.basis_set == "6-31+G**"
        assert gout.route_parameters == route
        assert gout.title == "L-cysteine neutral"
        assert gout.charge == 0
        assert gout.spin_multiplicity == 1

    def test_multiple_parameters_with_multiple_completed_lines(self):
        """
        This test makes sure that input files with multi-parameter keywords
        and route cards with multiple completed lines which are split by line break parse correctly.
        """
        filepath = f"{test_dir}/EC.log.gz"
        route_params = {
            "opt": {"loose": None, "maxcyc": "400"},
            "freq": None,
            "SCRF": "(SMD,READ)",
            "EmpiricalDispersion": "GD3BJ",
        }
        gout = GaussianOutput(filepath)
        assert gout.dieze_tag == "#N"
        assert gout.functional == "B3LYP"
        assert gout.basis_set == "6-311++G(2d,p)"
        assert gout.route_parameters == route_params
        assert gout.title == "H4 C3 O3"
        assert gout.charge == 0
        assert gout.spin_multiplicity == 1
