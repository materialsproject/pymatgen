from __future__ import annotations

import unittest

from pymatgen.core.structure import Molecule
from pymatgen.io.fiesta import FiestaInput, FiestaOutput
from pymatgen.util.testing import TEST_FILES_DIR


class TestFiestaInput(unittest.TestCase):
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
        self.cell_in = FiestaInput(
            mol,
            correlation_grid={"dE_grid": "0.500", "n_grid": "14"},
            exc_dft_option={"rdVxcpsi": "1"},
            cohsex_options={
                "eigMethod": "C",
                "mix_cohsex": "0.500",
                "nc_cohsex": "0",
                "nit_cohsex": "0",
                "nv_cohsex": "0",
                "resMethod": "V",
                "scf_cohsex_wf": "0",
            },
            gw_options={"nc_corr": "10", "nit_gw": "3", "nv_corr": "10"},
            bse_tddft_options={
                "do_bse": "1",
                "do_tddft": "0",
                "nc_bse": "382",
                "nit_bse": "50",
                "npsi_bse": "1",
                "nv_bse": "21",
            },
        )

    def test_init(self):
        mol = Molecule(["C", "H", "H", "H", "H"], self.coords)
        cell_in = FiestaInput(mol)
        assert cell_in.molecule.spin_multiplicity == 1

    def test_str_and_from_string(self):
        ans = (
            "# number of atoms and species\n   5    2\n# number of valence bands\n    5\n"
            "# number of points and spacing in eV for correlation grid\n    14    0.500\n"
            "# relire=1 ou recalculer=0 Exc DFT\n    1\n"
            "# number of COHSEX corrected occp and unoccp bands: C=COHSEX  H=HF\n    0    0   C\n"
            "# number of COHSEX iter, scf on wfns, mixing coeff; V=RI-V  I=RI-D\n    0   V       0       0.500\n"
            "# number of GW corrected occp and unoccp bands\n   10   10\n# number of GW iterations\n    3\n"
            "# dumping for BSE and TDDFT\n    1    0\n"
            "# number of occp. and virtual bands of BSE: nocore and up to 40 eVs\n    21   382\n"
            "# number of excitations needed and number of iterations\n    1   50\n"
            "# list of symbols in order\nC\nH\n"
            "# scaling factor\n    1.000\n# atoms x,y,z cartesian .. will be multiplied by scale\n 0.0 0.0 0.0 1\n"
            " 0.0 0.0 1.089 2\n 1.026719 0.0 -0.363 2\n -0.51336 -0.889165 -0.363 2\n -0.51336 0.889165 -0.363 2"
            "\n            "
        )
        assert str(self.cell_in) == ans
        cell_in = FiestaInput.from_str(ans)
        assert cell_in.GW_options["nc_corr"] == "10"
        assert cell_in.cohsex_options["eigMethod"] == "C"


class TestFiestaOutput(unittest.TestCase):
    def setUp(self):
        self.log_fiesta = FiestaOutput(f"{TEST_FILES_DIR}/log_fiesta")

    def test_props(self):
        out = self.log_fiesta
        assert out.data[0]["Gaps"]["Egap_QP_Linear"] == "10.4135"
        assert out.data[0]["HOMO"] == {
            "band": "HOMO",
            "eKS": "-7.3029",
            "eQP_Linear": "-9.5142",
            "eQP_SCF": "-8.9264",
            "eQP_old": "-7.7188",
            "eXX": "-15.9483",
            "sigma_c_Linear": "-0.4587",
            "sigma_c_SCF": "0.3900",
            "z": "0.87",
        }
