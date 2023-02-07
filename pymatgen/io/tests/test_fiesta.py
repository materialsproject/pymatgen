from __future__ import annotations

import os
import unittest

from pymatgen.core.structure import Molecule
from pymatgen.io.fiesta import FiestaInput, FiestaOutput
from pymatgen.util.testing import PymatgenTest


class FiestaInputTest(PymatgenTest):
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
        self.cellin = FiestaInput(
            mol,
            correlation_grid={"dE_grid": "0.500", "n_grid": "14"},
            Exc_DFT_option={"rdVxcpsi": "1"},
            COHSEX_options={
                "eigMethod": "C",
                "mix_cohsex": "0.500",
                "nc_cohsex": "0",
                "nit_cohsex": "0",
                "nv_cohsex": "0",
                "resMethod": "V",
                "scf_cohsex_wf": "0",
            },
            GW_options={"nc_corr": "10", "nit_gw": "3", "nv_corr": "10"},
            BSE_TDDFT_options={
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
        cellin = FiestaInput(mol)
        assert cellin.molecule.spin_multiplicity == 1

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
        assert str(self.cellin) == ans
        cellin = FiestaInput.from_string(ans)
        assert cellin.GW_options["nc_corr"] == "10"
        assert cellin.COHSEX_options["eigMethod"] == "C"


class FiestaOutputTest(PymatgenTest):
    def setUp(self):
        self.logfiesta = FiestaOutput(os.path.join(PymatgenTest.TEST_FILES_DIR, "log_fiesta"))

    def test_props(self):
        out = self.logfiesta
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


if __name__ == "__main__":
    unittest.main()
