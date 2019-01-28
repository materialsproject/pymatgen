# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import unittest
import os
from pymatgen.analysis.structure_prediction.volume_predictor import DLSVolumePredictor, \
    RLSVolumePredictor
from pymatgen.util.testing import PymatgenTest
from pymatgen.core import Structure
import warnings

dir_path = os.path.join(os.path.dirname(os.path.abspath(__file__)))


class RLSVolumePredictorTest(PymatgenTest):

    def setUp(self):
        warnings.filterwarnings("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_predict(self):
        s = PymatgenTest.get_structure("CsCl")
        nacl = PymatgenTest.get_structure("CsCl")
        nacl.replace_species({"Cs": "Na"})
        nacl.scale_lattice(184.384551033)
        p = RLSVolumePredictor(radii_type="ionic")
        self.assertAlmostEqual(p.predict(s, nacl), 342.84905395082535)
        p = RLSVolumePredictor(radii_type="atomic")
        self.assertAlmostEqual(p.predict(s, nacl), 391.884366481)
        lif = PymatgenTest.get_structure("CsCl")
        lif.replace_species({"Cs": "Li", "Cl": "F"})
        p = RLSVolumePredictor(radii_type="ionic")
        self.assertAlmostEqual(p.predict(lif, nacl), 74.268402413690467)
        p = RLSVolumePredictor(radii_type="atomic")
        self.assertAlmostEqual(p.predict(lif, nacl), 62.2808125839)

        lfpo = PymatgenTest.get_structure("LiFePO4")
        lmpo = PymatgenTest.get_structure("LiFePO4")
        lmpo.replace_species({"Fe": "Mn"})
        p = RLSVolumePredictor(radii_type="ionic")
        self.assertAlmostEqual(p.predict(lmpo, lfpo), 310.08253254420134)
        p = RLSVolumePredictor(radii_type="atomic")
        self.assertAlmostEqual(p.predict(lmpo, lfpo), 299.607967711)

        sto = PymatgenTest.get_structure("SrTiO3")
        scoo = PymatgenTest.get_structure("SrTiO3")
        scoo.replace_species({"Ti4+": "Co4+"})
        p = RLSVolumePredictor(radii_type="ionic")
        self.assertAlmostEqual(p.predict(scoo, sto), 56.162534974936463)
        p = RLSVolumePredictor(radii_type="atomic")
        self.assertAlmostEqual(p.predict(scoo, sto), 57.4777835108)

        # Use Ag7P3S11 as a test case:

        # (i) no oxidation states are assigned and CVP-atomic scheme is selected.
        aps = Structure.from_file(os.path.join(dir_path,
                                               "Ag7P3S11_mp-683910_primitive.cif"))
        apo = Structure.from_file(os.path.join(dir_path,
                                               "Ag7P3S11_mp-683910_primitive.cif"))
        apo.replace_species({"S": "O"})
        p = RLSVolumePredictor(radii_type="atomic", check_isostructural=False)
        self.assertAlmostEqual(p.predict(apo, aps), 1196.31384276)

        # (ii) Oxidation states are assigned.
        apo.add_oxidation_state_by_element({"Ag": 1, "P": 5, "O": -2})
        aps.add_oxidation_state_by_element({"Ag": 1, "P": 5, "S": -2})
        p = RLSVolumePredictor(radii_type="ionic")
        self.assertAlmostEqual(p.predict(apo, aps), 1165.23259079)
        p = RLSVolumePredictor(radii_type="atomic")
        self.assertAlmostEqual(p.predict(apo, aps), 1196.31384276)

    def test_modes(self):
        s = PymatgenTest.get_structure("CsCl")
        nacl = PymatgenTest.get_structure("CsCl")
        nacl.replace_species({"Cs": "Na"})
        nacl.scale_lattice(184.384551033)
        p = RLSVolumePredictor(radii_type="ionic", use_bv=False)
        self.assertRaises(ValueError, p.predict, s, nacl)
        p = RLSVolumePredictor(radii_type="ionic-atomic", use_bv=False)
        self.assertAlmostEqual(p.predict(s, nacl), 391.884366481)
        p = RLSVolumePredictor(radii_type="ionic-atomic", use_bv=True)
        self.assertAlmostEqual(p.predict(s, nacl), 342.84905395082535)


class DLSVolumePredictorTest(PymatgenTest):

    def test_predict(self):
        p = DLSVolumePredictor()
        p_fast = DLSVolumePredictor(cutoff=0.0)  # for speed on compressed cells

        fen = Structure.from_file(os.path.join(dir_path, "FeN_mp-6988.cif"))
        self.assertAlmostEqual(p.predict(fen), 18.2252568873)
        fen.scale_lattice(3.0)
        self.assertAlmostEqual(p.predict(fen), 18.2252568873)
        fen.scale_lattice(0.24)
        self.assertAlmostEqual(p.predict(fen), 18.2252568873)

        lfpo = PymatgenTest.get_structure("LiFePO4")
        lfpo.scale_lattice(10.1)
        self.assertAlmostEqual(p.predict(lfpo), 291.62094410192924)
        lfpo.scale_lattice(0.2)
        self.assertAlmostEqual(p_fast.predict(lfpo), 291.62094410192924)
        lmpo = PymatgenTest.get_structure("LiFePO4")
        lmpo.replace_species({"Fe": "Mn"})
        self.assertAlmostEqual(p.predict(lmpo), 290.795329052)


if __name__ == '__main__':
    unittest.main()
