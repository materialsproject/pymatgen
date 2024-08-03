from __future__ import annotations

import pytest
from pytest import approx

from pymatgen.analysis.structure_prediction.volume_predictor import DLSVolumePredictor, RLSVolumePredictor
from pymatgen.core import Structure
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest

TEST_DIR = f"{TEST_FILES_DIR}/analysis/structure_prediction"


class TestRLSVolumePredictor(PymatgenTest):
    def test_predict(self):
        struct = PymatgenTest.get_structure("CsCl")
        nacl = PymatgenTest.get_structure("CsCl")
        nacl.replace_species({"Cs": "Na"})
        nacl.scale_lattice(184.384551033)
        predictor = RLSVolumePredictor(radii_type="ionic")
        assert predictor.predict(struct, nacl) == approx(342.84905395082535)
        predictor = RLSVolumePredictor(radii_type="atomic")
        assert predictor.predict(struct, nacl) == approx(391.884366481)
        lif = PymatgenTest.get_structure("CsCl")
        lif.replace_species({"Cs": "Li", "Cl": "F"})
        predictor = RLSVolumePredictor(radii_type="ionic")
        assert predictor.predict(lif, nacl) == approx(74.268402413690467)
        predictor = RLSVolumePredictor(radii_type="atomic")
        assert predictor.predict(lif, nacl) == approx(62.2808125839)

        lfpo = PymatgenTest.get_structure("LiFePO4")
        lmpo = PymatgenTest.get_structure("LiFePO4")
        lmpo.replace_species({"Fe": "Mn"})
        predictor = RLSVolumePredictor(radii_type="ionic")
        assert predictor.predict(lmpo, lfpo) == approx(310.08253254420134)
        predictor = RLSVolumePredictor(radii_type="atomic")
        assert predictor.predict(lmpo, lfpo) == approx(299.607967711)

        sto = PymatgenTest.get_structure("SrTiO3")
        scoo = PymatgenTest.get_structure("SrTiO3")
        scoo.replace_species({"Ti4+": "Co4+"})
        predictor = RLSVolumePredictor(radii_type="ionic")
        assert predictor.predict(scoo, sto) == approx(56.162534974936463)
        predictor = RLSVolumePredictor(radii_type="atomic")
        assert predictor.predict(scoo, sto) == approx(57.4777835108)

        # Use Ag7P3S11 as a test case:

        # (i) no oxidation states are assigned and CVP-atomic scheme is selected.
        aps = Structure.from_file(f"{TEST_DIR}/Ag7P3S11_mp-683910_primitive.cif")
        apo = Structure.from_file(f"{TEST_DIR}/Ag7P3S11_mp-683910_primitive.cif")
        apo.replace_species({"S": "O"})
        predictor = RLSVolumePredictor(radii_type="atomic", check_isostructural=False)
        assert predictor.predict(apo, aps) == approx(1196.31384276)

        # (ii) Oxidation states are assigned.
        apo.add_oxidation_state_by_element({"Ag": 1, "P": 5, "O": -2})
        aps.add_oxidation_state_by_element({"Ag": 1, "P": 5, "S": -2})
        predictor = RLSVolumePredictor(radii_type="ionic")
        assert predictor.predict(apo, aps) == approx(1165.23259079)
        predictor = RLSVolumePredictor(radii_type="atomic")
        assert predictor.predict(apo, aps) == approx(1196.31384276)

    def test_modes(self):
        cs_cl = PymatgenTest.get_structure("CsCl")
        na_cl = PymatgenTest.get_structure("CsCl")
        na_cl.replace_species({"Cs": "Na"})
        na_cl.scale_lattice(184.384551033)
        vol_pred = RLSVolumePredictor(radii_type="ionic", use_bv=False)
        with pytest.raises(ValueError, match="Cannot find volume scaling based on radii choices specified"):
            vol_pred.predict(cs_cl, na_cl)
        vol_pred = RLSVolumePredictor(radii_type="ionic-atomic", use_bv=False)
        assert vol_pred.predict(cs_cl, na_cl) == approx(391.884366481)
        vol_pred = RLSVolumePredictor(radii_type="ionic-atomic", use_bv=True)
        assert vol_pred.predict(cs_cl, na_cl) == approx(342.84905395082535)


class TestDLSVolumePredictor(PymatgenTest):
    def test_predict(self):
        vol_pred = DLSVolumePredictor()
        p_fast = DLSVolumePredictor(cutoff=0.0)  # for speed on compressed cells
        p_nolimit = DLSVolumePredictor(min_scaling=None, max_scaling=None)  # no limits on scaling

        fen = Structure.from_file(f"{TEST_DIR}/FeN_mp-6988.cif")

        assert vol_pred.predict(fen) == approx(18.2252568873)
        fen.scale_lattice(fen.volume * 3.0)
        assert p_nolimit.predict(fen) == approx(18.2252568873)
        assert vol_pred.predict(fen) == approx(fen.volume * 0.5)
        fen.scale_lattice(fen.volume * 0.1)
        assert p_nolimit.predict(fen) == approx(18.2252568873)
        assert vol_pred.predict(fen) == approx(fen.volume * 1.5)
        assert p_fast.predict(fen) == approx(fen.volume * 1.5)

        lfpo = PymatgenTest.get_structure("LiFePO4")

        lfpo.scale_lattice(lfpo.volume * 3.0)
        assert p_nolimit.predict(lfpo) == approx(291.62094410192924)
        assert vol_pred.predict(lfpo) == approx(lfpo.volume * 0.5)
        lfpo.scale_lattice(lfpo.volume * 0.1)
        assert p_nolimit.predict(lfpo) == approx(291.62094410192924)
        assert vol_pred.predict(lfpo) == approx(lfpo.volume * 1.5)
        assert p_fast.predict(lfpo) == approx(lfpo.volume * 1.5)

        lmpo = PymatgenTest.get_structure("LiFePO4")
        lmpo.replace_species({"Fe": "Mn"})
        assert vol_pred.predict(lmpo) == approx(290.795329052)
