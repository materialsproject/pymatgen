from __future__ import annotations

import os
import unittest

import pytest
from pytest import approx

from pymatgen.apps.battery.analyzer import BatteryAnalyzer
from pymatgen.core.structure import Structure
from pymatgen.util.testing import PymatgenTest


class BatteryAnalyzerTest(PymatgenTest):
    def load_from_cif(self, filename, oxidations, working_ion="Li"):
        s = Structure.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, filename))
        s.add_oxidation_state_by_element(oxidations)
        return BatteryAnalyzer(s, working_ion)

    def load_from_internal(self, name, oxidations, working_ion="Li"):
        s = self.get_structure(name).copy()
        s.add_oxidation_state_by_element(oxidations)
        return BatteryAnalyzer(s, working_ion)

    def setUp(self):
        self.lifepo4 = self.load_from_internal("LiFePO4", {"Li": 1, "Fe": 2, "P": 5, "O": -2})
        self.nafepo4 = self.load_from_internal("NaFePO4", {"Na": 1, "Fe": 2, "P": 5, "O": -2}, working_ion="Na")
        self.la2coo4f = self.load_from_internal("La2CoO4F", {"La": 3, "Co": 3, "O": -2, "F": -1}, working_ion="F")
        self.fepo4 = self.load_from_cif("FePO4a.cif", {"Fe": 3, "P": 5, "O": -2})
        self.la2coo4 = self.load_from_cif("La2CoO4.cif", {"La": 3, "Co": 2, "O": -2}, working_ion="F")
        self.lifemnpo4 = self.load_from_cif("Li4Fe3Mn1(PO4)4.cif", {"Li": 1, "Fe": 2, "Mn": 2, "P": 5, "O": -2})
        self.li8nicofe208 = self.load_from_cif(
            "Li8Fe2NiCoO8.cif", {"Li": 1, "Fe": 2, "Mn": 2, "Co": 2, "Ni": 2, "O": -2}
        )
        self.li3v2p3o12 = self.load_from_internal("Li3V2(PO4)3", {"Li": 1, "V": 3, "O": -2, "P": 5})
        self.mgnif6 = self.load_from_cif("MgNiF6.cif", {"Mg": 2, "Ni": 4, "F": -1}, working_ion="F")

    def test_oxid_check(self):
        s = self.get_structure("LiFePO4")
        with pytest.raises(ValueError):
            BatteryAnalyzer(s, "Li")

    def test_capacitygrav_calculations(self):
        lifepo4_cap = 169.89053  # same as fepo4 cap
        nafepo4_cap = 154.20331
        la2coo4f_cap = 175.6564
        li3v2p3o12_cap_remove = 197.25339
        li3v2p3o12_cap_insert = 127.17129

        assert self.lifepo4.get_max_capgrav() == approx(lifepo4_cap, abs=1e-3)
        assert self.lifepo4.get_max_capgrav(remove=False) == 0
        assert self.lifepo4.get_max_capgrav(insert=False) == approx(lifepo4_cap, abs=1e-3)

        assert self.nafepo4.get_max_capgrav() == approx(nafepo4_cap, abs=1e-3)
        assert self.nafepo4.get_max_capgrav(remove=False) == 0

        assert self.fepo4.get_max_capgrav() == approx(lifepo4_cap, abs=1e-3)
        assert self.fepo4.get_max_capgrav(insert=False) == 0

        assert self.la2coo4f.get_max_capgrav() == approx(la2coo4f_cap, abs=1e-3)
        assert self.la2coo4.get_max_capgrav() == approx(la2coo4f_cap, abs=1e-3)
        assert self.la2coo4.get_max_capgrav(insert=False) == 0

        assert self.li3v2p3o12.get_max_capgrav(insert=False) == approx(li3v2p3o12_cap_remove, abs=1e-3)
        assert self.li3v2p3o12.get_max_capgrav(remove=False) == approx(li3v2p3o12_cap_insert, abs=1e-3)

    def test_capacityvol_calculations(self):
        lifepo4_cap = 594.17518
        nafepo4_cap = 542.86104

        fepo4_cap = 624.82289  # this is different than lifepo4 cap if lifepo4 volume not known

        assert self.lifepo4.get_max_capvol() == approx(lifepo4_cap, abs=1e-3)
        assert self.lifepo4.get_max_capvol(remove=False) == 0
        assert self.lifepo4.get_max_capvol(insert=False) == approx(lifepo4_cap, abs=1e-3)

        assert self.nafepo4.get_max_capvol() == approx(nafepo4_cap, abs=1e-3)
        assert self.nafepo4.get_max_capvol(remove=False) == 0
        assert self.nafepo4.get_max_capvol(insert=False) == approx(nafepo4_cap, abs=1e-3)

        assert self.fepo4.get_max_capvol() == approx(fepo4_cap, abs=1e-3)
        assert self.fepo4.get_max_capvol(remove=False) == approx(fepo4_cap, abs=1e-3)
        assert self.fepo4.get_max_capvol(insert=False) == 0

        # give the lifepo4 volume, should get lifepo4 capacity
        assert self.fepo4.get_max_capvol(volume=self.lifepo4.struc_oxid.volume) == approx(lifepo4_cap, abs=1e-3)

    def test_ion_removal(self):
        assert self.lifemnpo4.get_removals_int_oxid() == {1.0, 2.0, 3.0, 4.0}

        assert self.li8nicofe208.get_removals_int_oxid() == {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0}

        assert self.li3v2p3o12.get_removals_int_oxid() == {4.0, 6.0}

        assert self.mgnif6.get_removals_int_oxid() == {1.0, 2.0}


if __name__ == "__main__":
    unittest.main()
