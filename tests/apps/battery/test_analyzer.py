from __future__ import annotations

import os

import pytest
from pytest import approx

from pymatgen.apps.battery.analyzer import BatteryAnalyzer
from pymatgen.core.structure import Structure
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest


class TestBatteryAnalyzer(PymatgenTest):
    def load_from_cif(self, filename, oxidations, working_ion="Li"):
        struct = Structure.from_file(os.path.join(TEST_FILES_DIR, filename))
        struct.add_oxidation_state_by_element(oxidations)
        return BatteryAnalyzer(struct, working_ion)

    def load_from_internal(self, name, oxidations, working_ion="Li"):
        struct = self.get_structure(name).copy()
        struct.add_oxidation_state_by_element(oxidations)
        return BatteryAnalyzer(struct, working_ion)

    def setUp(self):
        self.li_fe_p_o4 = self.load_from_internal("LiFePO4", {"Li": 1, "Fe": 2, "P": 5, "O": -2})
        self.na_fe_p_o4 = self.load_from_internal("NaFePO4", {"Na": 1, "Fe": 2, "P": 5, "O": -2}, working_ion="Na")
        self.la2coo4f = self.load_from_internal("La2CoO4F", {"La": 3, "Co": 3, "O": -2, "F": -1}, working_ion="F")
        self.fe_p_o4 = self.load_from_cif("FePO4a.cif", {"Fe": 3, "P": 5, "O": -2})
        self.la2coo4 = self.load_from_cif("La2CoO4.cif", {"La": 3, "Co": 2, "O": -2}, working_ion="F")
        self.lifemnpo4 = self.load_from_cif("Li4Fe3Mn1(PO4)4.cif", {"Li": 1, "Fe": 2, "Mn": 2, "P": 5, "O": -2})
        self.li8nicofe208 = self.load_from_cif(
            "Li8Fe2NiCoO8.cif", {"Li": 1, "Fe": 2, "Mn": 2, "Co": 2, "Ni": 2, "O": -2}
        )
        self.li3v2p3o12 = self.load_from_internal("Li3V2(PO4)3", {"Li": 1, "V": 3, "O": -2, "P": 5})
        self.mgnif6 = self.load_from_cif("MgNiF6.cif", {"Mg": 2, "Ni": 4, "F": -1}, working_ion="F")

    def test_oxide_check(self):
        struct = self.get_structure("LiFePO4")
        with pytest.raises(ValueError, match="BatteryAnalyzer requires oxidation states assigned to structure"):
            BatteryAnalyzer(struct, "Li")

    def test_capacitygrav_calculations(self):
        li_fe_p_o4_cap = 169.89053  # same as fe_po4 cap
        na_fe_p_o4_cap = 154.20331
        la2_co_o4_f_cap = 175.6564
        li3_v2_p3_o12_cap_remove = 197.25339
        li3_v2_p3_o12_cap_insert = 127.17129

        assert self.li_fe_p_o4.get_max_capgrav() == approx(li_fe_p_o4_cap, abs=1e-3)
        assert self.li_fe_p_o4.get_max_capgrav(remove=False) == 0
        assert self.li_fe_p_o4.get_max_capgrav(insert=False) == approx(li_fe_p_o4_cap, abs=1e-3)

        assert self.na_fe_p_o4.get_max_capgrav() == approx(na_fe_p_o4_cap, abs=1e-3)
        assert self.na_fe_p_o4.get_max_capgrav(remove=False) == 0

        assert self.fe_p_o4.get_max_capgrav() == approx(li_fe_p_o4_cap, abs=1e-3)
        assert self.fe_p_o4.get_max_capgrav(insert=False) == 0

        assert self.la2coo4f.get_max_capgrav() == approx(la2_co_o4_f_cap, abs=1e-3)
        assert self.la2coo4.get_max_capgrav() == approx(la2_co_o4_f_cap, abs=1e-3)
        assert self.la2coo4.get_max_capgrav(insert=False) == 0

        assert self.li3v2p3o12.get_max_capgrav(insert=False) == approx(li3_v2_p3_o12_cap_remove, abs=1e-3)
        assert self.li3v2p3o12.get_max_capgrav(remove=False) == approx(li3_v2_p3_o12_cap_insert, abs=1e-3)

    def test_capacityvol_calculations(self):
        li_fe_p_o4_cap = 594.17518
        na_fe_p_o4_cap = 542.86104

        fe_p_o4_cap = 624.82289  # this is different than li_fe_p_o4_cap cap if li_fe_p_o4 volume not known

        assert self.li_fe_p_o4.get_max_capvol() == approx(li_fe_p_o4_cap, abs=1e-3)
        assert self.li_fe_p_o4.get_max_capvol(remove=False) == 0
        assert self.li_fe_p_o4.get_max_capvol(insert=False) == approx(li_fe_p_o4_cap, abs=1e-3)

        assert self.na_fe_p_o4.get_max_capvol() == approx(na_fe_p_o4_cap, abs=1e-3)
        assert self.na_fe_p_o4.get_max_capvol(remove=False) == 0
        assert self.na_fe_p_o4.get_max_capvol(insert=False) == approx(na_fe_p_o4_cap, abs=1e-3)

        assert self.fe_p_o4.get_max_capvol() == approx(fe_p_o4_cap, abs=1e-3)
        assert self.fe_p_o4.get_max_capvol(remove=False) == approx(fe_p_o4_cap, abs=1e-3)
        assert self.fe_p_o4.get_max_capvol(insert=False) == 0

        # give the lifepo4 volume, should get lifepo4 capacity
        assert self.fe_p_o4.get_max_capvol(volume=self.li_fe_p_o4.struc_oxid.volume) == approx(li_fe_p_o4_cap, abs=1e-3)

    def test_ion_removal(self):
        assert self.lifemnpo4.get_removals_int_oxid() == {1, 2, 3, 4}

        assert self.li8nicofe208.get_removals_int_oxid() == {1, 2, 3, 4, 5, 6, 7, 8}

        assert self.li3v2p3o12.get_removals_int_oxid() == {4, 6}

        assert self.mgnif6.get_removals_int_oxid() == {1, 2}
