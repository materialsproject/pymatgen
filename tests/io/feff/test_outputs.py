from __future__ import annotations

import unittest

from pymatgen.io.feff.outputs import LDos, Xmu
from pymatgen.util.testing import TEST_FILES_DIR

test_dir_reci = f"{TEST_FILES_DIR}/feff_reci_dos"


class TestFeffLdos(unittest.TestCase):
    filepath1 = f"{TEST_FILES_DIR}/feff.inp"
    filepath2 = f"{TEST_FILES_DIR}/ldos"
    ldos = LDos.from_file(filepath1, filepath2)

    reci_feffinp = f"{test_dir_reci}/feff.inp"
    reci_ldos = f"{test_dir_reci}/ldos"
    reci_dos = LDos.from_file(reci_feffinp, reci_ldos)

    def test_init(self):
        e_fermi = TestFeffLdos.ldos.complete_dos.efermi
        assert e_fermi == -11.430, "Did not read correct Fermi energy from ldos file"

    def test_complete_dos(self):
        complete_dos = TestFeffLdos.ldos.complete_dos
        assert (
            complete_dos.as_dict()["spd_dos"]["s"]["efermi"] == -11.430
        ), "Failed to construct complete_dos dict properly"

    def test_as_dict_and_from_dict(self):
        l2 = TestFeffLdos.ldos.charge_transfer_to_string()
        d = TestFeffLdos.ldos.as_dict()
        l3 = LDos.from_dict(d).charge_transfer_to_string()
        assert l2 == l3, "Feffldos to and from dict does not match"

    def test_reci_init(self):
        efermi = TestFeffLdos.reci_dos.complete_dos.efermi
        assert efermi == -9.672, "Did not read correct Fermi energy from ldos file"

    def test_reci_complete_dos(self):
        complete_dos = TestFeffLdos.reci_dos.complete_dos
        assert (
            complete_dos.as_dict()["spd_dos"]["s"]["efermi"] == -9.672
        ), "Failed to construct complete_dos dict properly"

    def test_reci_charge(self):
        charge_trans = TestFeffLdos.reci_dos.charge_transfer
        assert charge_trans["0"]["Na"]["s"] == 0.241
        assert charge_trans["1"]["O"]["tot"] == -0.594


class TestXmu(unittest.TestCase):
    def test_init(self):
        filepath1 = f"{TEST_FILES_DIR}/xmu.dat"
        filepath2 = f"{TEST_FILES_DIR}/feff.inp"
        x = Xmu.from_file(filepath1, filepath2)
        assert x.absorbing_atom == "O", "failed to read xmu.dat file properly"

    def test_as_dict_and_from_dict(self):
        filepath1 = f"{TEST_FILES_DIR}/xmu.dat"
        filepath2 = f"{TEST_FILES_DIR}/feff.inp"
        x = Xmu.from_file(filepath1, filepath2)
        data = x.data.tolist()
        d = x.as_dict()
        x2 = Xmu.from_dict(d)
        data2 = x2.data.tolist()
        assert data == data2, "Xmu to and from dict does not match"
