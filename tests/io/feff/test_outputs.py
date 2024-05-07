from __future__ import annotations

from unittest import TestCase

from pymatgen.io.feff.outputs import LDos, Xmu
from pymatgen.util.testing import TEST_FILES_DIR

FEFF_TEST_DIR = f"{TEST_FILES_DIR}/io/feff"


class TestFeffLdos(TestCase):
    filepath1 = f"{FEFF_TEST_DIR}/feff.inp"
    filepath2 = f"{FEFF_TEST_DIR}/ldos"
    ldos = LDos.from_file(filepath1, filepath2)

    reci_feffinp = f"{FEFF_TEST_DIR}/feff_reci_dos/feff.inp"
    reci_ldos = f"{FEFF_TEST_DIR}/feff_reci_dos/ldos"
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
        l2 = TestFeffLdos.ldos.charge_transfer_to_str()
        dct = TestFeffLdos.ldos.as_dict()
        l3 = LDos.from_dict(dct).charge_transfer_to_str()
        assert l2 == l3, "Feffldos to and from dict does not match"

    def test_reci_init(self):
        e_fermi = TestFeffLdos.reci_dos.complete_dos.efermi
        assert e_fermi == -9.672, "Did not read correct Fermi energy from ldos file"

    def test_reci_complete_dos(self):
        complete_dos = TestFeffLdos.reci_dos.complete_dos
        assert (
            complete_dos.as_dict()["spd_dos"]["s"]["efermi"] == -9.672
        ), "Failed to construct complete_dos dict properly"

    def test_reci_charge(self):
        charge_trans = TestFeffLdos.reci_dos.charge_transfer
        assert charge_trans["0"]["Na"]["s"] == 0.241
        assert charge_trans["1"]["O"]["tot"] == -0.594


class TestXmu(TestCase):
    def test_init(self):
        filepath1 = f"{FEFF_TEST_DIR}/xmu.dat"
        filepath2 = f"{FEFF_TEST_DIR}/feff.inp"
        x = Xmu.from_file(filepath1, filepath2)
        assert x.absorbing_atom == "O", "failed to read xmu.dat file properly"

    def test_as_dict_and_from_dict(self):
        filepath1 = f"{FEFF_TEST_DIR}/xmu.dat"
        filepath2 = f"{FEFF_TEST_DIR}/feff.inp"
        x = Xmu.from_file(filepath1, filepath2)
        data = x.data.tolist()
        dct = x.as_dict()
        x2 = Xmu.from_dict(dct)
        data2 = x2.data.tolist()
        assert data == data2, "Xmu to and from dict does not match"
