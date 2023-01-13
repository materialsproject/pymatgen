# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


from __future__ import annotations

import os
import unittest

from pymatgen.io.feff.outputs import LDos, Xmu
from pymatgen.util.testing import PymatgenTest

test_dir_reci = os.path.join(PymatgenTest.TEST_FILES_DIR, "feff_reci_dos")


class FeffLdosTest(unittest.TestCase):
    filepath1 = os.path.join(PymatgenTest.TEST_FILES_DIR, "feff.inp")
    filepath2 = os.path.join(PymatgenTest.TEST_FILES_DIR, "ldos")
    l = LDos.from_file(filepath1, filepath2)

    reci_feffinp = os.path.join(test_dir_reci, "feff.inp")
    reci_ldos = os.path.join(test_dir_reci, "ldos")
    reci_dos = LDos.from_file(reci_feffinp, reci_ldos)

    def test_init(self):
        efermi = FeffLdosTest.l.complete_dos.efermi
        self.assertEqual(efermi, -11.430, "Did not read correct Fermi energy from ldos file")

    def test_complete_dos(self):
        complete_dos = FeffLdosTest.l.complete_dos
        self.assertEqual(
            complete_dos.as_dict()["spd_dos"]["s"]["efermi"],
            -11.430,
            "Failed to construct complete_dos dict properly",
        )

    def test_as_dict_and_from_dict(self):
        l2 = FeffLdosTest.l.charge_transfer_to_string()
        d = FeffLdosTest.l.as_dict()
        l3 = LDos.from_dict(d).charge_transfer_to_string()
        self.assertEqual(l2, l3, "Feffldos to and from dict does not match")

    def test_reci_init(self):
        efermi = FeffLdosTest.reci_dos.complete_dos.efermi
        self.assertEqual(efermi, -9.672, "Did not read correct Fermi energy from ldos file")

    def test_reci_complete_dos(self):
        complete_dos = FeffLdosTest.reci_dos.complete_dos
        self.assertEqual(
            complete_dos.as_dict()["spd_dos"]["s"]["efermi"],
            -9.672,
            "Failed to construct complete_dos dict properly",
        )

    def test_reci_charge(self):
        charge_trans = FeffLdosTest.reci_dos.charge_transfer
        self.assertEqual(charge_trans["0"]["Na"]["s"], 0.241)
        self.assertEqual(charge_trans["1"]["O"]["tot"], -0.594)


class XmuTest(unittest.TestCase):
    def test_init(self):
        filepath1 = os.path.join(PymatgenTest.TEST_FILES_DIR, "xmu.dat")
        filepath2 = os.path.join(PymatgenTest.TEST_FILES_DIR, "feff.inp")
        x = Xmu.from_file(filepath1, filepath2)
        self.assertEqual(x.absorbing_atom, "O", "failed to read xmu.dat file properly")

    def test_as_dict_and_from_dict(self):
        filepath1 = os.path.join(PymatgenTest.TEST_FILES_DIR, "xmu.dat")
        filepath2 = os.path.join(PymatgenTest.TEST_FILES_DIR, "feff.inp")
        x = Xmu.from_file(filepath1, filepath2)
        data = x.data.tolist()
        d = x.as_dict()
        x2 = Xmu.from_dict(d)
        data2 = x2.data.tolist()
        self.assertEqual(data, data2, "Xmu to and from dict does not match")


if __name__ == "__main__":
    unittest.main()
