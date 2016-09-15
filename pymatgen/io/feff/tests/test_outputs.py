# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import unittest2 as unittest

import os

from pymatgen.io.feff.outputs import LDos, Xmu


test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        'test_files')


class FeffLdosTest(unittest.TestCase):

    filepath1 = os.path.join(test_dir, 'feff.inp')
    filepath2 = os.path.join(test_dir, 'ldos')
    l = LDos.from_file(filepath1, filepath2)

    def test_init(self):
        efermi = FeffLdosTest.l.complete_dos.efermi
        self.assertEqual(efermi, -11.430,
                         "Did not read correct Fermi energy from ldos file")

    def test_complete_dos(self):
        complete_dos = FeffLdosTest.l.complete_dos
        self.assertEqual(complete_dos.as_dict()['spd_dos']["s"]['efermi'],
                         - 11.430,
                         "Failed to construct complete_dos dict properly")

    def test_as_dict_and_from_dict(self):
        l2 = FeffLdosTest.l.charge_transfer_to_string()
        d = FeffLdosTest.l.as_dict()
        l3 = LDos.from_dict(d).charge_transfer_to_string()
        self.assertEqual(l2, l3, "Feffldos to and from dict does not match")


class XmuTest(unittest.TestCase):

    def test_init(self):
        filepath1 = os.path.join(test_dir, 'xmu.dat')
        filepath2 = os.path.join(test_dir, 'feff.inp')
        x = Xmu.from_file(filepath1, filepath2)
        self.assertEqual(x.absorbing_atom, 'O',
                         "failed to read xmu.dat file properly")

    def test_as_dict_and_from_dict(self):
        filepath1 = os.path.join(test_dir, 'xmu.dat')
        filepath2 = os.path.join(test_dir, 'feff.inp')
        x = Xmu.from_file(filepath1, filepath2)
        data=x.data.tolist()
        d=x.as_dict()
        x2 = Xmu.from_dict(d)
        data2= x2.data.tolist()
        self.assertEqual(data, data2, "Xmu to and from dict does not match")


if __name__ == '__main__':
    unittest.main()
