# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import os

from pymatgen.util.testing import PymatgenTest
from pymatgen.core.structure import Structure
from pymatgen.core.units import Ha_to_eV, bohr_to_ang
from pymatgen.io.abinit.abiobjects import *

import warnings

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        'test_files')


class LatticeFromAbivarsTest(PymatgenTest):
    def test_rprim_acell(self):
        l1 = lattice_from_abivars(acell=3 * [10], rprim=np.eye(3))
        self.assertAlmostEqual(l1.volume, bohr_to_ang ** 3 * 1000)
        assert l1.angles == (90, 90, 90)
        l2 = lattice_from_abivars(acell=3 * [10], angdeg=(90, 90, 90))
        assert l1 == l2

        l2 = lattice_from_abivars(acell=3 * [8], angdeg=(60, 60, 60))
        abi_rprimd = np.reshape([4.6188022, 0.0000000, 6.5319726,
                                 -2.3094011, 4.0000000, 6.5319726,
                                 -2.3094011, -4.0000000, 6.5319726], (3, 3)) * bohr_to_ang
        self.assertArrayAlmostEqual(l2.matrix, abi_rprimd)

        l3 = lattice_from_abivars(acell=[3, 6, 9], angdeg=(30, 40, 50))
        abi_rprimd = np.reshape([3.0000000, 0.0000000, 0.0000000,
                                 3.8567257, 4.5962667, 0.0000000,
                                 6.8944000, 4.3895544, 3.7681642], (3, 3)) * bohr_to_ang
        self.assertArrayAlmostEqual(l3.matrix, abi_rprimd)

        with self.assertRaises(ValueError):
            lattice_from_abivars(acell=[1, 1, 1], angdeg=(90, 90, 90), rprim=np.eye(3))
        with self.assertRaises(ValueError):
            lattice_from_abivars(acell=[1, 1, 1], angdeg=(-90, 90, 90))


class SpinModeTest(PymatgenTest):

    def test_base(self):
        polarized = SpinMode.as_spinmode("polarized")
        other_polarized = SpinMode.as_spinmode("polarized")
        unpolarized = SpinMode.as_spinmode("unpolarized")

        polarized.to_abivars()

        self.assertTrue(polarized is other_polarized)
        self.assertTrue(polarized == other_polarized)
        self.assertTrue(polarized != unpolarized)

        # Test pickle
        self.serialize_with_pickle(polarized)

        # Test dict methods
        self.assertMSONable(polarized)
        self.assertMSONable(unpolarized)


class SmearingTest(PymatgenTest):
    def test_base(self):
        fd1ev = Smearing.as_smearing("fermi_dirac:1 eV")
        fd1ev.to_abivars()

        self.assertTrue(fd1ev)

        same_fd = Smearing.as_smearing("fermi_dirac:" + str(1.0 / Ha_to_eV))

        self.assertTrue(same_fd == fd1ev)

        nosmear = Smearing.nosmearing()
        assert nosmear == Smearing.as_smearing("nosmearing")

        self.assertFalse(nosmear)
        self.assertTrue(nosmear != fd1ev)
        self.assertMSONable(nosmear)

        new_fd1ev = Smearing.from_dict(fd1ev.as_dict())
        self.assertTrue(new_fd1ev == fd1ev)

        # Test pickle
        self.serialize_with_pickle(fd1ev)

        # Test dict methods
        self.assertMSONable(fd1ev)


class ElectronsAlgorithmTest(PymatgenTest):
    def test_base(self):
        algo = ElectronsAlgorithm(nstep=70)
        abivars = algo.to_abivars()

        # Test pickle
        self.serialize_with_pickle(algo)

        # Test dict methods
        self.assertMSONable(algo)


class ElectronsTest(PymatgenTest):
    def test_base(self):
        default_electrons = Electrons()
        self.assertTrue(default_electrons.nsppol == 2)
        self.assertTrue(default_electrons.nspinor == 1)
        self.assertTrue(default_electrons.nspden == 2)

        abivars = default_electrons.to_abivars()

        # new = Electron.from_dict(default_electrons.as_dict())

        # Test pickle
        self.serialize_with_pickle(default_electrons, test_eq=False)

        custom_electrons = Electrons(spin_mode="unpolarized", smearing="marzari4:0.2 eV",
                                     algorithm=ElectronsAlgorithm(nstep=70), nband=10, charge=1.0,
                                     comment="Test comment")

        # Test dict methods
        self.assertMSONable(custom_electrons)


class KSamplingTest(PymatgenTest):

    def test_base(self):
        monkhorst = KSampling.monkhorst((3, 3, 3), (0.5, 0.5, 0.5), 0, False, False)
        gamma_centered = KSampling.gamma_centered((3, 3, 3), False, False)

        monkhorst.to_abivars()

        # Test dict methods
        self.assertMSONable(monkhorst)
        self.assertMSONable(gamma_centered)


class RelaxationTest(PymatgenTest):

    def test_base(self):
        atoms_and_cell = RelaxationMethod.atoms_and_cell()
        atoms_only = RelaxationMethod.atoms_only()

        atoms_and_cell.to_abivars()

        # Test dict methods
        self.assertMSONable(atoms_and_cell)
        self.assertMSONable(atoms_only)


class PPModelTest(PymatgenTest):

    def test_base(self):
        godby = PPModel.as_ppmodel("godby:12 eV")
        # print(godby)
        # print(repr(godby))
        godby.to_abivars()
        self.assertTrue(godby)

        same_godby = PPModel.as_ppmodel("godby:" + str(12.0 / Ha_to_eV))
        self.assertTrue(same_godby == godby)

        noppm = PPModel.get_noppmodel()

        self.assertFalse(noppm)
        self.assertTrue(noppm != godby)
        new_godby = PPModel.from_dict(godby.as_dict())
        self.assertTrue(new_godby == godby)

        # Test pickle
        self.serialize_with_pickle(godby)

        # Test dict methods
        self.assertMSONable(godby)


if __name__ == '__main__':
    import unittest

    unittest.main()
