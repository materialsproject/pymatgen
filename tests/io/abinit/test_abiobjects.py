from __future__ import annotations

import numpy as np
import pytest
from numpy.testing import assert_array_equal
from pytest import approx

from pymatgen.core.structure import Structure
from pymatgen.core.units import Ha_to_eV, bohr_to_ang
from pymatgen.io.abinit.abiobjects import (
    Electrons,
    ElectronsAlgorithm,
    KSampling,
    PPModel,
    RelaxationMethod,
    Smearing,
    SpinMode,
    lattice_from_abivars,
    species_by_znucl,
    structure_to_abivars,
)
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest


class TestLatticeFromAbivars(PymatgenTest):
    def test_rprim_acell(self):
        l1 = lattice_from_abivars(acell=3 * [10], rprim=np.eye(3))
        assert l1.volume == approx(bohr_to_ang**3 * 1000)
        assert l1.angles == (90, 90, 90)
        l2 = lattice_from_abivars(acell=3 * [10], angdeg=(90, 90, 90))
        assert l1 == l2

        l2 = lattice_from_abivars(acell=3 * [8], angdeg=(60, 60, 60))
        abi_rprimd = (
            np.reshape(
                [
                    4.6188022,
                    0.0000000,
                    6.5319726,
                    -2.3094011,
                    4.0000000,
                    6.5319726,
                    -2.3094011,
                    -4.0000000,
                    6.5319726,
                ],
                (3, 3),
            )
            * bohr_to_ang
        )
        assert np.allclose(l2.matrix, abi_rprimd)

        l3 = lattice_from_abivars(acell=[3, 6, 9], angdeg=(30, 40, 50))
        abi_rprimd = (
            np.reshape(
                [
                    3.0000000,
                    0.0000000,
                    0.0000000,
                    3.8567257,
                    4.5962667,
                    0.0000000,
                    6.8944000,
                    4.3895544,
                    3.7681642,
                ],
                (3, 3),
            )
            * bohr_to_ang
        )
        assert np.allclose(l3.matrix, abi_rprimd)

        with pytest.raises(ValueError, match="angdeg and rprimd are mutually exclusive"):
            lattice_from_abivars(acell=[1, 1, 1], angdeg=(90, 90, 90), rprim=np.eye(3))
        with pytest.raises(ValueError, match=r"Angles must be > 0 but got \[-90  90  90\]"):
            lattice_from_abivars(acell=[1, 1, 1], angdeg=(-90, 90, 90))

    def test_znucl_typat(self):
        """Test the order of typat and znucl in the Abinit input and enforce_typat, enforce_znucl."""
        # Ga  Ga1  1  0.33333333333333  0.666666666666667  0.500880  1.0
        # Ga  Ga2  1  0.66666666666667  0.333333333333333  0.000880  1.0
        # N  N3  1  0.333333333333333  0.666666666666667  0.124120  1.0
        # N  N4  1  0.666666666666667  0.333333333333333  0.624120  1.0
        gan = Structure.from_file(f"{TEST_FILES_DIR}/abinit/gan.cif")

        # By default, znucl is filled using the first new type found in sites.
        def_vars = structure_to_abivars(gan)
        def_znucl = def_vars["znucl"]
        assert_array_equal(def_znucl, [31, 7])
        def_typat = def_vars["typat"]
        assert_array_equal(def_typat, [1, 1, 2, 2])

        # But it's possible to enforce a particular value of typat and znucl.
        enforce_znucl = [7, 31]
        enforce_typat = [2, 2, 1, 1]
        enf_vars = structure_to_abivars(gan, enforce_znucl=enforce_znucl, enforce_typat=enforce_typat)
        assert_array_equal(enf_vars["znucl"], enforce_znucl)
        assert_array_equal(enf_vars["typat"], enforce_typat)
        assert_array_equal(def_vars["xred"], enf_vars["xred"])

        assert [s.symbol for s in species_by_znucl(gan)] == ["Ga", "N"]

        for itype1, itype2 in zip(def_typat, enforce_typat):
            assert def_znucl[itype1 - 1] == enforce_znucl[itype2 - 1]

        with pytest.raises(ValueError, match="Both enforce_znucl and enforce_typat are required"):
            structure_to_abivars(gan, enforce_znucl=enforce_znucl, enforce_typat=None)


class TestSpinMode(PymatgenTest):
    def test_base(self):
        polarized = SpinMode.as_spinmode("polarized")
        other_polarized = SpinMode.as_spinmode("polarized")
        unpolarized = SpinMode.as_spinmode("unpolarized")

        polarized.to_abivars()

        assert polarized is other_polarized
        assert polarized == other_polarized
        assert polarized != unpolarized

        # Test pickle
        self.serialize_with_pickle(polarized)

        # Test dict methods
        self.assert_msonable(polarized)
        self.assert_msonable(unpolarized)


class TestSmearing(PymatgenTest):
    def test_base(self):
        fd1ev = Smearing.as_smearing("fermi_dirac:1 eV")
        fd1ev.to_abivars()

        assert fd1ev

        same_fd = Smearing.as_smearing(f"fermi_dirac:{1.0 / Ha_to_eV}")

        assert same_fd == fd1ev

        no_smear = Smearing.nosmearing()
        assert no_smear == Smearing.as_smearing("nosmearing")

        assert not no_smear
        assert no_smear != fd1ev
        self.assert_msonable(no_smear)

        new_fd1ev = Smearing.from_dict(fd1ev.as_dict())
        assert new_fd1ev == fd1ev

        # Test pickle
        self.serialize_with_pickle(fd1ev)

        # Test dict methods
        self.assert_msonable(fd1ev)


class TestElectronsAlgorithm(PymatgenTest):
    def test_base(self):
        algo = ElectronsAlgorithm(nstep=70)
        _ = algo.to_abivars()

        # Test pickle
        self.serialize_with_pickle(algo)

        # Test dict methods
        self.assert_msonable(algo)


class TestElectrons(PymatgenTest):
    def test_base(self):
        default_electrons = Electrons()
        assert default_electrons.nsppol == 2
        assert default_electrons.nspinor == 1
        assert default_electrons.nspden == 2

        _ = default_electrons.to_abivars()

        # new = Electron.from_dict(default_electrons.as_dict())

        # Test pickle
        self.serialize_with_pickle(default_electrons, test_eq=False)

        custom_electrons = Electrons(
            spin_mode="unpolarized",
            smearing="marzari4:0.2 eV",
            algorithm=ElectronsAlgorithm(nstep=70),
            nband=10,
            charge=1.0,
            comment="Test comment",
        )

        # Test dict methods
        self.assert_msonable(custom_electrons)


class TestKSampling(PymatgenTest):
    def test_base(self):
        monkhorst = KSampling.monkhorst((3, 3, 3), (0.5, 0.5, 0.5), 0, False, False)
        gamma_centered = KSampling.gamma_centered((3, 3, 3), False, False)

        monkhorst.to_abivars()

        # Test dict methods
        self.assert_msonable(monkhorst)
        self.assert_msonable(gamma_centered)


class TestRelaxation(PymatgenTest):
    def test_base(self):
        atoms_and_cell = RelaxationMethod.atoms_and_cell()
        atoms_only = RelaxationMethod.atoms_only()

        out_vars = atoms_and_cell.to_abivars()
        assert {*out_vars} >= {"dilatmx", "ecutsm", "ionmov", "ntime", "optcell", "strfact"}

        # Test dict methods
        self.assert_msonable(atoms_and_cell)
        self.assert_msonable(atoms_only)


class TestPPModel(PymatgenTest):
    def test_base(self):
        godby = PPModel.as_ppmodel("godby:12 eV")
        godby.to_abivars()
        assert godby

        same_godby = PPModel.as_ppmodel(f"godby:{12.0 / Ha_to_eV}")
        assert same_godby == godby

        no_pp_model = PPModel.get_noppmodel()

        assert not no_pp_model
        assert no_pp_model != godby
        new_godby = PPModel.from_dict(godby.as_dict())
        assert new_godby == godby

        # Test pickle
        self.serialize_with_pickle(godby)

        # Test dict methods
        self.assert_msonable(godby)
