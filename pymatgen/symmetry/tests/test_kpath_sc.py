# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


from __future__ import annotations

import os
import unittest

import pytest

from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.symmetry.kpath import KPathSetyawanCurtarolo
from pymatgen.util.testing import PymatgenTest

test_dir_structs = os.path.join(PymatgenTest.TEST_FILES_DIR, "space_group_structs")


class BandStructureSCTest(PymatgenTest):
    def test_kpath_generation(self):
        triclinic = [1, 2]
        monoclinic = range(3, 16)
        orthorhombic = range(16, 75)
        tetragonal = range(75, 143)
        rhombohedral = range(143, 168)
        hexagonal = range(168, 195)
        cubic = range(195, 231)

        species = ["K", "La", "Ti"]
        coords = [[0.345, 5, 0.77298], [0.1345, 5.1, 0.77298], [0.7, 0.8, 0.9]]
        for i in range(230):
            sg_num = i + 1
            if sg_num in triclinic:
                lattice = Lattice(
                    [
                        [3.0233057319441246, 0, 0],
                        [0, 7.9850357844548681, 0],
                        [0, 0, 8.1136762279561818],
                    ]
                )
            elif sg_num in monoclinic:
                lattice = Lattice.monoclinic(2, 9, 1, 99)
            elif sg_num in orthorhombic:
                lattice = Lattice.orthorhombic(2, 9, 1)
            elif sg_num in tetragonal:
                lattice = Lattice.tetragonal(2, 9)
            elif sg_num in rhombohedral:
                lattice = Lattice.hexagonal(2, 95)
            elif sg_num in hexagonal:
                lattice = Lattice.hexagonal(2, 9)
            elif sg_num in cubic:
                lattice = Lattice.cubic(2)

            struct = Structure.from_spacegroup(sg_num, lattice, species, coords)
            _ = KPathSetyawanCurtarolo(struct)  # Throws error if something doesn't work, causing test to fail.

        struct_file_path = os.path.join(test_dir_structs, "ICSD_170.cif")
        struct = Structure.from_file(struct_file_path)
        hkp = KPathSetyawanCurtarolo(struct)
        assert hkp.name == "MCLC5"

    def test_kpath_acentered(self):
        species = ["K", "La", "Ti"]
        coords = [[0.345, 5, 0.77298], [0.1345, 5.1, 0.77298], [0.7, 0.8, 0.9]]
        lattice = Lattice.orthorhombic(2, 9, 1)
        struct = Structure.from_spacegroup(38, lattice, species, coords)
        kpath = KPathSetyawanCurtarolo(struct)

        kpoints = kpath._kpath["kpoints"]
        labels = list(kpoints)

        assert sorted(labels) == sorted(["\\Gamma", "A", "A_1", "R", "S", "T", "X", "X_1", "Y", "Z"])

        assert kpoints["\\Gamma"][0] == 0.00000000
        assert kpoints["\\Gamma"][1] == pytest.approx(0.00000000)
        assert kpoints["\\Gamma"][2] == pytest.approx(0.00000000)

        assert kpoints["A"][0] == pytest.approx(0.25308642)
        assert kpoints["A"][1] == pytest.approx(0.25308642)
        assert kpoints["A"][2] == pytest.approx(0.50000000)

        assert kpoints["A_1"][0] == pytest.approx(-0.25308642)
        assert kpoints["A_1"][1] == pytest.approx(0.74691358)
        assert kpoints["A_1"][2] == pytest.approx(0.50000000)

        assert kpoints["R"][0] == pytest.approx(0.00000000)
        assert kpoints["R"][1] == pytest.approx(0.50000000)
        assert kpoints["R"][2] == pytest.approx(0.50000000)

        assert kpoints["S"][0] == pytest.approx(0.00000000)
        assert kpoints["S"][1] == pytest.approx(0.50000000)
        assert kpoints["S"][2] == pytest.approx(0.00000000)

        assert kpoints["T"][0] == pytest.approx(-0.50000000)
        assert kpoints["T"][1] == pytest.approx(0.50000000)
        assert kpoints["T"][2] == pytest.approx(0.50000000)

        assert kpoints["X"][0] == pytest.approx(0.25308642)
        assert kpoints["X"][1] == pytest.approx(0.25308642)
        assert kpoints["X"][2] == pytest.approx(0.00000000)

        assert kpoints["X_1"][0] == pytest.approx(-0.25308642)
        assert kpoints["X_1"][1] == pytest.approx(0.74691358)
        assert kpoints["X_1"][2] == pytest.approx(0.00000000)

        assert kpoints["Y"][0] == pytest.approx(-0.50000000)
        assert kpoints["Y"][1] == pytest.approx(0.50000000)
        assert kpoints["Y"][2] == pytest.approx(0.00000000)

        assert kpoints["Z"][0] == pytest.approx(0.00000000)
        assert kpoints["Z"][1] == pytest.approx(0.00000000)
        assert kpoints["Z"][2] == pytest.approx(0.50000000)


if __name__ == "__main__":
    unittest.main()
