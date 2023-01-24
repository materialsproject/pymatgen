# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


from __future__ import annotations

import unittest

import pytest

from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.symmetry.kpath import KPathSeek
from pymatgen.util.testing import PymatgenTest

try:
    from seekpath import get_path
except ImportError:
    get_path = None


class KPathSeekTest(PymatgenTest):
    @unittest.skipIf(get_path is None, "No seek path present.")
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
                        [3.0233057319441246, 1, 0],
                        [0, 7.9850357844548681, 1],
                        [0, 1.2, 8.1136762279561818],
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
            kpath = KPathSeek(struct)  # Throws error if something doesn't work, causing test to fail.
            kpoints = kpath.get_kpoints()  # noqa: F841

    @unittest.skipIf(get_path is None, "No seek path present.")
    def test_kpath_acentered(self):
        species = ["K", "La", "Ti"]
        coords = [[0.345, 5, 0.77298], [0.1345, 5.1, 0.77298], [0.7, 0.8, 0.9]]
        lattice = Lattice.orthorhombic(2, 9, 1)
        struct = Structure.from_spacegroup(38, lattice, species, coords)
        kpath = KPathSeek(struct)

        kpoints = kpath._kpath["kpoints"]
        labels = list(kpoints)

        assert sorted(labels) == sorted(
            ["B_0", "B_2", "DELTA_0", "F_0", "GAMMA", "G_0", "G_2", "R", "R_2", "S", "T", "T_2", "Y", "Z", "Z_2"]
        )

        assert kpoints["GAMMA"][0] == pytest.approx(0.0)
        assert kpoints["GAMMA"][1] == pytest.approx(0.0)
        assert kpoints["GAMMA"][2] == pytest.approx(0.0)

        assert kpoints["Y"][0] == pytest.approx(0.5)
        assert kpoints["Y"][1] == pytest.approx(0.5)
        assert kpoints["Y"][2] == pytest.approx(0.0)

        assert kpoints["T"][0] == pytest.approx(0.5)
        assert kpoints["T"][1] == pytest.approx(0.5)
        assert kpoints["T"][2] == pytest.approx(0.5)

        assert kpoints["T_2"][0] == pytest.approx(0.5)
        assert kpoints["T_2"][1] == pytest.approx(0.5)
        assert kpoints["T_2"][2] == pytest.approx(-0.5)

        assert kpoints["Z"][0] == pytest.approx(0.0)
        assert kpoints["Z"][1] == pytest.approx(0.0)
        assert kpoints["Z"][2] == pytest.approx(0.5)

        assert kpoints["Z_2"][0] == pytest.approx(0.0)
        assert kpoints["Z_2"][1] == pytest.approx(0.0)
        assert kpoints["Z_2"][2] == pytest.approx(-0.5)

        assert kpoints["S"][0] == pytest.approx(0.0)
        assert kpoints["S"][1] == pytest.approx(0.5)
        assert kpoints["S"][2] == pytest.approx(0.0)

        assert kpoints["R"][0] == pytest.approx(0.0)
        assert kpoints["R"][1] == pytest.approx(0.5)
        assert kpoints["R"][2] == pytest.approx(0.5)

        assert kpoints["R_2"][0] == pytest.approx(0.0)
        assert kpoints["R_2"][1] == pytest.approx(0.5)
        assert kpoints["R_2"][2] == pytest.approx(-0.5)

        assert kpoints["DELTA_0"][0] == pytest.approx(-0.25308641975308643)
        assert kpoints["DELTA_0"][1] == pytest.approx(0.25308641975308643)
        assert kpoints["DELTA_0"][2] == pytest.approx(0.0)

        assert kpoints["F_0"][0] == pytest.approx(0.25308641975308643)
        assert kpoints["F_0"][1] == pytest.approx(0.7469135802469136)
        assert kpoints["F_0"][2] == pytest.approx(0.0)

        assert kpoints["B_0"][0] == pytest.approx(-0.25308641975308643)
        assert kpoints["B_0"][1] == pytest.approx(0.25308641975308643)
        assert kpoints["B_0"][2] == pytest.approx(0.5)

        assert kpoints["B_2"][0] == pytest.approx(-0.25308641975308643)
        assert kpoints["B_2"][1] == pytest.approx(0.25308641975308643)
        assert kpoints["B_2"][2] == pytest.approx(-0.5)

        assert kpoints["G_0"][0] == pytest.approx(0.25308641975308643)
        assert kpoints["G_0"][1] == pytest.approx(0.7469135802469136)
        assert kpoints["G_0"][2] == pytest.approx(0.5)

        assert kpoints["G_2"][0] == pytest.approx(0.25308641975308643)
        assert kpoints["G_2"][1] == pytest.approx(0.7469135802469136)
        assert kpoints["G_2"][2] == pytest.approx(-0.5)


if __name__ == "__main__":
    unittest.main()
