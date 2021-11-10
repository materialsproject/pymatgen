# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


import unittest

from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.symmetry.kpath import KPathSeek
from pymatgen.util.testing import PymatgenTest

try:
    from seekpath import get_path  # type: ignore
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
        labels = list(kpoints.keys())

        self.assertEqual(
            sorted(labels),
            sorted(
                [
                    "B_0",
                    "B_2",
                    "DELTA_0",
                    "F_0",
                    "GAMMA",
                    "G_0",
                    "G_2",
                    "R",
                    "R_2",
                    "S",
                    "T",
                    "T_2",
                    "Y",
                    "Z",
                    "Z_2",
                ]
            ),
        )

        self.assertAlmostEqual(kpoints["GAMMA"][0], 0.0)
        self.assertAlmostEqual(kpoints["GAMMA"][1], 0.0)
        self.assertAlmostEqual(kpoints["GAMMA"][2], 0.0)

        self.assertAlmostEqual(kpoints["Y"][0], 0.5)
        self.assertAlmostEqual(kpoints["Y"][1], 0.5)
        self.assertAlmostEqual(kpoints["Y"][2], 0.0)

        self.assertAlmostEqual(kpoints["T"][0], 0.5)
        self.assertAlmostEqual(kpoints["T"][1], 0.5)
        self.assertAlmostEqual(kpoints["T"][2], 0.5)

        self.assertAlmostEqual(kpoints["T_2"][0], 0.5)
        self.assertAlmostEqual(kpoints["T_2"][1], 0.5)
        self.assertAlmostEqual(kpoints["T_2"][2], -0.5)

        self.assertAlmostEqual(kpoints["Z"][0], 0.0)
        self.assertAlmostEqual(kpoints["Z"][1], 0.0)
        self.assertAlmostEqual(kpoints["Z"][2], 0.5)

        self.assertAlmostEqual(kpoints["Z_2"][0], 0.0)
        self.assertAlmostEqual(kpoints["Z_2"][1], 0.0)
        self.assertAlmostEqual(kpoints["Z_2"][2], -0.5)

        self.assertAlmostEqual(kpoints["S"][0], 0.0)
        self.assertAlmostEqual(kpoints["S"][1], 0.5)
        self.assertAlmostEqual(kpoints["S"][2], 0.0)

        self.assertAlmostEqual(kpoints["R"][0], 0.0)
        self.assertAlmostEqual(kpoints["R"][1], 0.5)
        self.assertAlmostEqual(kpoints["R"][2], 0.5)

        self.assertAlmostEqual(kpoints["R_2"][0], 0.0)
        self.assertAlmostEqual(kpoints["R_2"][1], 0.5)
        self.assertAlmostEqual(kpoints["R_2"][2], -0.5)

        self.assertAlmostEqual(kpoints["DELTA_0"][0], -0.25308641975308643)
        self.assertAlmostEqual(kpoints["DELTA_0"][1], 0.25308641975308643)
        self.assertAlmostEqual(kpoints["DELTA_0"][2], 0.0)

        self.assertAlmostEqual(kpoints["F_0"][0], 0.25308641975308643)
        self.assertAlmostEqual(kpoints["F_0"][1], 0.7469135802469136)
        self.assertAlmostEqual(kpoints["F_0"][2], 0.0)

        self.assertAlmostEqual(kpoints["B_0"][0], -0.25308641975308643)
        self.assertAlmostEqual(kpoints["B_0"][1], 0.25308641975308643)
        self.assertAlmostEqual(kpoints["B_0"][2], 0.5)

        self.assertAlmostEqual(kpoints["B_2"][0], -0.25308641975308643)
        self.assertAlmostEqual(kpoints["B_2"][1], 0.25308641975308643)
        self.assertAlmostEqual(kpoints["B_2"][2], -0.5)

        self.assertAlmostEqual(kpoints["G_0"][0], 0.25308641975308643)
        self.assertAlmostEqual(kpoints["G_0"][1], 0.7469135802469136)
        self.assertAlmostEqual(kpoints["G_0"][2], 0.5)

        self.assertAlmostEqual(kpoints["G_2"][0], 0.25308641975308643)
        self.assertAlmostEqual(kpoints["G_2"][1], 0.7469135802469136)
        self.assertAlmostEqual(kpoints["G_2"][2], -0.5)


if __name__ == "__main__":
    unittest.main()
