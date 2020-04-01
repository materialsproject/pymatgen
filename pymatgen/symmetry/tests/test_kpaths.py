# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


import unittest
import os

from pymatgen.util.testing import PymatgenTest
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.bandstructure import HighSymmKpath

from monty.serialization import loadfn

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "test_files")


class HighSymmKpathTest(PymatgenTest):
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
                    [[3.0233057319441246, 1, 0], [0, 7.9850357844548681, 1], [0, 1.2, 8.1136762279561818]]
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

            # Throws error if something doesn't work, causing test to fail.
            kpath = HighSymmKpath(struct, path_type="all")
            kpoints = kpath.get_kpoints()

    def test_continuous_kpath(self):
        bs = loadfn(os.path.join(test_dir, "Cu2O_361_bandstructure.json"))
        hskp = HighSymmKpath(bs.structure).get_continuous_path(bs)

        distance_map = [(3, False), (5, True), (1, True), (4, True), (3, True), (2, True), (1, True), (0, True)]
        labels = [
            ("\\Gamma", "R"),
            ("R", "M"),
            ("M", "X"),
            ("X", "R"),
            ("R", "\\Gamma"),
            ("\\Gamma", "M"),
            ("M", "X"),
            ("X", "\\Gamma"),
        ]

        self.assertEqual(hskp, (distance_map, labels))


if __name__ == "__main__":
    unittest.main()
