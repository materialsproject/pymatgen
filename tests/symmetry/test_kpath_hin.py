from __future__ import annotations

import unittest

from pytest import approx

from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.symmetry.kpath import KPathSeek
from pymatgen.util.testing import PymatgenTest

try:
    from seekpath import get_path
except ImportError:
    get_path = None


class TestKPathSeek(PymatgenTest):
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
        for idx in range(230):
            sg_num = idx + 1
            if sg_num in triclinic:
                lattice = Lattice([[3.02330573, 1, 0], [0, 7.98503578, 1], [0, 1.2, 8.11367622]])
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

        assert kpoints["GAMMA"] == approx([0.0, 0.0, 0.0])

        assert kpoints["Y"] == approx([0.5, 0.5, 0.0])

        assert kpoints["T"] == approx([0.5, 0.5, 0.5])

        assert kpoints["T_2"] == approx([0.5, 0.5, -0.5])

        assert kpoints["Z"] == approx([0.0, 0.0, 0.5])

        assert kpoints["Z_2"] == approx([0.0, 0.0, -0.5])

        assert kpoints["S"] == approx([0.0, 0.5, 0.0])

        assert kpoints["R"] == approx([0.0, 0.5, 0.5])

        assert kpoints["R_2"] == approx([0.0, 0.5, -0.5])

        assert kpoints["DELTA_0"] == approx([-0.25308641975308643, 0.25308641975308643, 0.0])

        assert kpoints["F_0"] == approx([0.25308641975308643, 0.7469135802469136, 0.0])

        assert kpoints["B_0"] == approx([-0.25308641975308643, 0.25308641975308643, 0.5])

        assert kpoints["B_2"] == approx([-0.25308641975308643, 0.25308641975308643, -0.5])

        assert kpoints["G_0"] == approx([0.25308641975308643, 0.7469135802469136, 0.5])

        assert kpoints["G_2"] == approx([0.25308641975308643, 0.7469135802469136, -0.5])
