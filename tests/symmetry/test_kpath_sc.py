from __future__ import annotations

from pytest import approx

from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.symmetry.kpath import KPathSetyawanCurtarolo
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest

test_dir_structs = f"{TEST_FILES_DIR}/space_group_structs"


class TestBandStructureSC(PymatgenTest):
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
                lattice = Lattice([[3.0233057319441246, 0, 0], [0, 7.9850357844548681, 0], [0, 0, 8.1136762279561818]])
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

        struct_file_path = f"{test_dir_structs}/ICSD_170.cif"
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

        assert list(kpoints["\\Gamma"]) == [0, 0, 0]

        assert kpoints["A"] == approx([0.25308642, 0.25308642, 0.5])

        assert kpoints["A_1"] == approx([-0.25308642, 0.74691358, 0.5])

        assert kpoints["R"] == approx([0, 0.5, 0.5])

        assert kpoints["S"] == approx([0, 0.5, 0])

        assert kpoints["T"] == approx([-0.5, 0.5, 0.5])

        assert kpoints["X"] == approx([0.25308642, 0.25308642, 0])

        assert kpoints["X_1"] == approx([-0.25308642, 0.74691358, 0])

        assert kpoints["Y"] == approx([-0.5, 0.5, 0])

        assert kpoints["Z"] == approx([0, 0, 0.5])
