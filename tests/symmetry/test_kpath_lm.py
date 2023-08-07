from __future__ import annotations

import numpy as np
from pytest import approx

from pymatgen.analysis.magnetism.analyzer import CollinearMagneticStructureAnalyzer
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.kpath import KPathLatimerMunro
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest

test_dir_structs = TEST_FILES_DIR


class TestKPathLatimerMunro(PymatgenTest):
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
        for sg_num in range(1, 231):
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
            _ = KPathLatimerMunro(struct)  # Throws error if something doesn't work, causing test to fail.

        struct_file_path = f"{test_dir_structs}/AgO_kpath_test.cif"
        struct = Structure.from_file(struct_file_path)
        _ = KPathLatimerMunro(struct)  # Throws error if something doesn't work, causing test to fail.

    def test_kpath_acentered(self):
        species = ["K", "La", "Ti"]
        coords = [[0.345, 5, 0.77298], [0.1345, 5.1, 0.77298], [0.7, 0.8, 0.9]]
        lattice = Lattice.orthorhombic(2, 9, 1)
        struct = Structure.from_spacegroup(38, lattice, species, coords)
        sga = SpacegroupAnalyzer(struct)
        struct_prim = sga.get_primitive_standard_structure(international_monoclinic=False)
        kpath = KPathLatimerMunro(struct_prim)

        kpoints = kpath._kpath["kpoints"]
        labels = list(kpoints)

        assert sorted(labels) == sorted(["a", "b", "c", "d", "d_{1}", "e", "f", "q", "q_{1}", "Γ"])

        assert kpoints["a"] == approx([0.0, 0.5, 0.0])

        assert kpoints["f"] == approx([-0.5, 0.5, 0.5])

        assert kpoints["c"] == approx([0.0, 0.0, 0.5])

        assert kpoints["b"] == approx([-0.5, 0.5, 0.0])

        assert kpoints["Γ"] == approx([0, 0, 0])

        assert kpoints["e"] == approx([0.0, 0.5, 0.5])

        assert kpoints["d_{1}"] == approx([0.2530864197530836, 0.25308641975308915, 0.0]) or kpoints["d"] == approx(
            [0.2530864197530836, 0.25308641975308915, 0.0]
        )

        assert kpoints["q_{1}"] == approx([0.2530864197530836, 0.25308641975308915, 0.5]) or kpoints["q"] == approx(
            [0.2530864197530836, 0.25308641975308915, 0.5]
        )

    def test_magnetic_kpath_generation(self):
        struct_file_path = f"{test_dir_structs}/LaMnO3_magnetic.mcif"
        struct = Structure.from_file(struct_file_path)
        mga = CollinearMagneticStructureAnalyzer(struct)
        col_spin_orig = mga.get_structure_with_spin()
        col_spin_orig.add_spin_by_site([0.0] * 20)
        col_spin_sym = SpacegroupAnalyzer(col_spin_orig)
        col_spin_prim = col_spin_sym.get_primitive_standard_structure(international_monoclinic=False)

        magmom_vec_list = [np.zeros(3) for site in col_spin_prim]
        magmom_vec_list[4:8] = [
            np.array([3.87, 3.87, 0.0]),
            np.array([3.87, 3.87, 0.0]),
            np.array([-3.87, -3.87, 0.0]),
            np.array([-3.87, -3.87, 0.0]),
        ]
        col_spin_prim.add_site_property("magmom", magmom_vec_list)

        kpath = KPathLatimerMunro(col_spin_prim, has_magmoms=True)

        kpoints = kpath._kpath["kpoints"]
        assert sorted(kpoints) == sorted(["a", "b", "c", "d", "d_{1}", "e", "f", "g", "g_{1}", "Γ"])

        assert kpoints["e"] == approx([-0.5, 0.0, 0.5])

        assert kpoints["g"] == approx([-0.5, -0.5, 0.5])

        assert kpoints["a"] == approx([-0.5, 0.0, 0.0])

        assert kpoints["g_{1}"] == approx([0.5, -0.5, 0.5])

        assert kpoints["f"] == approx([0.0, -0.5, 0.5])

        assert kpoints["c"] == approx([0.0, 0.0, 0.5])

        assert kpoints["b"] == approx([0.0, -0.5, 0.0])

        assert kpoints["Γ"] == approx([0, 0, 0])

        assert kpoints["d_{1}"] == approx([-0.5, -0.5, 0.0]) or kpoints["d"] == approx([-0.5, -0.5, 0.0])
