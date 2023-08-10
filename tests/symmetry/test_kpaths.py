from __future__ import annotations

import unittest

from monty.serialization import loadfn

from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest

try:
    from seekpath import get_path
except ImportError:
    get_path = None


class TestHighSymmKpath(PymatgenTest):
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

            # Throws error if something doesn't work, causing test to fail.
            kpath = HighSymmKpath(struct, path_type="all")
            _ = kpath.get_kpoints()

    def test_continuous_kpath(self):
        bs = loadfn(f"{TEST_FILES_DIR}/Cu2O_361_bandstructure.json")
        cont_bs = loadfn(f"{TEST_FILES_DIR}/Cu2O_361_bandstructure_continuous.json.gz")
        alt_bs = HighSymmKpath(bs.structure).get_continuous_path(bs)

        assert cont_bs.as_dict() == alt_bs.as_dict()
        assert alt_bs.kpoints[0].label == alt_bs.kpoints[-1].label
