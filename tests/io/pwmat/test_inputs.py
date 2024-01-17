from __future__ import annotations

from monty.io import zopen
from numpy.testing import assert_allclose

from pymatgen.core import Composition, Structure
from pymatgen.io.pwmat.inputs import ACExtractor, ACstrExtractor, AtomConfig, GenKpt, HighSymmetryPoint
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest


class TestACstrExtractor(PymatgenTest):
    def test_extract(self):
        filepath = f"{TEST_FILES_DIR}/pwmat/atom.config"
        ac_extractor = ACExtractor(file_path=filepath)
        with zopen(filepath, mode="rt") as f:
            ac_str_extractor = ACstrExtractor(atom_config_str="".join(f.readlines()))
        assert ac_extractor.n_atoms == ac_str_extractor.get_n_atoms()
        for idx in range(9):
            assert ac_extractor.lattice[idx] == ac_str_extractor.get_lattice()[idx]
        for idx in range(ac_extractor.n_atoms):
            assert ac_extractor.types[idx] == ac_str_extractor.get_types()[idx]
            assert ac_extractor.coords[idx * 3 + 0] == ac_str_extractor.get_coords()[idx * 3 + 0]
            assert ac_extractor.coords[idx * 3 + 1] == ac_str_extractor.get_coords()[idx * 3 + 1]
            assert ac_extractor.coords[idx * 3 + 2] == ac_str_extractor.get_coords()[idx * 3 + 2]
            assert ac_extractor.magmoms[idx] == ac_str_extractor.get_magmoms()[idx]


class TestAtomConfig(PymatgenTest):
    def test_init(self):
        filepath = f"{TEST_FILES_DIR}/pwmat/atom.config"
        structure = Structure.from_file(filepath)
        atom_config = AtomConfig(structure, sort_structure=False)
        assert atom_config.structure.composition == Composition("Cr2I6")

    def test_from_file(self):
        filepath = f"{TEST_FILES_DIR}/pwmat/atom.config"
        atom_config = AtomConfig.from_file(filename=filepath, mag=True)
        # assert Composition("Cr2I6").formula == atom_config.true_names
        for ii in range(8):
            assert "magmom" in atom_config.structure.sites[ii].properties

    def test_write_file(self):
        filepath = f"{TEST_FILES_DIR}/pwmat/atom.config"
        atom_config = AtomConfig.from_file(filepath)
        tmp_file = f"{self.tmp_path}/atom.config.testing.lzma"
        atom_config.write_file(tmp_file)
        tmp_atom_config = AtomConfig.from_file(filepath)
        assert_allclose(atom_config.structure.lattice.abc, tmp_atom_config.structure.lattice.abc, 5)


class TestGenKpt(PymatgenTest):
    def test_from_structure(self):
        filepath = f"{TEST_FILES_DIR}/pwmat/atom.config"
        structure = Structure.from_file(filepath)
        dim = 2
        density = 0.01
        gen_kpt = GenKpt.from_structure(structure, dim, density)
        assert hasattr(gen_kpt, "_kpath")
        assert hasattr(gen_kpt, "_reciprocal_lattice")
        assert hasattr(gen_kpt, "_density")

    def test_write_file(self):
        filepath = f"{TEST_FILES_DIR}/pwmat/atom.config"
        structure = Structure.from_file(filepath)
        dim = 2
        density = 0.01
        gen_kpt = GenKpt.from_structure(structure, dim, density)
        tmp_file = f"{self.tmp_path}/gen.kpt.testing.lzma"
        gen_kpt.write_file(tmp_file)
        tmp_gen_kpt_str = ""
        with zopen(tmp_file, mode="rt") as file:
            tmp_gen_kpt_str = file.read()
        assert gen_kpt.get_str() == tmp_gen_kpt_str


class TestHighSymmetryPoint(PymatgenTest):
    def test_from_structure(self):
        filepath = f"{TEST_FILES_DIR}/pwmat/atom.config"
        structure = Structure.from_file(filepath)
        dim = 2
        density = 0.01
        high_symmetry_points = HighSymmetryPoint.from_structure(structure, dim, density)
        assert hasattr(high_symmetry_points, "_reciprocal_lattice")
        assert hasattr(high_symmetry_points, "_kpath")
        assert hasattr(high_symmetry_points, "_density")

    def test_write_file(self):
        filepath = f"{TEST_FILES_DIR}/pwmat/atom.config"
        structure = Structure.from_file(filepath)
        dim = 2
        density = 0.01
        high_symmetry_points = HighSymmetryPoint.from_structure(structure, dim, density)
        tmp_filepath = f"{self.tmp_path}/HIGH_SYMMETRY_POINTS.testing.lzma"
        high_symmetry_points.write_file(tmp_filepath)
        tmp_high_symmetry_points_str = ""
        with zopen(tmp_filepath, "rt") as f:
            tmp_high_symmetry_points_str = f.read()
        assert tmp_high_symmetry_points_str == high_symmetry_points.get_str()



class TestStructure(PymatgenTest):
    def test_to_from_file(self):
        filepath = f"{TEST_FILES_DIR}/pwmat/atom.config"
        struct = Structure.from_file(filepath)
        tmp_filepath = f"{self.tmp_path}/atom.config.testing.lzma"
        struct.to(tmp_filepath, fmt="pwmat")
