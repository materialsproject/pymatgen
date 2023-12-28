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
        with zopen(filepath, "rt") as f:
            ac_str_extractor = ACstrExtractor(atom_config_str="".join(f.readlines()))
        assert (ac_extractor.num_atoms == ac_str_extractor.get_num_atoms())
        for ii in range(9):
            assert (ac_extractor.lattice[ii] == ac_str_extractor.get_lattice()[ii])
        for ii in range(ac_extractor.num_atoms):
            assert (ac_extractor.types[ii] == ac_str_extractor.get_types()[ii])
            assert (ac_extractor.coords[ii*3 + 0] == ac_str_extractor.get_coords()[ii*3 + 0])
            assert (ac_extractor.coords[ii*3 + 1] == ac_str_extractor.get_coords()[ii*3 + 1])
            assert (ac_extractor.coords[ii*3 + 2] == ac_str_extractor.get_coords()[ii*3 + 2])
            assert (ac_extractor.magmoms[ii] == ac_str_extractor.get_magmoms()[ii])


class TestAtomConfig(PymatgenTest):
    def test_init(self):
        filepath = f"{TEST_FILES_DIR}/pwmat/atom.config"
        structure = Structure.from_file(filepath)
        atom_config = AtomConfig(structure, sort_structure=False)
        assert atom_config.structure.composition == Composition("Cr2I6")

    def test_from_file(self):
        filepath = f"{TEST_FILES_DIR}/pwmat/atom.config"
        atom_config = AtomConfig.from_file(filename=filepath, mag=True)
        assert (Composition("Cr2I6").formula == atom_config.true_names)
        for ii in range(8):
            assert "magmom" in atom_config.structure.sites[ii].properties

    def test_write_file(self):
        filepath = f"{TEST_FILES_DIR}/pwmat/atom.config"
        atom_config = AtomConfig.from_file(filepath)
        tmp_file = f"{TEST_FILES_DIR}/pwmat/atom.config.testing"
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
        tmp_file = f"{TEST_FILES_DIR}/pwmat/gen.kpt.testing"
        gen_kpt.write_file(tmp_file)
        tmp_gen_kpt_str = ""
        with zopen(tmp_file) as f:
            tmp_gen_kpt_str = f.read()
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
        tmp_filepath = f"{TEST_FILES_DIR}/pwmat/HIGH_SYMMETRY_POINTS.testing"
        high_symmetry_points.write_file(tmp_filepath)
        tmp_high_symmetry_points_str = ""
        with zopen(tmp_high_symmetry_points_str, "rt") as f:
            tmp_high_symmetry_points_str = f.read()
        assert tmp_high_symmetry_points_str == high_symmetry_points.get_str()
