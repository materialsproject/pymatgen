from __future__ import annotations

import pytest
from monty.io import zopen
from numpy.testing import assert_allclose

from pymatgen.core import Composition, Structure
from pymatgen.io.pwmat.inputs import (
    ACExtractor,
    ACstrExtractor,
    AtomConfig,
    GenKpt,
    HighSymmetryPoint,
    LineLocator,
    ListLocator,
)
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest

TEST_DIR = f"{TEST_FILES_DIR}/io/pwmat"


@pytest.mark.parametrize(
    ("exclusion", "expected_idx"),
    [("", 1), ("average", 163), ("AVERAGE", 163)],
)
def test_line_locator(exclusion: str, expected_idx: int):
    filepath = f"{TEST_DIR}/MOVEMENT.lzma"
    aim_idx = LineLocator.locate_all_lines(file_path=filepath, content="FORCE", exclusion=exclusion)[0]
    assert aim_idx == expected_idx


@pytest.mark.parametrize(
    ("exclusion", "expected_idx"),
    [("", 0), ("average", 1), ("AVERAGE", 1)],
)
def test_list_locator(exclusion: str, expected_idx: int):
    strs_lst = ["Average Force=  0.12342E+01", "Force"]
    aim_idx = ListLocator.locate_all_lines(strs_lst=strs_lst, content="FORCE", exclusion=exclusion)[0]
    assert aim_idx == expected_idx


class TestACstrExtractor(PymatgenTest):
    def test_extract(self):
        filepath = f"{TEST_DIR}/atom.config"
        ac_extractor = ACExtractor(file_path=filepath)
        with zopen(filepath, mode="rt") as file:
            ac_str_extractor = ACstrExtractor(atom_config_str="".join(file.readlines()))
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
        filepath = f"{TEST_DIR}/atom.config"
        structure = Structure.from_file(filepath)
        atom_config = AtomConfig(structure, sort_structure=False)
        assert atom_config.structure.composition == Composition("Cr2I6")

    def test_from_file(self):
        filepath = f"{TEST_DIR}/atom.config"
        atom_config = AtomConfig.from_file(filename=filepath, mag=True)
        assert Composition("Cr2I6").to_pretty_string() == atom_config.true_names
        assert all("magmom" in site.properties for site in atom_config.structure)

    def test_write_file(self):
        filepath = f"{TEST_DIR}/atom.config"
        atom_config = AtomConfig.from_file(filepath)
        tmp_file = f"{self.tmp_path}/atom.config.testing.lzma"
        atom_config.write_file(tmp_file)
        tmp_atom_config = AtomConfig.from_file(filepath)
        assert_allclose(atom_config.structure.lattice.abc, tmp_atom_config.structure.lattice.abc, 5)


class TestGenKpt(PymatgenTest):
    def test_from_structure(self):
        pytest.importorskip("seekpath")
        filepath = f"{TEST_DIR}/atom.config"
        structure = Structure.from_file(filepath)
        gen_kpt = GenKpt.from_structure(structure, dim=2, density=0.01)
        assert gen_kpt.density == pytest.approx(0.0628318530)
        assert gen_kpt.reciprocal_lattice.shape == (3, 3)
        assert gen_kpt.kpath["path"] == [["GAMMA", "M", "K", "GAMMA"]]

    def test_write_file(self):
        pytest.importorskip("seekpath")
        filepath = f"{TEST_DIR}/atom.config"
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
        pytest.importorskip("seekpath")
        filepath = f"{TEST_DIR}/atom.config"
        structure = Structure.from_file(filepath)
        high_symmetry_points = HighSymmetryPoint.from_structure(structure, dim=2, density=0.01)
        assert list(high_symmetry_points.kpath) == ["kpoints", "path"]
        assert len(high_symmetry_points.kpath["path"]) == 1
        assert high_symmetry_points.density == pytest.approx(0.0628318530)
        assert high_symmetry_points.reciprocal_lattice.shape == (3, 3)

    def test_write_file(self):
        pytest.importorskip("seekpath")
        filepath = f"{TEST_DIR}/atom.config"
        structure = Structure.from_file(filepath)
        dim = 2
        density = 0.01
        high_symmetry_points = HighSymmetryPoint.from_structure(structure, dim, density)
        tmp_filepath = f"{self.tmp_path}/HIGH_SYMMETRY_POINTS.testing.lzma"
        high_symmetry_points.write_file(tmp_filepath)
        tmp_high_symmetry_points_str = ""
        with zopen(tmp_filepath, "rt") as file:
            tmp_high_symmetry_points_str = file.read()
        assert tmp_high_symmetry_points_str == high_symmetry_points.get_str()


def test_err_msg_on_seekpath_not_installed(monkeypatch):
    """Simulate and test error message when seekpath is not installed."""
    try:
        import seekpath  # noqa: F401
    except ImportError:
        with pytest.raises(RuntimeError, match="SeeK-path needs to be installed to use the convention of Hinuma et al"):
            GenKpt.from_structure(Structure.from_file(f"{TEST_DIR}/atom.config"), dim=2, density=0.01)
