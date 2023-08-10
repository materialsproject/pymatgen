# ruff: noqa: PT011
from __future__ import annotations

import os
from pathlib import Path
from shutil import which
from subprocess import TimeoutExpired

import pytest

from pymatgen.analysis.molecule_matcher import MoleculeMatcher
from pymatgen.core import Molecule
from pymatgen.io.packmol import PackmolBoxGen
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest

test_dir = f"{TEST_FILES_DIR}/packmol"


if which("packmol") is None:
    pytest.skip("packmol executable not present", allow_module_level=True)


ethanol_coords = [
    [0.00720, -0.56870, 0.00000],
    [-1.28540, 0.24990, 0.00000],
    [1.13040, 0.31470, 0.00000],
    [0.03920, -1.19720, 0.89000],
    [0.03920, -1.19720, -0.89000],
    [-1.31750, 0.87840, 0.89000],
    [-1.31750, 0.87840, -0.89000],
    [-2.14220, -0.42390, -0.00000],
    [1.98570, -0.13650, -0.00000],
]
ethanol_atoms = ["C", "C", "O", "H", "H", "H", "H", "H", "H"]

ethanol = Molecule(ethanol_atoms, ethanol_coords)


water_coords = [
    [9.626, 6.787, 12.673],
    [9.626, 8.420, 12.673],
    [10.203, 7.604, 12.673],
]
water_atoms = ["H", "H", "O"]

water = Molecule(water_atoms, water_coords)


class TestPackmolSet(PymatgenTest):
    def test_packmol_with_molecule(self):
        """Test coords input as Molecule."""
        pw = PackmolBoxGen().get_input_set(
            molecules=[
                {"name": "water", "number": 10, "coords": water},
                {"name": "ethanol", "number": 20, "coords": ethanol},
            ],
        )
        pw.write_input(self.tmp_path)
        pw.run(self.tmp_path)
        assert os.path.exists(os.path.join(self.tmp_path, "packmol_out.xyz"))
        out = Molecule.from_file(os.path.join(self.tmp_path, "packmol_out.xyz"))
        assert out.composition.num_atoms == 10 * 3 + 20 * 9

    def test_packmol_with_str(self):
        """Test coords input as strings."""
        pw = PackmolBoxGen().get_input_set(
            molecules=[
                {"name": "EMC", "number": 10, "coords": f"{test_dir}/subdir with spaces/EMC.xyz"},
                {"name": "LiTFSi", "number": 20, "coords": f"{test_dir}/LiTFSi.xyz"},
            ],
        )
        pw.write_input(self.tmp_path)
        pw.run(self.tmp_path)
        assert os.path.exists(os.path.join(self.tmp_path, "packmol_out.xyz"))
        out = Molecule.from_file(os.path.join(self.tmp_path, "packmol_out.xyz"))
        assert out.composition.num_atoms == 10 * 15 + 20 * 16

    def test_packmol_with_path(self):
        """Test coords input as Path. Use a subdirectory with spaces."""
        p1 = Path(f"{test_dir}/subdir with spaces/EMC.xyz")
        p2 = Path(f"{test_dir}/LiTFSi.xyz")
        pw = PackmolBoxGen().get_input_set(
            molecules=[
                {"name": "EMC", "number": 10, "coords": p1},
                {"name": "LiTFSi", "number": 20, "coords": p2},
            ],
        )
        pw.write_input(self.tmp_path)
        pw.run(self.tmp_path)
        assert os.path.exists(os.path.join(self.tmp_path, "packmol_out.xyz"))
        out = Molecule.from_file(os.path.join(self.tmp_path, "packmol_out.xyz"))
        assert out.composition.num_atoms == 10 * 15 + 20 * 16

    def test_control_params(self):
        """
        Check that custom control_params work and that ValueError
        is raised when 'ERROR' appears in stdout (even if return code is 0).
        """
        input_set = PackmolBoxGen(
            control_params={"maxit": 0, "nloop": 0},
        ).get_input_set(
            molecules=[
                {"name": "water", "number": 1000, "coords": water},
                {"name": "ethanol", "number": 2000, "coords": ethanol},
            ],
        )
        input_set.write_input(self.tmp_path)
        with open(os.path.join(self.tmp_path, "packmol.inp")) as f:
            input_string = f.read()
            assert "maxit 0" in input_string
            assert "nloop 0" in input_string
        with pytest.raises(ValueError):
            input_set.run(self.tmp_path)

    def test_timeout(self):
        """Check that the timeout works."""
        pw = PackmolBoxGen().get_input_set(
            molecules=[
                {"name": "water", "number": 1000, "coords": water},
                {"name": "ethanol", "number": 2000, "coords": ethanol},
            ],
        )
        pw.write_input(self.tmp_path)
        with pytest.raises(TimeoutExpired):
            pw.run(self.tmp_path, 1)

    def test_no_return_and_box(self):
        """
        Make sure the code raises an error if packmol doesn't
        exit cleanly. Also verify the box arg works properly.
        """
        pw = PackmolBoxGen().get_input_set(
            molecules=[
                {"name": "water", "number": 1000, "coords": water},
                {"name": "ethanol", "number": 2000, "coords": ethanol},
            ],
            box=[0, 0, 0, 2, 2, 2],
        )
        pw.write_input(self.tmp_path)
        with open(os.path.join(self.tmp_path, "packmol.inp")) as f:
            input_string = f.read()
            assert "inside box 0 0 0 2 2 2" in input_string
        with pytest.raises(ValueError):
            pw.run(self.tmp_path)

    def test_chdir_behavior(self):
        """
        Make sure the code returns to the starting directory whether
        or not packmol exits cleanly.
        """
        startdir = str(Path.cwd())
        # this one will not exit cleanly b/c the box is too small
        pw = PackmolBoxGen().get_input_set(
            molecules=[
                {"name": "water", "number": 1000, "coords": water},
                {"name": "ethanol", "number": 2000, "coords": ethanol},
            ],
            box=[0, 0, 0, 2, 2, 2],
        )
        pw.write_input(self.tmp_path)
        with pytest.raises(ValueError):
            pw.run(self.tmp_path)
        assert str(Path.cwd()) == startdir

        # this one will exit cleanly
        pw = PackmolBoxGen().get_input_set(
            molecules=[
                {"name": "water", "number": 1000, "coords": water},
                {"name": "ethanol", "number": 2000, "coords": ethanol},
            ],
        )
        pw.write_input(self.tmp_path)
        pw.run(self.tmp_path)
        assert str(Path.cwd()) == startdir

    def test_random_seed(self):
        """
        Confirm that seed = -1 generates random structures
        while seed = 1 is deterministic.
        """
        pytest.importorskip("openbabel")
        mm = MoleculeMatcher()

        # deterministic output
        pw = PackmolBoxGen(seed=1, inputfile="input.in", outputfile="output.xyz").get_input_set(
            # self.tmp_path,
            molecules=[
                {"name": "water", "number": 10, "coords": water},
                {"name": "ethanol", "number": 20, "coords": ethanol},
            ],
        )
        pw.write_input(self.tmp_path)
        pw.run(self.tmp_path)
        out1 = Molecule.from_file(os.path.join(self.tmp_path, "output.xyz"))
        pw.run(self.tmp_path)
        out2 = Molecule.from_file(os.path.join(self.tmp_path, "output.xyz"))
        assert mm.fit(out1, out2)

        # randomly generated structures
        pw = PackmolBoxGen(seed=-1, inputfile="input.in", outputfile="output.xyz").get_input_set(
            molecules=[
                {"name": "water", "number": 10, "coords": water},
                {"name": "ethanol", "number": 20, "coords": ethanol},
            ],
        )
        pw.write_input(self.tmp_path)
        pw.run(self.tmp_path)
        out1 = Molecule.from_file(os.path.join(self.tmp_path, "output.xyz"))
        pw.run(self.tmp_path)
        out2 = Molecule.from_file(os.path.join(self.tmp_path, "output.xyz"))
        assert not mm.fit(out1, out2)

    def test_arbitrary_filenames(self):
        """
        Make sure custom input and output filenames work.
        Use a subdirectory with spaces.
        """
        os.mkdir(os.path.join(self.tmp_path, "subdirectory with spaces"))
        pw = PackmolBoxGen(
            inputfile="input.in", outputfile=Path("output.xyz"), stdoutfile=Path("stdout.txt")
        ).get_input_set(
            molecules=[
                {"name": "water", "number": 10, "coords": water},
                {"name": "ethanol", "number": 20, "coords": ethanol},
            ],
        )
        pw.write_input(
            os.path.join(self.tmp_path, "subdirectory with spaces"),
        )
        assert os.path.exists(os.path.join(self.tmp_path, "subdirectory with spaces", "input.in"))
        pw.run(
            os.path.join(self.tmp_path, "subdirectory with spaces"),
        )
        assert os.path.exists(os.path.join(self.tmp_path, "subdirectory with spaces", "output.xyz"))
        assert os.path.exists(os.path.join(self.tmp_path, "subdirectory with spaces", "stdout.txt"))
        out = Molecule.from_file(os.path.join(self.tmp_path, "subdirectory with spaces", "output.xyz"))
        assert out.composition.num_atoms == 10 * 3 + 20 * 9
