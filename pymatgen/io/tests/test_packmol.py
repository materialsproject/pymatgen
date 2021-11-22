import os
import tempfile
from pathlib import Path
from subprocess import TimeoutExpired

import pytest
from pymatgen.analysis.molecule_matcher import MoleculeMatcher
from pymatgen.core import Molecule
from pymatgen.util.testing import PymatgenTest

from pymatgen.io.packmol import PackmolBoxGen

test_dir = os.path.join(PymatgenTest.TEST_FILES_DIR, "packmol")


@pytest.fixture
def ethanol():
    """
    Returns a Molecule of ethanol
    """
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

    return Molecule(ethanol_atoms, ethanol_coords)


@pytest.fixture
def water():
    """
    Returns a Molecule of water
    """
    water_coords = [
        [9.626, 6.787, 12.673],
        [9.626, 8.420, 12.673],
        [10.203, 7.604, 12.673],
    ]
    water_atoms = ["H", "H", "O"]

    return Molecule(water_atoms, water_coords)


class TestPackmolSet:
    def test_packmol_with_molecule(self, water, ethanol):
        """
        Test coords input as Molecule
        """
        with tempfile.TemporaryDirectory() as scratch_dir:
            pw = PackmolBoxGen().get_input_set(
                molecules=[
                    {"name": "water", "number": 10, "coords": water},
                    {"name": "ethanol", "number": 20, "coords": ethanol},
                ],
            )
            pw.write_input(scratch_dir)
            pw.run(scratch_dir)
            assert os.path.exists(os.path.join(scratch_dir, "packmol_out.xyz"))
            out = Molecule.from_file(os.path.join(scratch_dir, "packmol_out.xyz"))
            assert out.composition.num_atoms == 10 * 3 + 20 * 9

    def test_packmol_with_str(self):
        """
        Test coords input as strings
        """
        with tempfile.TemporaryDirectory() as scratch_dir:
            pw = PackmolBoxGen().get_input_set(
                molecules=[
                    {"name": "EMC", "number": 10, "coords": os.path.join(test_dir, "subdir with spaces", "EMC.xyz")},
                    {"name": "LiTFSi", "number": 20, "coords": os.path.join(test_dir, "LiTFSi.xyz")},
                ],
            )
            pw.write_input(scratch_dir)
            pw.run(scratch_dir)
            assert os.path.exists(os.path.join(scratch_dir, "packmol_out.xyz"))
            out = Molecule.from_file(os.path.join(scratch_dir, "packmol_out.xyz"))
            assert out.composition.num_atoms == 10 * 15 + 20 * 16

    def test_packmol_with_path(self):
        """
        Test coords input as Path. Use a subdirectory with spaces.
        """
        p1 = Path(os.path.join(test_dir, "subdir with spaces", "EMC.xyz"))
        p2 = Path(os.path.join(test_dir, "LiTFSi.xyz"))
        with tempfile.TemporaryDirectory() as scratch_dir:
            pw = PackmolBoxGen().get_input_set(
                molecules=[
                    {"name": "EMC", "number": 10, "coords": p1},
                    {"name": "LiTFSi", "number": 20, "coords": p2},
                ],
            )
            pw.write_input(scratch_dir)
            pw.run(scratch_dir)
            assert os.path.exists(os.path.join(scratch_dir, "packmol_out.xyz"))
            out = Molecule.from_file(os.path.join(scratch_dir, "packmol_out.xyz"))
            assert out.composition.num_atoms == 10 * 15 + 20 * 16

    def test_control_params(self, water, ethanol):
        """
        Check that custom control_params work and that ValueError
        is raised when 'ERROR' appears in stdout (even if return code is 0)
        """
        with tempfile.TemporaryDirectory() as scratch_dir:
            pw = PackmolBoxGen(control_params={"maxit": 0, "nloop": 0},).get_input_set(
                molecules=[
                    {"name": "water", "number": 1000, "coords": water},
                    {"name": "ethanol", "number": 2000, "coords": ethanol},
                ],
            )
            pw.write_input(scratch_dir)
            with open(os.path.join(scratch_dir, "packmol.inp"), "r") as f:
                input_string = f.read()
                assert "maxit 0" in input_string
                assert "nloop 0" in input_string
            with pytest.raises(ValueError):
                pw.run(scratch_dir)

    def test_timeout(self, water, ethanol):
        """
        Check that the timeout works
        """
        with tempfile.TemporaryDirectory() as scratch_dir:
            pw = PackmolBoxGen().get_input_set(
                molecules=[
                    {"name": "water", "number": 1000, "coords": water},
                    {"name": "ethanol", "number": 2000, "coords": ethanol},
                ],
            )
            pw.write_input(scratch_dir)
            with pytest.raises(TimeoutExpired):
                pw.run(scratch_dir, 1)

    def test_no_return_and_box(self, water, ethanol):
        """
        Make sure the code raises an error if packmol doesn't
        exit cleanly. Also verify the box arg works properly.
        """
        with tempfile.TemporaryDirectory() as scratch_dir:
            pw = PackmolBoxGen().get_input_set(
                molecules=[
                    {"name": "water", "number": 1000, "coords": water},
                    {"name": "ethanol", "number": 2000, "coords": ethanol},
                ],
                box=[0, 0, 0, 2, 2, 2],
            )
            pw.write_input(scratch_dir)
            with open(os.path.join(scratch_dir, "packmol.inp"), "r") as f:
                input_string = f.read()
                assert "inside box 0 0 0 2 2 2" in input_string
            with pytest.raises(ValueError):
                pw.run(scratch_dir)

    def test_chdir_behavior(self, water, ethanol):
        """
        Make sure the code returns to the starting directory whether
        or not packmol exits cleanly.
        """
        startdir = str(Path.cwd())
        with tempfile.TemporaryDirectory() as scratch_dir:
            # this one will not exit cleanly b/c the box is too small
            pw = PackmolBoxGen().get_input_set(
                molecules=[
                    {"name": "water", "number": 1000, "coords": water},
                    {"name": "ethanol", "number": 2000, "coords": ethanol},
                ],
                box=[0, 0, 0, 2, 2, 2],
            )
            pw.write_input(scratch_dir)
            with pytest.raises(ValueError):
                pw.run(scratch_dir)
            assert str(Path.cwd()) == startdir

        with tempfile.TemporaryDirectory() as scratch_dir:
            # this one will exit cleanly
            pw = PackmolBoxGen().get_input_set(
                molecules=[
                    {"name": "water", "number": 1000, "coords": water},
                    {"name": "ethanol", "number": 2000, "coords": ethanol},
                ],
            )
            pw.write_input(scratch_dir)
            pw.run(scratch_dir)
            assert str(Path.cwd()) == startdir

    def test_random_seed(self, water, ethanol):
        """
        Confirm that seed = -1 generates random structures
        while seed = 1 is deterministic
        """
        mm = MoleculeMatcher()

        # deterministic output
        with tempfile.TemporaryDirectory() as scratch_dir:
            pw = PackmolBoxGen(seed=1, inputfile="input.in", outputfile="output.xyz",).get_input_set(
                # scratch_dir,
                molecules=[
                    {"name": "water", "number": 10, "coords": water},
                    {"name": "ethanol", "number": 20, "coords": ethanol},
                ],
            )
            pw.write_input(scratch_dir)
            pw.run(scratch_dir)
            out1 = Molecule.from_file(os.path.join(scratch_dir, "output.xyz"))
            pw.run(scratch_dir)
            out2 = Molecule.from_file(os.path.join(scratch_dir, "output.xyz"))
            assert mm.fit(out1, out2)

        # randomly generated structures
        with tempfile.TemporaryDirectory() as scratch_dir:
            pw = PackmolBoxGen(seed=-1, inputfile="input.in", outputfile="output.xyz",).get_input_set(
                molecules=[
                    {"name": "water", "number": 10, "coords": water},
                    {"name": "ethanol", "number": 20, "coords": ethanol},
                ],
            )
            pw.write_input(scratch_dir)
            pw.run(scratch_dir)
            out1 = Molecule.from_file(os.path.join(scratch_dir, "output.xyz"))
            pw.run(scratch_dir)
            out2 = Molecule.from_file(os.path.join(scratch_dir, "output.xyz"))
            assert not mm.fit(out1, out2)

    def test_arbitrary_filenames(self, water, ethanol):
        """
        Make sure custom input and output filenames work.
        Use a subdirectory with spaces.
        """
        with tempfile.TemporaryDirectory() as scratch_dir:
            os.mkdir(os.path.join(scratch_dir, "subdirectory with spaces"))
            pw = PackmolBoxGen(
                inputfile="input.in", outputfile=Path("output.xyz"), stdoutfile=Path("stdout.txt")
            ).get_input_set(
                molecules=[
                    {"name": "water", "number": 10, "coords": water},
                    {"name": "ethanol", "number": 20, "coords": ethanol},
                ],
            )
            pw.write_input(
                os.path.join(scratch_dir, "subdirectory with spaces"),
            )
            assert os.path.exists(os.path.join(scratch_dir, "subdirectory with spaces", "input.in"))
            pw.run(
                os.path.join(scratch_dir, "subdirectory with spaces"),
            )
            assert os.path.exists(os.path.join(scratch_dir, "subdirectory with spaces", "output.xyz"))
            assert os.path.exists(os.path.join(scratch_dir, "subdirectory with spaces", "stdout.txt"))
            out = Molecule.from_file(os.path.join(scratch_dir, "subdirectory with spaces", "output.xyz"))
            assert out.composition.num_atoms == 10 * 3 + 20 * 9
