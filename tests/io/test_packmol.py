from __future__ import annotations

import os
from pathlib import Path
from shutil import which
from subprocess import TimeoutExpired

import numpy as np
import pytest

from pymatgen.analysis.molecule_matcher import MoleculeMatcher
from pymatgen.core import Molecule
from pymatgen.io.packmol import PackmolBoxGen
from pymatgen.util.testing import TEST_FILES_DIR, MatSciTest

TEST_DIR = f"{TEST_FILES_DIR}/io/packmol"
# error message is different in CI for unknown reasons (as of 2024-04-12)
# macOS: "Packmol failed with error code 173 and stderr: b'STOP 173\\n'"
# CI: "Packmol failed with return code 0 and stdout: Packmol was unable to
# put the molecules in the desired regions even without"
ERR_MSG_173 = "Packmol failed with "

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


helium = Molecule(["He"], [[0.0, 0.0, 0.0]])
xeon = Molecule(["Xe"], [[0.0, 0.0, 0.0]])


class TestPackmolSet(MatSciTest):
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
        assert os.path.isfile(f"{self.tmp_path}/packmol_out.xyz")
        out = Molecule.from_file(f"{self.tmp_path}/packmol_out.xyz")
        assert out.composition.num_atoms == 10 * 3 + 20 * 9

    def test_packmol_with_str(self):
        """Test coords input as strings."""
        pw = PackmolBoxGen().get_input_set(
            molecules=[
                {
                    "name": "EMC",
                    "number": 10,
                    "coords": f"{TEST_DIR}/subdir with spaces/EMC.xyz",
                },
                {"name": "LiTFSi", "number": 20, "coords": f"{TEST_DIR}/LiTFSi.xyz"},
            ],
        )
        pw.write_input(self.tmp_path)
        pw.run(self.tmp_path)
        assert os.path.isfile(f"{self.tmp_path}/packmol_out.xyz")
        out = Molecule.from_file(f"{self.tmp_path}/packmol_out.xyz")
        assert out.composition.num_atoms == 10 * 15 + 20 * 16

    def test_packmol_with_path(self):
        """Test coords input as Path. Use a subdirectory with spaces."""
        p1 = Path(f"{TEST_DIR}/subdir with spaces/EMC.xyz")
        p2 = Path(f"{TEST_DIR}/LiTFSi.xyz")
        pw = PackmolBoxGen().get_input_set(
            molecules=[
                {"name": "EMC", "number": 10, "coords": p1},
                {"name": "LiTFSi", "number": 20, "coords": p2},
            ],
        )
        # MatSciTest makes each test change to a temporary directory
        # Check here we can run in the current directory
        pw.write_input(".")
        pw.run(".")
        assert os.path.isfile(f"{self.tmp_path}/packmol_out.xyz")
        assert os.path.isfile(f"{self.tmp_path}/packmol.stdout")
        out = Molecule.from_file(f"{self.tmp_path}/packmol_out.xyz")
        assert out.composition.num_atoms == 10 * 15 + 20 * 16
        # MatSciTest makes each test change to a temporary directory
        # Check here we can run in a relative directory
        pw.write_input("somedir")
        pw.run("somedir")
        assert os.path.isfile(f"{self.tmp_path}/somedir/packmol_out.xyz")
        out = Molecule.from_file(f"{self.tmp_path}/somedir/packmol_out.xyz")
        assert os.path.isfile(f"{self.tmp_path}/somedir/packmol.stdout")
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
        with open(f"{self.tmp_path}/packmol.inp", encoding="utf-8") as file:
            input_string = file.read()
            assert "maxit 0" in input_string
            assert "nloop 0" in input_string
        with pytest.raises(ValueError, match=ERR_MSG_173):
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
        with open(f"{self.tmp_path}/packmol.inp", encoding="utf-8") as file:
            input_string = file.read()
            assert "inside box 0 0 0 2 2 2" in input_string
        with pytest.raises(ValueError, match=ERR_MSG_173):
            pw.run(self.tmp_path)

    def test_same_names(self):
        """Test failure when names are the same."""
        with pytest.raises(ValueError, match=r"Names of the molecules should all be different."):
            PackmolBoxGen().get_input_set(
                molecules=[
                    {"name": "water", "number": 10, "coords": water},
                    {"name": "water", "number": 20, "coords": ethanol},
                ],
            )

    def test_chdir_behavior(self):
        """
        Make sure the code returns to the starting directory whether
        or not packmol exits cleanly.
        """
        start_dir = str(Path.cwd())
        # this one will not exit cleanly b/c the box is too small
        pw = PackmolBoxGen().get_input_set(
            molecules=[
                {"name": "water", "number": 1000, "coords": water},
                {"name": "ethanol", "number": 2000, "coords": ethanol},
            ],
            box=[0, 0, 0, 2, 2, 2],
        )
        pw.write_input(self.tmp_path)

        with pytest.raises(ValueError, match=ERR_MSG_173):
            pw.run(self.tmp_path)
        assert str(Path.cwd()) == start_dir

        # this one will exit cleanly
        pw = PackmolBoxGen().get_input_set(
            molecules=[
                {"name": "water", "number": 1000, "coords": water},
                {"name": "ethanol", "number": 2000, "coords": ethanol},
            ],
        )
        pw.write_input(self.tmp_path)
        pw.run(self.tmp_path)
        assert str(Path.cwd()) == start_dir

    def test_random_seed(self):
        """
        Confirm that seed = -1 generates random structures
        while seed = 1 is deterministic.
        """
        pytest.importorskip("openbabel")
        mol_matcher = MoleculeMatcher()

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
        out1 = Molecule.from_file(f"{self.tmp_path}/output.xyz")
        pw.run(self.tmp_path)
        out2 = Molecule.from_file(f"{self.tmp_path}/output.xyz")
        assert mol_matcher.fit(out1, out2)

        # randomly generated structures
        pw = PackmolBoxGen(seed=-1, inputfile="input.in", outputfile="output.xyz").get_input_set(
            molecules=[
                {"name": "water", "number": 10, "coords": water},
                {"name": "ethanol", "number": 20, "coords": ethanol},
            ],
        )
        pw.write_input(self.tmp_path)
        pw.run(self.tmp_path)
        out1 = Molecule.from_file(f"{self.tmp_path}/output.xyz")
        pw.run(self.tmp_path)
        out2 = Molecule.from_file(f"{self.tmp_path}/output.xyz")
        assert not mol_matcher.fit(out1, out2)

    def test_arbitrary_filenames(self):
        """
        Make sure custom input and output filenames work.
        Use a subdirectory with spaces.
        """
        os.mkdir(f"{self.tmp_path}/subdirectory with spaces")
        pw = PackmolBoxGen(
            inputfile="input.in",
            outputfile=Path("output.xyz"),
            stdoutfile=Path("stdout.txt"),
        ).get_input_set(
            molecules=[
                {"name": "water", "number": 10, "coords": water},
                {"name": "ethanol", "number": 20, "coords": ethanol},
            ],
        )
        pw.write_input(
            f"{self.tmp_path}/subdirectory with spaces",
        )
        assert os.path.isfile(f"{self.tmp_path}/subdirectory with spaces/input.in")
        pw.run(
            f"{self.tmp_path}/subdirectory with spaces",
        )
        assert os.path.isfile(f"{self.tmp_path}/subdirectory with spaces/output.xyz")
        assert os.path.isfile(f"{self.tmp_path}/subdirectory with spaces/stdout.txt")
        out = Molecule.from_file(f"{self.tmp_path}/subdirectory with spaces/output.xyz")
        assert out.composition.num_atoms == 10 * 3 + 20 * 9

    def test_raise_individual_constraints(self):
        """
        Make sure it fails when providing individual constraints but not for all
        molecules.
        """
        with pytest.raises(
            ValueError,
            match=r"Either specify constraints for all molecules "
            r"or for none of them.",
        ):
            PackmolBoxGen(
                inputfile="input.in",
                outputfile=Path("output.xyz"),
                stdoutfile=Path("stdout.txt"),
            ).get_input_set(
                molecules=[
                    {"name": "water", "number": 10, "coords": water, "constraints": []},
                    {"name": "ethanol", "number": 20, "coords": ethanol, "constraints": ["inside sphere 1 2 3 4"]},
                ],
            )
        with pytest.raises(ValueError, match=r"Wrong packmol constraint: outsdie sphere 1 2 3 4"):
            PackmolBoxGen(
                inputfile="input.in",
                outputfile=Path("output.xyz"),
                stdoutfile=Path("stdout.txt"),
            ).get_input_set(
                molecules=[
                    {
                        "name": "water",
                        "number": 10,
                        "coords": water,
                        "constraints": ["outsdie sphere 1 2 3 4"],  # typo
                    },
                    {"name": "ethanol", "number": 20, "coords": ethanol, "constraints": ["inside sphere 1 2 3 4"]},
                ],
            )
        with pytest.raises(ValueError, match=r"Either specify individual constraints or box."):
            PackmolBoxGen(
                inputfile="input.in",
                outputfile=Path("output.xyz"),
                stdoutfile=Path("stdout.txt"),
            ).get_input_set(
                molecules=[
                    {"name": "water", "number": 10, "coords": water, "constraints": ["outside sphere 1 2 3 4"]},
                    {"name": "ethanol", "number": 20, "coords": ethanol, "constraints": ["inside sphere 1 2 3 4"]},
                ],
                box=[0, 0, 0, 5, 5, 5],
            )
        with pytest.raises(ValueError, match=r"You must use individual_constraints to use atoms_constraints."):
            PackmolBoxGen(
                inputfile="input.in",
                outputfile=Path("output.xyz"),
                stdoutfile=Path("stdout.txt"),
            ).get_input_set(
                molecules=[
                    {
                        "name": "water",
                        "number": 10,
                        "coords": water,
                        "atoms_constraints": [{"indices": [0, 2], "constraints": ["outside sphere 1 2 3 4"]}],
                    },
                    {"name": "ethanol", "number": 20, "coords": ethanol},
                ],
            )
        with pytest.raises(ValueError, match=r"Wrong packmol constraint for specific atoms: outsdie sphere 1 2 3 4"):
            PackmolBoxGen(
                inputfile="input.in",
                outputfile=Path("output.xyz"),
                stdoutfile=Path("stdout.txt"),
            ).get_input_set(
                molecules=[
                    {
                        "name": "water",
                        "number": 10,
                        "coords": water,
                        "constraints": ["inside cube 0 0 0 5"],
                        "atoms_constraints": [{"indices": [0, 2], "constraints": ["outsdie sphere 1 2 3 4"]}],
                    },
                    {"name": "ethanol", "number": 20, "coords": ethanol, "constraints": ["inside cube 0 0 0 5"]},
                ],
            )

    def test_individual_constraints(self):
        """
        Test individual constraints.
        """
        pw = PackmolBoxGen(
            inputfile="input.in",
            outputfile=Path("output.xyz"),
            stdoutfile=Path("stdout.txt"),
        ).get_input_set(
            molecules=[
                {"name": "helium", "number": 10, "coords": helium, "constraints": ["inside sphere 0 0 0 5"]},
                {
                    "name": "xeon",
                    "number": 100,
                    "coords": xeon,
                    "constraints": ["inside sphere 0 0 0 10", "outside sphere 0 0 0 5"],
                },
            ],
        )
        pw.write_input(f"{self.tmp_path}/with_constraints")
        assert os.path.isfile(f"{self.tmp_path}/with_constraints/input.in")
        with open(f"{self.tmp_path}/with_constraints/input.in") as f:
            inp = f.read()
        assert "  inside sphere 0 0 0 10\n  outside sphere 0 0 0 5" in inp
        pw.run(f"{self.tmp_path}/with_constraints")
        assert os.path.isfile(f"{self.tmp_path}/with_constraints/output.xyz")
        assert os.path.isfile(f"{self.tmp_path}/with_constraints/stdout.txt")
        out = Molecule.from_file(f"{self.tmp_path}/with_constraints/output.xyz")
        assert out.composition.num_atoms == 110
        for site in out:
            if site.specie.symbol == "He":
                assert np.linalg.norm(site.coords) <= 5.01
            else:
                assert 4.99 <= np.linalg.norm(site.coords) <= 10.01

    def test_atoms_constraints(self):
        """
        Test individual constraints.
        """
        pw = PackmolBoxGen(
            inputfile="input.in",
            outputfile=Path("output.xyz"),
            stdoutfile=Path("stdout.txt"),
            control_params={"precision": 0.001},
        ).get_input_set(
            molecules=[
                {
                    "name": "water",
                    "number": 5,
                    "coords": water,
                    "constraints": ["inside sphere 0 0 0 10"],
                    "atoms_constraints": [
                        {"indices": [2], "constraints": ["inside sphere 0 0 0 5"]},
                        {"indices": [0, 1], "constraints": ["outside sphere 0 0 0 5.5"]},
                    ],
                },
            ],
        )
        pw.write_input(f"{self.tmp_path}/with_atoms_constraints")
        assert os.path.isfile(f"{self.tmp_path}/with_atoms_constraints/input.in")
        with open(f"{self.tmp_path}/with_atoms_constraints/input.in") as f:
            inp = f.read()
        assert ("\n  atoms 3\n    inside sphere 0 0 0 5\n  end atoms\n") in inp
        assert ("\n  atoms 1 2\n    outside sphere 0 0 0 5.5\n  end atoms\n") in inp
        pw.run(f"{self.tmp_path}/with_atoms_constraints")
        assert os.path.isfile(f"{self.tmp_path}/with_atoms_constraints/output.xyz")
        assert os.path.isfile(f"{self.tmp_path}/with_atoms_constraints/stdout.txt")
        out = Molecule.from_file(f"{self.tmp_path}/with_atoms_constraints/output.xyz")
        assert out.composition.num_atoms == 15
        for site in out:
            if site.specie.symbol == "H":
                assert np.linalg.norm(site.coords) >= 5.49
            else:
                assert np.linalg.norm(site.coords) <= 5.01
