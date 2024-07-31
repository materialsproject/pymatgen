from __future__ import annotations

from pymatgen.core.structure import Molecule
from pymatgen.io.xtb.inputs import CRESTInput
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest

__author__ = "Alex Epstein"
__copyright__ = "Copyright 2020, The Materials Project"
__version__ = "0.1"

TEST_DIR = f"{TEST_FILES_DIR}/io/xtb/sample_CREST_output"
EXPECTED_DIR = f"{TEST_FILES_DIR}/io/xtb/expected_output"


class TestCRESTInput(PymatgenTest):
    """
    Checks that all attributes of CRESTInput match the expected values for
    sample inputs.
    """

    def test_coordinates_file(self):
        species = ["C", "O"]
        coords = [
            [-9.5782000000, 0.6241500000, 0.0000000000],
            [-7.5827400000, 0.5127000000, -0.0000000000],
        ]
        mol = Molecule(species=species, coords=coords)
        cin = CRESTInput(molecule=mol, coords_filename="crest_in.xyz")

        assert mol.as_dict() == cin.molecule.as_dict()
        assert cin.coords_filename == "crest_in.xyz"

    def test_constraints_file(self):
        constraints = {"atoms": [8, 1, 2], "force_constant": 0.5}
        mol = Molecule.from_file(f"{TEST_DIR}/crest_in.xyz")
        cin = CRESTInput(molecule=mol, constraints=constraints)
        with open(f"{EXPECTED_DIR}/expected_constrains.txt") as file:
            exp_con = file.read()
            assert (
                exp_con.strip()
                == cin.constrains_template(molecule=mol, reference_fnm="crest_in.xyz", constraints=constraints).strip()
            )
