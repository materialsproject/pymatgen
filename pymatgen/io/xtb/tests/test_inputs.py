# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


import os

from pymatgen.core.structure import Molecule
from pymatgen.io.xtb.inputs import CRESTInput
from pymatgen.util.testing import PymatgenTest

__author__ = "Alex Epstein"
__copyright__ = "Copyright 2020, The Materials Project"
__version__ = "0.1"

expected_output_dir = os.path.join(os.path.dirname(__file__), PymatgenTest.TEST_FILES_DIR, "xtb", "expected_output")
test_dir = os.path.join(os.path.dirname(__file__), PymatgenTest.TEST_FILES_DIR, "xtb", "sample_CREST_output")


class TestCRESTInput(PymatgenTest):
    """
    Checks that all attributes of CRESTInput match the expected values for
    sample inputs
    """

    def test_coordinates_file(self):
        species = ["C", "O"]
        coords = [
            [-9.5782000000, 0.6241500000, 0.0000000000],
            [-7.5827400000, 0.5127000000, -0.0000000000],
        ]
        mol = Molecule(species=species, coords=coords)
        cin = CRESTInput(molecule=mol, coords_filename="crest_in.xyz")

        self.assertDictEqual(mol.as_dict(), cin.molecule.as_dict())
        self.assertEqual("crest_in.xyz", cin.coords_filename)

    def test_constraints_file(self):
        constraints = {"atoms": [8, 1, 2], "force_constant": 0.5}
        mol = Molecule.from_file(os.path.join(test_dir, "crest_in.xyz"))
        cin = CRESTInput(molecule=mol, constraints=constraints)
        with open(os.path.join(expected_output_dir, "expected_constrains.txt")) as f:
            exp_con = f.read()
            self.assertEqual(
                exp_con.strip(),
                cin.constrains_template(molecule=mol, reference_fnm="crest_in.xyz", constraints=constraints).strip(),
            )
