from __future__ import annotations

import os

from pytest import approx

from pymatgen.core.structure import Molecule
from pymatgen.io.qchem.outputs import check_for_structure_changes
from pymatgen.io.xtb.outputs import CRESTOutput
from pymatgen.util.testing import TEST_FILES_DIR, MatSciTest

try:
    from openbabel import openbabel
except ImportError:
    openbabel = None
    print("OpenBabel not installed, parsed molecules structures will not be checked")

__author__ = "Alex Epstein"
__copyright__ = "Copyright 2020, The Materials Project"
__version__ = "0.1"

TEST_DIR = f"{TEST_FILES_DIR}/io/xtb/sample_CREST_output"
EXPECTED_DIR = f"{TEST_FILES_DIR}/io/xtb/expected_output"


class TestCRESTOutput(MatSciTest):
    """
    Checks that all attributes of CRESTOutput match the expected values for a
    sample CREST output directory.
    """

    def test_all(self):
        expected_cmd_options = {"g": "H2O", "c": "2"}
        expected_energies = [["-13.66580"] * 10, ["-13.66479"] * 27]
        expected_sorted_structures = [[], []]
        for filepath in os.listdir(EXPECTED_DIR):
            if filepath.endswith("xyz") and "_r" in filepath:
                n_conf = int(filepath.split("_")[0][-1])
                n_rot = int(filepath.split("_")[1].split(".")[0][-1])
                mol = Molecule.from_file(f"{EXPECTED_DIR}/{filepath}")
                expected_sorted_structures[n_conf].insert(n_rot, mol)

        crest_out = CRESTOutput(output_filename="crest_out.out", path=TEST_DIR)
        exp_best = Molecule.from_file(f"{EXPECTED_DIR}/expected_crest_best.xyz")
        for idx, c in enumerate(crest_out.sorted_structures_energies):
            for j, r in enumerate(c):
                if openbabel is not None:
                    assert check_for_structure_changes(r[0], expected_sorted_structures[idx][j]) == "no_change"
                assert float(r[1]) == approx(float(expected_energies[idx][j]), abs=1e-4)

        assert crest_out.properly_terminated
        if openbabel is not None:
            assert check_for_structure_changes(crest_out.lowest_energy_structure, exp_best) == "no_change"
        assert crest_out.cmd_options == expected_cmd_options
