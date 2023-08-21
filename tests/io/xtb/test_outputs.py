from __future__ import annotations

import os

from pytest import approx

from pymatgen.core.structure import Molecule
from pymatgen.io.qchem.outputs import check_for_structure_changes
from pymatgen.io.xtb.outputs import CRESTOutput
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest

try:
    from openbabel import openbabel
except ImportError:
    openbabel = None
    print("OpenBabel not found, parsed molecules structures will not be  checked")

__author__ = "Alex Epstein"
__copyright__ = "Copyright 2020, The Materials Project"
__version__ = "0.1"

test_dir = f"{TEST_FILES_DIR}/xtb/sample_CREST_output"
expected_dir = f"{TEST_FILES_DIR}/xtb/expected_output"


class TestCRESTOutput(PymatgenTest):
    """
    Checks that all attributes of CRESTOutput match the expected values for a
    sample CREST output directory.
    """

    def test_all(self):
        expected_cmd_options = {"g": "H2O", "c": "2"}
        expected_energies = [["-13.66580"] * 10, ["-13.66479"] * 27]
        expected_sorted_structures = [[], []]
        for filepath in os.listdir(expected_dir):
            if filepath.endswith("xyz") and "_r" in filepath:
                n_conf = int(filepath.split("_")[0][-1])
                n_rot = int(filepath.split("_")[1].split(".")[0][-1])
                m = Molecule.from_file(os.path.join(expected_dir, filepath))
                expected_sorted_structures[n_conf].insert(n_rot, m)

        crest_out = CRESTOutput(output_filename="crest_out.out", path=test_dir)
        exp_best = Molecule.from_file(f"{expected_dir}/expected_crest_best.xyz")
        for i, c in enumerate(crest_out.sorted_structures_energies):
            for j, r in enumerate(c):
                if openbabel:
                    assert check_for_structure_changes(r[0], expected_sorted_structures[i][j]) == "no_change"
                assert float(r[1]) == approx(float(expected_energies[i][j]), abs=1e-4)

        assert crest_out.properly_terminated
        if openbabel:
            assert check_for_structure_changes(crest_out.lowest_energy_structure, exp_best) == "no_change"
        assert crest_out.cmd_options == expected_cmd_options
