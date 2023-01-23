# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


from __future__ import annotations

import os
import unittest

from pytest import approx

from pymatgen.core.structure import Molecule
from pymatgen.io.qchem.outputs import check_for_structure_changes
from pymatgen.io.xtb.outputs import CRESTOutput
from pymatgen.util.testing import PymatgenTest

try:
    from openbabel import openbabel

    openbabel  # reference openbabel so it's not unused import
    have_babel = True
except ImportError:
    have_babel = False
    print("OpenBabel not found, parsed molecules structures will not be  checked")

__author__ = "Alex Epstein"
__copyright__ = "Copyright 2020, The Materials Project"
__version__ = "0.1"

test_dir = os.path.join(os.path.dirname(__file__), PymatgenTest.TEST_FILES_DIR, "xtb", "sample_CREST_output")
expected_output_dir = os.path.join(os.path.dirname(__file__), PymatgenTest.TEST_FILES_DIR, "xtb", "expected_output")


class TestCRESTOutput(PymatgenTest):
    """
    Checks that all attributes of CRESTOutput match the expected values for a
    sample CREST output directory.
    """

    def test_all(self):
        expected_cmd_options = {"g": "H2O", "c": "2"}
        expected_energies = [
            [
                "-13.66580",
                "-13.66580",
                "-13.66580",
                "-13.66580",
                "-13.66580",
                "-13.66580",
                "-13.66580",
                "-13.66580",
                "-13.66580",
                "-13.66580",
            ],
            [
                "-13.66479",
                "-13.66479",
                "-13.66479",
                "-13.66479",
                "-13.66479",
                "-13.66479",
                "-13.66479",
                "-13.66479",
                "-13.66479",
                "-13.66479",
                "-13.66479",
                "-13.66479",
                "-13.66479",
                "-13.66479",
                "-13.66479",
                "-13.66479",
                "-13.66479",
                "-13.66479",
                "-13.66479",
                "-13.66479",
                "-13.66479",
                "-13.66479",
                "-13.66479",
                "-13.66479",
                "-13.66479",
                "-13.66479",
                "-13.66479",
            ],
        ]
        expected_sorted_structures = [[], []]
        for f in os.listdir(expected_output_dir):
            if f.endswith("xyz") and "_r" in f:
                n_conf = int(f.split("_")[0][-1])
                n_rot = int(f.split("_")[1].split(".")[0][-1])
                m = Molecule.from_file(os.path.join(expected_output_dir, f))
                expected_sorted_structures[n_conf].insert(n_rot, m)

        cout = CRESTOutput(output_filename="crest_out.out", path=test_dir)
        exp_best = Molecule.from_file(os.path.join(expected_output_dir, "expected_crest_best.xyz"))
        for i, c in enumerate(cout.sorted_structures_energies):
            for j, r in enumerate(c):
                if have_babel:
                    assert check_for_structure_changes(r[0], expected_sorted_structures[i][j]) == "no_change"
                assert float(r[1]) == approx(float(expected_energies[i][j]), abs=1e-4)

        assert cout.properly_terminated is True
        if have_babel:
            assert check_for_structure_changes(cout.lowest_energy_structure, exp_best) == "no_change"
        assert cout.cmd_options == expected_cmd_options


if __name__ == "__main__":
    unittest.main()
