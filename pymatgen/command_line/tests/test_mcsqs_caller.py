import os
import unittest

from monty.os.path import which
from monty.serialization import loadfn

from pymatgen.command_line.mcsqs_caller import run_mcsqs
from pymatgen.util.testing import PymatgenTest

__author__ = "Handong Ling, Rachel Woods-Robinson"
__maintainer__ = "Handong Ling, Rachel Woods-Robinson"
__email__ = "handongling@berkeley.edu, rwoodsrobinson@lbl.gov"

test_dir = os.path.join(
    os.path.dirname(__file__), "..", "..", "..", "test_files", "mcsqs"
)


@unittest.skipIf(not which("mcsqs"), "mcsqs executable not present")
class McsqsCallerTest(PymatgenTest):
    def setUp(self):
        self.pztstructs = loadfn(os.path.join(test_dir, "pztstructs.json"))
        self.pztstructs2 = loadfn(os.path.join(test_dir, "pztstructs2.json"))
        self.struc = self.get_structure("Pb2TiZrO6")
        self.perfect_match_zzn_rs = loadfn(
            os.path.join(test_dir, "perfect_match_zzn_rs.json")
        )

    def test_mcsqs_caller_supercell(self):
        struc = self.struc.copy()
        struc.replace_species(
            {"Ti": {"Ti": 0.5, "Zr": 0.5}, "Zr": {"Ti": 0.5, "Zr": 0.5}}
        )
        sqs = run_mcsqs(
            struc, {2: 6, 3: 4}, scaling=[2, 1, 1], search_time=0.01, instances=1
        )

        matches = [sqs.bestsqs.matches(s) for s in self.pztstructs]
        self.assertIn(True, matches)

    def test_mcsqs_caller_total_atoms(self):
        struc = self.struc.copy()
        struc.replace_species(
            {"Ti": {"Ti": 0.5, "Zr": 0.5}, "Zr": {"Ti": 0.5, "Zr": 0.5}}
        )
        sqs = run_mcsqs(struc, {2: 6, 3: 4}, scaling=2, search_time=0.01, instances=1)

        matches = [sqs.bestsqs.matches(s) for s in self.pztstructs2]
        self.assertIn(True, matches)

    def test_mcsqs_caller_total_atoms_auto_instances(self):
        struc = self.struc.copy()
        struc.replace_species(
            {"Ti": {"Ti": 0.5, "Zr": 0.5}, "Zr": {"Ti": 0.5, "Zr": 0.5}}
        )
        sqs = run_mcsqs(
            struc, {2: 6, 3: 4}, scaling=2, search_time=0.01, instances=None
        )

        matches = [sqs.bestsqs.matches(s) for s in self.pztstructs2]
        self.assertIn(True, matches)

    def test_mcsqs_caller_parallel(self):
        # explicitly test with four instances

        struc = self.struc.copy()
        struc.replace_species(
            {"Ti": {"Ti": 0.5, "Zr": 0.5}, "Zr": {"Ti": 0.5, "Zr": 0.5}}
        )
        sqs = run_mcsqs(struc, {2: 6, 3: 4}, scaling=2, search_time=0.01, instances=4)

        matches = [sqs.bestsqs.matches(s) for s in self.pztstructs2]
        self.assertIn(True, matches)

    def test_mcsqs_perfect_match_error(self):

        scale = 32 / self.perfect_match_zzn_rs.num_sites
        sqs = run_mcsqs(
            self.perfect_match_zzn_rs,
            {2: 6, 3: 4},
            scaling=scale,
            search_time=1,
            instances=1,
        )

        self.assertEqual(sqs.objective_function, "Perfect_match")

    def test_mcsqs_perfect_match_error_parallel(self):

        scale = 32 / self.perfect_match_zzn_rs.num_sites
        sqs = run_mcsqs(
            self.perfect_match_zzn_rs,
            {2: 6, 3: 4},
            scaling=scale,
            search_time=1,
            instances=4,
        )

        self.assertEqual(sqs.objective_function, "Perfect_match")

    def test_mcsqs_caller_runtime_error(self):
        struc = self.struc.copy()
        struc.replace_species(
            {"Ti": {"Ti": 0.5, "Zr": 0.5}, "Zr": {"Ti": 0.5, "Zr": 0.5}}
        )
        struc.replace_species({"Pb": {"Ti": 0.2, "Pb": 0.8}})
        struc.replace_species({"O": {"F": 0.8, "O": 0.2}})
        self.assertRaises(RuntimeError, run_mcsqs, struc, {2: 6, 3: 4}, 10, 0.000001)
