from __future__ import annotations

from shutil import which

import pytest
from monty.serialization import loadfn

from pymatgen.command_line.mcsqs_caller import run_mcsqs
from pymatgen.core.structure import Structure
from pymatgen.util.testing import TEST_FILES_DIR, MatSciTest

__author__ = "Handong Ling, Rachel Woods-Robinson"
__maintainer__ = "Handong Ling, Rachel Woods-Robinson"
__email__ = "handongling@berkeley.edu, rwoodsrobinson@lbl.gov"


TEST_DIR = f"{TEST_FILES_DIR}/io/atat/mcsqs"


@pytest.mark.skipif(not (which("mcsqs") and which("str2cif")), reason="mcsqs executable not present")
class TestMcsqsCaller(MatSciTest):
    def setup_method(self):
        self.pzt_structs = loadfn(f"{TEST_DIR}/pzt-structs.json")
        self.pzt_structs2 = loadfn(f"{TEST_DIR}/pzt-structs-2.json")
        self.struct = self.get_structure("Pb2TiZrO6")
        self.perfect_match_zzn_rs = loadfn(f"{TEST_DIR}/perfect_match_zzn_rs.json")

    def test_mcsqs_caller_supercell(self):
        struct = self.struct.copy()
        struct.replace_species({"Ti": {"Ti": 0.5, "Zr": 0.5}, "Zr": {"Ti": 0.5, "Zr": 0.5}})
        sqs = run_mcsqs(struct, {2: 6, 3: 4}, scaling=[2, 1, 1], search_time=0.01, instances=1)

        matches = [sqs.bestsqs.matches(s) for s in self.pzt_structs]
        assert any(matches)

        assert isinstance(sqs.bestsqs, Structure)

        # ensures specific keys are present in cluster parsing for use in atomate
        assert set(sqs.clusters[0]) == {
            "multiplicity",
            "coordinates",
            "longest_pair_length",
            "num_points_in_cluster",
        }
        assert set(sqs.clusters[0]["coordinates"][0]) == {
            "cluster_function",
            "coordinates",
            "num_possible_species",
        }

    def test_mcsqs_caller_total_atoms(self):
        struct = self.struct.copy()
        struct.replace_species({"Ti": {"Ti": 0.5, "Zr": 0.5}, "Zr": {"Ti": 0.5, "Zr": 0.5}})
        sqs = run_mcsqs(struct, {2: 6, 3: 4}, scaling=2, search_time=0.01, instances=1)

        matches = [sqs.bestsqs.matches(s) for s in self.pzt_structs2]
        assert any(matches)

    def test_mcsqs_caller_total_atoms_auto_instances(self):
        struct = self.struct.copy()
        struct.replace_species({"Ti": {"Ti": 0.5, "Zr": 0.5}, "Zr": {"Ti": 0.5, "Zr": 0.5}})
        sqs = run_mcsqs(struct, {2: 6, 3: 4}, scaling=2, search_time=0.01, instances=None)

        matches = [sqs.bestsqs.matches(s) for s in self.pzt_structs2]
        assert any(matches)

    def test_mcsqs_caller_parallel(self):
        # explicitly test with four instances

        struct = self.struct.copy()
        struct.replace_species({"Ti": {"Ti": 0.5, "Zr": 0.5}, "Zr": {"Ti": 0.5, "Zr": 0.5}})
        sqs = run_mcsqs(struct, {2: 6, 3: 4}, scaling=2, search_time=0.01, instances=4)

        matches = [sqs.bestsqs.matches(s) for s in self.pzt_structs2]
        assert any(matches)

    def test_mcsqs_perfect_match_error(self):
        scale = 32 / len(self.perfect_match_zzn_rs)
        sqs = run_mcsqs(
            self.perfect_match_zzn_rs,
            {2: 6, 3: 4},
            scaling=scale,
            search_time=1,
            instances=1,
        )

        assert sqs.objective_function == "Perfect_match"

    def test_mcsqs_perfect_match_error_parallel(self):
        scale = 32 / len(self.perfect_match_zzn_rs)
        sqs = run_mcsqs(
            self.perfect_match_zzn_rs,
            {2: 6, 3: 4},
            scaling=scale,
            search_time=1,
            instances=4,
        )

        assert sqs.objective_function == "Perfect_match"

    @pytest.mark.xfail(reason="Fix later #4455")
    def test_mcsqs_caller_runtime_error(self):
        struct = self.struct.copy()
        struct.replace_species({"Ti": {"Ti": 0.5, "Zr": 0.5}, "Zr": {"Ti": 0.5, "Zr": 0.5}})
        struct.replace_species({"Pb": {"Ti": 0.2, "Pb": 0.8}})
        struct.replace_species({"O": {"F": 0.8, "O": 0.2}})
        with pytest.raises(RuntimeError, match="mcsqs exited before timeout reached"):
            run_mcsqs(struct, {2: 6, 3: 4}, 10, search_time=0.000001)
