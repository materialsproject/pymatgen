from __future__ import annotations

import filecmp
import os
import re
import shutil
import unittest

import pandas as pd
import pytest

from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.io.lammps.data import LammpsData
from pymatgen.io.lammps.inputs import LammpsInputFile, LammpsRun, LammpsTemplateGen, write_lammps_inputs
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest

test_dir = f"{TEST_FILES_DIR}/lammps"


class TestLammpsInputFile(PymatgenTest):
    @classmethod
    def setUpClass(cls):
        cls.filename = f"{test_dir}/lgps.in"

    def test_from_file(self):
        lmp_input = LammpsInputFile().from_file(self.filename)
        assert lmp_input.ncomments == 3
        assert lmp_input.nstages == 1
        assert lmp_input.stages_names == ["Stage 1"]
        assert lmp_input.stages[0]["commands"][0][1] == "metal"
        assert lmp_input.stages == [
            {
                "stage_name": "Stage 1",
                "commands": [
                    ("units", "metal"),
                    ("atom_style", "full"),
                    ("dimension", "3"),
                    ("pair_style", "hybrid/overlay morse 15 coul/long 15"),
                    ("kspace_style", "ewald 1e-4"),
                    ("boundary", "p p p"),
                    ("#", "2) System definition"),
                    ("read_data", "run_init.data"),
                    ("set", "type 1 charge 0.8803"),
                    ("set", "type 2 charge 1.2570"),
                    ("set", "type 3 charge 1.2580"),
                    ("set", "type 4 charge -1.048"),
                    ("neigh_modify", "every 1 delay 5 check yes"),
                    ("#", "3) Simulation settings"),
                    ("pair_coeff", "1 1 morse 0.0580 3.987 3.404"),
                    ("pair_coeff", "1 4 morse 0.0408 1.399 3.204"),
                    ("pair_coeff", "2 4 morse 0.3147 2.257 2.409"),
                    ("pair_coeff", "3 4 morse 0.4104 2.329 2.200"),
                    ("pair_coeff", "4 4 morse 0.0241 1.359 4.284"),
                    ("pair_coeff", "* * coul/long"),
                    ("#", "Part A : energy minimization"),
                    ("thermo", "1"),
                    ("thermo_style", "custom step lx ly lz press pxx pyy pzz pe"),
                    ("dump", "dmp all atom 5 run.dump"),
                    ("min_style", "cg"),
                    ("fix", "1 all box/relax iso 0.0 vmax 0.001"),
                    ("minimize", "1.0e-16 1.0e-16 5000 10000"),
                    ("write_data", "run.data"),
                ],
            }
        ]

    def test_get_string(self):
        lmp_input = LammpsInputFile().from_file(self.filename)
        string = lmp_input.get_string()
        assert "# LAMMPS input generated from LammpsInputFile with pymatgen v" in string
        assert "\nunits metal\natom_style full\ndimension 3\npair_style hybrid/overlay morse 15 coul/long" in string
        assert "15\nkspace_style ewald 1e-4\nboundary p p p\n# 2) System definition" in string
        assert "\nread_data run_init.data\nset type 1 charge 0.8803\nset type 2 charge" in string
        assert "1.2570\nset type 3 charge 1.2580\nset type 4 charge -1.048\nneigh_modify" in string
        assert "every 1 delay 5 check yes\n# 3) Simulation settings\npair_coeff 1 1" in string
        assert "morse 0.0580 3.987 3.404\npair_coeff 1 4 morse 0.0408 1.399 3.204\npair_coeff" in string
        assert "2 4 morse 0.3147 2.257 2.409\npair_coeff 3 4 morse 0.4104 2.329 2.200\npair_coeff" in string
        assert "4 4 morse 0.0241 1.359 4.284\npair_coeff * * coul/long\n# Part A : energy" in string
        assert "minimization\nthermo 1\nthermo_style custom step lx ly lz press pxx pyy pzz" in string
        assert "pe\ndump dmp all atom 5 run.dump\nmin_style cg\nfix 1 all" in string
        assert "box/relax iso 0.0 vmax 0.001\nminimize 1.0e-16 1.0e-16 5000 10000\nwrite_data run.data" in string

    def test_from_string(self):
        string = """# LGPS

# 1) Initialization
units metal
atom_style full
dimension 3
pair_style hybrid/overlay &
morse 15 coul/long 15
kspace_style ewald 1e-4
boundary p p p

# 2) System definition
read_data run_init.data
set type 1 charge  0.8803
set type 2 charge  1.2570
set type 3 charge  1.2580
set type 4 charge -1.048
neigh_modify every 1 delay 5 check yes

# 3) Simulation settings
pair_coeff 1 1 morse 0.0580 3.987 3.404
pair_coeff 1 4 morse 0.0408 1.399 3.204
pair_coeff 2 4 morse 0.3147 2.257 2.409
pair_coeff 3 4 morse 0.4104 2.329 2.200
pair_coeff 4 4 morse 0.0241 1.359 4.284
pair_coeff * * coul/long

# Part A : energy minimization
thermo 1
thermo_style custom step lx ly lz press pxx pyy pzz pe
dump dmp all atom 5 run.dump

min_style cg
fix 1 all box/relax iso 0.0 vmax 0.001
minimize 1.0e-16 1.0e-16 5000 10000
write_data run.data"""

        lmp_input = LammpsInputFile().from_str(string)
        assert lmp_input.stages == LammpsInputFile().from_file(self.filename).stages

    def test_stages_names(self):
        lmp_input = LammpsInputFile().from_file(self.filename)
        assert lmp_input.stages_names == ["Stage 1"]

        lmp_input = LammpsInputFile().from_file(self.filename, keep_stages=True)
        assert lmp_input.stages_names == [
            "Comment 1",
            "1) Initialization",
            "2) System definition",
            "3) Simulation settings",
            "Part A : energy minimization",
            "Stage 6",
        ]

        lmp_input.stages.append({"stage_name": "New stage", "commands": [("units", "metal")]})
        assert lmp_input.stages_names == [
            "Comment 1",
            "1) Initialization",
            "2) System definition",
            "3) Simulation settings",
            "Part A : energy minimization",
            "Stage 6",
            "New stage",
        ]

    def test_nstages(self):
        lmp_input = LammpsInputFile().from_file(self.filename)
        assert lmp_input.nstages == 1

        lmp_input = LammpsInputFile().from_file(self.filename, keep_stages=True)
        assert lmp_input.nstages == 6

    def test_ncomments(self):
        lmp_input = LammpsInputFile().from_file(self.filename)
        assert lmp_input.ncomments == 3

        lmp_input = LammpsInputFile().from_file(self.filename, keep_stages=True)
        assert lmp_input.ncomments == 1

    def test_get_args(self):
        lmp_input = LammpsInputFile().from_file(self.filename)
        units = lmp_input.get_args("units")
        assert units == "metal"

        sets = lmp_input.get_args("set")
        assert sets == ["type 1 charge 0.8803", "type 2 charge 1.2570", "type 3 charge 1.2580", "type 4 charge -1.048"]

        lmp_input = LammpsInputFile().from_file(self.filename, keep_stages=True)
        set1 = lmp_input.get_args("set", stage_name="2) System definition")
        assert set1[1] == "type 2 charge 1.2570"

    def test_contains_command(self):
        lmp_input = LammpsInputFile().from_file(self.filename)
        assert lmp_input.contains_command("set")
        assert not lmp_input.contains_command("sett")

        lmp_input = LammpsInputFile().from_file(self.filename, keep_stages=True)
        assert lmp_input.contains_command("units")
        assert lmp_input.contains_command("units", stage_name="1) Initialization")
        assert not lmp_input.contains_command("units", stage_name="Stage 6")

    def test_set_args(self):
        lmp_input = LammpsInputFile().from_file(self.filename)
        lmp_input.set_args(command="units", argument="atomic")
        assert lmp_input.get_args("units") == "atomic"

        lmp_input2 = lmp_input
        lmp_input2.set_args(command="set", argument="new set 2", how=2, stage_name="Stage 1")
        assert lmp_input.stages == lmp_input2.stages
        lmp_input2.set_args(command="set", argument="new set", how="first")
        assert lmp_input.get_args("set")[0] == "new set"
        assert lmp_input.get_args("set")[1] == "type 2 charge 1.2570"

        lmp_input3 = LammpsInputFile().from_file(self.filename, keep_stages=True)
        lmp_input4 = LammpsInputFile().from_file(self.filename, keep_stages=True)
        lmp_input4.set_args(command="set", argument="new set 2", how=2, stage_name="Stage 1")
        assert lmp_input3.stages == lmp_input4.stages
        lmp_input4.set_args(command="set", argument="new set 2", how="first", stage_name="2) System definition")
        assert lmp_input4.get_args("set", stage_name="2) System definition")[0] == "new set 2"
        assert lmp_input4.get_args("set", stage_name="2) System definition")[1] == "type 2 charge 1.2570"

    def test_add_stage(self):
        lmp_input = LammpsInputFile()
        lmp_input.add_stage(commands="units metal")
        assert lmp_input.stages == [{"stage_name": "Stage 1", "commands": [("units", "metal")]}]
        lmp_input.add_stage(commands={"pair_style": "eam"}, stage_name="Pair style")
        assert lmp_input.stages == [
            {"stage_name": "Stage 1", "commands": [("units", "metal")]},
            {"stage_name": "Pair style", "commands": [("pair_style", "eam")]},
        ]
        lmp_input.add_stage(commands=["boundary p p p", "atom_style full"], stage_name="Cell")
        assert lmp_input.stages == [
            {"stage_name": "Stage 1", "commands": [("units", "metal")]},
            {"stage_name": "Pair style", "commands": [("pair_style", "eam")]},
            {"stage_name": "Cell", "commands": [("boundary", "p p p"), ("atom_style", "full")]},
        ]

        lmp_input.add_stage(commands=["set type 2 charge 0.0"], after_stage="Stage 1")
        assert lmp_input.stages == [
            {"stage_name": "Stage 1", "commands": [("units", "metal")]},
            {"stage_name": "Stage 4", "commands": [("set", "type 2 charge 0.0")]},
            {"stage_name": "Pair style", "commands": [("pair_style", "eam")]},
            {"stage_name": "Cell", "commands": [("boundary", "p p p"), ("atom_style", "full")]},
        ]
        lmp_input.add_stage(stage_name="Empty stage")
        assert lmp_input.stages == [
            {"stage_name": "Stage 1", "commands": [("units", "metal")]},
            {"stage_name": "Stage 4", "commands": [("set", "type 2 charge 0.0")]},
            {"stage_name": "Pair style", "commands": [("pair_style", "eam")]},
            {"stage_name": "Cell", "commands": [("boundary", "p p p"), ("atom_style", "full")]},
            {"stage_name": "Empty stage", "commands": []},
        ]
        lmp_input.add_stage(
            {"stage_name": "New stage", "commands": [("newcommand", "newargs 1 2 3"), ("newcommand2", "newargs 4 5 6")]}
        )
        assert lmp_input.stages == [
            {"stage_name": "Stage 1", "commands": [("units", "metal")]},
            {"stage_name": "Stage 4", "commands": [("set", "type 2 charge 0.0")]},
            {"stage_name": "Pair style", "commands": [("pair_style", "eam")]},
            {"stage_name": "Cell", "commands": [("boundary", "p p p"), ("atom_style", "full")]},
            {"stage_name": "Empty stage", "commands": []},
            {
                "stage_name": "New stage",
                "commands": [("newcommand", "newargs 1 2 3"), ("newcommand2", "newargs 4 5 6")],
            },
        ]
        lmp_input.add_stage(
            {
                "stage_name": "New stage 2",
                "commands": [("newcommand", "newargs 1 2 3"), ("newcommand2", "newargs 4 5 6")],
            },
            after_stage="Stage 4",
        )
        assert lmp_input.stages == [
            {"stage_name": "Stage 1", "commands": [("units", "metal")]},
            {"stage_name": "Stage 4", "commands": [("set", "type 2 charge 0.0")]},
            {
                "stage_name": "New stage 2",
                "commands": [("newcommand", "newargs 1 2 3"), ("newcommand2", "newargs 4 5 6")],
            },
            {"stage_name": "Pair style", "commands": [("pair_style", "eam")]},
            {"stage_name": "Cell", "commands": [("boundary", "p p p"), ("atom_style", "full")]},
            {"stage_name": "Empty stage", "commands": []},
            {
                "stage_name": "New stage",
                "commands": [("newcommand", "newargs 1 2 3"), ("newcommand2", "newargs 4 5 6")],
            },
        ]

    def test_merge_stages(self):
        lmp_input = LammpsInputFile().from_file(self.filename, keep_stages=True)
        lmp_input.merge_stages(["1) Initialization", "3) Simulation settings"])
        assert lmp_input.stages_names == [
            "Comment 1",
            "Merge of: 1) Initialization, 3) Simulation settings",
            "2) System definition",
            "Part A : energy minimization",
            "Stage 6",
        ]
        assert lmp_input.stages == [
            {"stage_name": "Comment 1", "commands": [("#", "LGPS")]},
            {
                "stage_name": "Merge of: 1) Initialization, 3) Simulation settings",
                "commands": [
                    ("units", "metal"),
                    ("atom_style", "full"),
                    ("dimension", "3"),
                    ("pair_style", "hybrid/overlay morse 15 coul/long 15"),
                    ("kspace_style", "ewald 1e-4"),
                    ("boundary", "p p p"),
                    ("pair_coeff", "1 1 morse 0.0580 3.987 3.404"),
                    ("pair_coeff", "1 4 morse 0.0408 1.399 3.204"),
                    ("pair_coeff", "2 4 morse 0.3147 2.257 2.409"),
                    ("pair_coeff", "3 4 morse 0.4104 2.329 2.200"),
                    ("pair_coeff", "4 4 morse 0.0241 1.359 4.284"),
                    ("pair_coeff", "* * coul/long"),
                ],
            },
            {
                "stage_name": "2) System definition",
                "commands": [
                    ("read_data", "run_init.data"),
                    ("set", "type 1 charge 0.8803"),
                    ("set", "type 2 charge 1.2570"),
                    ("set", "type 3 charge 1.2580"),
                    ("set", "type 4 charge -1.048"),
                    ("neigh_modify", "every 1 delay 5 check yes"),
                ],
            },
            {
                "stage_name": "Part A : energy minimization",
                "commands": [
                    ("thermo", "1"),
                    ("thermo_style", "custom step lx ly lz press pxx pyy pzz pe"),
                    ("dump", "dmp all atom 5 run.dump"),
                ],
            },
            {
                "stage_name": "Stage 6",
                "commands": [
                    ("min_style", "cg"),
                    ("fix", "1 all box/relax iso 0.0 vmax 0.001"),
                    ("minimize", "1.0e-16 1.0e-16 5000 10000"),
                    ("write_data", "run.data"),
                ],
            },
        ]

        lmp_input = LammpsInputFile().from_file(self.filename, keep_stages=True)
        lmp_input.merge_stages(lmp_input.stages_names)
        assert len(lmp_input.stages_names) == 1
        assert lmp_input.stages == [
            {
                "stage_name": "Merge of: Comment 1, 1) Initialization, 2) System definition, 3) Simulation settings, "
                "Part A : "
                "energy minimization, Stage 6",
                "commands": [
                    ("#", "LGPS"),
                    ("units", "metal"),
                    ("atom_style", "full"),
                    ("dimension", "3"),
                    ("pair_style", "hybrid/overlay morse 15 coul/long 15"),
                    ("kspace_style", "ewald 1e-4"),
                    ("boundary", "p p p"),
                    ("read_data", "run_init.data"),
                    ("set", "type 1 charge 0.8803"),
                    ("set", "type 2 charge 1.2570"),
                    ("set", "type 3 charge 1.2580"),
                    ("set", "type 4 charge -1.048"),
                    ("neigh_modify", "every 1 delay 5 check yes"),
                    ("pair_coeff", "1 1 morse 0.0580 3.987 3.404"),
                    ("pair_coeff", "1 4 morse 0.0408 1.399 3.204"),
                    ("pair_coeff", "2 4 morse 0.3147 2.257 2.409"),
                    ("pair_coeff", "3 4 morse 0.4104 2.329 2.200"),
                    ("pair_coeff", "4 4 morse 0.0241 1.359 4.284"),
                    ("pair_coeff", "* * coul/long"),
                    ("thermo", "1"),
                    ("thermo_style", "custom step lx ly lz press pxx pyy pzz pe"),
                    ("dump", "dmp all atom 5 run.dump"),
                    ("min_style", "cg"),
                    ("fix", "1 all box/relax iso 0.0 vmax 0.001"),
                    ("minimize", "1.0e-16 1.0e-16 5000 10000"),
                    ("write_data", "run.data"),
                ],
            }
        ]

    def test_remove_stage(self):
        lmp_input = LammpsInputFile().from_file(self.filename)
        lmp_input.remove_stage("Stage 1")
        assert lmp_input.stages == []

        lmp_input = LammpsInputFile().from_file(self.filename, keep_stages=True)
        lmp_input.remove_stage("1) Initialization")
        assert lmp_input.stages_names == [
            "Comment 1",
            "2) System definition",
            "3) Simulation settings",
            "Part A : energy minimization",
            "Stage 6",
        ]

    def test_rename_stage(self):
        lmp_input = LammpsInputFile().from_file(self.filename)
        lmp_input.rename_stage("Stage 1", "Global stage")
        assert lmp_input.stages_names == ["Global stage"]

        lmp_input = LammpsInputFile().from_file(self.filename, keep_stages=True)
        lmp_input.rename_stage("Stage 6", "Final stage")
        assert lmp_input.stages_names == [
            "Comment 1",
            "1) Initialization",
            "2) System definition",
            "3) Simulation settings",
            "Part A : energy minimization",
            "Final stage",
        ]

    def test_remove_command(self):
        lmp_input = LammpsInputFile().from_file(self.filename)
        lmp_input.remove_command("set")
        assert lmp_input.stages == [
            {
                "stage_name": "Stage 1",
                "commands": [
                    ("units", "metal"),
                    ("atom_style", "full"),
                    ("dimension", "3"),
                    ("pair_style", "hybrid/overlay morse 15 coul/long 15"),
                    ("kspace_style", "ewald 1e-4"),
                    ("boundary", "p p p"),
                    ("#", "2) System definition"),
                    ("read_data", "run_init.data"),
                    ("neigh_modify", "every 1 delay 5 check yes"),
                    ("#", "3) Simulation settings"),
                    ("pair_coeff", "1 1 morse 0.0580 3.987 3.404"),
                    ("pair_coeff", "1 4 morse 0.0408 1.399 3.204"),
                    ("pair_coeff", "2 4 morse 0.3147 2.257 2.409"),
                    ("pair_coeff", "3 4 morse 0.4104 2.329 2.200"),
                    ("pair_coeff", "4 4 morse 0.0241 1.359 4.284"),
                    ("pair_coeff", "* * coul/long"),
                    ("#", "Part A : energy minimization"),
                    ("thermo", "1"),
                    ("thermo_style", "custom step lx ly lz press pxx pyy pzz pe"),
                    ("dump", "dmp all atom 5 run.dump"),
                    ("min_style", "cg"),
                    ("fix", "1 all box/relax iso 0.0 vmax 0.001"),
                    ("minimize", "1.0e-16 1.0e-16 5000 10000"),
                    ("write_data", "run.data"),
                ],
            }
        ]

        lmp_input = LammpsInputFile().from_file(self.filename, keep_stages=True)
        lmp_input._add_command(command="set type 1 charge 0.0", stage_name="1) Initialization")
        lmp_input.remove_command("set")
        assert lmp_input.stages == [
            {"stage_name": "Comment 1", "commands": [("#", "LGPS")]},
            {
                "stage_name": "1) Initialization",
                "commands": [
                    ("units", "metal"),
                    ("atom_style", "full"),
                    ("dimension", "3"),
                    ("pair_style", "hybrid/overlay morse 15 coul/long 15"),
                    ("kspace_style", "ewald 1e-4"),
                    ("boundary", "p p p"),
                ],
            },
            {
                "stage_name": "2) System definition",
                "commands": [("read_data", "run_init.data"), ("neigh_modify", "every 1 delay 5 check yes")],
            },
            {
                "stage_name": "3) Simulation settings",
                "commands": [
                    ("pair_coeff", "1 1 morse 0.0580 3.987 3.404"),
                    ("pair_coeff", "1 4 morse 0.0408 1.399 3.204"),
                    ("pair_coeff", "2 4 morse 0.3147 2.257 2.409"),
                    ("pair_coeff", "3 4 morse 0.4104 2.329 2.200"),
                    ("pair_coeff", "4 4 morse 0.0241 1.359 4.284"),
                    ("pair_coeff", "* * coul/long"),
                ],
            },
            {
                "stage_name": "Part A : energy minimization",
                "commands": [
                    ("thermo", "1"),
                    ("thermo_style", "custom step lx ly lz press pxx pyy pzz pe"),
                    ("dump", "dmp all atom 5 run.dump"),
                ],
            },
            {
                "stage_name": "Stage 6",
                "commands": [
                    ("min_style", "cg"),
                    ("fix", "1 all box/relax iso 0.0 vmax 0.001"),
                    ("minimize", "1.0e-16 1.0e-16 5000 10000"),
                    ("write_data", "run.data"),
                ],
            },
        ]

    def test_append(self):
        lmp_input = LammpsInputFile().from_file(self.filename)
        lmp_input2 = LammpsInputFile().from_file(self.filename)
        lmp_input.append(lmp_input2)
        assert lmp_input.nstages == 2
        assert lmp_input.ncomments == 6
        assert lmp_input.stages_names == ["Stage 1", "Stage 2"]

        lmp_input = LammpsInputFile().from_file(self.filename, keep_stages=True)
        lmp_input2 = LammpsInputFile().from_file(self.filename, keep_stages=True)
        lmp_input.append(lmp_input2)
        assert lmp_input.nstages == 12
        assert lmp_input.ncomments == 2
        assert lmp_input.stages_names == [
            "Comment 1",
            "1) Initialization",
            "2) System definition",
            "3) Simulation settings",
            "Part A : energy minimization",
            "Stage 6",
            "Comment 2",
            "Stage 8 (previously 1) Initialization)",
            "Stage 9 (previously 2) System definition)",
            "Stage 10 (previously 3) Simulation settings)",
            "Stage 11 (previously Part A : energy minimization)",
            "Stage 12",
        ]

    def test_add_command(self):
        lmp_input = LammpsInputFile()
        lmp_input.stages.append({"stage_name": "Init", "commands": [("units", "metal")]})
        lmp_input._add_command(stage_name="Init", command="boundary p p p")
        assert lmp_input.stages == [{"stage_name": "Init", "commands": [("units", "metal"), ("boundary", "p p p")]}]
        lmp_input._add_command(stage_name="Init", command="atom_style", args="full")
        assert lmp_input.stages == [
            {"stage_name": "Init", "commands": [("units", "metal"), ("boundary", "p p p"), ("atom_style", "full")]}
        ]

    def test_add_commands(self):
        lmp_input = LammpsInputFile()
        lmp_input.add_stage({"stage_name": "Init", "commands": [("units", "metal")]})
        lmp_input.add_commands(stage_name="Init", commands="boundary p p p")
        assert lmp_input.stages == [{"stage_name": "Init", "commands": [("units", "metal"), ("boundary", "p p p")]}]
        lmp_input.add_commands(stage_name="Init", commands=["atom_style full", "dimension 3"])
        assert lmp_input.stages == [
            {
                "stage_name": "Init",
                "commands": [("units", "metal"), ("boundary", "p p p"), ("atom_style", "full"), ("dimension", "3")],
            }
        ]
        lmp_input.add_commands(stage_name="Init", commands={"kspace_style": "ewald 1e-4"})
        assert lmp_input.stages == [
            {
                "stage_name": "Init",
                "commands": [
                    ("units", "metal"),
                    ("boundary", "p p p"),
                    ("atom_style", "full"),
                    ("dimension", "3"),
                    ("kspace_style", "ewald 1e-4"),
                ],
            }
        ]

    def test_add_comment(self):
        lmp_input = LammpsInputFile()
        lmp_input._add_comment(comment="First comment")
        assert lmp_input.stages == [{"stage_name": "Comment 1", "commands": [("#", "First comment")]}]
        lmp_input._add_comment(comment="Sub comment", stage_name="Comment 1")
        assert lmp_input.stages == [
            {"stage_name": "Comment 1", "commands": [("#", "First comment"), ("#", "Sub comment")]}
        ]
        lmp_input.stages.append({"stage_name": "Init", "commands": [("units", "metal")]})
        lmp_input._add_comment(comment="Inline comment", stage_name="Init")
        assert lmp_input.stages == [
            {"stage_name": "Comment 1", "commands": [("#", "First comment"), ("#", "Sub comment")]},
            {"stage_name": "Init", "commands": [("units", "metal"), ("#", "Inline comment")]},
        ]


class TestLammpsRun(unittest.TestCase):
    maxDiff = None

    def test_md(self):
        struct = Structure.from_spacegroup(225, Lattice.cubic(3.62126), ["Cu"], [[0, 0, 0]])
        ld = LammpsData.from_structure(struct, atom_style="atomic")
        ff = "pair_style eam\npair_coeff * * Cu_u3.eam"
        md = LammpsRun.md(data=ld, force_field=ff, temperature=1600.0, nsteps=10000)
        md.write_inputs(output_dir="md")
        with open(os.path.join("md", "in.md")) as f:
            md_script = f.read()
        script_string = """# Sample input script template for MD

# Initialization

units           metal
atom_style      atomic

# Atom definition

read_data       md.data
#read_restart    md.restart

# Force field settings (consult official document for detailed formats)

pair_style eam
pair_coeff * * Cu_u3.eam

# Create velocities
velocity        all create 1600.0 142857 mom yes rot yes dist gaussian

# Ensemble constraints
#fix             1 all nve
fix             1 all nvt temp 1600.0 1600.0 0.1
#fix             1 all npt temp 1600.0 1600.0 0.1 iso $pressure $pressure 1.0

# Various operations within timestepping
#fix             ...
#compute         ...

# Output settings
#thermo_style    custom ...  # control the thermo data type to output
thermo          100  # output thermo data every N steps
#dump            1 all atom 100 traj.*.gz  # dump a snapshot every 100 steps

# Actions
run             10000
"""
        assert md_script == script_string
        assert os.path.exists(os.path.join("md", "md.data"))

    @classmethod
    def tearDownClass(cls):
        temp_dirs = ["md"]
        for td in temp_dirs:
            if os.path.exists(td):
                shutil.rmtree(td)


class TestFunc(unittest.TestCase):
    def test_write_lammps_inputs(self):
        # script template
        with open(f"{TEST_FILES_DIR}/lammps/kappa.txt") as f:
            kappa_template = f.read()
        kappa_settings = {"method": "heat"}
        write_lammps_inputs(output_dir="heat", script_template=kappa_template, settings=kappa_settings)
        with open(os.path.join("heat", "in.lammps")) as f:
            kappa_script = f.read()
        fix_hot = re.search(r"fix\s+hot\s+all\s+([^\s]+)\s+", kappa_script)
        # placeholders supposed to be filled
        assert fix_hot.group(1) == "heat"
        fix_cold = re.search(r"fix\s+cold\s+all\s+([^\s]+)\s+", kappa_script)
        assert fix_cold.group(1) == "heat"
        lattice = re.search(r"lattice\s+fcc\s+(.*)\n", kappa_script)
        # parentheses not supposed to be filled
        assert lattice.group(1) == "${rho}"
        pair_style = re.search(r"pair_style\slj/cut\s+(.*)\n", kappa_script)
        assert pair_style.group(1) == "${rc}"

        with open(f"{TEST_FILES_DIR}/lammps/in.peptide") as f:
            peptide_script = f.read()
        # copy data file
        src = f"{TEST_FILES_DIR}/lammps/data.quartz"
        write_lammps_inputs(output_dir="path", script_template=peptide_script, data=src)
        dst = os.path.join("path", "data.peptide")
        assert filecmp.cmp(src, dst, shallow=False)
        # write data file from obj
        obj = LammpsData.from_file(src, atom_style="atomic")
        with pytest.warns(FutureWarning):
            write_lammps_inputs(output_dir="obj", script_template=peptide_script, data=obj)
            obj_read = LammpsData.from_file(os.path.join("obj", "data.peptide"), atom_style="atomic")
            pd.testing.assert_frame_equal(obj_read.masses, obj.masses)
            pd.testing.assert_frame_equal(obj_read.atoms, obj.atoms)

    @classmethod
    def tearDownClass(cls):
        temp_dirs = ["heat", "path", "obj"]
        for td in temp_dirs:
            if os.path.exists(td):
                shutil.rmtree(td)


class TestLammpsTemplateGen(PymatgenTest):
    def test_write_inputs(self):
        # simple script without data file
        lis = LammpsTemplateGen().get_input_set(
            script_template=f"{TEST_FILES_DIR}/lammps/kappa.txt",
            settings={"method": "heat"},
            data=None,
            data_filename="data.peptide",
        )
        assert len(lis) == 1
        lis.write_input(self.tmp_path / "heat")

        with open(self.tmp_path / "heat" / "in.lammps") as f:
            kappa_script = f.read()
        fix_hot = re.search(r"fix\s+hot\s+all\s+([^\s]+)\s+", kappa_script)
        # placeholders supposed to be filled
        assert fix_hot.group(1) == "heat"
        fix_cold = re.search(r"fix\s+cold\s+all\s+([^\s]+)\s+", kappa_script)
        assert fix_cold.group(1) == "heat"
        lattice = re.search(r"lattice\s+fcc\s+(.*)\n", kappa_script)
        # parentheses not supposed to be filled
        assert lattice.group(1) == "${rho}"
        pair_style = re.search(r"pair_style\slj/cut\s+(.*)\n", kappa_script)
        assert pair_style.group(1) == "${rc}"

        # script with data file
        obj = LammpsData.from_file(f"{TEST_FILES_DIR}/lammps/data.quartz", atom_style="atomic")
        lis = LammpsTemplateGen().get_input_set(
            script_template=f"{TEST_FILES_DIR}/lammps/in.peptide",
            settings=None,
            data=obj,
            data_filename="data.peptide",
        )
        assert len(lis) == 2
        assert isinstance(lis["data.peptide"], LammpsData)
        lis.write_input(self.tmp_path / "obj")

        obj_read = LammpsData.from_file(str(self.tmp_path / "obj" / "data.peptide"), atom_style="atomic")
        pd.testing.assert_frame_equal(obj_read.masses, obj.masses)
        pd.testing.assert_frame_equal(obj_read.atoms, obj.atoms)
