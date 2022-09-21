# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
import filecmp
import os
import re
import shutil
import tempfile
import unittest
from pathlib import Path

import pandas as pd
import pytest

from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.io.lammps.data import LammpsData
from pymatgen.io.lammps.inputs import (
    LammpsInputFile,
    LammpsRun,
    LammpsTemplateGen,
    write_lammps_inputs,
)
from pymatgen.util.testing import PymatgenTest

test_dir = os.path.join(PymatgenTest.TEST_FILES_DIR, "lammps")


class LammpsInputFileTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.filename = os.path.join(test_dir, "lgps.in")

    def test_from_file(self):
        lmp_input = LammpsInputFile().from_file(self.filename)
        self.assertEqual(lmp_input.ncomments, 1)
        self.assertEqual(lmp_input.nstages, 5)
        self.assertListEqual(
            list(lmp_input.input_settings.keys()),
            [
                "Comment 1",
                "1) Initialization ",
                "2) System definition ",
                "3) Simulation settings ",
                "Part A : energy minimization ",
                "Stage 5",
            ],
        )
        self.assertDictEqual(list(lmp_input.input_settings.values())[0], {"#": "LGPS"})
        self.assertDictEqual(
            list(lmp_input.input_settings.values())[1],
            {
                "units": "metal",
                "atom_style": "full",
                "dimension": "3",
                "pair_style": "hybrid/overlay morse 15 coul/long 15",
                "kspace_style": "ewald 1e-4",
                "boundary": "p p p",
            },
        )

    def test_get_string(self):
        lmp_input = LammpsInputFile().from_file(self.filename)
        string = lmp_input.get_string()
        self.assertIn("# LAMMPS input generated from LammpsInputFile\n\n# LGPS\n\n# 1) Initialization", string)
        self.assertIn(
            "\nunits metal\natom_style full\ndimension 3\npair_style hybrid/overlay morse 15 coul/long", string
        )
        self.assertIn("15\nkspace_style ewald 1e-4\nboundary p p p\n\n# 2) System definition", string)
        self.assertIn("\nread_data run_init.data\nset type 1 charge 0.8803\nset type 2 charge", string)
        self.assertIn("1.2570\nset type 3 charge 1.2580\nset type 4 charge -1.048\nneigh_modify", string)
        self.assertIn("every 1 delay 5 check yes\n\n# 3) Simulation settings \npair_coeff 1 1", string)
        self.assertIn("morse 0.0580 3.987 3.404\npair_coeff 1 4 morse 0.0408 1.399 3.204\npair_coeff", string)
        self.assertIn("2 4 morse 0.3147 2.257 2.409\npair_coeff 3 4 morse 0.4104 2.329 2.200\npair_coeff", string)
        self.assertIn("4 4 morse 0.0241 1.359 4.284\npair_coeff * * coul/long\n\n# Part A : energy", string)
        self.assertIn("minimization \nthermo 1\nthermo_style custom step lx ly lz press pxx pyy pzz", string)
        self.assertIn("pe\ndump dmp all atom 5 run.dump\n\n# Stage 5\nmin_style cg\nfix 1 all", string)
        self.assertIn(
            "box/relax iso 0.0 vmax 0.001\nminimize 1.0e-16 1.0e-16 5000 10000\nwrite_data run.data\n", string
        )

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

        lmp_input = LammpsInputFile().from_string(string)
        self.assertDictEqual(lmp_input.input_settings, LammpsInputFile().from_file(self.filename).input_settings)

    def test_add_stage(self):
        lmp_input = LammpsInputFile()
        lmp_input.add_stage(command="units metal")
        self.assertDictEqual(lmp_input.input_settings, {"Stage 1": {"units": "metal"}})
        lmp_input.add_stage(command={"pair_style": "eam"}, stage_name="Pair style")
        self.assertDictEqual(
            lmp_input.input_settings, {"Stage 1": {"units": "metal"}, "Pair style": {"pair_style": "eam"}}
        )
        lmp_input.add_stage(command=["boundary p p p", "atom_style full"], stage_name="Cell")
        self.assertDictEqual(
            lmp_input.input_settings,
            {
                "Stage 1": {"units": "metal"},
                "Pair style": {"pair_style": "eam"},
                "Cell": {"boundary": "p p p", "atom_style": "full"},
            },
        )

    def test_add_command(self):
        lmp_input = LammpsInputFile()
        lmp_input.add_stage(command="units metal", stage_name="Init")
        lmp_input.add_command(stage_name="Init", command="boundary p p p")
        self.assertDictEqual(lmp_input.input_settings, {"Init": {"units": "metal", "boundary": "p p p"}})
        lmp_input.add_command(stage_name="Init", command="atom_style", args="full")
        self.assertDictEqual(
            lmp_input.input_settings, {"Init": {"units": "metal", "boundary": "p p p", "atom_style": "full"}}
        )

    def test_add_comment(self):
        lmp_input = LammpsInputFile()
        lmp_input.add_comment(comment="First comment")
        self.assertDictEqual(lmp_input.input_settings, {"Comment 1": {"#": "First comment"}})
        lmp_input.add_comment(comment="Sub comment", stage_name="Comment 1")
        self.assertDictEqual(lmp_input.input_settings, {"Comment 1": {"#": "First comment\n# Sub comment"}})
        lmp_input.add_stage(command="units metal", stage_name="Init")
        lmp_input.add_comment(comment="Inline comment", stage_name="Init")
        self.assertDictEqual(
            lmp_input.input_settings,
            {"Comment 1": {"#": "First comment\n# Sub comment"}, "Init": {"units": "metal", "#": "Inline comment"}},
        )


class LammpsRunTest(unittest.TestCase):
    maxDiff = None

    def test_md(self):
        s = Structure.from_spacegroup(225, Lattice.cubic(3.62126), ["Cu"], [[0, 0, 0]])
        ld = LammpsData.from_structure(s, atom_style="atomic")
        ff = "\n".join(["pair_style eam", "pair_coeff * * Cu_u3.eam"])
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
        self.assertEqual(md_script, script_string)
        self.assertTrue(os.path.exists(os.path.join("md", "md.data")))

    @classmethod
    def tearDownClass(cls):
        temp_dirs = ["md"]
        for td in temp_dirs:
            if os.path.exists(td):
                shutil.rmtree(td)


class FuncTest(unittest.TestCase):
    def test_write_lammps_inputs(self):
        # script template
        with open(os.path.join(PymatgenTest.TEST_FILES_DIR, "lammps", "kappa.txt")) as f:
            kappa_template = f.read()
        kappa_settings = {"method": "heat"}
        write_lammps_inputs(output_dir="heat", script_template=kappa_template, settings=kappa_settings)
        with open(os.path.join("heat", "in.lammps")) as f:
            kappa_script = f.read()
        fix_hot = re.search(r"fix\s+hot\s+all\s+([^\s]+)\s+", kappa_script)
        # placeholders supposed to be filled
        self.assertEqual(fix_hot.group(1), "heat")
        fix_cold = re.search(r"fix\s+cold\s+all\s+([^\s]+)\s+", kappa_script)
        self.assertEqual(fix_cold.group(1), "heat")
        lattice = re.search(r"lattice\s+fcc\s+(.*)\n", kappa_script)
        # parentheses not supposed to be filled
        self.assertEqual(lattice.group(1), "${rho}")
        pair_style = re.search(r"pair_style\slj/cut\s+(.*)\n", kappa_script)
        self.assertEqual(pair_style.group(1), "${rc}")

        with open(os.path.join(PymatgenTest.TEST_FILES_DIR, "lammps", "in.peptide")) as f:
            peptide_script = f.read()
        # copy data file
        src = os.path.join(PymatgenTest.TEST_FILES_DIR, "lammps", "data.quartz")
        write_lammps_inputs(output_dir="path", script_template=peptide_script, data=src)
        dst = os.path.join("path", "data.peptide")
        self.assertTrue(filecmp.cmp(src, dst, shallow=False))
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


class LammpsTemplateGenTest(PymatgenTest):
    def test_write_inputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            # simple script without data file
            lis = LammpsTemplateGen().get_input_set(
                script_template=os.path.join(PymatgenTest.TEST_FILES_DIR, "lammps", "kappa.txt"),
                settings={"method": "heat"},
                data=None,
                data_filename="data.peptide",
            )
            tmpdir = Path(tmpdir)
            assert len(lis) == 1
            lis.write_input(tmpdir / "heat")

            with open(tmpdir / "heat" / "in.lammps") as f:
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
            obj = LammpsData.from_file(
                os.path.join(PymatgenTest.TEST_FILES_DIR, "lammps", "data.quartz"), atom_style="atomic"
            )
            lis = LammpsTemplateGen().get_input_set(
                script_template=os.path.join(PymatgenTest.TEST_FILES_DIR, "lammps", "in.peptide"),
                settings=None,
                data=obj,
                data_filename="data.peptide",
            )
            assert len(lis) == 2
            assert isinstance(lis["data.peptide"], LammpsData)
            lis.write_input(tmpdir / "obj")

            obj_read = LammpsData.from_file(str(tmpdir / "obj" / "data.peptide"), atom_style="atomic")
            pd.testing.assert_frame_equal(obj_read.masses, obj.masses)
            pd.testing.assert_frame_equal(obj_read.atoms, obj.atoms)


if __name__ == "__main__":
    unittest.main()
