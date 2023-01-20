# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
from __future__ import annotations

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
from pymatgen.io.lammps.inputs import LammpsRun, LammpsTemplateGen, write_lammps_inputs
from pymatgen.util.testing import PymatgenTest


class LammpsRunTest(unittest.TestCase):
    maxDiff = None

    def test_md(self):
        s = Structure.from_spacegroup(225, Lattice.cubic(3.62126), ["Cu"], [[0, 0, 0]])
        ld = LammpsData.from_structure(s, atom_style="atomic")
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
        assert fix_hot.group(1) == "heat"
        fix_cold = re.search(r"fix\s+cold\s+all\s+([^\s]+)\s+", kappa_script)
        assert fix_cold.group(1) == "heat"
        lattice = re.search(r"lattice\s+fcc\s+(.*)\n", kappa_script)
        # parentheses not supposed to be filled
        assert lattice.group(1) == "${rho}"
        pair_style = re.search(r"pair_style\slj/cut\s+(.*)\n", kappa_script)
        assert pair_style.group(1) == "${rc}"

        with open(os.path.join(PymatgenTest.TEST_FILES_DIR, "lammps", "in.peptide")) as f:
            peptide_script = f.read()
        # copy data file
        src = os.path.join(PymatgenTest.TEST_FILES_DIR, "lammps", "data.quartz")
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
