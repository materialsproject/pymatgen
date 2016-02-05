#!/usr/bin/env python

import os
import sys
import string
import json

import inspect

from pymatgen import Structure
from fireworks import Firework, Workflow, LaunchPad
from pymatgen.io.vasp.interfaces import VaspInputInterface, VaspFirework, VaspWorkflow

# get structure from Crystallographic Information File (CIF)
s = Structure.from_file('./mp-33088_Cr2FeO4.cif')

input=VaspInputInterface(s)
input.NEDOS=2000 # override default or add INCAR parameter

# Dump VASP Input into current directory for inspection
input.write_input()

# Complete definition of Firework Task(s) and add to
# Firework
task=VaspFirework(input)

# Save specification to yaml file for later inspection
# or manual add to launchpad with lpad script
task.to_file("simple_task.yaml")

# Adds single Firework to launchpad database
task.add_fw_to_launchpad()

