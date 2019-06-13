# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import unittest

from monty.os.path import which

import pymatgen.command_line.vampire_caller as vc

@unittest.skipIf(not which('vampire-serial'), "vampire executable not present")
class VampireCallerTest(unittest.TestCase):
    """Todo:
        * Put json files of magnetic orderings and energies in
        pymatgen/test_files with variable number of unique magnetic
        sites, sublattices, FM, FiM, and AFM ground state, etc.
        * Compare Tc calcs to known Tc
    """