# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import os
import unittest
import warnings

import pandas as pd
from monty.os.path import which

import pymatgen.command_line.vampire_caller as vampirecaller
from pymatgen.core.structure import Structure
from pymatgen.util.testing import PymatgenTest

test_dir = os.path.join(PymatgenTest.TEST_FILES_DIR, "magnetic_orderings")


@unittest.skipIf(not which("vampire-serial"), "vampire executable not present")
class VampireCallerTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):

        print("Testing with: ", which("vampire-serial"))

        cls.Mn3Al = pd.read_json(os.path.join(test_dir, "Mn3Al.json"))

        cls.compounds = [cls.Mn3Al]

        cls.structure_inputs = []
        cls.energy_inputs = []
        for c in cls.compounds:
            ordered_structures = list(c["structure"])
            ordered_structures = [Structure.from_dict(d) for d in ordered_structures]
            epa = list(c["energy_per_atom"])
            energies = [e * len(s) for (e, s) in zip(epa, ordered_structures)]

            cls.structure_inputs.append(ordered_structures)
            cls.energy_inputs.append(energies)

    def setUp(self):
        pass

    def tearDown(self):
        warnings.simplefilter("default")

    def test_vampire(self):
        for structs, energies in zip(self.structure_inputs, self.energy_inputs):
            settings = {"start_t": 0, "end_t": 500, "temp_increment": 50}
            vc = vampirecaller.VampireCaller(
                structs,
                energies,
                mc_box_size=3.0,
                equil_timesteps=1000,  # 1000
                mc_timesteps=2000,  # 2000
                user_input_settings=settings,
            )

            voutput = vc.output
            critical_temp = voutput.critical_temp
            self.assertAlmostEqual(400, critical_temp, delta=100)

        if os.path.exists("Mn3Al.mat"):
            os.remove("Mn3Al.mat")
        if os.path.exists("Mn3Al.ucf"):
            os.remove("Mn3Al.ucf")
        if os.path.exists("input"):
            os.remove("input")
        if os.path.exists("log"):
            os.remove("log")
        if os.path.exists("output"):
            os.remove("output")


if __name__ == "__main__":
    unittest.main()
